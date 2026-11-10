"""Vectorized BLS over many (period, duration) pairs for array libraries.

Algorithm (Kovacs et al. 2002, weighted): for a chunk of trial periods, phase-fold all points
at once -> (k, n), histogram weights w and w*y into nb phase bins with a scatter-add,
cumulative-sum (with wrap-around padding) and evaluate every box start and duration by
cumsum differences. With y mean-subtracted and w normalised to sum 1:
    r = sum_in w,  s = sum_in w*y,  depth = -s / (r (1 - r)),
    delta_chi2 = W_tot * s^2 / (r (1 - r))            (only s < 0, i.e. dips)
which matches astropy's "likelihood" objective up to a constant.

numpy is the tested reference; cupy (cuda) and mlx paths share the same code but are
UNTESTED here (no GPU / Apple Silicon on the CPU coordinator).
"""
from __future__ import annotations

import numpy as np


class _Ops:
    """Tiny adapter so one implementation runs on numpy, cupy and mlx.core."""

    def __init__(self, backend):
        self.backend = backend
        if backend == "numpy":
            self.xp = np
        elif backend == "cuda":
            import cupy as cp  # pip install cupy-cuda12x (Pascal sm_61 supported by CUDA 12)
            self.xp = cp
        elif backend == "mlx":
            import mlx.core as mx  # pip install mlx (Apple Silicon only)
            self.xp = mx
        else:
            raise ValueError(f"unknown vector backend {backend!r}")
        self.dtype = self.xp.float32 if backend == "mlx" else self.xp.float64

    def asarray(self, a, dtype=None):
        return self.xp.array(np.asarray(a), dtype=dtype or self.dtype)

    def scatter_add(self, idx, vals, size):
        xp = self.xp
        if self.backend == "mlx":
            return xp.zeros((size,), dtype=self.dtype).at[idx].add(vals)
        return xp.bincount(idx, weights=vals, minlength=size)

    def to_numpy(self, a):
        if self.backend == "cuda":
            return self.xp.asnumpy(a)
        return np.array(a)

    def take(self, a, idx):
        return self.xp.take_along_axis(a, idx, axis=1)


def _chunk_power(ops, t, y, w, wtot, P, dur, nb):
    xp = ops.xp
    k = P.shape[0]
    phase = xp.remainder(t[None, :] / P[:, None], 1.0)               # (k, n)
    b = xp.floor(phase * nb).astype(xp.int32)
    b = xp.minimum(b, nb - 1)
    flat = (xp.arange(k, dtype=xp.int32)[:, None] * nb + b).reshape(-1)
    W = ops.scatter_add(flat, xp.broadcast_to(w[None, :], (k, t.shape[0])).reshape(-1), k * nb).reshape(k, nb)
    Y = ops.scatter_add(flat, xp.broadcast_to((w * y)[None, :], (k, t.shape[0])).reshape(-1), k * nb).reshape(k, nb)
    # wrap-around padding then cumsum with a leading zero column
    W = xp.concatenate([xp.zeros((k, 1), dtype=W.dtype), W, W], axis=1)
    Y = xp.concatenate([xp.zeros((k, 1), dtype=Y.dtype), Y, Y], axis=1)
    CW = xp.cumsum(W, axis=1)
    CY = xp.cumsum(Y, axis=1)
    start = xp.arange(nb, dtype=xp.int32)[None, :]
    best = xp.zeros((k,), dtype=ops.dtype)
    best_j = xp.zeros((k,), dtype=xp.int32)
    best_d = xp.zeros((k,), dtype=xp.int32)
    for di in range(dur.shape[0]):
        nbin = xp.clip(xp.round(dur[di] / P * nb).astype(xp.int32), 1, nb // 2)    # (k,)
        end = start + nbin[:, None]
        r = ops.take(CW, end) - CW[:, :nb]
        s = ops.take(CY, end) - CY[:, :nb]
        denom = r * (1.0 - r)
        ok = (r > 1e-6) & (r < 1 - 1e-6) & (s < 0)
        pw = xp.where(ok, wtot * s * s / xp.where(ok, denom, 1.0), 0.0)
        j = xp.argmax(pw, axis=1).astype(xp.int32)
        m = xp.max(pw, axis=1)
        better = m > best
        best = xp.where(better, m, best)
        best_j = xp.where(better, j, best_j)
        best_d = xp.where(better, di, best_d)
    return best, best_j, best_d


def search(t, f, e, periods, durations, backend="numpy", chunk=128, max_bins=4096):
    ops = _Ops(backend)
    t_ref = float(t.min())
    tt = np.asarray(t, float) - t_ref          # float64 on host first: keeps float32 phase precise
    w = 1.0 / np.asarray(e, float) ** 2
    wtot = w.sum()
    w = w / wtot
    y = np.asarray(f, float) - np.sum(w * f)
    tx, yx, wx = ops.asarray(tt), ops.asarray(y), ops.asarray(w)
    dx = ops.asarray(durations)
    power = np.zeros(len(periods))
    jj = np.zeros(len(periods), int)
    dd = np.zeros(len(periods), int)
    nbs = np.zeros(len(periods), int)
    for a in range(0, len(periods), chunk):
        Pc = np.asarray(periods[a:a + chunk], float)
        nb = int(min(max_bins, max(64, np.ceil(2 * Pc.max() / durations.min()))))
        p, j, d = _chunk_power(ops, tx, yx, wx, wtot, ops.asarray(Pc), dx, nb)
        power[a:a + len(Pc)] = ops.to_numpy(p)
        jj[a:a + len(Pc)] = ops.to_numpy(j)
        dd[a:a + len(Pc)] = ops.to_numpy(d)
        nbs[a:a + len(Pc)] = nb
    i = int(np.argmax(power))
    P = float(periods[i])
    dur = float(durations[dd[i]])
    nbin = max(1, int(round(dur / P * nbs[i])))
    t0 = t_ref + (jj[i] + nbin / 2.0) / nbs[i] * P   # box centre
    return dict(period=P, t0=float(t0), duration=dur), power
