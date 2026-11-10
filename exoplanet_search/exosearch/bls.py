"""Box Least Squares search with selectable backend.

cpu   : astropy BoxLeastSquares (C, fast method)  -- tested, default on CPU workers
numpy : vectorized reference BLS (same algorithm as cuda/mlx) -- tested, for validation
cuda  : vectorized BLS on CuPy (Pascal sm_61 OK)  -- UNTESTED stub
mlx   : vectorized BLS on Apple MLX (mlx.core)     -- UNTESTED stub

All backends return the same dict; vetting stats are always recomputed on CPU with
astropy's compute_stats at the best (period, duration, t0) so results are comparable.
"""
from __future__ import annotations

import numpy as np
from astropy.timeseries import BoxLeastSquares
from scipy.ndimage import median_filter

DURATIONS = np.array([0.04, 0.06, 0.08, 0.10, 0.13, 0.17, 0.21, 0.25])  # days, ~1 h .. 6 h


def period_grid(t, pmin=0.5, pmax=None, min_transits=2, frequency_factor=1.0,
                durations=DURATIONS):
    baseline = float(t.max() - t.min())
    pmax = pmax or baseline / min_transits
    df = frequency_factor * durations.min() / baseline ** 2
    freqs = np.arange(1.0 / pmax, 1.0 / pmin, df)
    return np.sort(1.0 / freqs)


def sde(power, window=501):
    """Signal Detection Efficiency on a periodogram with its slow trend removed."""
    p = np.asarray(power, float)
    if len(p) > window:
        p = p - median_filter(p, size=window, mode="nearest")
    s = np.std(p)
    return float((p.max() - p.mean()) / s) if s > 0 else 0.0, p


def search(t, f, e, backend="cpu", periods=None, durations=DURATIONS, **kw):
    if periods is None:
        periods = period_grid(t, durations=durations)
    if backend == "cpu":
        bls = BoxLeastSquares(t, f, dy=e)
        res = bls.power(periods, durations, objective="likelihood", method="fast")
        power = np.nan_to_num(np.asarray(res.power), nan=0.0)
        i = int(np.argmax(power))
        out = dict(period=float(res.period[i]), t0=float(res.transit_time[i]),
                   duration=float(res.duration[i]))
    else:
        from . import bls_vector
        out, power = bls_vector.search(t, f, e, periods, durations, backend=backend, **kw)
    s, pdet = sde(power)
    out.update(periods=periods, power=power, power_detrended=pdet, sde=s, backend=backend)
    return out


def stats(t, f, e, period, duration, t0):
    bls = BoxLeastSquares(t, f, dy=e)
    return bls.compute_stats(period, duration, t0)
