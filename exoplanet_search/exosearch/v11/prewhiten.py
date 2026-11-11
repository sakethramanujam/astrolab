"""Pulsator prewhitening above fmin (default 5 c/d) with the comb guard (docs/v1.1_plan.md item 4).

Iterative sinusoid removal on the quality-masked, normalised flux before detrending:
  1. Lomb-Scargle amplitude spectrum of the residual on a fixed grid fmin + df * arange(n) up to the Nyquist
     frequency, df = 1 / (oversample * baseline).
  2. Highest peak f (comb-protected frequencies excluded), refined by parabolic interpolation; amplitude from a
     linear least-squares fit of (sin, cos, offset) at f.
  3. Stop when amplitude S/N < snr_stop (noise = mean amplitude within +-noise_window of f, the peak's own bins
     excluded) or after max_terms removals.
  4. Comb guard (RF Engineer): before removing f, fit the amplitudes at 2f and 3f. If both have S/N >= comb_min_snr
     and the ratios A(2f)/A(f), A(3f)/A(f) match the harmonic ladder of a periodic box within a factor comb_ratio_tol
     for some duration 1-6 h at P = m / f (m = 1, 2, 3), the peak is a transit comb: it is logged as `comb`, not
     removed, and the whole harmonic series k * f / m of every matching fundamental is excluded from later picks
     (a box at period P has power at every k / P, so protecting only f, 2f, 3f would let the other harmonics be removed);
     each protected line covers +-2 / baseline (main lobe and first sidelobes of the spectral window).
Box ladder: harmonic j of a box with duty cycle q has amplitude proportional to |sin(pi j q)| / j, so with f the m-th
harmonic, A(2f)/A(f) = |cos(pi m q)| and A(3f)/A(f) = |sin(3 pi m q)| / (3 |sin(pi m q)|), q = D f / m.
"""
from __future__ import annotations

import json

import numpy as np


def _fit(t, r, f):
    """Linear LSQ of r on (sin, cos, 1) at frequency f -> (amplitude, a_sin, b_cos, offset)."""
    w = 2 * np.pi * f * t
    X = np.stack([np.sin(w), np.cos(w), np.ones_like(t)], axis=1)
    coef = np.linalg.solve(X.T @ X, X.T @ r)
    return float(np.hypot(coef[0], coef[1])), coef


def _amp_spectrum(t, r, freqs):
    from astropy.timeseries import LombScargle
    p = LombScargle(t, r, fit_mean=True, center_data=True, normalization="psd").power(freqs, method="fast",
                                                                                         assume_regular_frequency=True)
    return np.sqrt(np.maximum(4.0 * p / t.size, 0.0))


def _noise(amp, freqs, i, window, df):
    m = (np.abs(freqs - freqs[i]) <= window) & (np.abs(freqs - freqs[i]) > 3 * df)
    return float(np.mean(amp[m])) if m.any() else np.inf


def box_ladder(m, q):
    s1 = np.abs(np.sin(np.pi * m * q))
    r2 = np.abs(np.cos(np.pi * m * q))
    r3 = np.abs(np.sin(3 * np.pi * m * q)) / (3 * np.maximum(s1, 1e-12))
    return r2, r3


def comb_multiples(f, a1, a2, a3, durations=(1 / 24, 6 / 24), multiples=(1, 2, 3), tol=2.0, n=61):
    """Multiples m for which (A(f), A(2f), A(3f)) match a box ladder with f the m-th harmonic of 1/P."""
    if not (a1 > 0 and a2 > 0 and a3 > 0):
        return []
    o2, o3 = a2 / a1, a3 / a1
    lt = np.log(tol)
    out = []
    for m in multiples:
        D = np.linspace(durations[0], durations[1], n)
        q = D * f / m
        q = q[q < 0.5]
        if not q.size:
            continue
        r2, r3 = box_ladder(m, q)
        ok = (np.abs(np.log(o2 / np.maximum(r2, 1e-3))) < lt) & (np.abs(np.log(o3 / np.maximum(r3, 1e-3))) < lt)
        if ok.any():
            out.append(m)
    return out


def is_comb(f, a1, a2, a3, durations=(1 / 24, 6 / 24), multiples=(1, 2, 3), tol=2.0, n=61):
    return bool(comb_multiples(f, a1, a2, a3, durations, multiples, tol, n))


def prewhiten(t, f, cfg, comb_guard=None):
    """-> (prewhitened flux, log dict). cfg = configs/v1.1.json['prewhiten']."""
    t = np.asarray(t, float)
    r = np.asarray(f, float).copy()
    guard = cfg["comb_guard"] if comb_guard is None else comb_guard
    cad = float(np.median(np.diff(t)))
    base = float(t.max() - t.min())
    fnyq = 0.5 / cad
    df = 1.0 / (cfg["oversample"] * base)
    fmin = float(cfg["fmin_cpd"])
    if fnyq <= fmin + 10 * df or t.size < 100:
        return r, dict(prewhiten_n=0, prewhiten_comb=0, prewhiten_terms="[]", prewhiten_comb_freqs="[]")
    n = int(np.floor((fnyq - fmin) / df))
    freqs = fmin + df * np.arange(n, dtype=float)
    win = float(cfg["noise_window_cpd"])
    terms, combs, excl = [], [], np.zeros(n, bool)
    for _ in range(int(cfg["max_terms"]) + 3 * int(cfg["max_terms"])):      # comb picks do not count as terms
        if len(terms) >= cfg["max_terms"]:
            break
        amp = _amp_spectrum(t, r, freqs)
        a_m = np.where(excl, -np.inf, amp)
        i = int(np.argmax(a_m))
        if not np.isfinite(a_m[i]):
            break
        fi = freqs[i]
        if 0 < i < n - 1 and np.isfinite(a_m[i - 1]) and np.isfinite(a_m[i + 1]):
            y0, y1, y2 = amp[i - 1] ** 2, amp[i] ** 2, amp[i + 1] ** 2
            den = y0 - 2 * y1 + y2
            if den < 0:
                fi = fi + df * 0.5 * (y0 - y2) / den
        a1, coef = _fit(t, r, fi)
        noise = _noise(amp, freqs, i, win, df)
        snr = a1 / noise if noise > 0 else np.inf
        if snr < cfg["snr_stop"]:
            break
        if guard:
            a2 = a3 = 0.0
            s2 = s3 = 0.0
            if 3 * fi < fnyq:
                a2, _ = _fit(t, r, 2 * fi)
                a3, _ = _fit(t, r, 3 * fi)
                i2 = min(int(round((2 * fi - fmin) / df)), n - 1)
                i3 = min(int(round((3 * fi - fmin) / df)), n - 1)
                s2, s3 = a2 / _noise(amp, freqs, i2, win, df), a3 / _noise(amp, freqs, i3, win, df)
            ms = comb_multiples(fi, a1, a2, a3, tuple(cfg["comb_durations_d"]), tuple(cfg["comb_multiples"]),
                                cfg["comb_ratio_tol"]) if (s2 >= cfg["comb_min_snr"] and s3 >= cfg["comb_min_snr"]) else []
            if ms:
                combs.append([float(fi), a1, float(a2), float(a3), ms])
                for m in ms:              # protect the whole harmonic series of every matching fundamental f / m
                    f0 = fi / m
                    for h in range(1, int(fnyq / f0) + 1):
                        excl |= np.abs(freqs - h * f0) <= 2.0 / base       # main lobe and first sidelobes
                continue
        w = 2 * np.pi * fi * t
        r = r - (coef[0] * np.sin(w) + coef[1] * np.cos(w))
        terms.append([float(fi), a1, float(snr)])
        excl |= np.abs(freqs - fi) <= 1.5 * df                          # do not re-pick the same bin
    return r, dict(prewhiten_n=len(terms), prewhiten_comb=len(combs), prewhiten_terms=json.dumps(terms),
                   prewhiten_comb_freqs=json.dumps(combs))
