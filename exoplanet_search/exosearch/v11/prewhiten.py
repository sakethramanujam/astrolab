"""Pulsator prewhitening above fmin (default 5 c/d) with the comb guard (docs/v1.1_plan.md item 4,
amended 2026-11-11 for Task T0).

Iterative sinusoid removal on the quality-masked, normalised flux before detrending:
  1. Star-level null calibration: check maximum Lomb-Scargle SNR on grid [fmin, fnyq]. If below
     min_detect_snr (calibrated on flipped null light curves to avoid noise false alarms on quiet stars,
     default 6.0), the star has no significant pulsator; no terms are removed.
  2. Lomb-Scargle amplitude spectrum of the residual on fixed grid fmin + df * arange(n),
     df = 1 / (oversample * baseline).
  3. Highest peak f (comb-protected frequencies excluded), parabolic interpolation, LSQ fit at f.
  4. Comb guard (Task T0 fix): before removing f, test whether f is an orbital harmonic m of a
     planetary transit with period P = m / f in [0.2, 5.0] days. For m in (1, 2, 3), check 2f and 3f
     against the box harmonic ladder. For m >= 2, check the fundamental f0 = f / m and whether
     the amplitude ratio A(f0)/A(f) matches the box ladder of a periodic transit of duration 1-6 h.
     If a transit comb is detected, log it in combs, exclude the entire harmonic series k * (f / m)
     from later picks (+-2/baseline), and do not remove f.
  5. Two-halves sector coherence: verify that candidate pulsation f is coherent across both halves
     of the observation baseline (significant in both halves, consistent amplitudes, and phase-coherent
     with the global fit). Non-coherent peaks (systematics, flares, or noise) are excluded.
  6. Stop when amplitude S/N < snr_stop or after max_terms removals.
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


def transit_comb_matches(t, r, fi, a1, noise_amp, fnyq, durations=(1 / 24, 6 / 24),
                         p_range=(0.2, 5.0), tol=2.5):
    """Detect if fi is an orbital harmonic m of a periodic transit comb (P = m / fi)."""
    ms = []
    # 1. Check multiples 1..3 with 2*fi and 3*fi if accessible
    if 3 * fi < fnyq:
        a2, _ = _fit(t, r, 2 * fi)
        a3, _ = _fit(t, r, 3 * fi)
        m123 = comb_multiples(fi, a1, a2, a3, durations, (1, 2, 3), tol)
        ms.extend(m123)

    # 2. Check candidate fundamentals f0 = fi / m for m >= 2 across planetary periods
    max_m = int(np.floor(p_range[1] * fi))
    min_m = max(2, int(np.ceil(p_range[0] * fi)))
    lt = np.log(max(tol, 2.5))

    for m in range(min_m, min_m + min(max_m - min_m + 1, 35)):
        if m in ms:
            continue
        P = m / fi
        f0 = fi / m
        af0, _ = _fit(t, r, f0)
        if af0 < 2.0 * noise_amp:
            continue
        # Duty cycle range for plausible planetary durations
        D = np.linspace(durations[0], min(durations[1], 0.35 * P), 31)
        q = D / P
        s1 = np.abs(np.sin(np.pi * q))
        sm = np.abs(np.sin(np.pi * m * q)) / m
        ratio_model = s1 / np.maximum(sm, 1e-6)
        obs_ratio = af0 / max(a1, 1e-9)
        if np.any(np.abs(np.log(obs_ratio / np.maximum(ratio_model, 1e-3))) < lt):
            ms.append(m)
    return sorted(list(set(ms)))


def is_coherent_sinusoid(t, r, fi, coef, noise_amp, m1, m2):
    """Check if sinusoid at fi is coherent across both halves of the observation baseline."""
    if np.sum(m1) < 50 or np.sum(m2) < 50:
        return True
    a1_half, c1 = _fit(t[m1], r[m1], fi)
    a2_half, c2 = _fit(t[m2], r[m2], fi)
    noise_half = noise_amp * np.sqrt(2.0)
    # Both halves must show significant amplitude
    if a1_half < 2.5 * noise_half or a2_half < 2.5 * noise_half:
        return False
    # Amplitudes across halves must be consistent
    if a1_half > 3.0 * a2_half or a2_half > 3.0 * a1_half:
        return False
    # Phase must be coherent across the halves relative to global fit
    phi_f = np.arctan2(coef[1], coef[0])
    phi1 = np.arctan2(c1[1], c1[0])
    phi2 = np.arctan2(c2[1], c2[0])
    if np.cos(phi1 - phi_f) < 0.3 or np.cos(phi2 - phi_f) < 0.3:
        return False
    return True


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

    # Star-level null calibration:
    # On null (flipped) light curves without pulsations, the maximum noise peak across
    # ~40,000 bins is ~4.3-4.5. Require initial peak SNR >= 6.0 before activating prewhitening.
    amp_init = _amp_spectrum(t, r, freqs)
    i_max = int(np.argmax(amp_init))
    noise_max = _noise(amp_init, freqs, i_max, win, df)
    snr_max = amp_init[i_max] / noise_max if noise_max > 0 else 0.0
    min_detect_snr = max(float(cfg.get("snr_detect", 6.0)), float(cfg["snr_stop"]))
    if snr_max < min_detect_snr:
        return r, dict(prewhiten_n=0, prewhiten_comb=0, prewhiten_terms="[]", prewhiten_comb_freqs="[]")

    t_mid = 0.5 * (t.min() + t.max())
    m1 = t < t_mid
    m2 = t >= t_mid
    has_halves = (np.sum(m1) >= 50 and np.sum(m2) >= 50)

    durations = tuple(cfg.get("comb_durations_d", [0.0416667, 0.25]))
    tol = float(cfg.get("comb_ratio_tol", 2.0))

    terms, combs, excl = [], [], np.zeros(n, bool)
    for _ in range(int(cfg["max_terms"]) + 4 * int(cfg["max_terms"])):      # comb picks do not count as terms
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
            ms = transit_comb_matches(t, r, fi, a1, noise, fnyq, durations=durations, tol=tol)
            if ms:
                a2, a3 = 0.0, 0.0
                if 2 * fi < fnyq:
                    a2, _ = _fit(t, r, 2 * fi)
                if 3 * fi < fnyq:
                    a3, _ = _fit(t, r, 3 * fi)
                combs.append([float(fi), a1, float(a2), float(a3), ms])
                for m in ms:              # protect the whole harmonic series of every matching fundamental f / m
                    f0 = fi / m
                    for h in range(1, int(fnyq / f0) + 1):
                        excl |= np.abs(freqs - h * f0) <= 2.0 / base       # main lobe and first sidelobes
                continue

        # Two-halves sector coherence: verify signal is stationary across the sector
        if has_halves and not is_coherent_sinusoid(t, r, fi, coef, noise, m1, m2):
            excl |= np.abs(freqs - fi) <= 1.5 * df
            continue

        w = 2 * np.pi * fi * t
        r = r - (coef[0] * np.sin(w) + coef[1] * np.cos(w))
        terms.append([float(fi), a1, float(snr)])
        excl |= np.abs(freqs - fi) <= 1.5 * df                          # do not re-pick the same bin
    return r, dict(prewhiten_n=len(terms), prewhiten_comb=len(combs), prewhiten_terms=json.dumps(terms),
                   prewhiten_comb_freqs=json.dumps(combs))
