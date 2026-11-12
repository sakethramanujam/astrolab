"""Pulsator prewhitening above fmin (default 5 c/d) with the comb guard (docs/v1.1_plan.md item 4,
amended 2026-11-11 for Task T0; see docs/v1.1_prewhitening_amendment.md).

Iterative sinusoid removal on the quality-masked, normalised flux before detrending:
  1. Star-level gate: check maximum Lomb-Scargle SNR on grid [fmin, fnyq]. If below
     snr_detect (configured in configs/v1.1.json["prewhiten"]), the star has no significant pulsator;
     no terms are removed.
  2. Lomb-Scargle amplitude spectrum of the residual on fixed grid fmin + df * arange(n),
     df = 1 / (oversample * baseline).
  3. Highest peak f (comb-protected frequencies excluded), parabolic interpolation, LSQ fit at f.
  4. Comb guard (Task T0 fix): before removing f, test whether f is an orbital harmonic m of a
     planetary transit with period P = m / f in comb_p_range_d (default [0.2, 5.0] days).
     For m in comb_multiples (1, 2, 3), check 2f and 3f against the box harmonic ladder (requiring
     SNR >= comb_min_snr, amplitude ratio tolerance comb_ratio_tol, and common-epoch phase consistency).
     For m >= 2:
       - Require local noise-calibrated fundamental f0 = f / m with SNR >= comb_min_snr.
       - Test box harmonic ratio A(f0)/A(f) against the box model R(m, q) across durations
         comb_durations_d, bounding m*q < 0.8 to prevent model pole divergence.
       - Require phase consistency between f and f0 pointing to a common epoch within comb_phase_tol.
       - Require at least one additional harmonic (e.g. k=2 or k=3) with SNR >= comb_min_snr and
         consistent epoch.
     If a transit comb is detected, log it in combs, exclude the harmonic series k * (f / m)
     from later picks (+-comb_exclusion_width/baseline) up to h_max = comb_max_harmonic_factor / (dmin * f0), capping
     total exclusion coverage to comb_max_coverage (default 20%), and do not remove f.
  5. Two-halves sector coherence (glitch and transient rejector): verify that candidate pulsation f is
     coherent across both halves of the observation baseline (significant in both halves >= half_min_snr,
     amplitude ratio <= half_amp_ratio, and phase-coherent cos(phi - phi_full) >= half_cos_phase).
     Non-coherent peaks are excluded without flux removal.
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
    p = LombScargle(t, r, fit_mean=True, center_data=True, normalization="psd").power(
        freqs, method="fast", assume_regular_frequency=True
    )
    return np.sqrt(np.maximum(4.0 * p / t.size, 0.0))


def _noise(amp, freqs, i, window, df):
    m = (np.abs(freqs - freqs[i]) <= window) & (np.abs(freqs - freqs[i]) > 3 * df)
    return float(np.mean(amp[m])) if m.any() else np.inf


def _local_noise(f, freqs, amp, window, df):
    m = (np.abs(freqs - f) <= window) & (np.abs(freqs - f) > 3 * df)
    return float(np.median(amp[m])) if np.any(m) else float(np.median(amp))


def box_ladder(m, q):
    s1 = np.abs(np.sin(np.pi * m * q))
    r2 = np.abs(np.cos(np.pi * m * q))
    r3 = np.abs(np.sin(3 * np.pi * m * q)) / (3 * np.maximum(s1, 1e-12))
    return r2, r3


def comb_multiples(f, a1, a2, a3, durations=(0.0208, 0.25), multiples=(1, 2, 3), tol=2.0, n=61):
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


def is_comb(f, a1, a2, a3, durations=(0.0208, 0.25), multiples=(1, 2, 3), tol=2.0, n=61):
    return bool(comb_multiples(f, a1, a2, a3, durations, multiples, tol, n))


def transit_comb_matches(t, r, fi, a1, noise_amp, fnyq, base=27.0, lf_freqs=None, lf_amp=None, cfg=None):
    """Detect if fi is an orbital harmonic m of a periodic transit comb (P = m / fi)."""
    cfg = cfg or {}
    durations = tuple(cfg.get("comb_durations_d", [0.0208, 0.25]))
    p_range = tuple(cfg.get("comb_p_range_d", [0.2, 5.0]))
    tol = float(cfg.get("comb_ratio_tol", 2.0))
    min_snr = float(cfg.get("comb_min_snr", 3.0))
    phase_tol = float(cfg.get("comb_phase_tol", 0.2))

    ms = []
    # Phase of fi
    w_i = 2 * np.pi * fi * t
    X_i = np.stack([np.sin(w_i), np.cos(w_i), np.ones_like(t)], axis=1)
    cf_i = np.linalg.solve(X_i.T @ X_i, X_i.T @ r)
    phi_i = float(np.arctan2(cf_i[0], cf_i[1]))
    t0_i = ((phi_i - np.pi) / (2 * np.pi * fi)) % (1.0 / fi)

    # 1. Candidate fundamentals f0 = fi / m for m >= 2 across planetary periods
    max_m = int(np.floor(p_range[1] * fi))
    min_m = max(2, int(np.ceil(p_range[0] * fi)))
    lt = np.log(tol)
    df_base = 1.0 / (5.0 * base)

    for m in range(min_m, max_m + 1):
        P = m / fi
        f0 = fi / m

        nz0 = _local_noise(f0, lf_freqs, lf_amp, 1.0, df_base) if lf_freqs is not None else noise_amp
        af0, cf0 = _fit(t, r, f0)
        if af0 < min_snr * nz0:
            continue

        D_grid = np.linspace(durations[0], min(durations[1], 0.35 * P), 31)
        valid = (m * D_grid / P < 0.8) & (D_grid / P < 0.5)
        if not np.any(valid):
            continue
        qv = D_grid[valid] / P
        s1 = np.abs(np.sin(np.pi * qv))
        sm = np.abs(np.sin(np.pi * m * qv)) / m
        ratio_model = s1 / np.maximum(sm, 1e-6)
        obs_ratio = af0 / max(a1, 1e-9)
        if not np.any(np.abs(np.log(obs_ratio / np.maximum(ratio_model, 1e-3))) < lt):
            continue

        phi_0 = float(np.arctan2(cf0[0], cf0[1]))
        t0_0 = ((phi_0 - np.pi) / (2 * np.pi * f0)) % P
        dt_i = abs((t0_0 - t0_i + 0.5 / fi) % (1.0 / fi) - 0.5 / fi)
        if dt_i > phase_tol / fi:
            continue

        # Agreement with another harmonic k * f0 (k=3 if m=2 else 2)
        k = 3 if m == 2 else 2
        fk = k * f0
        if fk < fnyq:
            nzk = _local_noise(fk, lf_freqs, lf_amp, 1.0, df_base) if (lf_freqs is not None and fk <= lf_freqs[-1]) else noise_amp
            afk, cfk = _fit(t, r, fk)
            if afk < min_snr * nzk:
                continue
            rk_model = np.abs(np.sin(np.pi * k * qv)) / (k * np.maximum(s1, 1e-6))
            obs_k = afk / max(af0, 1e-9)
            if not np.any((np.abs(np.log(obs_k / np.maximum(rk_model, 1e-3))) < lt) &
                          (np.abs(np.log(obs_ratio / np.maximum(ratio_model, 1e-3))) < lt)):
                continue
            phi_k = float(np.arctan2(cfk[0], cfk[1]))
            t0_k = ((phi_k - np.pi) / (2 * np.pi * fk)) % (1.0 / fk)
            dt_k = abs((t0_0 - t0_k + 0.5 / fk) % (1.0 / fk) - 0.5 / fk)
            if dt_k > phase_tol / fk:
                continue

        ms.append(m)

    # 2. Multiples 1..3 if no fundamental match found and fi period matches comb_p_range_d
    if not ms and 3 * fi < fnyq:
        p_fund = 1.0 / fi
        if p_range[0] <= p_fund <= p_range[1]:
            a2, c2 = _fit(t, r, 2 * fi)
            a3, c3 = _fit(t, r, 3 * fi)
            if a2 >= min_snr * noise_amp and a3 >= min_snr * noise_amp:
                m123 = comb_multiples(fi, a1, a2, a3, durations, tuple(cfg.get("comb_multiples", [1, 2, 3])), tol)
                for m in m123:
                    phi_2 = float(np.arctan2(c2[0], c2[1]))
                    t0_2 = ((phi_2 - np.pi) / (2 * np.pi * (2 * fi))) % (0.5 / fi)
                    dt_2 = abs((t0_i - t0_2 + 0.25 / fi) % (0.5 / fi) - 0.25 / fi)
                    if dt_2 <= phase_tol / (2 * fi):
                        ms.append(m)

    return sorted(list(set(ms)))


def is_coherent_sinusoid(t, r, fi, coef, noise_amp, m1, m2, cfg):
    """Check if sinusoid at fi is coherent across both halves of the observation baseline."""
    if np.sum(m1) < 50 or np.sum(m2) < 50:
        return True
    a1_half, c1 = _fit(t[m1], r[m1], fi)
    a2_half, c2 = _fit(t[m2], r[m2], fi)
    noise_half = noise_amp * np.sqrt(2.0)
    # Both halves must show significant amplitude
    half_min_snr = float(cfg.get("half_min_snr", 2.0))
    if a1_half < half_min_snr * noise_half or a2_half < half_min_snr * noise_half:
        return False
    # Amplitudes across halves must be consistent
    half_amp_ratio = float(cfg.get("half_amp_ratio", 5.0))
    if a1_half > half_amp_ratio * a2_half or a2_half > half_amp_ratio * a1_half:
        return False
    # Phase must be coherent across the halves relative to global fit
    half_cos_phase = float(cfg.get("half_cos_phase", -0.5))
    phi_f = np.arctan2(coef[1], coef[0])
    phi1 = np.arctan2(c1[1], c1[0])
    phi2 = np.arctan2(c2[1], c2[0])
    if np.cos(phi1 - phi_f) < half_cos_phase or np.cos(phi2 - phi_f) < half_cos_phase:
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

    # Star-level gate:
    amp_init = _amp_spectrum(t, r, freqs)
    i_max = int(np.argmax(amp_init))
    noise_max = _noise(amp_init, freqs, i_max, win, df)
    snr_max = amp_init[i_max] / noise_max if noise_max > 0 else 0.0
    min_detect_snr = max(float(cfg.get("snr_detect", 5.0)), float(cfg["snr_stop"]))
    if snr_max < min_detect_snr:
        return r, dict(prewhiten_n=0, prewhiten_comb=0, prewhiten_terms="[]", prewhiten_comb_freqs="[]")

    # Low-frequency spectrum for local noise estimation at f0 < fmin, extending high enough
    # for local-noise windows (window ~ 1.0 c/d) around f0 and its harmonics k*f0 up to ~10 c/d
    lf_freqs = np.arange(0.1, max(12.0, fmin + 2.0), df)
    lf_amp = _amp_spectrum(t, r, lf_freqs)

    t_mid = 0.5 * (t.min() + t.max())
    m1 = t < t_mid
    m2 = t >= t_mid
    has_halves = (np.sum(m1) >= 50 and np.sum(m2) >= 50)

    max_harm_factor = float(cfg.get("comb_max_harmonic_factor", 2.5))
    excl_width = float(cfg.get("comb_exclusion_width", 2.0))
    max_coverage = float(cfg.get("comb_max_coverage", 0.20))
    dmin = float(cfg.get("comb_durations_d", [0.0208, 0.25])[0])

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
            ms = transit_comb_matches(t, r, fi, a1, noise, fnyq, base, lf_freqs, lf_amp, cfg)
            if ms:
                a2, a3 = 0.0, 0.0
                if 2 * fi < fnyq:
                    a2, _ = _fit(t, r, 2 * fi)
                if 3 * fi < fnyq:
                    a3, _ = _fit(t, r, 3 * fi)
                combs.append([float(fi), a1, float(a2), float(a3), ms])
                for m in ms:
                    f0 = fi / m
                    h_max = min(int(fnyq / f0), int(np.ceil(max_harm_factor / (dmin * f0))))
                    for h in range(1, h_max + 1):
                        new_excl = excl | (np.abs(freqs - h * f0) <= excl_width / base)
                        if new_excl.mean() <= max_coverage:
                            excl = new_excl
                continue

        # Two-halves sector coherence: verify signal is stationary across the sector
        if has_halves and not is_coherent_sinusoid(t, r, fi, coef, noise, m1, m2, cfg):
            excl |= np.abs(freqs - fi) <= 1.5 * df
            continue

        w = 2 * np.pi * fi * t
        r = r - (coef[0] * np.sin(w) + coef[1] * np.cos(w))
        terms.append([float(fi), a1, float(snr)])
        excl |= np.abs(freqs - fi) <= 1.5 * df                          # do not re-pick the same bin

    return r, dict(prewhiten_n=len(terms), prewhiten_comb=len(combs), prewhiten_terms=json.dumps(terms),
                   prewhiten_comb_freqs=json.dumps(combs))
