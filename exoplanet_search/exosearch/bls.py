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


# ---------------------------------------------------------------- red-noise-aware detection statistic
def noise_profile(t, f, durations=DURATIONS):
    """Per-star noise at each trial duration D on the detrended, unbinned flux (upward outliers already
    clipped by detrend; gaps handled by requiring >= half the expected cadences in each window).
      sigma_w   = 1.4826 MAD(residuals)                           (white noise per cadence)
      CDPP_D    = 1.4826 MAD(running mean of residuals over D)
      beta_D    = CDPP_D sqrt(N_D) / sigma_w,  N_D = D / cadence   (red-noise factor; 1 for white noise)"""
    t = np.asarray(t, float)
    r = np.asarray(f, float) - np.median(f)
    sw = float(1.4826 * np.median(np.abs(r - np.median(r))))
    cad = float(np.median(np.diff(t)))
    cs = np.r_[0.0, np.cumsum(r)]
    out = {}
    for D in durations:
        j = np.searchsorted(t, t + D)
        n = j - np.arange(t.size)
        ok = n >= 0.5 * D / cad
        if ok.sum() < 10:
            out[float(D)] = (np.nan, np.nan)
            continue
        m = (cs[j[ok]] - cs[np.arange(t.size)[ok]]) / n[ok]
        cd = float(1.4826 * np.median(np.abs(m - np.median(m))))
        out[float(D)] = (cd, float(cd * np.sqrt(D / cad) / sw) if sw > 0 else np.nan)
    return sw, cad, out


def transit_count(t, period, t0, duration, bin_d):
    """Coverage-weighted number of transits with data: sum over epochs of min(1, n_in / n_expected)."""
    ph = ((t - t0 + 0.5 * period) % period) - 0.5 * period
    inn = np.abs(ph) < duration / 2
    if not inn.any():
        return 0.0
    ep = np.round((t[inn] - t0) / period).astype(int)
    _, cnt = np.unique(ep, return_counts=True)
    return float(np.minimum(1.0, cnt / max(duration / bin_d, 1.0)).sum())


def search_s(tb, fb, tflat, fflat, periods=None, durations=DURATIONS):
    """BLS scan thresholded on S = depth / (CDPP_D / sqrt(N_tr)) instead of SDE.

    For each D the binned points get a constant uncertainty sigma_D = CDPP_D sqrt(D / bin), so astropy's
    'snr' objective (depth / depth_err) equals depth sqrt(N_in / N_Db) / CDPP_D = S with a fractional N_tr.
    The star's best peak is the max over (P, D, T0); S is then recomputed there with the coverage-weighted
    transit count. SDE of the max-over-D S periodogram is logged but not thresholded."""
    if periods is None:
        periods = period_grid(tb, durations=durations)
    sw, cad, prof = noise_profile(tflat, fflat, durations)
    bin_d = float(np.median(np.diff(tb)))
    smax = np.full(len(periods), -np.inf)
    best = None
    for D in durations:
        cd, beta = prof[float(D)]
        if not (cd == cd and cd > 0):
            continue
        dy = np.full(tb.size, cd * np.sqrt(max(D / bin_d, 1.0)))
        r = BoxLeastSquares(tb, fb, dy=dy).power(periods, [D], objective="snr", method="fast")
        p = np.nan_to_num(np.asarray(r.power), nan=-np.inf)
        upd = p > smax
        smax[upd] = p[upd]
        i = int(np.argmax(p))
        if best is None or p[i] > best["S_scan"]:
            best = dict(period=float(r.period[i]), t0=float(r.transit_time[i]), duration=float(D),
                        depth_box=float(r.depth[i]), S_scan=float(p[i]), cdpp_d=cd, beta_d=beta)
    if best is None:
        raise ValueError("no valid noise estimate")
    ntr = transit_count(tb, best["period"], best["t0"], best["duration"], bin_d)
    s_val = best["depth_box"] * np.sqrt(ntr) / best["cdpp_d"] if ntr > 0 else 0.0
    smax[~np.isfinite(smax)] = 0.0
    sd, pdet = sde(smax)
    best.update(S=float(s_val), n_tr=ntr, sigma_w=sw, cadence_d=cad, sde=sd, periods=periods, power=smax,
                power_detrended=pdet, backend="cpu-S",
                cdpp_by_d={f"{k:.2f}": v[0] for k, v in prof.items()},
                beta_by_d={f"{k:.2f}": v[1] for k, v in prof.items()})
    return best
