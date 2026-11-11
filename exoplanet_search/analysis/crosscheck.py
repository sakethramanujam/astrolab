"""Cross-check of completeness, the null tail, ROC and DR25-style reliability.

Reads injection and null products (parquet under a products root, or the sample CSVs) and writes
strict JSON plus PNG figures. Holdout TICs are left out of every fit unless ``holdout_only`` is set.
The interim threshold config is an input, not something this module rewrites.

  python -m analysis --injections analysis/sample_injections.csv --null analysis/sample_null.csv --out /tmp/out
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import subprocess
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import optimize, stats

from . import figures

PACKAGE_ROOT = Path(__file__).resolve().parents[1]          # exoplanet_search/
REPO_ROOT = PACKAGE_ROOT.parent
P_EDGES = np.geomspace(0.5, 13.5, 9)
R_EDGES = np.geomspace(0.5, 20.0, 9)
# Rejecting vetting flags. ruwe_high and gaia_nss are warnings in exosearch.version, not rejects.
# low_s is the in-pipeline screening floor; detection is re-scored at S* instead.
WARNING_FLAGS = {"ruwe_high", "gaia_nss", "low_s"}
DEEP_REJECT = {"v_shape", "secondary_scan", "rp_gt_2rjup", "density_mismatch",
               "centroid_offset", "dip_off_target", "known_eb"}
PC_LABELS = {"pc", "planet", "candidate", "kp", "cp", "confirmed", "planet_candidate"}
FP_LABELS = {"fp", "false_positive", "false positive", "fa", "false_alarm", "eb"}
WILSON_Z = 1.959963984540054  # 95%


# ---------------------------------------------------------------- JSON
def json_safe(o):
    """Recursively replace NaN/Inf and numpy scalars so the result is strict JSON."""
    if isinstance(o, dict):
        return {str(k): json_safe(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [json_safe(v) for v in o]
    if isinstance(o, np.ndarray):
        return json_safe(o.tolist())
    if isinstance(o, np.generic):
        o = o.item()
    if isinstance(o, float) and not math.isfinite(o):
        return None
    if o is pd.NA:
        return None
    return o


def dump_json(path, obj):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    text = json.dumps(json_safe(obj), indent=2, allow_nan=False) + "\n"
    path.write_text(text)
    return text


# ---------------------------------------------------------------- small stats
def wilson_interval(k, n, z=WILSON_Z):
    """Wilson score interval. Returns (None, None) when n is 0."""
    k = float(k)
    n = float(n)
    if n <= 0:
        return (None, None)
    p = k / n
    z2 = z * z
    den = 1.0 + z2 / n
    center = (p + z2 / (2.0 * n)) / den
    half = z * math.sqrt(p * (1.0 - p) / n + z2 / (4.0 * n * n)) / den
    return (max(0.0, center - half), min(1.0, center + half))


def expected_snr(depth, n_tr, cdpp):
    """SCHEMA expected S: depth_inj / cdpp_dur * sqrt(n_tr). None when the ratio is undefined."""
    try:
        depth, n_tr, cdpp = float(depth), float(n_tr), float(cdpp)
    except (TypeError, ValueError):
        return None
    if not (math.isfinite(depth) and math.isfinite(n_tr) and math.isfinite(cdpp)):
        return None
    if n_tr <= 0 or cdpp <= 0 or depth < 0:
        return None
    return depth * math.sqrt(n_tr) / cdpp


def _ci(values, lo=0.025, hi=0.975, min_n=20):
    a = np.asarray(list(values), float)
    a = a[np.isfinite(a)]
    if a.size < min_n:
        return None
    return [float(np.quantile(a, lo)), float(np.quantile(a, hi))]


def _star_groups(tics: np.ndarray):
    order = np.argsort(tics, kind="mergesort")
    sorted_tics = tics[order]
    if sorted_tics.size == 0:
        return []
    cuts = np.flatnonzero(sorted_tics[1:] != sorted_tics[:-1]) + 1
    return np.split(order, cuts)


def _draw_groups(groups, rng):
    choice = rng.integers(0, len(groups), size=len(groups))
    return np.concatenate([groups[int(i)] for i in choice])


def mann_whitney_auc(pos_scores, neg_scores):
    """AUC of the score, higher is more positive. Ties use midranks. None if either class is empty."""
    pos = np.asarray(pos_scores, float)
    neg = np.asarray(neg_scores, float)
    pos = pos[np.isfinite(pos)]
    neg = neg[np.isfinite(neg)]
    n_pos, n_neg = pos.size, neg.size
    if n_pos == 0 or n_neg == 0:
        return None
    scores = np.concatenate([pos, neg])
    ranks = stats.rankdata(scores, method="average")
    return float((ranks[:n_pos].sum() - n_pos * (n_pos + 1) / 2.0) / (n_pos * n_neg))


def roc_curve(pos_scores, neg_scores, max_points=200):
    pos = np.asarray(pos_scores, float)
    neg = np.asarray(neg_scores, float)
    pos = pos[np.isfinite(pos)]
    neg = neg[np.isfinite(neg)]
    auc = mann_whitney_auc(pos, neg)
    if pos.size == 0 or neg.size == 0:
        return {"fpr": [], "tpr": [], "auc": None, "n_pos": int(pos.size), "n_neg": int(neg.size)}
    scores = np.concatenate([pos, neg])
    labels = np.concatenate([np.ones(pos.size), np.zeros(neg.size)])
    order = np.argsort(-scores, kind="mergesort")
    labels = labels[order]
    tps = np.cumsum(labels)
    fps = np.cumsum(1.0 - labels)
    fpr = np.concatenate([[0.0], fps / neg.size])
    tpr = np.concatenate([[0.0], tps / pos.size])
    if fpr.size > max_points:
        idx = np.linspace(0, fpr.size - 1, max_points).astype(int)
        fpr, tpr = fpr[idx], tpr[idx]
    return {"fpr": fpr.tolist(), "tpr": tpr.tolist(), "auc": auc,
            "n_pos": int(pos.size), "n_neg": int(neg.size)}


# ---------------------------------------------------------------- gamma CDF
def _gamma_nll(theta, snr, detected):
    shape = math.exp(theta[0])
    scale = math.exp(theta[1])
    p = stats.gamma.cdf(snr, shape, scale=scale)
    p = np.clip(p, 1e-12, 1.0 - 1e-12)
    return float(-np.sum(detected * np.log(p) + (1.0 - detected) * np.log(1.0 - p)))


def _fit_gamma_once(snr, detected, start=None):
    snr = np.asarray(snr, float)
    detected = np.asarray(detected, float)
    starts = [np.log([2.0, 5.0]), np.log([1.2, 10.0]), np.log([6.0, 2.0])]
    if start is not None:
        starts = [np.asarray(start, float)]
    best = None
    bounds = [(math.log(0.15), math.log(40.0)), (math.log(0.15), math.log(300.0))]
    for st in starts:
        res = optimize.minimize(_gamma_nll, st, args=(snr, detected), method="L-BFGS-B", bounds=bounds)
        if best is None or res.fun < best.fun:
            best = res
    shape, scale = (float(math.exp(best.x[0])), float(math.exp(best.x[1])))
    on_bound = (
        abs(best.x[0] - bounds[0][0]) < 1e-4 or abs(best.x[0] - bounds[0][1]) < 1e-4
        or abs(best.x[1] - bounds[1][0]) < 1e-4 or abs(best.x[1] - bounds[1][1]) < 1e-4
    )
    return shape, scale, bool(on_bound), best.x


def fit_gamma_cdf(snr, detected, tics, *, n_boot, seed):
    """Unbinned Bernoulli MLE of a two-parameter gamma CDF, with a star-level bootstrap CI."""
    snr = np.asarray(snr, float)
    detected = np.asarray(detected, float)
    tics = np.asarray(tics)
    ok = np.isfinite(snr) & np.isfinite(detected) & (snr >= 0)
    snr, detected, tics = snr[ok], detected[ok], tics[ok]
    n_det = int(detected.sum())
    n_miss = int(len(detected) - n_det)
    out = {
        "n": int(len(snr)), "n_detected": n_det, "n_missed": n_miss,
        "shape": None, "scale": None, "shape_ci95": None, "scale_ci95": None,
        "on_boundary": None, "n_boot": int(n_boot), "n_boot_ok": 0, "seed": int(seed),
        "curve": None, "meaningful": False,
        "note": None,
    }
    if len(snr) < 8 or n_det < 2 or n_miss < 2:
        out["note"] = "too few injections with both detections and non-detections to fit a gamma CDF"
        return out
    shape, scale, on_bound, theta = _fit_gamma_once(snr, detected)
    out.update(shape=shape, scale=scale, on_boundary=on_bound)
    rng = np.random.default_rng(seed)
    groups = _star_groups(tics)
    shapes, scales = [], []
    for _ in range(int(n_boot)):
        idx = _draw_groups(groups, rng)
        y = detected[idx]
        if y.sum() < 2 or (len(y) - y.sum()) < 2:
            continue
        try:
            a, sc, bound, _ = _fit_gamma_once(snr[idx], y, start=theta)
        except Exception:
            continue
        if bound or not (0.2 < a < 35 and 0.2 < sc < 250):
            continue
        shapes.append(a)
        scales.append(sc)
    out["n_boot_ok"] = int(len(shapes))
    out["shape_ci95"] = _ci(shapes)
    out["scale_ci95"] = _ci(scales)
    grid = np.linspace(0.0, max(40.0, float(np.quantile(snr, 0.99)) if len(snr) else 40.0), 61)
    p = stats.gamma.cdf(grid, shape, scale=scale)
    if shapes:
        band = np.vstack([stats.gamma.cdf(grid, a, scale=sc) for a, sc in zip(shapes, scales)])
        p_lo = np.quantile(band, 0.025, axis=0)
        p_hi = np.quantile(band, 0.975, axis=0)
    else:
        p_lo = p_hi = None
    out["curve"] = {
        "snr": grid.tolist(),
        "p": np.asarray(p, float).tolist(),
        "p_lo": None if p_lo is None else np.asarray(p_lo, float).tolist(),
        "p_hi": None if p_hi is None else np.asarray(p_hi, float).tolist(),
    }
    out["meaningful"] = bool(len(snr) >= 200 and min(n_det, n_miss) >= 30 and not on_bound and out["shape_ci95"])
    if not out["meaningful"]:
        out["note"] = "gamma CDF was fit, but the sample is too small (or the MLE sits on a bound) for the parameters to be a measurement"
    return out


# ---------------------------------------------------------------- GPD / Gumbel
def _pwm_gpd(y):
    y = np.sort(np.asarray(y, float))
    n = y.size
    if n < 2:
        return None
    a0 = float(y.mean())
    j = np.arange(1, n + 1)
    a1 = float(np.mean((n - j) / (n - 1.0) * y))
    den = a0 - 2.0 * a1
    if den == 0 or not math.isfinite(den):
        return None
    xi = 2.0 - a0 / den
    sigma = 2.0 * a0 * a1 / den
    if not (math.isfinite(xi) and math.isfinite(sigma) and sigma > 0 and -0.9 < xi < 1.5):
        return None
    return xi, sigma


def _fit_gpd_excesses(y, min_k=15):
    y = np.asarray(y, float)
    y = y[np.isfinite(y) & (y > 0)]
    if y.size < min_k:
        return None
    xi = sigma = None
    method = None
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        try:
            c, _loc, scale = stats.genpareto.fit(y, floc=0)
            xi, sigma = float(c), float(scale)
            method = "mle"
        except Exception:
            xi = sigma = None
    if xi is None or not (math.isfinite(xi) and math.isfinite(sigma) and sigma > 0 and -0.9 < xi < 1.5):
        pwm = _pwm_gpd(y)
        if pwm is None:
            return None
        xi, sigma = pwm
        method = "pwm"
    return {"xi": float(xi), "sigma": float(sigma), "k": int(y.size), "method": method}


def gpd_s_star(xi, sigma, u, zeta, p):
    """S* with P(S > S*) = p, using a GPD for the excess over u and zeta = P(S > u)."""
    if sigma <= 0 or zeta <= 0 or p <= 0 or p >= zeta:
        return None
    q = p / zeta
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        excess = float(stats.genpareto.isf(q, xi, loc=0.0, scale=sigma))
    if not math.isfinite(excess) or excess < 0:
        if xi < 0:
            excess = -sigma / xi
        else:
            return None
    return float(u + excess)


def _gpd_from_sample(S, n_stars, percentile):
    S = np.asarray(S, float)
    S = S[np.isfinite(S)]
    n = int(S.size)
    p = 1.0 / float(n_stars) if n_stars else None
    if n < 20 or p is None or p <= 0:
        return None
    u = float(np.quantile(S, percentile / 100.0))
    excess = S[S > u] - u
    fit = _fit_gpd_excesses(excess)
    zeta = float(excess.size / n)
    max_obs = float(S.max())
    rate_in_tail = bool(zeta > 0 and p < zeta)
    xi = sigma = method = None
    if fit is not None:
        xi, sigma, method = fit["xi"], fit["sigma"], fit["method"]
    empirical = float(np.quantile(S, max(0.0, 1.0 - p))) if 0 < p < 1 else None
    if xi is not None and rate_in_tail:
        s_star = gpd_s_star(xi, sigma, u, zeta, p)
        how = "gpd"
    elif not rate_in_tail:
        s_star = empirical
        how = "empirical_quantile_because_rate_is_not_in_the_tail"
    else:
        # The target rate sits in the tail, but there are not enough exceedances to fit it.
        s_star = None
        how = "gpd_not_fit"
    extrapolated = None if s_star is None else bool(s_star > max_obs)
    return {
        "u": u, "k_tail": int(excess.size), "n": n, "zeta": zeta, "max_observed": max_obs,
        "n_stars": int(n_stars), "false_alarm_rate": p, "xi": xi, "sigma": sigma, "method": method,
        "s_star": s_star, "s_star_how": how, "empirical_quantile": empirical,
        "rate_in_tail": rate_in_tail, "extrapolated": extrapolated,
        "percentile": float(percentile),
    }


def mean_residual_life(S, n_grid=24):
    S = np.asarray(S, float)
    S = np.sort(S[np.isfinite(S)])
    if S.size < 10:
        return {"u": [], "mean_excess": [], "se": [], "k": []}
    grid = np.linspace(float(np.quantile(S, 0.50)), float(np.quantile(S, 0.95)), n_grid)
    u, mean, se, k = [], [], [], []
    for thr in grid:
        exc = S[S > thr] - thr
        if exc.size < 5:
            continue
        u.append(float(thr))
        mean.append(float(exc.mean()))
        se.append(float(exc.std(ddof=1) / math.sqrt(exc.size)) if exc.size > 1 else None)
        k.append(int(exc.size))
    return {"u": u, "mean_excess": mean, "se": se, "k": k}


def threshold_stability(S, n_stars, percentile_grid=(80, 85, 90, 92, 95, 97, 98)):
    S = np.asarray(S, float)
    S = S[np.isfinite(S)]
    rows = []
    for pct in percentile_grid:
        if S.size < 20:
            break
        u = float(np.quantile(S, pct / 100.0))
        excess = S[S > u] - u
        fit = _fit_gpd_excesses(excess, min_k=15)
        zeta = float(excess.size / S.size) if S.size else 0.0
        p = 1.0 / float(n_stars)
        s_star = None if fit is None else gpd_s_star(fit["xi"], fit["sigma"], u, zeta, p)
        rows.append({
            "percentile": float(pct), "u": u, "k": int(excess.size),
            "xi": None if fit is None else fit["xi"],
            "sigma": None if fit is None else fit["sigma"],
            "s_star": s_star,
        })
    usable = [r["s_star"] for r in rows if r["s_star"] is not None and r["k"] >= 30 and r["xi"] is not None]
    assessed = len(usable) >= 3
    unstable = False
    if assessed:
        med = float(np.median(usable))
        unstable = bool(med > 0 and (max(usable) - min(usable)) / med > 0.5)
    return {
        "thresholds": rows,
        "assessed": assessed,
        "s_star_unstable": unstable if assessed else None,
        "note": None if assessed else "stability not assessed: fewer than three thresholds with at least 30 exceedances",
    }


def fit_tail(S, tics, *, n_stars, n_boot, seed, percentile=95.0):
    """GPD above the given percentile of per-star max S, with a star-level bootstrap CI on S*."""
    S = np.asarray(S, float)
    tics = np.asarray(tics)
    ok = np.isfinite(S)
    S, tics = S[ok], tics[ok]
    point = _gpd_from_sample(S, n_stars, percentile)
    out = {
        "seed": int(seed), "n_boot": int(n_boot), "n_boot_ok": 0, "percentile": float(percentile),
        "xi_ci95": None, "sigma_ci95": None, "s_star_ci95": None, "meaningful": False,
    }
    if point is None:
        out.update(n=int(S.size), n_stars=int(n_stars), xi=None, sigma=None, s_star=None,
                   rate_in_tail=None, extrapolated=None, max_observed=None if S.size == 0 else float(np.nanmax(S)),
                   note="too few null stars to fit a tail")
        out["stability"] = {"thresholds": [], "assessed": False, "s_star_unstable": None,
                            "note": "stability not assessed"}
        out["mrl"] = mean_residual_life(S)
        return out
    out.update(point)
    rng = np.random.default_rng(seed)
    groups = _star_groups(tics)
    xis, sigmas, sstars = [], [], []
    for _ in range(int(n_boot)):
        draw = _gpd_from_sample(S[_draw_groups(groups, rng)], n_stars, percentile)
        if draw is None or draw["xi"] is None or draw["s_star"] is None:
            continue
        if not math.isfinite(draw["s_star"]):
            continue
        xis.append(draw["xi"])
        sigmas.append(draw["sigma"])
        sstars.append(draw["s_star"])
    out["n_boot_ok"] = int(len(sstars))
    out["xi_ci95"] = _ci(xis)
    out["sigma_ci95"] = _ci(sigmas)
    out["s_star_ci95"] = _ci(sstars)
    out["stability"] = threshold_stability(S, n_stars)
    out["mrl"] = mean_residual_life(S)
    k = out.get("k_tail") or 0
    out["meaningful"] = bool(k >= 50 and S.size >= 1000 and out["s_star_ci95"] and out["rate_in_tail"])
    if out["xi"] is None:
        noun = "exceedance" if k == 1 else "exceedances"
        out["note"] = (f"{k} {noun} above the {percentile:g}th percentile; "
                       "a generalized Pareto fit needs at least 15 and is not a measurement below 50.")
    elif out.get("extrapolated"):
        out["note"] = "S* lies beyond the largest observed null value; the false-alarm rate is an extrapolation"
    elif not out["meaningful"]:
        out["note"] = "GPD was fit, but the tail sample is too small for S* to be a measurement"
    else:
        out["note"] = None
    return out


def fit_gumbel_maxima(S, *, n_stars, n_boot, seed, tics=None):
    """Gumbel fit to per-star maxima. N_eff = exp(loc/scale) for an exponential single trial."""
    S = np.asarray(S, float)
    ok = np.isfinite(S)
    S = S[ok]
    if tics is None:
        tics = np.arange(S.size)
    else:
        tics = np.asarray(tics)[ok]
    assumption = ("N_eff = exp(location / scale). This is the number of i.i.d. exponential trials, "
                  "with the same scale as the fitted Gumbel, whose maximum has that location. "
                  "It is not a count of BLS frequency samples.")
    out = {
        "n": int(S.size), "loc": None, "scale": None, "n_eff": None, "s_star": None,
        "loc_ci95": None, "scale_ci95": None, "n_eff_ci95": None, "s_star_ci95": None,
        "n_boot": int(n_boot), "n_boot_ok": 0, "seed": int(seed), "assumption": assumption,
        "note": None,
    }
    if S.size < 20 or n_stars <= 1:
        out["note"] = "too few per-star maxima to fit a Gumbel"
        return out

    def _one(sample):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            loc, scale = stats.gumbel_r.fit(sample)
        loc, scale = float(loc), float(scale)
        if not (math.isfinite(loc) and math.isfinite(scale) and scale > 0):
            return None
        p = 1.0 / float(n_stars)
        s_star = float(stats.gumbel_r.ppf(1.0 - p, loc, scale))
        return loc, scale, float(math.exp(loc / scale)), s_star

    point = _one(S)
    if point is None:
        out["note"] = "Gumbel fit failed"
        return out
    out["loc"], out["scale"], out["n_eff"], out["s_star"] = point
    out["extrapolated"] = bool(out["s_star"] > float(np.max(S)))
    notes = []
    if S.size < 200:
        notes.append("Too few per-star maxima for N_eff or this S* to be a measurement.")
    if out["extrapolated"]:
        notes.append("Gumbel S* lies beyond the largest observed null maximum.")
    out["note"] = " ".join(notes) or None
    rng = np.random.default_rng(seed)
    groups = _star_groups(tics)
    locs, scales, neffs, sstars = [], [], [], []
    for _ in range(int(n_boot)):
        got = _one(S[_draw_groups(groups, rng)])
        if got is None:
            continue
        locs.append(got[0])
        scales.append(got[1])
        neffs.append(got[2])
        sstars.append(got[3])
    out["n_boot_ok"] = int(len(neffs))
    out["loc_ci95"] = _ci(locs)
    out["scale_ci95"] = _ci(scales)
    out["n_eff_ci95"] = _ci(neffs)
    out["s_star_ci95"] = _ci(sstars)
    return out


# ---------------------------------------------------------------- reliability formula
def reliability_from_counts(n_fp, n_pc, n_rejected, n_sim):
    """Thompson et al. 2018 reliability. E = n_rejected / n_sim from simulations only.

    R = 1 - (N_FP,obs / N_PC,obs) * (1 - E) / E. The Wilson interval on E is transformed
    through that expression. Observed counts are treated as fixed. R is null when the
    observed counts are absent, when N_PC is 0, or when E is 0.
    """
    n_sim = int(n_sim)
    n_rejected = int(n_rejected)
    out = {
        "n_sim": n_sim, "n_rejected": n_rejected, "n_fp_obs": n_fp, "n_pc_obs": n_pc,
        "E": None, "E_ci95": [None, None], "R": None, "R_ci95": [None, None],
        "R_floored_at_0": None, "note": None,
    }
    if n_sim <= 0:
        out["note"] = "no simulated false alarms with S >= S* and vetting logged in this bin; E is undefined"
        return out
    e = n_rejected / n_sim
    lo, hi = wilson_interval(n_rejected, n_sim)
    out["E"] = float(e)
    out["E_ci95"] = [lo, hi]
    if n_fp is None or n_pc is None:
        out["note"] = "E is measured from the simulated nulls; R is left null until observed N_FP and N_PC are supplied"
        if n_sim < 30:
            noun = "false alarm" if n_sim == 1 else "false alarms"
            out["note"] += f" Only {n_sim} simulated {noun} entered this count; E is not a measurement."
        return out
    n_fp, n_pc = int(n_fp), int(n_pc)
    out["n_fp_obs"], out["n_pc_obs"] = n_fp, n_pc
    if n_pc <= 0 or e <= 0:
        out["note"] = "R is undefined because N_PC,obs is 0 or E is 0"
        return out
    ratio = n_fp / n_pc

    def r_of(eff):
        if eff is None or eff <= 0:
            return None
        return float(1.0 - ratio * (1.0 - eff) / eff)

    r = r_of(e)
    out["R"] = r
    out["R_ci95"] = [r_of(lo), r_of(hi)]
    out["R_floored_at_0"] = None if r is None else max(0.0, r)
    notes = []
    if r is not None and r < 0:
        notes.append("the formula is negative: implied false-alarm leakage exceeds N_PC,obs; R_floored_at_0 is 0")
    if n_sim < 30:
        noun = "false alarm" if n_sim == 1 else "false alarms"
        notes.append(f"only {n_sim} simulated {noun} entered this count; E is not a measurement")
    if notes:
        out["note"] = " ".join(notes)
    return out


# ---------------------------------------------------------------- tables
def read_table(path: Path) -> pd.DataFrame:
    path = Path(path)
    if path.suffix == ".parquet":
        return pd.read_parquet(path)
    if path.suffix in {".csv", ".txt"}:
        return pd.read_csv(path)
    try:
        return pd.read_parquet(path)
    except Exception:
        return pd.read_csv(path)


def discover_products(root: Path):
    """Parquet products under ``root`` or ``root/products/v*/`` following SCHEMA.md."""
    root = Path(root)
    bases = [root]
    if (root / "products").is_dir():
        bases.extend(sorted((root / "products").glob("v*")))
    inj, null = [], []
    for base in bases:
        inj.extend(sorted(base.glob("injections/**/*.parquet")))
        null.extend(sorted(base.glob("nulltests/**/*.parquet")))
    finals = [p for p in null if p.name.startswith("null_") and "maxS" in p.name]
    if finals:
        null = finals
    return inj, null


def _infer_mode(path: Path) -> str:
    text = str(path).lower()
    if "scramble" in text:
        return "scramble"
    if "invert" in text:
        return "invert"
    return "unknown"


def load_holdout_tics(path) -> set[int]:
    if path is None or not Path(path).exists():
        return set()
    tics = set()
    for line in Path(path).read_text().splitlines():
        line = line.split("#", 1)[0].strip()
        if line:
            tics.add(int(float(line)))
    return tics


def _as_bool(series: pd.Series) -> pd.Series:
    if series.dtype == bool:
        return series.fillna(False)
    return series.astype(str).str.strip().str.lower().isin(["1", "true", "t", "yes"])


def holdout_mask(df: pd.DataFrame, holdout_tics: set[int]) -> pd.Series:
    if "holdout" in df.columns:
        col = _as_bool(df["holdout"])
    else:
        col = pd.Series(False, index=df.index)
    if "tic" in df.columns and holdout_tics:
        col = col | df["tic"].astype(int).isin(holdout_tics)
    return col


def rejecting_mask(df: pd.DataFrame) -> pd.Series:
    """True when a rejecting vetting check failed. Ignores low_s and warning-only flags."""
    if "failed_checks" in df.columns:
        def bad(val):
            if val is None or (isinstance(val, float) and not math.isfinite(val)):
                return False
            parts = [p.strip() for p in str(val).split(",")]
            parts = [p for p in parts if p and p.lower() != "nan" and p not in WARNING_FLAGS]
            return len(parts) > 0
        return df["failed_checks"].map(bad)
    bad = pd.Series(False, index=df.index)
    for col in df.columns:
        if col.startswith("fail_"):
            name = col[len("fail_"):]
            if name in WARNING_FLAGS:
                continue
            bad = bad | _as_bool(df[col]).fillna(False)
    if "deep_flags" in df.columns:
        def _flags(val):
            if val is None or (isinstance(val, float) and not math.isfinite(val)):
                return False
            return bool({p for p in str(val).split(",") if p and p.lower() != "nan"} & DEEP_REJECT)
        bad = bad | df["deep_flags"].map(_flags)
    return bad


def rescore_injections(df: pd.DataFrame, s_star_of=None) -> pd.DataFrame:
    """Detection is matched and S >= S*. The stored ``detected`` column (S >= s_min) is not used."""
    out = df.copy()
    if "s_star" not in out.columns:
        if s_star_of is None:
            raise ValueError("s_star column or s_star_of is required")
        cams = out["camera"] if "camera" in out.columns else pd.Series([np.nan] * len(out), index=out.index)
        out["s_star"] = [
            s_star_of(int(s), None if pd.isna(c) else int(c))
            for s, c in zip(out["sector"], cams)
        ]
    matched = _as_bool(out["matched"]) if "matched" in out.columns else pd.Series(False, index=out.index)
    status = out["status"].astype(str).str.lower().eq("ok") if "status" in out.columns else True
    S = pd.to_numeric(out["S"], errors="coerce")
    thr = pd.to_numeric(out["s_star"], errors="coerce")
    out["detected_star"] = status & matched & S.notna() & thr.notna() & (S >= thr)
    out["passed_star"] = out["detected_star"] & ~rejecting_mask(out)
    return out


def _binned_map(df, pcol, rcol, det_col, pas_col):
    n = np.zeros((8, 8), dtype=int)
    det = np.zeros((8, 8), dtype=int)
    pas = np.zeros((8, 8), dtype=int)
    if len(df):
        pb = np.digitize(pd.to_numeric(df[pcol], errors="coerce"), P_EDGES) - 1
        rb = np.digitize(pd.to_numeric(df[rcol], errors="coerce"), R_EDGES) - 1
        for i in range(8):
            for j in range(8):
                m = (rb == i) & (pb == j)
                n[i, j] = int(m.sum())
                if n[i, j]:
                    det[i, j] = int(df.loc[df.index[m], det_col].sum())
                    pas[i, j] = int(df.loc[df.index[m], pas_col].sum())
    cells = []
    for i in range(8):
        row = []
        for j in range(8):
            lo, hi = wilson_interval(int(det[i, j]), int(n[i, j]))
            row.append({
                "n": int(n[i, j]), "detected": int(det[i, j]), "passed": int(pas[i, j]),
                "fraction": None if n[i, j] == 0 else float(det[i, j] / n[i, j]),
                "ci95_wilson": [lo, hi],
            })
        cells.append(row)
    fraction = np.divide(det, n, out=np.full(n.shape, np.nan, dtype=float), where=n > 0)
    return {"cells": cells, "fraction": fraction, "counts": n, "detected": det, "passed": pas}


def _summarize(S: np.ndarray):
    S = np.asarray(S, float)
    S = S[np.isfinite(S)]
    if S.size == 0:
        return None
    edges = np.linspace(0, 40, 81)
    counts, _ = np.histogram(np.clip(S, 0, 40), bins=edges)
    return {
        "n": int(S.size),
        "mean": float(S.mean()),
        "std": float(S.std(ddof=1)) if S.size > 1 else None,
        "min": float(S.min()),
        "max": float(S.max()),
        "quantiles": {str(q): float(np.quantile(S, q)) for q in (0.5, 0.9, 0.95, 0.99)},
        "histogram": {"edges": edges.tolist(), "counts": counts.astype(int).tolist()},
    }


def _two_sample(a, b, seed):
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    a, b = a[np.isfinite(a)], b[np.isfinite(b)]
    if a.size < 5 or b.size < 5:
        return None
    ks = stats.ks_2samp(a, b, method="auto")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ad = stats.anderson_ksamp([a, b])
    return {
        "n_invert": int(a.size), "n_scramble": int(b.size),
        "ks_statistic": float(ks.statistic), "ks_pvalue": float(ks.pvalue),
        "anderson_darling_statistic": float(ad.statistic),
        "anderson_darling_pvalue": None if ad.pvalue is None else float(ad.pvalue),
        "note": "Anderson-Darling p-values from scipy are capped at 0.25 unless a permutation method is used.",
    }


# ---------------------------------------------------------------- threshold config
def load_threshold(path: Path):
    raw = json.loads(Path(path).read_text())
    table = {}
    for sector, rec in (raw.get("sectors") or {}).items():
        if rec.get("s_star") is not None:
            table[(int(sector), None)] = float(rec["s_star"])
        for cam, crec in (rec.get("per_camera") or {}).items():
            if crec.get("s_star") is not None:
                table[(int(sector), int(cam))] = float(crec["s_star"])
    digest = hashlib.sha256(Path(path).read_bytes()).hexdigest()
    return table, raw.get("population") or {}, raw, digest


def resolve_s_star(sector, camera, table, override):
    if override is not None:
        return float(override), "cli"
    sector = int(sector)
    if camera is not None and not (isinstance(camera, float) and math.isnan(camera)):
        camera = int(camera)
        if (sector, camera) in table:
            return table[(sector, camera)], "per_camera"
    else:
        camera = None
    if (sector, None) in table:
        return table[(sector, None)], "sector"
    if (70, None) in table:
        return table[(70, None)], "borrowed_sector_70"
    return None, "missing"


def _attach_threshold(df, table, override):
    out = df.copy()
    vals, srcs = [], []
    cams = out["camera"] if "camera" in out.columns else pd.Series([np.nan] * len(out), index=out.index)
    for sector, camera in zip(out["sector"], cams):
        cam = None if pd.isna(camera) else int(camera)
        v, src = resolve_s_star(sector, cam, table, override)
        vals.append(v)
        srcs.append(src)
    out["s_star"] = vals
    out["s_star_source"] = srcs
    return out


def _quality_threshold(path: Path):
    raw = json.loads(Path(path).read_text())
    return float(raw["threshold"]), raw.get("rule")


def apply_quality(df, quality, threshold):
    """Drop stars with beta_max > threshold. Stars missing from the table are dropped."""
    q = quality.copy()
    q["tic"] = q["tic"].astype(int)
    q["sector"] = q["sector"].astype(int)
    if "drop" not in q.columns:
        q["drop"] = pd.to_numeric(q["beta_max"], errors="coerce") > threshold
    keep_cols = ["tic", "sector", "drop"]
    base = df.copy()
    base["tic"] = base["tic"].astype(int)
    base["sector"] = base["sector"].astype(int)
    merged = base.merge(q[keep_cols], on=["tic", "sector"], how="left")
    merged["drop"] = merged["drop"].fillna(True).astype(bool)
    return merged.loc[~merged["drop"]].drop(columns=["drop"]), int(merged["drop"].sum())


def _prepare(df, holdout_tics, holdout_only, quality, q_threshold):
    n_in = int(len(df))
    if "status" in df.columns:
        ok = df["status"].astype(str).str.lower().eq("ok")
        n_status = int((~ok).sum())
        df = df.loc[ok].copy()
    else:
        n_status = 0
        df = df.copy()
    n_quality = 0
    if quality is not None:
        df, n_quality = apply_quality(df, quality, q_threshold)
    held = holdout_mask(df, holdout_tics)
    if holdout_only:
        used = df.loc[held].copy()
        n_excl = int((~held).sum())
    else:
        used = df.loc[~held].copy()
        n_excl = int(held.sum())
    return used, {"n_input": n_in, "n_dropped_status": n_status, "n_dropped_quality": n_quality,
                  "n_excluded_holdout": n_excl, "n_used": int(len(used))}


def _column_snr(df):
    if not {"depth_inj", "n_tr", "cdpp_dur"}.issubset(df.columns):
        return pd.Series(np.nan, index=df.index, dtype=float)
    depth = pd.to_numeric(df["depth_inj"], errors="coerce")
    n_tr = pd.to_numeric(df["n_tr"], errors="coerce")
    cdpp = pd.to_numeric(df["cdpp_dur"], errors="coerce")
    snr = depth * np.sqrt(n_tr) / cdpp
    return snr.where((n_tr > 0) & (cdpp > 0) & np.isfinite(depth))


def _uniq(df, col):
    if df is None or col not in df.columns:
        return []
    return sorted({str(v) for v in df[col].dropna().unique()})


def _git_head():
    try:
        return subprocess.check_output(
            ["git", "-C", str(REPO_ROOT), "rev-parse", "--short=10", "HEAD"], text=True).strip()
    except Exception:
        return None


def _one(values):
    return values[0] if len(values) == 1 else None


def _bin_index(values, edges):
    return np.digitize(pd.to_numeric(values, errors="coerce"), edges) - 1


def _period_radius(df):
    pcol = next((c for c in ("period", "P", "P_inj") if c in df.columns), None)
    rcol = next((c for c in ("rp_earth", "rp", "rp_inj", "radius") if c in df.columns), None)
    return pcol, rcol


# ---------------------------------------------------------------- orchestration
def run_crosscheck(*, injections=None, nulls=None, products_root=None, candidates=None, quality=None,
                   out_dir, threshold_config=None, s_star=None, holdout_file=None, holdout_only=False,
                   quality_config=None, n_boot=2000, seed=20261110, n_stars=None, make_figures=True):
    out_dir = Path(out_dir)
    threshold_config = Path(threshold_config or PACKAGE_ROOT / "configs" / "interim_threshold_v1.json")
    quality_config = Path(quality_config or PACKAGE_ROOT / "configs" / "star_quality_v1.json")
    holdout_file = Path(holdout_file) if holdout_file else PACKAGE_ROOT / "holdout_tics.txt"
    table, population, thr_raw, thr_sha = load_threshold(threshold_config)
    q_threshold, q_rule = _quality_threshold(quality_config) if quality_config.exists() else (3.0, None)

    inj_paths, null_paths = [], []
    if products_root:
        found_i, found_n = discover_products(products_root)
        inj_paths, null_paths = found_i, found_n
    if injections:
        inj_paths = [Path(injections)]
    if nulls:
        null_paths = [Path(p) for p in nulls]
    if not inj_paths or not null_paths:
        raise FileNotFoundError("need injections and null products (paths or --products-root)")

    inj = pd.concat([read_table(p) for p in inj_paths], ignore_index=True)
    null_frames = []
    for p in null_paths:
        part = read_table(p)
        if "mode" not in part.columns or part["mode"].isna().all():
            part = part.copy()
            part["mode"] = _infer_mode(p)
        null_frames.append(part)
    null = pd.concat(null_frames, ignore_index=True)

    qdf = read_table(quality) if quality else None
    hold_tics = load_holdout_tics(holdout_file)
    inj_used, inj_counts = _prepare(inj, hold_tics, holdout_only, qdf, q_threshold)
    null_used, null_counts = _prepare(null, hold_tics, holdout_only, qdf, q_threshold)
    inj_used = _attach_threshold(inj_used, table, s_star)
    null_used = _attach_threshold(null_used, table, s_star) if len(null_used) else null_used
    if len(inj_used):
        inj_used = rescore_injections(inj_used)
        inj_used["snr_exp_recomputed"] = _column_snr(inj_used)
    else:
        inj_used["snr_exp_recomputed"] = []
        inj_used["detected_star"] = []
        inj_used["passed_star"] = []

    # --- completeness
    snr_diff = None
    if len(inj_used) and "snr_exp" in inj_used.columns:
        both = pd.to_numeric(inj_used["snr_exp"], errors="coerce")
        rec = inj_used["snr_exp_recomputed"]
        m = both.notna() & rec.notna()
        if m.any():
            snr_diff = float((both[m] - rec[m]).abs().max())
    gamma = {"shape": None, "scale": None, "note": "no injections in this split", "meaningful": False,
             "seed": int(seed), "n": 0}
    if len(inj_used):
        gamma = fit_gamma_cdf(
            inj_used["snr_exp_recomputed"].to_numpy(),
            inj_used["detected_star"].to_numpy().astype(float),
            inj_used["tic"].astype(int).to_numpy() if "tic" in inj_used.columns else np.arange(len(inj_used)),
            n_boot=n_boot, seed=seed,
        )
    maps = {}
    map_fraction = np.full((8, 8), np.nan)
    map_counts = np.zeros((8, 8))
    if len(inj_used) and {"P_inj", "rp_inj"}.issubset(inj_used.columns):
        overall = _binned_map(inj_used, "P_inj", "rp_inj", "detected_star", "passed_star")
        maps["all"] = {"cells": overall["cells"]}
        map_fraction, map_counts = overall["fraction"], overall["counts"]
        for sector, sub in inj_used.groupby(inj_used["sector"].astype(int)):
            maps[f"sector:{int(sector)}"] = {"cells": _binned_map(sub, "P_inj", "rp_inj", "detected_star", "passed_star")["cells"]}
    n_detected = int(inj_used["detected_star"].sum()) if len(inj_used) else 0
    n_passed = int(inj_used["passed_star"].sum()) if len(inj_used) else 0
    if len(inj_used) and "detected" in inj_used.columns:
        n_flag = int(_as_bool(inj_used["detected"]).sum())
    else:
        n_flag = None
    completeness = {
        **inj_counts,
        "n_detected": n_detected,
        "n_passed": n_passed,
        "n_detected_pipeline_flag": n_flag,
        "snr_exp_column_max_abs_diff": snr_diff,
        "gamma": gamma,
        "maps": maps,
        "p_edges": P_EDGES.tolist(),
        "r_edges": R_EDGES.tolist(),
        "interval": "Wilson score, 95%, on the detection fraction in each period-radius cell",
        "detection_rule": "matched period and epoch, and S >= S*(sector, camera). Aliases are not matches. The stored detected flag is S >= s_min and is not used.",
    }

    # --- null distributions
    def _values(frame):
        if frame is None or len(frame) == 0 or "S" not in frame.columns:
            return np.array([])
        return pd.to_numeric(frame["S"], errors="coerce").to_numpy()

    modes = sorted({str(m) for m in null_used["mode"].dropna().unique()}) if len(null_used) and "mode" in null_used.columns else []
    by_mode = {}
    collapsed = {}
    for mode in modes:
        sub = null_used[null_used["mode"].astype(str) == mode]
        by_mode[mode] = _summarize(_values(sub))
        if by_mode[mode] is not None:
            by_mode[mode]["n_stars"] = int(sub["tic"].astype(int).nunique()) if "tic" in sub.columns else int(len(sub))
            by_mode[mode]["n_realizations"] = int(sub["realization"].nunique()) if "realization" in sub.columns else 1
        if "tic" in sub.columns:
            collapsed[mode] = sub.groupby(sub["tic"].astype(int))["S"].mean().to_numpy()
        else:
            collapsed[mode] = _values(sub)
    for required in ("invert", "scramble"):
        by_mode.setdefault(required, None)
    combined_S = _values(null_used)
    comparison = {"invert_vs_scramble": None,
                  "note": "inverted and scrambled nulls are both required for the two-sample comparison"}
    if collapsed.get("invert") is not None and collapsed.get("scramble") is not None \
            and len(collapsed.get("invert", [])) and len(collapsed.get("scramble", [])):
        comparison = {
            "invert_vs_scramble": {
                "star_collapsed": {
                    **(_two_sample(collapsed["invert"], collapsed["scramble"], seed) or {}),
                    "definition": "one value per TIC: the mean of that star's per-realization maximum S",
                },
                "pooled_rows": {
                    **(_two_sample(_values(null_used[null_used["mode"].astype(str) == "invert"]),
                                   _values(null_used[null_used["mode"].astype(str) == "scramble"]), seed) or {}),
                    "definition": "every realization is a row; scramble realizations of one star are not independent",
                },
            },
            "note": "The star-collapsed test is the one to quote. The pooled-row test double-counts stars that have several scramble realizations.",
        }
    null_distribution = {
        **null_counts,
        "by_mode": by_mode,
        "combined": _summarize(combined_S),
        "comparison": comparison,
        "definition": "Each row is already the maximum S of one star in one null realization. Modes are whatever values the mode column (or the null directory name) identifies.",
    }

    # --- tail, per mode and combined, per sector
    def _n_stars(sector, frame):
        if n_stars is not None:
            return int(n_stars), "cli"
        rec = population.get(str(int(sector))) or {}
        if rec.get("kept"):
            return int(rec["kept"]), "interim_threshold_v1.population.kept"
        if len(frame) and "tic" in frame.columns:
            return int(frame["tic"].astype(int).nunique()), "null_sample_unique_tics"
        return max(int(len(frame)), 1), "null_sample_rows"

    def _fit_frame(frame, sector):
        ns, src = _n_stars(sector, frame)
        if len(frame) == 0:
            fit = fit_tail(np.array([]), np.array([]), n_stars=ns, n_boot=n_boot, seed=seed + 1, percentile=95.0)
        else:
            tics = frame["tic"].astype(int).to_numpy() if "tic" in frame.columns else np.arange(len(frame))
            fit = fit_tail(pd.to_numeric(frame["S"], errors="coerce").to_numpy(), tics,
                           n_stars=ns, n_boot=n_boot, seed=seed + 1, percentile=95.0)
        gumbel = fit_gumbel_maxima(
            pd.to_numeric(frame["S"], errors="coerce").to_numpy() if len(frame) else np.array([]),
            tics=frame["tic"].astype(int).to_numpy() if len(frame) and "tic" in frame.columns else None,
            n_stars=ns, n_boot=n_boot, seed=seed + 2)
        interim = None
        rec = (thr_raw.get("sectors") or {}).get(str(int(sector))) or {}
        if rec.get("s_star") is not None:
            interim = float(rec["s_star"])
        fit["gumbel"] = gumbel
        fit["n_stars_source"] = src
        fit["reference_interim_s_star"] = interim
        fit["reference_note"] = ("reference_interim_s_star is the interim layer's post-vetting GPD at u=6 on vetted "
                                 "inverted peaks. The fit in this file is a pre-vetting GPD above the 95th percentile "
                                 "of per-star max S. The two numbers answer different questions.")
        return fit

    tail_modes = {}
    frames_for_tail = {m: null_used[null_used["mode"].astype(str) == m] for m in modes} if len(null_used) else {}
    if len(null_used):
        frames_for_tail["combined"] = null_used
    for mode, frame in frames_for_tail.items():
        by_sector = {}
        if len(frame) == 0:
            tail_modes[mode] = None
            continue
        for sector, sub in frame.groupby(frame["sector"].astype(int)):
            by_sector[str(int(sector))] = _fit_frame(sub, int(sector))
        tail_modes[mode] = {"by_sector": by_sector}
    for required in ("invert", "scramble", "combined"):
        tail_modes.setdefault(required, None)
    tail = {
        "definition": ("Generalized Pareto fit to per-star maximum S above that sample's 95th percentile, "
                       "separately for each null mode and for the pooled rows. S* solves P(S > S*) = 1/N_stars "
                       "per sector. This is pre-vetting. It is not the interim post-vetting fit above u = 6."),
        "by_mode": tail_modes,
    }

    # --- ROC
    def _roc_block(pos_df, neg_df):
        if len(pos_df) == 0 or len(neg_df) == 0:
            return {"auc": None, "auc_ci95": None, "n_pos": int(len(pos_df)), "n_neg": int(len(neg_df)),
                    "fpr": [], "tpr": [], "note": "need both injections and null maxima"}
        matched = _as_bool(pos_df["matched"]) if "matched" in pos_df.columns else pd.Series(False, index=pos_df.index)
        pos_score = np.where(matched.to_numpy(), pd.to_numeric(pos_df["S"], errors="coerce").to_numpy(), 0.0)
        neg_score = pd.to_numeric(neg_df["S"], errors="coerce").to_numpy()
        pos_tic = pos_df["tic"].astype(int).to_numpy() if "tic" in pos_df.columns else np.arange(len(pos_df))
        neg_tic = neg_df["tic"].astype(int).to_numpy() if "tic" in neg_df.columns else np.arange(len(neg_df))
        curve = roc_curve(pos_score, neg_score)
        rng = np.random.default_rng(seed + 3)
        pg, ng = _star_groups(pos_tic), _star_groups(neg_tic)
        aucs = []
        for _ in range(int(n_boot)):
            if not pg or not ng:
                break
            ps = pos_score[_draw_groups(pg, rng)]
            ns = neg_score[_draw_groups(ng, rng)]
            a = mann_whitney_auc(ps, ns)
            if a is not None:
                aucs.append(a)
        curve["auc_ci95"] = _ci(aucs)
        curve["n_boot"] = int(n_boot)
        curve["n_boot_ok"] = int(len(aucs))
        curve["seed"] = int(seed + 3)
        matched_scores = pos_score[matched.to_numpy() & np.isfinite(pos_score)]
        curve["n_matched"] = int(matched_scores.size)
        curve["auc_matched_recoveries_only"] = mann_whitney_auc(matched_scores, neg_score[np.isfinite(neg_score)])
        curve["score"] = ("injection score is S when the recovery is matched, else 0; a negative is a null maximum S. "
                          "auc_matched_recoveries_only drops the unmatched injections instead of scoring them 0.")
        return curve

    roc = {"overall": _roc_block(inj_used, null_used) if len(null_used) else _roc_block(inj_used, null_used),
           "by_sector": {}, "by_sector_camera": {}}
    if len(inj_used) and len(null_used):
        sectors = sorted(set(inj_used["sector"].astype(int)).union(set(null_used["sector"].astype(int))))
        for sector in sectors:
            roc["by_sector"][str(sector)] = _roc_block(
                inj_used[inj_used["sector"].astype(int) == sector],
                null_used[null_used["sector"].astype(int) == sector])
            if "camera" in inj_used.columns and "camera" in null_used.columns:
                cams = sorted(set(inj_used.loc[inj_used["sector"].astype(int) == sector, "camera"].dropna().astype(int))
                              .union(null_used.loc[null_used["sector"].astype(int) == sector, "camera"].dropna().astype(int)))
                for cam in cams:
                    roc["by_sector_camera"][f"{sector}:{cam}"] = _roc_block(
                        inj_used[(inj_used["sector"].astype(int) == sector) & (inj_used["camera"].astype(int) == cam)],
                        null_used[(null_used["sector"].astype(int) == sector) & (null_used["camera"].astype(int) == cam)])

    # --- reliability
    cand_df = None
    if candidates:
        cand_df = read_table(candidates)
        if "tic" in cand_df.columns:
            held = holdout_mask(cand_df, hold_tics)
            cand_df = cand_df.loc[held if holdout_only else ~held].copy()

    def _obs_counts(p_lo, p_hi, r_lo, r_hi, sector=None):
        if cand_df is None or len(cand_df) == 0:
            return None, None
        pcol, rcol = _period_radius(cand_df)
        disp_col = next((c for c in ("disposition", "label", "disp") if c in cand_df.columns), None)
        if pcol is None or rcol is None or disp_col is None:
            return None, None
        sub = cand_df
        if sector is not None and "sector" in sub.columns:
            sub = sub[sub["sector"].astype(int) == int(sector)]
        if p_lo is None:
            m = pd.Series(True, index=sub.index)
        else:
            p = pd.to_numeric(sub[pcol], errors="coerce")
            r = pd.to_numeric(sub[rcol], errors="coerce")
            m = (p >= p_lo) & (p < p_hi) & (r >= r_lo) & (r < r_hi)
        labels = sub.loc[m, disp_col].astype(str).str.strip().str.lower()
        return int(labels.isin(FP_LABELS).sum()), int(labels.isin(PC_LABELS).sum())

    def _sim_counts(frame):
        if frame is None or len(frame) == 0:
            return 0, 0
        S = pd.to_numeric(frame["S"], errors="coerce")
        thr = pd.to_numeric(frame["s_star"], errors="coerce")
        above = S.notna() & thr.notna() & (S >= thr)
        if "checks_logged" in frame.columns:
            logged = _as_bool(frame["checks_logged"])
        else:
            logged = pd.Series(True, index=frame.index)
        evaluable = above & logged
        rejected = evaluable & rejecting_mask(frame)
        return int(rejected.sum()), int(evaluable.sum())

    rel_rows = []
    n_rej, n_sim = _sim_counts(null_used)
    n_fp, n_pc = _obs_counts(None, None, None, None)
    overall = reliability_from_counts(n_fp, n_pc, n_rej, n_sim)
    # per period-radius cell, counted when the null has a radius
    if len(null_used):
        pcol, rcol = _period_radius(null_used)
    else:
        pcol = rcol = None
    if pcol and rcol and len(null_used):
        pb = _bin_index(null_used[pcol], P_EDGES)
        rb = _bin_index(null_used[rcol], R_EDGES)
        for i in range(8):
            for j in range(8):
                m = (rb == i) & (pb == j)
                if not np.any(m):
                    continue
                sub = null_used.loc[null_used.index[m]]
                kr, ns = _sim_counts(sub)
                if ns == 0 and cand_df is None:
                    continue
                fp, pc = _obs_counts(float(P_EDGES[j]), float(P_EDGES[j + 1]), float(R_EDGES[i]), float(R_EDGES[i + 1]))
                cell = reliability_from_counts(fp, pc, kr, ns)
                cell["period"] = [float(P_EDGES[j]), float(P_EDGES[j + 1])]
                cell["rp_earth"] = [float(R_EDGES[i]), float(R_EDGES[i + 1])]
                if ns or (fp or pc):
                    rel_rows.append(cell)
    reliability = {
        "formula": "R = 1 - (N_FP,obs / N_PC,obs) * (1 - E) / E",
        "citation": "Thompson, S. E., et al. 2018, ApJS, 235, 38",
        "E_definition": ("E is the fraction of simulated nulls with S >= S* and checks logged that fail a rejecting "
                         "vetting check. It is computed only from inverted/scrambled runs. low_s is ignored because "
                         "the threshold is re-applied as S*. Warning flags (ruwe_high, gaia_nss) are not failures."),
        "candidates_path": None if not candidates else str(candidates),
        "overall": overall,
        "bins": rel_rows,
        "p_edges": P_EDGES.tolist(),
        "r_edges": R_EDGES.tolist(),
        "uncertainty": "Wilson 95% interval on E, transformed through R. N_FP,obs and N_PC,obs are treated as fixed counts.",
    }

    # --- provenance and the small-sample judgement
    pv = sorted(set(_uniq(inj, "pipeline_version") + _uniq(null, "pipeline_version")))
    ch = sorted(set(_uniq(inj, "config_hash") + _uniq(null, "config_hash")))
    gh = sorted(set(_uniq(inj, "git_hash") + _uniq(null, "git_hash")))
    k_tail = 0
    for mode_block in tail_modes.values():
        if not mode_block:
            continue
        for fit in mode_block["by_sector"].values():
            k_tail = max(k_tail, int(fit.get("k_tail") or 0))
    reasons = []
    if inj_counts["n_used"] < 200:
        reasons.append(f"only {inj_counts['n_used']} injections in the fit (a gamma CDF wants several hundred, with both outcomes well populated)")
    if null_counts["n_used"] < 1000:
        reasons.append(f"only {null_counts['n_used']} null rows (the 95th-percentile tail is then a handful of points, not a measurement)")
    if k_tail < 50:
        reasons.append(f"the largest GPD tail here has {k_tail} exceedances; 50 is a bare minimum for a quoted shape")
    if "scramble" not in modes or "invert" not in modes:
        reasons.append("one of the two null experiments is absent, so the inverted-versus-scrambled comparison cannot be made")
    small = bool(reasons)
    if small:
        small_note = ("These inputs are too small for a meaningful measurement; the numbers below are a demonstration "
                      "that the code path runs, not a measurement of the survey. " + " ".join(reasons))
    else:
        small_note = "Sample size meets the module's minimum for quoting the gamma CDF, the GPD tail, and both null experiments."

    used_thresholds = []
    if len(inj_used) and "s_star" in inj_used.columns:
        for (sector, camera, source), sub in inj_used.groupby(["sector", "camera", "s_star_source"], dropna=False):
            used_thresholds.append({
                "sector": int(sector),
                "camera": None if pd.isna(camera) else int(camera),
                "s_star": None if pd.isna(sub["s_star"].iloc[0]) else float(sub["s_star"].iloc[0]),
                "source": str(source),
                "n_injections": int(len(sub)),
            })

    assumptions = [
        "Expected S is recomputed as depth_inj * sqrt(n_tr) / cdpp_dur. The stored snr_exp column is only a cross-check.",
        "A detection is a period-and-epoch match with S >= S*(sector, camera). S* is the per-camera value in interim_threshold_v1.json when that value is present, else the sector value, else the sector-70 value. --s-star overrides all of them. s_min is not a detection threshold.",
        "Holdout is the union of the holdout column and holdout_tics.txt. The default split is the complement. --holdout-only fits the held-out TICs and nothing else.",
        "The star-quality cut is applied only when a quality table is passed. The threshold is read from star_quality_v1.json and is not re-tuned. Stars absent from that table are excluded.",
        "The gamma CDF is two-parameter and asymptotes to 1. The period-radius grid is the fixed 8x8 log grid from exosearch.reliability, with Wilson 95% intervals.",
        "The GPD is fit above the sample 95th percentile of pre-vetting per-star max S. N_stars is population.kept from the interim config when that sector is listed, so the target rate is one false alarm per searched star, not one per null row.",
        "N_eff assumes a single BLS trial is exponential with the Gumbel scale. It is a cross-check, not a frequency count.",
        "E uses only simulated nulls. R uses observed PC and FP counts from --candidates when that file is present.",
    ]
    payload = {
        "seed": int(seed),
        "component_seeds": {"gamma": int(seed), "gpd": int(seed + 1), "gumbel": int(seed + 2), "roc": int(seed + 3)},
        "pipeline_version": _one(pv),
        "config_hash": _one(ch),
        "git_hash": _one(gh),
        "pipeline_version_values": pv,
        "config_hash_values": ch,
        "git_hash_values": gh,
        "analysis_git_hash": _git_head(),
        "holdout": {"mode": "holdout" if holdout_only else "tuning", "file": str(holdout_file),
                    "n_file": len(hold_tics)},
        "small_sample": small,
        "small_sample_note": small_note,
        "quality_cut": {"applied": qdf is not None, "threshold": q_threshold, "rule": q_rule,
                        "config": str(quality_config)},
        "threshold": {"config": str(threshold_config), "sha256": thr_sha,
                      "override": s_star, "used": used_thresholds},
        "inputs": {
            "injections": [str(p) for p in inj_paths],
            "nulls": [str(p) for p in null_paths],
            "products_root": None if not products_root else str(products_root),
            "null_modes": modes,
            "candidates": None if not candidates else str(candidates),
        },
        "assumptions": assumptions,
        "completeness": completeness,
        "null_distribution": null_distribution,
        "tail": tail,
        "roc": roc,
        "reliability": reliability,
    }

    # figure bundle: the stability panel follows the invert sample when it exists, else the pool
    headline = None
    for key in ("invert", "combined"):
        block = tail_modes.get(key)
        if block and block.get("by_sector"):
            headline = max(block["by_sector"].values(), key=lambda f: f.get("n") or 0)
            break
    roc_curves = []
    if roc["overall"].get("fpr"):
        roc_curves.append({"name": "all", "fpr": roc["overall"]["fpr"], "tpr": roc["overall"]["tpr"],
                           "auc": roc["overall"].get("auc")})
    if len(roc["by_sector"]) > 1:
        for key, block in list(roc["by_sector"].items())[:6]:
            if block.get("fpr"):
                roc_curves.append({"name": f"S{key}", "fpr": block["fpr"], "tpr": block["tpr"], "auc": block.get("auc")})
    rel_plot = []
    if overall.get("E") is not None:
        rel_plot.append({"label": "all simulated nulls", "E": overall["E"],
                         "lo": overall["E_ci95"][0], "hi": overall["E_ci95"][1]})
    for cell in rel_rows[:12]:
        if cell.get("E") is None:
            continue
        rel_plot.append({"label": f"P {cell['period'][0]:.2g}–{cell['period'][1]:.2g}, R {cell['rp_earth'][0]:.2g}–{cell['rp_earth'][1]:.2g}",
                         "E": cell["E"], "lo": cell["E_ci95"][0], "hi": cell["E_ci95"][1]})
    bundle = {
        "snr": inj_used["snr_exp_recomputed"].to_numpy() if len(inj_used) else [],
        "detected": inj_used["detected_star"].to_numpy().astype(float) if len(inj_used) else [],
        "gamma": gamma,
        "map_fraction": map_fraction,
        "map_counts": map_counts,
        "p_edges": P_EDGES,
        "r_edges": R_EDGES,
        "null_by_mode": (
            {m: _values(null_used[null_used["mode"].astype(str) == m]) if len(null_used) else [] for m in modes}
            | ({"combined": combined_S} if len(modes) > 1 else {})
        ),
        "stability": (headline or {}).get("stability", {}).get("thresholds"),
        "mrl": (headline or {}).get("mrl"),
        "stability_title": "Threshold stability of the largest null sample",
        "roc_curves": roc_curves,
        "reliability_rows": rel_plot,
    }
    prov_keys = ["seed", "pipeline_version", "config_hash", "git_hash", "analysis_git_hash",
                 "pipeline_version_values", "config_hash_values", "git_hash_values", "small_sample", "small_sample_note"]
    prov = {k: payload[k] for k in prov_keys}
    dump_json(out_dir / "crosscheck.json", payload)
    dump_json(out_dir / "completeness.json", {**prov, "completeness": completeness})
    dump_json(out_dir / "null_distribution.json", {**prov, "null_distribution": null_distribution})
    dump_json(out_dir / "tail.json", {**prov, "tail": tail})
    dump_json(out_dir / "roc.json", {**prov, "roc": roc})
    dump_json(out_dir / "reliability.json", {**prov, "reliability": reliability})
    if make_figures:
        figures.write_all(out_dir, bundle)
    return payload


def build_parser():
    p = argparse.ArgumentParser(
        prog="analysis",
        description="Independent cross-check of exosearch completeness, null tails, ROC and reliability.",
    )
    p.add_argument("--injections", type=Path, help="injection CSV or parquet (one file)")
    p.add_argument("--null", dest="nulls", type=Path, action="append", help="null max-S CSV or parquet; repeatable")
    p.add_argument("--products-root", type=Path, help="directory that contains products/v<version>/ or the version directory itself")
    p.add_argument("--candidates", type=Path, help="observed vetting outcomes; R is left null when this is omitted")
    p.add_argument("--quality", type=Path, help="optional star-quality table with tic, sector, and beta_max or drop")
    p.add_argument("--threshold-config", type=Path, default=PACKAGE_ROOT / "configs" / "interim_threshold_v1.json")
    p.add_argument("--s-star", type=float, default=None, help="use this S* for every sector and camera")
    p.add_argument("--holdout-file", type=Path, default=PACKAGE_ROOT / "holdout_tics.txt")
    p.add_argument("--holdout-only", action="store_true", help="fit and evaluate the held-out TICs only")
    p.add_argument("--quality-config", type=Path, default=PACKAGE_ROOT / "configs" / "star_quality_v1.json")
    p.add_argument("--n-stars", type=int, default=None, help="override N_stars in the false-alarm rate 1/N_stars")
    p.add_argument("--n-boot", type=int, default=2000)
    p.add_argument("--seed", type=int, default=20261110)
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--no-plots", action="store_true")
    return p


def main(argv=None):
    args = build_parser().parse_args(argv)
    result = run_crosscheck(
        injections=args.injections, nulls=args.nulls, products_root=args.products_root,
        candidates=args.candidates, quality=args.quality, out_dir=args.out,
        threshold_config=args.threshold_config, s_star=args.s_star,
        holdout_file=args.holdout_file, holdout_only=args.holdout_only,
        quality_config=args.quality_config, n_boot=args.n_boot, seed=args.seed,
        n_stars=args.n_stars, make_figures=not args.no_plots,
    )
    gamma = result["completeness"]["gamma"]
    print(f"wrote {args.out}")
    print(f"split={result['holdout']['mode']} injections={result['completeness']['n_used']} "
          f"detected={result['completeness']['n_detected']} nulls={result['null_distribution']['n_used']}")
    print(f"gamma shape={gamma.get('shape')} scale={gamma.get('scale')}")
    print(f"small_sample={result['small_sample']}")
    print(result["small_sample_note"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
