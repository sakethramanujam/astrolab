"""Second-pass prewhitening validation (v1.1 amendment 2026-11-11, Section 2b).

New seeds: base seed 20261112, injection runs 12-15, hosts redrawn with the new base seed.
Disjoint from validate.py SEED 20261110 / run 0.

Arms:
- quiet stars with prewhitening and without (runs 12-15)
- the 25 Gaia DR3 stars in docs/preregistration/dsct_s70_2026-11-11/dsct_s70.csv:
    * real search with prewhitening
    * injections with prewhitening

Classify each Gaia star by dominant Lomb-Scargle frequency of quality-masked normalised
sector-70 flux before prewhitening:
  > 5 c/d: delta Scuti
  <= 5 c/d: gamma Doradus

Pass criteria (all three must hold):
  (i) at least 90% of Gaia stars have top search peak that is not a pulsation frequency or harmonic;
  (ii) injection recovery on Gaia stars within 5 pp of quiet stars at same S (bin by expected S; report per bin and overall with binomial CIs);
  (iii) quiet-star recovery with prewhitening drops by no more than 2 pp vs without.
Report per group (delta Scuti, gamma Dor) and total.
Output JSON + CSV.
"""
from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMBA_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[_v] = "1"

import argparse
import hashlib
import json
from pathlib import Path
from typing import Dict, Any, List, Optional, Tuple, Union

import numpy as np
import pandas as pd
from scipy import stats

from .. import batch, lightcurve, targets
from . import run, settings

BASE_SEED = 20261112
INJECTION_RUNS = [12, 13, 14, 15]
EXCLUDED_BASE_SEED = 20261110
EXCLUDED_RUNS = [0]
GAIA_CSV_PATH = Path("docs/preregistration/dsct_s70_2026-11-11/dsct_s70.csv")

S_REF = 10.0
DEFAULT_S_BINS = [0.0, 10.0, 15.0, 25.0, 50.0, 500.0]


def seed_for_pw2(base_seed: int, sector: int, run: int, tic: int) -> int:
    return int.from_bytes(hashlib.sha256(f"{base_seed}:{sector}:{run}:{tic}".encode()).digest()[:4], "big")


def classify_pulsator_dominant_frequency(t: np.ndarray, f_norm: np.ndarray) -> Tuple[str, float]:
    """Classify star by dominant Lomb-Scargle frequency of quality-masked normalised flux before prewhitening:
    above 5 c/d is delta Scuti, 5 c/d or below is gamma Dor.
    """
    from astropy.timeseries import LombScargle
    baseline = float(np.max(t) - np.min(t))
    df = 1.0 / (5.0 * baseline)
    fmax = 1.0 / (2.0 * float(np.median(np.diff(t))))
    freqs = np.arange(0.1, fmax, df)
    power = LombScargle(t, f_norm, fit_mean=True, center_data=True, normalization="psd").power(
        freqs, method="fast", assume_regular_frequency=True
    )
    idx_max = int(np.argmax(power))
    f_dom = float(freqs[idx_max])
    grp = "delta_scuti" if f_dom > 5.0 else "gamma_dor"
    return grp, f_dom


def find_pulsation_frequencies(
    t: np.ndarray,
    f_norm: np.ndarray,
    snr_stop: float = 4.0,
    max_terms: int = 20,
    window_cpd: float = 2.0,
    fmin_cpd: float = 0.1,
) -> List[float]:
    """Extract prominent pulsation frequencies before prewhitening using Lomb-Scargle iterative search."""
    from astropy.timeseries import LombScargle
    baseline = float(np.max(t) - np.min(t))
    df = 1.0 / (5.0 * baseline)
    fnyq = 1.0 / (2.0 * float(np.median(np.diff(t))))
    freqs = np.arange(fmin_cpd, fnyq, df)

    r = f_norm.copy()
    peaks = []
    for _ in range(max_terms):
        p = LombScargle(t, r, fit_mean=True, center_data=True, normalization="psd").power(
            freqs, method="fast", assume_regular_frequency=True
        )
        amp = np.sqrt(np.maximum(4.0 * p / t.size, 0.0))
        imax = int(np.argmax(amp))
        f_cand = freqs[imax]
        # Noise estimate
        m = (np.abs(freqs - f_cand) <= window_cpd) & (np.abs(freqs - f_cand) > 3 * df)
        noise = float(np.mean(amp[m])) if m.any() else np.median(amp)
        snr = amp[imax] / max(noise, 1e-12)
        if snr < snr_stop:
            break
        peaks.append(float(f_cand))
        # Subtract sinusoid
        w = 2 * np.pi * f_cand * t
        X = np.stack([np.sin(w), np.cos(w), np.ones_like(t)], axis=1)
        try:
            coef = np.linalg.solve(X.T @ X, X.T @ r)
            r -= (coef[0] * np.sin(w) + coef[1] * np.cos(w))
        except Exception:
            break
    return peaks


def is_pulsation_or_harmonic(
    f_found: float,
    pulsation_freqs: List[float],
    freq_tol: float = 0.01,
    harmonic_tol: float = 0.01,
    harmonics: Tuple[float, ...] = (1.0, 2.0, 3.0, 0.5, 0.3333333333333333),
) -> bool:
    """Check if search peak frequency is a known pulsation frequency or harmonic."""
    if f_found <= 0:
        return False
    for f_puls in pulsation_freqs:
        if f_puls <= 0:
            continue
        # Direct match
        if abs(f_found / f_puls - 1.0) < freq_tol:
            return True
        # Harmonic match
        for h in harmonics:
            if abs(f_found / (h * f_puls) - 1.0) < harmonic_tol:
                return True
    return False


def binomial_ci_wilson(k: int, n: int, confidence: float = 0.95) -> Tuple[float, float]:
    """Wilson score binomial confidence interval."""
    if n == 0:
        return (0.0, 0.0)
    z = stats.norm.ppf(1.0 - (1.0 - confidence) / 2.0)
    p = k / n
    denominator = 1.0 + z ** 2 / n
    centre = (p + z ** 2 / (2.0 * n)) / denominator
    half_width = z * np.sqrt((p * (1.0 - p) + z ** 2 / (4.0 * n)) / n) / denominator
    return (float(max(0.0, centre - half_width)), float(min(1.0, centre + half_width)))


def evaluate_pw2_report(
    gaia_results: List[Dict[str, Any]],
    gaia_injections: List[Dict[str, Any]],
    quiet_injections: List[Dict[str, Any]],
    cfg: Optional[Dict[str, Any]] = None,
    s_bins: Optional[List[float]] = None,
    unclean_reasons: Optional[Dict[int, str]] = None,
) -> Dict[str, Any]:
    """Evaluate all three pass criteria:
    (i) at least 90% of Gaia stars have top search peak that is not a pulsation frequency or harmonic;
    (ii) injection recovery on Gaia stars within 5 percentage points of quiet stars at same S;
    (iii) quiet-star recovery with prewhitening drops by no more than 2 points vs without.
    Report per group (delta Scuti, gamma Dor, unclassified) and total, with binomial CIs.
    """
    if cfg is None:
        cfg = settings.cfg()
    pw2_cfg = cfg.get("validation_pw2", {})
    freq_tol = float(pw2_cfg.get("frequency_match_tol", 0.01))
    harm_tol = float(pw2_cfg.get("harmonic_match_tol", 0.01))
    harmonics = tuple(pw2_cfg.get("harmonics", [1, 2, 3, 0.5, 0.3333333333333333]))

    if s_bins is None:
        s_bins = DEFAULT_S_BINS

    # Criterion (i): Real search on Gaia stars
    gaia_rows = []
    for g in gaia_results:
        tic = int(g["tic"])
        grp = g.get("group", "unclassified")
        f_dom = float(g.get("f_dom", 0.0))
        p_found = float(g.get("period", 0.0))
        f_found = 1.0 / p_found if p_found > 0 else 0.0
        puls_freqs = g.get("pulsation_frequencies", [])
        
        # Check if already marked unclean or missing
        reason = g.get("unclean_reason")
        if unclean_reasons and tic in unclean_reasons:
            reason = unclean_reasons[tic]

        if reason is not None:
            clean_peak = False
            is_puls = (reason == "pulsation_harmonic")
        else:
            is_puls = is_pulsation_or_harmonic(f_found, puls_freqs, freq_tol, harm_tol, harmonics)
            clean_peak = not is_puls
            if not clean_peak:
                reason = "pulsation_harmonic"

        gaia_rows.append(dict(
            tic=tic,
            group=grp,
            f_dom=f_dom,
            P_found=p_found,
            f_found=f_found,
            pulsation_frequencies=puls_freqs,
            is_pulsation=is_puls,
            clean_peak=clean_peak,
            unclean_reason=reason if not clean_peak else "",
        ))

    if gaia_rows:
        df_gaia_real = pd.DataFrame(gaia_rows)
    else:
        df_gaia_real = pd.DataFrame(columns=[
            "tic", "group", "f_dom", "P_found", "f_found", "pulsation_frequencies", "is_pulsation", "clean_peak", "unclean_reason"
        ])

    crit1_groups = {}
    if not df_gaia_real.empty and "group" in df_gaia_real.columns:
        for grp_name, sub in df_gaia_real.groupby("group"):
            n_tot = len(sub)
            n_clean = int(sub.clean_peak.sum())
            rate = float(n_clean / n_tot) if n_tot > 0 else 0.0
            crit1_groups[grp_name] = {
                "n_total": n_tot,
                "n_clean": n_clean,
                "clean_rate": rate,
                "ci95": binomial_ci_wilson(n_clean, n_tot),
            }

    total_gaia = len(df_gaia_real)
    total_clean = int(df_gaia_real.clean_peak.sum()) if total_gaia > 0 else 0
    total_clean_rate = float(total_clean / total_gaia) if total_gaia > 0 else 0.0
    # At least 90% of the full registered 25 Gaia stars must be clean; fails if empty/incomplete
    crit1_pass = bool(total_gaia >= 25 and total_clean_rate >= 0.90)

    # Criterion (ii): Injection recovery on Gaia stars vs quiet stars at same expected S
    df_gaia_inj = pd.DataFrame(gaia_injections) if gaia_injections else pd.DataFrame(columns=["tic", "expected_S", "recovered", "group"])
    df_quiet_inj = pd.DataFrame(quiet_injections) if quiet_injections else pd.DataFrame(columns=["tic", "expected_S", "recovered", "prewhiten"])

    if not df_quiet_inj.empty and "prewhiten" in df_quiet_inj.columns:
        df_quiet_pw = df_quiet_inj[df_quiet_inj.prewhiten == True].copy()
        df_quiet_nopw = df_quiet_inj[df_quiet_inj.prewhiten == False].copy()
    else:
        df_quiet_pw = df_quiet_inj.copy()
        df_quiet_nopw = pd.DataFrame()

    crit2_bins = []
    crit2_pass = False

    # Helper function to compute bin comparison
    def _compute_bins_comparison(df_target, df_ref):
        bins_out = []
        df_t = df_target.copy()
        df_r = df_ref.copy()
        df_t["s_bin"] = pd.cut(df_t["expected_S"], s_bins, right=False)
        df_r["s_bin"] = pd.cut(df_r["expected_S"], s_bins, right=False)
        all_bins = sorted(
            set(df_t["s_bin"].dropna().unique()).union(set(df_r["s_bin"].dropna().unique())),
            key=lambda x: x.left,
        )
        b_pass = bool(len(all_bins) > 0)
        for b in all_bins:
            sub_t = df_t[df_t.s_bin == b]
            sub_r = df_r[df_r.s_bin == b]
            nt, nr = len(sub_t), len(sub_r)
            rt = float(sub_t.recovered.mean()) if nt > 0 else 0.0
            rr = float(sub_r.recovered.mean()) if nr > 0 else 0.0
            diff_pp = rt - rr
            ci_t = binomial_ci_wilson(int(sub_t.recovered.sum()), nt) if nt > 0 else (0.0, 0.0)
            ci_r = binomial_ci_wilson(int(sub_r.recovered.sum()), nr) if nr > 0 else (0.0, 0.0)
            passed = bool(diff_pp >= -0.05) if (nt >= 5 and nr >= 5) else True
            if nt >= 5 and nr >= 5 and not passed:
                b_pass = False
            bins_out.append(dict(
                bin=str(b),
                n_target=nt,
                rec_target=rt,
                ci95_target=ci_t,
                n_quiet=nr,
                rec_quiet=rr,
                ci95_quiet=ci_r,
                diff_pp=diff_pp,
                pass_=passed,
            ))
        return bins_out, b_pass

    crit2_groups = {}
    if not df_gaia_inj.empty and not df_quiet_pw.empty:
        # Total Gaia injection recovery vs quiet pw
        crit2_bins, crit2_pass = _compute_bins_comparison(df_gaia_inj, df_quiet_pw)

        # Per-group breakdown
        if "group" in df_gaia_inj.columns:
            for grp_name, sub_g in df_gaia_inj.groupby("group"):
                grp_bins, grp_b_pass = _compute_bins_comparison(sub_g, df_quiet_pw)
                g_sub_s10 = sub_g[sub_g.expected_S >= S_REF]
                q_sub_s10 = df_quiet_pw[df_quiet_pw.expected_S >= S_REF]
                g_grp_rec = float(g_sub_s10.recovered.mean()) if len(g_sub_s10) > 0 else 0.0
                q_grp_rec = float(q_sub_s10.recovered.mean()) if len(q_sub_s10) > 0 else 0.0
                grp_diff = g_grp_rec - q_grp_rec
                grp_ci = binomial_ci_wilson(int(g_sub_s10.recovered.sum()), len(g_sub_s10)) if len(g_sub_s10) > 0 else (0.0, 0.0)
                crit2_groups[grp_name] = {
                    "overall_rec": g_grp_rec,
                    "overall_ci95": grp_ci,
                    "overall_diff_pp": grp_diff,
                    "overall_pass": bool(len(g_sub_s10) > 0 and len(q_sub_s10) > 0 and grp_diff >= -0.05),
                    "bins": grp_bins,
                }

        # Overall S >= 10 comparison
        g_s10 = df_gaia_inj[df_gaia_inj.expected_S >= S_REF]
        q_s10 = df_quiet_pw[df_quiet_pw.expected_S >= S_REF]
        ng_s10 = len(g_s10)
        nq_s10 = len(q_s10)
        g_rec_tot = float(g_s10.recovered.mean()) if ng_s10 > 0 else 0.0
        q_rec_tot = float(q_s10.recovered.mean()) if nq_s10 > 0 else 0.0
        g_ci_tot = binomial_ci_wilson(int(g_s10.recovered.sum()), ng_s10) if ng_s10 > 0 else (0.0, 0.0)
        q_ci_tot = binomial_ci_wilson(int(q_s10.recovered.sum()), nq_s10) if nq_s10 > 0 else (0.0, 0.0)
        diff_tot = g_rec_tot - q_rec_tot
        crit2_overall_pass = bool(ng_s10 > 0 and nq_s10 > 0 and diff_tot >= -0.05)
        crit2_pass = bool(crit2_pass and crit2_overall_pass)
    else:
        crit2_overall_pass = False
        diff_tot = 0.0
        g_rec_tot, q_rec_tot = 0.0, 0.0
        g_ci_tot, q_ci_tot = (0.0, 0.0), (0.0, 0.0)
        ng_s10, nq_s10 = 0, 0

    # Criterion (iii): Quiet-star recovery pw vs nopw drops by <= 2 points
    if not df_quiet_pw.empty and not df_quiet_nopw.empty:
        q_pw_s10 = df_quiet_pw[df_quiet_pw.expected_S >= S_REF]
        q_nopw_s10 = df_quiet_nopw[df_quiet_nopw.expected_S >= S_REF]
        n_pw = len(q_pw_s10)
        n_nopw = len(q_nopw_s10)
        r_pw = float(q_pw_s10.recovered.mean()) if n_pw > 0 else 0.0
        r_nopw = float(q_nopw_s10.recovered.mean()) if n_nopw > 0 else 0.0
        ci_pw = binomial_ci_wilson(int(q_pw_s10.recovered.sum()), n_pw) if n_pw > 0 else (0.0, 0.0)
        ci_nopw = binomial_ci_wilson(int(q_nopw_s10.recovered.sum()), n_nopw) if n_nopw > 0 else (0.0, 0.0)
        quiet_drop = r_nopw - r_pw
        crit3_pass = bool(n_pw > 0 and n_nopw > 0 and quiet_drop <= 0.02)
    else:
        r_pw, r_nopw, quiet_drop = 0.0, 0.0, 0.0
        ci_pw, ci_nopw = (0.0, 0.0), (0.0, 0.0)
        n_pw, n_nopw = 0, 0
        crit3_pass = False

    overall_pass = bool(crit1_pass and crit2_pass and crit3_pass)

    return {
        "pass": overall_pass,
        "note": "25 Gaia stars is coarse (one star = 4 percentage points).",
        "criterion_1_real_gaia_search": {
            "pass": crit1_pass,
            "threshold": 0.90,
            "total_stars": total_gaia,
            "total_clean": total_clean,
            "total_clean_rate": total_clean_rate,
            "total_clean_ci95": binomial_ci_wilson(total_clean, total_gaia),
            "groups": crit1_groups,
            "real_gaia_table": df_gaia_real.to_dict(orient="records"),
        },
        "criterion_2_gaia_vs_quiet_injections": {
            "pass": crit2_pass,
            "max_drop_pp": 0.05,
            "overall_diff_pp": diff_tot,
            "overall_pass": crit2_overall_pass,
            "overall_gaia_rec": g_rec_tot,
            "overall_gaia_ci95": g_ci_tot,
            "overall_quiet_rec": q_rec_tot,
            "overall_quiet_ci95": q_ci_tot,
            "groups": crit2_groups,
            "bins": crit2_bins,
        },
        "criterion_3_quiet_prewhiten_drop": {
            "pass": crit3_pass,
            "max_allowed_drop_pp": 0.02,
            "recovery_pw": r_pw,
            "ci95_pw": ci_pw,
            "recovery_nopw": r_nopw,
            "ci95_nopw": ci_nopw,
            "drop_pp": quiet_drop,
        },
    }


def classify_gaia_star_from_file(lc_path: Union[str, Path]) -> Dict[str, Any]:
    """Load light curve, quality-mask, normalize flux, classify dominant frequency and find pulsation terms."""
    p = Path(lc_path)
    if not p.exists():
        raise FileNotFoundError(f"Light curve not found: {p}")
    t, f, e, q, meta = lightcurve.load_spoc_raw(p)
    good = lightcurve.good_mask(f, e, q)
    if np.sum(good) < 500:
        raise ValueError(f"Too few good cadences: {np.sum(good)}")
    tg = t[good]
    med = np.nanmedian(f[good])
    fg = f[good] / med

    grp, f_dom = classify_pulsator_dominant_frequency(tg, fg)
    puls_freqs = find_pulsation_frequencies(tg, fg)
    return {
        "tic": int(meta.get("tic") or 0),
        "group": grp,
        "f_dom": f_dom,
        "pulsation_frequencies": puls_freqs,
    }


def build_pw2_manifests(
    weights: List[Tuple[str, float]],
    out: Union[str, Path],
    sector: int = 70,
    base_seed: int = BASE_SEED,
    injection_runs: List[int] = INJECTION_RUNS,
    gaia_csv: Path = GAIA_CSV_PATH,
) -> List[Dict[str, Any]]:
    """Build task manifests for prewhitening validation on new seeds."""
    from .. import inject, store
    con = store.connect()
    cfg_inj = json.loads((settings.ROOT / "configs" / "injection_v1.1.json").read_text())

    # Redraw host selection with new base seed
    cfg_host_draw = dict(cfg_inj, seed=base_seed)
    quiet_hosts_df = inject.select_hosts(sector, cfg_host_draw, con)
    quiet_hosts = set(quiet_hosts_df.tic.astype(int))

    # Read the 25 Gaia DR3 stars
    gaia_path = settings.ROOT / gaia_csv
    if gaia_path.exists():
        df_gaia = pd.read_csv(gaia_path)
        gaia_tics = set(df_gaia["tic"].dropna().astype(int))
    else:
        gaia_tics = set()

    all_tics = quiet_hosts.union(gaia_tics)
    jobs = run.build([sector], weights, out, production=False, tics_filter=all_tics, tasks_override=[], write=False)

    from .. import holdout

    keep = []
    for j in jobs:
        tic = int(j["tic"])
        # Blinding: keep holdout stars OUT of every arm
        if holdout.is_holdout(tic):
            continue

        rm = run.host_rm(j)
        if rm:
            j.update(rstar=rm[0], mstar=rm[1])

        tasks = []

        if tic in quiet_hosts and rm:
            # Quiet star arms: injections with prewhitening and without across runs 12-15
            for r in injection_runs:
                seed = seed_for_pw2(base_seed, sector, r, tic)
                base = dict(mode="inject", run=r, seed=seed, wrong_density=False)
                tasks.append(dict(base, tag=f"quiet_pw_r{r}", prewhiten=True))
                tasks.append(dict(base, tag=f"quiet_nopw_r{r}", prewhiten=False))

        if tic in gaia_tics:
            # Real search with prewhitening (even without stellar parameters)
            tasks.append(dict(mode="real", tag="gaia_real_pw", prewhiten=True))
            # Injections with prewhitening (only if stellar parameters exist)
            if rm:
                for r in injection_runs:
                    seed = seed_for_pw2(base_seed, sector, r, tic)
                    base = dict(mode="inject", run=r, seed=seed, wrong_density=False)
                    tasks.append(dict(base, tag=f"gaia_inj_pw_r{r}", prewhiten=True))

        if tasks:
            j["tasks"] = tasks
            keep.append(j)

    cfg = settings.cfg()
    run.write_manifests(keep, out, cfg, cfg_inj, weights, name="validation_pw2")
    return keep


def run_pw2_report(
    dirs: List[Union[str, Path]],
    out_dir: Union[str, Path],
    sector: int = 70,
    gaia_csv: Path = GAIA_CSV_PATH,
    cfg: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """Generate validation report from completed task outputs."""
    if cfg is None:
        cfg = settings.cfg()

    out_p = Path(out_dir)
    out_p.mkdir(parents=True, exist_ok=True)

    df_rows = pd.concat([run.load_rows(d) for d in dirs], ignore_index=True) if dirs else pd.DataFrame()
    if "status" in df_rows.columns:
        df_rows = df_rows[df_rows.status == "ok"].copy()

    gaia_path = settings.ROOT / gaia_csv
    gaia_tics = set()
    if gaia_path.exists():
        df_gaia_meta = pd.read_csv(gaia_path)
        gaia_tics = set(df_gaia_meta["tic"].dropna().astype(int))

    # Pre-extract pulsation frequencies and classifications for Gaia stars
    gaia_class_map = {}
    man = None
    try:
        man = targets.fetch_manifest(sector).set_index("tic")
    except Exception:
        pass

    for tic in gaia_tics:
        try:
            fname = man.loc[tic].fname if (man is not None and tic in man.index) else None
            p = batch.lc_path(sector, fname, download=False) if fname else None
            if p and p.exists():
                cinfo = classify_gaia_star_from_file(p)
                gaia_class_map[tic] = cinfo
        except Exception:
            pass

    # Extract gaia_results (real search with pw)
    gaia_results = []
    df_real = df_rows[df_rows.task == "gaia_real_pw"] if "task" in df_rows.columns else pd.DataFrame()
    evaluated_tics = set()

    for r in df_real.itertuples():
        tic = int(r.tic)
        evaluated_tics.add(tic)
        cinfo = gaia_class_map.get(tic, {})
        grp = cinfo.get("group", getattr(r, "group", "unclassified"))
        f_dom = cinfo.get("f_dom", 0.0)
        puls_freqs = cinfo.get("pulsation_frequencies", [])
        period = getattr(r, "period", getattr(r, "P_grid", 0.0))
        if pd.isna(period):
            period = getattr(r, "P_grid", 0.0)
        
        # Determine unclean reason if any
        reason = None
        if not cinfo or f_dom <= 0:
            reason = "no dominant frequency"

        gaia_results.append({
            "tic": tic,
            "group": grp,
            "f_dom": f_dom,
            "period": float(period) if pd.notna(period) else 0.0,
            "pulsation_frequencies": puls_freqs,
            "unclean_reason": reason,
        })

    # Account for all 25 registered stars (F2 BLOCKER)
    # Any star not evaluated counts as NOT clean and is listed with its reason
    from .. import holdout
    stellar_params_map = run.star_params(sector, list(gaia_tics)) if gaia_tics else {}

    for tic in gaia_tics:
        if tic in evaluated_tics:
            continue
        cinfo = gaia_class_map.get(tic, {})
        grp = cinfo.get("group", "unclassified")
        f_dom = cinfo.get("f_dom", 0.0)
        puls_freqs = cinfo.get("pulsation_frequencies", [])

        # Diagnose reason
        fname = man.loc[tic].fname if (man is not None and tic in man.index) else None
        p = batch.lc_path(sector, fname, download=False) if fname else None
        has_lc = bool(p and p.exists())
        
        sp = stellar_params_map.get(tic, {})
        j = run.base_job(sector, tic, "dummy", None, None, sp)
        has_rm = bool(run.host_rm(j))

        if holdout.is_holdout(tic):
            reason = "holdout star"
        elif not has_lc:
            reason = "missing light curve"
        elif not has_rm:
            reason = "missing stellar parameters"
        elif f_dom <= 0:
            reason = "no dominant frequency"
        else:
            reason = "task error"

        gaia_results.append({
            "tic": tic,
            "group": grp,
            "f_dom": f_dom,
            "period": 0.0,
            "pulsation_frequencies": puls_freqs,
            "unclean_reason": reason,
        })

    # Extract gaia_injections
    gaia_injections = []
    df_g_inj = df_rows[df_rows.task.str.startswith("gaia_inj_pw_r")] if "task" in df_rows.columns else pd.DataFrame()
    for r in df_g_inj.itertuples():
        tic = int(r.tic)
        grp = gaia_class_map.get(tic, {}).get("group", "delta_scuti")
        exp_s = getattr(r, "snr_exp", getattr(r, "expected_S", getattr(r, "S", 0.0)))
        exp_s = float(exp_s) if pd.notna(exp_s) else 0.0
        rec = bool(getattr(r, "matched", False) and (getattr(r, "S", 0.0) >= S_REF))
        gaia_injections.append({
            "tic": tic,
            "group": grp,
            "expected_S": exp_s,
            "recovered": rec,
        })

    # Extract quiet_injections
    quiet_injections = []
    df_q_pw = df_rows[df_rows.task.str.startswith("quiet_pw_r")] if "task" in df_rows.columns else pd.DataFrame()
    for r in df_q_pw.itertuples():
        tic = int(r.tic)
        exp_s = getattr(r, "snr_exp", getattr(r, "expected_S", getattr(r, "S", 0.0)))
        exp_s = float(exp_s) if pd.notna(exp_s) else 0.0
        rec = bool(getattr(r, "matched", False) and (getattr(r, "S", 0.0) >= S_REF))
        quiet_injections.append({
            "tic": tic,
            "expected_S": exp_s,
            "recovered": rec,
            "prewhiten": True,
        })

    df_q_nopw = df_rows[df_rows.task.str.startswith("quiet_nopw_r")] if "task" in df_rows.columns else pd.DataFrame()
    for r in df_q_nopw.itertuples():
        tic = int(r.tic)
        exp_s = getattr(r, "snr_exp", getattr(r, "expected_S", getattr(r, "S", 0.0)))
        exp_s = float(exp_s) if pd.notna(exp_s) else 0.0
        rec = bool(getattr(r, "matched", False) and (getattr(r, "S", 0.0) >= S_REF))
        quiet_injections.append({
            "tic": tic,
            "expected_S": exp_s,
            "recovered": rec,
            "prewhiten": False,
        })

    summary = evaluate_pw2_report(gaia_results, gaia_injections, quiet_injections, cfg=cfg)

    # Output JSON
    json_path = out_p / "validation_pw2.json"
    with open(json_path, "w") as fp:
        json.dump(summary, fp, indent=2, default=float)

    # Output CSVs
    # 1. Real search on Gaia stars
    c1_records = summary["criterion_1_real_gaia_search"]["real_gaia_table"]
    df_c1 = pd.DataFrame(c1_records)
    if "pulsation_frequencies" in df_c1.columns:
        df_c1["pulsation_frequencies"] = df_c1["pulsation_frequencies"].apply(lambda freqs: ";".join(f"{x:.4f}" for x in freqs))
    df_c1.to_csv(out_p / "gaia_stars_pw_evaluation.csv", index=False)

    # 2. Injections summary table by bin
    crit2_bins = summary["criterion_2_gaia_vs_quiet_injections"]["bins"]
    pd.DataFrame(crit2_bins).to_csv(out_p / "gaia_vs_quiet_injections_bins.csv", index=False)

    return summary


def main():
    ap = argparse.ArgumentParser(description="v1.1 Prewhitening validation (T2)")
    sp = ap.add_subparsers(dest="cmd", required=True)
    b = sp.add_parser("build")
    b.add_argument("--weights", required=True)
    b.add_argument("--out", required=True)
    b.add_argument("--sector", type=int, default=70)

    r = sp.add_parser("report")
    r.add_argument("dirs", nargs="+")
    r.add_argument("--out", required=True)
    r.add_argument("--sector", type=int, default=70)
    a = ap.parse_args()

    if a.cmd == "build":
        weights = [(h, float(x)) for h, x in (p.split(":") for p in a.weights.split(","))]
        build_pw2_manifests(weights, a.out, sector=a.sector)
    else:
        run_pw2_report(a.dirs, a.out, sector=a.sector)


if __name__ == "__main__":
    main()
