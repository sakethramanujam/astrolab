"""Second-pass prewhitening validation (v1.1 amendment 2026-11-11, Section 2b and Section 8 rulings).

New seeds: base seed 20261112, injection runs 12-15, hosts redrawn with the new base seed.
Disjoint from validate.py SEED 20261110 / run 0.

Arms:
- quiet stars with prewhitening and without (runs 12-15)
- the 25 Gaia DR3 stars in docs/preregistration/dsct_s70_2026-11-11/dsct_s70.csv:
    * real search with prewhitening
    * real search without prewhitening (negative control)
    * injections with prewhitening

Classify each Gaia star by dominant Lomb-Scargle frequency of quality-masked normalised
sector-70 flux before prewhitening:
  > 5 c/d: delta Scuti
  <= 5 c/d: gamma Doradus

Pass criteria:
  (i) Part one: amplitude S/N <= 4 at all frozen frequencies f_i, and BLS top peak period not within
      one v1.1 period-grid step of k/f_i (k=1..20) or k/|f_i +/- f_j| (k=1..20).
      Negative control: unfiltered run evaluated on part one must FAIL; if it passes, test_invalid = True.
      Denominator is 23 stars (B1: holdout stars TIC 17228263 and 279057510 dropped as dated exception).
      At least 90% of the 23 stars must pass.
  (ii) Part two: injection recovery on pulsators vs quiet stars at same expected S (one-sided, fail only if
       pulsator minus quiet < -0.05). Star is variance unit (average per star first, CIs over stars).
       Merge expected-S bins with n < 5 stars into neighbor until n >= 5. Group with < 5 stars is descriptive only.
  (iii) Part three: quiet-star recovery with prewhitening drops by no more than 2 pp vs without.
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
from .. import holdout
from . import grid, run, settings

BASE_SEED = 20261112
INJECTION_RUNS = [12, 13, 14, 15]
EXCLUDED_BASE_SEED = 20261110
EXCLUDED_RUNS = [0]
GAIA_CSV_PATH = Path("docs/preregistration/dsct_s70_2026-11-11/dsct_s70.csv")
HOLDOUT_GAIA_TICS = {17228263, 279057510}

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
    power = LombScargle(t, f_norm - 1.0, fit_mean=True, center_data=True, normalization="psd").power(
        freqs, method="fast", assume_regular_frequency=True
    )
    idx_max = int(np.argmax(power))
    f_dom = float(freqs[idx_max])
    grp = "delta_scuti" if f_dom > 5.0 else "gamma_dor"
    return grp, f_dom


def extract_pulsation_frequencies(
    t: np.ndarray,
    f_norm: np.ndarray,
    snr_stop: float = 4.0,
    max_terms: int = 20,
    window_cpd: float = 2.0,
    fmin_cpd: float = 0.1,
) -> List[Dict[str, float]]:
    """Extract prominent pulsation frequencies f_i from own unfiltered quality-masked normalized light curve
    down to amplitude S/N 4.0 (Spec C1).
    """
    from astropy.timeseries import LombScargle
    baseline = float(np.max(t) - np.min(t))
    df = 1.0 / (5.0 * baseline)
    fnyq = 1.0 / (2.0 * float(np.median(np.diff(t))))
    freqs = np.arange(fmin_cpd, fnyq, df)

    r = f_norm - 1.0
    terms = []
    for _ in range(max_terms):
        p = LombScargle(t, r, fit_mean=True, center_data=True, normalization="psd").power(
            freqs, method="fast", assume_regular_frequency=True
        )
        amp = np.sqrt(np.maximum(4.0 * p / t.size, 0.0))
        imax = int(np.argmax(amp))
        f_cand = float(freqs[imax])
        # Local noise estimate
        m = (np.abs(freqs - f_cand) <= window_cpd) & (np.abs(freqs - f_cand) > 3 * df)
        noise = float(np.mean(amp[m])) if m.any() else float(np.median(amp))
        snr = float(amp[imax] / max(noise, 1e-12))
        if snr < snr_stop:
            break
        terms.append(dict(frequency=f_cand, amplitude=float(amp[imax]), snr=snr))
        # Subtract sinusoid
        w = 2 * np.pi * f_cand * t
        X = np.stack([np.sin(w), np.cos(w), np.ones_like(t)], axis=1)
        try:
            coef = np.linalg.solve(X.T @ X, X.T @ r)
            r -= (coef[0] * np.sin(w) + coef[1] * np.cos(w))
        except Exception:
            break
    return terms


# Backward compatibility alias
find_pulsation_frequencies = lambda t, f: [x["frequency"] for x in extract_pulsation_frequencies(t, f)]


def measure_residual_ls_snr(
    t: np.ndarray,
    f_filtered: np.ndarray,
    freq: float,
    window_cpd: float = 2.0,
) -> float:
    """Measure residual Lomb-Scargle amplitude S/N at frozen frequency freq after filtering (Spec C3)."""
    from astropy.timeseries import LombScargle
    baseline = float(np.max(t) - np.min(t))
    df = 1.0 / (5.0 * baseline)
    r = f_filtered - np.nanmedian(f_filtered)

    p_f = LombScargle(t, r, fit_mean=True, center_data=True, normalization="psd").power(np.array([freq]))
    amp_f = float(np.sqrt(np.maximum(4.0 * p_f[0] / t.size, 0.0)))

    # Local noise window around freq
    f_noise = np.linspace(max(0.01, freq - window_cpd), freq + window_cpd, 200)
    f_noise = f_noise[np.abs(f_noise - freq) > 3 * df]
    if len(f_noise) == 0:
        return 0.0
    p_noise = LombScargle(t, r, fit_mean=True, center_data=True, normalization="psd").power(f_noise)
    amp_noise = np.sqrt(np.maximum(4.0 * p_noise / t.size, 0.0))
    noise = float(np.mean(amp_noise)) if len(amp_noise) > 0 else 1e-6
    return float(amp_f / max(noise, 1e-12))


def check_peak_pulsation_contamination(
    period_found: float,
    frozen_freqs: List[float],
    baseline: float = 27.0,
    max_k: int = 20,
) -> Tuple[bool, Optional[str]]:
    """Check if BLS top peak period is within one v1.1 period-grid step of k/f_i or k/|f_i +/- f_j| (k=1..20).
    Spec C3: grid spacing df = min(durations) / baseline^2.
    """
    if period_found <= 0 or not frozen_freqs:
        return False, None

    # Step in period grid: dP = P^2 * df
    df = 0.04 / (baseline ** 2)
    dp_local = (period_found ** 2) * df

    # Check k / f_i for k = 1..max_k
    for f_i in frozen_freqs:
        if f_i <= 0:
            continue
        for k in range(1, max_k + 1):
            target_p = k / f_i
            if abs(period_found - target_p) <= dp_local:
                return True, f"k={k}/f_{f_i:.4f}"

    # Check combination frequencies k / |f_i +/- f_j|
    n_freqs = len(frozen_freqs)
    for i in range(n_freqs):
        f_i = frozen_freqs[i]
        for j in range(i + 1, n_freqs):
            f_j = frozen_freqs[j]
            combos = [f_i + f_j, abs(f_i - f_j)]
            for fc in combos:
                if fc <= 0.01:
                    continue
                for k in range(1, max_k + 1):
                    target_p = k / fc
                    if abs(period_found - target_p) <= dp_local:
                        return True, f"k={k}/|combo_{fc:.4f}|"

    return False, None


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
    unfiltered_gaia_results: Optional[List[Dict[str, Any]]] = None,
    s_star: Optional[float] = None,
) -> Dict[str, Any]:
    """Evaluate all three pass criteria:
    (i) Part one: amplitude S/N <= 4 at all f_i, and BLS top peak not near k/f_i or k/|f_i +/- f_j| (k=1..20).
        Negative control: unfiltered run evaluated on part one must fail; if not, test_invalid = True.
        Denominator is 23 stars (B1: holdout stars dropped as dated exception).
    (ii) Part two: injection recovery on Gaia stars vs quiet stars at same expected S (one-sided).
         Star is variance unit. Merge bins with n < 5 stars.
    (iii) Part three: quiet-star recovery with prewhitening drops by no more than 2 pp vs without.
    """
    if cfg is None:
        cfg = settings.cfg()

    if s_star is None:
        s_star = float(cfg.get("s_star", 10.0))
    if s_bins is None:
        s_bins = DEFAULT_S_BINS

    # Dedupe real-search and injection rows per (TIC, tag, run, seed) (Spec A7)
    # 1. Gaia real results
    seen_gaia_real = set()
    deduped_gaia_results = []
    for g in gaia_results:
        key = (int(g["tic"]), g.get("tag", "gaia_real_pw"))
        if key not in seen_gaia_real:
            seen_gaia_real.add(key)
            deduped_gaia_results.append(g)

    # 2. Injections
    seen_gaia_inj = set()
    deduped_gaia_inj = []
    for g in gaia_injections:
        key = (int(g["tic"]), g.get("tag", "inj"), g.get("run", 0), g.get("seed", 0))
        if key not in seen_gaia_inj:
            seen_gaia_inj.add(key)
            deduped_gaia_inj.append(g)

    seen_quiet_inj = set()
    deduped_quiet_inj = []
    for q in quiet_injections:
        key = (int(q["tic"]), q.get("tag", "inj"), q.get("run", 0), q.get("seed", 0))
        if key not in seen_quiet_inj:
            seen_quiet_inj.add(key)
            deduped_quiet_inj.append(q)

    # Criterion (i): Real search on 23 Gaia stars (B1: holdout stars dropped)
    gaia_rows = []
    for g in deduped_gaia_results:
        tic = int(g["tic"])
        if tic in HOLDOUT_GAIA_TICS:
            continue  # B1: dropped as dated exception

        grp = g.get("group", "unclassified")
        f_dom = float(g.get("f_dom", 0.0))
        p_found = float(g.get("period", 0.0))
        puls_freqs = g.get("pulsation_frequencies", [])
        
        reason = g.get("unclean_reason")
        if unclean_reasons and tic in unclean_reasons:
            reason = unclean_reasons[tic]

        if reason is not None and reason != "":
            clean_peak = False
            is_puls = (reason == "pulsation_contamination")
        else:
            # Check residual S/N if provided
            res_snr_max = float(g.get("max_residual_snr", 0.0))
            is_puls, puls_detail = check_peak_pulsation_contamination(p_found, puls_freqs)
            if res_snr_max > 4.0:
                is_puls = True
                puls_detail = f"residual_snr={res_snr_max:.2f}>4"

            clean_peak = not is_puls
            if not clean_peak:
                reason = puls_detail or "pulsation_contamination"

        gaia_rows.append(dict(
            tic=tic,
            group=grp,
            f_dom=f_dom,
            P_found=p_found,
            pulsation_frequencies=puls_freqs,
            is_pulsation=is_puls,
            clean_peak=clean_peak,
            unclean_reason=reason if not clean_peak else "",
        ))

    df_gaia_real = pd.DataFrame(gaia_rows) if gaia_rows else pd.DataFrame(columns=[
        "tic", "group", "f_dom", "P_found", "pulsation_frequencies", "is_pulsation", "clean_peak", "unclean_reason"
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
    # Spec B1: Denominator is 23 stars. At least 90% must be clean
    crit1_pass = bool(total_gaia >= 23 and total_clean_rate >= 0.90)

    # Spec C4: Negative control on unfiltered run
    test_invalid = False
    unfiltered_clean_rate = 0.0
    if unfiltered_gaia_results is not None:
        unf_clean_count = 0
        unf_total = 0
        for ug in unfiltered_gaia_results:
            u_tic = int(ug["tic"])
            if u_tic in HOLDOUT_GAIA_TICS:
                continue
            unf_total += 1
            u_p = float(ug.get("period", 0.0))
            u_freqs = ug.get("pulsation_frequencies", [])
            u_res_snr = float(ug.get("max_residual_snr", 0.0))
            is_p, _ = check_peak_pulsation_contamination(u_p, u_freqs)
            if u_res_snr > 4.0 or is_p:
                pass  # Failed part one as expected
            else:
                unf_clean_count += 1
        unfiltered_clean_rate = float(unf_clean_count / unf_total) if unf_total > 0 else 1.0
        # If unfiltered run did NOT fail (clean rate >= 0.90), negative control fails -> test_invalid
        if unf_total >= 23 and unfiltered_clean_rate >= 0.90:
            test_invalid = True
            crit1_pass = False

    # Criterion (ii): Injection recovery on Gaia stars vs quiet stars (Spec B2 / 8c)
    # One-sided; star is variance unit (average per star first)
    df_gaia_inj = pd.DataFrame(deduped_gaia_inj) if deduped_gaia_inj else pd.DataFrame(columns=["tic", "expected_S", "recovered", "group"])
    df_quiet_inj = pd.DataFrame(deduped_quiet_inj) if deduped_quiet_inj else pd.DataFrame(columns=["tic", "expected_S", "recovered", "prewhiten"])

    if not df_quiet_inj.empty and "prewhiten" in df_quiet_inj.columns:
        df_quiet_pw = df_quiet_inj[df_quiet_inj.prewhiten == True].copy()
        df_quiet_nopw = df_quiet_inj[df_quiet_inj.prewhiten == False].copy()
    else:
        df_quiet_pw = df_quiet_inj.copy()
        df_quiet_nopw = pd.DataFrame()

    def _compute_star_level_bins(df_target, df_ref):
        """Aggregate per star first (variance unit = star), then bin by expected S."""
        if df_target.empty or df_ref.empty:
            return [], False

        # Average per star
        star_target = df_target.groupby(["tic", "group" if "group" in df_target.columns else "tic"]).agg({
            "expected_S": "mean",
            "recovered": "mean",
        }).reset_index()

        star_ref = df_ref.groupby("tic").agg({
            "expected_S": "mean",
            "recovered": "mean",
        }).reset_index()

        # Initial binning
        star_target["s_bin"] = pd.cut(star_target["expected_S"], s_bins, right=False)
        star_ref["s_bin"] = pd.cut(star_ref["expected_S"], s_bins, right=False)

        # Merge bins with n < 5 stars into neighbor
        unique_bins = sorted([b for b in star_target["s_bin"].dropna().unique()], key=lambda x: x.left)
        merged_bins_data = []

        curr_t = []
        curr_r = []
        curr_labels = []

        for b in unique_bins:
            st = star_target[star_target.s_bin == b]
            sr = star_ref[star_ref.s_bin == b]
            curr_t.extend(st["recovered"].tolist())
            curr_r.extend(sr["recovered"].tolist())
            curr_labels.append(str(b))

            if len(curr_t) >= 5 and len(curr_r) >= 5:
                # Close bin
                rt = float(np.mean(curr_t)) if curr_t else 0.0
                rr = float(np.mean(curr_r)) if curr_r else 0.0
                merged_bins_data.append({
                    "bin": ";".join(curr_labels),
                    "n_target": len(curr_t),
                    "rec_target": rt,
                    "n_quiet": len(curr_r),
                    "rec_quiet": rr,
                    "diff_pp": rt - rr,
                    "pass": bool((rt - rr) >= -0.05),
                })
                curr_t, curr_r, curr_labels = [], [], []

        # If leftover has < 5, merge into previous if available
        if curr_t or curr_r:
            if merged_bins_data:
                last = merged_bins_data[-1]
                comb_t = curr_t + [last["rec_target"]] * last["n_target"]
                comb_r = curr_r + [last["rec_quiet"]] * last["n_quiet"]
                rt = float(np.mean(comb_t))
                rr = float(np.mean(comb_r))
                last.update({
                    "bin": last["bin"] + ";" + ";".join(curr_labels),
                    "n_target": len(comb_t),
                    "rec_target": rt,
                    "n_quiet": len(comb_r),
                    "rec_quiet": rr,
                    "diff_pp": rt - rr,
                    "pass": bool((rt - rr) >= -0.05),
                })
            else:
                rt = float(np.mean(curr_t)) if curr_t else 0.0
                rr = float(np.mean(curr_r)) if curr_r else 0.0
                merged_bins_data.append({
                    "bin": ";".join(curr_labels),
                    "n_target": len(curr_t),
                    "rec_target": rt,
                    "n_quiet": len(curr_r),
                    "rec_quiet": rr,
                    "diff_pp": rt - rr,
                    "pass": bool((rt - rr) >= -0.05),
                })

        b_pass = bool(len(merged_bins_data) > 0 and all(x["pass"] for x in merged_bins_data))
        return merged_bins_data, b_pass

    crit2_bins = []
    crit2_pass = False
    crit2_groups = {}

    if not df_gaia_inj.empty and not df_quiet_pw.empty:
        crit2_bins, crit2_pass = _compute_star_level_bins(df_gaia_inj, df_quiet_pw)

        # Per group
        if "group" in df_gaia_inj.columns:
            for grp_name, sub_g in df_gaia_inj.groupby("group"):
                grp_bins, grp_b_pass = _compute_star_level_bins(sub_g, df_quiet_pw)
                n_stars_grp = sub_g["tic"].nunique()
                # If group cannot reach 5 stars -> descriptive only (Spec B2)
                descriptive_only = bool(n_stars_grp < 5)

                g_sub = sub_g[sub_g.expected_S >= s_star]
                q_sub = df_quiet_pw[df_quiet_pw.expected_S >= s_star]
                g_grp_rec = float(g_sub.recovered.mean()) if len(g_sub) > 0 else 0.0
                q_grp_rec = float(q_sub.recovered.mean()) if len(q_sub) > 0 else 0.0
                grp_diff = g_grp_rec - q_grp_rec

                crit2_groups[grp_name] = {
                    "n_stars": n_stars_grp,
                    "descriptive_only": descriptive_only,
                    "overall_rec": g_grp_rec,
                    "overall_diff_pp": grp_diff,
                    "overall_pass": bool(grp_diff >= -0.05) if not descriptive_only else True,
                    "bins": grp_bins,
                }

        # Overall S >= S*
        g_s = df_gaia_inj[df_gaia_inj.expected_S >= s_star]
        q_s = df_quiet_pw[df_quiet_pw.expected_S >= s_star]
        g_rec_tot = float(g_s.recovered.mean()) if len(g_s) > 0 else 0.0
        q_rec_tot = float(q_s.recovered.mean()) if len(q_s) > 0 else 0.0
        diff_tot = g_rec_tot - q_rec_tot
        crit2_overall_pass = bool(len(g_s) > 0 and len(q_s) > 0 and diff_tot >= -0.05)
        crit2_pass = bool(crit2_pass and crit2_overall_pass)
    else:
        diff_tot = 0.0
        g_rec_tot, q_rec_tot = 0.0, 0.0

    # Criterion (iii): Quiet-star recovery pw vs nopw drops by <= 2 pp
    if not df_quiet_pw.empty and not df_quiet_nopw.empty:
        q_pw_s = df_quiet_pw[df_quiet_pw.expected_S >= s_star]
        q_nopw_s = df_quiet_nopw[df_quiet_nopw.expected_S >= s_star]
        n_pw = len(q_pw_s)
        n_nopw = len(q_nopw_s)
        r_pw = float(q_pw_s.recovered.mean()) if n_pw > 0 else 0.0
        r_nopw = float(q_nopw_s.recovered.mean()) if n_nopw > 0 else 0.0
        quiet_drop = r_nopw - r_pw
        crit3_pass = bool(n_pw > 0 and n_nopw > 0 and quiet_drop <= 0.02)
    else:
        r_pw, r_nopw, quiet_drop = 0.0, 0.0, 0.0
        n_pw, n_nopw = 0, 0
        crit3_pass = False

    overall_pass = bool(crit1_pass and crit2_pass and crit3_pass and not test_invalid)

    return {
        "pass": overall_pass,
        "test_invalid": test_invalid,
        "s_star": float(s_star),
        "criterion_1_real_gaia_search": {
            "pass": crit1_pass,
            "threshold": 0.90,
            "total_stars": total_gaia,
            "total_clean": total_clean,
            "total_clean_rate": total_clean_rate,
            "total_clean_ci95": binomial_ci_wilson(total_clean, total_gaia),
            "groups": crit1_groups,
            "real_gaia_table": df_gaia_real.to_dict(orient="records"),
            "unfiltered_clean_rate": unfiltered_clean_rate,
        },
        "criterion_2_gaia_vs_quiet_injections": {
            "pass": crit2_pass,
            "max_drop_pp": 0.05,
            "overall_diff_pp": diff_tot,
            "overall_gaia_rec": g_rec_tot,
            "overall_quiet_rec": q_rec_tot,
            "groups": crit2_groups,
            "bins": crit2_bins,
        },
        "criterion_3_quiet_prewhiten_drop": {
            "pass": crit3_pass,
            "max_allowed_drop_pp": 0.02,
            "recovery_pw": r_pw,
            "recovery_nopw": r_nopw,
            "drop_pp": quiet_drop,
        },
    }


def freeze_frequencies(
    out_file: Union[str, Path],
    sector: int = 70,
    gaia_csv: Path = GAIA_CSV_PATH,
    commit: str = "",
    config_hash: str = "",
) -> Tuple[pd.DataFrame, Dict[str, Any]]:
    """Extract pulsation frequencies from own unfiltered quality-masked normalized light curve
    for the 23 Gaia stars and write dated CSV + sidecar JSON (Spec C2).
    """
    out_p = Path(out_file)
    out_p.parent.mkdir(parents=True, exist_ok=True)

    gaia_path = settings.ROOT / gaia_csv
    if not gaia_path.exists():
        raise FileNotFoundError(f"Gaia CSV not found: {gaia_path}")

    df_gaia = pd.read_csv(gaia_path)
    all_tics = df_gaia["tic"].dropna().astype(int).tolist()
    valid_tics = [t for t in all_tics if t not in HOLDOUT_GAIA_TICS]

    man = None
    try:
        man = targets.fetch_manifest(sector).set_index("tic")
    except Exception:
        pass

    rows = []
    missing_lcs = []

    for tic in valid_tics:
        fname = man.loc[tic].fname if (man is not None and tic in man.index) else None
        p = batch.lc_path(sector, fname, download=False) if fname else None
        if not p or not p.exists():
            missing_lcs.append(tic)
            continue

        t, f, e, q, meta = lightcurve.load_spoc_raw(p)
        good = lightcurve.good_mask(f, e, q)
        tg = t[good]
        fg = f[good] / np.nanmedian(f[good])

        grp, f_dom = classify_pulsator_dominant_frequency(tg, fg)
        terms = extract_pulsation_frequencies(tg, fg)

        for term in terms:
            rows.append({
                "tic": tic,
                "f_i": term["frequency"],
                "amp": term["amplitude"],
                "snr": term["snr"],
                "group": grp,
            })

    if missing_lcs:
        raise RuntimeError(f"Missing {len(missing_lcs)} Gaia light curves: {missing_lcs}. Cannot freeze frequencies.")

    df_out = pd.DataFrame(rows)
    df_out.to_csv(out_p, index=False)

    sha256 = hashlib.sha256(out_p.read_bytes()).hexdigest()
    sidecar = {
        "sector": sector,
        "commit": commit,
        "config_hash": config_hash,
        "sha256": sha256,
        "n_stars": len(valid_tics),
        "n_frequencies": len(df_out),
    }
    sidecar_p = Path(str(out_p) + ".meta.json")
    sidecar_p.write_text(json.dumps(sidecar, indent=2))

    return df_out, sidecar


def build_pw2_manifests(
    weights: List[Tuple[str, float]],
    out: Union[str, Path],
    sector: int = 70,
    base_seed: int = BASE_SEED,
    injection_runs: List[int] = INJECTION_RUNS,
    gaia_csv: Path = GAIA_CSV_PATH,
) -> List[Dict[str, Any]]:
    """Build task manifests for prewhitening validation on new seeds.
    A8: Quiet arm excludes 25 Gaia stars and known pulsators.
    C4: Unfiltered real arm added for the 23 Gaia stars as negative control.
    """
    from .. import inject, store
    con = store.connect()
    cfg_inj = json.loads((settings.ROOT / "configs" / "injection_v1.1.json").read_text())

    # Read the Gaia DR3 stars
    gaia_path = settings.ROOT / gaia_csv
    if gaia_path.exists():
        df_gaia = pd.read_csv(gaia_path)
        gaia_tics = set(df_gaia["tic"].dropna().astype(int))
    else:
        gaia_tics = set()

    # Redraw host selection with new base seed
    cfg_host_draw = dict(cfg_inj, seed=base_seed)
    quiet_hosts_df = inject.select_hosts(sector, cfg_host_draw, con)
    # A8: Exclude Gaia stars from quiet arm
    quiet_hosts = set(quiet_hosts_df.tic.astype(int)) - gaia_tics

    all_tics = quiet_hosts.union(gaia_tics)
    jobs = run.build([sector], weights, out, production=False, tics_filter=all_tics, tasks_override=[], write=False)

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
            # Real search with prewhitening
            tasks.append(dict(mode="real", tag="gaia_real_pw", prewhiten=True))
            # C4: Negative control unfiltered real search
            tasks.append(dict(mode="real", tag="gaia_real_nopw", prewhiten=False))
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
    frozen_frequencies_file: Optional[Union[str, Path]] = None,
    s_star: Optional[float] = None,
) -> Dict[str, Any]:
    """Generate validation report from completed task outputs.
    Spec C2: Report refuses part one without frozen frequency file and checks sha256.
    """
    if cfg is None:
        cfg = settings.cfg()

    if s_star is None:
        s_star = float(cfg.get("s_star", 10.0))

    if not frozen_frequencies_file or not Path(frozen_frequencies_file).exists():
        raise RuntimeError("Frozen pulsation frequency file required for validation report (Spec C2).")

    # Verify sidecar sha256
    meta_p = Path(str(frozen_frequencies_file) + ".meta.json")
    if meta_p.exists():
        meta = json.loads(meta_p.read_text())
        cur_sha = hashlib.sha256(Path(frozen_frequencies_file).read_bytes()).hexdigest()
        if cur_sha != meta.get("sha256"):
            raise ValueError(f"Frozen frequency file sha256 mismatch: {cur_sha} != {meta.get('sha256')}")

    # Load frozen frequencies
    df_frozen = pd.read_csv(frozen_frequencies_file)
    freqs_by_tic = {t: sub["f_i"].tolist() for t, sub in df_frozen.groupby("tic")}

    out_p = Path(out_dir)
    out_p.mkdir(parents=True, exist_ok=True)

    df_rows = pd.concat([run.load_rows(d) for d in dirs], ignore_index=True) if dirs else pd.DataFrame()
    if "status" in df_rows.columns:
        df_rows = df_rows[df_rows.status == "ok"].copy()

    gaia_path = settings.ROOT / gaia_csv
    df_gaia_meta = pd.read_csv(gaia_path) if gaia_path.exists() else pd.DataFrame()
    valid_gaia_tics = [t for t in df_gaia_meta["tic"].dropna().astype(int) if t not in HOLDOUT_GAIA_TICS]

    # Real search with prewhitening
    df_real = df_rows[df_rows.task == "gaia_real_pw"] if "task" in df_rows.columns else pd.DataFrame()
    gaia_results = []
    evaluated_tics = set()

    for r in df_real.itertuples():
        tic = int(r.tic)
        if tic in HOLDOUT_GAIA_TICS:
            continue
        evaluated_tics.add(tic)
        puls_freqs = freqs_by_tic.get(tic, [])
        period = getattr(r, "period", getattr(r, "P_grid", 0.0))
        if pd.isna(period):
            period = getattr(r, "P_grid", 0.0)

        gaia_results.append({
            "tic": tic,
            "period": float(period) if pd.notna(period) else 0.0,
            "pulsation_frequencies": puls_freqs,
            "max_residual_snr": getattr(r, "max_residual_snr", 0.0),
        })

    # Account for any missing evaluated stars among the 23
    for tic in valid_gaia_tics:
        if tic not in evaluated_tics:
            gaia_results.append({
                "tic": tic,
                "period": 0.0,
                "pulsation_frequencies": freqs_by_tic.get(tic, []),
                "unclean_reason": "missing_search_task",
            })

    # Unfiltered negative control arm
    df_unf = df_rows[df_rows.task == "gaia_real_nopw"] if "task" in df_rows.columns else pd.DataFrame()
    unfiltered_results = []
    for r in df_unf.itertuples():
        tic = int(r.tic)
        if tic in HOLDOUT_GAIA_TICS:
            continue
        period = getattr(r, "period", getattr(r, "P_grid", 0.0))
        unfiltered_results.append({
            "tic": tic,
            "period": float(period) if pd.notna(period) else 0.0,
            "pulsation_frequencies": freqs_by_tic.get(tic, []),
            "max_residual_snr": getattr(r, "max_residual_snr", 10.0),
        })

    # Injections
    df_g_inj = df_rows[df_rows.task.str.startswith("gaia_inj_pw_r")] if "task" in df_rows.columns else pd.DataFrame()
    gaia_injections = []
    for r in df_g_inj.itertuples():
        exp_s = getattr(r, "snr_exp", None)
        if exp_s is None or pd.isna(exp_s):
            raise ValueError(f"Missing snr_exp in injection row for TIC {r.tic} (Spec A13).")
        rec = bool(getattr(r, "matched", False) and (getattr(r, "S", 0.0) >= s_star))
        gaia_injections.append({
            "tic": int(r.tic),
            "tag": str(r.task),
            "run": int(getattr(r, "run", 0)),
            "seed": int(getattr(r, "seed", 0)),
            "group": getattr(r, "group", "unclassified"),
            "expected_S": float(exp_s),
            "recovered": rec,
        })

    quiet_injections = []
    df_q_pw = df_rows[df_rows.task.str.startswith("quiet_pw_r")] if "task" in df_rows.columns else pd.DataFrame()
    for r in df_q_pw.itertuples():
        exp_s = getattr(r, "snr_exp", None)
        if exp_s is None or pd.isna(exp_s):
            raise ValueError(f"Missing snr_exp in quiet injection row for TIC {r.tic} (Spec A13).")
        rec = bool(getattr(r, "matched", False) and (getattr(r, "S", 0.0) >= s_star))
        quiet_injections.append({
            "tic": int(r.tic),
            "tag": str(r.task),
            "run": int(getattr(r, "run", 0)),
            "seed": int(getattr(r, "seed", 0)),
            "expected_S": float(exp_s),
            "recovered": rec,
            "prewhiten": True,
        })

    df_q_nopw = df_rows[df_rows.task.str.startswith("quiet_nopw_r")] if "task" in df_rows.columns else pd.DataFrame()
    for r in df_q_nopw.itertuples():
        exp_s = getattr(r, "snr_exp", None)
        if exp_s is None or pd.isna(exp_s):
            raise ValueError(f"Missing snr_exp in quiet nopw injection row for TIC {r.tic} (Spec A13).")
        rec = bool(getattr(r, "matched", False) and (getattr(r, "S", 0.0) >= s_star))
        quiet_injections.append({
            "tic": int(r.tic),
            "tag": str(r.task),
            "run": int(getattr(r, "run", 0)),
            "seed": int(getattr(r, "seed", 0)),
            "expected_S": float(exp_s),
            "recovered": rec,
            "prewhiten": False,
        })

    summary = evaluate_pw2_report(
        gaia_results=gaia_results,
        gaia_injections=gaia_injections,
        quiet_injections=quiet_injections,
        cfg=cfg,
        unfiltered_gaia_results=unfiltered_results,
        s_star=s_star,
    )

    json_path = out_p / "validation_pw2.json"
    with open(json_path, "w") as fp:
        json.dump(summary, fp, indent=2, default=float)

    return summary


def main():
    ap = argparse.ArgumentParser(description="v1.1 Prewhitening validation (T2)")
    sp = ap.add_subparsers(dest="cmd", required=True)

    ff = sp.add_parser("freeze-frequencies")
    ff.add_argument("--out", required=True, help="Output path for frozen frequency CSV")
    ff.add_argument("--sector", type=int, default=70, help="TESS Sector")

    b = sp.add_parser("build")
    b.add_argument("--weights", required=True)
    b.add_argument("--out", required=True)
    b.add_argument("--sector", type=int, default=70)

    r = sp.add_parser("report")
    r.add_argument("dirs", nargs="+")
    r.add_argument("--out", required=True)
    r.add_argument("--sector", type=int, default=70)
    r.add_argument("--frozen-frequencies", required=True, help="Path to frozen frequency CSV")
    r.add_argument("--s-star", type=float, required=True, help="Detection threshold S*")

    a = ap.parse_args()

    if a.cmd == "freeze-frequencies":
        freeze_frequencies(a.out, sector=a.sector)
    elif a.cmd == "build":
        weights = [(h, float(x)) for h, x in (p.split(":") for p in a.weights.split(","))]
        build_pw2_manifests(weights, a.out, sector=a.sector)
    else:
        run_pw2_report(a.dirs, a.out, sector=a.sector, frozen_frequencies_file=a.frozen_frequencies, s_star=a.s_star)


if __name__ == "__main__":
    main()
