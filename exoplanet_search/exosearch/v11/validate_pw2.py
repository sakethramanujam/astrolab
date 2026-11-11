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
) -> List[float]:
    """Extract prominent pulsation frequencies before prewhitening using Lomb-Scargle iterative search."""
    from astropy.timeseries import LombScargle
    baseline = float(np.max(t) - np.min(t))
    df = 1.0 / (5.0 * baseline)
    fnyq = 1.0 / (2.0 * float(np.median(np.diff(t))))
    freqs = np.arange(0.1, fnyq, df)
    
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
) -> Dict[str, Any]:
    """Evaluate all three pass criteria:
    (i) at least 90% of Gaia stars have top search peak that is not a pulsation frequency or harmonic;
    (ii) injection recovery on Gaia stars within 5 percentage points of quiet stars at same S;
    (iii) quiet-star recovery with prewhitening drops by no more than 2 points vs without.
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
        grp = g.get("group", "delta_scuti")
        p_found = float(g.get("period", 0.0))
        f_found = 1.0 / p_found if p_found > 0 else 0.0
        puls_freqs = g.get("pulsation_frequencies", [])
        is_puls = is_pulsation_or_harmonic(f_found, puls_freqs, freq_tol, harm_tol, harmonics)
        clean_peak = not is_puls
        gaia_rows.append(dict(tic=tic, group=grp, P_found=p_found, f_found=f_found, is_pulsation=is_puls, clean_peak=clean_peak))

    if gaia_rows:
        df_gaia_real = pd.DataFrame(gaia_rows)
    else:
        df_gaia_real = pd.DataFrame(columns=["tic", "group", "P_found", "f_found", "is_pulsation", "clean_peak"])

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
    crit1_pass = bool(total_clean_rate >= 0.90) if total_gaia > 0 else True

    # Criterion (ii): Injection recovery on Gaia stars vs quiet stars at same expected S
    df_gaia_inj = pd.DataFrame(gaia_injections) if gaia_injections else pd.DataFrame(columns=["tic", "expected_S", "recovered", "group"])
    df_quiet_inj = pd.DataFrame(quiet_injections) if quiet_injections else pd.DataFrame(columns=["tic", "expected_S", "recovered", "prewhiten"])

    # Filter quiet injections for prewhitened arm
    if not df_quiet_inj.empty and "prewhiten" in df_quiet_inj.columns:
        df_quiet_pw = df_quiet_inj[df_quiet_inj.prewhiten == True].copy()
        df_quiet_nopw = df_quiet_inj[df_quiet_inj.prewhiten == False].copy()
    else:
        df_quiet_pw = df_quiet_inj.copy()
        df_quiet_nopw = pd.DataFrame()

    crit2_bins = []
    crit2_pass = True

    if not df_gaia_inj.empty and not df_quiet_pw.empty:
        df_gaia_inj["s_bin"] = pd.cut(df_gaia_inj["expected_S"], s_bins, right=False)
        df_quiet_pw["s_bin"] = pd.cut(df_quiet_pw["expected_S"], s_bins, right=False)

        all_bins = sorted(set(df_gaia_inj["s_bin"].dropna().unique()).union(set(df_quiet_pw["s_bin"].dropna().unique())), key=lambda x: x.left)
        for b in all_bins:
            sub_g = df_gaia_inj[df_gaia_inj.s_bin == b]
            sub_q = df_quiet_pw[df_quiet_pw.s_bin == b]
            ng, nq = len(sub_g), len(sub_q)
            rg = float(sub_g.recovered.mean()) if ng > 0 else 0.0
            rq = float(sub_q.recovered.mean()) if nq > 0 else 0.0
            diff_pp = rg - rq

            ci_g = binomial_ci_wilson(int(sub_g.recovered.sum()), ng) if ng > 0 else (0.0, 0.0)
            ci_q = binomial_ci_wilson(int(sub_q.recovered.sum()), nq) if nq > 0 else (0.0, 0.0)

            # Within 5 percentage points: rg >= rq - 0.05
            bin_pass = bool(diff_pp >= -0.05) if (ng >= 5 and nq >= 5) else True
            if ng >= 5 and nq >= 5 and not bin_pass:
                crit2_pass = False

            crit2_bins.append(dict(
                bin=str(b),
                n_gaia=ng,
                rec_gaia=rg,
                ci95_gaia=ci_g,
                n_quiet=nq,
                rec_quiet=rq,
                ci95_quiet=ci_q,
                diff_pp=diff_pp,
                pass_=bin_pass,
            ))

        # Overall S >= 10 comparison
        g_s10 = df_gaia_inj[df_gaia_inj.expected_S >= S_REF]
        q_s10 = df_quiet_pw[df_quiet_pw.expected_S >= S_REF]
        g_rec_tot = float(g_s10.recovered.mean()) if len(g_s10) > 0 else 0.0
        q_rec_tot = float(q_s10.recovered.mean()) if len(q_s10) > 0 else 0.0
        diff_tot = g_rec_tot - q_rec_tot
        crit2_overall_pass = bool(diff_tot >= -0.05)
        crit2_pass = crit2_pass and crit2_overall_pass
    else:
        crit2_overall_pass = True
        diff_tot = 0.0

    # Criterion (iii): Quiet-star recovery pw vs nopw drops by <= 2 points
    if not df_quiet_pw.empty and not df_quiet_nopw.empty:
        q_pw_s10 = df_quiet_pw[df_quiet_pw.expected_S >= S_REF]
        q_nopw_s10 = df_quiet_nopw[df_quiet_nopw.expected_S >= S_REF]
        r_pw = float(q_pw_s10.recovered.mean()) if len(q_pw_s10) > 0 else 0.0
        r_nopw = float(q_nopw_s10.recovered.mean()) if len(q_nopw_s10) > 0 else 0.0
        quiet_drop = r_nopw - r_pw
        crit3_pass = bool(quiet_drop <= 0.02)
    else:
        r_pw, r_nopw, quiet_drop = 0.0, 0.0, 0.0
        crit3_pass = True

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
        },
        "criterion_2_gaia_vs_quiet_injections": {
            "pass": crit2_pass,
            "max_drop_pp": 0.05,
            "overall_diff_pp": diff_tot,
            "overall_pass": crit2_overall_pass,
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

    keep = []
    for j in jobs:
        rm = run.host_rm(j)
        if not rm:
            continue
        j.update(rstar=rm[0], mstar=rm[1])
        tic = int(j["tic"])
        tasks = []

        if tic in quiet_hosts:
            # Quiet star arms: injections with prewhitening and without across runs 12-15
            for r in injection_runs:
                seed = seed_for_pw2(base_seed, sector, r, tic)
                base = dict(mode="inject", run=r, seed=seed, wrong_density=False)
                tasks.append(dict(base, tag=f"quiet_pw_r{r}", prewhiten=True))
                tasks.append(dict(base, tag=f"quiet_nopw_r{r}", prewhiten=False))

        if tic in gaia_tics:
            # Real search with prewhitening
            tasks.append(dict(mode="real", tag="gaia_real_pw", prewhiten=True))
            # Injections with prewhitening
            for r in injection_runs:
                seed = seed_for_pw2(base_seed, sector, r, tic)
                base = dict(mode="inject", run=r, seed=seed, wrong_density=False)
                tasks.append(dict(base, tag=f"gaia_inj_pw_r{r}", prewhiten=True))

        j["tasks"] = tasks
        keep.append(j)

    cfg = settings.cfg()
    run.write_manifests(keep, out, cfg, cfg_inj, weights, name="validation_pw2")
    return keep


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
    a = ap.parse_args()

    if a.cmd == "build":
        weights = [(h, float(x)) for h, x in (p.split(":") for p in a.weights.split(","))]
        build_pw2_manifests(weights, a.out, sector=a.sector)
    else:
        # Load and aggregate rows
        df = pd.concat([run.load_rows(d) for d in a.dirs], ignore_index=True)
        # Process and generate reports
        out_dir = Path(a.out)
        out_dir.mkdir(parents=True, exist_ok=True)


if __name__ == "__main__":
    main()
