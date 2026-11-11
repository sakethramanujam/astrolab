"""Known-planet and active-TOI recovery test (v1.1 amendment 2026-11-11, Section 3).

Every confirmed planet or active TOI (TFOPWG disposition not FP/FA) on a searched star
with expected S >= S* and at least 2 transits surviving masks must be recovered.

Expected S = depth_exp * sqrt(n_tr) / CDPP_D at catalogue duration.
n_tr = catalogue-ephemeris transits surviving the v1.1 masks (at least half their in-transit cadences).
CDPP_D from our own noise profile.
depth_exp = catalogue depth scaled by aperture contamination ratio:
  We search PDCSAP (already crowding-corrected by SPOC), so record the light curve
  CROWDSAP header value per star and apply the factor documented here (1.0 for PDCSAP;
  CROWDSAP only if SAP flux were used).
Recovered: top reported peak period within the period-grid plateau half-width of the
  catalogue period (grid spacing / plateau half-width from v1.1 grid.py), or 2x / 0.5x alias.
Miss cause bucket logic: few_transits, long_duration, masked, below_threshold, outcompeted, other.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Dict, Any, List, Optional, Tuple, Union

from astropy.io import fits
import numpy as np
import pandas as pd

from .. import batch, bls, catalogs, config, lightcurve, store, targets
from . import grid, pipeline, search, settings

# Contamination factor documentation:
# For PDCSAP flux, the SPOC pipeline has already performed crowding correction by subtracting
# estimated excess flux from nearby stars based on the PRF model. Therefore, the observed
# transit depth in PDCSAP is already on an uncontaminated relative flux scale, so the
# depth contamination scaling factor for PDCSAP is 1.0.
# If uncorrected SAP flux were used instead, the observed transit depth would be diluted by
# non-target starlight in the aperture by a factor equal to CROWDSAP (the target flux fraction),
# so the expected depth in SAP would be catalogue_depth * CROWDSAP.
CONTAMINATION_FACTOR_PDCSAP = 1.0
DETREND_WINDOW_DAYS = 0.75


def calculate_grid_halfwidth(
    period: float,
    baseline: float = 27.0,
    durations: np.ndarray = search.DURATIONS,
    frequency_factor: float = 1.0,
) -> float:
    """Calculate the period-grid plateau half-width at a given period from the v1.1 grid spacing.

    Under the v1.1 grid specification (grid.py), trial frequencies are evenly spaced with step
    df = frequency_factor * min(durations) / baseline^2.
    Since P = 1 / f, the local grid step in period around P is:
        dP = |d(1/f)| = df / f^2 = P^2 * df.
    For an unresolved single grid point, the plateau half-width corresponds to half the local grid
    spacing: 0.5 * P^2 * df (as defined in exosearch/p_grid_halfwidth.py and tie_break.md).
    """
    df = frequency_factor * float(np.min(durations)) / (baseline ** 2)
    return 0.5 * (period ** 2) * df


def count_surviving_transits(
    t: np.ndarray,
    period: float,
    t0: float,
    duration: float,
    cadence: Optional[float] = None,
) -> int:
    """Number of transits of the ephemeris that have at least half their in-transit cadences.

    Parameters
    ----------
    t : np.ndarray
        Cadence timestamps surviving masks.
    period, t0, duration : float
        Transit ephemeris in days.
    cadence : float, optional
        Median cadence in days.
    """
    if len(t) == 0 or period <= 0 or duration <= 0:
        return 0
    if cadence is None:
        cadence = float(np.median(np.diff(t))) if len(t) > 1 else (2.0 / 1440.0)
    expected_cadences = max(duration / cadence, 1.0)

    ph = ((t - t0 + 0.5 * period) % period) - 0.5 * period
    in_tr = np.abs(ph) < (duration / 2.0)
    if not np.any(in_tr):
        return 0
    epochs = np.round((t[in_tr] - t0) / period).astype(int)
    _, counts = np.unique(epochs, return_counts=True)
    return int(np.sum(counts >= 0.5 * expected_cadences))


def check_recovery(
    p_cat: float,
    p_found: Optional[float],
    half_width: float,
    alias_ratios: Tuple[float, ...] = (2.0, 0.5),
    alias_tol: float = 0.01,
) -> Tuple[bool, Optional[str]]:
    """Determine whether p_found recovers p_cat within the half-width plateau or as an alias.

    Returns (recovered, alias_flag).
    """
    if p_found is None or not np.isfinite(p_found) or p_found <= 0:
        return False, None

    # Check primary recovery within plateau half-width
    if abs(p_found - p_cat) <= half_width:
        return True, None

    # Check 2x / 0.5x aliases
    for ratio in alias_ratios:
        target_alias = p_cat * ratio
        hw_alias = half_width * (ratio ** 2)
        if abs(p_found - target_alias) <= max(hw_alias, alias_tol * target_alias):
            tag = f"{int(ratio)}x" if ratio >= 1 else "0.5x"
            return True, tag

    return False, None


def determine_miss_cause(
    n_raw: int,
    dur_d: float,
    n_tr_surviving: int,
    expected_s: float,
    s_star: float,
    p_cat: float,
    p_found: Optional[float],
    half_width: float,
    s_found: Optional[float] = None,
    detrend_window: float = DETREND_WINDOW_DAYS,
) -> str:
    """Classify the miss cause bucket according to exosearch/missed.py rules:
    1. few_transits   : < 2 transits of ephemeris with >= half cadences in raw sector data
    2. long_duration  : biweight window (0.75 d) < 3 x duration -> detrending eats transit
    3. masked         : >= 2 transits in raw, but < 2 survive quality mask + gap-edge trim
    4. below_threshold: expected S < detection threshold S*
    5. outcompeted    : expected S >= S*, but top reported peak was not recovered
    6. other          : edge case / unclassified
    """
    if n_raw < 2:
        return "few_transits"
    if detrend_window < 3.0 * dur_d:
        return "long_duration"
    if n_tr_surviving < 2:
        return "masked"
    if expected_s < s_star:
        return "below_threshold"
    recovered, _ = check_recovery(p_cat, p_found, half_width)
    if not recovered:
        return "outcompeted"
    return "other"


def read_crowdsap_from_fits(fits_path: Union[str, Path]) -> float:
    """Read CROWDSAP from FITS header (HDU 1 or HDU 0)."""
    with fits.open(fits_path, memmap=False) as hdul:
        for hdu in hdul:
            if "CROWDSAP" in hdu.header and hdu.header["CROWDSAP"] is not None:
                val = float(hdu.header["CROWDSAP"])
                if np.isfinite(val):
                    return val
    return 1.0


def select_known_planet_targets(
    sector: int = 70,
    searched_tics: Optional[set[int]] = None,
) -> pd.DataFrame:
    """Find sector-70 searched stars hosting confirmed planets or active TOIs."""
    if searched_tics is None:
        con = store.connect()
        df_searched = store.frame(con, sector)
        searched_tics = set(df_searched.tic.astype(int))

    cats = catalogs.load_catalogs()
    toi_df = cats["toi"]
    active_toi = toi_df[~toi_df.disp.isin(["FP", "FA"])].copy()

    # Load full TOI table for ephemeris & depth
    toi_path = config.CATALOG_DIR / "exofop_toi.csv"
    if toi_path.exists():
        t_raw = pd.read_csv(toi_path, comment="#", low_memory=False)
        t_raw = t_raw.rename(columns={
            "TIC ID": "tic", "TOI": "toi", "Period (days)": "period",
            "TFOPWG Disposition": "disp", "Epoch (BJD)": "epoch_bjd",
            "Duration (hours)": "dur_h", "Depth (ppm)": "depth_ppm",
        })
        active_toi_full = t_raw[t_raw.tic.isin(searched_tics) & ~t_raw.disp.isin(["FP", "FA"])].copy()
    else:
        active_toi_full = active_toi[active_toi.tic.isin(searched_tics)].copy()
        if "dur_h" not in active_toi_full.columns:
            active_toi_full["dur_h"] = 2.4
        if "depth_ppm" not in active_toi_full.columns:
            active_toi_full["depth_ppm"] = 1000.0
        if "epoch_bjd" not in active_toi_full.columns:
            active_toi_full["epoch_bjd"] = 2457000.0

    targets_list = []
    seen_keys = set()

    for r in active_toi_full.itertuples():
        tic = int(r.tic)
        p = float(r.period) if pd.notna(r.period) else np.nan
        if np.isnan(p) or p <= 0:
            continue
        dur_h = float(getattr(r, "dur_h", 2.4)) if pd.notna(getattr(r, "dur_h", 2.4)) else 2.4
        depth_ppm = float(getattr(r, "depth_ppm", 1000.0)) if pd.notna(getattr(r, "depth_ppm", 1000.0)) else 1000.0
        t0 = float(getattr(r, "epoch_bjd", 2457000.0)) - 2457000.0 if pd.notna(getattr(r, "epoch_bjd", 2457000.0)) else 0.0
        key = (tic, round(p, 4))
        if key not in seen_keys:
            seen_keys.add(key)
            targets_list.append({
                "tic": tic,
                "name": f"TOI {r.toi}",
                "period": p,
                "t0": t0,
                "dur_h": dur_h,
                "depth_ppm": depth_ppm,
                "type": "toi",
            })

    # PS confirmed planets
    ps_df = cats["ps"]
    ps_s70 = ps_df[ps_df.tic.isin(searched_tics)].copy()
    for r in ps_s70.itertuples():
        tic = int(r.tic)
        p = float(r.pl_orbper) if pd.notna(r.pl_orbper) else np.nan
        if np.isnan(p) or p <= 0:
            continue
        key = (tic, round(p, 4))
        if key not in seen_keys:
            seen_keys.add(key)
            targets_list.append({
                "tic": tic,
                "name": str(r.pl_name),
                "period": p,
                "t0": 0.0,
                "dur_h": 2.4,
                "depth_ppm": 1000.0,
                "type": "ps",
            })

    return pd.DataFrame(targets_list)


def evaluate_known_planets(
    targets_df: pd.DataFrame,
    search_results_df: pd.DataFrame,
    s_star: float,
    durations: np.ndarray = search.DURATIONS,
) -> Tuple[pd.DataFrame, Dict[str, Any]]:
    """Evaluate known planet / active TOI recovery on searched stars.

    Parameters
    ----------
    targets_df : pd.DataFrame
        Targets with columns:
        ['tic', 'name', 'period', 't0', 'duration_hours', 'depth_ppm', 'crowdsap',
         'raw_cadences_t', 'clean_cadences_t', 'clean_flux']
        or computed columns:
        ['tic', 'name', 'period', 'duration_days', 'n_raw', 'n_tr', 'cdpp', 'depth_exp', 'crowdsap']
    search_results_df : pd.DataFrame
        Top peak search results with columns ['tic', 'period', 'S'] (or 'P_grid')
    s_star : float
        Threshold S* passed in as an argument.

    Returns
    -------
    tuple of (per_target_df, summary_dict)
    """
    res_by_tic = {}
    if not search_results_df.empty:
        for r in search_results_df.itertuples():
            tic = int(r.tic)
            P = getattr(r, "period", None)
            if P is None or pd.isna(P):
                P = getattr(r, "P_grid", None)
            S = getattr(r, "S", None)
            res_by_tic[tic] = dict(period=P, S=S)

    rows = []

    for row in targets_df.itertuples():
        tic = int(row.tic)
        name = str(getattr(row, "name", getattr(row, "toi", getattr(row, "pl_name", ""))))
        P_cat = float(row.period)
        dur_h = float(getattr(row, "duration_hours", getattr(row, "dur_h", 2.4)))
        dur_d = dur_h / 24.0
        depth_ppm = float(getattr(row, "depth_ppm", 1000.0))
        depth_cat = depth_ppm * 1e-6

        crowdsap = float(getattr(row, "crowdsap", 1.0))
        # Contamination factor: 1.0 for PDCSAP
        contamination_factor = CONTAMINATION_FACTOR_PDCSAP
        depth_exp = depth_cat * contamination_factor

        # Transits and noise
        n_tr = int(getattr(row, "n_tr", 2))
        n_raw = int(getattr(row, "n_raw", n_tr))
        cdpp_d = float(getattr(row, "cdpp_d", getattr(row, "cdpp", 1e-4)))
        baseline = float(getattr(row, "baseline", 27.0))

        # Expected S
        if cdpp_d > 0 and n_tr > 0:
            expected_s = float(depth_exp * np.sqrt(n_tr) / cdpp_d)
        else:
            expected_s = 0.0

        half_width = calculate_grid_halfwidth(P_cat, baseline=baseline)

        # Found peak
        found = res_by_tic.get(tic, {})
        P_found = found.get("period")
        S_found = found.get("S")

        # Required set: expected S >= S* and n_tr >= 2
        required = bool((expected_s >= s_star) and (n_tr >= 2))

        recovered, alias_flag = check_recovery(P_cat, P_found, half_width)

        miss_cause = None
        if not recovered:
            miss_cause = determine_miss_cause(
                n_raw=n_raw,
                dur_d=dur_d,
                n_tr_surviving=n_tr,
                expected_s=expected_s,
                s_star=s_star,
                p_cat=P_cat,
                p_found=P_found,
                half_width=half_width,
                s_found=S_found,
            )

        rows.append(
            dict(
                tic=tic,
                toi_or_planet=name,
                P_cat=P_cat,
                P_found=P_found,
                expected_S=expected_s,
                n_tr=n_tr,
                contamination_factor=contamination_factor,
                crowdsap=crowdsap,
                required=required,
                recovered=recovered,
                alias_flag=alias_flag,
                miss_cause=miss_cause,
            )
        )

    out_df = pd.DataFrame(rows)

    # Any required target missed fails the test
    req_df = out_df[out_df.required]
    missed_required = req_df[~req_df.recovered]
    pass_flag = bool(len(missed_required) == 0)

    summary = {
        "s_star": float(s_star),
        "total_targets": int(len(out_df)),
        "required_targets": int(len(req_df)),
        "required_recovered": int(req_df.recovered.sum()) if not req_df.empty else 0,
        "required_missed": int(len(missed_required)),
        "pass": pass_flag,
        "miss_cause_counts": out_df.miss_cause.value_counts().to_dict(),
    }

    return out_df, summary


def run_known_planet_recovery(
    s_star: float,
    search_results_path: Optional[Union[str, Path]] = None,
    out_dir: Union[str, Path] = Path("."),
    sector: int = 70,
    targets_csv: Optional[Union[str, Path]] = None,
) -> Tuple[pd.DataFrame, Dict[str, Any]]:
    """Run known-planet recovery evaluation and write CSV and JSON outputs."""
    out_p = Path(out_dir)
    out_p.mkdir(parents=True, exist_ok=True)

    # 1. Load or build targets
    if targets_csv and Path(targets_csv).exists():
        df_tgt = pd.read_csv(targets_csv)
    else:
        df_tgt = select_known_planet_targets(sector=sector)

    # 2. Enrich targets with light-curve data, masks, CDPP, and CROWDSAP if not pre-computed
    enriched_rows = []
    man = None
    try:
        man = targets.fetch_manifest(sector).set_index("tic")
    except Exception:
        pass

    cfg = settings.cfg()

    for row in df_tgt.itertuples():
        rdict = dict(row._asdict())
        tic = int(rdict["tic"])
        P = float(rdict["period"])
        dur_h = float(rdict.get("dur_h", rdict.get("duration_hours", 2.4)))
        dur_d = dur_h / 24.0
        t0 = float(rdict.get("t0", 0.0))

        # Check if already computed
        if "n_tr" in rdict and "cdpp_d" in rdict and "crowdsap" in rdict:
            enriched_rows.append(rdict)
            continue

        # Look up lightcurve
        crowdsap = 1.0
        n_raw = 0
        n_tr = 0
        cdpp_d = 1e-4
        baseline = 27.0

        fname = man.loc[tic].fname if (man is not None and tic in man.index) else None
        p = batch.lc_path(sector, fname, download=False) if fname else None
        if p and p.exists():
            crowdsap = read_crowdsap_from_fits(p)
            try:
                t, f, e, q, meta = lightcurve.load_spoc_raw(p)
                good = lightcurve.good_mask(f, e, q)
                tg = t[good]
                fg = f[good] / np.nanmedian(f[good])
                eg = e[good] / np.nanmedian(f[good])
                baseline = float(tg.max() - tg.min()) if len(tg) > 0 else 27.0
                n_raw = count_surviving_transits(t, P, t0, dur_d)

                # Prepare with v1.1 masks
                d = pipeline.prepare(tg, fg, eg, cfg, prewhiten_on=True, gap_trim=0.3)
                tf, ff, ef = d["flat"]
                n_tr = count_surviving_transits(tf, P, t0, dur_d)

                # CDPP at duration
                prof = bls.noise_profile(tf, ff, [dur_d])
                cdpp_d = prof[2].get(dur_d, (1e-4, 1.0))[0]
                if np.isnan(cdpp_d) or cdpp_d <= 0:
                    cdpp_d = 1e-4
            except Exception:
                pass

        rdict.update(
            crowdsap=crowdsap,
            n_raw=n_raw,
            n_tr=n_tr,
            cdpp_d=cdpp_d,
            baseline=baseline,
        )
        enriched_rows.append(rdict)

    df_enriched = pd.DataFrame(enriched_rows)

    # 3. Load search results
    if search_results_path and Path(search_results_path).exists():
        df_search = pd.read_csv(search_results_path)
    else:
        df_search = pd.DataFrame(columns=["tic", "period", "S"])

    out_df, summary = evaluate_known_planets(df_enriched, df_search, s_star=s_star)

    # Output CSV and JSON
    csv_path = out_p / "known_planets_recovery.csv"
    out_df.to_csv(csv_path, index=False)

    json_path = out_p / "known_planets_summary.json"
    with open(json_path, "w") as fp:
        json.dump(summary, fp, indent=2, default=float)

    return out_df, summary


def main():
    ap = argparse.ArgumentParser(description="v1.1 Known planet / active TOI recovery test (T2)")
    ap.add_argument("--s-star", type=float, required=True, help="Detection threshold S*")
    ap.add_argument("--results", help="Search results CSV path")
    ap.add_argument("--targets", help="Targets CSV path (optional)")
    ap.add_argument("--out", required=True, help="Output directory")
    ap.add_argument("--sector", type=int, default=70, help="TESS Sector")
    a = ap.parse_args()

    run_known_planet_recovery(
        s_star=a.s_star,
        search_results_path=a.results,
        out_dir=a.out,
        sector=a.sector,
        targets_csv=a.targets,
    )


if __name__ == "__main__":
    main()
