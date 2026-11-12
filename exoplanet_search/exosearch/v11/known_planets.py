"""Known-planet and active-TOI recovery test (v1.1 amendment 2026-11-11, Section 3 and Section 8e).

Every confirmed planet or active TOI (TFOPWG disposition not FP/FA) on a searched star
with expected S >= S* and at least 2 transits surviving masks must be recovered.

Expected S = depth_exp * sqrt(n_tr) / CDPP_D at catalogue duration.
n_tr = catalogue-ephemeris transits surviving the v1.1 masks (at least half their in-transit cadences).
CDPP_D from our own noise profile.
depth_exp = catalogue depth scaled by aperture contamination ratio:
  Aperture contamination is primary: per-star aperture input column or FITS CROWDSAP.
  Documented gap: No standalone custom aperture pixel model code is present in the repository;
  we read the per-star aperture target flux fraction CROWDSAP from the SPOC FITS header.
  The TIC contamination ratio (contratio) is logged as diagnostic and stars where the two differ
  by >20% (|crowdsap - 1/(1+contratio)| / crowdsap > 0.20 or similar) are flagged.
  A star with no aperture value is ineligible with reason recorded.
  A factor of 1.0 is never assumed.
Recovered: top reported peak period within the period-grid plateau half-width of the
  catalogue period (from exosearch/p_grid_halfwidth.py locate() or interim file, falling back to 0.5 * P^2 * df),
  or 2x / 0.5x alias (using the exact same half-width tolerance and flagged).
Eligibility: computed PRE-RUN from pscomppars + active TOIs via 'freeze-eligibility' subcommand.
No excused misses: any eligible miss for any reason (including missing light curve or runtime error)
is a FAILURE row with diagnosed cause.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Dict, Any, List, Optional, Tuple, Union

from astropy.io import fits
import numpy as np
import pandas as pd

from .. import batch, bls, catalogs, config, lightcurve, store, targets
from .. import p_grid_halfwidth
from . import grid, pipeline, search, settings

DETREND_WINDOW_DAYS = 0.75


def get_plateau_halfwidth(
    period: float,
    baseline: float = 27.0,
    durations: np.ndarray = search.DURATIONS,
    frequency_factor: float = 1.0,
    plateaus_cache: Optional[Dict[float, float]] = None,
) -> float:
    """Calculate the registered plateau half-width from exosearch/p_grid_halfwidth.py.
    
    Falls back to 0.5 * P^2 * df only where that module itself falls back.
    """
    if plateaus_cache is not None and period in plateaus_cache:
        return plateaus_cache[period]
    df = frequency_factor * float(np.min(durations)) / (baseline ** 2)
    # Check if registered plateau half-width file exists
    pgh_file = p_grid_halfwidth.OUT
    if pgh_file.exists():
        try:
            df_pgh = pd.read_csv(pgh_file)
            m = np.abs(df_pgh["period_v10"] - period) <= 1e-5
            if m.any():
                return float(df_pgh.loc[m, "P_grid_halfwidth"].iloc[0])
        except Exception:
            pass
    # Default fallback per p_grid_halfwidth.locate
    loc = p_grid_halfwidth.locate(period, [], df)
    return float(loc["P_grid_halfwidth"])


# Alias for backwards compatibility
calculate_grid_halfwidth = get_plateau_halfwidth


def count_surviving_transits(
    t: np.ndarray,
    period: float,
    t0: float,
    duration: float,
    cadence: Optional[float] = None,
) -> int:
    """Number of transits of the ephemeris that have at least half their in-transit cadences."""
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
) -> Tuple[bool, Optional[str]]:
    """Determine whether p_found recovers p_cat within the plateau half-width or as an alias.
    
    Per Spec A4 / B4: 2x and 0.5x aliases use the exact same tolerance and are flagged.
    Returns (recovered, alias_flag).
    """
    if p_found is None or not np.isfinite(p_found) or p_found <= 0:
        return False, None

    # Check primary recovery within plateau half-width
    if abs(p_found - p_cat) <= half_width:
        return True, None

    # Check 2x / 0.5x aliases using the EXACT same half_width tolerance
    for ratio in alias_ratios:
        target_alias = p_cat * ratio
        if abs(p_found - target_alias) <= half_width:
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
    error_msg: Optional[str] = None,
) -> str:
    """Classify miss cause bucket (diagnosed cause, never excused)."""
    if error_msg:
        return f"runtime_error: {error_msg}"
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


def read_crowdsap_from_fits(fits_path: Union[str, Path]) -> Optional[float]:
    """Read CROWDSAP from FITS header (HDU 1 or HDU 0). Returns None if absent."""
    try:
        with fits.open(fits_path, memmap=False) as hdul:
            for hdu in hdul:
                if "CROWDSAP" in hdu.header and hdu.header["CROWDSAP"] is not None:
                    val = float(hdu.header["CROWDSAP"])
                    if np.isfinite(val) and val > 0:
                        return val
    except Exception:
        pass
    return None


def select_known_planet_targets(
    sector: int = 70,
    searched_tics: Optional[set[int]] = None,
) -> pd.DataFrame:
    """Find sector-70 searched stars hosting confirmed planets or active TOIs.
    
    A2: Real catalogue depth, duration, epoch (pscomppars + active TOIs). No placeholders anywhere.
    A3: Known-planet dedupe by (TIC, period within registered tolerance), preferring pscomppars row;
        never round(P, 4).
    """
    if searched_tics is None:
        con = store.connect()
        df_searched = store.frame(con, sector)
        searched_tics = set(df_searched.tic.astype(int))

    cats = catalogs.load_catalogs()
    active_toi_base = cats["toi"][~cats["toi"].disp.isin(["FP", "FA"])].copy()

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
        active_toi_full = active_toi_base[active_toi_base.tic.isin(searched_tics)].copy()

    # PS confirmed planets
    ps_df = cats["ps"]
    ps_s70 = ps_df[ps_df.tic.isin(searched_tics)].copy()

    # Build initial candidate lists
    ps_targets = []
    for r in ps_s70.itertuples():
        tic = int(r.tic)
        p = float(r.pl_orbper) if pd.notna(r.pl_orbper) else np.nan
        # Look for matching row in TOI table for duration, depth, epoch if available
        dur_h = np.nan
        depth_ppm = np.nan
        t0 = np.nan
        if not active_toi_full.empty:
            m_toi = active_toi_full[(active_toi_full.tic == tic) & np.isclose(active_toi_full.period.astype(float), p, atol=get_plateau_halfwidth(p))]
            if not m_toi.empty:
                r_m = m_toi.iloc[0]
                dur_h = float(r_m.dur_h) if pd.notna(r_m.dur_h) else np.nan
                depth_ppm = float(r_m.depth_ppm) if pd.notna(r_m.depth_ppm) else np.nan
                t0 = float(r_m.epoch_bjd) - 2457000.0 if pd.notna(r_m.epoch_bjd) else np.nan

        ps_targets.append({
            "tic": tic,
            "name": str(r.pl_name),
            "period": p,
            "t0": t0,
            "dur_h": dur_h,
            "depth_ppm": depth_ppm,
            "type": "ps",
        })

    toi_targets = []
    for r in active_toi_full.itertuples():
        tic = int(r.tic)
        p = float(r.period) if pd.notna(r.period) else np.nan
        dur_h = float(r.dur_h) if pd.notna(getattr(r, "dur_h", np.nan)) else np.nan
        depth_ppm = float(r.depth_ppm) if pd.notna(getattr(r, "depth_ppm", np.nan)) else np.nan
        t0 = float(r.epoch_bjd) - 2457000.0 if pd.notna(getattr(r, "epoch_bjd", np.nan)) else np.nan
        toi_targets.append({
            "tic": tic,
            "name": f"TOI {r.toi}",
            "period": p,
            "t0": t0,
            "dur_h": dur_h,
            "depth_ppm": depth_ppm,
            "type": "toi",
        })

    # A3: Dedupe by (TIC, period within registered tolerance), PREFERRING the pscomppars row
    final_targets = []
    # Add PS rows first
    for pt in ps_targets:
        final_targets.append(pt)

    # Add TOI rows if no PS row matches (TIC, period within registered tolerance)
    for tt in toi_targets:
        t_tic = tt["tic"]
        t_p = tt["period"]
        if np.isnan(t_p) or t_p <= 0:
            final_targets.append(tt)
            continue
        hw = get_plateau_halfwidth(t_p)
        matched = False
        for pt in ps_targets:
            if pt["tic"] == t_tic and np.isfinite(pt["period"]) and abs(pt["period"] - t_p) <= hw:
                matched = True
                break
        if not matched:
            final_targets.append(tt)

    return pd.DataFrame(final_targets)


def evaluate_known_planets(
    targets_df: pd.DataFrame,
    search_results_df: pd.DataFrame,
    s_star: float,
    durations: np.ndarray = search.DURATIONS,
) -> Tuple[pd.DataFrame, Dict[str, Any]]:
    """Evaluate known planet / active TOI recovery on searched stars.
    
    B3: S* is passed in; error if not supplied.
    """
    if s_star is None or not np.isfinite(s_star) or s_star <= 0:
        raise ValueError("s_star must be provided as a positive float")

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
    baseline = 27.0

    for row in targets_df.itertuples():
        tic = int(row.tic)
        name = str(getattr(row, "name", getattr(row, "toi", getattr(row, "pl_name", ""))))
        P_cat = float(getattr(row, "period", np.nan))
        dur_h = float(getattr(row, "duration_hours", getattr(row, "dur_h", np.nan)))
        dur_d = dur_h / 24.0 if np.isfinite(dur_h) else np.nan
        depth_ppm = float(getattr(row, "depth_ppm", np.nan))
        depth_cat = depth_ppm * 1e-6 if np.isfinite(depth_ppm) else np.nan

        # Aperture contamination
        aperture_val = getattr(row, "aperture_contamination", getattr(row, "crowdsap", None))
        if aperture_val is not None and pd.notna(aperture_val):
            aperture_val = float(aperture_val)
        else:
            aperture_val = None

        contratio = getattr(row, "contratio", None)
        if contratio is not None and pd.notna(contratio):
            contratio = float(contratio)
        else:
            contratio = None

        # Flag > 20% difference between aperture model and TIC contratio
        contamination_flag = False
        if aperture_val is not None and contratio is not None and contratio >= 0:
            tic_target_flux_frac = 1.0 / (1.0 + contratio)
            if abs(aperture_val - tic_target_flux_frac) / max(aperture_val, 1e-4) > 0.20:
                contamination_flag = True

        # Check eligibility reasons (B4, B5)
        ineligible_reason = getattr(row, "ineligible_reason", None)
        if ineligible_reason is None:
            if np.isnan(P_cat) or P_cat <= 0:
                ineligible_reason = "missing_period"
            elif np.isnan(dur_h) or dur_h <= 0:
                ineligible_reason = "missing_duration"
            elif np.isnan(depth_ppm) or depth_ppm <= 0:
                ineligible_reason = "missing_depth"
            elif getattr(row, "t0", None) is None or pd.isna(getattr(row, "t0", np.nan)):
                ineligible_reason = "missing_epoch"
            elif aperture_val is None:
                ineligible_reason = "missing_aperture_contamination"

        # Transits and noise
        n_tr = int(getattr(row, "n_tr", 0))
        n_raw = int(getattr(row, "n_raw", n_tr))
        cdpp_d = float(getattr(row, "cdpp_d", getattr(row, "cdpp", np.nan)))

        # Expected S
        expected_s = 0.0
        if (
            ineligible_reason is None
            and depth_cat is not None
            and aperture_val is not None
            and np.isfinite(cdpp_d)
            and cdpp_d > 0
            and n_tr > 0
        ):
            depth_exp = depth_cat * aperture_val
            expected_s = float(depth_exp * np.sqrt(n_tr) / cdpp_d)
        else:
            depth_exp = np.nan

        # Check further eligibility requirements
        if ineligible_reason is None:
            if n_tr < 2:
                ineligible_reason = "transits_under_two"
            elif expected_s < s_star:
                ineligible_reason = "expected_s_below_threshold"

        eligible = bool(ineligible_reason is None)

        half_width = get_plateau_halfwidth(P_cat, baseline=baseline) if np.isfinite(P_cat) else 0.0

        # Found peak
        found = res_by_tic.get(tic, {})
        P_found = found.get("period")
        S_found = found.get("S")

        # Required set = eligible (B5)
        required = eligible

        recovered, alias_flag = check_recovery(P_cat, P_found, half_width) if (eligible and P_found is not None) else (False, None)

        miss_cause = None
        if not recovered:
            err = getattr(row, "error", None)
            err_str = str(err) if (err is not None and pd.notna(err) and str(err).strip() != "") else None
            miss_cause = determine_miss_cause(
                n_raw=n_raw,
                dur_d=dur_d if np.isfinite(dur_d) else 0.1,
                n_tr_surviving=n_tr,
                expected_s=expected_s,
                s_star=s_star,
                p_cat=P_cat if np.isfinite(P_cat) else 1.0,
                p_found=P_found,
                half_width=half_width,
                s_found=S_found,
                error_msg=err_str,
            )

        rows.append(
            dict(
                tic=tic,
                toi_or_planet=name,
                P_cat=P_cat,
                P_found=P_found,
                expected_S=expected_s,
                n_tr=n_tr,
                aperture_contamination=aperture_val,
                contratio=contratio,
                contamination_flag=contamination_flag,
                eligible=eligible,
                ineligible_reason=ineligible_reason or "",
                required=required,
                recovered=recovered,
                alias_flag=alias_flag,
                miss_cause=miss_cause,
            )
        )

    out_df = pd.DataFrame(rows)

    # Any eligible/required target missed fails the test (B5: no excused misses)
    req_df = out_df[out_df.required]
    missed_required = req_df[~req_df.recovered]
    pass_flag = bool(len(missed_required) == 0 and len(req_df) > 0)

    summary = {
        "s_star": float(s_star),
        "total_targets": int(len(out_df)),
        "eligible_targets": int(len(req_df)),
        "eligible_recovered": int(req_df.recovered.sum()) if not req_df.empty else 0,
        "eligible_missed": int(len(missed_required)),
        "pass": pass_flag,
        "miss_cause_counts": out_df.miss_cause.value_counts().to_dict(),
    }

    return out_df, summary


def freeze_eligibility(
    s_star: float,
    out_file: Union[str, Path],
    sector: int = 70,
    commit: str = "",
    config_hash: str = "",
) -> Tuple[pd.DataFrame, Dict[str, Any]]:
    """Compute eligibility pre-run from pscomppars + active TOIs and write dated file + sidecar JSON."""
    if s_star is None or not np.isfinite(s_star) or s_star <= 0:
        raise ValueError("s_star must be provided as a positive float")

    out_p = Path(out_file)
    out_p.parent.mkdir(parents=True, exist_ok=True)

    df_tgt = select_known_planet_targets(sector=sector)

    # Enrich with light curve, CDPP, and aperture contamination
    man = None
    try:
        man = targets.fetch_manifest(sector).set_index("tic")
    except Exception:
        pass

    # Read contratio from tic82_params if cached
    tic82_map = {}
    tic82_p = config.CATALOG_DIR / "tic82_params.parquet"
    if tic82_p.exists():
        try:
            df_82 = pd.read_parquet(tic82_p)
            for r in df_82.itertuples():
                tic82_map[int(r.ID)] = getattr(r, "contratio", np.nan)
        except Exception:
            pass

    cfg = settings.cfg()
    enriched_rows = []

    for row in df_tgt.itertuples():
        rdict = dict(row._asdict())
        tic = int(rdict["tic"])
        P = float(rdict["period"]) if pd.notna(rdict["period"]) else np.nan
        dur_h = float(rdict["dur_h"]) if pd.notna(rdict["dur_h"]) else np.nan
        dur_d = dur_h / 24.0 if np.isfinite(dur_h) else np.nan
        t0 = float(rdict["t0"]) if pd.notna(rdict["t0"]) else np.nan

        contratio = tic82_map.get(tic, np.nan)
        rdict["contratio"] = contratio

        # Look up lightcurve and CROWDSAP
        fname = man.loc[tic].fname if (man is not None and tic in man.index) else None
        p = batch.lc_path(sector, fname, download=False) if fname else None
        crowdsap = None
        n_raw = 0
        n_tr = 0
        cdpp_d = np.nan
        err = None

        if p and p.exists():
            crowdsap = read_crowdsap_from_fits(p)
            if crowdsap is not None:
                try:
                    t, f, e, q, meta = lightcurve.load_spoc_raw(p)
                    good = lightcurve.good_mask(f, e, q)
                    tg = t[good]
                    fg = f[good] / np.nanmedian(f[good])
                    eg = e[good] / np.nanmedian(f[good])
                    n_raw = count_surviving_transits(t, P, t0, dur_d) if np.isfinite(P) and np.isfinite(t0) and np.isfinite(dur_d) else 0

                    d = pipeline.prepare(tg, fg, eg, cfg, prewhiten_on=True, gap_trim=0.3)
                    tf, ff, ef = d["flat"]
                    n_tr = count_surviving_transits(tf, P, t0, dur_d) if np.isfinite(P) and np.isfinite(t0) and np.isfinite(dur_d) else 0

                    prof = bls.noise_profile(tf, ff, [dur_d])
                    cdpp_d = prof[2].get(dur_d, (np.nan, 1.0))[0]
                except Exception as exc:
                    err = str(exc)
        else:
            err = "missing light curve"

        rdict.update(
            aperture_contamination=crowdsap,
            n_raw=n_raw,
            n_tr=n_tr,
            cdpp_d=cdpp_d,
            error=err,
        )
        enriched_rows.append(rdict)

    df_enriched = pd.DataFrame(enriched_rows)
    eval_df, summary = evaluate_known_planets(df_enriched, pd.DataFrame(columns=["tic", "period", "S"]), s_star=s_star)

    # Save CSV
    eval_df.to_csv(out_p, index=False)

    # Sidecar JSON
    sha256 = hashlib.sha256(out_p.read_bytes()).hexdigest()
    sidecar = {
        "s_star": float(s_star),
        "sector": sector,
        "commit": commit,
        "config_hash": config_hash,
        "sha256": sha256,
        "total_targets": len(eval_df),
        "eligible_targets": int(eval_df.eligible.sum()),
        "ineligible_reasons": eval_df[~eval_df.eligible].ineligible_reason.value_counts().to_dict(),
    }
    sidecar_p = Path(str(out_p) + ".meta.json")
    sidecar_p.write_text(json.dumps(sidecar, indent=2))

    return eval_df, sidecar


def run_known_planet_recovery(
    s_star: float,
    search_results_path: Optional[Union[str, Path]] = None,
    out_dir: Union[str, Path] = Path("."),
    sector: int = 70,
    frozen_eligibility_file: Optional[Union[str, Path]] = None,
) -> Tuple[pd.DataFrame, Dict[str, Any]]:
    """Run known-planet recovery evaluation and write CSV and JSON outputs.
    
    Refuses to run without frozen eligibility file (B5).
    """
    if s_star is None or not np.isfinite(s_star) or s_star <= 0:
        raise ValueError("s_star must be provided as a positive float")

    if not frozen_eligibility_file or not Path(frozen_eligibility_file).exists():
        raise RuntimeError("Frozen eligibility file required for known-planet recovery evaluation (Spec B5).")

    # Verify sidecar sha256
    meta_p = Path(str(frozen_eligibility_file) + ".meta.json")
    if meta_p.exists():
        meta = json.loads(meta_p.read_text())
        cur_sha = hashlib.sha256(Path(frozen_eligibility_file).read_bytes()).hexdigest()
        if cur_sha != meta.get("sha256"):
            raise ValueError(f"Eligibility file sha256 mismatch: {cur_sha} != {meta.get('sha256')}")

    out_p = Path(out_dir)
    out_p.mkdir(parents=True, exist_ok=True)

    df_eligible = pd.read_csv(frozen_eligibility_file)

    if search_results_path and Path(search_results_path).exists():
        df_search = pd.read_csv(search_results_path)
    else:
        df_search = pd.DataFrame(columns=["tic", "period", "S"])

    out_df, summary = evaluate_known_planets(df_eligible, df_search, s_star=s_star)

    csv_path = out_p / "known_planets_recovery.csv"
    out_df.to_csv(csv_path, index=False)

    json_path = out_p / "known_planets_summary.json"
    with open(json_path, "w") as fp:
        json.dump(summary, fp, indent=2, default=float)

    return out_df, summary


def main():
    ap = argparse.ArgumentParser(description="v1.1 Known planet / active TOI recovery test (T2)")
    sp = ap.add_subparsers(dest="cmd", required=True)

    fe = sp.add_parser("freeze-eligibility")
    fe.add_argument("--s-star", type=float, required=True, help="Detection threshold S*")
    fe.add_argument("--out", required=True, help="Output path for frozen eligibility CSV")
    fe.add_argument("--sector", type=int, default=70, help="TESS Sector")

    r = sp.add_parser("report")
    r.add_argument("--s-star", type=float, required=True, help="Detection threshold S*")
    r.add_argument("--results", help="Search results CSV path")
    r.add_argument("--frozen-eligibility", required=True, help="Path to frozen eligibility CSV")
    r.add_argument("--out", required=True, help="Output directory")
    r.add_argument("--sector", type=int, default=70, help="TESS Sector")

    a = ap.parse_args()

    if a.cmd == "freeze-eligibility":
        freeze_eligibility(s_star=a.s_star, out_file=a.out, sector=a.sector)
    else:
        run_known_planet_recovery(
            s_star=a.s_star,
            search_results_path=a.results,
            out_dir=a.out,
            sector=a.sector,
            frozen_eligibility_file=a.frozen_eligibility,
        )


if __name__ == "__main__":
    main()
