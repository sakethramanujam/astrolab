"""Known-planet and active-TOI recovery test (v1.1 amendment 2026-11-11, Section 3 and Section 8e).

Every confirmed planet or active TOI (TFOPWG disposition not FP/FA) on a searched star
with expected S >= S* and at least 2 transits surviving masks must be recovered.

Expected S = depth_exp * sqrt(n_tr) / CDPP_D at catalogue duration (computed ONCE, pre-run, by
'freeze-eligibility'; the report never recomputes it).
n_tr = catalogue-ephemeris transits surviving the v1.1 masks (at least half their in-transit cadences).
CDPP_D from our own noise profile.
depth_exp = catalogue depth scaled by aperture contamination ratio:
  Aperture contamination is primary: per-star aperture input column or FITS CROWDSAP.
  Documented gap: No standalone custom aperture pixel model code is present in the repository;
  we read the per-star aperture target flux fraction CROWDSAP from the SPOC FITS header.
  The TIC contamination ratio (contratio) is logged as diagnostic and stars where the two differ
  by >20% are flagged. A star with no aperture value is ineligible with reason recorded.
  A factor of 1.0 is never assumed.
Catalogue parameters: confirmed planets from the NASA Exoplanet Archive pscomppars table (fetched over TAP with
  pl_trandep, pl_trandur, pl_tranmid and cached under a dated filename with a sha256 sidecar); active TOIs from the
  ExoFOP TOI table (Depth (ppm), Duration (hours), Epoch (BJD)). Only transiting pscomppars rows (tran_flag = 1). A
  pscomppars row missing a transit parameter takes it from the matching TOI row, else from its own catalogue values
  (depth = (Rp/R*)^2; T14 from P, a/R*, b, Rp/R*); the source of every parameter is recorded. A TOI within 1% in
  period of a confirmed planet on the same star is that planet (dropped). eph_sigma_h: catalogue ephemeris
  uncertainty propagated to mid-S70 (diagnostic).
Recovered: the star's reported peak period lies within the star's own plateau half-width of the catalogue period
  (p_grid_halfwidth.locate on that star's search plateaus and grid spacing; falling back to 0.5 * P^2 * df only where
  locate itself falls back), with S >= S*. 2x and 0.5x aliases use the same tolerance and are reported as aliases;
  they are NOT recoveries.
No excused misses: any eligible miss for any reason (including missing light curve or runtime error)
is a FAILURE row with diagnosed cause. A light-curve or runtime error during freeze-eligibility is an error row
(required, failing), never an ineligible row.
"""
from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import io
import json
from pathlib import Path
from typing import Dict, Any, List, Optional, Tuple, Union

from astropy.io import fits
import numpy as np
import pandas as pd

from .. import batch, bls, config, lightcurve, store, targets
from .. import p_grid_halfwidth
from . import pipeline, search, settings

DETREND_WINDOW_DAYS = 0.75

PS_TAP_URL = "https://exoplanetarchive.ipac.caltech.edu/TAP/sync"
PS_COLUMNS = ("pl_name", "hostname", "tic_id", "tran_flag", "pl_orbper", "pl_orbpererr1", "pl_trandep", "pl_trandur",
              "pl_tranmid", "pl_tranmiderr1", "pl_ratror", "pl_ratdor", "pl_imppar", "disc_facility")
PS_QUERY = "select " + ",".join(PS_COLUMNS) + " from pscomppars where tic_id is not null"
PS_PREFIX = "nasa_pscomppars_full_"
TOI_FILE = "exofop_toi.csv"
TOI_RENAME = {"TIC ID": "tic", "TOI": "toi", "Period (days)": "period", "Period (days) err": "period_err",
              "TFOPWG Disposition": "disp", "Epoch (BJD)": "epoch_bjd", "Epoch (BJD) err": "epoch_err",
              "Duration (hours)": "dur_h", "Depth (ppm)": "depth_ppm"}
S70_MID_BTJD = 3221.1          # mid-sector reference for ephemeris propagation uncertainty (diagnostic only)
CATALOG_MATCH_REL = 0.01       # PS-vs-TOI same-planet match: |dP| / P (two planets cannot orbit 1% apart)
BTJD_OFFSET = 2457000.0

# Frozen eligibility schema (one row per target; written once pre-run, read back by the report)
FROZEN_STR = ("name", "type", "disp", "depth_src", "dur_src", "t0_src", "base_ineligible_reason", "error",
              "eligibility", "ineligible_reason")
FROZEN_FLOAT = ("period", "period_err", "t0", "dur_h", "depth_ppm", "eph_sigma_h", "aperture_contamination", "contratio",
                "cdpp_d", "t_span_d", "expected_S_raw", "expected_S")
FROZEN_INT = ("tic", "n_raw", "n_tr")
FROZEN_BOOL = ("contamination_flag", "required")
FROZEN_COLUMNS = ("tic", "name", "type", "disp", "period", "period_err", "t0", "dur_h", "depth_ppm", "depth_src", "dur_src",
                  "t0_src", "eph_sigma_h", "aperture_contamination", "contratio", "contamination_flag", "n_raw", "n_tr",
                  "cdpp_d", "t_span_d", "expected_S_raw", "expected_S", "base_ineligible_reason", "error", "eligibility",
                  "ineligible_reason", "required")

STAGES = ("error", "masked", "detrending", "below_S*", "outranked", "vetting")


def _sha256(p: Union[str, Path]) -> str:
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def _s(x) -> str:
    """String cell -> '' for None / NaN / blank."""
    if x is None:
        return ""
    if isinstance(x, float) and not np.isfinite(x):
        return ""
    s = str(x).strip()
    return "" if s.lower() == "nan" else s


def _f(x) -> float:
    try:
        v = float(x)
    except (TypeError, ValueError):
        return np.nan
    return v


# ---------------------------------------------------------------- catalogues (Spec A2)
def pscomppars_path(cache_dir: Optional[Union[str, Path]] = None) -> Optional[Path]:
    """Newest dated full-column pscomppars cache whose sidecar sha256 matches the file, or None."""
    d = Path(cache_dir or config.CATALOG_DIR)
    for p in sorted(d.glob(PS_PREFIX + "*.csv"), reverse=True):
        meta = Path(str(p) + ".meta.json")
        if meta.exists() and json.loads(meta.read_text()).get("sha256") == _sha256(p):
            return p
    return None


def fetch_pscomppars(cache_dir: Optional[Union[str, Path]] = None, date_tag: Optional[str] = None) -> Path:
    """Query pscomppars over the public NASA Exoplanet Archive TAP service with the transit columns and cache it as
    nasa_pscomppars_full_<YYYYMMDD>.csv plus a .meta.json sidecar (UTC fetch time, query, rows, sha256)."""
    import requests
    d = Path(cache_dir or config.CATALOG_DIR)
    d.mkdir(parents=True, exist_ok=True)
    now = dt.datetime.now(dt.timezone.utc)
    tag = date_tag or now.strftime("%Y%m%d")
    r = requests.get(PS_TAP_URL, params=dict(query=PS_QUERY, format="csv"), timeout=300)
    r.raise_for_status()
    if b"<html" in r.content[:500].lower():
        raise IOError("pscomppars TAP returned an HTML page instead of CSV")
    df = pd.read_csv(io.BytesIO(r.content), comment="#", low_memory=False)
    missing = [c for c in PS_COLUMNS if c not in df.columns]
    if df.empty or missing:
        raise IOError(f"pscomppars TAP result unusable (rows={len(df)}, missing columns={missing})")
    out = d / f"{PS_PREFIX}{tag}.csv"
    out.write_bytes(r.content)
    meta = dict(fetched_at_utc=now.isoformat(timespec="seconds"), source=PS_TAP_URL, query=PS_QUERY, rows=int(len(df)),
                columns=list(PS_COLUMNS), sha256=_sha256(out))
    Path(str(out) + ".meta.json").write_text(json.dumps(meta, indent=2))
    return out


def load_pscomppars(path: Optional[Union[str, Path]] = None,
                    cache_dir: Optional[Union[str, Path]] = None) -> Tuple[pd.DataFrame, Path]:
    """Dated full-column pscomppars cache -> (DataFrame with int 'tic', path). Never the 4-column v1.0 cache."""
    p = Path(path) if path else pscomppars_path(cache_dir)
    if p is None or not p.exists():
        raise FileNotFoundError("no dated full-column pscomppars cache (nasa_pscomppars_full_<date>.csv); "
                                "run 'python -m exosearch.v11.known_planets fetch-pscomppars' (Spec A2)")
    meta_p = Path(str(p) + ".meta.json")
    if not meta_p.exists() or json.loads(meta_p.read_text()).get("sha256") != _sha256(p):
        raise ValueError(f"pscomppars cache {p.name}: sidecar missing or sha256 mismatch")
    df = pd.read_csv(p, comment="#", low_memory=False)
    missing = [c for c in PS_COLUMNS if c not in df.columns]
    if missing:
        raise ValueError(f"pscomppars cache {p.name} lacks columns {missing}")
    df["tic"] = pd.to_numeric(df["tic_id"].astype(str).str.replace("TIC", "").str.strip(), errors="coerce")
    df = df[df.tic.notna()].copy()
    df["tic"] = df.tic.astype(np.int64)
    return df, p


def load_toi_table(path: Optional[Union[str, Path]] = None) -> Tuple[pd.DataFrame, Path]:
    """ExoFOP TOI table (period, epoch, duration, depth from ExoFOP's own columns) -> (DataFrame, path)."""
    p = Path(path) if path else config.CATALOG_DIR / TOI_FILE
    if not p.exists():
        raise FileNotFoundError(f"ExoFOP TOI table not found: {p}")
    raw = pd.read_csv(p, comment="#", low_memory=False)
    missing = [c for c in TOI_RENAME if c not in raw.columns]
    if missing:
        raise ValueError(f"ExoFOP TOI table lacks columns {missing}")
    df = raw.rename(columns=TOI_RENAME)
    df["tic"] = pd.to_numeric(df["tic"], errors="coerce")
    df = df[df.tic.notna()].copy()
    df["tic"] = df.tic.astype(np.int64)
    df["disp"] = df["disp"].map(_s)
    return df, p


# ---------------------------------------------------------------- tolerances and recovery (Spec A4)
def get_plateau_halfwidth(
    period: float,
    baseline: float = 27.0,
    durations: np.ndarray = search.DURATIONS,
    frequency_factor: float = 1.0,
    plateaus_list: Optional[Union[str, List[Any]]] = None,
    grid_df: Optional[float] = None,
    tic: Optional[int] = None,
    sector: Optional[int] = None,
    pgh_table: Optional[pd.DataFrame] = None,
) -> float:
    """The star's own registered plateau half-width at `period` (exosearch/p_grid_halfwidth.locate).

    plateaus_list / grid_df: that star's search plateaus ([[P_lo, P_hi, width, hw], ...], JSON allowed) and grid
    spacing in 1/d. Without them, the (tic, sector) row of p_grid_halfwidth_v1.csv is used when the period lies on
    that row's in-band plateau. Another star's row is never used. Otherwise locate's own fallback, 0.5 * P^2 * df.
    """
    if not np.isfinite(period) or period <= 0:
        return 0.0
    df = float(grid_df) if (grid_df is not None and np.isfinite(_f(grid_df)) and float(grid_df) > 0) else \
        frequency_factor * float(np.min(durations)) / (baseline ** 2)
    if plateaus_list is not None:
        pl = json.loads(plateaus_list) if isinstance(plateaus_list, str) else plateaus_list
        return float(p_grid_halfwidth.locate(period, pl, df)["P_grid_halfwidth"])
    if tic is not None and sector is not None:
        tab = pgh_table if pgh_table is not None else p_grid_halfwidth.load()
        if tab is not None and len(tab):
            m = tab[(tab.tic.astype(np.int64) == int(tic)) & (tab.sector.astype(np.int64) == int(sector))]
            for r in m.itertuples():
                in_band = str(getattr(r, "P_grid_halfwidth_in_band", True)).lower() in ("true", "1")
                if in_band and abs(float(r.period_v10) - period) <= float(r.P_grid_halfwidth):
                    return float(r.P_grid_halfwidth)
    return float(p_grid_halfwidth.locate(period, [], df)["P_grid_halfwidth"])


# Alias for backwards compatibility
calculate_grid_halfwidth = get_plateau_halfwidth


def bls_peak_tolerance(period: float, dur_h: float, t_span_d: float, period_err: float = np.nan,
                       nsigma: float = 3.0) -> float:
    """Known-gate amendment 2026-11-12, item 1: BLS peak width D * P / T_span (catalogue duration, the star's own data
    span) plus nsigma times the catalogue period error, in quadrature. 0.0 when an input is missing."""
    P, D, T = _f(period), _f(dur_h) / 24.0, _f(t_span_d)
    if not (P > 0 and D > 0 and T > 0):
        return 0.0
    sp = _f(period_err)
    sp = abs(sp) if np.isfinite(sp) else 0.0
    return float(np.hypot(D * P / T, nsigma * sp))



def count_surviving_transits(
    t: np.ndarray,
    period: float,
    t0: float,
    duration: float,
    cadence: Optional[float] = None,
) -> int:
    """Number of transits of the ephemeris that have at least half their in-transit cadences."""
    if len(t) == 0 or not (period > 0) or not (duration > 0) or not np.isfinite(t0):
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
    s_found: Optional[float] = None,
    s_star: Optional[float] = None,
    alias_scale: bool = False,
    alias_recovered: bool = False,
) -> Tuple[bool, str]:
    """-> (recovered, match_type) with match_type in {'direct', '2x', '0.5x', 'none'}.

    Recovered only for a direct match (|P_found - P_cat| <= half-width) with S_found >= S* (when S* is given).
    2x / 0.5x aliases use the exact same half-width, are reported as match types, and never count as recovered.
    Known-gate amendment 2026-11-12, item 5: alias_scale scales the tolerance with the alias period (ratio x
    half_width, as the BLS peak width and the period error both scale with P) and alias_recovered counts a flagged
    alias as recovered.
    """
    pf = _f(p_found)
    if not np.isfinite(pf) or pf <= 0 or not np.isfinite(p_cat):
        return False, "none"
    mt = "none"
    if abs(pf - p_cat) <= half_width:
        mt = "direct"
    else:
        for ratio in alias_ratios:
            if abs(pf - p_cat * ratio) <= (ratio * half_width if alias_scale else half_width):
                mt = f"{int(ratio)}x" if ratio >= 1 else "0.5x"
                break
    s_ok = True if s_star is None else bool(np.isfinite(_f(s_found)) and _f(s_found) >= s_star)
    return bool((mt == "direct" or (alias_recovered and mt != "none")) and s_ok), mt


def losing_stage(
    match_type: str,
    s_found: Optional[float],
    s_star: float,
    vetted: Optional[bool],
    error: str = "",
    s_cat: Optional[float] = None,
    n_tr_cat: Optional[float] = None,
    depth_ratio_cat: Optional[float] = None,
    detrend_loss_ratio: float = 0.5,
) -> str:
    """Stage that lost an eligible planet: error, masked, detrending, below_S*, outranked, vetting ('' if recovered).

    s_cat / n_tr_cat / depth_ratio_cat: the v1.1 statistic, surviving transit count and measured/expected depth at
    the catalogue ephemeris on the star's own searched light curve (diagnostics; optional)."""
    if error:
        return "error"
    if np.isfinite(_f(n_tr_cat)) and _f(n_tr_cat) < 2:
        return "masked"
    weak = "detrending" if (np.isfinite(_f(depth_ratio_cat)) and _f(depth_ratio_cat) < detrend_loss_ratio) else "below_S*"
    sf = _f(s_found)
    if match_type == "direct":
        if not (np.isfinite(sf) and sf >= s_star):
            return weak
        if vetted is not None and not bool(vetted):
            return "vetting"
        return ""
    sc = _f(s_cat)
    if np.isfinite(sc):
        return "outranked" if sc >= s_star else weak
    return "outranked" if (np.isfinite(sf) and sf >= s_star) else "below_S*"


def read_crowdsap_from_fits(fits_path: Union[str, Path]) -> Optional[float]:
    """Read CROWDSAP from the FITS headers. Returns None if absent (an unreadable file raises)."""
    with fits.open(fits_path, memmap=False) as hdul:
        for hdu in hdul:
            if "CROWDSAP" in hdu.header and hdu.header["CROWDSAP"] is not None:
                val = float(hdu.header["CROWDSAP"])
                if np.isfinite(val) and val > 0:
                    return val
    return None


# ---------------------------------------------------------------- targets (Spec A2, A3)
def select_known_planet_targets(
    sector: int = 70,
    searched_tics: Optional[set] = None,
    ps_df: Optional[pd.DataFrame] = None,
    toi_df: Optional[pd.DataFrame] = None,
    pgh_table: Optional[pd.DataFrame] = None,
) -> pd.DataFrame:
    """Searched stars hosting confirmed planets (pscomppars) or active TOIs (TFOPWG not FP/FA).

    A2: real catalogue depth, duration, epoch; pscomppars first, a missing pscomppars value from the matching TOI
        row (source recorded per parameter). No placeholders.
    A3: dedupe by (TIC, period within the registered tolerance), preferring the pscomppars row; never round(P, 4).
    """
    if searched_tics is None:
        con = store.connect()
        searched_tics = set(store.frame(con, sector).tic.astype(int))
    searched_tics = {int(t) for t in searched_tics}
    if ps_df is None:
        ps_df, _ = load_pscomppars()
    if toi_df is None:
        toi_df, _ = load_toi_table()
    if pgh_table is None:
        pgh_table = p_grid_halfwidth.load()

    def tol(p, tic):
        return get_plateau_halfwidth(p, tic=tic, sector=sector, pgh_table=pgh_table)

    toi = toi_df[toi_df.tic.astype(np.int64).isin(searched_tics)].copy()
    toi = toi[~toi.disp.map(_s).isin(["FP", "FA"])]
    ps = ps_df[ps_df.tic.astype(np.int64).isin(searched_tics)].copy()
    if "tran_flag" in ps.columns:                       # transiting planets only (RV-only planets have no transit)
        ps = ps[pd.to_numeric(ps.tran_flag, errors="coerce") == 1]

    def eph_sigma_h(P, sP, sT0, t0):
        if not (np.isfinite(P) and P > 0 and np.isfinite(t0)):
            return np.nan
        n = abs(S70_MID_BTJD - t0) / P
        return float(24.0 * np.hypot(_f(sT0) if np.isfinite(_f(sT0)) else 0.0, n * (_f(sP) if np.isfinite(_f(sP)) else 0.0)))

    ps_targets = []
    for r in ps.itertuples():
        tic = int(r.tic)
        p = _f(getattr(r, "pl_orbper", np.nan))
        depth_pct = _f(getattr(r, "pl_trandep", np.nan))                  # pscomppars: percent
        vals = dict(depth_ppm=depth_pct * 1e4 if np.isfinite(depth_pct) else np.nan,
                    dur_h=_f(getattr(r, "pl_trandur", np.nan)),            # hours
                    t0=_f(getattr(r, "pl_tranmid", np.nan)) - BTJD_OFFSET)  # BJD -> BTJD
        src = {k: ("pscomppars" if np.isfinite(v) else "") for k, v in vals.items()}
        if any(not np.isfinite(v) for v in vals.values()) and np.isfinite(p):
            m = toi[(toi.tic == tic) & (np.abs(toi.period.astype(float) - p) <= max(tol(p, tic), CATALOG_MATCH_REL * p))]
            if len(m):
                rm = m.iloc[0]
                fill = dict(depth_ppm=_f(rm.depth_ppm), dur_h=_f(rm.dur_h), t0=_f(rm.epoch_bjd) - BTJD_OFFSET)
                for k in vals:
                    if not np.isfinite(vals[k]) and np.isfinite(fill[k]):
                        vals[k] = fill[k]
                        src[k] = f"exofop_toi {rm.toi}"
        # last resort, still catalogue values: depth = (Rp/R*)^2, T14 from P, a/R*, b, Rp/R* (circular orbit)
        k_rr, a_r, b_imp = _f(getattr(r, "pl_ratror", np.nan)), _f(getattr(r, "pl_ratdor", np.nan)), \
            _f(getattr(r, "pl_imppar", np.nan))
        if not np.isfinite(vals["depth_ppm"]) and k_rr > 0:
            vals["depth_ppm"], src["depth_ppm"] = 1e6 * k_rr ** 2, "pscomppars pl_ratror^2"
        if not np.isfinite(vals["dur_h"]) and k_rr > 0 and a_r > 1 and np.isfinite(b_imp) and p > 0:
            arg = np.sqrt(max((1 + k_rr) ** 2 - b_imp ** 2, 0.0)) / (a_r * np.sqrt(max(1 - (b_imp / a_r) ** 2, 1e-12)))
            if 0 < arg <= 1:
                vals["dur_h"], src["dur_h"] = 24.0 * p / np.pi * float(np.arcsin(arg)), "pscomppars T14(P,a/R*,b,Rp/R*)"
        sig = eph_sigma_h(p, getattr(r, "pl_orbpererr1", np.nan), getattr(r, "pl_tranmiderr1", np.nan), vals["t0"]) \
            if src["t0"] == "pscomppars" else np.nan
        ps_targets.append(dict(tic=tic, name=_s(getattr(r, "pl_name", "")), type="ps", disp="CONFIRMED", period=p,
                               period_err=_f(getattr(r, "pl_orbpererr1", np.nan)), **vals,
                               depth_src=src["depth_ppm"], dur_src=src["dur_h"], t0_src=src["t0"], eph_sigma_h=sig))

    toi_targets = []
    for r in toi.itertuples():
        vals = dict(depth_ppm=_f(r.depth_ppm), dur_h=_f(r.dur_h), t0=_f(r.epoch_bjd) - BTJD_OFFSET)
        src = {k: ("exofop_toi" if np.isfinite(v) else "") for k, v in vals.items()}
        sig = eph_sigma_h(_f(r.period), getattr(r, "period_err", np.nan), getattr(r, "epoch_err", np.nan), vals["t0"])
        toi_targets.append(dict(tic=int(r.tic), name=f"TOI {r.toi}", type="toi", disp=_s(r.disp), period=_f(r.period),
                                period_err=_f(getattr(r, "period_err", np.nan)), **vals,
                                depth_src=src["depth_ppm"], dur_src=src["dur_h"], t0_src=src["t0"], eph_sigma_h=sig))

    def same(a, b, rel=0.0):
        return (a["tic"] == b["tic"] and np.isfinite(a["period"]) and np.isfinite(b["period"])
                and abs(a["period"] - b["period"]) <= max(tol(b["period"], b["tic"]), rel * b["period"]))

    final = []
    for pt in ps_targets:
        if not any(same(ex, pt) for ex in final):
            final.append(pt)
    for tt in toi_targets:
        # a TOI is the same planet as a confirmed planet on the star within CATALOG_MATCH_REL (catalogue periods of
        # one planet differ by far more than one grid step); TOIs among themselves use the registered tolerance
        if not any(same(ex, tt, CATALOG_MATCH_REL if ex["type"] == "ps" else 0.0) for ex in final):
            final.append(tt)
    cols = ["tic", "name", "type", "disp", "period", "period_err", "t0", "dur_h", "depth_ppm", "depth_src", "dur_src",
            "t0_src", "eph_sigma_h"]
    return pd.DataFrame(final, columns=cols)


# ---------------------------------------------------------------- expected S and eligibility (computed once)
def derive_expected(rd: Dict[str, Any], cfg: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
    """Threshold-independent part of eligibility plus expected S from the frozen inputs of one target."""
    cfg = cfg if cfg is not None else settings.cfg()
    kp = cfg.get("known_planets", {})
    contam_tol = float(kp.get("contamination_diff_tol", 0.20))
    min_tr = int(kp.get("min_transits", 2))
    s_scale = float(kp.get("expected_s_scale", 1.0))      # amendment 2026-11-12 item 4 (injection obs/exp ratio)
    window = float(kp.get("detrend_window_d", DETREND_WINDOW_DAYS))
    P, dur_h, depth_ppm, t0 = _f(rd.get("period")), _f(rd.get("dur_h")), _f(rd.get("depth_ppm")), _f(rd.get("t0"))
    ap = _f(rd.get("aperture_contamination"))
    cr = _f(rd.get("contratio"))
    err = _s(rd.get("error"))
    n_raw = int(_f(rd.get("n_raw")) if np.isfinite(_f(rd.get("n_raw"))) else 0)
    n_tr = int(_f(rd.get("n_tr")) if np.isfinite(_f(rd.get("n_tr"))) else 0)
    cdpp = _f(rd.get("cdpp_d"))
    flag = bool(np.isfinite(ap) and ap > 0 and np.isfinite(cr) and cr >= 0
                and abs(ap - 1.0 / (1.0 + cr)) / max(ap, 1e-4) > contam_tol)
    reason = ""
    if not (P > 0):
        reason = "missing_period"
    elif not (dur_h > 0):
        reason = "missing_duration"
    elif not (depth_ppm > 0):
        reason = "missing_depth"
    elif not np.isfinite(t0):
        reason = "missing_epoch"
    elif not err and not (ap > 0):
        reason = "missing_aperture_contamination"
    exp_s = exp_raw = np.nan
    if not reason and not err:
        if n_raw < min_tr:
            reason = "transits_under_min"
        elif n_tr < min_tr or window < 3.0 * dur_h / 24.0:
            reason = "in_declared_mask"
        if np.isfinite(cdpp) and cdpp > 0 and n_tr > 0:
            exp_raw = float(depth_ppm * 1e-6 * ap * np.sqrt(n_tr) / cdpp)
            exp_s = exp_raw * s_scale
        elif not reason:
            err = "no finite CDPP at catalogue duration"
    return dict(contamination_flag=flag, n_raw=n_raw, n_tr=n_tr, expected_S_raw=exp_raw, expected_S=exp_s,
                base_ineligible_reason=reason, error=err)


def apply_threshold(df: pd.DataFrame, s_star: float, required_types: Optional[Tuple[str, ...]] = None) -> pd.DataFrame:
    """Eligibility at S* from the frozen expected S: 'eligible' | 'ineligible' | 'error'; required = eligible or error.
    An error row (missing light curve, runtime error) is never ineligible."""
    if s_star is None or not np.isfinite(s_star) or s_star <= 0:
        raise ValueError("s_star must be provided as a positive float")
    out = df.copy()
    elig, reason = [], []
    for r in out.itertuples():
        base, err = _s(r.base_ineligible_reason), _s(r.error)
        if base:
            elig.append("ineligible"); reason.append(base)
        elif err:
            elig.append("error"); reason.append("")
        elif not (_f(r.expected_S) >= s_star):
            elig.append("ineligible"); reason.append("expected_s_below_threshold")
        else:
            elig.append("eligible"); reason.append("")
    out["eligibility"] = elig
    out["ineligible_reason"] = reason
    out["required"] = out.eligibility.isin(["eligible", "error"])
    if required_types is not None:       # amendment 2026-11-12 item 3: TOIs are scored, decided by ExoFOP disposition
        out["required"] = out["required"] & out["type"].map(_s).isin(list(required_types))
    return out


def read_frozen(path: Union[str, Path]) -> pd.DataFrame:
    """Frozen eligibility CSV -> typed DataFrame (round trip of freeze_eligibility; blank strings stay '')."""
    df = pd.read_csv(path, dtype={c: str for c in FROZEN_STR}, keep_default_na=False, na_values={
        c: ["", "nan", "NaN"] for c in FROZEN_FLOAT + FROZEN_INT})
    missing = [c for c in FROZEN_COLUMNS if c not in df.columns]
    if missing:
        raise ValueError(f"frozen eligibility file {Path(path).name} lacks columns {missing}")
    for c in FROZEN_STR:
        df[c] = df[c].map(_s)
    for c in FROZEN_FLOAT:
        df[c] = pd.to_numeric(df[c], errors="raise").astype(float)
    for c in FROZEN_INT:
        df[c] = pd.to_numeric(df[c], errors="raise").astype(np.int64)
    for c in FROZEN_BOOL:
        df[c] = df[c].map(lambda x: str(x).strip().lower() in ("true", "1"))
    return df


def measure_target(rd: Dict[str, Any], path: Optional[Path], cfg: Dict[str, Any]) -> Dict[str, Any]:
    """Light-curve inputs for one target: CROWDSAP, raw and surviving transit counts, CDPP at catalogue duration.
    A missing light curve or any exception is returned as `error`, never swallowed."""
    out = dict(aperture_contamination=np.nan, n_raw=0, n_tr=0, cdpp_d=np.nan, t_span_d=np.nan, error="")
    if path is None or not Path(path).exists():
        out["error"] = "missing light curve"
        return out
    P, t0 = _f(rd.get("period")), _f(rd.get("t0"))
    dur_d = _f(rd.get("dur_h")) / 24.0
    try:
        cs = read_crowdsap_from_fits(path)
        out["aperture_contamination"] = cs if cs is not None else np.nan
        if not (P > 0 and dur_d > 0 and np.isfinite(t0)):
            return out
        t, f, e, q, _ = lightcurve.load_spoc_raw(path)
        good = lightcurve.good_mask(f, e, q)
        med = np.nanmedian(f[good])
        tg, fg, eg = t[good], f[good] / med, e[good] / med
        out["n_raw"] = count_surviving_transits(tg, P, t0, dur_d)
        d = pipeline.prepare(tg, fg, eg, cfg, prewhiten_on=True, gap_trim=0.3)
        tf, ff, _ = d["flat"]
        out["n_tr"] = count_surviving_transits(tf, P, t0, dur_d)
        out["t_span_d"] = float(np.max(tf) - np.min(tf)) if len(tf) else np.nan
        prof = bls.noise_profile(tf, ff, [dur_d])
        out["cdpp_d"] = float(prof[2][float(dur_d)][0])
    except Exception as exc:                                    # recorded as an error row, never as ineligible
        out["error"] = f"{type(exc).__name__}: {exc}"[:300]
    return out


def freeze_eligibility(
    s_star: float,
    out_file: Union[str, Path],
    sector: int = 70,
    commit: str = "",
    config_hash: str = "",
    targets_df: Optional[pd.DataFrame] = None,
    lc_lookup=None,
    contratio_map: Optional[Dict[int, float]] = None,
    cfg: Optional[Dict[str, Any]] = None,
    searched_tics: Optional[set] = None,
) -> Tuple[pd.DataFrame, Dict[str, Any]]:
    """Compute expected S once, pre-run, for every target and write the frozen CSV + sidecar JSON (sha256).

    lc_lookup(tic) -> light-curve path or None (default: the sector's SPOC manifest, local files only)."""
    if s_star is None or not np.isfinite(s_star) or s_star <= 0:
        raise ValueError("s_star must be provided as a positive float")
    cfg = cfg if cfg is not None else settings.cfg()
    out_p = Path(out_file)
    out_p.parent.mkdir(parents=True, exist_ok=True)
    prov = {}
    if targets_df is None:
        ps_df, ps_path = load_pscomppars()
        toi_df, toi_path = load_toi_table()
        targets_df = select_known_planet_targets(sector=sector, searched_tics=searched_tics, ps_df=ps_df, toi_df=toi_df)
        prov = dict(pscomppars_file=ps_path.name, pscomppars_sha256=_sha256(ps_path),
                    toi_file=toi_path.name, toi_sha256=_sha256(toi_path),
                    n_searched_stars=None if searched_tics is None else len(searched_tics))
    if lc_lookup is None:
        man = targets.fetch_manifest(sector).set_index("tic")

        def lc_lookup(tic):
            return batch.lc_path(sector, man.loc[tic].fname, download=False) if tic in man.index else None
    if contratio_map is None:
        contratio_map = {}
        p82 = config.CATALOG_DIR / "tic82_params.parquet"
        if p82.exists():
            d82 = pd.read_parquet(p82)
            if "contratio" in d82.columns:
                contratio_map = {int(r.ID): _f(r.contratio) for r in d82.itertuples()}

    rows = []
    for r in targets_df.to_dict(orient="records"):
        tic = int(r["tic"])
        rd = {c: r.get(c, "") for c in ("tic", "name", "type", "disp", "period", "t0", "dur_h", "depth_ppm",
                                        "depth_src", "dur_src", "t0_src")}
        rd["period_err"] = _f(r.get("period_err", np.nan))
        rd["eph_sigma_h"] = _f(r.get("eph_sigma_h", np.nan))
        rd["tic"] = tic
        rd["contratio"] = contratio_map.get(tic, np.nan)
        rd.update(measure_target(rd, lc_lookup(tic), cfg))
        rd.update(derive_expected(rd, cfg))
        rows.append(rd)
    df = apply_threshold(pd.DataFrame(rows, columns=[c for c in FROZEN_COLUMNS if c not in
                                                     ("eligibility", "ineligible_reason", "required")]), s_star,
                         _required_types(cfg))
    df = df[list(FROZEN_COLUMNS)]
    df.to_csv(out_p, index=False)

    sidecar = dict(
        timestamp_utc=dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
        s_star=float(s_star), sector=int(sector),
        commit=commit or settings.commit(), config_hash=config_hash or settings.config_hash(),
        sha256=_sha256(out_p), total_targets=int(len(df)),
        eligibility_counts=df.eligibility.value_counts().to_dict(),
        ineligible_reasons=df[df.eligibility == "ineligible"].ineligible_reason.value_counts().to_dict(),
        **prov,
    )
    Path(str(out_p) + ".meta.json").write_text(json.dumps(sidecar, indent=2))
    return df, sidecar


# ---------------------------------------------------------------- evaluation
def _required_types(cfg: Dict[str, Any]) -> Optional[Tuple[str, ...]]:
    rt = cfg.get("known_planets", {}).get("required_types")
    return tuple(rt) if rt else None


def _search_index(search_results_df: Optional[pd.DataFrame]) -> Dict[int, Dict[str, Any]]:
    res = {}
    if search_results_df is None or search_results_df.empty:
        return res
    d = search_results_df
    if "mode" in d.columns:
        d = d[d["mode"] == "real"]
    for r in d.to_dict(orient="records"):
        tic = int(r["tic"])
        P = r.get("period")
        if P is None or not np.isfinite(_f(P)):
            P = r.get("P_grid")
        rec = dict(period=_f(P), S=_f(r.get("S")), status=_s(r.get("status")) or "ok", error=_s(r.get("error")),
                   plateaus=r.get("plateaus") if _s(r.get("plateaus")) else None, grid_df=_f(r.get("grid_df")),
                   vetted=(None if "passed_lc_checks" not in r or not np.isfinite(_f(r.get("passed_lc_checks")))
                           else bool(_f(r.get("passed_lc_checks")))))
        if tic in res:      # A7: exactly one real-search row per TIC
            prev = res[tic]
            if not (np.isclose(prev["period"], rec["period"], equal_nan=True) and
                    np.isclose(prev["S"], rec["S"], equal_nan=True) and prev["status"] == rec["status"]):
                raise ValueError(f"Conflicting duplicate real-search rows for TIC {tic}")
        res[tic] = rec
    return res


def evaluate_known_planets(
    targets_df: pd.DataFrame,
    search_results_df: Optional[pd.DataFrame],
    s_star: float,
    durations: np.ndarray = search.DURATIONS,
    cfg: Optional[Dict[str, Any]] = None,
    diagnostics_df: Optional[pd.DataFrame] = None,
    sector: int = 70,
) -> Tuple[pd.DataFrame, Dict[str, Any]]:
    """Evaluate recovery of the frozen targets against the real-search rows (one per TIC).

    targets_df: the frozen eligibility table (read_frozen), or raw target rows, for which the frozen quantities are
    derived here with the same function freeze_eligibility uses. diagnostics_df (optional, keyed by tic + name):
    S_cat, n_tr_cat, depth_ratio_cat at the catalogue ephemeris, used only to name the losing stage.
    B3: S* is passed in; error if not supplied.
    """
    if s_star is None or not np.isfinite(s_star) or s_star <= 0:
        raise ValueError("s_star must be provided as a positive float")
    cfg = cfg if cfg is not None else settings.cfg()
    kp = cfg.get("known_planets", {})
    alias_ratios = tuple(float(x) for x in kp.get("alias_ratios", [2.0, 0.5]))
    tol_mode = kp.get("tolerance", "plateau")            # "bls_peak": amendment 2026-11-12 item 1
    nsig_p = float(kp.get("period_err_nsigma", 3.0))
    alias_rec = bool(kp.get("alias_recovered", False))   # amendment 2026-11-12 item 5

    tdf = targets_df.copy()
    if "period" not in tdf.columns and "P_cat" in tdf.columns:
        tdf = tdf.rename(columns={"P_cat": "period"})
    if "expected_S" not in tdf.columns or "base_ineligible_reason" not in tdf.columns:
        recs = []
        for r in tdf.to_dict(orient="records"):
            r.update(derive_expected(r, cfg))
            recs.append(r)
        tdf = pd.DataFrame(recs)
    for c in ("name", "error", "type", "disp"):
        if c not in tdf.columns:
            tdf[c] = ""
    tdf = apply_threshold(tdf, s_star, _required_types(cfg))

    res_by_tic = _search_index(search_results_df)
    diag = {}
    if diagnostics_df is not None and not diagnostics_df.empty:
        for r in diagnostics_df.to_dict(orient="records"):
            diag[(int(r["tic"]), _s(r.get("name")))] = r

    rows = []
    for r in tdf.to_dict(orient="records"):
        tic, name = int(r["tic"]), _s(r.get("name"))
        P_cat = _f(r.get("period"))
        found = res_by_tic.get(tic)
        err = _s(r.get("error"))
        if found is None:
            search_err = "missing search result"
        elif found["status"] != "ok":
            search_err = f"search error: {found['error'] or found['status']}"
        else:
            search_err = ""
        P_found = found["period"] if found else np.nan
        S_found = found["S"] if found else np.nan
        if tol_mode == "bls_peak":
            hw = bls_peak_tolerance(P_cat, r.get("dur_h"), r.get("t_span_d"), r.get("period_err"), nsig_p)
        else:
            hw = get_plateau_halfwidth(P_cat, plateaus_list=found["plateaus"] if found else None,
                                       grid_df=found["grid_df"] if found else None, tic=tic, sector=sector) \
                if np.isfinite(P_cat) else 0.0
        rec, mt = check_recovery(P_cat, P_found, hw, alias_ratios=alias_ratios, s_found=S_found, s_star=s_star,
                                 alias_scale=tol_mode == "bls_peak", alias_recovered=alias_rec) \
            if not search_err else (False, "none")
        vetted = found["vetted"] if found else None
        recovered = bool(rec and (vetted is None or vetted))
        stage = ""
        if (r["required"] or r["eligibility"] == "eligible") and not recovered:
            dg = diag.get((tic, name), {})
            stage = losing_stage(mt, S_found, s_star, vetted, error=err or search_err, s_cat=dg.get("S_cat"),
                                 n_tr_cat=dg.get("n_tr_cat"), depth_ratio_cat=dg.get("depth_ratio_cat"))
        rows.append(dict(
            tic=tic, toi_or_planet=name, type=_s(r.get("type")), disp=_s(r.get("disp")), P_cat=P_cat,
            period_err=_f(r.get("period_err")), t_span_d=_f(r.get("t_span_d")), expected_S_raw=_f(r.get("expected_S_raw")), depth_ppm=_f(r.get("depth_ppm")),
            dur_h=_f(r.get("dur_h")), t0=_f(r.get("t0")), eph_sigma_h=_f(r.get("eph_sigma_h")), n_tr=int(_f(r.get("n_tr")) or 0) if np.isfinite(_f(r.get("n_tr"))) else 0,
            expected_S=_f(r.get("expected_S")), aperture_contamination=_f(r.get("aperture_contamination")),
            contratio=_f(r.get("contratio")), contamination_flag=bool(r.get("contamination_flag", False)),
            eligibility=r["eligibility"], ineligible_reason=r["ineligible_reason"], error=err or search_err,
            required=bool(r["required"]), P_found=P_found, S_found=S_found, half_width=hw, match_type=mt,
            vetted=vetted, recovered=recovered if r["required"] else bool(rec and (vetted is None or vetted)),
            stage=stage,
        ))
    out_df = pd.DataFrame(rows)
    req = out_df[out_df.required] if len(out_df) else out_df
    missed = req[~req.recovered] if len(req) else req
    summary = dict(
        s_star=float(s_star),
        total_targets=int(len(out_df)),
        eligible_targets=int((out_df.eligibility == "eligible").sum()) if len(out_df) else 0,
        error_targets=int((out_df.eligibility == "error").sum()) if len(out_df) else 0,
        required_targets=int(len(req)),
        eligible_recovered=int(req.recovered.sum()) if len(req) else 0,
        eligible_missed=int(len(missed)),
        required_aliases=int(req.match_type.isin(["2x", "0.5x"]).sum()) if len(req) else 0,
        stage_counts=missed.stage.value_counts().to_dict() if len(missed) else {},
        tolerance=tol_mode, alias_recovered=alias_rec, required_types=list(_required_types(cfg) or []),
        pass_=bool(len(missed) == 0 and len(req) > 0),
    )
    summary["pass"] = summary.pop("pass_")
    return out_df, summary


def _read_table(p: Union[str, Path]) -> pd.DataFrame:
    p = Path(p)
    if p.is_dir():
        from . import run
        return run.load_rows(p)
    return pd.read_parquet(p) if p.suffix == ".parquet" else pd.read_csv(p)


def run_known_planet_recovery(
    s_star: float,
    search_results_path: Optional[Union[str, Path]] = None,
    out_dir: Union[str, Path] = Path("."),
    sector: int = 70,
    frozen_eligibility_file: Optional[Union[str, Path]] = None,
    cfg: Optional[Dict[str, Any]] = None,
    diagnostics_path: Optional[Union[str, Path]] = None,
) -> Tuple[pd.DataFrame, Dict[str, Any]]:
    """Known-planet recovery report from the frozen eligibility file and the real-search rows.

    Refuses to run without the frozen eligibility file or with a sha256 mismatch (B5); a missing search-results
    path is an error, not an empty result."""
    if s_star is None or not np.isfinite(s_star) or s_star <= 0:
        raise ValueError("s_star must be provided as a positive float")
    cfg = cfg if cfg is not None else settings.cfg()
    if not frozen_eligibility_file:
        frozen_eligibility_file = cfg.get("known_planets", {}).get("frozen_eligibility_file")
    if not frozen_eligibility_file or not Path(frozen_eligibility_file).exists():
        raise RuntimeError("Frozen eligibility file required for known-planet recovery evaluation (Spec B5).")
    meta_p = Path(str(frozen_eligibility_file) + ".meta.json")
    if not meta_p.exists():
        raise FileNotFoundError(f"Missing sidecar metadata file: {meta_p.name} (Spec B5)")
    meta = json.loads(meta_p.read_text())
    cur = _sha256(frozen_eligibility_file)
    if cur != meta.get("sha256"):
        raise ValueError(f"Eligibility file sha256 mismatch: {cur} != {meta.get('sha256')}")
    if not search_results_path or not Path(search_results_path).exists():
        raise FileNotFoundError("search results required (real-search rows, CSV / parquet / run output dir)")

    df_frozen = read_frozen(frozen_eligibility_file)
    df_search = _read_table(search_results_path)
    df_diag = _read_table(diagnostics_path) if diagnostics_path else None
    out_df, summary = evaluate_known_planets(df_frozen, df_search, s_star=s_star, cfg=cfg, diagnostics_df=df_diag,
                                             sector=sector)
    summary["frozen_eligibility_sha256"] = cur
    out_p = Path(out_dir)
    out_p.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(out_p / "known_planets_recovery.csv", index=False)
    (out_p / "known_planets_summary.json").write_text(json.dumps(summary, indent=2, default=float))
    return out_df, summary


def main():
    ap = argparse.ArgumentParser(description="v1.1 Known planet / active TOI recovery test (T2)")
    sp = ap.add_subparsers(dest="cmd", required=True)
    fp = sp.add_parser("fetch-pscomppars", help="fetch pscomppars with transit columns (dated cache + sha256)")
    fp.add_argument("--date-tag")
    fe = sp.add_parser("freeze-eligibility")
    fe.add_argument("--s-star", type=float, required=True, help="Detection threshold S*")
    fe.add_argument("--out", required=True, help="Output path for frozen eligibility CSV")
    fe.add_argument("--sector", type=int, default=70, help="TESS Sector")
    fe.add_argument("--searched-tics", help="JSON list of searched TICs (default: the results store for the sector)")
    r = sp.add_parser("report")
    r.add_argument("--s-star", type=float, required=True, help="Detection threshold S*")
    r.add_argument("--results", required=True, help="Real-search rows (CSV, parquet or run output dir)")
    r.add_argument("--frozen-eligibility", help="Path to frozen eligibility CSV")
    r.add_argument("--diagnostics", help="Optional catalogue-ephemeris diagnostics (tic, name, S_cat, ...)")
    r.add_argument("--out", required=True, help="Output directory")
    r.add_argument("--sector", type=int, default=70, help="TESS Sector")
    a = ap.parse_args()
    if a.cmd == "fetch-pscomppars":
        print(fetch_pscomppars(date_tag=a.date_tag))
    elif a.cmd == "freeze-eligibility":
        st = set(int(x) for x in json.loads(Path(a.searched_tics).read_text())) if a.searched_tics else None
        freeze_eligibility(s_star=a.s_star, out_file=a.out, sector=a.sector, searched_tics=st)
    else:
        _, s = run_known_planet_recovery(s_star=a.s_star, search_results_path=a.results, out_dir=a.out,
                                         sector=a.sector, frozen_eligibility_file=a.frozen_eligibility,
                                         diagnostics_path=a.diagnostics)
        print(json.dumps(s, indent=1, default=float))


if __name__ == "__main__":
    main()
