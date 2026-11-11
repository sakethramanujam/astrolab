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
from typing import Dict, Any, List, Optional, Tuple

import numpy as np
import pandas as pd

from .. import batch, bls, config, lightcurve, store, targets
from . import grid, search, settings

# Contamination factor documentation:
# For PDCSAP flux, the SPOC pipeline has already performed crowding correction by subtracting
# estimated excess flux from nearby stars based on the PRF model. Therefore, the observed
# transit depth in PDCSAP is already on an uncontaminated relative flux scale, so the
# depth contamination scaling factor for PDCSAP is 1.0.
# If uncorrected SAP flux were used instead, the observed transit depth would be diluted by
# non-target starlight in the aperture by a factor equal to CROWDSAP (the target flux fraction),
# so the expected depth in SAP would be catalogue_depth * CROWDSAP.
CONTAMINATION_FACTOR_PDCSAP = 1.0


def calculate_grid_halfwidth(
    period: float,
    baseline: float = 27.0,
    durations: np.ndarray = search.DURATIONS,
    min_transits: int = 2,
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
        # Check alias within half-width of the alias period or within tolerance
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
    detrend_window: float = 0.75,
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
            P = getattr(r, "period", None) or getattr(r, "P_grid", None)
            S = getattr(r, "S", None)
            res_by_tic[tic] = dict(period=P, S=S)

    rows = []
    baseline = 27.0

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
