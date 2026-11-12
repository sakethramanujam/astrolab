"""Registration checks for v1.1 amendment (2026-11-11, Section 5).

Asserts that no excluded TIC appears in the paired (intersection) set of flipped
and scrambled null tables for a sector, and reports the paired count.
NOT part of beta_select.py.
"""
from __future__ import annotations

from pathlib import Path
from typing import Dict, Any, Union
import pandas as pd


class PopulationGuardError(RuntimeError):
    """The flipped and scrambled null populations are not the registered same population. Registration stops."""


def check_registration_population(
    flipped_table: Union[pd.DataFrame, str, Path],
    scrambled_table: Union[pd.DataFrame, str, Path],
    exclusion_file: Union[pd.DataFrame, str, Path],
    sector: int | None = None,
) -> Dict[str, Any]:
    """Given flipped and scrambled null tables and the exclusion CSV:
    asserts that no excluded TIC appears in the paired (intersection) set
    and returns the paired count and summary stats.

    Parameters
    ----------
    flipped_table : pd.DataFrame or path
        Results dataframe or path for the flipped (inverted) null.
    scrambled_table : pd.DataFrame or path
        Results dataframe or path for the scrambled null.
    exclusion_file : pd.DataFrame or path
        Exclusion dataframe or path (configs/scramble_exclude_s0070_2026-11-11.csv).
    sector : int, optional
        Sector to filter if present in dataframes.

    Returns
    -------
    dict
        {
            "paired_count": int,
            "flipped_count": int,
            "scrambled_count": int,
            "excluded_count": int,
            "excluded_in_paired": 0,
            "status": "ok",
        }

    Raises
    ------
    PopulationGuardError
        If any excluded TIC appears in the paired set, the exclusion set is empty (also after the sector
        filter), has no tic column, or the paired set is empty. Never a bare assert (survives python -O).
    """
    if isinstance(flipped_table, (str, Path)):
        df_flip = pd.read_csv(flipped_table)
    else:
        df_flip = flipped_table.copy()

    if isinstance(scrambled_table, (str, Path)):
        df_scr = pd.read_csv(scrambled_table)
    else:
        df_scr = scrambled_table.copy()

    if isinstance(exclusion_file, (str, Path)):
        df_ex = pd.read_csv(exclusion_file)
    else:
        df_ex = exclusion_file.copy()

    ex_col = "tic" if "tic" in df_ex.columns else ("TIC" if "TIC" in df_ex.columns else None)
    if ex_col is None:
        raise PopulationGuardError("exclusion table has no tic column")
    if len(df_ex) == 0:
        raise PopulationGuardError("exclusion set is empty: the guard would pass trivially")

    if sector is not None:
        if "sector" in df_flip.columns:
            df_flip = df_flip[df_flip["sector"] == sector]
        if "sector" in df_scr.columns:
            df_scr = df_scr[df_scr["sector"] == sector]
        if "sector" in df_ex.columns:
            df_ex = df_ex[df_ex["sector"] == sector]
            if len(df_ex) == 0:
                raise PopulationGuardError(f"sector filter {sector} empties the exclusion set")

    # Status / quality filter if status column present
    if "status" in df_flip.columns:
        df_flip = df_flip[df_flip["status"] == "ok"]
    if "status" in df_scr.columns:
        df_scr = df_scr[df_scr["status"] == "ok"]

    flip_tics = set(df_flip["tic"].astype(int))
    scr_tics = set(df_scr["tic"].astype(int))
    ex_tics = set(df_ex[ex_col].astype(int))

    paired_tics = flip_tics.intersection(scr_tics)
    if not paired_tics:
        raise PopulationGuardError("paired (flipped and scrambled) set is empty")
    overlap = paired_tics.intersection(ex_tics)

    if overlap:
        raise PopulationGuardError(
            f"Population guard failed: {len(overlap)} excluded TIC(s) appear in the paired "
            f"null set (e.g. {sorted(overlap)[:5]})")

    return {
        "paired_count": len(paired_tics),
        "flipped_count": len(flip_tics),
        "scrambled_count": len(scr_tics),
        "excluded_count": len(ex_tics),
        "excluded_in_paired": 0,
        "status": "ok",
    }


def guard_registered_nulls(flip_tag="invert", scr_tag="scramble_r0", sectors=(70, 106)):
    """Load the same null tables beta_select loads and run the guard per sector with the registered exclusion CSV.

    Sectors with no registered exclusion file are not exempt silently: the guard is only run for sectors that
    have one and the others are reported as 'no_exclusion_registered' (the paired set is not constrained).
    """
    from .. import interim
    from . import settings
    flip, scr = interim.load_null(flip_tag), interim.load_null(scr_tag)
    files = settings.exclusion_files()
    out = {}
    for sec in sectors:
        rel = files.get(str(sec), files.get("*"))
        if rel is None:
            out[sec] = dict(status="no_exclusion_registered")
            continue
        ex = pd.read_csv(settings.ROOT / rel)
        out[sec] = check_registration_population(flip, scr, ex, sector=sec)
    return out


if __name__ == "__main__":
    import argparse
    import json
    ap = argparse.ArgumentParser(description="Same-population guard (run before beta_select)")
    ap.add_argument("--flip", default="invert")
    ap.add_argument("--scr", default="scramble_r0")
    a = ap.parse_args()
    try:
        print(json.dumps(guard_registered_nulls(a.flip, a.scr), indent=1, default=int))
    except PopulationGuardError as e:
        raise SystemExit(f"POPULATION GUARD FAILED: {e}")
