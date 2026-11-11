"""Registration checks for v1.1 amendment (2026-11-11, Section 5).

Asserts that no excluded TIC appears in the paired (intersection) set of flipped
and scrambled null tables for a sector, and reports the paired count.
NOT part of beta_select.py.
"""
from __future__ import annotations

from pathlib import Path
from typing import Dict, Any, Union
import pandas as pd


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
    AssertionError
        If any excluded TIC appears in the paired set.
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

    if sector is not None:
        if "sector" in df_flip.columns:
            df_flip = df_flip[df_flip["sector"] == sector]
        if "sector" in df_scr.columns:
            df_scr = df_scr[df_scr["sector"] == sector]
        if "sector" in df_ex.columns:
            df_ex = df_ex[df_ex["sector"] == sector]

    # Status / quality filter if status column present
    if "status" in df_flip.columns:
        df_flip = df_flip[df_flip["status"] == "ok"]
    if "status" in df_scr.columns:
        df_scr = df_scr[df_scr["status"] == "ok"]

    flip_tics = set(df_flip["tic"].astype(int))
    scr_tics = set(df_scr["tic"].astype(int))
    ex_col = "tic" if "tic" in df_ex.columns else "TIC"
    ex_tics = set(df_ex[ex_col].astype(int))

    paired_tics = flip_tics.intersection(scr_tics)
    overlap = paired_tics.intersection(ex_tics)

    assert len(overlap) == 0, (
        f"Population guard failed: {len(overlap)} excluded TIC(s) appear in the paired "
        f"null set (e.g. {sorted(overlap)[:5]})"
    )

    return {
        "paired_count": len(paired_tics),
        "flipped_count": len(flip_tics),
        "scrambled_count": len(scr_tics),
        "excluded_count": len(ex_tics),
        "excluded_in_paired": 0,
        "status": "ok",
    }
