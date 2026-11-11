"""v1.1 scrambled-null transform: per-segment cyclic rotation with fallback.

REGISTERED v1.1 SCRAMBLE RULE (2026-11-11, room decision):
(1) Input is the quality-masked light curve; apply lightcurve.gap_edge_mask
    (same gap=0.5 d, trim as the detrender uses) BEFORE scrambling and drop those cadences.
(2) Split into contiguous segments at gaps > 0.5 d.
(3) Rotate each segment's flux and flux_err (together) by its own independent random
    cyclic offset drawn from a numpy default_rng seeded by the task seed; time stamps
    stay fixed; samples never cross a segment/gap boundary.
(4) Fallback only if a segment is too short for a meaningful rotation (shorter than
    3x the longest trial duration in the v1.1 config grid): shuffle within that segment
    in contiguous blocks of length >= 3x the longest trial duration (if the segment is
    shorter than one block, leave it unrotated and record scheme='unrotated').
(5) Return the per-segment offsets/scheme so the pipeline can record them in the output row
    (e.g. scr_offsets, scr_scheme).
(6) Put any parameters in configs/v1.1.json under a new 'scramble' section.
"""
from __future__ import annotations

import json
from typing import Any, Dict, List, Tuple

import numpy as np

from .. import lightcurve


def segment_rotate(
    t: np.ndarray,
    f: np.ndarray,
    e: np.ndarray,
    seed: int,
    cfg: Dict[str, Any],
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, Dict[str, Any]]:
    """Apply the v1.1 scrambled-null transform.

    Parameters
    ----------
    t : np.ndarray
        Time stamps (quality-masked).
    f : np.ndarray
        Flux values.
    e : np.ndarray
        Flux error values.
    seed : int
        Seed for the random number generator.
    cfg : dict
        Configuration dictionary (contains 'scramble' section).

    Returns
    -------
    t_out : np.ndarray
        Time stamps after gap-edge cadences are dropped.
    f_out : np.ndarray
        Scrambled flux array.
    e_out : np.ndarray
        Scrambled flux error array.
    info : dict
        Dictionary containing 'scr_offsets' and 'scr_scheme' (and segment metadata).
    """
    scr_cfg = cfg.get("scramble", {})
    gap_d = scr_cfg.get("gap_d", 0.5)
    gap_trim_d = scr_cfg.get("gap_trim_d", 0.3)
    block_mult = scr_cfg.get("block_duration_multiplier", 3.0)
    max_d = scr_cfg.get("longest_trial_duration_d", 0.25)
    min_seg_d = block_mult * max_d

    # Rule (1): Quality-masked input; apply gap_edge_mask BEFORE scrambling and drop cadences
    keep = lightcurve.gap_edge_mask(t, gap=gap_d, trim=gap_trim_d)
    t = t[keep]
    f = f[keep]
    e = e[keep]

    if t.size == 0:
        return t, f, e, {"scr_offsets": [], "scr_scheme": "empty", "n_segments": 0}

    # Rule (2): Split into contiguous segments at gaps > 0.5 d
    dt = np.diff(t)
    split_indices = np.where(dt > gap_d)[0] + 1
    seg_indices = np.split(np.arange(t.size), split_indices)

    rng = np.random.default_rng(seed)

    f_out = np.zeros_like(f)
    e_out = np.zeros_like(e)

    schemes: List[str] = []
    offsets: List[int] = []

    for seg_idx in seg_indices:
        n_seg = len(seg_idx)
        if n_seg == 0:
            continue
        seg_t = t[seg_idx]
        seg_span = float(seg_t[-1] - seg_t[0])
        seg_f = f[seg_idx]
        seg_e = e[seg_idx]

        cad = float(np.median(np.diff(seg_t))) if n_seg > 1 else 0.0
        n_block = max(int(np.ceil(min_seg_d / cad)), 1) if cad > 0 else 1

        # Check rotation condition vs fallback
        if seg_span >= min_seg_d:
            # Rule (3): Independent random cyclic offset
            # offset in [0, n_seg)
            offset = int(rng.integers(0, n_seg))
            f_rot = np.roll(seg_f, offset)
            e_rot = np.roll(seg_e, offset)
            f_out[seg_idx] = f_rot
            e_out[seg_idx] = e_rot
            schemes.append("rotate")
            offsets.append(offset)
        elif n_seg >= n_block and n_seg > 1:
            # Rule (4): Fallback for short segment: shuffle contiguous blocks of length >= min_seg_d
            blocks = [np.arange(i, min(i + n_block, n_seg)) for i in range(0, n_seg, n_block)]
            order = rng.permutation(len(blocks))
            idx_perm = np.concatenate([blocks[i] for i in order])
            f_out[seg_idx] = seg_f[idx_perm]
            e_out[seg_idx] = seg_e[idx_perm]
            schemes.append("block_shuffle")
            offsets.append(-1)
        else:
            # Rule (4): Shorter than one block, leave unrotated
            f_out[seg_idx] = seg_f
            e_out[seg_idx] = seg_e
            schemes.append("unrotated")
            offsets.append(0)

    # Summarize scheme across segments
    unique_schemes = sorted(set(schemes))
    if len(unique_schemes) == 1:
        overall_scheme = unique_schemes[0]
    else:
        overall_scheme = "+".join(unique_schemes)

    info = {
        "scr_offsets": offsets,
        "scr_scheme": overall_scheme,
        "scr_schemes": schemes,
        "n_segments": len(seg_indices),
    }

    return t, f_out, e_out, info
