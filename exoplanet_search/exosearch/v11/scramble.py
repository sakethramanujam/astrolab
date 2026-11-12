"""v1.1 scrambled-null transform: per-segment cyclic rotation rule.

REVISED v1.1 SCRAMBLE RULE (Scientist ruling 2026-11-11, supersedes the fallback):
(1) Input is the quality-masked light curve; apply lightcurve.gap_edge_mask
    (gap=0.5 d, trim=0.3 d) BEFORE scrambling and drop those cadences.
(2) Split into contiguous segments at gaps > 0.5 d.
(3) Rotate EVERY segment, regardless of length, by its own independent random
    cyclic offset from a numpy default_rng seeded by the task seed; flux and flux_err
    move together, time stamps fixed, no sample crosses a gap.
(4) There is NO block-shuffle fallback and NO leave-unrotated path: both deleted.
    A 1-sample segment has no non-zero rotation and stays as is (the only unrotated case);
    for longer segments the offset is drawn from 1..n_seg-1, i.e. a zero offset is redrawn.
(5) Return per-segment offsets so the pipeline records them (scr_offsets,
    scr_scheme='segment_rotation').
(6) Parameters in configs/v1.1.json 'scramble' section.
"""
from __future__ import annotations

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
    """Apply the v1.1 scrambled-null transform (segment rotation rule).

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

    # Rule (1): Quality-masked input; apply gap_edge_mask BEFORE scrambling and drop cadences
    keep = lightcurve.gap_edge_mask(t, gap=gap_d, trim=gap_trim_d)
    t = t[keep]
    f = f[keep]
    e = e[keep]

    if t.size == 0:
        return t, f, e, {"scr_offsets": [], "scr_scheme": "segment_rotation", "n_segments": 0}

    # Rule (2): Split into contiguous segments at gaps > 0.5 d
    dt = np.diff(t)
    split_indices = np.where(dt > gap_d)[0] + 1
    seg_indices = np.split(np.arange(t.size), split_indices)

    rng = np.random.default_rng(seed)

    f_out = np.zeros_like(f)
    e_out = np.zeros_like(e)

    offsets: List[int] = []

    # Rule (3): Rotate EVERY segment, regardless of length, by its own independent random cyclic offset.
    # Rule (4): a 1-sample segment cannot rotate (offset 0, documented); otherwise offsets are 1..n_seg-1.
    for seg_idx in seg_indices:
        n_seg = len(seg_idx)
        if n_seg == 0:
            continue
        offset = int(rng.integers(1, n_seg)) if n_seg > 1 else 0      # zero is redrawn: offsets are in 1..n_seg-1
        f_out[seg_idx] = np.roll(f[seg_idx], offset)
        e_out[seg_idx] = np.roll(e[seg_idx], offset)
        offsets.append(offset)

    info = {
        "scr_offsets": offsets,
        "scr_scheme": "segment_rotation",
        "n_segments": len(seg_indices),
    }

    return t, f_out, e_out, info
