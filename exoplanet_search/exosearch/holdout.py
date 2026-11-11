"""Held-out stars: a seeded, rule-based random 20% of TIC IDs.

A TIC is held out when sha256("<SEED>:<tic>") / 2**256 < FRACTION. The rule extends to any future
sector without re-drawing, and `holdout_tics.txt` (committed) lists the held-out IDs among every star
searched so far. Held-out stars are excluded from all tuning, threshold choices and cut choices;
their metrics are reported separately per sector group.
"""
from __future__ import annotations

import hashlib
from pathlib import Path

SEED = 20261110
FRACTION = 0.20
FILE = Path(__file__).resolve().parents[1] / "holdout_tics.txt"


def is_holdout(tic) -> bool:
    h = hashlib.sha256(f"{SEED}:{int(tic)}".encode()).digest()
    return int.from_bytes(h, "big") / 2 ** 256 < FRACTION


def split(tics):
    """-> (tuning, holdout) lists."""
    tics = [int(t) for t in tics]
    hold = [t for t in tics if is_holdout(t)]
    hs = set(hold)
    return [t for t in tics if t not in hs], hold


def write(tics):
    hold = sorted({int(t) for t in tics if is_holdout(t)})
    head = [f"# held-out TIC IDs (never used for tuning). rule: sha256('{SEED}:<tic>')/2**256 < {FRACTION}",
            f"# seed={SEED} fraction={FRACTION} n_searched={len(set(int(t) for t in tics))} n_holdout={len(hold)}"]
    FILE.write_text("\n".join(head + [str(t) for t in hold]) + "\n")
    return hold
