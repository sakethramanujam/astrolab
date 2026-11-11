"""Refuse to start a worker on a dirty checkout (v1.1 rule; disclosure docs/v1.0_disclosure_ties_and_numba_cache.md).

Tracked files must match HEAD (untracked outputs such as site data are allowed). Production entry points call
`require_clean()`; dry runs may pass allow_dirty=True explicitly, and the stamp then records it.
"""
from __future__ import annotations

import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def dirty_files(root=ROOT):
    out = subprocess.run(["git", "status", "--porcelain", "--untracked-files=no"], cwd=root,
                         capture_output=True, text=True, check=True).stdout
    return [l[3:] for l in out.splitlines() if l.strip()]


def require_clean(allow_dirty=False, root=ROOT):
    files = dirty_files(root)
    if files and not allow_dirty:
        raise SystemExit(f"refusing to start: {len(files)} tracked file(s) differ from HEAD (commit first): "
                         + ", ".join(files[:10]))
    return files
