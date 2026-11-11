"""Helpers for the reliability batches (injection, null tests, labelled-set checks).

* worker claims: write N into $CLAIM_DIR/<pid>.workers; scripts/continuous.sh shrinks the runner so the
  runner plus batches stay under COMBINED_CAP (8) workers
* light-curve paths that survive the runner's FITS pruning (keep/ hard links, else re-download)
* resumable process-pool map that appends parquet parts
"""
from __future__ import annotations

import contextlib
import os
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import pandas as pd

from . import config, targets

CLAIM_DIR = Path(os.environ.get("CLAIM_DIR", "/tmp/exosearch-claims"))
MAST_FILE = "https://mast.stsci.edu/api/v0.1/Download/file/?uri=mast:TESS/product/"


@contextlib.contextmanager
def claim(n):
    CLAIM_DIR.mkdir(parents=True, exist_ok=True)
    f = CLAIM_DIR / f"{os.getpid()}.workers"
    f.write_text(str(int(n)))
    try:
        yield
    finally:
        f.unlink(missing_ok=True)


def keep_dir(sector):
    return config.DATA_ROOT / "keep" / f"s{int(sector):04d}" / "lc"


def lc_path(sector, fname, download=True):
    for d in (config.lc_dir(sector), keep_dir(sector)):
        p = d / fname
        if p.exists() and p.stat().st_size > 0:
            return p
    if not download:
        return None
    return targets.download_one(fname, MAST_FILE + fname, keep_dir(sector))


def pmap(fn, tasks, out_dir, workers, key="key", task_key="key", part_size=500, log=print):
    """Run fn(task) -> dict|list[dict] over tasks (dicts with a unique `task_key`), resumably; result rows
    carry the task's id in column `key`.
    Results are appended as out_dir/part-XXXXX.parquet; tasks whose key is already present are skipped."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    done = set()
    parts = sorted(out_dir.glob("part-*.parquet"))
    for p in parts:
        done |= set(pd.read_parquet(p, columns=[key])[key])
    todo = [t for t in tasks if t[task_key] not in done]
    log(f"{len(tasks)} tasks, {len(done)} done, {len(todo)} to run on {workers} workers")
    buf, n, t0 = [], len(parts), time.time()
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futs = [ex.submit(fn, t) for t in todo]
        for i, fu in enumerate(as_completed(futs), 1):
            r = fu.result()
            buf.extend(r if isinstance(r, list) else [r])
            if len(buf) >= part_size or i == len(futs):
                pd.DataFrame(buf).to_parquet(out_dir / f"part-{n:05d}.parquet", index=False)
                n += 1
                buf = []
                el = time.time() - t0
                log(f"  {i}/{len(todo)}  {el / i * workers:.2f} s/task/worker  eta {(len(todo) - i) * el / i / 60:.0f} min")
    return pd.concat([pd.read_parquet(p) for p in sorted(out_dir.glob("part-*.parquet"))], ignore_index=True) \
        if any(out_dir.glob("part-*.parquet")) else pd.DataFrame()
