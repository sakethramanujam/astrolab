"""Sector manifests (MAST bulk-download scripts) and resumable, idempotent downloads."""
from __future__ import annotations

import os
import re
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import pandas as pd
import requests

from . import config

MANIFEST_URL = ("https://archive.stsci.edu/missions/tess/download_scripts/sector/"
                "tesscurl_sector_{sector}_lc.sh")
_LINE = re.compile(r"-o\s+(\S+_lc\.fits)\s+(\S+)")


def fetch_manifest(sector: int, refresh: bool = False) -> pd.DataFrame:
    """All SPOC 2-min light curves in a sector -> DataFrame[tic, fname, url]. Cached."""
    sdir = config.sector_dir(sector)
    sdir.mkdir(parents=True, exist_ok=True)
    cache = sdir / f"tesscurl_sector_{sector}_lc.sh"
    if refresh or not cache.exists() or cache.stat().st_size == 0:
        r = requests.get(MANIFEST_URL.format(sector=sector), timeout=120)
        r.raise_for_status()
        tmp = cache.with_suffix(".part")
        tmp.write_text(r.text)
        os.replace(tmp, cache)
    rows = []
    for line in cache.read_text().splitlines():
        m = _LINE.search(line)
        if m:
            fname, url = m.groups()
            rows.append((int(fname.split("-")[2]), fname, url))
    return pd.DataFrame(rows, columns=["tic", "fname", "url"]).drop_duplicates("tic")


def select_targets(manifest: pd.DataFrame, n: int, seed: int = 42) -> pd.DataFrame:
    """Deterministic random subset, so re-runs pick the same stars."""
    n = min(n, len(manifest))
    return manifest.sample(n=n, random_state=seed).sort_values("tic").reset_index(drop=True)


def _valid_fits(path: Path) -> bool:
    try:
        size = path.stat().st_size
        if size == 0 or size % 2880:
            return False
        with open(path, "rb") as fh:
            return fh.read(6) == b"SIMPLE"
    except OSError:
        return False


def download_one(fname: str, url: str, dest_dir: Path, retries: int = 4) -> Path:
    """Download if missing/corrupt. Writes .part then renames, so it is safe to interrupt."""
    dest_dir.mkdir(parents=True, exist_ok=True)
    path = dest_dir / fname
    if _valid_fits(path):
        return path
    tmp = path.with_name(path.name + ".part")
    last = None
    for attempt in range(retries):
        try:
            with requests.get(url, stream=True, timeout=120) as r:
                r.raise_for_status()
                with open(tmp, "wb") as fh:
                    for chunk in r.iter_content(1 << 20):
                        fh.write(chunk)
            if not _valid_fits(tmp):
                raise IOError(f"invalid FITS from {url}")
            os.replace(tmp, path)
            return path
        except Exception as exc:  # network flake -> back off and retry
            last = exc
            time.sleep(2 ** attempt)
    tmp.unlink(missing_ok=True)
    raise RuntimeError(f"download failed for {fname}: {last}")


def download_many(targets: pd.DataFrame, sector: int, threads: int = 8, log=print) -> dict:
    dest = config.lc_dir(sector)
    ok, failed = {}, {}
    with ThreadPoolExecutor(max_workers=threads) as pool:
        futs = {pool.submit(download_one, r.fname, r.url, dest): r.tic for r in targets.itertuples()}
        for i, fut in enumerate(as_completed(futs), 1):
            tic = futs[fut]
            try:
                ok[tic] = fut.result()
            except Exception as exc:
                failed[tic] = str(exc)
            if i % 100 == 0 or i == len(futs):
                log(f"  downloads: {i}/{len(futs)} ({len(failed)} failed)")
    return {"ok": ok, "failed": failed}
