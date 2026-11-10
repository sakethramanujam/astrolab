"""Sector queue: work through every TESS 2-min sector over time, resuming where it left off.

Progress is derived from the results DB (stars with status ok vs. the sector manifest), so the
queue never drifts from what was actually searched.  Order: unfinished sectors already started
first, then newest sectors first (least mined by others), then older ones.

    python -m exosearch.sectorqueue status
    python -m exosearch.sectorqueue next [--count 3]
"""
from __future__ import annotations

import argparse
import json
import re
import time

import requests

from . import config, store, targets

INDEX_URL = "https://archive.stsci.edu/missions/tess/download_scripts/sector/"
_RX = re.compile(r"tesscurl_sector_(\d+)_lc\.sh")


def available_sectors(max_age_days=7):
    cache = config.DATA_ROOT / "sectors_available.json"
    if cache.exists() and time.time() - cache.stat().st_mtime < max_age_days * 86400:
        return json.loads(cache.read_text())
    html = requests.get(INDEX_URL, timeout=60).text
    sectors = sorted({int(m) for m in _RX.findall(html) if 0 < int(m) < 500})
    config.DATA_ROOT.mkdir(parents=True, exist_ok=True)
    cache.write_text(json.dumps(sectors))
    return sectors


def manifest_size(sector):
    cache = config.sector_dir(sector) / f"tesscurl_sector_{sector}_lc.sh"
    if not cache.exists():
        return None
    return len(targets.fetch_manifest(sector))


def status():
    con = store.connect()
    done = dict(con.execute("SELECT sector, COUNT(*) FROM results WHERE status='ok' GROUP BY sector"))
    errs = dict(con.execute("SELECT sector, COUNT(*) FROM results WHERE status!='ok' GROUP BY sector"))
    out = []
    for s in available_sectors():
        n = manifest_size(s)
        d = done.get(s, 0)
        e = errs.get(s, 0)
        if n and d + e >= n:
            st = "done"
        elif d:
            st = "in_progress"
        else:
            st = "pending"
        out.append(dict(sector=s, n_targets=n, n_searched=d, n_errors=e, status=st))
    return out


def next_sectors(count=1):
    st = status()
    started = [r["sector"] for r in st if r["status"] == "in_progress"]
    pending = sorted((r["sector"] for r in st if r["status"] == "pending"), reverse=True)
    return (started + pending)[:count]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["status", "next"])
    ap.add_argument("--count", type=int, default=1)
    a = ap.parse_args()
    if a.cmd == "next":
        print(" ".join(map(str, next_sectors(a.count))))
    else:
        for r in status():
            if r["status"] != "pending":
                print(r)
        print(sum(r["status"] == "pending" for r in status()), "sectors pending")


if __name__ == "__main__":
    main()
