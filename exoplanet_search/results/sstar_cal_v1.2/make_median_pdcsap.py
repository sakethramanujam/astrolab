"""Build median_pdcsap_s59-64.csv for the 2026-11-24 brightness cut (registration_v1.2.md; ea3fa7c).

Per S59-S64 manifest job: median of PDCSAP_FLUX over cadences with QUALITY == 0 and finite flux, read
read-only from the job's local SPOC 2-min file (exosearch.batch.lc_path(..., download=False)).

  python make_median_pdcsap.py MANIFEST.jsonl.gz OUT.csv [WORKERS]
"""
import csv
import gzip
import json
import sys
from concurrent.futures import ProcessPoolExecutor

import numpy as np

SECTORS = (59, 60, 61, 62, 63, 64)


def one(t):
    s, tic, fname, path = t
    if path is None:
        return (s, tic, fname, 0, 0, "", "lc_not_found")
    try:
        from astropy.io import fits
        with fits.open(path, memmap=True, lazy_load_hdus=True) as h:
            d = h[1].data
            f = np.asarray(d["PDCSAP_FLUX"], dtype=float)
            q = np.asarray(d["QUALITY"])
        g = (q == 0) & np.isfinite(f)
        n = int(g.sum())
        med = float(np.median(f[g])) if n else float("nan")
        return (s, tic, fname, 1, n, repr(med), "" if n else "no_good_cadences")
    except Exception as e:
        return (s, tic, fname, 1, 0, "", f"read_error:{type(e).__name__}:{e}")


def main():
    from exosearch import batch
    man, out = sys.argv[1], sys.argv[2]
    workers = int(sys.argv[3]) if len(sys.argv) > 3 else 8
    with gzip.open(man, "rt") as fh:
        fh.readline()
        jobs = [json.loads(line) for line in fh if line.strip()]
    tasks = []
    for j in jobs:
        s, tic, fn = int(j["sector"]), int(j["tic"]), j.get("fname")
        if s not in SECTORS:
            continue
        p = batch.lc_path(s, fn, download=False) if fn else None
        tasks.append((s, tic, fn, str(p) if p else None))
    with ProcessPoolExecutor(max_workers=workers) as ex:
        rows = sorted(ex.map(one, tasks, chunksize=64))
    with open(out, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["sector", "tic", "fname", "lc_found", "n_good", "median_pdcsap", "error"])
        w.writerows(rows)


if __name__ == "__main__":
    main()
