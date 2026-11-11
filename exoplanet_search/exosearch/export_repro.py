"""Export the raw per-row result tables as gzipped CSV for a clean-room reproduction from SCHEMA.md.

  python -m exosearch.export_repro          # -> data/reproduction/*.csv.gz + MANIFEST.json (in the repo)

Snapshot of products that are still growing (lanes running); MANIFEST.json records row counts, sha256 and the
pipeline_version / config_hash / git_hash values present. Columns follow SCHEMA.md; anything that could carry a
host-specific path (light-curve paths, error strings, profiling) is dropped and every string cell is checked.
"""
from __future__ import annotations

import datetime as dt
import gzip
import hashlib
import io
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd

from . import interim, nulltest

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "reproduction"
DROP = {"error", "path", "runtime_s", "load_s", "detrend_s", "search_s", "job_path", "fname"}
BAD = re.compile(r"(/home/|/media/|/Users/|\\\\|[A-Za-z]:\\\\|192\.168\.|(^|[^0-9.])10\.\d+\.\d+\.\d+|(^|[^0-9.])100\.\d+\.\d+\.\d+)")


def clean(d):
    d = d[[c for c in d.columns if c not in DROP]].copy()
    for c in d.columns:
        if d[c].dtype == object:
            s = d[c].dropna().astype(str)
            hit = s[s.str.contains(BAD)]
            if len(hit):
                raise ValueError(f"column {c} carries a path/host-like value, e.g. {hit.iloc[0][:40]!r}")
    return d


def write(d, name):
    buf = io.StringIO()
    d.to_csv(buf, index=False)
    raw = buf.getvalue().encode()
    p = OUT / name
    with open(p, "wb") as fh, gzip.GzipFile(fileobj=fh, mode="wb", mtime=0, filename="") as gz:   # deterministic gzip
        gz.write(raw)
    return dict(file=name, rows=int(len(d)), columns=list(d.columns), bytes_gz=p.stat().st_size,
                sha256_gz=hashlib.sha256(p.read_bytes()).hexdigest(),
                pipeline_version=sorted(map(str, d.pipeline_version.dropna().unique())) if "pipeline_version" in d else None,
                config_hash=sorted(map(str, d.config_hash.dropna().unique())) if "config_hash" in d else None,
                git_hash=sorted(map(str, d.git_hash.dropna().unique())) if "git_hash" in d else None)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    P = interim.PROD
    files = []
    fin = P / "nulltests" / "invert" / "null_invert_maxS.parquet"
    null = pd.read_parquet(fin) if fin.exists() else nulltest.finalize(interim.read_parts(P / "nulltests" / "invert" / "parts" / "part-*.parquet"))
    null["S_exact"] = null.depth * np.sqrt(null.n_tr) / null.cdpp_d          # interim-layer derived (depth audit)
    files.append(write(clean(null.sort_values(["sector", "tic"])), "null_invert_maxS.csv.gz"))
    inj = interim.read_parts(P / "injections" / "*" / "part-*.parquet")
    inj["depth_box_implied"] = inj.S * inj.cdpp_d_rec / np.sqrt(inj.n_tr_rec)   # not stored by v1.0.0 injections
    inj["S_exact"] = inj.depth_rec * np.sqrt(inj.n_tr_rec) / inj.cdpp_d_rec
    files.append(write(clean(inj.sort_values(["sector", "run", "tic"])), "injections_v1.csv.gz"))
    q = interim.read_parts(P / "analysis_layer" / "star_quality_s*_parts" / "part-*.parquet")
    files.append(write(clean(q.sort_values(["sector", "tic"])), "star_quality_v1.csv.gz"))
    from . import sigma_p
    sp = sigma_p.load()
    if sp is not None:                          # plateau half-width at the stored P of each candidate-table row
        files.append(write(clean(sp.sort_values(["sector", "tic"])), "sigma_p_v1.csv.gz"))
    man = dict(snapshot=dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%MZ"),
               note="snapshot of products still being extended by running batches; re-export later for the full sets",
               schema="SCHEMA.md", configs=["configs/star_quality_v1.json", "configs/interim_threshold_v1.json",
                                            "configs/injection_v1.json", "configs/detection_threshold.json"],
               derived_columns={"null_invert_maxS.csv.gz": ["S_exact"],
                                "injections_v1.csv.gz": ["depth_box_implied", "S_exact"],
                                "sigma_p_v1.csv.gz": ["sigma_P", "sigma_P_in_band", "plateau_width", "plateau_centre"]},
               sigma_P="half-width (days) of the in-band plateau containing the stored v1.0.0 period, from the interim "
                       "fixed-grid search (docs/preregistration/tie_break.md amendment 1); half a grid step when the "
                       "stored period is not in the tie band (sigma_P_in_band false)",
               depth_audit="S uses depth_box (BLS fast-method depth); depth is the compute_stats box depth. "
                           "S_exact = depth*sqrt(n_tr)/cdpp_d is the interim-layer statistic; see SCHEMA.md",
               files=files)
    (OUT / "MANIFEST.json").write_text(json.dumps(man, indent=1))
    print(json.dumps([{k: f[k] for k in ("file", "rows", "bytes_gz")} for f in files]))


if __name__ == "__main__":
    main()
