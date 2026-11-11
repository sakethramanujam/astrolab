"""Star-quality table for the pre-registered cut in configs/star_quality_v1.json (analysis layer over v1.0.0).

For every searched star: sigma_w, CDPP_D and beta_D at each trial duration, measured on the star's own
quality-masked PDCSAP light curve after the frozen detrend. `drop` = beta_max > threshold.
  python -m exosearch.quality --sectors 70,106 --workers 2
Output: products/v<ver>/analysis_layer/star_quality_s<sector>.parquet
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from . import batch, bls, config, holdout, lightcurve, store, targets, task
from .version import stamp

CFG = Path(__file__).resolve().parents[1] / "configs" / "star_quality_v1.json"


def rule():
    return json.loads(CFG.read_text())


def out_dir():
    return config.DATA_ROOT / "products" / f"v{stamp()['pipeline_version']}" / "analysis_layer"


def run_one(job):
    row = dict(key=job["key"], tic=job["tic"], sector=job["sector"], holdout=holdout.is_holdout(job["tic"]))
    try:
        t, f, e, meta = lightcurve.load_spoc(job["path"])
        d = task.prepare_arrays(t, f, e, meta)
        tf, ff, _ = d["flat"]
        sw, cad, prof = bls.noise_profile(tf, ff)
        betas = [v[1] for v in prof.values() if v[1] == v[1]]
        row.update(tmag=meta.get("tmag"), camera=meta.get("camera"), ccd=meta.get("ccd"), sigma_w=sw,
                   beta_max=float(max(betas)) if betas else np.nan, n_cad=int(tf.size), status="ok")
        for k, (cd, b) in prof.items():
            row[f"cdpp_{k:.2f}"] = cd
            row[f"beta_{k:.2f}"] = b
    except Exception as exc:
        row.update(status="error", error=f"{type(exc).__name__}: {exc}"[:300])
    return row


def load(sectors=(70, 106)):
    """Quality table with the cut applied (`drop`); stars that failed to load are dropped too."""
    thr = rule()["threshold"]
    fs = [out_dir() / f"star_quality_s{int(s):04d}.parquet" for s in sectors]
    q = pd.concat([pd.read_parquet(f) for f in fs if f.exists()], ignore_index=True)
    q["drop"] = (q.status != "ok") | ~(q.beta_max <= thr)
    return q


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sectors", default="70,106")
    ap.add_argument("--workers", type=int, default=2)
    a = ap.parse_args()
    con = store.connect()
    od = out_dir()
    for s in [int(x) for x in a.sectors.split(",")]:
        df = store.frame(con, s)
        df = df[df.status == "ok"]
        man = targets.fetch_manifest(s).set_index("tic")
        jobs = []
        for tic in df.tic:
            tic = int(tic)
            if tic in man.index:
                p = batch.lc_path(s, man.loc[tic].fname, download=False)
                if p is not None:
                    jobs.append(dict(key=f"q:{s}:{tic}", tic=tic, sector=s, path=str(p)))
        with batch.claim(a.workers):
            res = batch.pmap(run_one, jobs, od / f"star_quality_s{s:04d}_parts", a.workers, part_size=1000)
        res["n_searched"] = len(df)
        res.to_parquet(od / f"star_quality_s{s:04d}.parquet", index=False)
        thr = rule()["threshold"]
        print(s, "searched", len(df), "with LC", len(jobs), "ok", int((res.status == "ok").sum()),
              "dropped", int(((res.status != "ok") | ~(res.beta_max <= thr)).sum()))


if __name__ == "__main__":
    main()
