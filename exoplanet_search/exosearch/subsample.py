"""Random-subsample runs through the full v1.0.0 chain (analysis layer; never writes the results DB or an S*).

  python -m exosearch.subsample --sector 105 --mode invert --n 2000 --workers 3   # flipped null, S105 subsample
  python -m exosearch.subsample --sector 70  --mode real   --n 2000 --workers 3   # v1 real search, S70 subsample

Stars: searched (status ok), tuning (non-holdout); for --mode real on a sector with a star-quality table the draw is
from quality-passing stars (the population the interim S* uses). Seeded draw. Each star goes through
nulltest.run_one (same search, screens and deep checks as the nulls); mode 'real' leaves the flux untouched.
Output: products/v<ver>/analysis_layer/subsample/s<sector>_<mode>_n<n>/{parts/, rows.parquet, summary.json}.
Every output is labelled "subsample, not an S*".
"""
from __future__ import annotations

import argparse
import json

import numpy as np
import pandas as pd

from . import batch, config, holdout, nulltest, quality, stellar, store, targets
from .version import stamp

LABEL = "subsample, not an S*"
_transform = nulltest.transform


def _t(mode, t, f, e, meta, seed):
    return (t, f, e) if mode == "real" else _transform(mode, t, f, e, meta, seed)


nulltest.transform = _t          # analysis-layer hook; the frozen module is unchanged


def draw(sector, n, mode, seed):
    df = store.frame(store.connect(), sector)
    df = df[df.status == "ok"]
    tics = [int(t) for t in df.tic if not holdout.is_holdout(int(t))]
    pop = "searched, tuning"
    if mode == "real":
        try:
            q = quality.load()
            keep = set(q[(q.sector == sector) & ~q["drop"]].tic.astype(int))
            if keep:
                tics, pop = [t for t in tics if t in keep], "searched, tuning, star_quality_v1 kept"
        except Exception:
            pass
    rng = np.random.default_rng(seed)
    return sorted(int(x) for x in rng.choice(tics, size=min(n, len(tics)), replace=False)), pop, len(tics)


def summary(d, sector, mode, pop, n_pop):
    ok = d[d.status == "ok"]
    fl = ok["flags"].fillna("").str.split(",").map(lambda x: set(x) - {"", "low_s"})
    out = dict(label=LABEL, sector=sector, mode=mode, population=pop, population_size=n_pop, n=int(len(ok)),
               n_error=int((d.status != "ok").sum()),
               max_S_quantiles={str(p): float(ok.S.quantile(p)) for p in (.1, .25, .5, .75, .9, .95, .99)},
               median_beta_d=float(ok.beta_d.median()), frac_beta_lt_1=float((ok.beta_d < 1).mean()),
               frac_S_ge_12p4=float((ok.S >= 12.4).mean()), n_S_ge_12p4=int((ok.S >= 12.4).sum()),
               frac_S_ge_12p4_screens_passed=float(((ok.S >= 12.4) & (fl.map(len) == 0)).mean()),
               hist_edges=list(range(0, 41)), hist=np.histogram(ok.S.clip(0, 39.99), bins=range(0, 41))[0].tolist())
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sector", type=int, required=True)
    ap.add_argument("--mode", choices=["invert", "real"], required=True)
    ap.add_argument("--n", type=int, default=2000)
    ap.add_argument("--workers", type=int, default=3)
    ap.add_argument("--seed", type=int, default=20261110)
    a = ap.parse_args()
    st = stamp()
    out = config.DATA_ROOT / "products" / f"v{st['pipeline_version']}" / "analysis_layer" / "subsample" / f"s{a.sector:04d}_{a.mode}_n{a.n}"
    out.mkdir(parents=True, exist_ok=True)
    tics, pop, n_pop = draw(a.sector, a.n, a.mode, a.seed)
    (out / "provenance.json").write_text(json.dumps(dict(st, label=LABEL, sector=a.sector, mode=a.mode, n=a.n, seed=a.seed,
                                                         population=pop, population_size=n_pop), indent=1))
    man = targets.fetch_manifest(a.sector).set_index("tic")
    tp = stellar.tic_params(tics)
    jobs = []
    for tic in tics:
        if tic not in man.index:
            continue
        x = tp.loc[tic].to_dict() if tic in tp.index else {}
        rho, erho = stellar.tic_rho(x) if x else (None, None)
        jobs.append(dict(key=f"sub:{a.mode}:{a.sector}:{tic}", tic=tic, sector=a.sector, mode=a.mode, realization=0,
                         seed=nulltest.seed_for("invert", a.sector, tic, 0), fname=man.loc[tic].fname, stamp=st,
                         rho=rho, rho_err=erho, teff=x.get("Teff"), logg=x.get("logg"), ruwe=None, nss=None))
    with batch.claim(a.workers):
        d = batch.pmap(run_one, jobs, out / "parts", a.workers, part_size=100)
    d["label"] = LABEL
    d.to_parquet(out / "rows.parquet", index=False)
    s = summary(d, a.sector, a.mode, pop, n_pop)
    (out / "summary.json").write_text(json.dumps(s, indent=1))
    print(json.dumps({k: v for k, v in s.items() if k not in ("hist", "hist_edges")}, indent=1))


def run_one(job):
    job = dict(job)
    try:
        job["path"] = str(batch.lc_path(job["sector"], job.pop("fname"), download=True))
    except Exception as exc:
        return dict(key=job["key"], tic=job["tic"], sector=job["sector"], mode=job["mode"], status="error",
                    error=f"download {type(exc).__name__}"[:200])
    return nulltest.run_one(job)


if __name__ == "__main__":
    main()
