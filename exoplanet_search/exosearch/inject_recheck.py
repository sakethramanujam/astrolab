"""Do v1.0.0 injection rows keep their recovered/vetted status when the search kernels are compiled fresh?

The stored injection rows ran with wotan's cached numba kernels. On hosts where a fresh JIT is used, exact BLS
ties can resolve to a different grid point (docs/v1.0_disclosure_ties_and_numba_cache.md). This recomputes a
seeded sample of stored S70 injection rows, plus every stored row on a gate star with an in-band tie, with the
unchanged v1.0.0 code and a per-process empty NUMBA_CACHE_DIR, logs the tie band / plateaus of each search, and
compares matched / detected / passed_all and the interim S_exact >= S* status.

  python -m exosearch.inject_recheck prepare --n 300 --out /tmp/injrecheck_jobs.pkl     # needs the results DB
  python -m exosearch.inject_recheck run --jobs /tmp/injrecheck_jobs.pkl --workers 6 --cache-dir /tmp --out rc.csv.gz
  python -m exosearch.inject_recheck compare rc.csv.gz --summary docs/host_gate/injection_recheck.json
"""
from __future__ import annotations

import argparse
import json
import os
import pickle
import time
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SEED = 20261110
TIE_REL = 1e-9
FLAGS = ("matched", "detected", "passed_all")


def stored():
    from .inject import PRODUCTS
    d = PRODUCTS / "v1.0.0" / "injections" / "injection_v1"
    x = pd.concat(pd.read_parquet(p) for p in sorted(d.glob("part-*.parquet")))
    return x[(x.sector == 70) & (x.status == "ok")]


def prepare(n, out, gate_csv=None):
    from . import inject, targets
    cfg = json.loads((ROOT / "configs" / "injection_v1.json").read_text())
    x = stored()
    pick = x.sample(min(n, len(x)), random_state=SEED)
    if gate_csv:                                       # every stored row on a gate star whose search had a tie
        g = pd.read_csv(gate_csv)
        tied = set(g.loc[g.n_tie_band.fillna(1) > 1, "tic"].astype(int))
        pick = pd.concat([pick, x[x.tic.isin(tied)]]).drop_duplicates("key")
    want = {}
    for r in pick.itertuples():
        want.setdefault(int(r.tic), set()).add(int(r.run))
    runs = sorted({int(r) for r in pick.run})
    man = targets.fetch_manifest(70).set_index("tic")
    jobs = []
    for j in inject.build_jobs(cfg, runs, [70]):
        if j["tic"] in want:
            keep = [(r, s) for r, s in zip(j["runs"], j["seeds"]) if r in want[j["tic"]]]
            j.update(runs=[r for r, _ in keep], seeds=[s for _, s in keep], path=None, fname=man.loc[j["tic"]].fname)
            jobs.append(j)
    Path(out).write_bytes(pickle.dumps(jobs))
    print(f"{sum(len(j['runs']) for j in jobs)} rows on {len(jobs)} hosts -> {out}")


_TIE = {}


def _init(cache_dir):
    import tempfile
    os.environ["NUMBA_CACHE_DIR"] = tempfile.mkdtemp(prefix="nbc-", dir=cache_dir)
    from . import bls
    orig = bls.search_s

    def search_s(*a, **k):                             # unchanged v1.0.0 search; only logs the tie band
        o = orig(*a, **k)
        p = np.asarray(o["power"], float)
        top = float(np.nanmax(p))
        idx = np.flatnonzero(p >= top - TIE_REL * abs(top))
        runs = np.split(idx, np.flatnonzero(np.diff(idx) > 1) + 1)
        _TIE.update(n_tie_periods=int(idx.size), n_plateaus=len(runs),
                    plateau_width_max=int(max(len(r) for r in runs)))
        return o
    bls.search_s = search_s


def _star(job):
    from . import batch, inject
    job = dict(job, path=str(batch.lc_path(job["sector"], job["fname"], download=True)))
    rows = []
    t, f, e, q, meta = inject.lightcurve.load_spoc_raw(job["path"])
    good = inject.lightcurve.good_mask(f, e, q)
    med = np.nanmedian(f[good])
    lc = (t, f / med, e / med, good, meta)
    for r, sd in zip(job["runs"], job["seeds"]):
        _TIE.clear()
        row = inject.run_one(job, lc, r, sd)
        row.update(_TIE)
        rows.append({k: row.get(k) for k in ("key", "tic", "run", "seed", "status", "error", "P_inj", "P_rec", "t0_rec",
                                             "dur_rec", "depth_rec", "S", "n_tr_rec", "cdpp_d_rec", "s_min",
                                             *FLAGS, "failed_checks", "n_tie_periods", "n_plateaus",
                                             "plateau_width_max")})
    return rows


def run(jobs_file, workers, cache_dir, out):
    from multiprocessing import get_context
    jobs = pickle.loads(Path(jobs_file).read_bytes())
    t0 = time.time()
    with get_context("spawn").Pool(workers, initializer=_init, initargs=(cache_dir,)) as pool:
        rows = [r for rs in pool.imap_unordered(_star, jobs) for r in rs]
    pd.DataFrame(rows).to_csv(out, index=False)
    print(f"{len(rows)} rows in {time.time() - t0:.0f}s -> {out}")


def compare(files, summary, s_star=None):
    if s_star is None:                                 # provisional interim S70 S* (configs/interim_threshold_v1.json)
        s_star = json.loads((ROOT / "configs" / "interim_threshold_v1.json").read_text())["sectors"]["70"]["s_star"]
    new = pd.concat(pd.read_csv(f) for f in files).set_index("key")
    old = stored().set_index("key").loc[new.index]
    se = lambda d, dep, n, c: d[dep] * np.sqrt(d[n]) / d[c]
    out = dict(n_rows=int(len(new)), n_hosts=int(new.tic.nunique()), s_star=s_star,
               tie_band_gt1=int((new.n_tie_periods > 1).sum()), plateau_rows=int((new.plateau_width_max > 1).sum()),
               disjoint_rows=int((new.n_plateaus > 1).sum()),
               P_rec_changed=int((new.P_rec != old.P_rec).sum()), dur_changed=int((new.dur_rec != old.dur_rec).sum()),
               max_abs_dS=float((new.S - old.S).abs().max()))
    for f in FLAGS:
        out[f"{f}_changed"] = int((new[f].astype(bool) != old[f].astype(bool)).sum())
    sn, so = se(new, "depth_rec", "n_tr_rec", "cdpp_d_rec"), se(old, "depth_rec", "n_tr_rec", "cdpp_d_rec")
    rec_new, rec_old = new.matched.astype(bool) & (sn >= s_star), old.matched.astype(bool) & (so >= s_star)
    vet_new = rec_new & new.failed_checks.fillna("").eq("")
    vet_old = rec_old & old.failed_checks.fillna("").eq("")
    out.update(recovered_interim_changed=int((rec_new != rec_old).sum()), vetted_interim_changed=int((vet_new != vet_old).sum()),
               max_abs_dS_exact=float((sn - so).abs().max()))
    ch = new.index[(new.P_rec != old.P_rec) | (rec_new != rec_old) | (vet_new != vet_old)
                   | np.any([new[f].astype(bool) != old[f].astype(bool) for f in FLAGS], axis=0)]
    out["changed_rows"] = [dict(key=k, P_old=float(old.P_rec[k]), P_new=float(new.P_rec[k]),
                                matched=[bool(old.matched[k]), bool(new.matched[k])],
                                passed_all=[bool(old.passed_all[k]), bool(new.passed_all[k])],
                                n_tie_periods=int(new.n_tie_periods[k])) for k in ch]
    out["completeness_unaffected"] = not any(out[k] for k in ("matched_changed", "detected_changed", "passed_all_changed",
                                                               "recovered_interim_changed", "vetted_interim_changed"))
    if summary:
        Path(summary).write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps({k: v for k, v in out.items() if k != "changed_rows"}, indent=1))
    return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    sp = ap.add_subparsers(dest="cmd", required=True)
    p = sp.add_parser("prepare"); p.add_argument("--n", type=int, default=300); p.add_argument("--out", required=True)
    p.add_argument("--gate-csv", default=None)
    r = sp.add_parser("run"); r.add_argument("--jobs", required=True); r.add_argument("--workers", type=int, default=2)
    r.add_argument("--cache-dir", default="/tmp"); r.add_argument("--out", required=True)
    c = sp.add_parser("compare"); c.add_argument("files", nargs="+"); c.add_argument("--summary", default=None)
    c.add_argument("--s-star", type=float, default=None)
    a = ap.parse_args()
    if a.cmd == "prepare":
        prepare(a.n, a.out, a.gate_csv)
    elif a.cmd == "run":
        run(a.jobs, a.workers, a.cache_dir, a.out)
    else:
        compare(a.files, a.summary, a.s_star)
