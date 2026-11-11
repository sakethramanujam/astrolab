"""Cross-host reproducibility gate (dry run on the current interim code; no production rows are written).

Every host searches the same ~200 S70 stars (`configs/host_gate_tics.txt`) with the frozen v1.0.0 search path
(task.prepare + task.search_prepared, float64 CPU, 1 thread per worker) and records S, S_exact, best (P, D, T0).
Gate: for every star, |S_host - S_ref| <= 1e-6 and |S_exact_host - S_exact_ref| <= 1e-6 and an identical best
(P, D) across hosts. Rows carry only an opaque host_id (h1, h2, ...), arch and library versions.

  python -m exosearch.host_gate pick                         # writes configs/host_gate_tics.txt (seeded)
  python -m exosearch.host_gate run --host-id h2 --workers 12 --out /tmp/gate_h2.csv.gz
  python -m exosearch.host_gate compare gate_h1.csv.gz gate_h2.csv.gz ... --summary docs/host_gate/summary.json
"""
from __future__ import annotations

import argparse
import json
import os
import platform
import time

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMBA_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np
import pandas as pd

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TICS = ROOT / "configs" / "host_gate_tics.txt"
SECTOR, N, SEED = 70, 200, 20261110
TOL = 1e-6
LIBS = ("numpy", "scipy", "astropy", "numba", "llvmlite", "wotan", "lightkurve", "pandas")


def lib_versions():
    import importlib.metadata as md
    v = {k: md.version(k) for k in LIBS}
    try:
        cfg = np.show_config(mode="dicts")
        v["blas"] = f'{cfg["Build Dependencies"]["blas"]["name"]} {cfg["Build Dependencies"]["blas"].get("version", "")}'.strip()
    except Exception:
        v["blas"] = "unknown"
    v["python"] = platform.python_version()
    return v


def pick():
    from . import targets
    man = targets.fetch_manifest(SECTOR)
    tics = man.tic.astype(int).sort_values().sample(N, random_state=SEED).sort_values()
    TICS.write_text("\n".join(map(str, tics)) + "\n")
    print(f"{len(tics)} TICs -> {TICS.relative_to(ROOT)}")


def _one(args):
    tic, fname = args
    from . import batch, task
    t0 = time.perf_counter()
    try:
        p = batch.lc_path(SECTOR, fname, download=True)
        dl = time.perf_counter() - t0
        r = task.search_prepared(task.prepare(p))
        se = r["depth"] * np.sqrt(r["n_tr"]) / r["cdpp_d"]
        return dict(tic=tic, status="ok", S=r["S"], S_exact=se, period=r["period"], duration=r["duration"], t0=r["t0"],
                    depth_box=r["depth_box"], depth=r["depth"], n_tr=r["n_tr"], cdpp_d=r["cdpp_d"], beta_d=r["beta_d"],
                    download_s=dl, search_s=time.perf_counter() - t0 - dl)
    except Exception as exc:                                      # error text kept short and path-free
        return dict(tic=tic, status="error", error=type(exc).__name__)


def run(host_id, workers, out):
    from multiprocessing import get_context
    from . import targets
    man = targets.fetch_manifest(SECTOR).set_index("tic")
    tics = [int(x) for x in TICS.read_text().split()]
    jobs = [(t, man.loc[t].fname) for t in tics]
    t0 = time.time()
    with get_context("spawn").Pool(workers) as pool:
        rows = list(pool.imap_unordered(_one, jobs))
    wall = time.time() - t0
    d = pd.DataFrame(rows).sort_values("tic")
    meta = dict(host_id=host_id, arch=platform.machine(), system=platform.system(), lib_versions=lib_versions(),
                workers=workers, wall_s=wall)
    d["host_id"], d["arch"], d["lib_versions"] = host_id, meta["arch"], json.dumps(meta["lib_versions"], sort_keys=True)
    d.to_csv(out, index=False)
    ok = d[d.status == "ok"]
    meta.update(n_ok=int(len(ok)), n_err=int((d.status != "ok").sum()),
                search_s_per_star=float(ok.search_s.median()) if len(ok) else None,
                stars_per_hour=float(3600 * workers / ok.search_s.mean()) if len(ok) else None)
    Path(str(out) + ".meta.json").write_text(json.dumps(meta, indent=1))
    print(json.dumps(meta))


def compare(files, summary):
    ds = [pd.read_csv(f) for f in files]
    metas = [json.loads(Path(str(f) + ".meta.json").read_text()) for f in files]
    ref = ds[0][ds[0].status == "ok"].set_index("tic")
    hosts = []
    for d, m in zip(ds, metas):
        x = d[d.status == "ok"].set_index("tic")
        common = ref.index.intersection(x.index)
        a, b = ref.loc[common], x.loc[common]
        dS = (a.S - b.S).abs()
        dSx = (a.S_exact - b.S_exact).abs()
        same_pd = (a.period == b.period) & (a.duration == b.duration)
        bad = common[(dS > TOL) | (dSx > TOL) | ~same_pd]
        hosts.append(dict(host_id=m["host_id"], arch=m["arch"], system=m["system"], lib_versions=m["lib_versions"],
                          workers=m["workers"], n_ok=m["n_ok"], n_err=m["n_err"], n_compared=int(len(common)),
                          max_abs_dS=float(dS.max()), max_abs_dS_exact=float(dSx.max()),
                          n_PD_mismatch=int((~same_pd).sum()), n_fail=int(len(bad)),
                          fail_tics=[int(t) for t in bad[:20]], stars_per_hour=m.get("stars_per_hour"),
                          search_s_per_star_median=m.get("search_s_per_star"),
                          passes=bool(len(bad) == 0 and m["n_err"] == 0 and len(common) == len(ref))))
    res = dict(gate="cross-host dry run, v1.0.0 search path + interim S_exact, sector 70",
               reference=metas[0]["host_id"], tolerance=TOL, n_stars=int(len(ref)), hosts=hosts,
               all_pass=all(h["passes"] for h in hosts))
    if summary:
        Path(summary).write_text(json.dumps(res, indent=1))
    print(json.dumps(res, indent=1))


def main():
    ap = argparse.ArgumentParser()
    sp = ap.add_subparsers(dest="cmd", required=True)
    sp.add_parser("pick")
    r = sp.add_parser("run"); r.add_argument("--host-id", required=True); r.add_argument("--workers", type=int, default=4); r.add_argument("--out", required=True)
    c = sp.add_parser("compare"); c.add_argument("files", nargs="+"); c.add_argument("--summary")
    a = ap.parse_args()
    if a.cmd == "pick":
        pick()
    elif a.cmd == "run":
        run(a.host_id, a.workers, a.out)
    else:
        compare(a.files, a.summary)


if __name__ == "__main__":
    main()
