"""Cross-host reproducibility gate (dry run on the current interim code; no production rows are written).

Every host searches the same ~200 S70 stars (`configs/host_gate_tics.txt`) with the frozen v1.0.0 search path
(task.prepare + task.search_prepared, float64 CPU, 1 thread per worker) and records S, S_exact, best (P, D, T0).
Gate: for every star, |S_host - S_ref| <= 1e-6 and |S_exact_host - S_exact_ref| <= 1e-6 and an identical best
(P, D) across hosts. Rows carry only an opaque host_id (h1, h2, ...), arch and library versions.

  python -m exosearch.host_gate pick                         # writes configs/host_gate_tics.txt (seeded)
  python -m exosearch.host_gate run --host-id h2 --workers 12 --out /tmp/gate_h2.csv.gz
  python -m exosearch.host_gate trace 26016053                # per-stage hashes, to localise a divergence
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


TIE_REL = 1e-9          # registered tie band: |S_i - S_max| <= TIE_REL * S_max (docs/preregistration/tie_break.md)


def apply_fixes():
    """Dry-run fixes for cross-host reproducibility (interim code only; v1.0.0 untouched):
    1. period grid = fmin + df * arange(n) (no fused multiply-add in np.arange's fill loop);
    2. ties within TIE_REL * S_max -> lowest period, then shortest duration; band size logged per star;
    3. wotan's numba kernels recompiled with fastmath=False, cache=False (and parallel=False)."""
    import numba
    import wotan
    import importlib
    import pkgutil
    from numba.core.registry import CPUDispatcher
    from astropy.timeseries import BoxLeastSquares
    from . import bls
    mods = [importlib.import_module(f"wotan.{m.name}") for m in pkgutil.iter_modules(wotan.__path__)]
    new = {}
    for m in mods:                                   # one replacement per original kernel, shared by every module
        for k, v in list(vars(m).items()):
            if isinstance(v, CPUDispatcher):
                if id(v) not in new:
                    opts = dict(v.targetoptions); opts.update(fastmath=False, cache=False, parallel=False)
                    new[id(v)] = numba.jit(**opts)(v.py_func)
                setattr(m, k, new[id(v)])

    def period_grid(t, pmin=0.5, pmax=None, min_transits=2, frequency_factor=1.0, durations=bls.DURATIONS):
        baseline = float(t.max() - t.min())
        pmax = pmax or baseline / min_transits
        df = frequency_factor * durations.min() / baseline ** 2
        fmin = 1.0 / pmax
        n = int(np.ceil((1.0 / pmin - fmin) / df))
        freqs = fmin + df * np.arange(n, dtype=float)
        freqs = freqs[freqs < 1.0 / pmin]
        return np.sort(1.0 / freqs)

    def search_s(tb, fb, tflat, fflat, periods=None, durations=bls.DURATIONS):
        if periods is None:
            periods = period_grid(tb, durations=durations)
        sw, cad, prof = bls.noise_profile(tflat, fflat, durations)
        bin_d = float(np.median(np.diff(tb)))
        smax = np.full(len(periods), -np.inf)
        per_d = []
        for D in durations:
            cd, beta = prof[float(D)]
            if not (cd == cd and cd > 0):
                continue
            dy = np.full(tb.size, cd * np.sqrt(max(D / bin_d, 1.0)))
            r = BoxLeastSquares(tb, fb, dy=dy).power(periods, [D], objective="snr", method="fast")
            p = np.nan_to_num(np.asarray(r.power), nan=-np.inf)
            smax = np.maximum(smax, p)
            per_d.append((float(D), cd, beta, r, p))
        if not per_d:
            raise ValueError("no valid noise estimate")
        top = max(float(p.max()) for *_, p in per_d)
        band = top - TIE_REL * abs(top)
        cands = []
        for D, cd, beta, r, p in per_d:
            for i in np.flatnonzero(p >= band):
                cands.append((float(r.period[i]), D, i, cd, beta, r, p))
        cands.sort(key=lambda c: (c[0], c[1]))                      # lowest period, then shortest duration
        P, D, i, cd, beta, r, p = cands[0]
        best = dict(period=P, t0=float(r.transit_time[i]), duration=D, depth_box=float(r.depth[i]), S_scan=float(p[i]),
                    cdpp_d=cd, beta_d=beta, n_tie_band=len(cands))
        ntr = bls.transit_count(tb, best["period"], best["t0"], best["duration"], bin_d)
        s_val = best["depth_box"] * np.sqrt(ntr) / best["cdpp_d"] if ntr > 0 else 0.0
        smax[~np.isfinite(smax)] = 0.0
        sd, pdet = bls.sde(smax)
        best.update(S=float(s_val), n_tr=ntr, sigma_w=sw, cadence_d=cad, sde=sd, periods=periods, power=smax,
                    power_detrended=pdet, backend="cpu-S-fixed")
        return best

    bls.period_grid, bls.search_s = period_grid, search_s
    return len(new)


def pick():
    from . import targets
    man = targets.fetch_manifest(SECTOR)
    tics = man.tic.astype(int).sort_values().sample(N, random_state=SEED).sort_values()
    TICS.write_text("\n".join(map(str, tics)) + "\n")
    print(f"{len(tics)} TICs -> {TICS.relative_to(ROOT)}")


_FIXED = {"on": False}
_TIE = {}


def _init(fixed, cache_dir):
    if cache_dir:                                   # per-process empty numba cache: every kernel compiled fresh
        import tempfile
        os.environ["NUMBA_CACHE_DIR"] = tempfile.mkdtemp(prefix="nbc-", dir=cache_dir)
    if fixed:
        apply_fixes()
        _FIXED["on"] = True


def _one(args):
    tic, fname = args[:2]
    sector = args[2] if len(args) > 2 else SECTOR
    from . import batch, task
    t0 = time.perf_counter()
    try:
        p = batch.lc_path(sector, fname, download=True)
        dl = time.perf_counter() - t0
        d = task.prepare(p)
        if _FIXED["on"]:
            from . import bls
            o = bls.search_s(d["binned"][0], d["binned"][1], d["flat"][0], d["flat"][1])
            st = bls.stats(*d["binned"], o["period"], o["duration"], o["t0"])
            dep = float(st["depth"][0]) if np.ndim(st["depth"]) else float(st["depth"])
            r = dict(S=o["S"], depth=dep, n_tr=o["n_tr"], cdpp_d=o["cdpp_d"], period=o["period"], duration=o["duration"],
                     t0=o["t0"], depth_box=o["depth_box"], beta_d=o["beta_d"], n_tie_band=o["n_tie_band"])
        else:
            r = task.search_prepared(d)
        se = r["depth"] * np.sqrt(r["n_tr"]) / r["cdpp_d"]
        return dict(tic=tic, status="ok", S=r["S"], S_exact=se, period=r["period"], duration=r["duration"], t0=r["t0"],
                    depth_box=r["depth_box"], depth=r["depth"], n_tr=r["n_tr"], cdpp_d=r["cdpp_d"], beta_d=r["beta_d"],
                    n_tie_band=r.get("n_tie_band"),
                    download_s=dl, search_s=time.perf_counter() - t0 - dl)
    except Exception as exc:                                      # error text kept short and path-free
        return dict(tic=tic, status="error", error=type(exc).__name__)


def run(host_id, workers, out, fixed=False, fresh_cache_dir=None, jobs=None, sector=SECTOR):
    from multiprocessing import get_context
    from . import targets
    if jobs is None:
        man = targets.fetch_manifest(SECTOR).set_index("tic")
        tics = [int(x) for x in TICS.read_text().split()]
        jobs = [(t, man.loc[t].fname, sector) for t in tics]
    t0 = time.time()
    with get_context("spawn").Pool(workers, initializer=_init, initargs=(fixed, fresh_cache_dir)) as pool:
        rows = list(pool.imap_unordered(_one, jobs))
    wall = time.time() - t0
    d = pd.DataFrame(rows).sort_values("tic")
    meta = dict(host_id=host_id, arch=platform.machine(), system=platform.system(), lib_versions=lib_versions(),
                workers=workers, wall_s=wall, fixes_applied=bool(fixed), fresh_numba_cache=bool(fresh_cache_dir))
    d["host_id"], d["arch"], d["lib_versions"] = host_id, meta["arch"], json.dumps(meta["lib_versions"], sort_keys=True)
    d.to_csv(out, index=False)
    ok = d[d.status == "ok"]
    meta.update(n_ok=int(len(ok)), n_err=int((d.status != "ok").sum()),
                search_s_per_star=float(ok.search_s.median()) if len(ok) else None,
                stars_per_hour=float(3600 * workers / ok.search_s.mean()) if len(ok) else None)
    Path(str(out) + ".meta.json").write_text(json.dumps(meta, indent=1))
    print(json.dumps(meta))


def v10check(workers, out, cache_dir, n=N, seed=SEED, sector=105):
    """v1.0.0 exposure: do freshly compiled wotan kernels reproduce the stored v1.0.0 rows (mostly cached code)?"""
    from . import store, targets
    d = store.frame(store.connect(), sector)
    d = d[(d.status == "ok") & (d.pipeline_version == "1.0.0") & d.S.notna()]
    d = d.sample(min(n, len(d)), random_state=seed).sort_values("tic")
    man = targets.fetch_manifest(sector).set_index("tic")
    jobs = [(int(t), man.loc[int(t)].fname, sector) for t in d.tic]
    run("h1", workers, out, fixed=False, fresh_cache_dir=cache_dir, jobs=jobs, sector=sector)
    x = pd.read_csv(out).set_index("tic")
    ref = d.set_index("tic").loc[x.index]
    ok = x.status == "ok"
    same_pd = (x.period == ref.period) & (x.duration == ref.duration)
    dS = (x.S - ref.S).abs()
    res = dict(check="v1.0.0 rows (S105 real search) recomputed with per-process empty numba caches (fresh JIT)",
               n=int(len(x)), n_ok=int(ok.sum()), n_PD_mismatch=int((ok & ~same_pd).sum()),
               n_dS_gt_tol=int((ok & (dS > TOL)).sum()), max_abs_dS=float(dS[ok].max()),
               git_hashes_in_sample=sorted(map(str, ref.git_hash.dropna().unique())),
               mismatch_tics=[int(t) for t in x.index[ok & (~same_pd | (dS > TOL))][:30]])
    Path(str(out) + ".check.json").write_text(json.dumps(res, indent=1))
    print(json.dumps(res, indent=1))


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
        relP = (a.period / b.period - 1).abs()
        ulp = ~same_pd & (a.duration == b.duration) & (relP < 1e-12)          # same grid point, last-bit difference
        flip = ~same_pd & ~ulp                                                  # a different best peak
        flips = [dict(tic=int(t), P_ref=float(a.period[t]), P_host=float(b.period[t]), D_ref=float(a.duration[t]),
                      D_host=float(b.duration[t]), S_ref=float(a.S[t]), dS=float(b.S[t] - a.S[t])) for t in common[flip]]
        hosts.append(dict(host_id=m["host_id"], arch=m["arch"], system=m["system"], lib_versions=m["lib_versions"],
                          workers=m["workers"], n_ok=m["n_ok"], n_err=m["n_err"], n_compared=int(len(common)),
                          max_abs_dS=float(dS.max()), max_abs_dS_exact=float(dSx.max()),
                          n_PD_mismatch=int((~same_pd).sum()), n_P_last_bit_only=int(ulp.sum()), n_peak_flip=int(flip.sum()),
                          peak_flips=flips, max_rel_dP_same_grid_point=float(relP[same_pd | ulp].max()), n_fail=int(len(bad)),
                          fail_tics=[int(t) for t in bad[:20]], stars_per_hour=m.get("stars_per_hour"),
                          search_s_per_star_median=m.get("search_s_per_star"),
                          fixes_applied=m.get("fixes_applied", False),
                          tie_band_gt1=None if "n_tie_band" not in x or x.n_tie_band.isna().all() else
                          dict(n=int((x.n_tie_band > 1).sum()), frac=round(float((x.n_tie_band > 1).mean()), 4),
                               max=int(x.n_tie_band.max())),
                          passes=bool(len(bad) == 0 and m["n_err"] == 0 and len(common) == len(ref))))
    res = dict(gate="cross-host dry run, v1.0.0 search path + interim S_exact, sector 70",
               reference=metas[0]["host_id"], tolerance=TOL, n_stars=int(len(ref)), hosts=hosts,
               all_pass=all(h["passes"] for h in hosts))
    if summary:
        Path(summary).write_text(json.dumps(res, indent=1))
    print(json.dumps(res, indent=1))


def trace(tics):
    """Per-stage array hashes for a star (raw, trend, flat, binned, period grid, CDPP) to localise a divergence."""
    import hashlib
    from . import batch, bls, targets, task
    h = lambda a: hashlib.sha256(np.ascontiguousarray(np.asarray(a, dtype=float)).tobytes()).hexdigest()[:12]
    man = targets.fetch_manifest(SECTOR).set_index("tic")
    for tic in tics:
        d = task.prepare(batch.lc_path(SECTOR, man.loc[tic].fname))
        t, f, _ = d["raw"]; _, ff, _ = d["flat"]; tb, fb, _ = d["binned"]
        _, _, prof = bls.noise_profile(d["flat"][0], ff)
        print(json.dumps(dict(tic=tic, raw_t=h(t), raw_f=h(f), trend=h(d["trend"]), flat=h(ff), binned_t=h(tb), binned_f=h(fb),
                              grid=h(bls.period_grid(tb)), cdpp=h([prof[k][0] for k in sorted(prof)]))))


def main():
    ap = argparse.ArgumentParser()
    sp = ap.add_subparsers(dest="cmd", required=True)
    sp.add_parser("pick")
    r = sp.add_parser("run"); r.add_argument("--host-id", required=True); r.add_argument("--workers", type=int, default=4); r.add_argument("--out", required=True)
    r.add_argument("--fixed", action="store_true", help="apply the dry-run reproducibility fixes (grid, tie-break, wotan no fastmath/cache)")
    r.add_argument("--fresh-cache-dir", default=None, help="parent dir for per-process empty numba caches")
    v = sp.add_parser("v10check", help="recompute ~200 stored v1.0.0 S105 rows with fresh numba compiles; compare to the DB")
    v.add_argument("--workers", type=int, default=2); v.add_argument("--out", required=True); v.add_argument("--cache-dir", default="/tmp")
    tr = sp.add_parser("trace"); tr.add_argument("tics", nargs="+", type=int)
    c = sp.add_parser("compare"); c.add_argument("files", nargs="+"); c.add_argument("--summary")
    a = ap.parse_args()
    if a.cmd == "pick":
        pick()
    elif a.cmd == "trace":
        trace(a.tics)
    elif a.cmd == "run":
        run(a.host_id, a.workers, a.out, fixed=a.fixed, fresh_cache_dir=a.fresh_cache_dir)
    elif a.cmd == "v10check":
        v10check(a.workers, a.out, a.cache_dir)
    else:
        compare(a.files, a.summary)


if __name__ == "__main__":
    main()
