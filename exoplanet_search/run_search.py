#!/usr/bin/env python3
"""TESS transit search: manifest -> download -> detrend+BLS per star (Dask) -> vet -> cross-match -> plots.

Examples
  # local CPU run, 8 polite workers
  nice -n 19 .venv/bin/python run_search.py --sector 70 --n-targets 500 --workers 8
  # use an external scheduler that CPU / CUDA / MLX workers have joined
  nice -n 19 .venv/bin/python run_search.py --sector 70 --n-targets 2000 --scheduler tcp://<scheduler-host>:8786
"""
from __future__ import annotations

import argparse
import os
import shutil
import sys
import tempfile
import time
from pathlib import Path

# one thread per worker process; keeps CPU use predictable on a shared box
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

sys.path.insert(0, str(Path(__file__).resolve().parent))
from exosearch import catalogs, config, plots, store, targets, vet  # noqa: E402
from exosearch.task import process_star  # noqa: E402

RESOURCES = {"cuda": {"GPU": 1}, "mlx": {"MLX": 1}}


def log(*a):
    print(time.strftime("%H:%M:%S"), *a, flush=True)


def get_client(args):
    from dask.distributed import Client, LocalCluster
    if args.scheduler:
        client = Client(args.scheduler)
        if args.upload_code:  # ship exosearch/ to workers that don't have the repo checked out
            tmp = Path(tempfile.mkdtemp(prefix="exosearch_", dir="/tmp"))
            z = shutil.make_archive(str(tmp / "exosearch"), "zip", Path(__file__).parent, "exosearch")
            client.upload_file(z)
            shutil.rmtree(tmp, ignore_errors=True)
        return client, None
    cluster = LocalCluster(n_workers=args.workers, threads_per_worker=1, processes=True,
                           host=args.listen, scheduler_port=args.port,
                           dashboard_address=f"{args.listen}:{args.dashboard_port}")
    return Client(cluster), cluster


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sector", type=int, default=70)
    ap.add_argument("--n-targets", type=int, default=500)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--backend", choices=["auto", "cpu", "numpy", "cuda", "mlx"], default="auto",
                    help="auto = each worker picks by its resource tag (GPU=1 -> cuda, MLX=1 -> mlx, else cpu)")
    ap.add_argument("--scheduler", help="tcp://host:8786 of an existing Dask scheduler")
    ap.add_argument("--upload-code", action="store_true", help="upload exosearch/ to remote workers")
    ap.add_argument("--workers", type=int, default=8, help="LocalCluster worker processes (no --scheduler)")
    ap.add_argument("--listen", default="127.0.0.1", help="LocalCluster bind address (LAN/Tailscale IP to let others join)")
    ap.add_argument("--port", type=int, default=0)
    ap.add_argument("--dashboard-port", type=int, default=0, help="0 = random free port")
    ap.add_argument("--download-threads", type=int, default=8)
    ap.add_argument("--redo", action="store_true", help="re-run stars already in the DB")
    ap.add_argument("--top", type=int, default=15)
    ap.add_argument("--plots", type=int, default=10)
    args = ap.parse_args()

    t_all = time.perf_counter()
    man = targets.fetch_manifest(args.sector)
    sel = targets.select_targets(man, args.n_targets, args.seed)
    log(f"sector {args.sector}: {len(man)} SPOC 2-min targets, selected {len(sel)} (seed {args.seed})")

    t0 = time.perf_counter()
    dl = targets.download_many(sel, args.sector, threads=args.download_threads, log=log)
    log(f"download/cache check: {len(dl['ok'])} ok, {len(dl['failed'])} failed in {time.perf_counter()-t0:.0f}s")

    con = store.connect()
    done = set() if args.redo else store.done_tics(con, args.sector)
    todo = sel[sel.tic.isin(dl["ok"].keys()) & ~sel.tic.isin(done)]
    log(f"{len(done & set(sel.tic))} already done, {len(todo)} to search, backend={args.backend}")

    t_search = time.perf_counter()
    if len(todo):
        from dask.distributed import as_completed
        client, cluster = get_client(args)
        log(f"dask: {client.dashboard_link}  workers={len(client.scheduler_info()['workers'])}")
        res = RESOURCES.get(args.backend)
        futs = [client.submit(process_star, r.tic, args.sector, r.fname, r.url, args.backend,
                              key=f"star-{args.sector}-{r.tic}-{args.backend}", resources=res, pure=False)
                for r in todo.itertuples()]
        n_err = 0
        for i, fut in enumerate(as_completed(futs), 1):
            row = fut.result()
            n_err += row["status"] != "ok"
            store.upsert(con, row)
            if i % 25 == 0 or i == len(futs):
                el = time.perf_counter() - t_search
                log(f"  searched {i}/{len(futs)}  errors={n_err}  wall {el:.0f}s ({el/i:.2f}s/star wall)")
        client.close()
        if cluster:
            cluster.close()
    search_wall = time.perf_counter() - t_search

    log("cross-matching against ExoFOP TOI/CTOI + NASA Exoplanet Archive")
    cats = catalogs.load_catalogs()
    df = store.frame(con, args.sector)
    for r in df[df.status == "ok"].itertuples():
        row = r._asdict()
        row.pop("Index", None)
        row.update(catalogs.crossmatch(r.tic, r.period, cats))
        store.upsert(con, row)
    df = vet.flag_common_periods(store.frame(con, args.sector))
    for r in df.itertuples():
        row = r._asdict()
        row.pop("Index", None)
        store.upsert(con, row)
    out_dir = config.RESULTS_DIR
    df.to_parquet(out_dir / f"results_s{args.sector:04d}.parquet", index=False)

    ok = df[df.status == "ok"].copy()
    ok["known"] = ok.apply(lambda r: ("TOI " + str(r.known_toi)) if r.known_toi else
                           (r.known_planet or (("CTOI " + str(r.known_ctoi)) if r.known_ctoi else "")), axis=1)
    top = ok.sort_values(["passed", "sde"], ascending=[False, False]).head(args.top)
    cols = ["tic", "tmag", "period", "duration", "depth", "sde", "snr", "odd_even_sigma",
            "secondary_sigma", "rp_earth", "passed", "flags", "known", "toi_period_match", "planet_period_match"]
    import pandas as pd
    with pd.option_context("display.width", 250, "display.max_columns", 30):
        print(top[cols].to_string(index=False, float_format=lambda x: f"{x:.4g}"))
    top[cols].to_csv(out_dir / f"top_s{args.sector:04d}.csv", index=False)

    for r in top.head(args.plots).itertuples():
        row = r._asdict()
        p = config.lc_dir(args.sector) / sel.set_index("tic").fname.get(r.tic, "")
        if p.is_file():
            png = plots.diagnostic(p, row, config.plots_dir(args.sector) / f"TIC{r.tic}_s{args.sector:04d}.png")
            log(f"  plot {png}")

    rt = ok.runtime_s
    log(f"stars ok={len(ok)} errors={int((df.status!='ok').sum())}  per-star task time "
        f"median {rt.median():.2f}s mean {rt.mean():.2f}s  search wall {search_wall:.0f}s  "
        f"total wall {time.perf_counter()-t_all:.0f}s")


if __name__ == "__main__":
    main()
