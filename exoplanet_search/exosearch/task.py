"""One star = one task. Runs on any Dask worker (CPU, CUDA or MLX)."""
from __future__ import annotations

import datetime as dt
import os
import socket
import time

import numpy as np

from . import bls, config, lightcurve, targets, vet


def resolve_backend(requested: str) -> str:
    """'auto' picks by the worker's resource tags: GPU=1 -> cuda, MLX=1 -> mlx, else cpu."""
    if requested != "auto":
        return requested
    try:
        from distributed import get_worker
        w = get_worker()
        res = getattr(getattr(w, "state", None), "total_resources", None) or getattr(w, "total_resources", {}) or {}
        if res.get("GPU"):
            return "cuda"
        if res.get("MLX"):
            return "mlx"
    except (ImportError, ValueError):
        pass
    return "cpu"


def prepare(path, window=0.75, bin_minutes=10.0):
    t, f, e, meta = lightcurve.load_spoc(path)
    if len(t) < 500:
        raise ValueError(f"too few good cadences ({len(t)})")
    tf, ff, ef, trend, _ = lightcurve.detrend(t, f, e, window=window)
    tb, fb, eb = lightcurve.bin_lc(tf, ff, ef, minutes=bin_minutes) if bin_minutes else (tf, ff, ef)
    return dict(raw=(t, f, e), trend=trend, flat=(tf, ff, ef), binned=(tb, fb, eb), meta=meta)


def process_star(tic, sector, fname, url, backend="auto", window=0.75, bin_minutes=10.0):
    t_start = time.perf_counter()
    be = resolve_backend(backend)
    row = dict(tic=int(tic), sector=int(sector), backend=be, host=socket.gethostname(),
               processed_at=dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"))
    try:
        # Remote workers without the shared disk download into their own EXO_DATA_ROOT.
        path = targets.download_one(fname, url, config.lc_dir(sector))
        d = prepare(path, window, bin_minutes)
        tb, fb, eb = d["binned"]
        out = bls.search(tb, fb, eb, backend=be)
        st = bls.stats(tb, fb, eb, out["period"], out["duration"], out["t0"])
        row.update({k: d["meta"].get(k) for k in ("tmag", "teff", "radius", "ra", "dec")})
        row.update(n_points=int(len(d["raw"][0])), period=out["period"], t0=out["t0"],
                   duration=out["duration"], sde=out["sde"])
        row.update(vet.vet(st, out, d["meta"], tb, fb))
        row["status"] = "ok"
    except Exception as exc:
        row.update(status="error", error=f"{type(exc).__name__}: {exc}"[:500])
    row["runtime_s"] = time.perf_counter() - t_start
    for k, v in list(row.items()):  # plain Python types for SQLite
        if isinstance(v, (np.floating, np.integer)):
            row[k] = v.item()
    return row
