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


def prepare_arrays(t, f, e, meta, window=0.75, bin_minutes=10.0):
    """Detrend + bin already-loaded (normalised) flux. Shared by the search, injection and null tests."""
    if len(t) < 500:
        raise ValueError(f"too few good cadences ({len(t)})")
    tf, ff, ef, trend, _ = lightcurve.detrend(t, f, e, window=window)
    tb, fb, eb = lightcurve.bin_lc(tf, ff, ef, minutes=bin_minutes) if bin_minutes else (tf, ff, ef)
    return dict(raw=(t, f, e), trend=trend, flat=(tf, ff, ef), binned=(tb, fb, eb), meta=meta)


def prepare(path, window=0.75, bin_minutes=10.0):
    t, f, e, meta = lightcurve.load_spoc(path)
    return prepare_arrays(t, f, e, meta, window, bin_minutes)


def search_prepared(d, backend="cpu", s_min=None):
    """BLS (S statistic) + per-star vetting on a prepared light curve -> result fields (as stored in the DB).
    Only the CPU path computes S; the GPU/MLX backends are exploratory and not used for v1 results."""
    from .version import stamp
    tb, fb, eb = d["binned"]
    tf, ff, _ = d["flat"]
    out = bls.search_s(tb, fb, tf, ff)
    st = bls.stats(tb, fb, eb, out["period"], out["duration"], out["t0"])
    row = {k: d["meta"].get(k) for k in ("tmag", "teff", "radius", "ra", "dec", "camera", "ccd")}
    row.update(n_points=int(len(d["raw"][0])), period=out["period"], t0=out["t0"],
               duration=out["duration"], sde=out["sde"], S=out["S"], S_scan=out["S_scan"], n_tr=out["n_tr"],
               cdpp_d=out["cdpp_d"], beta_d=out["beta_d"], sigma_w=out["sigma_w"], depth_box=out["depth_box"],
               s_min=s_min if s_min is not None else vet.s_star(d["meta"].get("sector"), d["meta"].get("camera"))[0])
    row.update(vet.vet(st, out, d["meta"], tb, fb, s_min=row["s_min"]))
    row.update(stamp())
    return row


def process_star(tic, sector, fname, url, backend="auto", window=0.75, bin_minutes=10.0):
    t_start = time.perf_counter()
    be = resolve_backend(backend)
    row = dict(tic=int(tic), sector=int(sector), backend=be, host=socket.gethostname(),
               processed_at=dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"))
    try:
        # Remote workers without the shared disk download into their own EXO_DATA_ROOT.
        path = targets.download_one(fname, url, config.lc_dir(sector))
        d = prepare(path, window, bin_minutes)
        row.update(search_prepared(d, be))
        row["status"] = "ok"
    except Exception as exc:
        row.update(status="error", error=f"{type(exc).__name__}: {exc}"[:500])
    row["runtime_s"] = time.perf_counter() - t_start
    for k, v in list(row.items()):  # plain Python types for SQLite
        if isinstance(v, (np.floating, np.integer)):
            row[k] = v.item()
    return row
