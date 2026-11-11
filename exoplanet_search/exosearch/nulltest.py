"""Noise-only searches for the empirical false-alarm rate (Kepler DR25-style).

invert   : mask flares (upward outliers of the detrended flux, +-1 cadence), then flux -> 2 - flux.
           Real dips become brightenings that the search cannot fit, so whatever the frozen search
           still finds is noise/systematics. This is the main single-sector false-alarm test.
scramble : cut the quality-masked light curve into ~1-day blocks of cadences and permute them
           (seeded), keeping the time stamps. Destroys periodicity but keeps local noise; with one
           sector there are only ~25 blocks, so treat this as limited until multi-sector searches.
Each star gets the identical search + screens (+ light-curve deep checks for screen survivors).
The per-star max-S table is the product; analysis/ fits its tail (GPD / Gumbel).

  python -m exosearch.nulltest --mode invert --sectors 70,106 --workers 12
"""
from __future__ import annotations

import argparse
import hashlib
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd

from . import batch, config, deepvet, density, holdout, lightcurve, stellar, store, targets, task, vet
from .version import packages, stamp

BLOCK_D = 1.0
S_FLOOR = 6.0      # every null peak with S >= S_FLOOR gets every check logged (for vetting effectiveness E)
SCREEN_CHECKS = ["few_transits", "odd_even", "secondary", "too_deep", "too_long", "variable", "rp_too_large", "no_dip"]
DEEP_CHECKS = ["v_shape", "secondary_scan", "rp_gt_2rjup", "density_mismatch", "ruwe_high", "gaia_nss"]


def seed_for(mode, sector, tic, realization=0, base=20261110):
    return int.from_bytes(hashlib.sha256(f"{base}:{mode}:{realization}:{sector}:{tic}".encode()).digest()[:4], "big")


def flare_mask(t, f, e, meta, clip=3.0):
    """True for cadences to keep: drops upward outliers (> clip sigma above the biweight trend) and neighbours."""
    d = task.prepare_arrays(t, f, e, meta)
    trend = d["trend"]
    r = f / trend - 1
    ok = np.isfinite(r)
    s = 1.4826 * np.nanmedian(np.abs(r[ok] - np.nanmedian(r[ok])))
    up = ok & (r > clip * max(s, 1e-6))
    up = up | np.r_[up[1:], False] | np.r_[False, up[:-1]]
    return ~up


def transform(mode, t, f, e, meta, seed):
    if mode == "invert":
        keep = flare_mask(t, f, e, meta)
        return t[keep], 2.0 - f[keep], e[keep]
    if mode == "scramble":
        rng = np.random.default_rng(seed)
        cad = float(np.median(np.diff(t)))
        n = max(int(BLOCK_D / cad), 10)
        blocks = [np.arange(i, min(i + n, t.size)) for i in range(0, t.size, n)]
        order = rng.permutation(len(blocks))
        idx = np.concatenate([blocks[i] for i in order])
        return t, f[idx], e[idx]
    raise ValueError(mode)


def run_one(job):
    t_start = time.perf_counter()
    row = dict(key=job["key"], tic=job["tic"], sector=job["sector"], mode=job["mode"], seed=job["seed"],
               realization=job["realization"],
               holdout=holdout.is_holdout(job["tic"]), **job["stamp"])
    try:
        t, f, e, meta = lightcurve.load_spoc(job["path"])
        row.update(camera=meta.get("camera"), ccd=meta.get("ccd"), tmag=meta.get("tmag"))
        t2, f2, e2 = transform(job["mode"], t, f, e, meta, job["seed"])
        d = task.prepare_arrays(t2, f2, e2, meta)
        res = task.search_prepared(d)
        keep = ("period", "t0", "duration", "S", "S_scan", "n_tr", "cdpp_d", "beta_d", "sigma_w", "depth_box",
                "sde", "snr", "depth", "n_transits", "rp_earth", "flags", "passed", "s_min", "radius")
        row.update({k: res.get(k) for k in keep})
        row["deep_flags"] = None
        fl = set((res["flags"] or "").split(","))
        for c in SCREEN_CHECKS:
            row[f"fail_{c}"] = c in fl
        row["checks_logged"] = bool(res["S"] >= S_FLOOR or res["passed"])
        if row["checks_logged"]:
            tf, ff, ef = d["flat"]
            sm = deepvet.shape_metrics(tf, ff, ef, res["period"], res["t0"])
            dv = dict(sm, rp_rjup=float(np.sqrt(max(res["depth"], 0)) * (res["radius"] or np.nan) / deepvet.R_JUP_RSUN),
                      diff_peak_snr=None, n_epochs=0, ruwe=job.get("ruwe"), gaia_nss=job.get("nss"))
            u = stellar.quad_ld(job.get("teff"), job.get("logg"))
            fit = density.fit(tf, ff, ef, res["period"], sm["t0_fine"], sm["duration_fine"], res["depth"], u, job.get("rho") or 1.0)
            if fit:
                dv["rho_z"] = density.density_z(fit["rho_transit"], fit["rho_transit_logerr"], job.get("rho"), job.get("rho_err"))
            row["rho_z"] = dv.get("rho_z")
            row.update(v_ratio=dv.get("v_ratio"), v_sigma=dv.get("v_sigma"), sec_scan_sigma=dv.get("sec_scan_sigma"),
                       rp_rjup=dv.get("rp_rjup"))
            dfl = deepvet.flags_for(dv)
            row["deep_flags"] = ",".join(dfl)
            for c in DEEP_CHECKS:
                row[f"fail_{c}"] = c in dfl
        row["status"] = "ok"
    except Exception as exc:
        row.update(status="error", error=f"{type(exc).__name__}: {exc}"[:300])
    row["runtime_s"] = time.perf_counter() - t_start
    return row


def finalize(df):
    """Common-period flag among the null peaks themselves, then pass/fail through all light-curve vetting."""
    out = []
    for (s, m), g in df[df.status == "ok"].groupby(["sector", "mode"]):
        g = vet.flag_common_periods(g.assign(status="ok"))
        g["fail_common_period"] = g["flags"].fillna("").str.contains("common_period")
        dfl = g.deep_flags.fillna("").str.split(",")
        g["passed_all"] = (g.passed == 1) & dfl.map(lambda x: not (set(x) & deepvet.REJECT))
        out.append(g)
    return pd.concat(out, ignore_index=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=["invert", "scramble"], required=True)
    ap.add_argument("--sectors", default="70,106")
    ap.add_argument("--workers", type=int, default=12)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--realization", type=int, default=0, help="scramble realization (own seed); invert has one")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    st = stamp()
    tag = a.mode if a.mode == "invert" else f"scramble_r{a.realization}"
    out = Path(a.out) if a.out else config.DATA_ROOT / "products" / f"v{st['pipeline_version']}" / "nulltests" / tag
    out.mkdir(parents=True, exist_ok=True)
    (out / "provenance.json").write_text(json.dumps(dict(st, packages=packages(), mode=a.mode, realization=a.realization,
                                                         block_d=BLOCK_D, s_floor=S_FLOOR), indent=1))
    con = store.connect()
    jobs = []
    for s in [int(x) for x in a.sectors.split(",")]:
        df = store.frame(con, s)
        df = df[df.status == "ok"]
        man = targets.fetch_manifest(s).set_index("tic")
        tp = stellar.tic_params(df.tic.tolist())
        gd = stellar.gaia_dr3(tp.GAIA.tolist())
        for tic in df.tic:
            tic = int(tic)
            if tic not in man.index:
                continue
            p = batch.lc_path(s, man.loc[tic].fname, download=False)
            if p is None:
                continue
            x = tp.loc[tic].to_dict() if tic in tp.index else {}
            rho, erho = stellar.tic_rho(x) if x else (None, None)
            g = str(x.get("GAIA", ""))
            gr = gd.loc[int(g)] if g.isdigit() and int(g) in gd.index else None
            jobs.append(dict(key=f"{tag}:{s}:{tic}", tic=tic, sector=s, mode=a.mode, realization=a.realization,
                             seed=seed_for(a.mode, s, tic, a.realization),
                             path=str(p), stamp=st, rho=rho, rho_err=erho, teff=x.get("Teff"), logg=x.get("logg"),
                             ruwe=None if gr is None or pd.isna(gr.ruwe) else float(gr.ruwe),
                             nss=None if gr is None or pd.isna(gr.non_single_star) else int(gr.non_single_star)))
    if a.limit:
        jobs = jobs[:a.limit]
    with batch.claim(a.workers):
        df = batch.pmap(run_one, jobs, out / "parts", a.workers)
    fin = finalize(df)
    fin.to_parquet(out / f"null_{tag}_maxS.parquet", index=False)
    print(fin.groupby("sector").agg(n=("S", "size"), passed=("passed", "sum"), passed_all=("passed_all", "sum")).to_string())


if __name__ == "__main__":
    main()
