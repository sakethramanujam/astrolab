"""Labelled signals for tuning and testing the vetting checks.

planet : the star hosts a TOI dispositioned KP or CP and our best BLS peak matches its period
         (ratio 1, 2, 1/2 within 1%)
fp     : the best peak matches a TOI dispositioned FP or FA, or the Prsa et al. (2022) TESS EB period
PCs/APCs are unlabelled. Each row carries `holdout` (see holdout.py).
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from . import config, holdout, store, vet

RATIOS = (1.0, 2.0, 0.5)


def _match(p, ref, tol=0.01):
    return ref == ref and ref > 0 and any(abs(p * k / ref - 1) < tol for k in RATIOS)


def toi_table():
    t = pd.read_csv(config.CATALOG_DIR / "exofop_toi.csv", comment="#", low_memory=False)
    return t.rename(columns={"TIC ID": "tic", "TOI": "toi", "Period (days)": "period",
                             "TFOPWG Disposition": "disp", "Epoch (BJD)": "epoch_bjd",
                             "Duration (hours)": "dur_h", "Depth (ppm)": "depth_ppm",
                             "Planet Radius (R_Earth)": "rp_toi"})


def labelled(sector, con=None):
    con = con or store.connect()
    df = vet.flag_common_periods(store.frame(con, sector))
    df = df[df.status == "ok"]
    toi = toi_table()
    eb = pd.read_csv(config.CATALOG_DIR / "tess_ebs_prsa2022.csv")
    eb["Per"] = pd.to_numeric(eb.Per, errors="coerce")
    toi = toi[toi.tic.isin(df.tic)]
    eb = eb[eb.TIC.isin(df.tic)]
    rows = []
    for r in df[df.tic.isin(set(toi.tic) | set(eb.TIC))].itertuples():
        lab, src = None, None
        for x in toi[toi.tic == r.tic].itertuples():
            if _match(r.period, x.period):
                if x.disp in ("KP", "CP"):
                    lab, src = "planet", f"TOI {x.toi} {x.disp}"
                    break
                if x.disp in ("FP", "FA") and lab is None:
                    lab, src = "fp", f"TOI {x.toi} {x.disp}"
        if lab is None:
            for x in eb[eb.TIC == r.tic].itertuples():
                if _match(r.period, x.Per):
                    lab, src = "fp", "TESS EB"
        if lab:
            rows.append(dict(tic=int(r.tic), sector=int(sector), label=lab, source=src))
    out = pd.DataFrame(rows)
    out["holdout"] = out.tic.map(holdout.is_holdout)
    return out.merge(df, on=["tic", "sector"], how="left")


def lc_checks(job):
    """Light-curve-level deep checks (shape, secondary scan, density, RUWE) for one labelled signal."""
    from types import SimpleNamespace
    from . import deepvet, task
    r = SimpleNamespace(**job)
    out = dict(key=job["key"], tic=r.tic, sector=r.sector)
    try:
        d = task.prepare(job["path"])
        t, f, e = d["flat"]
        sm = deepvet.shape_metrics(t, f, e, float(r.period), float(r.t0))
        out.update(sm)
        out.update(deepvet.stellar_checks(r, job["path"], sm["t0_fine"], sm["duration_fine"]))
        if r.radius == r.radius and r.radius and r.depth > 0:
            out["rp_rjup"] = float(np.sqrt(r.depth) * r.radius / deepvet.R_JUP_RSUN)
        out["status"] = "ok"
    except Exception as exc:
        out.update(status="error", error=f"{type(exc).__name__}: {exc}"[:300])
    return out


def label_peak(tic, period, toi, eb):
    """(label, source) for a peak on `tic`, or (None, None)."""
    lab, src = None, None
    for x in toi[toi.tic == tic].itertuples():
        if _match(period, x.period):
            if x.disp in ("KP", "CP"):
                return "planet", f"TOI {x.toi} {x.disp}"
            if x.disp in ("FP", "FA") and lab is None:
                lab, src = "fp", f"TOI {x.toi} {x.disp}"
    if lab is None:
        for x in eb[eb.TIC == tic].itertuples():
            if _match(period, x.Per):
                return "fp", "TESS EB"
    return lab, src


def v1_job(job):
    """Frozen v1 search + every check (screens, light-curve deep checks, density, RUWE, and the pixel-level
    centroid / Gaia / EB-catalogue tests) on one TOI or EB host."""
    from types import SimpleNamespace
    from . import deepvet, task
    from .version import stamp
    out = dict(key=job["key"], tic=job["tic"], sector=job["sector"], holdout=holdout.is_holdout(job["tic"]), **stamp())
    try:
        d = task.prepare(job["path"])
        res = task.search_prepared(d)
        out.update({k: res.get(k) for k in ("period", "t0", "duration", "depth", "S", "S_scan", "n_tr", "cdpp_d",
                                              "beta_d", "sigma_w", "sde", "snr", "flags", "passed", "rp_earth",
                                              "radius", "tmag", "ra", "dec", "camera", "s_min", "odd_even_sigma",
                                              "secondary_sigma", "dur_ratio", "n_transits")})
        lab, src = label_peak(job["tic"], res["period"], job["toi"], job["eb"])
        out.update(label=lab, source=src)
        if lab is None:
            out["status"] = "unlabelled"
            return out
        r = SimpleNamespace(tic=job["tic"], sector=job["sector"], period=res["period"], t0=res["t0"],
                            duration=res["duration"], depth=res["depth"], radius=res["radius"], ra=res["ra"],
                            dec=res["dec"], tmag=res["tmag"])
        dv = deepvet.vet_star(r, job["fname"])
        out.update({k: v for k, v in dv.items() if k not in ("tic", "sector", "updated_at")})
        out["status"] = "ok"
    except Exception as exc:
        out.update(status="error", error=f"{type(exc).__name__}: {exc}"[:300])
    return out


def main_v1(sectors, workers):
    from . import batch, targets
    from .version import stamp
    st = stamp()
    out = config.DATA_ROOT / "products" / f"v{st['pipeline_version']}" / "labelled_v1"
    toi = toi_table()
    eb = pd.read_csv(config.CATALOG_DIR / "tess_ebs_prsa2022.csv")
    eb["Per"] = pd.to_numeric(eb.Per, errors="coerce")
    jobs = []
    for s in sectors:
        man = targets.fetch_manifest(s).set_index("tic")
        hosts = sorted((set(toi.tic) | set(eb.TIC)) & set(man.index))
        for tic in hosts:
            p = batch.lc_path(s, man.loc[tic].fname)
            jobs.append(dict(key=f"{s}:{tic}", tic=int(tic), sector=s, path=str(p), fname=man.loc[tic].fname,
                             toi=toi[toi.tic == tic], eb=eb[eb.TIC == tic]))
    with batch.claim(workers):
        res = batch.pmap(v1_job, jobs, out / "parts", workers)
    res.to_parquet(out / "labelled_v1.parquet", index=False)
    print(res.groupby(["sector", "label"]).size().to_string())


def main():
    import argparse
    from . import batch, targets
    from .version import stamp
    ap = argparse.ArgumentParser()
    ap.add_argument("--sectors", default="70,106,107")
    ap.add_argument("--workers", type=int, default=12)
    ap.add_argument("--v1", action="store_true", help="re-search TOI/EB hosts with the frozen v1 search + all checks")
    a = ap.parse_args()
    if a.v1:
        return main_v1([int(x) for x in a.sectors.split(",")], a.workers)
    st = stamp()
    out = config.DATA_ROOT / "products" / f"v{st['pipeline_version']}" / "labelled"
    out.mkdir(parents=True, exist_ok=True)
    con = store.connect()
    L = pd.concat([labelled(int(s), con) for s in a.sectors.split(",")], ignore_index=True)
    jobs = []
    for s, g in L.groupby("sector"):
        man = targets.fetch_manifest(int(s)).set_index("tic")
        for r in g.itertuples():
            p = batch.lc_path(int(s), man.loc[r.tic].fname)
            jobs.append(dict(key=f"{s}:{r.tic}", tic=int(r.tic), sector=int(s), period=r.period, t0=r.t0,
                             duration=r.duration, depth=r.depth, radius=r.radius, path=str(p)))
    L.to_parquet(out / "labelled_signals.parquet", index=False)
    with batch.claim(a.workers):
        res = batch.pmap(lc_checks, jobs, out / "lc_checks", a.workers)
    print(res.status.value_counts().to_dict())


if __name__ == "__main__":
    main()
