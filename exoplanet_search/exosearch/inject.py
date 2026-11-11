"""Injection-recovery through the full search + vetting chain (Scientist spec; configs/injection_v1.json).

Per host light curve and run: draw (P, Rp, b, T0), build a batman transit (circular, Claret quadratic LD,
a/R* from Kepler's third law with TIC M*, R*), multiply it into the normalised PDCSAP flux of every
cadence BEFORE quality masking and detrending, then run exactly the production code path:
quality mask -> biweight detrend -> bin -> BLS -> screens -> common-period list -> light-curve deep
checks (V-shape, secondary scan, Rp > 2 RJup) -> transit fit / density check. Pixel-level checks
(centroid, dip location, known EB, Gaia neighbours) cannot be exercised by a light-curve injection and
are counted as passed for on-target signals (hosts with catalogue EBs are excluded).

  python -m exosearch.inject --config configs/injection_v1.json --runs 0 --workers 4
"""
from __future__ import annotations

import argparse
import hashlib
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd

from . import batch, bls, config, nulltest, deepvet, density, holdout, labels, lightcurve, stellar, store, targets, task, vet
from .version import packages, stamp

R_EARTH_RSUN = 1 / 109.1
G, MSUN, RSUN, DAY = 6.674e-8, 1.989e33, 6.957e10, 86400.0
PRODUCTS = config.DATA_ROOT / "products"


def seed_for(base, sector, run, tic):
    return int.from_bytes(hashlib.sha256(f"{base}:{sector}:{run}:{tic}".encode()).digest()[:4], "big")


def select_hosts(sector, cfg, con=None):
    """Searched stars, stratified by Tmag x R*, without TOI hosts or catalogue EBs. Deterministic."""
    con = con or store.connect()
    df = store.frame(con, sector)
    df = df[(df.status == "ok") & df.radius.notna() & (df.radius > 0) & df.tmag.notna()]
    toi = labels.toi_table()
    eb = pd.read_csv(config.CATALOG_DIR / "tess_ebs_prsa2022.csv")
    df = df[~df.tic.isin(toi.tic) & ~df.tic.isin(eb.TIC)]
    te, re_ = cfg["strata"]["tmag_edges"], cfg["strata"]["rstar_edges"]
    df = df.assign(tbin=np.digitize(df.tmag, te[1:-1]), rbin=np.digitize(df.radius, re_[1:-1]))
    rng = np.random.default_rng(seed_for(cfg["seed"], sector, -1, 0))
    n = cfg["hosts_per_sector"]
    cells = [g for _, g in df.groupby(["tbin", "rbin"])]
    quota = n // len(cells)
    pick = [g.sample(min(quota, len(g)), random_state=int(rng.integers(2 ** 31))) for g in cells]
    chosen = pd.concat(pick)
    rest = df[~df.tic.isin(chosen.tic)]
    if len(chosen) < n and len(rest):
        chosen = pd.concat([chosen, rest.sample(min(n - len(chosen), len(rest)), random_state=int(rng.integers(2 ** 31)))])
    return chosen.sort_values("tic")


def common_periods(sector, con=None):
    d = vet.flag_common_periods(store.frame(con or store.connect(), sector))
    return np.sort(d.loc[d["flags"].fillna("").str.contains("common_period"), "period"].to_numpy(float))


def t14(P, aR, k, b):
    inc = np.arccos(b / aR)
    arg = np.sqrt(max((1 + k) ** 2 - b ** 2, 0)) / (aR * np.sin(inc))
    return P / np.pi * np.arcsin(min(arg, 1.0))


def cdpp_at(t, f, width, intransit):
    """CDPP_D at D = injected duration, by the same estimator the search uses, with injected transits masked."""
    return bls.noise_profile(t[~intransit], f[~intransit], [width])[2][float(width)][0]


def run_star(job):
    """All runs for one host light curve: load + quality-mask once, then one seeded injection per run."""
    t0 = time.perf_counter()
    try:
        t, f, e, q, meta = lightcurve.load_spoc_raw(job["path"])
        good = lightcurve.good_mask(f, e, q)
        med = np.nanmedian(f[good])
        lc = (t, f / med, e / med, good, meta)
    except Exception as exc:
        return [dict(job=job["key"], key=f"{job['key']}:{r}", tic=job["tic"], sector=job["sector"], run=r,
                     status="error", error=f"load: {type(exc).__name__}: {exc}"[:300]) for r in job["runs"]]
    load_s = time.perf_counter() - t0
    rows = []
    for r, sd in zip(job["runs"], job["seeds"]):
        row = run_one(job, lc, r, sd)
        row.update(job=job["key"], load_s=load_s)
        rows.append(row)
    return rows


def run_one(job, lc, run, seed):
    """One injection into one (already loaded) light curve -> result row."""
    t_start = time.perf_counter()
    cfg = job["cfg"]
    row = dict(key=f"{job['sector']}:{run}:{job['tic']}", tic=job["tic"], sectors=str(job["sector"]),
               sector=job["sector"], run=run, seed=seed, holdout=holdout.is_holdout(job["tic"]), **job["stamp"])
    rng = np.random.default_rng(seed)
    try:
        t, f, e, good, meta = lc
        tg = t[good]
        baseline = float(tg.max() - tg.min())
        R, M = job["rstar"], job["mstar"]
        u = stellar.quad_ld(job["teff"], job["logg"])
        P = float(np.exp(rng.uniform(np.log(0.5), np.log(baseline / 2))))
        rp = float(np.exp(rng.uniform(np.log(0.5), np.log(20.0))))
        b = float(rng.uniform(0.0, 1.0))
        T0 = float(tg.min() + rng.uniform(0, P))
        k = rp * R_EARTH_RSUN / R
        aR = float((G * M * MSUN * (P * DAY) ** 2 / (4 * np.pi ** 2)) ** (1 / 3) / (R * RSUN))
        row.update(tmag=job["tmag"], rstar=R, mstar=M, teff=job["teff"], logg=job["logg"], u1=u[0], u2=u[1],
                   P_inj=P, rp_inj=rp, b_inj=b, t0_inj=T0, k_inj=k, aR_inj=aR, baseline=baseline)
        if aR <= 1 + k:
            row.update(status="unphysical")
            return row
        mdl = density.model(t, T0, P, k, aR, b, u)
        dur = t14(P, aR, k, b)
        fi = f * mdl
        row.update(depth_inj=float(1 - mdl.min()), dur_inj=float(dur))
        # transits with >= half their expected cadences surviving the quality mask
        ph = ((tg - T0 + 0.5 * P) % P) - 0.5 * P
        intr = np.abs(ph) < dur / 2
        ep = np.round((tg[intr] - T0) / P).astype(int)
        cad = np.median(np.diff(tg))
        n_tr = int(sum(np.sum(ep == x) >= 0.5 * dur / cad for x in np.unique(ep))) if dur > 0 else 0
        row["n_tr"] = n_tr
        t_prep = time.perf_counter()
        d = task.prepare_arrays(tg, fi[good], e[good], meta)
        row["detrend_s"] = time.perf_counter() - t_prep
        tf, ff, ef = d["flat"]
        phf = ((tf - T0 + 0.5 * P) % P) - 0.5 * P
        w = max(dur, 1 / 48)
        c = cdpp_at(tf, ff, w, np.abs(phf) < w)
        row.update(cdpp_dur=c, snr_exp=float(row["depth_inj"] / c * np.sqrt(n_tr)) if c and c > 0 else np.nan)
        t_s = time.perf_counter()
        res = task.search_prepared(d)
        row["search_s"] = time.perf_counter() - t_s
        Pr, t0r = res["period"], res["t0"]
        dt0 = ((t0r - T0 + 0.5 * P) % P) - 0.5 * P
        rc = cfg["recovery"]
        matched = abs(Pr / P - 1) < rc["period_tol"] and abs(dt0) < rc["t0_tol_durations"] * max(dur, 1e-3)
        alias = None
        for a in rc["aliases"]:
            if abs(Pr / (P * a) - 1) < rc["period_tol"]:
                alias = {2: "2x", 0.5: "1/2x", 3: "3x"}.get(a, "1/3x")
        failed = [x for x in (res["flags"] or "").split(",") if x]
        if len(job["common"]) and np.any(np.abs(job["common"] / Pr - 1) < vet.THRESH["common_period_tol"]):
            failed.append("common_period")
        detected = bool(matched and res["S"] >= res["s_min"])
        warn = []
        fl0 = set(failed)
        for c in nulltest.SCREEN_CHECKS + ["common_period"]:
            row[f"fail_{c}"] = c in fl0
        row["checks_logged"] = bool(matched or res["S"] >= nulltest.S_FLOOR)
        if row["checks_logged"]:                      # deep light-curve checks + density, as for real candidates
            sm = deepvet.shape_metrics(tf, ff, ef, Pr, t0r)
            dv = dict(sm, rp_rjup=float(np.sqrt(max(res["depth"], 0)) * R / deepvet.R_JUP_RSUN))
            fit = density.fit(tf, ff, ef, Pr, sm["t0_fine"], sm["duration_fine"], res["depth"], u, job["rho"] or 1.0)
            if fit:
                dv.update(fit)
                dv["rho_z"] = density.density_z(fit["rho_transit"], fit["rho_transit_logerr"], job["rho"], job["rho_err"])
            dv.update(ruwe=job["ruwe"], gaia_nss=job["nss"], diff_peak_snr=None, n_epochs=0)
            fl = deepvet.flags_for(dv)
            for c in nulltest.DEEP_CHECKS:
                row[f"fail_{c}"] = c in fl
            failed += [x for x in fl if x in deepvet.REJECT]
            warn = [x for x in fl if x not in deepvet.REJECT and x != "centroid_inconclusive"]
            row.update(v_ratio=dv.get("v_ratio"), sec_scan_sigma=dv.get("sec_scan_sigma"), rho_z=dv.get("rho_z"),
                       fit_ar=dv.get("fit_ar"), fit_b=dv.get("fit_b"))
        row.update(status="ok", P_rec=Pr, t0_rec=t0r, dur_rec=res["duration"], depth_rec=res["depth"],
                   sde=res["sde"], snr=res["snr"], S=res["S"], S_scan=res["S_scan"], n_tr_rec=res["n_tr"],
                   cdpp_d_rec=res["cdpp_d"], beta_d=res["beta_d"], sigma_w=res["sigma_w"], s_min=res["s_min"],
                   camera=res.get("camera"), n_transits_rec=res["n_transits"], matched=bool(matched),
                   alias=alias, detected=detected, passed_all=bool(detected and not failed),
                   failed_checks=",".join(failed), warnings=",".join(warn))
    except Exception as exc:
        row.update(status="error", error=f"{type(exc).__name__}: {exc}"[:300])
    row["runtime_s"] = time.perf_counter() - t_start
    return row


def build_jobs(cfg, runs, sectors=None, log=print):
    con = store.connect()
    st = stamp()
    jobs = []
    for s in sectors or cfg["sectors"]:
        hosts = select_hosts(s, cfg, con)
        tp = stellar.tic_params(hosts.tic.tolist(), log=log)
        gd = stellar.gaia_dr3(tp.GAIA.tolist(), log=log)
        man = targets.fetch_manifest(s).set_index("tic")
        com = common_periods(s, con)
        for h in hosts.itertuples():
            tic = int(h.tic)
            if tic not in tp.index or tic not in man.index:
                continue
            x = tp.loc[tic]
            R = float(x.rad) if x.rad == x.rad else float(h.radius)
            if x.mass == x.mass and x.mass > 0:
                M = float(x.mass)
            elif x.logg == x.logg:
                M = float(10 ** (x.logg - 4.438) * R ** 2)
            else:
                continue
            if not (R > 0 and 0.05 < M < 10):
                continue
            path = batch.lc_path(s, man.loc[tic].fname, download=False)
            if path is None:
                continue
            rho, erho = stellar.tic_rho(x.to_dict())
            g = str(x.GAIA)
            gr = gd.loc[int(g)] if g.isdigit() and int(g) in gd.index else None
            jobs.append(dict(key=f"{s}:{tic}:{runs[0]}-{runs[-1]}", tic=tic, sector=s, runs=list(runs),
                                 seeds=[seed_for(cfg["seed"], s, r, tic) for r in runs], path=str(path), rstar=R, mstar=M, teff=float(x.Teff) if x.Teff == x.Teff else None,
                                 logg=float(x.logg) if x.logg == x.logg else None, tmag=float(h.tmag), rho=rho, rho_err=erho,
                                 ruwe=None if gr is None or pd.isna(gr.ruwe) else float(gr.ruwe),
                                 nss=None if gr is None or pd.isna(gr.non_single_star) else int(gr.non_single_star),
                                 common=com, cfg=cfg, stamp=st))
    return jobs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="configs/injection_v1.json")
    ap.add_argument("--runs", default="0", help="e.g. 0 or 0-4")
    ap.add_argument("--sectors", default=None)
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--limit", type=int, default=0, help="first N jobs only (timing)")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    cfg = json.loads(Path(a.config).read_text())
    lo, _, hi = a.runs.partition("-")
    runs = list(range(int(lo), int(hi or lo) + 1))
    sectors = [int(x) for x in a.sectors.split(",")] if a.sectors else None
    st = stamp()
    out = Path(a.out or PRODUCTS / f"v{st['pipeline_version']}" / "injections" / cfg["name"])
    out.mkdir(parents=True, exist_ok=True)
    (out / "provenance.json").write_text(json.dumps(dict(st, packages=packages(), config=cfg), indent=1))
    jobs = build_jobs(cfg, runs, sectors)
    if a.limit:
        jobs = jobs[:a.limit]
    with batch.claim(a.workers):
        df = batch.pmap(run_star, jobs, out, a.workers, key="job")
    ok = df[df.status == "ok"]
    print(f"{len(df)} rows, ok {len(ok)}, detected {ok.detected.mean():.3f}, passed_all {ok.passed_all.mean():.3f}")


if __name__ == "__main__":
    main()
