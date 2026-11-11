"""v1.1 per-star processing for every population: real search, flipped (inverted) null, scrambled nulls, injections.

One job = one star with all of its tasks, so a star's real, null and injection runs always run on the same host.
Chain per task: load SPOC 2-min light curve -> (null transform | transit injection) -> quality mask -> prewhitening
above 5 c/d with the comb guard -> biweight detrend (no-fastmath kernels) -> 10-min bins -> v1.1 search
(exosearch.v11.search) -> screens (exosearch.v11.vet11) -> light-curve deep checks and density fit when S >= 6 or the
screens pass (as v1.0.0). Injection draws reuse the v1.0.0 seeds and formulas (exosearch.inject), so v1.1 injections
are the same signals in the same light curves; a registered subset is drawn from "true" stellar densities that differ
from TIC (configs/injection_v1.1.json, "wrong_density") while the search keeps the TIC density for the cap.
"""
from __future__ import annotations

import hashlib
import json
import time

import numpy as np

from .. import batch, deepvet, density, lightcurve, nulltest, stellar
from ..inject import G, MSUN, RSUN, DAY, R_EARTH_RSUN, t14, cdpp_at
from . import kernels, scramble, search, settings, vet11
from .prewhiten import prewhiten

KEEP_META = ("tmag", "teff", "radius", "ra", "dec", "camera", "ccd")


def load(job, download=True):
    p = batch.lc_path(job["sector"], job["fname"], download=download)
    if p is None:
        raise FileNotFoundError("light curve not available")
    t, f, e, q, meta = lightcurve.load_spoc_raw(p)
    good = lightcurve.good_mask(f, e, q)
    med = np.nanmedian(f[good])
    dumps = vet11.dump_times(t, q, job["cfg"]["momentum_dump"]["quality_bit"])
    return t, f / med, e / med, good, meta, dumps


def prepare(t, f, e, cfg, prewhiten_on=True, comb_guard=None, gap_trim=0.3):
    kernels.install()
    if len(t) < 500:
        raise ValueError(f"too few good cadences ({len(t)})")
    if prewhiten_on:
        f, pwlog = prewhiten(t, f, cfg["prewhiten"], comb_guard)
    else:
        pwlog = dict(prewhiten_n=0, prewhiten_comb=0, prewhiten_terms="[]", prewhiten_comb_freqs="[]", prewhiten_off=True)
    tf, ff, ef, _, _ = lightcurve.detrend(t, f, e, window=0.75, gap_trim=gap_trim)
    tb, fb, eb = lightcurve.bin_lc(tf, ff, ef, minutes=10.0)
    return dict(flat=(tf, ff, ef), binned=(tb, fb, eb), pw=pwlog)


def analyse(d, job, meta, dumps, primary, cap=True, deep=True):
    cfg = job["cfg"]
    tb, fb, eb = d["binned"]
    tf, ff, ef = d["flat"]
    rho = job.get("rho")
    best, sens = search.search(tb, fb, eb, tf, ff, rho=rho, dumps=dumps, cfg=cfg, primary=primary,
                               prior=job.get("beta_prior"), cap=cap)
    from astropy.timeseries import BoxLeastSquares
    st = BoxLeastSquares(tb, fb, dy=eb).compute_stats(best["P_grid"], best["duration"], best["t0"])
    row = dict(best)
    row.update(vet11.screens(st, best, dict(meta, radius=job.get("radius") or meta.get("radius")), cfg))
    for o, v in sens.items():
        for k in ("S", "period", "P_grid", "duration", "t0", "P_grid_halfwidth"):
            row[f"{k}_opt_{o}"] = None if v is None else v[k]
    row["checks_logged"] = bool(best["S"] >= cfg["s_floor_checks"] or row["passed_screens"])
    row["deep_flags"] = None
    if deep and row["checks_logged"]:
        P, t0 = best["P_grid"], best["t0"]
        sm = deepvet.shape_metrics(tf, ff, ef, P, t0)
        R = job.get("radius") or meta.get("radius") or np.nan
        dv = dict(sm, rp_rjup=float(np.sqrt(max(row["depth"], 0)) * R / deepvet.R_JUP_RSUN) if row["depth"] == row["depth"] else None,
                  diff_peak_snr=None, n_epochs=0, ruwe=job.get("ruwe"), gaia_nss=job.get("nss"))
        u = stellar.quad_ld(job.get("teff"), job.get("logg"))
        fit = density.fit(tf, ff, ef, P, sm["t0_fine"], sm["duration_fine"], max(row["depth"], 1e-6) if row["depth"] == row["depth"] else 1e-6,
                          u, rho or 1.0)
        if fit:
            dv["rho_z"] = density.density_z(fit["rho_transit"], fit["rho_transit_logerr"], rho, job.get("rho_err"))
        fl = deepvet.flags_for(dv)
        row.update(v_ratio=dv.get("v_ratio"), v_sigma=dv.get("v_sigma"), sec_scan_sigma=dv.get("sec_scan_sigma"),
                   rp_rjup=dv.get("rp_rjup"), rho_z=dv.get("rho_z"), deep_flags=",".join(fl))
    rej = set((row["flags"] or "").split(",")) - {""}
    rej |= set((row["deep_flags"] or "").split(",")) & deepvet.REJECT
    row["passed_lc_checks"] = bool(not rej)          # every light-curve check; common_period is applied at analysis
    for k in ("cdpp_by_d", "beta_by_d"):
        row[k] = best[k]
    return row


# ---------------------------------------------------------------- injections
def dens_seed(key):
    return int.from_bytes(hashlib.sha256(f"v1.1-density:{key}".encode()).digest()[:4], "big")


def draw_density(job, run, cfg_inj, force=False):
    """-> (class, rho_true / rho_TIC, R_true / R_TIC). Registered wrong-density subset (>= 20% of injections);
    force=True (validation only) always draws a wrong density."""
    wd = cfg_inj.get("wrong_density")
    if not wd:
        return "tic", 1.0, 1.0
    rng = np.random.default_rng(dens_seed(f"{job['sector']}:{run}:{job['tic']}"))
    if rng.uniform() >= (1.0 if force else wd["fraction"]):
        return "tic", 1.0, 1.0
    cls = wd["classes"][int(rng.integers(len(wd["classes"])))]
    if cls == "binary":           # unresolved companion: TIC radius inflated by up to sqrt(2); true star denser
        infl = float(rng.uniform(*wd["binary_radius_inflation"]))
        return cls, infl ** 3, 1.0 / infl
    if cls == "subgiant":         # misclassified subgiant: true density 3-10x below TIC
        return cls, 1.0 / float(rng.uniform(*wd["subgiant_density_factor"])), 1.0
    lo, hi = wd["loguniform_density_factor"]                 # spread: log-uniform, both directions, beyond 27x
    fac = float(np.exp(rng.uniform(np.log(lo), np.log(hi))))
    return cls, fac, fac ** (-1.0 / 3.0)


def inject_one(job, lc, run, seed, primary, prewhiten_on=True, comb_guard=None, cap=True, wrong_density=True,
               draws=None, force_density=False):
    t, f, e, good, meta, dumps = lc
    cfg_inj = job["cfg_inj"]
    rng = np.random.default_rng(seed)
    tg = t[good]
    baseline = float(tg.max() - tg.min())
    R, M = job["rstar"], job["mstar"]
    u = stellar.quad_ld(job.get("teff"), job.get("logg"))
    P = float(np.exp(rng.uniform(np.log(0.5), np.log(baseline / 2))))          # same draw order as v1.0.0
    rp = float(np.exp(rng.uniform(np.log(0.5), np.log(20.0))))
    b = float(rng.uniform(0.0, 1.0))
    T0 = float(tg.min() + rng.uniform(0, P))
    if draws:                                                                  # validation overrides (USP / short-D sets)
        P = draws.get("P", P)
        T0 = float(tg.min() + rng.uniform(0, P)) if "P" in draws else T0
        rp, b = draws.get("rp", rp), draws.get("b", b)
    cls, rho_fac, r_fac = draw_density(job, run, cfg_inj, force_density) if wrong_density else ("tic", 1.0, 1.0)
    R_true = R * r_fac
    aR_tic = float((G * M * MSUN * (P * DAY) ** 2 / (4 * np.pi ** 2)) ** (1 / 3) / (R * RSUN))
    aR = aR_tic * rho_fac ** (1.0 / 3.0)                                       # a/R* scales as rho^(1/3)
    k = rp * R_EARTH_RSUN / R_true
    row = dict(P_inj=P, rp_inj=rp, b_inj=b, t0_inj=T0, k_inj=k, aR_inj=aR, aR_tic=aR_tic, dens_class=cls,
               dens_factor=rho_fac, rstar=R, mstar=M, baseline=baseline)
    if aR <= 1 + k:
        row["status"] = "unphysical"
        return row
    mdl = density.model(t, T0, P, k, aR, b, u)
    dur = t14(P, aR, k, b)
    fi = f * mdl
    row.update(depth_inj=float(1 - mdl.min()), dur_inj=float(dur))
    ph = ((tg - T0 + 0.5 * P) % P) - 0.5 * P
    intr = np.abs(ph) < dur / 2
    ep = np.round((tg[intr] - T0) / P).astype(int)
    cad = np.median(np.diff(tg))
    row["n_tr_inj"] = int(sum(np.sum(ep == x) >= 0.5 * dur / cad for x in np.unique(ep))) if dur > 0 else 0
    d = prepare(tg, fi[good], e[good], job["cfg"], prewhiten_on, comb_guard)
    tf, ff, _ = d["flat"]
    phf = ((tf - T0 + 0.5 * P) % P) - 0.5 * P
    w = max(dur, 1 / 48)
    c = cdpp_at(tf, ff, w, np.abs(phf) < w)
    row.update(cdpp_dur=c, snr_exp=float(row["depth_inj"] / c * np.sqrt(row["n_tr_inj"])) if c and c > 0 else np.nan)
    res = analyse(d, job, meta, dumps, primary, cap=cap)
    rc = cfg_inj["recovery"]

    def match(Pr, t0r):
        if Pr is None or t0r is None:
            return False, None
        dt0 = ((t0r - T0 + 0.5 * P) % P) - 0.5 * P
        m = abs(Pr / P - 1) < rc["period_tol"] and abs(dt0) < rc["t0_tol_durations"] * max(dur, 1e-3)
        al = None
        for a in rc["aliases"]:
            if abs(Pr / (P * a) - 1) < rc["period_tol"]:
                al = {2: "2x", 0.5: "1/2x", 3: "3x"}.get(a, "1/3x")
        return bool(m), al

    row.update(res)
    row.update({f"pw_{k}": v for k, v in d["pw"].items()})
    row["matched"], row["alias"] = match(res["P_grid"], res["t0"])
    for o in settings.BETA_OPTIONS:
        row[f"matched_opt_{o}"] = match(res.get(f"P_grid_opt_{o}"), res.get(f"t0_opt_{o}"))[0]
    row["depth_ratio"] = float(res["depth"] / row["depth_inj"]) if row["depth_inj"] > 0 else np.nan
    row["status"] = "ok"
    return row


# ---------------------------------------------------------------- job runner
def run_job(job):
    """All tasks of one star -> list of rows. job['tasks']: [{'mode': 'real'|'invert'|'scramble'|'inject', ...}]."""
    kernels.install()                                  # before any detrend (the inverted null's flare mask detrends too)
    t_load = time.perf_counter()
    base = dict(job=job["key"], tic=job["tic"], sector=job["sector"], holdout=job["holdout"], **job["prov"])
    try:
        lc = load(job)
    except Exception as exc:
        return [dict(base, key=f"{job['key']}:{tk['tag']}", task=tk["tag"], mode=tk["mode"], status="error",
                     error=f"load: {type(exc).__name__}: {exc}"[:300]) for tk in job["tasks"]]
    load_s = time.perf_counter() - t_load
    t, f, e, good, meta, dumps = lc
    rows = []
    for tk in job["tasks"]:
        t0 = time.perf_counter()
        row = dict(base, key=f"{job['key']}:{tk['tag']}", task=tk["tag"], mode=tk["mode"], realization=tk.get("realization"),
                   run=tk.get("run"), seed=tk.get("seed"), load_s=load_s, n_dumps=int(len(dumps)),
                   **{k: meta.get(k) for k in KEEP_META})
        try:
            pw_on = tk.get("prewhiten", True)
            if tk["mode"] == "inject":
                row.update(inject_one(job, lc, tk["run"], tk["seed"], job["primary"], pw_on, tk.get("comb_guard"),
                                      tk.get("cap", True), tk.get("wrong_density", True), tk.get("draws"),
                                      tk.get("force_density", False)))
            else:
                tg, fg, eg = t[good], f[good], e[good]
                if tk["mode"] == "invert":
                    tg, fg, eg = nulltest.transform(tk["mode"], tg, fg, eg, meta, tk["seed"])
                elif tk["mode"] == "scramble":
                    tg, fg, eg, scr_info = scramble.segment_rotate(tg, fg, eg, tk["seed"], job["cfg"])
                    row.update(scr_info)
                elif tk["mode"] != "real":
                    raise ValueError(tk["mode"])
                gap_trim = 0.0 if tk["mode"] == "scramble" else 0.3
                d = prepare(tg, fg, eg, job["cfg"], pw_on, tk.get("comb_guard"), gap_trim=gap_trim)
                row.update(analyse(d, job, meta, dumps, job["primary"], cap=tk.get("cap", True)))
                row.update({f"pw_{k}": v for k, v in d["pw"].items()})
                row["status"] = "ok"
        except Exception as exc:
            row.update(status="error", error=f"{type(exc).__name__}: {exc}"[:300])
        row["runtime_s"] = time.perf_counter() - t0
        for k, v in list(row.items()):
            if isinstance(v, (np.floating, np.integer, np.bool_)):
                row[k] = v.item()
        rows.append(row)
    return rows
