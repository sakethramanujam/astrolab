"""Interim analysis layer over the frozen v1.0.0 search (search code, config_hash and stored S unchanged).

1. Star-quality cut (configs/star_quality_v1.json, frozen; registered after unblinding) applied identically to the
   real search population, the inverted/scrambled nulls and the injection hosts.
2. Interim detection threshold S*: generalized-Pareto tail fit to the max S of *vetted* inverted-null peaks
   (every check passes except the S threshold) on tuning (non-holdout) quality-passing stars, solved for a
   post-vetting false-alarm rate of 1 / N_stars per sector; star-level bootstrap for the CI.
3. Completeness re-scored consistently: detected = matched and S >= S*; passed = detected and every check passes.
The analysis PR (analysis/) supersedes this with its own fits; this module only keeps the site honest meanwhile.
  python -m exosearch.interim            # writes configs/interim_threshold_v1.json + products/.../analysis_layer/
"""
from __future__ import annotations

import glob
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from . import config, deepvet, nulltest, quality, store, vet
from .version import PIPELINE_VERSION, stamp

ROOT = Path(__file__).resolve().parents[1]
THR_FILE = ROOT / "configs" / "interim_threshold_v1.json"
PROD = config.DATA_ROOT / "products" / f"v{PIPELINE_VERSION}"
U_DEFAULT = 6.0           # = nulltest.S_FLOOR: every check is logged above it
U_SENS = (6.0, 7.0, 8.0)
N_BOOT = 2000
SEED = 20261110
REJECT_ALL = set(nulltest.SCREEN_CHECKS) | {"common_period"} | deepvet.REJECT


def read_parts(pattern):
    fs = sorted(glob.glob(str(pattern)))
    return pd.concat([pd.read_parquet(f) for f in fs], ignore_index=True) if fs else pd.DataFrame()


def load_null(tag="invert"):
    fin = PROD / "nulltests" / tag / f"null_{tag}_maxS.parquet"
    if fin.exists():
        d = pd.read_parquet(fin)
    else:
        d = read_parts(PROD / "nulltests" / tag / "parts" / "part-*.parquet")
        d = nulltest.finalize(d)          # adds fail_common_period among the null peaks themselves
    d = d[d.status == "ok"].copy()
    fl = d["flags"].fillna("").str.split(",").map(lambda x: set(x) - {"", "low_s"})
    dfl = d.deep_flags.fillna("").str.split(",").map(set)
    d["vetted"] = (fl.map(len) == 0) & ~d.fail_common_period.fillna(False).astype(bool) \
        & dfl.map(lambda x: not (x & deepvet.REJECT)) & d.checks_logged.fillna(False).astype(bool)
    return d


def with_quality(d, q):
    m = d.merge(q[["tic", "sector", "drop", "beta_max"]], on=["tic", "sector"], how="left")
    m["drop"] = m["drop"].fillna(True).astype(bool)      # no quality measurement -> excluded
    return m


def gpd_threshold(S_vetted, n_null, p, u):
    """S* with n_null * P(vetted, S > S*) = p * n_null; GPD on vetted exceedances over u."""
    y = np.sort(S_vetted[S_vetted > u] - u)
    k = len(y)
    if k < 5:
        return dict(u=u, k=k, s_star=None, xi=None, sigma=None, note="too few vetted exceedances")
    xi, _, sig = stats.genpareto.fit(y, floc=0)
    q = p * n_null / k                     # conditional exceedance probability above u
    s = u + float(stats.genpareto.isf(q, xi, 0, sig)) if q < 1 else float(np.quantile(S_vetted, 1 - p * n_null / max(len(S_vetted), 1)))
    return dict(u=u, k=int(k), xi=float(xi), sigma=float(sig), s_star=s)


def fit_sector(x, n_stars, u=U_DEFAULT, n_boot=N_BOOT, seed=SEED):
    """x: quality-passing tuning null rows for one sector. Returns point estimate, bootstrap CI, empirical check."""
    p = 1.0 / n_stars
    Sv = np.where(x.vetted.values, x.S.values, -np.inf)
    n = len(x)
    pt = gpd_threshold(Sv[np.isfinite(Sv)], n, p, u)
    rng = np.random.default_rng(seed)
    boots = []
    for _ in range(n_boot):
        b = Sv[rng.integers(0, n, n)]
        r = gpd_threshold(b[np.isfinite(b)], n, p, u)
        if r["s_star"] is not None and np.isfinite(r["s_star"]):
            boots.append(r["s_star"])
    boots = np.array(boots)
    emp_k = {s: int((x.vetted & (x.S >= s)).sum()) for s in (7, 8, 10, 12, 15, 20)}
    ci = (lambda a, lo, hi: [float(np.quantile(a, lo)), float(np.quantile(a, hi))] if len(a) > 20 else None)
    return dict(n_null=int(n), n_stars=int(n_stars), target_rate=p, u=u, k_tail=pt["k"], xi=pt["xi"], sigma=pt["sigma"],
                s_star=pt["s_star"], ci68=ci(boots, .16, .84), ci95=ci(boots, .025, .975), n_boot_ok=int(len(boots)),
                n_vetted=int(x.vetted.sum()), vetted_exceed=emp_k, max_vetted=float(x[x.vetted].S.max()) if x.vetted.any() else None,
                sensitivity={str(uu): gpd_threshold(Sv[np.isfinite(Sv)], n, p, uu)["s_star"] for uu in U_SENS})


def search_population(sectors, q):
    con = store.connect()
    out = {}
    for s in sectors:
        df = store.frame(con, s)
        df = df[df.status == "ok"]
        qs = q[q.sector == s]
        keep = set(qs[~qs["drop"]].tic.astype(int))
        out[int(s)] = dict(searched=int(len(df)), measured=int(len(qs)), dropped=int(len(df) - len(keep & set(df.tic.astype(int)))),
                           dropped_beta=int(((qs.status == "ok") & ~(qs.beta_max <= quality.rule()["threshold"])).sum()),
                           dropped_unmeasured=int(len(df) - (qs.status == "ok").sum()), kept=int(len(keep & set(df.tic.astype(int)))))
    return out


def rescore_injections(inj, q, s_star_by_sector):
    d = inj[inj.status == "ok"].copy()
    d = d.merge(q[["tic", "sector", "drop"]], on=["tic", "sector"], how="left")
    d["drop"] = d["drop"].fillna(True).astype(bool)
    d = d[~d["drop"]].copy()
    thr = d.sector.map(lambda s: s_star_by_sector.get(int(s), np.nan))
    failed = d.failed_checks.fillna("").str.split(",").map(lambda x: set(x) - {"", "low_s"})
    d["detected_v10"], d["passed_all_v10"] = d.detected, d.passed_all
    d["detected"] = d.matched.astype(bool) & (d.S >= thr)
    d["passed_all"] = d.detected & (failed.map(len) == 0)
    d["s_star"] = thr
    return d


def main():
    st = stamp()
    rule = quality.rule()
    q = quality.load()
    null = with_quality(load_null("invert"), q)
    sectors = sorted(int(s) for s in null.sector.unique())
    pop = search_population(sorted(set(sectors) | {70, 106}), q)
    from . import interim_review as rv                     # lazy: rv imports this module
    res = {}
    for s in sectors:
        x = null[(null.sector == s) & ~null["drop"] & ~null.holdout.astype(bool)]
        _, u0, stable = rv.gpd_diagnostics(x, pop[s]["kept"], s, plot=False)
        r = fit_sector(x, pop[s]["kept"], u=u0)
        sysv = [d["s_star"] for d in stable if d.get("s_star")]
        r.update(u_rule=f"lowest u >= {U_DEFAULT} with k >= {rv.MIN_K} whose xi is within 1 bootstrap SE of xi at every higher u with k >= 30",
                 u_stable=[d["u"] for d in stable],
                 sys_range=[float(min(sysv)), float(max(sysv))] if sysv else None,
                 sys_halfwidth=float((max(sysv) - min(sysv)) / 2) if sysv else None)
        r["dropped_null_quality"] = int(((null.sector == s) & null["drop"]).sum())
        r["held_out_excluded"] = int(((null.sector == s) & null.holdout.astype(bool) & ~null["drop"]).sum())
        res[str(s)] = r
        cams = {}
        for cam, xc in x.groupby("camera"):
            k_c = int((xc.vetted & (xc.S > u0)).sum())
            if k_c >= 30:
                c = fit_sector(xc, pop[s]["kept"] * len(xc) / len(x), u=u0, n_boot=500)      # own alpha = 1/N_cam
                Svc = xc.S.values[xc.vetted.values]
                c.update(alpha="own (1/N_cam)", s_star_sector_alpha=gpd_threshold(Svc, len(xc), 1.0 / pop[s]["kept"], u0)["s_star"])
                cams[str(int(cam))] = c
            else:
                cams[str(int(cam))] = dict(n_null=int(len(xc)), k_tail=k_c,
                                           s_star=None, note="fewer than 30 vetted exceedances; sector value applies")
        r["per_camera"] = cams
    out = dict(name="interim_threshold_v1", layer="interim analysis layer over frozen exosearch v1.0.0",
               pipeline_version=st["pipeline_version"], config_hash=st["config_hash"], git_hash=st["git_hash"],
               statistic="S = depth_box / (CDPP_D / sqrt(N_tr))", interim=True,
               definition=("S* solves N_null^-1 * #(vetted inverted-null max-S > S*) = 1/N_stars per sector, with a GPD tail fit "
                           "above u chosen per sector by the parameter-stability rule (u >= 6: all checks logged above it), quoted as value (95% stat CI) "
                           "+- systematic half-width over the stable u range, on tuning (non-holdout) stars passing {rule['name']}; "
                           f"CI from {N_BOOT} star-level bootstrap resamples. False-alarm rate is post-vetting. "
                           "Partial null: the sample grows as the inverted run proceeds."),
               quality_cut=rule["rule"], population=pop, sectors=res)
    THR_FILE.write_text(json.dumps(out, indent=1, default=float))
    (PROD / "analysis_layer").mkdir(parents=True, exist_ok=True)
    (PROD / "analysis_layer" / "interim_threshold_v1.json").write_text(THR_FILE.read_text())
    print(json.dumps({s: {k: r[k] for k in ("n_null", "n_stars", "u", "k_tail", "xi", "s_star", "ci95", "sys_range", "sys_halfwidth", "vetted_exceed")}
                      for s, r in res.items()}, indent=1, default=float))
    print(json.dumps(pop))


def s_star_by_sector():
    if not THR_FILE.exists():
        return {}
    j = json.loads(THR_FILE.read_text())
    return {int(s): r["s_star"] for s, r in j["sectors"].items() if r.get("s_star")}


if __name__ == "__main__":
    main()
