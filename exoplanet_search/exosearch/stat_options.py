"""Interim-layer variants of the detection statistic, applied post hoc to stored peaks (v1.0.0 unchanged).

Depth audit. v1.0.0 computes S = depth_box * sqrt(n_tr) / cdpp_d, where depth_box is the depth astropy's BLS
'fast' method reports at the S-scan peak. At short periods / long durations that value can exceed the depth of any
box at the reported (P, D, T0) (checked against compute_stats, an unweighted box mean and a 'slow' BLS refit), so
S is inflated for some non-transit peaks. `depth` (compute_stats at the same P, D, T0: in-transit vs out-of-transit
mean on the binned flux) is the box depth, and is what Rp, too_deep, the site and the injection depth_rec already use.
  S_exact = depth * sqrt(n_tr) / cdpp_d        (injections: depth_rec, n_tr_rec, cdpp_d_rec)

Red-noise treatments on top of S_exact, beta' replacing beta (S' = S_exact * beta / beta'):
  a  none                beta' = beta
  b  hard floor          beta' = max(1, beta)
  c  significance floor  beta' = max(beta, 1 - 2 sigma_beta), sigma_beta = 1 / sqrt(2 N_w), N_w = n_cad * cadence / D
  d  empirical Bayes     beta' = posterior mean of beta under a normal prior per duration D (prior from the
                         star-quality beta_D of every kept star at that D; measurement variance sigma_beta^2)
The peak is not re-chosen, so every S' is a lower bound on that variant's per-star max and each S* is a lower bound.
Identical rule on the flipped null and on the injections.
  python -m exosearch.stat_options [--no-cut]  -> products/v<ver>/analysis_layer/stat_options/stat_options[_nocut].json
  (--no-cut: no star-quality cut on the null or the injections; alpha = 1 / all measured searched stars)
"""
from __future__ import annotations

import json

import numpy as np
import pandas as pd

from . import bls, interim, interim_review as rv, quality

CAD = 2.0 / 1440
EDGES = (0, 10, 15, 25, 50, 100, 1e9)


def prior(q):
    """Per duration: (median beta, between-star variance tau^2) from kept stars' own-light-curve beta_D."""
    out = {}
    k = q[(q.status == "ok") & ~q["drop"]]
    for D in bls.DURATIONS:
        b = k[f"beta_{D:.2f}"].dropna()
        nw = (k.loc[b.index, "n_cad"] * CAD / D).clip(lower=1)
        s2 = 1 / (2 * nw)
        out[f"{D:.2f}"] = (float(b.median()), float(max(b.var() - s2.mean(), 1e-4)))
    return out


def variants(S, beta, D, n_cad, pri):
    nw = np.maximum(np.asarray(n_cad, float) * CAD / np.asarray(D, float), 1.0)
    sb = 1 / np.sqrt(2 * nw)
    m = np.array([pri[f"{d:.2f}"][0] if f"{d:.2f}" in pri else 1.0 for d in D])
    t2 = np.array([pri[f"{d:.2f}"][1] if f"{d:.2f}" in pri else 1.0 for d in D])
    bp = {"a": beta,
          "b": np.maximum(1.0, beta),
          "c": np.maximum(beta, 1 - 2 * sb),
          "d": (beta / sb ** 2 + m / t2) / (1 / sb ** 2 + 1 / t2)}
    return {k: S * beta / v for k, v in bp.items()}, sb


def main(sector=70, cut=True):
    q = quality.load()
    pri = prior(q[q.sector == sector])
    null = interim.with_quality(interim.load_null("invert"), q)
    x = null[(null.sector == sector) & (~null["drop"] if cut else True) & ~null.holdout.astype(bool)].copy()
    x = x.merge(q[["tic", "sector", "n_cad"]], on=["tic", "sector"], how="left")
    x["S_exact"] = x.depth * np.sqrt(x.n_tr) / x.cdpp_d
    thr = json.loads(interim.THR_FILE.read_text())["sectors"][str(sector)]
    N = thr["n_stars"] if cut else int((q[q.sector == sector].status == "ok").sum())   # no cut: all measured searched stars
    inj = interim.read_parts(interim.PROD / "injections" / "*" / "part-*.parquet")
    inj = inj[(inj.status == "ok") & (inj.sector == sector) & ~inj.holdout.astype(bool)]
    inj = inj.merge(q[["tic", "sector", "beta_max", "status", "n_cad"]].rename(columns={"status": "q_status"}), on=["tic", "sector"], how="left")
    inj = inj[inj.q_status.eq("ok") & ((inj.beta_max <= quality.rule()["threshold"]) if cut else True)].copy()
    inj["S_exact"] = inj.depth_rec * np.sqrt(inj.n_tr_rec) / inj.cdpp_d_rec
    inj["allowed"] = (inj.depth_inj < rv.MAX_DEPTH) & (inj.rp_inj < rv.MAX_RP)
    failed = inj.failed_checks.fillna("").str.split(",").map(lambda f: set(f) - {"", "low_s"})
    clean = failed.map(len) == 0
    vn, sbn = variants(x.S_exact.values, x.beta_d.values, x.duration.values, x.n_cad.fillna(x.n_cad.median()).values, pri)
    vi, _ = variants(inj.S_exact.values, inj.beta_d.values, inj.dur_rec.values, inj.n_cad.fillna(x.n_cad.median()).values, pri)
    opts = {"v1.0 S (stored)": (x.S_v10.values, inj.S.values), "a: exact depth, no floor": (vn["a"], vi["a"]),
            "b: hard floor": (vn["b"], vi["b"]), "c: 1-2 sigma_beta floor": (vn["c"], vi["c"]), "d: empirical Bayes": (vn["d"], vi["d"])}
    res = dict(sector=sector, quality_cut=bool(cut), n_stars=int(N), n_null=int(len(x)), n_inj=int(len(inj)), prior_by_D=pri,
               sigma_beta_median=float(np.median(sbn)), note="S* values are lower bounds (stored peaks, not re-searched)",
               options={})
    for name, (sn, si) in opts.items():
        y = x.copy(); y["S"] = sn
        _, u0, stable = rv.gpd_diagnostics(y, N, sector, plot=False)
        f = interim.fit_sector(y, N, u=u0, n_boot=500)
        sysv = [d["s_star"] for d in stable if d.get("s_star")]
        det = inj.matched.astype(bool).values & (si >= f["s_star"])
        vet = det & clean.values
        al = inj.allowed.values
        b = pd.cut(inj.snr_exp, list(EDGES), right=False)
        comp = {}
        for k_, idx in inj.groupby(b, observed=True).groups.items():
            ii = inj.index.get_indexer(idx); ii = ii[al[ii]]
            comp[str(k_)] = dict(n=int(len(ii)), detected=round(float(det[ii].mean()), 3) if len(ii) else None,
                                 vetted=round(float(vet[ii].mean()), 3) if len(ii) else None)
        res["options"][name] = dict(u=u0, k_tail=f["k_tail"], xi=f["xi"], s_star=f["s_star"], ci95=f["ci95"],
                                    sys_range=[min(sysv), max(sysv)] if sysv else None, max_vetted=f["max_vetted"],
                                    null_changed_gt1pct_vs_v10=int((np.abs(sn / x.S_v10.values - 1) > 0.01).sum()),
                                    null_lowered_gt1pct_vs_exact=int((sn < vn["a"] * 0.99).sum()),
                                    null_raised_gt1pct_vs_exact=int((sn > vn["a"] * 1.01).sum()),
                                    inj_changed_gt1pct_vs_exact=int((np.abs(si / vi["a"] - 1) > 0.01).sum()),
                                    completeness_allowed=comp)
        print(name, json.dumps({k: v for k, v in res["options"][name].items() if k != "completeness_allowed"}, default=float), flush=True)
    out = interim.PROD / "analysis_layer" / "stat_options"
    out.mkdir(parents=True, exist_ok=True)
    (out / ("stat_options.json" if cut else "stat_options_nocut.json")).write_text(json.dumps(res, indent=1, default=float))


if __name__ == "__main__":
    import sys
    main(cut="--no-cut" not in sys.argv)
