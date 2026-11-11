"""Headline numbers with and without the star-quality cut (interim layer, statistic S_exact).

Per sector with its own flipped null: S* (u by the stability rule, star-level bootstrap), expected post-vetting
false alarms per sector at both thresholds, vetted completeness at expected S 15-25, and the S105 counts
(v1.0.0 screens passed and S_exact >= S*; S105 has no own null, so these are 'not calibrated' counts against the
S70 threshold).
  python -m exosearch.cut_sensitivity   -> products/v<ver>/analysis_layer/cut_sensitivity/cut_sensitivity.json
"""
from __future__ import annotations

import json

import numpy as np
import pandas as pd

from . import interim, interim_review as rv, quality, store


def _fit(x, n):
    _, u, _ = rv.gpd_diagnostics(x, n, 0, plot=False)
    f = interim.fit_sector(x, n, u=u, n_boot=500)
    return dict(s_star=f["s_star"], ci95=f["ci95"], u=u, k_tail=f["k_tail"], n_null=int(len(x)), n_stars=int(n))


def _efa(x, n_stars, u, s):
    Sv = x.S.values[x.vetted.values]
    return float(rv.tail_rate(Sv, len(x), u, s) * n_stars)


def main():
    q = quality.load((70, 105, 106))
    thr = json.loads(interim.THR_FILE.read_text())["sectors"]
    null = interim.with_quality(interim.load_null("invert"), q)
    inj = interim.read_parts(interim.PROD / "injections" / "*" / "part-*.parquet")
    res = dict(statistic="S_exact", note="cut = star_quality_v1 (beta_max <= 3); no cut = every measured searched star; "
                                         "alpha = 1/N_stars of that population", sectors={})
    for sec, r in thr.items():
        s = int(sec)
        qs = q[q.sector == s]
        tun = null[(null.sector == s) & ~null.holdout.astype(bool)]
        meas = tun.merge(qs[["tic", "status"]].rename(columns={"status": "q_status"}), on="tic", how="left")
        pops = {"with_cut": (tun[~tun["drop"]], int(r["n_stars"]), True),
                "without_cut": (meas[meas.q_status.eq("ok")], int((qs.status == "ok").sum()), False)}
        out = {}
        for k, (x, n, cut) in pops.items():
            f = _fit(x, n) if k == "without_cut" else dict(s_star=r["s_star"], ci95=r["ci95"], u=r["u"], k_tail=r["k_tail"],
                                                           n_null=r["n_null"], n_stars=n)
            c = rv.completeness(inj, q, quality.rule()["threshold"] if cut else None, f["s_star"], sector=s)
            c = c[c.allowed]
            m = (c.snr_exp >= 15) & (c.snr_exp < 25)
            f["completeness_15_25"] = dict(n=int(m.sum()), detected=round(float(c[m].det.mean()), 3), vetted=round(float(c[m].vet.mean()), 3))
            out[k] = (f, x)
        a, b = out["with_cut"][0], out["without_cut"][0]
        a["efa_per_sector_at_own_s_star"] = round(_efa(out["with_cut"][1], a["n_stars"], a["u"], a["s_star"]), 3)
        b["efa_per_sector_at_own_s_star"] = round(_efa(out["without_cut"][1], b["n_stars"], b["u"], b["s_star"]), 3)
        b["efa_per_sector_at_with_cut_s_star"] = round(_efa(out["without_cut"][1], b["n_stars"], b["u"], a["s_star"]), 3)
        res["sectors"][sec] = dict(with_cut=a, without_cut=b)
    # S105: real search, not calibrated (no own null); counts against the S70 thresholds
    d = store.frame(store.connect(), 105)
    d = d[(d.status == "ok") & d.S.notna()].copy()
    d["S_exact"] = d.depth * np.sqrt(d.n_tr) / d.cdpp_d
    q105 = q[q.sector == 105]
    ref = res["sectors"].get("70")
    if ref:
        scr = d.passed == 1
        row = dict(rows=int(len(d)), passed_screens_v10=int(scr.sum()),
                   ge_s_star_without_cut=int((scr & (d.S_exact >= ref["without_cut"]["s_star"])).sum()),
                   ge_s_star_with_cut_threshold_all_stars=int((scr & (d.S_exact >= ref["with_cut"]["s_star"])).sum()),
                   s_star_ref=dict(with_cut=ref["with_cut"]["s_star"], without_cut=ref["without_cut"]["s_star"]))
        if len(q105):
            keep = set(q105[~q105["drop"]].tic.astype(int))
            k = d.tic.astype(int).isin(keep)
            row["ge_s_star_with_cut"] = int((scr & k & (d.S_exact >= ref["with_cut"]["s_star"])).sum())
            row["quality_measured"] = int(len(q105))
        else:
            row["ge_s_star_with_cut"] = None
            row["quality_note"] = "star quality not yet measured for S105; with-cut count pending"
        res["s105_not_calibrated"] = row
    out = interim.PROD / "analysis_layer" / "cut_sensitivity"
    out.mkdir(parents=True, exist_ok=True)
    (out / "cut_sensitivity.json").write_text(json.dumps(res, indent=1, default=float))
    print(json.dumps(res, indent=1, default=float))


if __name__ == "__main__":
    main()
