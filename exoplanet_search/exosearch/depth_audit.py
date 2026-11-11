"""Depth-column audit (interim layer; v1.0.0 unchanged).

Compares depth_box (BLS 'fast' depth at the S peak; enters S) with depth (compute_stats box depth at the same
P, D, T0; enters Rp, too_deep, no_dip, odd/even, secondary, rp_too_large, site depth_ppm, injection depth_rec,
rp_rjup and the density fit) on every population, and counts how many rows cross S thresholds when S is recomputed
as S_exact = depth * sqrt(n_tr) / cdpp_d.
  python -m exosearch.depth_audit   -> products/v<ver>/analysis_layer/depth_audit/depth_audit.json
"""
from __future__ import annotations

import json

import numpy as np
import pandas as pd

from . import interim, store

DP_BINS = (0, 0.02, 0.05, 0.1, 0.2, 0.6)


def ratio_stats(r):
    r = pd.Series(r).replace([np.inf, -np.inf], np.nan)
    bad = (r > 2) | (r < 0.5) | r.isna()
    q = r.dropna()
    return dict(n=int(len(r)), median=round(float(q.median()), 3) if len(q) else None,
                p10=round(float(q.quantile(.1)), 3) if len(q) else None, p90=round(float(q.quantile(.9)), 3) if len(q) else None,
                p99=round(float(q.quantile(.99)), 3) if len(q) else None, n_gt2x_or_lt_half_or_nonpos=int(bad.sum()),
                frac_bad=round(float(bad.mean()), 4) if len(r) else None)


def block(d, depth_a, depth_b, S, s_exact, P, D, thr=(6, 7.3, 10, 12.39)):
    r = np.where(depth_b > 0, depth_a / depth_b, np.nan)
    dp = D / P
    out = dict(all=ratio_stats(r), by_D_over_P={}, by_S={},
               S_exact_over_S_median=round(float(np.nanmedian(s_exact / S)), 3),
               crossings={})
    for lo, hi in zip(DP_BINS[:-1], DP_BINS[1:]):
        m = (dp >= lo) & (dp < hi)
        out["by_D_over_P"][f"{lo}-{hi}"] = ratio_stats(r[m])
    for lo, hi in ((0, 6), (6, 10), (10, 12.39), (12.39, 1e9)):
        m = (S >= lo) & (S < hi)
        out["by_S"][f"{lo}-{hi}"] = ratio_stats(r[m])
    for t in thr:
        out["crossings"][str(t)] = dict(S_ge=int((S >= t).sum()), S_exact_ge=int((s_exact >= t).sum()),
                                        lost=int(((S >= t) & (s_exact < t)).sum()), gained=int(((S < t) & (s_exact >= t)).sum()))
    return out


def main():
    res = dict(note="ratio = depth_box / depth (injections: implied depth_box = S*cdpp_d_rec/sqrt(n_tr_rec) over depth_rec); "
                    "'bad' = ratio > 2, < 0.5, or depth <= 0", populations={})
    con = store.connect()
    for s in (70, 105, 106):
        d = store.frame(con, s)
        d = d[(d.status == "ok") & d.depth_box.notna()] if "status" in d else d[d.depth_box.notna()]
        se = d.depth * np.sqrt(d.n_tr) / d.cdpp_d
        res["populations"][f"real_s{s:04d}"] = block(d, d.depth_box.values, d.depth.values, d.S.values, se.values, d.period.values, d.duration.values)
    n = interim.load_null("invert")
    n = n[(n.status == "ok") & n.depth_box.notna()]
    for s, g in n.groupby("sector"):
        se = g.depth * np.sqrt(g.n_tr) / g.cdpp_d
        res["populations"][f"null_invert_s{int(s):04d}"] = block(g, g.depth_box.values, g.depth.values, g.S.values, se.values, g.period.values, g.duration.values)
    i = interim.read_parts(interim.PROD / "injections" / "*" / "part-*.parquet")
    i = i[(i.status == "ok") & i.matched.astype(bool)]
    imp = i.S * i.cdpp_d_rec / np.sqrt(i.n_tr_rec)
    se = i.depth_rec * np.sqrt(i.n_tr_rec) / i.cdpp_d_rec
    res["populations"]["injections_matched"] = block(i, imp.values, i.depth_rec.values, i.S.values, se.values, i.P_rec.values, i.dur_rec.values)
    res["populations"]["injections_matched"]["depth_rec_over_depth_inj"] = ratio_stats((i.depth_rec / i.depth_inj).values)
    out = interim.PROD / "analysis_layer" / "depth_audit"
    out.mkdir(parents=True, exist_ok=True)
    (out / "depth_audit.json").write_text(json.dumps(res, indent=1))
    for k, v in res["populations"].items():
        print(k, json.dumps(v["all"]), "S_exact/S", v["S_exact_over_S_median"], json.dumps(v["crossings"]))


if __name__ == "__main__":
    main()
