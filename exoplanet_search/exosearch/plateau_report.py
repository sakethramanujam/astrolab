"""Plateau statistics under the registered tie-break (amendment 1): disjoint peaks vs plateau-only, plateau width
against the number of transits, and the distance between the lowest-edge period and the plateau centre.
Widths are reported only; no grid change is proposed here (any change follows Ofir 2014 and is registered first).

  python -m exosearch.plateau_report --set S70_gate=gate.csv.gz --set S105_v10=plat105.csv.gz \
      --v10check v10check.csv.gz --out docs/host_gate/plateau_report.json
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

NTR_BINS = [(0, 2.5, "<=2"), (2.5, 4.5, "3-4"), (4.5, 8.5, "5-8"), (8.5, 1e9, ">=9")]


def q(x):
    x = pd.Series(x).dropna()
    return None if x.empty else dict(n=int(len(x)), min=float(x.min()), median=float(x.median()),
                                     p90=float(x.quantile(0.9)), max=float(x.max()))


def one(d):
    d = d[d.status == "ok"].copy()
    tie = d.n_tie_band.fillna(1) > 1
    plat = d.plateau_width > 1
    out = dict(n=int(len(d)), tie_band_gt1=int(tie.sum()),
               single_point=int(((d.n_plateaus == 1) & ~plat).sum()),
               plateau_only=int(((d.n_plateaus == 1) & plat).sum()),
               disjoint=int(d.disjoint.astype(bool).sum()),
               disjoint_with_chosen_plateau=int((d.disjoint.astype(bool) & plat).sum()))
    out["frac_disjoint"] = out["disjoint"] / out["n"]
    out["frac_plateau_only"] = out["plateau_only"] / out["n"]
    out["frac_of_tied_disjoint"] = out["disjoint"] / max(out["tie_band_gt1"], 1)
    out["frac_of_tied_plateau_only"] = out["plateau_only"] / max(out["tie_band_gt1"], 1)
    out["plateau_width_steps_chosen"] = q(d.plateau_width)
    out["plateau_width_steps_max"] = q(d.plateau_width_max)
    nt = d.n_tr.astype(float)
    out["width_vs_n_tr"] = [dict(n_tr=lab, n=int(m.sum()), n_width_gt1=int((m & plat).sum()),
                                 width_chosen=q(d.plateau_width[m]), width_max=q(d.plateau_width_max[m]))
                            for lo, hi, lab in NTR_BINS for m in [(nt > lo) & (nt <= hi)]]
    rho = d[["n_tr", "plateau_width_max"]].astype(float).corr(method="spearman").iloc[0, 1]
    out["spearman_n_tr_vs_width_max"] = None if rho != rho else float(rho)
    off = (d.period - d.P_lowest_edge)
    out["lowest_edge_to_centre"] = dict(days=q(off[plat]), sigma_P=q((off / d.sigma_P)[plat]),
                                        note="plateau stars only (width > 1); 0 for single points by construction")
    out["sigma_P_days"] = dict(all=q(d.sigma_P), plateau=q(d.sigma_P[plat]))
    return out


def flips(d, v10check, stored_period):
    """Do the v1.0.0 cached-vs-fresh (P, D) flips lie on one plateau (benign) or on disjoint peaks?"""
    c = pd.read_csv(v10check).set_index("tic")
    d = d.set_index("tic")
    res = []
    for tic in c.index[(c.period != stored_period.reindex(c.index))]:
        if tic not in d.index:
            continue
        r = d.loc[tic]
        pl = json.loads(r.plateaus)
        h = lambda P: 0.5 * P * P * float(r.grid_df)

        def where(P):
            for i, p in enumerate(pl):
                if p[0] - h(P) * 1.0001 <= P <= p[1] + h(P) * 1.0001:
                    return i
            return None
        a, b = where(float(stored_period[tic])), where(float(c.period[tic]))
        res.append(dict(tic=int(tic), P_cached=float(stored_period[tic]), P_fresh=float(c.period[tic]),
                        plateau_cached=a, plateau_fresh=b, n_plateaus=int(r.n_plateaus),
                        same_plateau=a is not None and a == b))
    return dict(n=len(res), same_plateau=sum(x["same_plateau"] for x in res),
                disjoint_peaks=sum(x["plateau_cached"] is not None and x["plateau_fresh"] is not None
                                   and not x["same_plateau"] for x in res),
                not_located=sum(x["plateau_cached"] is None or x["plateau_fresh"] is None for x in res), rows=res)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--set", action="append", required=True, help="name=csv")
    ap.add_argument("--v10check", default=None, help="fresh-JIT v1.0 rows of the S105 set")
    ap.add_argument("--flip-set", default="S105_v10")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    out = {}
    for s in a.set:
        name, f = s.split("=", 1)
        out[name] = one(pd.read_csv(f))
        if a.v10check and name == a.flip_set:
            from . import store
            st = store.frame(store.connect(), 105).set_index("tic").period
            out[name]["v10_cache_flips"] = flips(pd.read_csv(f), a.v10check, st)
    Path(a.out).write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps({k: {kk: vv for kk, vv in v.items() if kk not in ("width_vs_n_tr",)} for k, v in out.items()}, indent=1)[:6000])
