"""Post-hoc CDPP floor on stored null peaks (analysis layer; v1.0.0 and stored S unchanged).

CDPP_D' = max(CDPP_D, sigma_w / sqrt(N_D)), N_D = D / cadence. Since beta_D = CDPP_D sqrt(N_D) / sigma_w, this is
CDPP_D' = CDPP_D * max(1, 1 / beta_D) and S' = S * min(1, beta_D) at the stored peak. The peak itself is not re-chosen,
so S' is a lower bound on the floored per-star max (another trial could become the max); the floor only lowers S.
  python -m exosearch.beta_floor        # -> products/v<ver>/analysis_layer/beta_floor/beta_floor.json
"""
from __future__ import annotations

import json

import numpy as np

from . import interim, interim_review as rv, quality


def main(sector=70):
    q = quality.load()
    x = interim.with_quality(interim.load_null("invert"), q)
    x = x[(x.sector == sector) & ~x["drop"] & ~x.holdout.astype(bool)].copy()
    thr = json.loads(interim.THR_FILE.read_text())["sectors"][str(sector)]
    N, s_star = thr["n_stars"], thr["s_star"]
    low = x.beta_d < 1
    v = x.vetted.astype(bool)
    res = dict(sector=sector, n_null=int(len(x)), beta_lt_1=int(low.sum()), beta_lt_1_frac=float(low.mean()),
               beta_lt_0p8=int((x.beta_d < 0.8).sum()), beta_lt_0p5=int((x.beta_d < 0.5).sum()),
               vetted_tail_gt6=dict(n=int((v & (x.S > 6)).sum()), beta_lt_1=int((v & (x.S > 6) & low).sum())),
               vetted_tail_ge_s_star=dict(s_star=s_star, n=int((v & (x.S >= s_star)).sum()), beta_lt_1=int((v & (x.S >= s_star) & low).sum())),
               all_tail_gt6=dict(n=int((x.S > 6).sum()), beta_lt_1=int(((x.S > 6) & low).sum())),
               beta_median=float(x.beta_d.median()), beta_median_vetted_tail=float(x[v & (x.S > 6)].beta_d.median()))
    y = x.copy()
    y["S"] = x.S * np.minimum(1.0, x.beta_d.fillna(1.0))
    _, u0, stable = rv.gpd_diagnostics(y, N, sector, plot=False)
    f = interim.fit_sector(y, N, u=u0)
    sysv = [d["s_star"] for d in stable if d.get("s_star")]
    res["floored"] = dict(u=u0, k_tail=f["k_tail"], xi=f["xi"], s_star=f["s_star"], ci95=f["ci95"],
                          sys_range=[min(sysv), max(sysv)] if sysv else None, max_vetted=f["max_vetted"])
    res["committed_interim"] = dict(n_null=thr["n_null"], u=thr["u"], k_tail=thr["k_tail"], xi=thr["xi"], s_star=s_star, ci95=thr["ci95"],
                                    sys_range=thr.get("sys_range"))
    _, u1, stable1 = rv.gpd_diagnostics(x, N, sector, plot=False)          # same sample, no floor: like-for-like
    g = interim.fit_sector(x, N, u=u1)
    sys1 = [d["s_star"] for d in stable1 if d.get("s_star")]
    res["unfloored_same_sample"] = dict(u=u1, k_tail=g["k_tail"], xi=g["xi"], s_star=g["s_star"], ci95=g["ci95"],
                                        sys_range=[min(sys1), max(sys1)] if sys1 else None)
    res["shift_same_sample"] = f["s_star"] - g["s_star"]
    res["beta_hist"] = np.histogram(x.beta_d.dropna().clip(0, 2), bins=np.linspace(0, 2, 21))[0].tolist()
    out = interim.PROD / "analysis_layer" / "beta_floor"
    out.mkdir(parents=True, exist_ok=True)
    (out / "beta_floor.json").write_text(json.dumps(res, indent=1, default=float))
    print(json.dumps(res, indent=1, default=float))


if __name__ == "__main__":
    main()
