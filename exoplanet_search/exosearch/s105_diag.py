"""Sector-105 diagnostics for signals that crossed the (borrowed, not applicable) sector-70 S* (analysis layer).

  python -m exosearch.s105_diag --candidates <progress.json with the listed rows> [--sector 105]

For each listed signal: fraction of its in-data epochs within 0.1 d of a momentum dump (QUALITY bit 32), of a
gap / orbit edge, or inside a scattered-light window (QUALITY bits 2048 | 4096), against a control of the same
(P, D) with uniformly random epochs; P near 13.7/n d (TESS orbit aliases, 1%) or near the dump spacing / n;
clusters of signals sharing P (0.5%) and phase (0.1 d). Population comparison: median beta_D at the peak and
median max S of the sector's real v1 search against the sector-70 flipped null (same estimator, flipped flux).
Output: products/v<ver>/analysis_layer/s<sector>/{s105_diag.json, s105_survivors.csv, s105_timing.png}
"""
from __future__ import annotations

import argparse
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from . import batch, interim, lightcurve, quality, store, targets

DESAT, SCAT = 32, 2048 | 4096
WIN, GAP_D = 0.1, 0.5


def events(sector, df):
    man = targets.fetch_manifest(sector).set_index("tic")
    out = {}
    for cam, g in df.groupby("camera"):
        for tic in g.sort_values("tmag").tic:                                    # a bright star per camera
            if tic not in man.index:
                continue
            p = batch.lc_path(sector, man.loc[tic].fname, download=True)
            if p is None:
                continue
            t, f, e, q, meta = lightcurve.load_spoc_raw(p)
            qi = q.astype(int)
            good = t[lightcurve.good_mask(f, e, q)]
            br = np.where(np.diff(good) > GAP_D)[0]
            out[int(cam)] = dict(dumps=np.unique(np.round(t[(qi & DESAT) > 0], 4)), scat=t[(qi & SCAT) > 0],
                                 edges=np.sort(np.r_[good[0], good[br], good[br + 1], good[-1]]), good=good,
                                 tmin=good[0], tmax=good[-1], ref_tic=int(tic))
            break
    return out


def epoch_metrics(P, t0, ev, D):
    k = np.arange(np.ceil((ev["tmin"] - t0) / P), np.floor((ev["tmax"] - t0) / P) + 1)
    tc = t0 + P * k
    g = ev["good"]
    tc = tc[np.abs(g[np.clip(np.searchsorted(g, tc), 0, g.size - 1)] - tc) < max(D, 0.05)]
    if tc.size == 0:
        return None
    near = lambda arr: (np.min(np.abs(tc[:, None] - arr[None, :]), axis=1) < WIN) if arr.size else np.zeros(tc.size, bool)
    return dict(n=int(tc.size), dump=float(near(ev["dumps"]).mean()), edge=float(near(ev["edges"]).mean()),
                scat=float(near(ev["scat"]).mean()))


def alias_expected(P, rng, n_max=30, trials=300):
    """Expected number of 13.7/n hits (1%) if each period were jittered uniformly in log P by +-10% (control)."""
    hit = lambda x: any(abs(x / (13.7 / n) - 1) < 0.01 for n in range(1, n_max))
    return float(np.mean([sum(hit(p * np.exp(rng.uniform(-.1, .1))) for p in P) for _ in range(trials)]))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--candidates", required=True)
    ap.add_argument("--sector", type=int, default=105)
    a = ap.parse_args()
    out = interim.PROD / "analysis_layer" / f"s{a.sector:04d}"
    out.mkdir(parents=True, exist_ok=True)
    c = pd.DataFrame([x for x in json.loads(open(a.candidates).read())["candidates"]
                      if x["status"] == "needs_review" and x["sector"] == a.sector])
    df = store.frame(store.connect(), a.sector)
    df = df[df.status == "ok"]
    c = c.merge(df[["tic", "camera", "duration", "beta_d", "n_tr"]], on="tic", how="left")
    ev = events(a.sector, df)
    rng = np.random.default_rng(20261110)
    rows, ctl = [], []
    for r in c.itertuples():
        e = ev.get(int(r.camera))
        m = epoch_metrics(r.period, r.t0_btjd, e, r.duration) if e else None
        cm = [epoch_metrics(r.period, r.t0_btjd + rng.uniform(0, r.period), e, r.duration) for _ in range(50)] if e else []
        cm = [x for x in cm if x]
        alias = [n for n in range(1, 30) if abs(r.period / (13.7 / n) - 1) < 0.01]
        rows.append(dict(tic=int(r.tic), camera=int(r.camera), period=r.period, t0=r.t0_btjd, depth_ppm=r.depth_ppm, S=r.S,
                         duration_h=r.duration_h, n_tr=r.n_tr, beta_d=r.beta_d, tmag=r.tmag, flags=r.flags,
                         orbit_alias_n=alias[0] if alias else None,
                         **({f"ep_{k}": v for k, v in m.items()} if m else {}),
                         **({f"ctl_{k}": float(np.mean([x[k] for x in cm])) for k in ("dump", "edge", "scat")} if cm else {})))
    s = pd.DataFrame(rows)
    # clusters in (P, phase)
    s = s.sort_values("period").reset_index(drop=True)
    grp = np.zeros(len(s), int); gid = 0
    for i in range(len(s)):
        if grp[i]:
            continue
        gid += 1; grp[i] = gid
        for j in range(i + 1, len(s)):
            if not grp[j] and abs(s.period[j] / s.period[i] - 1) < 0.005:
                dph = ((s.t0[j] - s.t0[i] + 0.5 * s.period[i]) % s.period[i]) - 0.5 * s.period[i]
                if abs(dph) < 0.1:
                    grp[j] = gid
    s["cluster"] = grp
    cl = s.groupby("cluster").filter(lambda g: len(g) > 1)
    w = lambda col: float(np.average(s[col].dropna(), weights=s.loc[s[col].notna(), "ep_n"])) if col in s and s[col].notna().any() else None
    dumps = {c_: (np.diff(v["dumps"][np.r_[True, np.diff(v["dumps"]) > 0.05]]).tolist()) for c_, v in ev.items()}
    sp = np.median(np.concatenate([np.array(x) for x in dumps.values() if len(x)])) if any(len(x) for x in dumps.values()) else None
    dump_alias = [int(r.tic) for r in s.itertuples() if sp and any(abs(r.period / (sp / n) - 1) < 0.01 for n in range(1, 10))]
    # population comparison
    null = interim.load_null("invert")
    q = quality.load()
    pop = {f"S{a.sector}_real_v1": dict(n=int(len(df)), median_S=float(df.S.median()), p99_S=float(df.S.quantile(.99)),
                                         median_beta_d=float(df.beta_d.median()), frac_S_ge_12p4=float((df.S >= 12.39).mean()))}
    for sec in (70, 106):
        x = null[null.sector == sec]
        qq = q[(q.sector == sec) & (q.status == "ok")]
        pop[f"S{sec}"] = dict(flipped_null_n=int(len(x)), flipped_median_S=float(x.S.median()) if len(x) else None,
                              flipped_p99_S=float(x.S.quantile(.99)) if len(x) else None,
                              flipped_median_beta_d=float(x.beta_d.median()) if len(x) else None,
                              flipped_frac_S_ge_12p4=float((x.S >= 12.39).mean()) if len(x) else None,
                              quality_median_beta_max=float(qq.beta_max.median()), real_v1_S="not available (v0 rows)")
    res = dict(sector=a.sector, n=int(len(s)), reference_lcs={k: v["ref_tic"] for k, v in ev.items()},
               n_dumps={k: int(v["dumps"].size) for k, v in ev.items()}, dump_spacing_d=sp,
               epoch_frac_near_dump=w("ep_dump"), control_near_dump=w("ctl_dump"),
               epoch_frac_near_edge=w("ep_edge"), control_near_edge=w("ctl_edge"),
               epoch_frac_scattered=w("ep_scat"), control_scattered=w("ctl_scat"),
               n_orbit_alias=int(s.orbit_alias_n.notna().sum()),
               n_orbit_alias_expected=alias_expected(s.period.values, rng), orbit_alias_tics=s[s.orbit_alias_n.notna()].tic.tolist(),
               n_dump_alias=len(dump_alias), dump_alias_tics=dump_alias,
               n_signals_with_majority_epochs_near_dump=int((s.ep_dump > 0.5).sum()),
               n_signals_le2_epochs=int((s.ep_n <= 2).sum()),
               clusters=[dict(tics=g.tic.tolist(), period=float(g.period.median())) for _, g in cl.groupby("cluster")],
               period_quartiles=s.period.quantile([.25, .5, .75]).tolist(), population=pop)
    s.to_csv(out / "s105_survivors.csv", index=False)
    (out / "s105_diag.json").write_text(json.dumps(res, indent=1, default=float))
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.4), dpi=100)
    ax[0].hist(s.period, bins=np.geomspace(0.5, 15, 30)); ax[0].set_xscale("log")
    for n in range(1, 15):
        ax[0].axvline(13.7 / n, color="#e45756", lw=0.5, ls=":")
    ax[0].set_xlabel("period (d); dotted: 13.7/n"); ax[0].set_ylabel("signals")
    ax[1].scatter(s.period, s.S, s=10); ax[1].set_xscale("log"); ax[1].set_xlabel("period (d)"); ax[1].set_ylabel("S")
    fig.suptitle(f"S{a.sector}: {len(s)} signals above the borrowed S70 S* (not applicable to S{a.sector})", fontsize=9)
    fig.tight_layout(); fig.savefig(out / "s105_timing.png"); plt.close(fig)
    print(json.dumps(res, indent=1, default=float))
    print(s[["tic", "period", "t0", "depth_ppm", "S", "ep_n", "ep_dump", "ctl_dump", "ep_edge", "ep_scat", "orbit_alias_n", "beta_d", "flags"]].round(4).to_string())


if __name__ == "__main__":
    main()
