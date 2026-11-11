"""Diagnostics of the inverted-null tail for the methods paper (analysis layer over frozen v1.0.0).

* per-star max null S against the star's red-noise factor beta_max and against Tmag (Spearman rho),
* where null peak epochs fall relative to momentum dumps (QUALITY bit 32, reaction-wheel desaturation)
  and to data-gap / orbit edges, against a control with the same (P, D) and a uniformly random epoch.
  python -m exosearch.nulldiag
Writes products/v<ver>/analysis_layer/diag/{null_S_vs_beta.png, null_S_vs_tmag.png, null_peak_timing.png, diag.json}.
"""
from __future__ import annotations

import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

from . import batch, interim, lightcurve, quality, targets

DESAT_BIT = 32
GAP_D = 0.5
U = interim.U_DEFAULT


def sector_events(sector, tics):
    """Momentum-dump times and gap/orbit edges from reference light curves (one per camera)."""
    man = targets.fetch_manifest(sector).set_index("tic")
    out = {}
    for cam, tic in tics.items():
        p = batch.lc_path(sector, man.loc[tic].fname, download=False)
        t, f, e, q, meta = lightcurve.load_spoc_raw(p)
        dumps = t[(q.astype(int) & DESAT_BIT) > 0]
        good = t[lightcurve.good_mask(f, e, q)]
        br = np.where(np.diff(good) > GAP_D)[0]
        edges = np.r_[good[0], good[br], good[br + 1], good[-1]]
        out[int(cam)] = dict(dumps=np.unique(np.round(dumps, 4)), edges=np.sort(edges), tmin=good[0], tmax=good[-1], good=good)
    return out


def epoch_distances(P, t0, ev):
    tc = t0 + P * np.arange(np.ceil((ev["tmin"] - t0) / P), np.floor((ev["tmax"] - t0) / P) + 1)
    g = ev["good"]
    tc = tc[np.abs(g[np.clip(np.searchsorted(g, tc), 0, g.size - 1)] - tc) < 0.05]     # epochs with data
    if tc.size == 0 or ev["dumps"].size == 0:
        return np.array([]), np.array([])
    dd = np.min(np.abs(tc[:, None] - ev["dumps"][None, :]), axis=1)
    de = np.min(np.abs(tc[:, None] - ev["edges"][None, :]), axis=1)
    return dd, de


def main():
    out_dir = interim.PROD / "analysis_layer" / "diag"
    out_dir.mkdir(parents=True, exist_ok=True)
    q = quality.load()
    d = interim.with_quality(interim.load_null("invert"), q)
    d = d.merge(q[["tic", "sector", "tmag", "sigma_w"]].rename(columns={"tmag": "tmag_q", "sigma_w": "sigma_w_q"}),
                on=["tic", "sector"], how="left")
    res = {"n_null": int(len(d)), "n_quality_pass": int((~d["drop"]).sum())}
    sp = {}
    for name, sel in (("all", d.index == d.index), ("quality_pass", ~d["drop"]), ("quality_pass_S_gt_u", ~d["drop"] & (d.S > U)),
                      ("vetted_quality_pass", ~d["drop"] & d.vetted)):
        x = d[sel]
        sp[name] = dict(n=int(len(x)),
                        rho_S_beta_max=float(stats.spearmanr(x.S, x.beta_max, nan_policy="omit")[0]),
                        rho_S_beta_peak=float(stats.spearmanr(x.S, x.beta_d, nan_policy="omit")[0]),
                        rho_S_tmag=float(stats.spearmanr(x.S, x.tmag, nan_policy="omit")[0]),
                        rho_S_sigma_w=float(stats.spearmanr(x.S, x.sigma_w_q, nan_policy="omit")[0]))
    res["spearman"] = sp
    tail = d[~d["drop"] & (d.S > U)]
    res["tail_vs_bulk"] = dict(
        beta_max_median_tail=float(tail.beta_max.median()), beta_max_median_bulk=float(d[~d["drop"] & (d.S <= U)].beta_max.median()),
        frac_single_epoch_tail=float((tail.n_tr <= 1.5).mean()), frac_duration_max_tail=float((tail.duration >= 0.25).mean()),
        frac_duration_min_tail=float((tail.duration <= 0.04).mean()),
        flags_tail={k: int(v) for k, v in tail["flags"].fillna("").str.split(",").explode().value_counts().items() if k})
    # plots: S vs beta, S vs Tmag
    for col, fname, xl, logx in (("beta_max", "null_S_vs_beta.png", "star red-noise factor  beta_max", True),
                                 ("tmag", "null_S_vs_tmag.png", "TESS magnitude", False)):
        fig, ax = plt.subplots(figsize=(6.4, 4.2), dpi=110)
        for lab, sel, c, a in (("dropped by quality cut", d["drop"], "#bbbbbb", .5), ("kept", ~d["drop"] & ~d.vetted, "#4c78a8", .45),
                               ("kept, passes all vetting", ~d["drop"] & d.vetted, "#e45756", .9)):
            x = d[sel]
            ax.scatter(x[col], x.S.clip(lower=3), s=5, c=c, alpha=a, label=f"{lab} ({len(x)})", lw=0)
        ax.set_yscale("log"); ax.set_xscale("log" if logx else "linear")
        if col == "beta_max":
            ax.axvline(quality.rule()["threshold"], color="k", ls=":", lw=1)
        ax.set_xlabel(xl); ax.set_ylabel("per-star max S (inverted light curve)")
        r = sp["quality_pass"][f"rho_S_{'beta_max' if col == 'beta_max' else 'tmag'}"]
        ax.set_title(f"Inverted null, sector 70 (partial, N={len(d)}); Spearman rho (kept) = {r:.2f}", fontsize=9)
        ax.legend(fontsize=7, loc="upper left", markerscale=2)
        fig.tight_layout(); fig.savefig(out_dir / fname); plt.close(fig)
    # timing vs momentum dumps and gap edges
    timing = {}
    rng = np.random.default_rng(20261110)
    fig, axs = plt.subplots(1, 2, figsize=(9, 3.6), dpi=110)
    for s, ds in d.groupby("sector"):
        refs = ds.groupby("camera").tic.first().to_dict()
        ev = sector_events(int(s), refs)
        for lab, sel in (("tail S>u (kept)", ~ds["drop"] & (ds.S > U)), ("bulk S<=u (kept)", ~ds["drop"] & (ds.S <= U))):
            obs_d, obs_e, ctl_d, ctl_e = [], [], [], []
            for r in ds[sel].itertuples():
                e = ev.get(int(r.camera))
                if e is None:
                    continue
                a, b = epoch_distances(r.period, r.t0, e)
                obs_d += list(a); obs_e += list(b)
                for _ in range(5):
                    a, b = epoch_distances(r.period, e["tmin"] + rng.uniform(0, r.period), e)
                    ctl_d += list(a); ctl_e += list(b)
            obs_d, obs_e, ctl_d, ctl_e = map(np.asarray, (obs_d, obs_e, ctl_d, ctl_e))
            win = 0.1
            timing[f"{int(s)}|{lab}"] = dict(
                n_epochs=int(obs_d.size), frac_within_0p1d_dump=float((obs_d < win).mean()) if obs_d.size else None,
                control_frac_dump=float((ctl_d < win).mean()) if ctl_d.size else None,
                frac_within_0p5d_edge=float((obs_e < .5).mean()) if obs_e.size else None,
                control_frac_edge=float((ctl_e < .5).mean()) if ctl_e.size else None,
                n_dumps={int(c): int(v["dumps"].size) for c, v in ev.items()})
            if lab.startswith("tail"):
                bins = np.linspace(0, 2, 41)
                axs[0].hist(ctl_d, bins, density=True, histtype="step", color="#999999", label="control (random epoch)")
                axs[0].hist(obs_d, bins, density=True, histtype="step", color="#e45756", lw=1.5, label=f"null tail epochs (n={obs_d.size})")
                axs[1].hist(ctl_e, bins, density=True, histtype="step", color="#999999", label="control")
                axs[1].hist(obs_e, bins, density=True, histtype="step", color="#e45756", lw=1.5, label="null tail epochs")
    axs[0].set_xlabel("epoch distance to nearest momentum dump (d)"); axs[1].set_xlabel("epoch distance to nearest gap / orbit edge (d)")
    for a in axs:
        a.legend(fontsize=7); a.set_ylabel("density")
    fig.suptitle(f"Inverted-null peaks with S > {U:g}: timing vs spacecraft events", fontsize=9)
    fig.tight_layout(); fig.savefig(out_dir / "null_peak_timing.png"); plt.close(fig)
    res["timing"] = timing
    (out_dir / "diag.json").write_text(json.dumps(res, indent=1, default=float))
    print(json.dumps(res, indent=1, default=float))


if __name__ == "__main__":
    main()
