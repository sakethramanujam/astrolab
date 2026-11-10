"""Figures for the site's Learn page, made from our own sector-70 data.

    python -m exosearch.learn_figs [--sector 70] [--tic 301289516]
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from . import bls, config, lightcurve, store, targets
from .task import prepare

OUT = Path(__file__).resolve().parents[1] / "site" / "learn"
BG, FG, ACC, ACC2, MUTED = "#0b1020", "#e6e9f5", "#7cc4ff", "#ffb86b", "#5b6585"

plt.rcParams.update({
    "figure.facecolor": BG, "axes.facecolor": BG, "savefig.facecolor": BG, "axes.edgecolor": MUTED,
    "axes.labelcolor": FG, "xtick.color": FG, "ytick.color": FG, "text.color": FG,
    "axes.grid": True, "grid.color": "#1c2540", "grid.linewidth": 0.6, "font.size": 10,
    "axes.titleweight": "bold", "axes.spines.top": False, "axes.spines.right": False})


def _fits(sector, tic):
    man = targets.fetch_manifest(sector).set_index("tic")
    m = man.loc[tic]
    return targets.download_one(m.fname, m.url, config.lc_dir(sector))


def _save(fig, name):
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT / name, dpi=90, bbox_inches="tight")
    plt.close(fig)
    print("wrote", OUT / name)


def _fold(t, P, t0):
    return ((t - t0 + 0.5 * P) % P) - 0.5 * P


def _bin(x, y, nb, lo, hi):
    e = np.linspace(lo, hi, nb + 1)
    i = np.digitize(x, e) - 1
    ok = (i >= 0) & (i < nb)
    s = np.bincount(i[ok], y[ok], nb)
    n = np.bincount(i[ok], minlength=nb)
    m = n > 0
    return (0.5 * (e[1:] + e[:-1]))[m], s[m] / n[m]


def star_figures(sector, tic, row):
    d = prepare(_fits(sector, tic))
    t, f, _ = d["raw"]
    tf, ff, ef = d["flat"]
    tb, fb, eb = d["binned"]
    P, t0, dur = row["period"], row["t0"], row["duration"]

    fig, ax = plt.subplots(figsize=(10, 3.2))
    ax.plot(t, f, ".", ms=1, color=MUTED, alpha=0.7, label="PDCSAP flux (2-min)")
    ax.plot(t, d["trend"], color=ACC2, lw=1.4, label="wotan biweight trend (0.75 d window)")
    ax.set(xlabel="time (BTJD = BJD - 2457000)", ylabel="relative brightness",
           title=f"1  A TESS light curve: TIC {tic}, sector {sector}")
    ax.legend(loc="lower left", fontsize=8, frameon=False)
    _save(fig, "01_lightcurve.png")

    fig, ax = plt.subplots(figsize=(10, 3.2))
    ax.plot(tf, ff, ".", ms=1, color=MUTED, alpha=0.6, label="detrended (flux / trend)")
    ax.plot(tb, fb, ".", ms=2.5, color=ACC, label="binned to 10 min")
    n = np.arange(np.ceil((tf.min() - t0) / P), np.floor((tf.max() - t0) / P) + 1)
    for x in t0 + n * P:
        ax.axvline(x, color=ACC2, alpha=0.25, lw=1)
    ax.set(xlabel="time (BTJD)", ylabel="relative brightness",
           title="2  After detrending: flat, with transits (orange lines) left intact")
    ax.legend(loc="lower left", fontsize=8, frameon=False)
    _save(fig, "02_detrended.png")

    res = bls.search(tb, fb, eb, backend="cpu")
    p = res["power_detrended"]
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 3.4), gridspec_kw=dict(width_ratios=[3, 1]))
    a1.plot(res["periods"], p, color=ACC, lw=0.6)
    a1.axhline(p.mean(), color=FG, lw=0.8, ls="--")
    a1.axhline(p.mean() + 9 * p.std(), color=ACC2, lw=0.8, ls=":")
    a1.annotate(f"peak  P = {res['period']:.4f} d\nSDE = {res['sde']:.1f}", (res["period"], p.max()),
                xytext=(20, -30), textcoords="offset points", color=ACC2,
                arrowprops=dict(arrowstyle="->", color=ACC2))
    a1.set(xscale="log", xlabel="trial period (days)", ylabel="BLS power",
           title="3  BLS periodogram: how well a box-shaped dip fits at each period")
    a2.hist(p, bins=80, orientation="horizontal", color=MUTED)
    a2.axhline(p.mean(), color=FG, lw=0.8, ls="--", label="mean")
    a2.axhline(p.mean() + 9 * p.std(), color=ACC2, lw=0.8, ls=":", label="mean + 9 sigma")
    a2.axhline(p.max(), color=ACC, lw=1.2, label="peak")
    a2.set(xscale="log", xlabel="count", title="SDE = (peak - mean) / sigma")
    a2.legend(fontsize=7, frameon=False)
    _save(fig, "03_periodogram.png")

    ph = _fold(tf, P, t0) * 24
    w = max(4 * dur * 24, 4)
    m = np.abs(ph) < w
    fig, ax = plt.subplots(figsize=(7, 3.4))
    ax.plot(ph[m], (ff[m] - 1) * 1e6, ".", ms=1.5, color=MUTED, alpha=0.6)
    x, y = _bin(ph[m], (ff[m] - 1) * 1e6, 60, -w, w)
    ax.plot(x, y, "o", ms=4, color=ACC)
    box = np.where(np.abs(np.linspace(-w, w, 400)) < dur * 12, -row["depth"] * 1e6, 0)
    ax.plot(np.linspace(-w, w, 400), box, color=ACC2, lw=1.5, label="best-fit box")
    ax.set(xlabel="hours from mid-transit", ylabel="ppm", ylim=(-4 * row["depth"] * 1e6, 3 * row["depth"] * 1e6),
           title=f"4  Folded at {P:.4f} d: depth {row['depth']*1e6:.0f} ppm, {dur*24:.1f} h")
    ax.legend(fontsize=8, frameon=False)
    _save(fig, "04_fold.png")


def eb_figure(sector, row):
    d = prepare(_fits(sector, int(row["tic"])))
    tf, ff, _ = d["flat"]
    P, t0, dur = row["period"], row["t0"], row["duration"]
    w = max(3 * dur * 24, 4)
    fig, axs = plt.subplots(1, 2, figsize=(11, 3.4))
    ph = _fold(tf, P, t0) * 24
    ep = np.round((tf - t0) / P).astype(int)
    for par, c, lab in ((0, ACC, "even transits"), (1, ACC2, "odd transits")):
        m = (np.abs(ph) < w) & (ep % 2 == par)
        axs[0].plot(*_bin(ph[m], (ff[m] - 1) * 1e6, 50, -w, w), "o-", ms=3, color=c, label=lab)
    axs[0].set(xlabel="hours from mid-transit", ylabel="ppm",
               title=f"5a  Odd vs even: {row['odd_even_sigma']:.1f} sigma apart")
    axs[0].legend(fontsize=8, frameon=False)
    ph2 = _fold(tf, P, t0 + 0.5 * P) * 24
    m = np.abs(ph2) < w
    axs[1].plot(ph2[m], (ff[m] - 1) * 1e6, ".", ms=1, color=MUTED, alpha=0.5)
    axs[1].plot(*_bin(ph2[m], (ff[m] - 1) * 1e6, 50, -w, w), "o", ms=3, color=ACC2)
    axs[1].set(xlabel="hours from phase 0.5", title=f"5b  Halfway round the orbit: {row['secondary_sigma']:.1f} sigma dip")
    fig.suptitle(f"A rejected signal (TIC {int(row['tic'])}): flags = {row['flags']}", fontsize=9, color=MUTED, y=-0.02)
    _save(fig, "05_odd_even_secondary.png")


def sector_figures(df):
    ok = df[df.status == "ok"]
    fig, ax = plt.subplots(figsize=(10, 3.2))
    ax.hist(ok.period, bins=np.geomspace(0.5, 14, 80), color=ACC)
    ax.set(xscale="log", xlabel="best BLS period (days)", ylabel="stars",
           title="6  Best periods across every star: pile-ups are the spacecraft, not planets")
    for x in (6.4, 12.0):
        ax.axvline(x, color=ACC2, ls=":", lw=1)
    _save(fig, "06_systematics.png")

    fig, ax = plt.subplots(figsize=(10, 3.2))
    ax.hist(ok.sde.clip(0, 40), bins=80, color=MUTED, label="all stars")
    ax.hist(ok[ok.passed == 1].sde.clip(0, 40), bins=80, range=(0, 40), color=ACC2, label="passed every screen")
    ax.axvline(9, color=ACC, ls="--", label="SDE threshold 9")
    ax.set(yscale="log", xlabel="SDE", ylabel="stars", title="7  Most stars show nothing; a few stand out")
    ax.legend(fontsize=8, frameon=False)
    _save(fig, "07_sde.png")


def diff_image_figure(sector, df):
    """08: out-of-transit image vs difference image for an on-target and an off-target signal."""
    from . import deepvet
    dv = deepvet.frame(deepvet.connect(), sector)
    if dv.empty:
        return
    m = df.merge(dv, on=["tic", "sector"])
    m = m[m.diff_peak_snr.fillna(0) > 5]
    off = m[m.deep_flags.fillna("").str.contains("dip_off_target|centroid_offset")].sort_values("diff_peak_snr", ascending=False)
    on = m[(m.toi_period_match == 1) & (m.deep_flags.fillna("") == "")].sort_values("diff_peak_snr", ascending=False)
    picks = [("on target", on.iloc[0])] if len(on) else []
    if len(off):
        picks.append(("off target", off.iloc[0]))
    if not picks:
        return
    man = targets.fetch_manifest(sector).set_index("tic")
    fig, axs = plt.subplots(len(picks), 2, figsize=(9, 4.2 * len(picks)), squeeze=False)
    for row_i, (lab, r) in enumerate(picks):
        tp = deepvet.fetch_tpf(sector, man.loc[int(r.tic)].fname)
        res = deepvet.centroid_test(tp, r.period, r.t0, r.duration, images=True)
        for j, (img, c, title) in enumerate(((res["_oot"], res["_oot_c"], "out of transit"),
                                             (res["_diff"], res["_diff_c"], "light lost in transit"))):
            ax = axs[row_i, j]
            ax.imshow(img, origin="lower", cmap="magma")
            ax.plot(res["_oot_c"][1], res["_oot_c"][0], "+", color=ACC, ms=16, mew=2, label="target")
            if j == 1:
                ax.plot(c[1], c[0], "x", color=ACC2, ms=12, mew=2, label="dip centroid")
                ax.legend(fontsize=7, loc="upper right", frameon=False, labelcolor=FG)
            ax.set(title=f"TIC {int(r.tic)} ({lab}): {title}", xticks=[], yticks=[])
            ax.title.set_fontsize(9)
            ax.grid(False)
        axs[row_i, 1].text(0.02, 0.02, f"offset {res['centroid_offset_px']:.2f} px ({res['centroid_offset_sigma']:.1f}σ)",
                           transform=axs[row_i, 1].transAxes, color=FG, fontsize=8)
    _save(fig, "08_difference_image.png")
    for p in (config.sector_dir(sector) / "tp").glob("*.fits"):
        p.unlink()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sector", type=int, default=70)
    ap.add_argument("--tic", type=int, default=301289516)   # TOI-4517.01, a known planet
    a = ap.parse_args()
    df = store.frame(store.connect(), a.sector)
    row = df[df.tic == a.tic].iloc[0].to_dict()
    star_figures(a.sector, a.tic, row)
    eb = df[(df.status == "ok") & df["flags"].fillna("").str.contains("odd_even") &
            df["flags"].fillna("").str.contains("secondary") & (df.depth < 0.05) & (df.snr > 20)]
    if len(eb):
        eb_figure(a.sector, eb.sort_values("snr", ascending=False).iloc[0].to_dict())
    sector_figures(df)
    diff_image_figure(a.sector, df)


if __name__ == "__main__":
    main()
