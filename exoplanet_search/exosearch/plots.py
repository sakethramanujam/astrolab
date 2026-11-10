"""Diagnostic sheet: raw+trend, periodogram, phase fold, odd/even, secondary."""
from __future__ import annotations

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from . import bls
from .task import prepare


def _fold(t, P, t0):
    return ((t - t0 + 0.5 * P) % P) - 0.5 * P


def _binned(x, y, nb=80, lo=None, hi=None):
    lo = x.min() if lo is None else lo
    hi = x.max() if hi is None else hi
    edges = np.linspace(lo, hi, nb + 1)
    idx = np.digitize(x, edges) - 1
    ok = (idx >= 0) & (idx < nb)
    s = np.bincount(idx[ok], y[ok], nb)
    n = np.bincount(idx[ok], minlength=nb)
    c = 0.5 * (edges[1:] + edges[:-1])
    m = n > 0
    return c[m], s[m] / n[m]


def diagnostic(path, row, out_png, dpi=110):
    d = prepare(path)
    t, f, _ = d["raw"]
    tf, ff, ef = d["flat"]
    tb, fb, eb = d["binned"]
    res = bls.search(tb, fb, eb, backend="cpu")
    P, t0, dur = row["period"], row["t0"], row["duration"]
    fig = plt.figure(figsize=(13, 9))
    gs = fig.add_gridspec(3, 3, hspace=0.45, wspace=0.3)
    ax = fig.add_subplot(gs[0, :])
    ax.plot(t, f, ",", color="0.5")
    ax.plot(t, d["trend"], "-", color="C1", lw=1)
    n = np.arange(np.ceil((t.min() - t0) / P), np.floor((t.max() - t0) / P) + 1)
    for tt in t0 + n * P:
        ax.axvline(tt, color="C3", alpha=0.3, lw=0.8)
    ax.set(xlabel="BTJD", ylabel="PDCSAP (norm)", title=(
        f"TIC {row['tic']}  S{row['sector']}  Tmag={row.get('tmag') or float('nan'):.2f}  "
        f"P={P:.5f} d  depth={row['depth']*1e6:.0f} ppm  SDE={row['sde']:.1f}  SNR={row['snr']:.1f}  "
        f"known={row.get('known_toi') or row.get('known_planet') or 'none'}"))
    ax = fig.add_subplot(gs[1, :])
    ax.plot(res["periods"], res["power_detrended"], lw=0.6)
    ax.axvline(P, color="C3", alpha=0.4)
    ax.set(xscale="log", xlabel="Period (d)", ylabel="BLS power (detrended)")
    ph = _fold(tf, P, t0) * 24
    w = max(3 * dur * 24, 3)
    ax = fig.add_subplot(gs[2, 0])
    m = np.abs(ph) < w
    ax.plot(ph[m], ff[m], ".", ms=1.5, color="0.6")
    ax.plot(*_binned(ph[m], ff[m], 60, -w, w), "o-", ms=3, color="C0")
    ax.set(xlabel="hours from mid-transit", ylabel="flux", title="phase fold")
    ax = fig.add_subplot(gs[2, 1])
    ep = np.round((tf - t0) / P).astype(int)
    for parity, c in ((0, "C0"), (1, "C2")):
        mm = m & (ep % 2 == parity)
        if mm.sum() > 5:
            ax.plot(*_binned(ph[mm], ff[mm], 40, -w, w), "o-", ms=3, color=c,
                    label="even" if parity == 0 else "odd")
    ax.legend(fontsize=8)
    ax.set(xlabel="hours", title=f"odd/even  {row['odd_even_sigma']:.1f}σ")
    ax = fig.add_subplot(gs[2, 2])
    ph2 = _fold(tf, P, t0 + 0.5 * P) * 24
    m2 = np.abs(ph2) < w
    ax.plot(ph2[m2], ff[m2], ".", ms=1.5, color="0.6")
    if m2.sum() > 5:
        ax.plot(*_binned(ph2[m2], ff[m2], 60, -w, w), "o-", ms=3, color="C1")
    ax.set(xlabel="hours from phase 0.5", title=f"secondary  {row['secondary_sigma']:.1f}σ")
    fig.text(0.01, 0.005, f"flags: {row.get('flags') or 'none'}   (screening output, not a validated planet)",
             fontsize=8)
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return out_png
