"""PNG figures for the cross-check. Always writes a file, including when a fit was not possible."""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams.update({
    "figure.dpi": 120,
    "savefig.dpi": 140,
    "font.size": 10,
    "axes.grid": True,
    "grid.alpha": 0.3,
    "axes.axisbelow": True,
})

_BLUE = "#0072B2"
_ORANGE = "#E69F00"
_GREEN = "#009E73"
_VERM = "#D55E00"


def _save(fig, path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def _empty(path: Path, message: str):
    fig, ax = plt.subplots(figsize=(6.2, 3.6))
    ax.axis("off")
    ax.text(0.5, 0.5, message, ha="center", va="center", wrap=True)
    _save(fig, path)


def completeness_gamma(path: Path, bundle: dict):
    snr = np.asarray([] if bundle.get("snr") is None else bundle.get("snr"), float)
    det = np.asarray([] if bundle.get("detected") is None else bundle.get("detected"), float)
    gamma = bundle.get("gamma") or {}
    if snr.size == 0:
        _empty(path, "No injections in this split.\nGamma CDF not fit.")
        return
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    edges = np.array([0, 2, 4, 6, 8, 10, 12, 15, 20, 30, 45, 80], float)
    for i in range(len(edges) - 1):
        m = (snr >= edges[i]) & (snr < edges[i + 1]) & np.isfinite(snr)
        n = int(m.sum())
        if n == 0:
            continue
        k = float(det[m].sum())
        f = k / n
        # Wilson half-width is already in the JSON; a simple binomial SE is enough on the figure.
        se = np.sqrt(max(f * (1 - f), 0.25 / n) / n)
        xc = np.sqrt(max(edges[i], 1e-3) * edges[i + 1])
        ax.errorbar(xc, f, yerr=se, fmt="o", color=_BLUE, ms=5, lw=1)
    curve = gamma.get("curve") or {}
    xs = np.asarray(curve.get("snr") or [], float)
    if xs.size:
        ax.plot(xs, curve["p"], color=_ORANGE, lw=2, label="gamma CDF")
        if curve.get("p_lo"):
            ax.fill_between(xs, curve["p_lo"], curve["p_hi"], color=_ORANGE, alpha=0.25, lw=0, label="95% bootstrap")
    ax.set_xlim(0, 40)
    ax.set_ylim(-0.05, 1.05)
    ax.set_xlabel("expected S  (depth × √n_tr / CDPP)")
    ax.set_ylabel("detection fraction")
    ax.set_title("Completeness vs expected SNR")
    if gamma.get("shape") is not None:
        ax.legend(frameon=False)
    else:
        ax.text(0.98, 0.05, "gamma CDF not fit", transform=ax.transAxes, ha="right", color="#444")
    _save(fig, path)


def completeness_map(path: Path, bundle: dict):
    fraction = np.asarray(bundle.get("map_fraction") if bundle.get("map_fraction") is not None else np.full((8, 8), np.nan), float)
    counts = np.asarray(bundle.get("map_counts") if bundle.get("map_counts") is not None else np.zeros((8, 8)), float)
    if fraction.size == 0:
        _empty(path, "No period–radius map.")
        return
    fig, ax = plt.subplots(figsize=(6.2, 5.2))
    shown = np.ma.array(fraction, mask=~np.isfinite(fraction) | (counts <= 0))
    im = ax.imshow(shown, origin="lower", aspect="auto", vmin=0, vmax=1, cmap="viridis")
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="detected / injected")
    p_edges = np.asarray(bundle.get("p_edges"), float)
    r_edges = np.asarray(bundle.get("r_edges"), float)
    ax.set_xticks(range(len(p_edges) - 1))
    ax.set_xticklabels([f"{v:.2g}" for v in p_edges[:-1]], rotation=45, ha="right")
    ax.set_yticks(range(len(r_edges) - 1))
    ax.set_yticklabels([f"{v:.2g}" for v in r_edges[:-1]])
    ax.set_xlabel("period bin lower edge (days)")
    ax.set_ylabel("radius bin lower edge (R⊕)")
    ax.set_title("Detection completeness")
    _save(fig, path)


def null_distributions(path: Path, bundle: dict):
    by_mode = bundle.get("null_by_mode") or {}
    present = {k: np.asarray(v, float) for k, v in by_mode.items() if v is not None and len(v)}
    if not present:
        _empty(path, "No null maxima in this split.")
        return
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    edges = np.linspace(0, 30, 31)
    colors = {"invert": _BLUE, "scramble": _ORANGE, "combined": _GREEN}
    for name, sample in present.items():
        ax.hist(sample, bins=edges, histtype="step", density=True, lw=1.8,
                color=colors.get(name, _VERM), label=f"{name} (n={sample.size})")
    ax.set_xlabel("per-star maximum S")
    ax.set_ylabel("density")
    ax.set_title("Null distribution of maximum S")
    ax.legend(frameon=False)
    if "invert" in present and "scramble" not in present:
        ax.text(0.98, 0.95, "no scrambled null in this file", transform=ax.transAxes, ha="right", va="top", color="#444")
    _save(fig, path)


def tail_stability(path: Path, bundle: dict):
    rows = bundle.get("stability") or []
    mrl = bundle.get("mrl") or {}
    if not rows and not mrl.get("u"):
        _empty(path, "Tail stability was not assessed.\nToo few null exceedances to fit a generalized Pareto.")
        return
    fig, axes = plt.subplots(1, 3, figsize=(9.6, 3.6))
    u = [r.get("u") for r in rows if r.get("xi") is not None]
    xi = [r.get("xi") for r in rows if r.get("xi") is not None]
    ss = [r.get("s_star") for r in rows if r.get("s_star") is not None and r.get("xi") is not None]
    us = [r.get("u") for r in rows if r.get("s_star") is not None and r.get("xi") is not None]
    axes[0].plot(mrl.get("u") or [], mrl.get("mean_excess") or [], color=_BLUE)
    if mrl.get("se"):
        uu = np.asarray(mrl["u"], float)
        ee = np.asarray(mrl["mean_excess"], float)
        se = np.asarray(mrl["se"], float)
        axes[0].fill_between(uu, ee - se, ee + se, color=_BLUE, alpha=0.2, lw=0)
    axes[0].set_title("mean residual life")
    axes[0].set_xlabel("threshold u")
    axes[0].set_ylabel("mean excess")
    if xi:
        axes[1].plot(u, xi, "o-", color=_ORANGE)
        axes[1].set_xlabel("threshold u")
        axes[1].set_ylabel("ξ")
    else:
        axes[1].axis("off")
        axes[1].text(0.5, 0.5, "GPD shape not fit\n(fewer than 15 exceedances)", ha="center", va="center")
    axes[1].set_title("GPD shape vs u")
    if ss:
        axes[2].plot(us, ss, "o-", color=_VERM)
        axes[2].set_xlabel("threshold u")
        axes[2].set_ylabel("S*")
    else:
        axes[2].axis("off")
        axes[2].text(0.5, 0.5, "S* not fit from the GPD", ha="center", va="center")
    axes[2].set_title("S* vs u")
    fig.suptitle(bundle.get("stability_title") or "Threshold stability", fontsize=11)
    _save(fig, path)


def roc(path: Path, bundle: dict):
    curves = bundle.get("roc_curves") or []
    curves = [c for c in curves if c.get("fpr")]
    if not curves:
        _empty(path, "ROC not computed.\nNeed injections and null maxima in this split.")
        return
    fig, ax = plt.subplots(figsize=(5.6, 5.2))
    ax.plot([0, 1], [0, 1], ls="--", color="#888", lw=1)
    for c in curves:
        label = c["name"]
        if c.get("auc") is not None:
            label += f"  AUC {c['auc']:.2f}"
        ax.plot(c["fpr"], c["tpr"], lw=1.6, label=label)
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.02, 1.02)
    ax.set_xlabel("false-alarm rate (null maxima)")
    ax.set_ylabel("true-positive rate (injections)")
    ax.set_title("ROC on S")
    ax.legend(frameon=False, fontsize=8)
    _save(fig, path)


def reliability(path: Path, bundle: dict):
    rows = bundle.get("reliability_rows") or []
    if not rows:
        _empty(path, "Vetting effectiveness E is undefined.\nNo simulated null cleared S* with checks logged.")
        return
    fig, ax = plt.subplots(figsize=(6.8, 4.0))
    y = np.arange(len(rows))
    e = [r["E"] if r["E"] is not None else np.nan for r in rows]
    lo = [r["E"] - r["lo"] if r["E"] is not None and r["lo"] is not None else np.nan for r in rows]
    hi = [r["hi"] - r["E"] if r["E"] is not None and r["hi"] is not None else np.nan for r in rows]
    ax.barh(y, e, color=_BLUE, alpha=0.85, xerr=[lo, hi], capsize=3)
    ax.set_yticks(y)
    ax.set_yticklabels([r["label"] for r in rows], fontsize=8)
    ax.set_xlim(0, 1.05)
    ax.set_xlabel("vetting effectiveness E")
    ax.set_title("E from simulated nulls (Wilson 95%)")
    _save(fig, path)


def write_all(out_dir: Path, bundle: dict):
    fig = Path(out_dir) / "figures"
    completeness_gamma(fig / "completeness_gamma.png", bundle)
    completeness_map(fig / "completeness_map.png", bundle)
    null_distributions(fig / "null_distributions.png", bundle)
    tail_stability(fig / "tail_stability.png", bundle)
    roc(fig / "roc.png", bundle)
    reliability(fig / "reliability_effectiveness.png", bundle)
