"""v1.2 step 4: T14 retention diagnostic for the standard biweight detrend.  REPORTS ONLY; it is not a cut.

The detrend (`lightcurve.detrend`, wotan biweight, window 0.75 d, gap_trim 0.3 d) can partly absorb long transits.
This diagnostic plants box transits of known T14 and depth, runs the standard detrend unchanged and reports how
much of the injected depth survives, per T14 bin.  It returns numbers only: no pass/fail, and it never touches a
vetting flag, `final_verdict`, S or S*.

Retention estimator (per planted light curve)
    flat = lightcurve.detrend(t, f, e, window=0.75, gap_trim=0.3)   (cadences dropped by the detrend are gone)
    in-transit: |((t - t0 + P/2) mod P) - P/2| < T14/2  (known epoch t0, period P, true T14)
    depth_rec = mean_w(flat out of transit) - mean_w(flat in transit), inverse-variance weights 1/e^2
    retention = depth_rec / depth_injected          (1 = fully kept, 0 = fully absorbed; may exceed 1 with noise)
A planted transit is "kept" when retention >= `tolerance` (default 0.5).  A case with no in-transit (or no
out-of-transit) cadence after the detrend has NaN retention, is counted in `n` but not `n_valid`, and is not kept.
Bins (true T14, hours, edges 3 and 6): `T14 < 3 h`, `3 <= T14 <= 6 h`, `T14 > 6 h`.
"""
from __future__ import annotations

from typing import Any, Dict, Iterable, Optional, Sequence

import numpy as np

from .lightcurve import detrend

DETREND_WINDOW_D = 0.75
DEFAULT_EDGES_H = (3.0, 6.0)
DEFAULT_TOLERANCE = 0.5


def t14_cfg(cfg: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
    """Resolve bin edges / tolerance from `cfg["t14_retention"]`; missing keys give the documented defaults."""
    out = dict(bin_edges_h=list(DEFAULT_EDGES_H), tolerance=DEFAULT_TOLERANCE)
    blk = cfg.get("t14_retention") if isinstance(cfg, dict) else None
    if isinstance(blk, dict):
        if blk.get("bin_edges_h") is not None:
            out["bin_edges_h"] = [float(x) for x in blk["bin_edges_h"]]
        if blk.get("tolerance") is not None:
            out["tolerance"] = float(blk["tolerance"])
    return out


def bin_labels(edges_h: Sequence[float] = DEFAULT_EDGES_H):
    lo, hi = edges_h
    f = lambda x: f"{x:g}"
    return (f"T14 < {f(lo)} h", f"{f(lo)} <= T14 <= {f(hi)} h", f"T14 > {f(hi)} h")


def t14_bin(t14_h: float, edges_h: Sequence[float] = DEFAULT_EDGES_H) -> str:
    """Bin label for a true T14 in hours: `< lo`, `lo <= T14 <= hi`, `> hi`."""
    lo, hi = edges_h
    labels = bin_labels(edges_h)
    return labels[0] if t14_h < lo else (labels[2] if t14_h > hi else labels[1])


def plant_box(t, f, t0: float, period: float, t14_h: float, depth: float):
    """Copy of `f` with a periodic box transit (subtracts `depth` in transit) of duration T14."""
    t = np.asarray(t, float)
    out = np.array(f, float, copy=True)
    ph = ((t - t0 + 0.5 * period) % period) - 0.5 * period
    out[np.abs(ph) < t14_h / 48.0] -= depth
    return out


def box_depth(t, flat, e, t0: float, period: float, t14_h: float) -> float:
    """Inverse-variance out-of-transit minus in-transit mean of the detrended flux (NaN if either is empty)."""
    t, flat, e = (np.asarray(x, float) for x in (t, flat, e))
    ph = ((t - t0 + 0.5 * period) % period) - 0.5 * period
    inn = np.abs(ph) < t14_h / 48.0
    if not inn.any() or inn.all():
        return float("nan")
    w = 1.0 / e ** 2
    mean = lambda m: float(np.sum(w[m] * flat[m]) / np.sum(w[m]))
    return mean(~inn) - mean(inn)


def retention_one(t, f, e, t0: float, period: float, t14_h: float, depth: float) -> float:
    """Retention (recovered / injected depth) of one planted light curve after the standard detrend."""
    tt, flat, ee, _, _ = detrend(np.asarray(t, float), np.asarray(f, float), np.asarray(e, float),
                                 window=DETREND_WINDOW_D)
    return box_depth(tt, flat, ee, t0, period, t14_h) / depth


def t14_retention_diagnostic(cases: Iterable[Dict[str, Any]], tolerance: float = DEFAULT_TOLERANCE,
                             edges_h: Sequence[float] = DEFAULT_EDGES_H) -> Dict[str, Any]:
    """Per-T14-bin retention breakdown for planted transits.

    `cases`: dicts with t, f (flux that already contains the planted box), e, t0, period, t14_h, depth.
    Returns dict(labels, edges_h, tolerance, estimator, bins={label: dict(n, n_valid, n_kept, frac_kept,
    mean_retention)}, retentions=[(t14_h, label, retention)]).  Diagnostic only: there is no pass/fail field.
    """
    labels = bin_labels(edges_h)
    acc = {lab: [] for lab in labels}
    rows = []
    for c in cases:
        r = retention_one(c["t"], c["f"], c["e"], c["t0"], c["period"], c["t14_h"], c["depth"])
        lab = t14_bin(c["t14_h"], edges_h)
        acc[lab].append(r)
        rows.append((float(c["t14_h"]), lab, float(r)))
    bins = {}
    for lab in labels:
        r = np.array(acc[lab], float)
        ok = np.isfinite(r)
        kept = int(np.sum(ok & (r >= tolerance)))
        bins[lab] = dict(n=int(len(r)), n_valid=int(ok.sum()), n_kept=kept,
                         frac_kept=(kept / int(ok.sum())) if ok.any() else None,
                         mean_retention=float(r[ok].mean()) if ok.any() else None)
    return dict(labels=list(labels), edges_h=[float(x) for x in edges_h], tolerance=float(tolerance),
                estimator="retention = (IVW mean OOT - IVW mean in-transit of detrended flux) / injected depth",
                bins=bins, retentions=rows)
