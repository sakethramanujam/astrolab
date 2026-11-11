"""v1.1 light-curve screens. Same checks as v1.0.0 (exosearch/vet.py) except:

* odd_even: flagged only when |d_odd - d_even| / hypot(e_odd, e_even) > 3 AND
  |d_odd - d_even| / mean(d_odd, d_even) > 10%   (docs/v1.1_odd_even_secondary.md, agreed with the Scientist);
* secondary: flagged only when d_sec / e_sec > 3 AND d_sec / d_primary > 10%;
* momentum_dump: a transit epoch whose centre lies within 0.1 d of a reaction-wheel momentum dump (QUALITY bit 32,
  "desaturation event") is removed from depth and N_tr before S (exosearch/v11/search.py). The signal is rejected
  (`momentum_dump`) when it has <= 2 epochs with data and any of them is dump-adjacent; otherwise a warning
  (`momentum_dump_epoch`). QUALITY bit 1024 is not a dump flag in SPOC products (collateral cosmic ray) and is not used;
* period_ambiguous (tie_break.md amendment 1a): number of transits with data <= 2;
* the depth used for Rp, too_deep and no_dip is the exact box depth of the statistic (dump epochs removed).
There is no S threshold here: v1.1 S* is not frozen, so `low_s` is never set and nothing is called detected.
"""
from __future__ import annotations

import numpy as np

from .. import vet as vet10

R_EARTH_PER_RSUN = vet10.R_EARTH_PER_RSUN


def dump_times(t_all, quality, bit=32, gap=0.05):
    """Centres of contiguous runs of cadences flagged with the momentum-dump bit."""
    t_all = np.asarray(t_all, float)
    m = ((np.asarray(quality, int) & bit) > 0) & np.isfinite(t_all)
    td = np.sort(t_all[m])
    if td.size == 0:
        return np.array([])
    groups = np.split(td, np.flatnonzero(np.diff(td) > gap) + 1)
    return np.array([float(g.mean()) for g in groups])


def screens(st, best, meta, cfg):
    v = cfg["vet"]
    th = vet10.THRESH
    depth = float(best["depth_exact"])
    _, depth_err = map(float, st["depth"])
    d_odd, e_odd = map(float, st["depth_odd"])
    d_even, e_even = map(float, st["depth_even"])
    d_sec, e_sec = map(float, st["depth_phased"])
    eo = np.hypot(e_odd, e_even)
    oe = abs(d_odd - d_even) / eo if eo > 0 else 0.0
    mean_oe = 0.5 * (d_odd + d_even)
    oe_rel = abs(d_odd - d_even) / abs(mean_oe) if mean_oe != 0 else np.inf
    sec_sig = d_sec / e_sec if e_sec > 0 else 0.0
    sec_rel = d_sec / depth if depth > 0 else np.inf
    n_tr = int(np.sum(np.asarray(st["per_transit_count"]) > 0))
    P, dur = best["P_grid"], best["duration"]
    t_exp = vet10.expected_duration(P, meta.get("radius"))
    dur_ratio = dur / t_exp if t_exp > 0 else np.nan
    rstar = meta.get("radius") or np.nan
    rp = np.sqrt(max(depth, 0)) * rstar * R_EARTH_PER_RSUN if (rstar == rstar and depth == depth) else np.nan
    flags, warns = [], []
    if n_tr < th["min_transits"]:
        flags.append("few_transits")
    if oe > v["odd_even_sigma"] and oe_rel > v["odd_even_rel"]:
        flags.append("odd_even")
    if sec_sig > v["secondary_sigma"] and sec_rel > v["secondary_rel"]:
        flags.append("secondary")
    if depth > th["max_depth"]:
        flags.append("too_deep")
    if dur_ratio == dur_ratio and dur_ratio > th["max_dur_ratio"]:
        flags.append("too_long")
    if sec_sig < -th["variable_sigma"]:
        flags.append("variable")
    if rp == rp and rp > th["max_rp_earth"]:
        flags.append("rp_too_large")
    if not (depth > 0):
        flags.append("no_dip")
    if best.get("n_dump_epochs", 0) > 0:
        if best.get("n_epochs", 0) <= cfg["momentum_dump"]["reject_max_epochs"]:
            flags.append("momentum_dump")
        else:
            warns.append("momentum_dump_epoch")
    return dict(depth=depth, depth_err=depth_err, odd_even_sigma=float(oe), odd_even_rel=float(oe_rel),
                depth_odd=d_odd, depth_even=d_even, secondary_depth=d_sec, secondary_sigma=float(sec_sig),
                secondary_rel=float(sec_rel), n_transits=n_tr, period_ambiguous=bool(n_tr <= 2),
                dur_expected=float(t_exp), dur_ratio=float(dur_ratio), rp_earth=float(rp),
                flags=",".join(flags), warnings=",".join(warns), passed_screens=int(not flags))
