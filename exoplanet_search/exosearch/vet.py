"""Simple, transparent vetting metrics. These are screens, not validation."""
from __future__ import annotations

import numpy as np

R_EARTH_PER_RSUN = 109.1

THRESH = dict(sde=9.0, snr=7.0, min_transits=2, odd_even_sigma=3.0, secondary_sigma=3.0,
              max_depth=0.05, max_dur_ratio=1.5)


def expected_duration(period_d, radius_rsun):
    """Central-transit T14 (days) for a star of given radius; M ~ R on the main sequence."""
    r = radius_rsun if radius_rsun and np.isfinite(radius_rsun) and radius_rsun > 0 else 1.0
    m = min(max(r, 0.1), 1.5) if r < 2 else 1.2      # crude mass guess; giants ~1.2 Msun
    rho = m / r ** 3
    return (13.0 / 24.0) * (period_d / 365.25) ** (1 / 3) * rho ** (-1 / 3)


def vet(st, bls_out, meta, t, f):
    depth, depth_err = map(float, st["depth"])
    d_odd, e_odd = map(float, st["depth_odd"])
    d_even, e_even = map(float, st["depth_even"])
    d_sec, e_sec = map(float, st["depth_phased"])   # depth at phase 0.5 = secondary eclipse
    snr = depth / depth_err if depth_err > 0 else 0.0
    oe = abs(d_odd - d_even) / np.hypot(e_odd, e_even) if np.hypot(e_odd, e_even) > 0 else 0.0
    sec_sig = d_sec / e_sec if e_sec > 0 else 0.0
    n_tr = int(np.sum(np.asarray(st["per_transit_count"]) > 0))
    P, dur = bls_out["period"], bls_out["duration"]
    t_exp = expected_duration(P, meta.get("radius"))
    dur_ratio = dur / t_exp if t_exp > 0 else np.nan
    rstar = meta.get("radius") or np.nan
    rp = np.sqrt(max(depth, 0)) * rstar * R_EARTH_PER_RSUN if rstar == rstar else np.nan

    flags = []
    if bls_out["sde"] < THRESH["sde"]:
        flags.append("low_sde")
    if snr < THRESH["snr"]:
        flags.append("low_snr")
    if n_tr < THRESH["min_transits"]:
        flags.append("few_transits")
    if oe > THRESH["odd_even_sigma"]:
        flags.append("odd_even")
    if sec_sig > THRESH["secondary_sigma"]:
        flags.append("secondary")
    if depth > THRESH["max_depth"]:
        flags.append("too_deep")
    if dur_ratio == dur_ratio and dur_ratio > THRESH["max_dur_ratio"]:
        flags.append("too_long")
    if depth <= 0:
        flags.append("no_dip")
    return dict(depth=depth, depth_err=depth_err, snr=float(snr), odd_even_sigma=float(oe),
                depth_odd=d_odd, depth_even=d_even, secondary_depth=d_sec, secondary_sigma=float(sec_sig),
                n_transits=n_tr, dur_expected=float(t_exp), dur_ratio=float(dur_ratio),
                rp_earth=float(rp), flags=",".join(flags), passed=int(not flags))
