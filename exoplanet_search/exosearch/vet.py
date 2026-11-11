"""Simple, transparent vetting metrics. These are screens, not validation."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

R_EARTH_PER_RSUN = 109.1

THRESH = dict(min_transits=3, odd_even_sigma=3.0, secondary_sigma=3.0,
              max_depth=0.05, max_dur_ratio=1.5, max_rp_earth=25.0, variable_sigma=5.0,
              common_period_tol=0.005, common_period_min=4)


def expected_duration(period_d, radius_rsun):
    """Central-transit T14 (days) for a star of given radius; M ~ R on the main sequence."""
    r = radius_rsun if radius_rsun and np.isfinite(radius_rsun) and radius_rsun > 0 else 1.0
    m = min(max(r, 0.1), 1.5) if r < 2 else 1.2      # crude mass guess; giants ~1.2 Msun
    rho = m / r ** 3
    return (13.0 / 24.0) * (period_d / 365.25) ** (1 / 3) * rho ** (-1 / 3)


THRESHOLD_FILE = Path(__file__).resolve().parents[1] / "configs" / "detection_threshold.json"


def s_star(sector=None, camera=None):
    """Detection threshold on S for (sector, camera); derived from inversions, so it lives outside THRESH
    (and outside the frozen config hash). Returns (value, provisional)."""
    cfg = json.loads(THRESHOLD_FILE.read_text())
    v = cfg.get("per_sector_camera", {}).get(f"{sector}:{camera}", cfg["default"])
    return float(v), bool(cfg.get("provisional", True))


def vet(st, bls_out, meta, t, f, s_min=None):
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
    if s_min is None:
        s_min = s_star(meta.get("sector"), meta.get("camera"))[0]
    if bls_out.get("S", 0.0) < s_min:
        flags.append("low_s")               # detection threshold (SDE and white-noise SNR are logged only)
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
    if sec_sig < -THRESH["variable_sigma"]:
        flags.append("variable")            # brightening at phase 0.5: sinusoidal/rotational
    if rp == rp and rp > THRESH["max_rp_earth"]:
        flags.append("rp_too_large")        # > ~2.2 R_Jup: stellar companion or giant host
    if depth <= 0:
        flags.append("no_dip")
    return dict(depth=depth, depth_err=depth_err, snr=float(snr), odd_even_sigma=float(oe),
                depth_odd=d_odd, depth_even=d_even, secondary_depth=d_sec, secondary_sigma=float(sec_sig),
                n_transits=n_tr, dur_expected=float(t_exp), dur_ratio=float(dur_ratio),
                rp_earth=float(rp), flags=",".join(flags), passed=int(not flags))


def flag_common_periods(df):
    """Flag periods shared by many unrelated stars in a sector (spacecraft systematics, e.g.
    momentum dumps / orbit-gap aliases). Returns a copy with 'common_period' added to flags."""
    import pandas as pd  # noqa: F401
    df = df.copy()
    ok = df.status == "ok"
    P = df.loc[ok, "period"].to_numpy()
    if len(P) < 50:
        return df
    tol = THRESH["common_period_tol"]
    near = np.array([(np.abs(P / p - 1) < tol).sum() - 1 for p in P])
    broad = np.array([(np.abs(P / p - 1) < 20 * tol).sum() - 1 for p in P])
    expected = (broad - near) / 19.0             # background rate from the surrounding 20x window
    lim = np.maximum(THRESH["common_period_min"], 4 * expected + 3)
    hit = df.index[ok][near >= lim]
    for i in hit:
        f = df.at[i, "flags"] or ""
        if "common_period" not in f:
            df.at[i, "flags"] = ",".join(x for x in (f, "common_period") if x)
        df.at[i, "passed"] = 0
    return df
