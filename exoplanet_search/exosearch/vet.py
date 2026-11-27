"""Simple, transparent vetting metrics. These are screens, not validation."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

R_EARTH_PER_RSUN = 109.1

THRESH = dict(min_transits=3, odd_even_sigma=3.0, secondary_sigma=3.0,
              max_depth=0.05, max_dur_ratio=1.5, max_rp_earth=25.0, variable_sigma=5.0,
              common_period_tol=0.005, common_period_min=4)

# --- v1.2 step 5c: the 5% depth screen as a setting (configs/v1.1.json "vet"."too_deep") ---
TOO_DEEP = "too_deep"
TOO_DEEP_MODES = ("cut", "flag", "off")
# Defaults reproduce v1.1 exactly: depth > 0.05 is a rejection (the strict comparison keeps depth == 0.05).
DEFAULT_TOO_DEEP = dict(mode="cut", max_depth=THRESH["max_depth"])


def too_deep_cfg(cfg=None, too_deep=None):
    """Resolve the too_deep depth-screen settings from `cfg["vet"]["too_deep"]` (and optional override).

    Keys:
      mode:       'cut'  - depth > max_depth is a rejection (default; exactly v1.1);
                  'flag' - the screen is evaluated and recorded, but it is never a rejection;
                  'off'  - the screen is not applied (and a stored too_deep flag is not a rejection).
      max_depth:  depth threshold (default 0.05, vet.THRESH["max_depth"], the v1.1 5% cap).

    `too_deep` may be a mode string or a dict with `mode` / `max_depth`.  Configs that omit the keys
    resolve to DEFAULT_TOO_DEEP ('cut', 0.05), so v1.1 rows and existing configs are unchanged.
    The comparison is strict (`depth > max_depth`): depth == max_depth is kept.
    """
    out = dict(DEFAULT_TOO_DEEP)
    v = cfg.get("vet") if isinstance(cfg, dict) else None
    sources = []
    if isinstance(v, dict) and v.get("too_deep") is not None:
        sources.append(v["too_deep"])
    if too_deep is not None:
        sources.append(too_deep)
    for src in sources:
        if isinstance(src, str):
            out["mode"] = src
        elif isinstance(src, dict):
            if src.get("mode") is not None:
                out["mode"] = str(src["mode"])
            if src.get("max_depth") is not None:
                out["max_depth"] = float(src["max_depth"])
        else:
            raise ValueError(f"too_deep must be a mode string or a dict, got {src!r}")
    out["mode"] = str(out["mode"]).strip().lower()
    if out["mode"] not in TOO_DEEP_MODES:
        raise ValueError(f"too_deep mode must be one of {TOO_DEEP_MODES}, got {out['mode']!r}")
    if not np.isfinite(out["max_depth"]) or out["max_depth"] < 0:
        raise ValueError(f"too_deep max_depth must be finite and >= 0, got {out['max_depth']!r}")
    return out


def too_deep(depth, cfg=None, too_deep=None):
    """Evaluate the depth screen for one depth: -> dict(mode, max_depth, flagged, reject).

    `flagged` is the strict `depth > max_depth` result (NaN / missing depth -> False); `reject` is
    `flagged` only in mode 'cut'.  Mode 'flag' still reports `flagged` (so it can be recorded) but
    never rejects; mode 'off' reports the screen as not applied through `mode` and sets both False.
    This is the ONE implementation, shared by the search-time screen (vet11.screens) and the final
    verdict (vet_final), which re-evaluates from the row's depth so a rescore may use another
    threshold or mode than the search-time flag.
    """
    c = too_deep_cfg(cfg, too_deep)
    try:
        d = float(depth)
    except (TypeError, ValueError):
        d = np.nan
    flagged = bool(np.isfinite(d) and d > c["max_depth"])
    return dict(mode=c["mode"], max_depth=c["max_depth"],
                flagged=bool(flagged and c["mode"] != "off"),
                reject=bool(flagged and c["mode"] == "cut"))


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
