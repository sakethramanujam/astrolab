"""Amended light-curve vetting, applied after the search (v1.1 vetter amendment 2026-11-12).

docs/preregistration/v1.1_amendment_vetter_2026-11-12.md. Three rule changes from the Kepler DR25 Robovetter
(Thompson et al. 2018, ApJS 235, 38); every other check and threshold is unchanged:
  (a) secondary-eclipse exemption: a secondary (the phase-0.5 screen 'secondary' or the deep 'secondary_scan') does
      not reject when the geometric albedo it implies, A_g = depth_sec / (Rp/a)^2, is below 1 and Rp < 30 R_Earth;
  (b) radius cap 30 R_Earth, replacing both earlier caps (deep 'rp_gt_2rjup', screen 'rp_too_large' at 25 R_Earth);
  (c) stellar-density mismatch is a flag, not a rejection;
  (d) leave-one-orbit-out secondary test: a 'secondary_scan' rejection stands only if the scan still clears its 5 sigma
      threshold after dropping any single orbit (minimum over the leave-one-orbit-out scans >= 5).
The search code and its rows are unchanged (search-path hash unchanged): this layer reads the logged row fields and,
only for rows with a secondary flag, rebuilds the light curve the pipeline vetted (same load, null transform and
prepare) to measure the secondary-scan depth at the scan's own phase. Output column: passed_lc_checks_v2.
"""
from __future__ import annotations

from typing import Any, Dict, Optional

import numpy as np
import pandas as pd

from .. import deepvet

AMENDMENT = "docs/preregistration/v1.1_amendment_vetter_2026-11-12.md"
R_EARTH_PER_RJUP = 11.209
RP_CAP_REARTH = 30.0
ALBEDO_MAX = 1.0
G_SI = 6.674e-11
RHO_SUN_KG_M3 = 1408.0
LOO_SIGMA = deepvet.THRESH["sec_scan_sigma"]                         # (d): same 5 sigma as the scan itself
DROPPED_DEEP = {"rp_gt_2rjup", "density_mismatch"}                  # (b) replaced, (c) flag only
REJECT_V2 = (set(deepvet.REJECT) - DROPPED_DEEP) | {"rp_gt_30re"}
SECONDARY = {"secondary", "secondary_scan"}


def _f(x) -> float:
    try:
        v = float(x)
    except (TypeError, ValueError):
        return np.nan
    return v


def _set(s) -> set:
    return set(str(s).split(",")) - {"", "nan", "None"} if s is not None and s == s else set()


def a_over_rstar(P: float, rho: Optional[float], dur: Optional[float] = None):
    """a/R* from the TIC stellar density (solar units, as stellar.tic_rho; circular Kepler orbit); else
    P / (pi * duration). -> (value, source)."""
    r = _f(rho)
    if r > 0 and P > 0:
        rho_si = r * RHO_SUN_KG_M3
        return float((G_SI * rho_si * (P * 86400.0) ** 2 / (3.0 * np.pi)) ** (1.0 / 3.0)), "tic_density"
    d = _f(dur)
    if d > 0 and P > 0:
        return float(P / (np.pi * d)), "duration"
    return np.nan, "none"


def implied_albedo(depth_sec: float, depth: float, a_rs: float) -> float:
    """A_g = depth_sec / (Rp/a)^2 with (Rp/a)^2 = depth / (a/R*)^2."""
    if not (depth > 0 and a_rs > 0) or not np.isfinite(depth_sec):
        return np.nan
    return float(depth_sec / (depth / a_rs ** 2))


def rp_earth(row: Dict[str, Any]) -> float:
    rj = _f(row.get("rp_rjup"))
    return rj * R_EARTH_PER_RJUP if np.isfinite(rj) else _f(row.get("rp_earth"))


def scan_depth(t, f, e, P: float, t0: float) -> Dict[str, float]:
    """deepvet.shape_metrics plus the box depth at the secondary scan's best phase (same base, box and phase)."""
    sm = deepvet.shape_metrics(t, f, e, P, t0)
    out = dict(sec_scan_sigma=sm["sec_scan_sigma"], sec_scan_phase=sm["sec_scan_phase"], sec_scan_depth=np.nan,
               duration_fine=sm["duration_fine"], t0_fine=sm["t0_fine"])
    p0 = sm["sec_scan_phase"]
    if p0 is None:
        return out
    dur, t0f = sm["duration_fine"], sm["t0_fine"]
    x = (((t - t0f + 0.5 * P) % P) - 0.5 * P) / dur
    oot = (np.abs(x) > 1.0) & (np.abs(x) < 3.0)
    base = np.median(f[oot]) if oot.sum() > 10 else 1.0
    ph = ((t - t0f) / P) % 1.0
    m = np.abs(ph - p0) < (dur / P) / 2
    if m.any():
        out["sec_scan_depth"] = float(base - np.mean(f[m]))
    return out


def loo_scan_min(t, f, e, P: float, t0: float, t0_ref: float) -> Dict[str, float]:
    """(d): the secondary scan (deepvet.shape_metrics, unchanged) re-run with each single orbit left out; orbit k is
    t0_ref + k*P <= t < t0_ref + (k+1)*P (t0_ref: the full scan's refined epoch). -> min sigma and the orbit giving it."""
    k = np.floor((t - t0_ref) / P).astype(int)
    best, which, n = np.inf, None, 0
    for o in np.unique(k):
        m = k != o
        if m.sum() < 100:
            continue
        s = float(deepvet.shape_metrics(t[m], f[m], e[m], P, t0)["sec_scan_sigma"])
        n += 1
        if s < best:
            best, which = s, int(o)
    return dict(loo_min_sigma=float(best) if n else np.nan, loo_min_orbit=which, loo_n_orbits=n)


def rebuild_flat(job: Dict[str, Any], row: Dict[str, Any]):
    """The detrended light curve the pipeline vetted for this row's task (exosearch.v11.pipeline.run_job)."""
    from .. import nulltest
    from . import pipeline, scramble
    t, f, e, good, meta, dumps = pipeline.load(job, download=False)
    tg, fg, eg = t[good], f[good], e[good]
    mode = row.get("mode", "real")
    if mode == "invert":
        tg, fg, eg = nulltest.transform(mode, tg, fg, eg, meta, int(row.get("seed") or 0))
    elif mode == "scramble":
        tg, fg, eg, _ = scramble.segment_rotate(tg, fg, eg, int(row.get("seed")), job["cfg"])
    elif mode != "real":
        raise ValueError(f"re-vetting not defined for mode {mode}")
    d = pipeline.prepare(tg, fg, eg, job["cfg"], True, None, gap_trim=0.0 if mode == "scramble" else 0.3)
    return d["flat"]


def amended(row: Dict[str, Any], scan: Optional[Dict[str, float]] = None, rho: Optional[float] = None) -> Dict[str, Any]:
    """Amended verdict for one row. scan: scan_depth() output (needed only when 'secondary_scan' fired)."""
    flags, deep = _set(row.get("flags")), _set(row.get("deep_flags"))
    P, depth = _f(row.get("P_grid")) if np.isfinite(_f(row.get("P_grid"))) else _f(row.get("period")), _f(row.get("depth"))
    rp = rp_earth(row)
    dur = (scan or {}).get("duration_fine", _f(row.get("duration")))
    a_rs, a_src = a_over_rstar(P, rho, dur)
    small = bool(np.isfinite(rp) and rp < RP_CAP_REARTH)
    alb_scr = implied_albedo(_f(row.get("secondary_depth")), depth, a_rs)
    alb_scan = implied_albedo((scan or {}).get("sec_scan_depth", np.nan), depth, a_rs)
    rej = set(flags) - {"rp_too_large"}
    if "secondary" in rej and small and np.isfinite(alb_scr) and alb_scr < ALBEDO_MAX:
        rej.discard("secondary")
    if np.isfinite(rp) and rp > RP_CAP_REARTH:
        rej.add("rp_gt_30re")
    drej = deep & REJECT_V2
    loo = _f((scan or {}).get("loo_min_sigma"))
    if "secondary_scan" in drej and small and np.isfinite(alb_scan) and alb_scan < ALBEDO_MAX:
        drej.discard("secondary_scan")                                     # (a)
    if "secondary_scan" in drej and np.isfinite(loo) and loo < LOO_SIGMA:
        drej.discard("secondary_scan")                                     # (d)
    rej |= drej
    warn = (deep - drej) | ({"density_mismatch"} & deep)
    return dict(passed_lc_checks_v2=bool(not rej), reject_v2=",".join(sorted(rej)), warn_v2=",".join(sorted(warn)),
                rp_earth_v2=rp, a_rs_v2=a_rs, a_rs_src_v2=a_src, albedo_secondary_v2=alb_scr,
                albedo_scan_v2=alb_scan, sec_scan_depth_v2=(scan or {}).get("sec_scan_depth", np.nan),
                sec_scan_phase_v2=(scan or {}).get("sec_scan_phase", np.nan), loo_min_sigma_v2=loo,
                loo_min_orbit_v2=(scan or {}).get("loo_min_orbit"))


def revet_one(r: Dict[str, Any], j: Dict[str, Any], cfg: Dict[str, Any]) -> Dict[str, Any]:
    """Amended vetting for one search row. j: its manifest job (light-curve file, TIC density); the null seed is the
    row's. A light curve is rebuilt only when the deep checks raised 'secondary_scan'."""
    scan = None
    if "secondary_scan" in _set(r.get("deep_flags")):
        t, f, e = rebuild_flat(dict(j, cfg=cfg), r)
        P, t0 = float(r["P_grid"]), float(r["t0"])
        scan = scan_depth(t, f, e, P, t0)
        if not np.isclose(scan["sec_scan_sigma"], _f(r.get("sec_scan_sigma")), rtol=1e-6, atol=1e-6):
            raise RuntimeError(f"TIC {r['tic']}: rebuilt secondary scan {scan['sec_scan_sigma']} != row "
                               f"{r.get('sec_scan_sigma')}")
        scan.update(loo_scan_min(t, f, e, P, t0, scan["t0_fine"]))
    return dict(r, **amended(r, scan, j.get("rho")))


def revet(rows: pd.DataFrame, jobs: Dict[int, Dict[str, Any]], cfg: Dict[str, Any], workers: int = 1) -> pd.DataFrame:
    """Amended vetting for search rows (jobs: manifest jobs by TIC)."""
    recs = rows.to_dict(orient="records")
    args = [(r, jobs.get(int(r["tic"]), {}), cfg) for r in recs]
    if workers > 1:
        from concurrent.futures import ProcessPoolExecutor
        with ProcessPoolExecutor(workers) as ex:
            out = list(ex.map(_revet_star, args, chunksize=4))
    else:
        out = [_revet_star(a) for a in args]
    return pd.DataFrame(out)


def _revet_star(a):
    return revet_one(*a)
