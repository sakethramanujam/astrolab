"""Transit fit (batman, circular orbit, fixed period, Claret LD) -> transit-implied stellar density.

rho_star = 3 pi (a/R*)^3 / (G P^2)   (Seager & Mallen-Ornelas 2003; Rp << R*, e = 0)
Compared with the TIC v8.2 density as z = (log rho_tr - log rho_TIC) / sigma_log. Large |z| suggests
a blend (diluted, V-shaped -> low a/R*), a wrong host, or an eccentric orbit, so the flag is a reject
only beyond a threshold tuned on known planets vs false positives in the tuning stars.
"""
from __future__ import annotations

import numpy as np
from scipy.optimize import least_squares

from .stellar import G_CGS, RHO_SUN

DAY = 86400.0


def model(t, t0, per, k, ar, b, u):
    import batman
    p = batman.TransitParams()
    p.t0, p.per, p.rp, p.a = t0, per, k, ar
    p.inc = float(np.degrees(np.arccos(np.clip(b / ar, 0, 1))))
    p.ecc, p.w, p.limb_dark, p.u = 0.0, 90.0, "quadratic", list(u)
    return batman.TransitModel(p, t).light_curve(p)


def rho_from_ar(ar, per_d):
    return 3 * np.pi * ar ** 3 / (G_CGS * (per_d * DAY) ** 2) / RHO_SUN       # solar units


def ar_from_rho(rho_sun, per_d):
    return (rho_sun * RHO_SUN * G_CGS * (per_d * DAY) ** 2 / (3 * np.pi)) ** (1 / 3)


def fit(t, f, e, per, t0, dur, depth, u, rho_guess=1.0):
    """Least-squares fit of (t0, k, a/R*, b) on the detrended flux near transit. Returns dict or None."""
    ph = ((t - t0 + 0.5 * per) % per) - 0.5 * per
    w = np.abs(ph) < max(2.5 * dur, 0.1)
    if w.sum() < 30:
        return None
    tt, ff, ee = t[w], f[w], e[w]
    k0 = float(np.sqrt(max(depth, 1e-5)))
    ar0 = float(np.clip(ar_from_rho(rho_guess, per), 1.5, 200))
    x0 = np.array([t0, min(k0, 0.4), ar0, 0.3])
    lo, hi = [t0 - 0.5 * dur, 1e-3, 1.2, 0.0], [t0 + 0.5 * dur, 0.6, 400.0, 1.2]
    res_fn = lambda x: (model(tt, x[0], per, x[1], x[2], min(x[3], 1 + x[1]), u) - ff) / ee
    best = None
    for b0 in (0.1, 0.5, 0.85):                      # a/R*-b degeneracy: try several starts, keep the best
        x0[3] = b0
        try:
            r = least_squares(res_fn, x0, bounds=(lo, hi), x_scale="jac", max_nfev=400)
        except Exception:
            continue
        if best is None or r.cost < best.cost:
            best = r
    if best is None:
        return None
    J = best.jac
    dof = max(len(ff) - 4, 1)
    s2 = 2 * best.cost / dof
    try:
        cov = np.linalg.pinv(J.T @ J) * max(s2, 1.0)
        sig = np.sqrt(np.clip(np.diag(cov), 0, None))
    except Exception:
        sig = np.full(4, np.nan)
    t0f, k, ar, b = best.x
    rho = rho_from_ar(ar, per)
    return dict(fit_t0=float(t0f), fit_k=float(k), fit_ar=float(ar), fit_ar_err=float(sig[2]), fit_b=float(b),
                fit_chi2r=float(s2), rho_transit=float(rho),
                rho_transit_logerr=float(3 * sig[2] / ar / np.log(10)) if ar > 0 else np.nan)


def density_z(rho_tr, logerr_tr, rho_tic, err_tic):
    if not (rho_tr and rho_tic and rho_tr > 0 and rho_tic > 0):
        return None
    le_tic = err_tic / rho_tic / np.log(10) if err_tic else 0.13
    s = np.hypot(logerr_tr if np.isfinite(logerr_tr) else 0.3, le_tic)
    return float((np.log10(rho_tr) - np.log10(rho_tic)) / max(s, 0.05))
