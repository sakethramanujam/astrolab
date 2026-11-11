"""v1.1 search on one prepared light curve: capped BLS scan, beta treatment, plateau tie-break, exact-depth S.

Per trial duration D (the v1.0.0 set, 1-6 h) and period P on the fixed grid, inside the duration cap:
  scan:   astropy BLS 'snr' objective (fast method) with sigma_D = CDPP_D sqrt(D/bin) on the binned flux, i.e. the
          box contrast over its own in-and-out error, times the beta treatment factor beta_D / beta'_D.
  peak:   the registered plateau tie-break on the scan (tie_break.md amendment 1/1a) gives the first candidate; up to
          refine_top_k distinct local peaks (|dP/P| > refine_min_sep_rel) are refined, and the star's peak is the
          candidate with the largest S (ties within TIE_REL: lowest period).
  S:      depth_exact * sqrt(N_tr) / CDPP_D * sqrt(1 - D/P) * beta_D / beta'_D
          depth_exact = weighted (1/err^2) out-of-transit minus in-transit mean on the binned flux at
          (P_grid, t0, D) (the compute_stats box depth), with in-transit points of momentum-dump epochs removed;
          N_tr = coverage-weighted transit count over the remaining epochs; sqrt(1 - D/P) restores the
          out-of-transit term (white-noise calibration, docs/v1.1_odd_even_secondary.md).
Every beta option (a-d) is searched (its own peak) and stored for sensitivity rows; the configured option is the
statistic. S_v10 (v1.0.0 definition: BLS fast-method depth, all epochs, no out-of-transit term, no beta treatment)
and S_exact (exact depth, no out-of-transit term, no beta treatment) are kept at the chosen peak.
"""
from __future__ import annotations

import json

import numpy as np
from astropy.timeseries import BoxLeastSquares

from .. import bls as bls10
from . import grid

DURATIONS = bls10.DURATIONS            # 0.04 .. 0.25 d, unchanged from v1.0.0
TIE_REL = 1e-9                          # registered tie band


# ---------------------------------------------------------------- beta treatment (beta_treatment.md)
def sigma_beta(n_cad, cadence_d, D):
    nw = max(float(n_cad) * float(cadence_d) / float(D), 1.0)
    return 1.0 / np.sqrt(2.0 * nw)


def beta_prime(option, beta, sb, prior_D=None):
    if option == "a":
        return beta
    if option == "b":
        return max(1.0, beta)
    if option == "c":
        return max(beta, 1.0 - 2.0 * sb)
    if option == "d":
        if prior_D is None:
            raise ValueError("option d needs the registered beta prior (configs/v1.1_beta_prior.json)")
        m, t2 = prior_D
        return (beta / sb ** 2 + m / t2) / (1.0 / sb ** 2 + 1.0 / t2)
    raise ValueError(f"unknown beta option {option!r}")


# ---------------------------------------------------------------- exact box depth with momentum-dump epochs removed
def in_transit(t, P, t0, D):
    return np.abs(((t - t0 + 0.5 * P) % P) - 0.5 * P) < 0.5 * D


def dump_epochs(P, t0, epochs, dumps, window):
    if dumps is None or len(dumps) == 0 or len(epochs) == 0:
        return np.zeros(len(epochs), bool)
    centres = t0 + np.asarray(epochs, float) * P
    return np.min(np.abs(centres[:, None] - np.asarray(dumps, float)[None, :]), axis=1) < window


def exact_box(tb, fb, eb, P, t0, D, dumps=None, window=0.1):
    """-> dict(depth, n_tr, n_epochs, n_dump_epochs, n_tr_all). Weighted in/out means (compute_stats depth)."""
    w = 1.0 / np.asarray(eb, float) ** 2
    inn = in_transit(tb, P, t0, D)
    ep = np.round((tb - t0) / P).astype(np.int64)
    ue = np.unique(ep[inn])
    bad = dump_epochs(P, t0, ue, dumps, window)
    keep_in = inn & ~np.isin(ep, ue[bad])
    out = ~inn
    bin_d = float(np.median(np.diff(tb)))
    n_all = bls10.transit_count(tb, P, t0, D, bin_d)
    if keep_in.sum() == 0 or out.sum() == 0:
        return dict(depth=np.nan, n_tr=0.0, n_epochs=int(len(ue)), n_dump_epochs=int(bad.sum()), n_tr_all=n_all)
    m_in = float(np.sum(w[keep_in] * fb[keep_in]) / np.sum(w[keep_in]))
    m_out = float(np.sum(w[out] * fb[out]) / np.sum(w[out]))
    tk = tb[keep_in | out]
    n_tr = bls10.transit_count(tk, P, t0, D, bin_d)
    return dict(depth=m_out - m_in, n_tr=n_tr, n_epochs=int(len(ue)), n_dump_epochs=int(bad.sum()), n_tr_all=n_all)


# ---------------------------------------------------------------- plateaus (tie_break.md amendment 1/1a)
def _step(periods, i):
    n = len(periods)
    return float(periods[min(i + 1, n - 1)] - periods[max(i - 1, 0)]) / (2 if 0 < i < n - 1 else 1)


def _plateau(periods, lo, hi):
    c = lo + (hi - lo) // 2
    width_days = float(periods[hi] - periods[lo]) + 0.5 * (_step(periods, lo) + _step(periods, hi))
    return dict(lo=lo, hi=hi, centre=c, width=int(hi - lo + 1), P_lo=float(periods[lo]), P_hi=float(periods[hi]),
                hw=width_days / 2)


def plateaus_in_band(per_d, periods):
    """Registered rule on the global maximum: in-band grid points over all durations, split into contiguous plateaus
    and disjoint peaks; returns (plateaus sorted by period, {index: [durations in band]})."""
    top = max(float(np.max(x["p"])) for x in per_d)
    band = top - TIE_REL * abs(top)
    inb = {}
    for x in per_d:
        for i in np.flatnonzero(x["p"] >= band):
            inb.setdefault(int(i), []).append(x["D"])
    idx = np.array(sorted(inb))
    runs = np.split(idx, np.flatnonzero(np.diff(idx) > 1) + 1)
    return [_plateau(periods, int(r[0]), int(r[-1])) for r in runs], inb


def local_plateau(smax, periods, i):
    """Contiguous run around index i where the max-over-D scan is within TIE_REL of smax[i]."""
    v = smax[i]
    lo = hi = i
    while lo > 0 and smax[lo - 1] >= v - TIE_REL * abs(v):
        lo -= 1
    while hi < len(smax) - 1 and smax[hi + 1] >= v - TIE_REL * abs(v):
        hi += 1
    return _plateau(periods, lo, hi)


def candidates(per_d, periods, k, sep):
    """First: the registered tie-break choice. Then distinct local maxima of the max-over-D scan, descending."""
    plats, inb = plateaus_in_band(per_d, periods)
    ch = plats[0]
    c = ch["centre"]
    Ds = sorted(inb.get(c, [])) or sorted({d for i in range(ch["lo"], ch["hi"] + 1) for d in inb[i]})
    first = dict(plat=ch, idx=c, D=Ds[0], global_peak=True)
    out = [first]
    P = np.stack([x["p"] for x in per_d])
    smax = P.max(axis=0)
    if k <= 1:
        return out, plats, inb, smax
    loc = np.flatnonzero(np.isfinite(smax) & (smax >= np.r_[-np.inf, smax[:-1]]) & (smax >= np.r_[smax[1:], -np.inf]))
    loc = loc[np.argsort(-smax[loc], kind="stable")]
    chosen_P = [0.5 * (ch["P_lo"] + ch["P_hi"])]
    for i in loc:
        if len(out) >= k:
            break
        if any(abs(periods[i] / p - 1) <= sep for p in chosen_P):
            continue
        pl = local_plateau(smax, periods, int(i))
        j = pl["centre"]
        col = P[:, j]
        jd = int(np.flatnonzero(col >= col.max() - TIE_REL * abs(col.max()))[0])    # shortest D in band
        out.append(dict(plat=pl, idx=j, D=per_d[jd]["D"], global_peak=False))
        chosen_P.append(0.5 * (pl["P_lo"] + pl["P_hi"]))
    return out, plats, inb, smax


# ---------------------------------------------------------------- the search
def search(tb, fb, eb, tflat, fflat, rho=None, dumps=None, cfg=None, options=("a", "b", "c", "d"), primary=None,
           prior=None, cap=True):
    """-> (primary result dict, {option: summary}). prior: {"0.04": [median, tau2], ...} for option d."""
    cfg = cfg or {}
    k = int(cfg.get("refine_top_k", 5))
    sep = float(cfg.get("refine_min_sep_rel", 0.01))
    capc = cfg.get("duration_cap", {"margin": 3.0, "q_max_no_density": 0.1})
    mdc = cfg.get("momentum_dump", {"window_d": 0.1})
    periods = grid.period_grid(tb, DURATIONS)
    sw, cad, prof = bls10.noise_profile(tflat, fflat, DURATIONS)
    n_cad = int(np.asarray(tflat).size)
    bin_d = float(np.median(np.diff(tb)))
    scans = []
    for D in DURATIONS:
        D = float(D)
        cd, beta = prof[D]
        if not (cd == cd and cd > 0 and beta == beta and beta > 0):
            continue
        al = grid.allowed(periods, D, rho, capc["margin"], capc["q_max_no_density"]) if cap else (D < periods)
        if not al.any():
            continue
        dy = np.full(tb.size, cd * np.sqrt(max(D / bin_d, 1.0)))
        r = BoxLeastSquares(tb, fb, dy=dy).power(periods[al], [D], objective="snr", method="fast")
        p = np.full(periods.size, -np.inf)
        t0 = np.full(periods.size, np.nan)
        dep = np.full(periods.size, np.nan)
        p[al] = np.nan_to_num(np.asarray(r.power, float), nan=-np.inf)
        t0[al] = np.asarray(r.transit_time, float)
        dep[al] = np.asarray(r.depth, float)
        sb = sigma_beta(n_cad, cad, D)
        scans.append(dict(D=D, cd=cd, beta=beta, sb=sb, p=p, t0=t0, depth_box=dep))
    if not scans:
        raise ValueError("no valid noise estimate inside the duration cap")
    res = {}
    for o in options:
        fac = {}
        for x in scans:
            pD = None if prior is None else prior.get(f"{x['D']:.2f}")
            fac[x["D"]] = x["beta"] / beta_prime(o, x["beta"], x["sb"], pD) if (o != "d" or pD is not None) else np.nan
        if o == "d" and any(not np.isfinite(v) for v in fac.values()):
            res[o] = None
            continue
        per_d = [dict(x, p=x["p"] * fac[x["D"]]) for x in scans]
        cands, plats, inb, smax = candidates(per_d, periods, k, sep)
        best = None
        for cnd in cands:
            x = next(y for y in per_d if y["D"] == cnd["D"])
            i = cnd["idx"]
            Pg, t0, D = float(periods[i]), float(x["t0"][i]), cnd["D"]
            eb_ = exact_box(tb, fb, eb, Pg, t0, D, dumps, mdc.get("window_d", 0.1))
            s_ex = eb_["depth"] * np.sqrt(eb_["n_tr"]) / x["cd"] if eb_["n_tr"] > 0 and eb_["depth"] == eb_["depth"] else 0.0
            s_oot = s_ex * np.sqrt(max(1.0 - D / Pg, 0.0))
            S = float(s_oot * fac[D])
            pl = cnd["plat"]
            cand = dict(S=S, S_exact=float(s_ex), S_oot=float(s_oot), period=0.5 * (pl["P_lo"] + pl["P_hi"]), P_grid=Pg,
                        t0=t0, duration=D, depth_exact=float(eb_["depth"]), depth_box=float(x["depth_box"][i]),
                        S_scan=float(x["p"][i]), cdpp_d=float(x["cd"]), beta_d=float(x["beta"]),
                        beta_prime=float(x["beta"] / fac[D]), sigma_beta=float(x["sb"]), n_tr=float(eb_["n_tr"]),
                        n_tr_all=float(eb_["n_tr_all"]), n_epochs=eb_["n_epochs"], n_dump_epochs=eb_["n_dump_epochs"],
                        P_grid_halfwidth=pl["hw"], plateau_width=pl["width"], P_lowest_edge=pl["P_lo"],
                        refine_rank=cands.index(cnd), global_scan_peak=cnd["global_peak"])
            cand["S_v10"] = float(cand["depth_box"] * np.sqrt(cand["n_tr_all"]) / cand["cdpp_d"]) if cand["n_tr_all"] > 0 else 0.0
            if best is None or cand["S"] > best["S"] + TIE_REL * abs(best["S"]) or (
                    abs(cand["S"] - best["S"]) <= TIE_REL * abs(best["S"]) and cand["P_grid"] < best["P_grid"]):
                best = cand
        g = plats
        best.update(n_tie_band=int(sum(len(v) for v in inb.values())), n_plateaus=len(g), disjoint=len(g) > 1,
                    plateau_width_max=max(q["width"] for q in g),
                    plateaus=json.dumps([[q["P_lo"], q["P_hi"], q["width"], q["hw"]] for q in g]),
                    n_candidates=len(cands))
        res[o] = best
    prim = primary or options[0]
    if res.get(prim) is None:
        raise ValueError(f"primary beta option {prim} unavailable")
    out = dict(res[prim])
    smax_all = np.max(np.stack([x["p"] for x in scans]), axis=0)
    smax_all[~np.isfinite(smax_all)] = 0.0
    out.update(sigma_w=float(sw), cadence_d=float(cad), n_cad=n_cad, grid_n=int(periods.size),
               grid_df=float(1.0 / periods[0] - 1.0 / periods[1]), grid_sha=grid.grid_sha256(periods),
               n_durations_used=len(scans), sde=float(bls10.sde(smax_all)[0]),
               cdpp_by_d=json.dumps({f"{x['D']:.2f}": x["cd"] for x in scans}),
               beta_by_d=json.dumps({f"{x['D']:.2f}": x["beta"] for x in scans}))
    sens = {o: (None if v is None else {kk: v[kk] for kk in ("S", "period", "P_grid", "duration", "t0", "P_grid_halfwidth")})
            for o, v in res.items()}
    return out, sens
