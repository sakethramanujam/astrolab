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

ODD_EVEN = "odd_even"
# Distinct flag when fewer than 2 odd or fewer than 2 even transits remain after the odd-even transit mask
# (n_odd < 2 or n_even < 2). A row with `odd_even_too_few` does NOT pass light-curve vetting checks, and is
# NOT counted as an `odd_even` depth-comparison failure.
ODD_EVEN_TOO_FEW = "odd_even_too_few"

_EPS = 1e-9

DEFAULT_ODD_EVEN_MASK = dict(
    enabled=False,
    min_coverage=0.75,
    min_gap_dump_dist_d=0.375,
    min_gap_dump_dist_h=9.0,
    gap_d=0.5,
    oot_baseline_t14=1.0,
    oot_one_side=True,
)


def odd_even_mask_cfg(cfg=None, odd_even_mask=None):
    """Resolve odd-even transit mask settings from `cfg` (and optional override `odd_even_mask`).

    Configs that omit the mask keys behave as mask-off (`enabled=False`), reproducing v1.1 behaviour.
    Distance threshold `Y` defaults to 9 hours = 0.375 days (half the 0.75 d biweight detrend window).
    """
    out = dict(DEFAULT_ODD_EVEN_MASK)
    v = (cfg or {}).get("vet") if isinstance(cfg, dict) else None
    sources = []
    if isinstance(cfg, dict):
        if "odd_even_mask" in cfg:
            sources.append(cfg["odd_even_mask"])
        if isinstance(v, dict):
            for k_flat, k_dst in (
                ("odd_even_mask_enabled", "enabled"),
                ("odd_even_min_coverage", "min_coverage"),
                ("odd_even_min_gap_dump_dist_d", "min_gap_dump_dist_d"),
                ("odd_even_min_gap_dump_dist_h", "min_gap_dump_dist_h"),
                ("odd_even_gap_d", "gap_d"),
                ("odd_even_oot_baseline_t14", "oot_baseline_t14"),
                ("odd_even_oot_one_side", "oot_one_side"),
            ):
                if k_flat in v and v[k_flat] is not None:
                    out[k_dst] = v[k_flat]
            if "odd_even_mask" in v:
                sources.append(v["odd_even_mask"])
        scr_gap = (cfg.get("scramble") or {}).get("gap_d") if isinstance(cfg.get("scramble"), dict) else None
        if scr_gap is not None and "gap_d" not in ((v or {}).get("odd_even_mask") or {}):
            out["gap_d"] = float(scr_gap)
    if odd_even_mask is not None:
        sources.append(odd_even_mask)

    for src in sources:
        if isinstance(src, bool):
            out["enabled"] = bool(src)
        elif isinstance(src, dict):
            d_set = "min_gap_dump_dist_d" in src or "Y_d" in src or "y_d" in src
            h_set = "min_gap_dump_dist_h" in src or "Y_h" in src or "y_h" in src
            for k_alias, k_dst in (
                ("enabled", "enabled"),
                ("mask_enabled", "enabled"),
                ("min_coverage", "min_coverage"),
                ("X", "min_coverage"),
                ("x", "min_coverage"),
                ("min_gap_dump_dist_d", "min_gap_dump_dist_d"),
                ("Y_d", "min_gap_dump_dist_d"),
                ("y_d", "min_gap_dump_dist_d"),
                ("min_gap_dump_dist_h", "min_gap_dump_dist_h"),
                ("Y_h", "min_gap_dump_dist_h"),
                ("y_h", "min_gap_dump_dist_h"),
                ("gap_d", "gap_d"),
                ("oot_baseline_t14", "oot_baseline_t14"),
                ("oot_t14", "oot_baseline_t14"),
                ("oot_one_side", "oot_one_side"),
                ("one_side", "oot_one_side"),
            ):
                if k_alias in src and src[k_alias] is not None:
                    out[k_dst] = src[k_alias]
            if h_set and not d_set:
                out["min_gap_dump_dist_d"] = float(out["min_gap_dump_dist_h"]) / 24.0
            elif d_set and not h_set:
                out["min_gap_dump_dist_h"] = float(out["min_gap_dump_dist_d"]) * 24.0

    out["enabled"] = bool(out["enabled"])
    out["min_coverage"] = float(out["min_coverage"])
    out["min_gap_dump_dist_d"] = float(out["min_gap_dump_dist_d"])
    out["min_gap_dump_dist_h"] = float(out["min_gap_dump_dist_d"]) * 24.0
    out["gap_d"] = float(out["gap_d"])
    out["oot_baseline_t14"] = float(out["oot_baseline_t14"])
    out["oot_one_side"] = bool(out["oot_one_side"])
    return out


def odd_even_mask_enabled(cfg=None, odd_even_mask=None) -> bool:
    """True when the v1.2 odd-even transit mask is enabled in `cfg` / `odd_even_mask`."""
    return bool(odd_even_mask_cfg(cfg, odd_even_mask)["enabled"])


def dump_times(t_all, quality, bit=32, gap=0.05):
    """Centres of contiguous runs of cadences flagged with the momentum-dump bit."""
    t_all = np.asarray(t_all, float)
    m = ((np.asarray(quality, int) & bit) > 0) & np.isfinite(t_all)
    td = np.sort(t_all[m])
    if td.size == 0:
        return np.array([])
    groups = np.split(td, np.flatnonzero(np.diff(td) > gap) + 1)
    return np.array([float(g.mean()) for g in groups])


def gap_edges(t, gap_d: float = 0.5) -> np.ndarray:
    """Gap edge timestamps in a light-curve time array `t`.

    Data gap definition: reuses the existing light-curve gap threshold `gap_d = 0.5` d from
    `exosearch.lightcurve.gap_edge_mask` (and `configs/v1.1.json` `scramble.gap_d = 0.5`). Any interval
    between consecutive valid cadences `t[i]` and `t[i + 1]` with `t[i + 1] - t[i] > gap_d` is a data gap,
    and both `t[i]` (last cadence before the gap) and `t[i + 1]` (first cadence after the gap) are gap edges.
    """
    t = np.asarray(t, dtype=float)
    t = np.sort(t[np.isfinite(t)])
    if t.size < 2:
        return np.array([], dtype=float)
    idx = np.flatnonzero(np.diff(t) > float(gap_d))
    if idx.size == 0:
        return np.array([], dtype=float)
    return np.sort(np.unique(np.concatenate([t[idx], t[idx + 1]])))


def odd_even_transit_mask(
    t,
    P: float,
    t0: float,
    duration: float,
    *,
    cfg=None,
    dumps=None,
    gap_edges_t=None,
    odd_even_mask=None,
    min_coverage: float | None = None,
    min_gap_dump_dist_d: float | None = None,
    gap_d: float | None = None,
    oot_baseline_t14: float | None = None,
    oot_one_side: bool | None = None,
):
    """Evaluate the fixed up-front v1.2 odd-even transit mask for each transit epoch in `t`.

    The mask is decided purely from the cadence timestamps `t`, candidate `(P, t0, duration)` (where
    `duration` is the candidate transit duration T14), data gaps, and momentum dumps — never from the
    observed transit depths.

    A transit with midpoint `t_mid = t0 + k * P`, ingress `t_ing = t_mid - 0.5 * T14`, and egress
    `t_egr = t_mid + 0.5 * T14` is kept iff all three conditions hold:
      1. In-transit coverage >= `X` (default `0.75`):
         Measured on the same time series `t` used for the odd-even depths. With median cadence step
         `dt = median(diff(t))`, the expected number of cadences in the in-transit window
         `[t_mid - 0.5 * T14, t_mid + 0.5 * T14)` is `n_expected = T14 / dt`, and
         `coverage = min(1.0, n_in / (T14 / dt))`, where `n_in` is the number of cadences in
         `[t_mid - 0.5 * T14, t_mid + 0.5 * T14)`. Kept when `coverage >= X`.
      2. Distance from data gaps and momentum dumps >= `Y` (default `9` hours = `0.375` days):
         `min(|t_mid - t_gap_edge|) >= Y` for all gap edges `t_gap_edge` (last cadence before and first
         cadence after any step `diff(t) > gap_d`, default `gap_d = 0.5` d matching
         `lightcurve.gap_edge_mask`), and `min(|t_mid - t_dump|) >= Y` for all momentum-dump times `t_dump`
         (`vet11.dump_times` / QUALITY bit 32).
      3. Out-of-transit (OOT) baseline >= `B * T14` (default `B = 1.0`):
         Immediately before ingress `[t_ing - B * T14, t_ing)` and/or immediately after egress
         `[t_egr, t_egr + B * T14)` (at least one side when `oot_one_side=True`, the default; both sides
         when `oot_one_side=False`), unbroken by data gaps (`> gap_d`) and with full grid coverage over
         `B * T14`.

    Epoch parity follows `astropy.timeseries.BoxLeastSquares.compute_stats`: epoch index
    `k = round((t_mid - t0) / P)` is odd when `k % 2 != 0` and even when `k % 2 == 0`.
    """
    mcfg = odd_even_mask_cfg(cfg, odd_even_mask)
    if min_coverage is not None:
        mcfg["min_coverage"] = float(min_coverage)
    if min_gap_dump_dist_d is not None:
        mcfg["min_gap_dump_dist_d"] = float(min_gap_dump_dist_d)
        mcfg["min_gap_dump_dist_h"] = float(min_gap_dump_dist_d) * 24.0
    if gap_d is not None:
        mcfg["gap_d"] = float(gap_d)
    if oot_baseline_t14 is not None:
        mcfg["oot_baseline_t14"] = float(oot_baseline_t14)
    if oot_one_side is not None:
        mcfg["oot_one_side"] = bool(oot_one_side)

    X = mcfg["min_coverage"]
    Y = mcfg["min_gap_dump_dist_d"]
    gd = mcfg["gap_d"]
    B = mcfg["oot_baseline_t14"]
    one_side = mcfg["oot_one_side"]

    t = np.asarray(t, dtype=float)
    t = np.sort(t[np.isfinite(t)])
    P = float(P)
    t0 = float(t0)
    D = float(duration)

    if t.size < 2 or not (P > 0 and D > 0 and np.isfinite(P) and np.isfinite(t0) and np.isfinite(D)):
        empty_i = np.array([], dtype=np.int64)
        empty_f = np.array([], dtype=float)
        empty_b = np.array([], dtype=bool)
        return dict(
            epochs=empty_i, t_mid=empty_f, is_odd=empty_b, is_even=empty_b,
            n_in=empty_i, n_expected=0.0, dt=np.nan, coverage=empty_f, pass_coverage=empty_b,
            gap_dist=empty_f, pass_gap=empty_b, dump_dist=empty_f, pass_dump=empty_b,
            oot_before_t14=empty_f, oot_after_t14=empty_f,
            pass_oot_before=empty_b, pass_oot_after=empty_b, pass_oot=empty_b,
            keep=empty_b, kept_epochs=empty_i, kept_odd_epochs=empty_i, kept_even_epochs=empty_i,
            n_odd=0, n_even=0, too_few=True,
            gap_edges=np.array([], dtype=float), dumps=np.array([], dtype=float), params=mcfg,
        )

    dt = float(np.median(np.diff(t)))
    n_exp = float(D / dt) if dt > 0 else 0.0
    ge = gap_edges(t, gap_d=gd) if gap_edges_t is None else np.sort(np.asarray(gap_edges_t, dtype=float))
    ge = ge[np.isfinite(ge)]
    dm = np.asarray(dumps if dumps is not None else [], dtype=float)
    dm = np.sort(dm[np.isfinite(dm)])

    k_min = int(np.ceil((t[0] - 0.5 * D - t0) / P - _EPS))
    k_max = int(np.floor((t[-1] + 0.5 * D - t0) / P + _EPS))
    epochs = np.arange(k_min, k_max + 1, dtype=np.int64) if k_max >= k_min else np.array([], dtype=np.int64)
    n_ep = epochs.size

    t_mid_arr = t0 + epochs.astype(float) * P
    is_odd = (epochs % 2) != 0
    is_even = ~is_odd

    n_in = np.zeros(n_ep, dtype=np.int64)
    coverage = np.zeros(n_ep, dtype=float)
    pass_coverage = np.zeros(n_ep, dtype=bool)
    gap_dist = np.full(n_ep, np.inf, dtype=float)
    pass_gap = np.ones(n_ep, dtype=bool)
    dump_dist = np.full(n_ep, np.inf, dtype=float)
    pass_dump = np.ones(n_ep, dtype=bool)
    oot_before_t14 = np.zeros(n_ep, dtype=float)
    oot_after_t14 = np.zeros(n_ep, dtype=float)
    pass_oot_before = np.zeros(n_ep, dtype=bool)
    pass_oot_after = np.zeros(n_ep, dtype=bool)

    L_req = B * D
    n_req_oot = max(1, int(np.floor(L_req / dt + 1e-6))) if (L_req > 0 and dt > 0) else 0

    for i, t_mid in enumerate(t_mid_arr):
        t_ing = t_mid - 0.5 * D
        t_egr = t_mid + 0.5 * D

        # 1. In-transit coverage on [t_mid - 0.5 * D, t_mid + 0.5 * D)
        m_in_i = (t >= t_ing - _EPS) & (t < t_egr - _EPS)
        cnt = int(np.sum(m_in_i))
        cov = float(min(1.0, cnt / n_exp)) if n_exp > 0 else 0.0
        n_in[i] = cnt
        coverage[i] = cov
        pass_coverage[i] = bool(cov >= X - _EPS)

        # 2. Distance from midpoint to nearest gap edge and momentum dump
        if ge.size > 0:
            gd_i = float(np.min(np.abs(ge - t_mid)))
            gap_dist[i] = gd_i
            pass_gap[i] = bool(gd_i >= Y - _EPS)
        if dm.size > 0:
            dd_i = float(np.min(np.abs(dm - t_mid)))
            dump_dist[i] = dd_i
            pass_dump[i] = bool(dd_i >= Y - _EPS)

        # 3. Out-of-transit baseline immediately before ingress and/or after egress
        if B <= 0:
            oot_before_t14[i] = np.inf
            oot_after_t14[i] = np.inf
            pass_oot_before[i] = True
            pass_oot_after[i] = True
            continue

        # Pre-ingress window [t_ing - L_req, t_ing)
        has_gap_pre = bool(ge.size > 0 and np.any((ge > t_ing - L_req + _EPS) & (ge < t_ing + _EPS)))
        t_pre_all = t[t < t_ing - _EPS]
        if t_pre_all.size > 0 and not has_gap_pre and (t_ing - t_pre_all[-1]) <= gd + _EPS:
            # Restrict to contiguous segment after the last gap <= t_ing
            if ge.size > 0:
                ge_before = ge[ge <= t_ing + _EPS]
                if ge_before.size > 0:
                    t_pre_all = t_pre_all[t_pre_all >= ge_before[-1] - _EPS]
            t_pre_win = t_pre_all[t_pre_all >= t_ing - L_req - _EPS]
            n_pre = t_pre_win.size
            span_pre = float(t_ing - t_pre_all[0]) if t_pre_all.size > 0 else 0.0
            if n_pre >= n_req_oot and (span_pre >= L_req - _EPS or (n_pre * dt) >= L_req - _EPS):
                oot_before_t14[i] = float(max(n_pre * dt, min(span_pre, L_req)) / D)
                pass_oot_before[i] = True
            else:
                oot_before_t14[i] = float((n_pre * dt) / D)
                pass_oot_before[i] = False

        # Post-egress window [t_egr, t_egr + L_req) (or (t_egr, t_egr + L_req] when grid-aligned)
        has_gap_post = bool(ge.size > 0 and np.any((ge > t_egr - _EPS) & (ge < t_egr + L_req - _EPS)))
        t_post_all = t[t >= t_egr - _EPS]
        if t_post_all.size > 0 and not has_gap_post and (t_post_all[0] - t_egr) <= gd + _EPS:
            if ge.size > 0:
                ge_after = ge[ge >= t_egr - _EPS]
                if ge_after.size > 0:
                    t_post_all = t_post_all[t_post_all <= ge_after[0] + _EPS]
            n_post_lo = int(np.sum((t_post_all >= t_egr - _EPS) & (t_post_all < t_egr + L_req - _EPS)))
            n_post_hi = int(np.sum((t_post_all > t_egr + _EPS) & (t_post_all <= t_egr + L_req + _EPS)))
            n_post = max(n_post_lo, n_post_hi)
            span_post = float(t_post_all[-1] - t_egr) if t_post_all.size > 0 else 0.0
            if (
                n_post >= n_req_oot
                and (t_post_all[-1] > t_egr + _EPS or n_req_oot == 1)
                and (span_post >= L_req - _EPS or (n_post * dt) >= L_req - _EPS)
            ):
                oot_after_t14[i] = float(max(n_post * dt, min(span_post, L_req)) / D)
                pass_oot_after[i] = True
            else:
                oot_after_t14[i] = float(
                    (n_post_hi if t_post_all[-1] <= t_egr + _EPS else n_post) * dt / D
                )
                pass_oot_after[i] = False

    pass_oot = (pass_oot_before | pass_oot_after) if one_side else (pass_oot_before & pass_oot_after)
    keep = pass_coverage & pass_gap & pass_dump & pass_oot
    kept_epochs = epochs[keep]
    kept_odd_epochs = epochs[keep & is_odd]
    kept_even_epochs = epochs[keep & is_even]
    n_odd = int(kept_odd_epochs.size)
    n_even = int(kept_even_epochs.size)

    return dict(
        epochs=epochs,
        t_mid=t_mid_arr,
        is_odd=is_odd,
        is_even=is_even,
        n_in=n_in,
        n_expected=n_exp,
        dt=dt,
        coverage=coverage,
        pass_coverage=pass_coverage,
        gap_dist=gap_dist,
        pass_gap=pass_gap,
        dump_dist=dump_dist,
        pass_dump=pass_dump,
        oot_before_t14=oot_before_t14,
        oot_after_t14=oot_after_t14,
        pass_oot_before=pass_oot_before,
        pass_oot_after=pass_oot_after,
        pass_oot=pass_oot,
        keep=keep,
        kept_epochs=kept_epochs,
        kept_odd_epochs=kept_odd_epochs,
        kept_even_epochs=kept_even_epochs,
        n_odd=n_odd,
        n_even=n_even,
        too_few=bool(n_odd < 2 or n_even < 2),
        gap_edges=ge,
        dumps=dm,
        params=mcfg,
    )


def _weighted_box_depth(m_in_sub, y_out: float, var_out: float, f: np.ndarray, w: np.ndarray):
    if np.any(m_in_sub) and np.isfinite(var_out):
        w_sum = float(np.sum(w[m_in_sub]))
        if w_sum > 0:
            var_m = 1.0 / w_sum
            y_m = float(np.sum(f[m_in_sub] * w[m_in_sub])) * var_m
            return float(y_out - y_m), float(np.sqrt(var_m + var_out))
    return np.nan, np.inf


def masked_odd_even(
    t,
    f,
    e,
    P: float,
    t0: float,
    duration: float,
    *,
    cfg=None,
    dumps=None,
    gap_edges_t=None,
    odd_even_mask=None,
    **kw,
):
    """Compute odd/even transit depths and verdict on the transits that pass `odd_even_transit_mask`.

    Only transits passing the mask are included in `depth_odd` and `depth_even`. The out-of-transit
    baseline mean `y_out` and variance `var_out` use all cadences outside any transit window (`~m_in`),
    matching `astropy.timeseries.BoxLeastSquares.compute_stats`. Rejected transits are excluded from
    the odd/even depths only; their in-transit cadences are still excluded from `y_out`.

    Too-few rule:
      If fewer than 2 odd or fewer than 2 even transits remain after the mask (`n_odd < 2` or `n_even < 2`),
      `flag` is set to `ODD_EVEN_TOO_FEW` (`"odd_even_too_few"`), `odd_even_fail` is `False` (it is NOT
      counted as an `odd_even` depth-comparison failure), and `passed` is `False` (the row does not pass
      light-curve checks).
    Otherwise (`n_odd >= 2` and `n_even >= 2`):
      Runs the standard v1.1 3-sigma and 10% test on `(depth_odd, depth_even)` of the masked transit set.
    """
    t = np.asarray(t, dtype=float)
    f = np.asarray(f, dtype=float)
    e = np.asarray(e, dtype=float)
    valid = np.isfinite(t) & np.isfinite(f) & np.isfinite(e) & (e > 0)
    order = np.argsort(t[valid])
    t, f, e = t[valid][order], f[valid][order], e[valid][order]

    v = (cfg or {}).get("vet") if isinstance(cfg, dict) else None
    thresh_sigma = float((v or {}).get("odd_even_sigma", vet10.THRESH["odd_even_sigma"]))
    thresh_rel = float((v or {}).get("odd_even_rel", 0.1))

    mask = odd_even_transit_mask(
        t, P, t0, duration, cfg=cfg, dumps=dumps, gap_edges_t=gap_edges_t, odd_even_mask=odd_even_mask, **kw
    )

    P = float(P)
    t0 = float(t0)
    D = float(duration)
    w = 1.0 / (e ** 2)
    ep = np.round((t - t0) / P).astype(np.int64)
    m_in_win = (t >= (t0 + ep * P - 0.5 * D - _EPS)) & (t < (t0 + ep * P + 0.5 * D - _EPS))
    m_out = ~m_in_win

    if np.any(m_out):
        w_out = float(np.sum(w[m_out]))
        var_out = 1.0 / w_out if w_out > 0 else np.inf
        y_out = float(np.sum(f[m_out] * w[m_out])) * var_out if w_out > 0 else np.nan
    else:
        var_out, y_out = np.inf, np.nan

    def _epoch_cadences(kept_eps):
        if kept_eps.size == 0:
            return np.zeros(t.size, dtype=bool)
        return m_in_win & np.isin(ep, kept_eps)

    m_odd = _epoch_cadences(mask["kept_odd_epochs"])
    m_even = _epoch_cadences(mask["kept_even_epochs"])

    d_odd, e_odd = _weighted_box_depth(m_odd, y_out, var_out, f, w)
    d_even, e_even = _weighted_box_depth(m_even, y_out, var_out, f, w)

    eo = float(np.hypot(e_odd, e_even)) if (np.isfinite(e_odd) and np.isfinite(e_even)) else np.inf
    if eo > 0 and np.isfinite(eo) and np.isfinite(d_odd) and np.isfinite(d_even):
        oe_sigma = float(abs(d_odd - d_even) / eo)
    else:
        oe_sigma = 0.0
    mean_oe = 0.5 * (d_odd + d_even) if (np.isfinite(d_odd) and np.isfinite(d_even)) else np.nan
    if np.isfinite(mean_oe) and mean_oe != 0:
        oe_rel = float(abs(d_odd - d_even) / abs(mean_oe))
    else:
        oe_rel = np.inf

    # Optional LOTO diagnostic (logged only; never used as a pass/fail rule)
    loto_sigmas = []
    if not mask["too_few"]:
        for drop_k in mask["kept_epochs"]:
            ko = mask["kept_odd_epochs"][mask["kept_odd_epochs"] != drop_k]
            ke = mask["kept_even_epochs"][mask["kept_even_epochs"] != drop_k]
            if ko.size == 0 or ke.size == 0:
                continue
            do_k, eo_k = _weighted_box_depth(_epoch_cadences(ko), y_out, var_out, f, w)
            de_k, ee_k = _weighted_box_depth(_epoch_cadences(ke), y_out, var_out, f, w)
            err_k = float(np.hypot(eo_k, ee_k))
            if err_k > 0 and np.isfinite(err_k) and np.isfinite(do_k) and np.isfinite(de_k):
                loto_sigmas.append(float(abs(do_k - de_k) / err_k))
    loto_min_sigma = float(min(loto_sigmas)) if loto_sigmas else np.nan
    loto_max_sigma = float(max(loto_sigmas)) if loto_sigmas else np.nan

    if mask["too_few"]:
        flag = ODD_EVEN_TOO_FEW
        odd_even_fail = False
        too_few = True
        passed = False
    else:
        too_few = False
        odd_even_fail = bool(oe_sigma > thresh_sigma and oe_rel > thresh_rel)
        flag = ODD_EVEN if odd_even_fail else None
        passed = not odd_even_fail

    return dict(
        depth_odd=float(d_odd),
        depth_odd_err=float(e_odd),
        depth_even=float(d_even),
        depth_even_err=float(e_even),
        odd_even_sigma=float(oe_sigma),
        odd_even_rel=float(oe_rel),
        n_odd=mask["n_odd"],
        n_even=mask["n_even"],
        too_few=too_few,
        odd_even_fail=odd_even_fail,
        flag=flag,
        passed=passed,
        loto_min_sigma=loto_min_sigma,
        loto_max_sigma=loto_max_sigma,
        mask=mask,
    )


def screens(st, best, meta, cfg, lc=None, dumps=None):
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
    if odd_even_mask_enabled(cfg):
        if lc is None:
            raise ValueError("odd_even_mask is enabled but no light curve (lc) was supplied to screens()")
        tb, fb, eb = lc[0], lc[1], lc[2]
        oe_res = masked_odd_even(tb, fb, eb, P, best.get("t0", 0.0), dur, cfg=cfg, dumps=dumps)
        d_odd, d_even = oe_res["depth_odd"], oe_res["depth_even"]
        oe, oe_rel = oe_res["odd_even_sigma"], oe_res["odd_even_rel"]
        if oe_res["flag"] is not None:
            flags.append(oe_res["flag"])
    else:
        if oe > v["odd_even_sigma"] and oe_rel > v["odd_even_rel"]:
            flags.append("odd_even")
    if sec_sig > v["secondary_sigma"] and sec_rel > v["secondary_rel"]:
        flags.append("secondary")
    td = vet10.too_deep(depth, cfg)
    if td["mode"] == "cut" and td["flagged"]:
        flags.append(vet10.TOO_DEEP)
    elif td["mode"] == "flag" and td["flagged"]:
        warns.append(vet10.TOO_DEEP)          # recorded, not a screen rejection (vet_final is authoritative)
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
