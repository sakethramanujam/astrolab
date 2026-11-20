"""v1.2 catalogue cross-check, amended rules (exoplanet_search/docs/xcheck_catalogue_v1.2_amendment_2026-11-20.md).

Same scope, data and read discipline as xcheck.py (its pass-1 / pass-2 loaders are reused unchanged).
Changes: period tol = P_cat*T14/(2*T_span); t0 tol = max(T14/2, 3*sigma_pred); pre-score labels outside_grid and
outside_period_range; sibling -> found_sibling_instead (checked before aliases); aliases 1:2 2:1 1:3 3:1 2:3 3:2.
Outputs: *_v2.csv / xcheck_summary_v2.json next to the originals (originals untouched).
"""
from __future__ import annotations

import json
from fractions import Fraction
from pathlib import Path

import numpy as np
import pandas as pd

from exosearch import bls, config, lightcurve
import xcheck as x1          # reuse pass1_scan, pass2_load_science, epoch-error lookups, catalogue selection

D = Path(__file__).resolve().parent
DUR = np.asarray(bls.DURATIONS, float)                 # days
DMIN = float(DUR.min())
ALIASES = [Fraction(1, 2), Fraction(2, 1), Fraction(1, 3), Fraction(3, 1), Fraction(2, 3), Fraction(3, 2)]


def t_span(tic: int, sector: int) -> float:
    """Last minus first cadence passing lightcurve.good_mask (finite flux/err, QUALITY==0, err>0), days (BTJD)."""
    hits = sorted(config.lc_dir(sector).glob(f"*-{int(tic):016d}-*_lc.fits"))
    if len(hits) != 1:
        return np.nan
    t, f, e, q, _ = lightcurve.load_spoc_raw(hits[0])
    t = t[lightcurve.good_mask(f, e, q)]
    return float(t.max() - t.min()) if t.size >= 2 else np.nan


def ptol(p, t14, span):
    return p * t14 / (2.0 * span)


def main():
    parts = sorted((D.parent / "rows" / "parts").glob("*.parquet"))
    pass1, _ = x1.pass1_scan(parts)
    tics = set(int(t) for t in pass1.tic.unique())
    ps_df, ps_path = x1.kp.load_pscomppars()
    toi_df, toi_path = x1.kp.load_toi_table()
    assert x1.p_grid_halfwidth.load() is None          # sector-independent selection, as in xcheck.py
    tg = x1.kp.select_known_planet_targets(sector=59, searched_tics=tics, ps_df=ps_df, toi_df=toi_df, pgh_table=None)
    cp = tg[(tg.type == "ps") | ((tg.type == "toi") & tg.disp.isin(["CP", "KP"]))].reset_index(drop=True)
    inscope = pass1[pass1.tic.isin(set(cp.tic.astype(int)))].copy()
    rows, _ = x1.pass2_load_science(inscope)            # asserts: no S69, no holdout, real only, keys == in-scope
    ps_err, toi_n, toi_s = x1.build_epoch_error_lookups(ps_df, toi_df)

    spans = {(int(r.tic), int(r.sector)): t_span(r.tic, r.sector) for r in rows.itertuples()}

    pl = []
    for r in cp.to_dict("records"):
        s0, s0src = x1.lookup_sigma_t0(int(r["tic"]), x1._to_str(r["name"]), x1._to_str(r["t0_src"]), ps_err, toi_n, toi_s)
        pl.append(dict(tic=int(r["tic"]), name=x1._to_str(r["name"]), type=r["type"], disp=x1._to_str(r["disp"]),
                       P=float(r["period"]), sP=x1._to_float(r["period_err"]), T0=float(r["t0"]),   # days, BTJD
                       sT0=s0, sT0_src=s0src, T14=x1._to_float(r["dur_h"]) / 24.0,                    # hours -> days
                       dep=x1._to_float(r["depth_ppm"]) * 1e-6))                                     # ppm -> fraction
    by_tic = {}
    for p in pl:
        by_tic.setdefault(p["tic"], []).append(p)

    out = []
    for p in pl:
        for sr in rows[rows.tic == p["tic"]].sort_values("sector").to_dict("records"):
            span = spans[(p["tic"], int(sr["sector"]))]
            pmax = np.sqrt(DMIN / sr["grid_df"]) / 2.0 if sr["grid_df"] == sr["grid_df"] else np.nan  # grid.py: baseline/2
            o = dict(tic=p["tic"], name=p["name"], type=p["type"], disp=p["disp"], sector=int(sr["sector"]), key=sr["key"],
                     status=sr["status"], P_cat_d=p["P"], period_err_d=p["sP"], T0_cat_btjd=p["T0"], sigma_t0_d=p["sT0"],
                     sigma_t0_src=p["sT0_src"], T14_cat_d=p["T14"], depth_cat_frac=p["dep"], T_span_d=span,
                     grid_baseline_d=2 * pmax, pmax_d=pmax, P_grid_d=sr["P_grid"], t0_ours_btjd=sr["t0"],
                     duration_ours_d=sr["duration"], depth_ours_frac=sr["depth"], S=sr["S"],
                     period_pass=None, t0_pass=None, duration_pass=None, depth_pass=None, full_pass=False, scored=False,
                     alias_ratio="", sibling="", P_tol_d=np.nan, dP_d=np.nan, dt_h=np.nan, sigma_pred_h=np.nan,
                     t0_tol_h=np.nan, depth_ratio=np.nan, n_pass=0, reason_class="")
            complete = all(np.isfinite(v) and v > 0 for v in (p["P"], p["T14"], p["dep"])) and np.isfinite(p["T0"])
            if np.isfinite(p["T14"]) and p["T14"] > 0.25:
                o["reason_class"] = "outside_grid"
            elif np.isfinite(p["P"]) and (p["P"] < 0.5 or p["P"] > pmax):
                o["reason_class"] = "outside_period_range"
            elif not complete:
                o["reason_class"] = "catalogue_incomplete"
            elif sr["status"] != "ok":
                o.update(scored=True, reason_class="error_row", period_pass=False, t0_pass=False, duration_pass=False,
                         depth_pass=False)
            elif not np.isfinite(span) or span <= 0:
                o.update(reason_class="no_span")
            else:
                o["scored"] = True
                P, Pg, T14 = p["P"], float(sr["P_grid"]), p["T14"]
                tol = ptol(P, T14, span)
                dP = Pg - P
                per = abs(dP) <= tol
                n = int(round((sr["t0"] - p["T0"]) / P))
                dt = sr["t0"] - (p["T0"] + n * P)
                sig = np.hypot(p["sT0"] if np.isfinite(p["sT0"]) else 0.0, abs(n) * (p["sP"] if np.isfinite(p["sP"]) else 0.0))
                t0tol = max(T14 / 2.0, 3.0 * sig)
                t0ok = abs(dt) <= t0tol
                m = np.flatnonzero(np.isclose(DUR, sr["duration"]))
                inear = int(np.argmin(np.abs(DUR - T14)))
                durok = len(m) == 1 and abs(int(m[0]) - inear) <= 1
                ratio = sr["depth"] / p["dep"]
                depok = bool(np.isfinite(ratio) and 0.5 <= ratio <= 1.5)
                sib = ""
                ali = ""
                if not per:
                    for q in by_tic[p["tic"]]:
                        if q is not p and np.isfinite(q["P"]) and np.isfinite(q["T14"]) and q["T14"] > 0 and \
                                abs(Pg - q["P"]) <= ptol(q["P"], q["T14"], span):
                            sib = q["name"]
                            break
                    for a in ALIASES:
                        Pa = float(a) * P
                        if abs(Pg - Pa) <= ptol(Pa, T14, span):
                            ali = f"{a.numerator}:{a.denominator}"
                            break
                full = bool(per and t0ok and durok and depok)
                if per:
                    cls = "full_pass" if full else "+".join(k for k, v in (("t0", t0ok), ("duration", durok), ("depth", depok)) if not v)
                elif sib:
                    cls = "found_sibling_instead"
                elif ali:
                    cls = "alias_" + ali
                else:
                    cls = "period_mismatch"
                o.update(period_pass=bool(per), t0_pass=bool(t0ok), duration_pass=bool(durok), depth_pass=depok, full_pass=full,
                         alias_ratio=ali, sibling=sib, P_tol_d=tol, dP_d=dP, dt_h=24 * dt, sigma_pred_h=24 * sig,
                         t0_tol_h=24 * t0tol, depth_ratio=ratio, n_pass=int(per) + int(t0ok) + int(durok) + int(depok),
                         reason_class=cls)
            out.append(o)
    s = pd.DataFrame(out)
    # planet level: full_pass if any sector passes, else best row (most criteria, then lowest sector).
    # A planet's label is outside_period_range only if every sector row is (per-star pmax differs by sector).
    rank = s.assign(_fp=-s.full_pass.astype(int), _sc=-s.scored.astype(int), _np=-s.n_pass)
    best = rank.sort_values(["tic", "name", "_fp", "_sc", "_np", "sector"]).groupby(["tic", "name"], sort=True).head(1)
    best = best.drop(columns=["_fp", "_sc", "_np"]).reset_index(drop=True)
    nsec = s.groupby(["tic", "name"]).agg(n_sectors=("sector", "size"), n_sectors_full_pass=("full_pass", "sum")).reset_index()
    p = best.merge(nsec, on=["tic", "name"])
    assert len(p) == len(cp)
    # siblings credited: target planets whose best row found a sibling; the sibling counts as recovered by that row
    sib_rows = s[s.reason_class == "found_sibling_instead"]
    sib_names = set(zip(sib_rows.tic, sib_rows.sibling))
    p["credited_via_sibling"] = [(t, n) in sib_names for t, n in zip(p.tic, p.name)]
    s.to_csv(D / "xcheck_planet_sector_v2.csv", index=False)
    p.to_csv(D / "xcheck_planet_v2.csv", index=False)
    p[~p.reason_class.isin(["full_pass", "outside_grid", "outside_period_range"])].to_csv(D / "xcheck_mismatches_v2.csv", index=False)
    sc, ss = p[p.scored], s[s.scored]
    cnt = lambda c: int(c.fillna(False).astype(bool).sum())
    summ = dict(
        n_planets_in_scope=len(p), n_planet_sector_rows=len(s),
        labels_planet_level={k: int((p.reason_class == k).sum()) for k in ("outside_grid", "outside_period_range", "catalogue_incomplete", "no_span")},
        n_scored=len(sc), n_scored_planet_sector_rows=len(ss),
        planet_level_pass={k: cnt(sc[k + "_pass"]) for k in ("period", "t0", "duration", "depth")} | {"full_pass": cnt(sc.full_pass)},
        planet_sector_level_pass={k: cnt(ss[k + "_pass"]) for k in ("period", "t0", "duration", "depth")} | {"full_pass": cnt(ss.full_pass)},
        counts_by_class_planet=p.reason_class.value_counts().sort_index().to_dict(),
        counts_by_class_planet_sector=s.reason_class.value_counts().sort_index().to_dict(),
        aliases_by_ratio_planet={k: int((p.alias_ratio == k).sum()) for k in ("1:2", "2:1", "1:3", "3:1", "2:3", "3:2")},
        aliases_by_ratio_planet_sector={k: int((s.alias_ratio == k).sum()) for k in ("1:2", "2:1", "1:3", "3:1", "2:3", "3:2")},
        n_found_sibling_instead_planet=int((p.reason_class == "found_sibling_instead").sum()),
        n_sibling_planets_credited=len(sib_names), n_sibling_credited_not_otherwise_full_pass=int((p.credited_via_sibling & ~p.full_pass).sum()),
        t_span_vs_grid_baseline_median_diff_d=float((s.T_span_d - s.grid_baseline_d).median()),
        n_rows_loaded=len(rows), pscomppars_file=ps_path.name, asserts_passed=True)
    (D / "xcheck_summary_v2.json").write_text(json.dumps(summ, indent=1, default=int) + "\n")
    print(json.dumps(summ, indent=1, default=int))


if __name__ == "__main__":
    main()
