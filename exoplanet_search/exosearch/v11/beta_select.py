"""Apply the registered beta-treatment selection rule (docs/preregistration/beta_treatment.md) exactly as written.

Inputs (nulls and injections only; no candidate row): the v1.0 interim flipped null and scrambled realization r0 for
S70 and S106, the v1.0 interim injections (S70 runs 0-1, S106 run 0), star_quality_v1. Per sector: tuning
(non-holdout) stars, star-quality cut applied, statistic S_exact, options a-d applied at the stored peak
(exosearch.stat_options.variants; S' = S_exact * beta / beta').

1. Gate: the 95% star-level bootstrap CI (paired over stars present in both nulls, 2,000 resamples, seed 20261110)
   of S*_flip - S*_scr contains 0, in every sector. Options failing in any sector are out.
2. Among passing options: highest vetted completeness at expected S 15-25 (allowed injection range; S70 and S106
   injections pooled, each scored at its own sector's S* = max(S*_flip, S*_scr)).
3. Differences within 2 percentage points are ties; ties go to the option that modifies the fewest null peaks by more
   than 1% (flipped + scrambled tuning/quality null peaks of both sectors, relative to the unmodified S_exact).
4. No option passes: b, and the failure is reported.
Disclosed alongside (not part of the rule as written): an equivalence (TOST) reading of step 1 (90% CI of the
difference inside +-0.5) and a paired-bootstrap tie test for step 3.

  python -m exosearch.v11.beta_select --out DIR [--dry-run-flip-only]
"""
from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMBA_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[_v] = "1"          # before numpy is imported: one thread per worker (reproducible sums, no oversubscription)

import argparse
import json
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

from .. import interim, interim_review as rv, quality, stat_options

SECTORS = (70, 106)
OPTIONS = ("a", "b", "c", "d")
N_BOOT = 2000
SEED = 20261110
TIE_PP = 0.02
TOST_MARGIN = 0.5
SNR_BIN = (15, 25)
ORDER_SMALLEST = None          # decided by the counts (rule 3), not by a fixed order


def null_frame(tag, q, sector):
    d = interim.with_quality(interim.load_null(tag), q)
    x = d[(d.sector == sector) & ~d["drop"] & ~d.holdout.astype(bool)].copy()
    x = x.merge(q[["tic", "sector", "n_cad"]], on=["tic", "sector"], how="left")
    x["S_exact"] = x.depth * np.sqrt(x.n_tr) / x.cdpp_d
    return x.drop_duplicates("tic").set_index("tic")


def option_S(x, pri, med_ncad):
    v, _ = stat_options.variants(x.S_exact.values, x.beta_d.values, x.duration.values,
                                 x.n_cad.fillna(med_ncad).values, pri)
    return v


def _sstar(S, vetted, n, p, u):
    Sv = np.where(vetted, S, -np.inf)
    Sv = Sv[np.isfinite(Sv)]
    return interim.gpd_threshold(Sv, n, p, u)["s_star"]


def sector_option(args):
    sector, opt, flip_tag, scr_tag = args
    q = quality.load()
    pri = stat_options.prior(q[q.sector == sector])
    N = interim.search_population([sector], q)[sector]["kept"]
    p = 1.0 / N
    xf, xs = null_frame(flip_tag, q, sector), null_frame(scr_tag, q, sector)
    common = xf.index.intersection(xs.index)
    xf, xs = xf.loc[common], xs.loc[common]
    med = float(q[q.sector == sector].n_cad.median())
    Sf, Ss = option_S(xf, pri, med)[opt], option_S(xs, pri, med)[opt]
    Sa_f, Sa_s = option_S(xf, pri, med)["a"], option_S(xs, pri, med)["a"]
    vf, vs = xf.vetted.values.astype(bool), xs.vetted.values.astype(bool)
    _, uf, _ = rv.gpd_diagnostics(pd.DataFrame(dict(S=Sf, vetted=vf)), N, sector, plot=False)
    _, us, _ = rv.gpd_diagnostics(pd.DataFrame(dict(S=Ss, vetted=vs)), N, sector, plot=False)
    n = len(common)
    s_f, s_s = _sstar(Sf, vf, n, p, uf), _sstar(Ss, vs, n, p, us)
    rng = np.random.default_rng(SEED)
    diffs, bf, bs = [], [], []
    for _ in range(N_BOOT):
        i = rng.integers(0, n, n)
        a, b = _sstar(Sf[i], vf[i], n, p, uf), _sstar(Ss[i], vs[i], n, p, us)
        if a is not None and b is not None and np.isfinite(a) and np.isfinite(b):
            diffs.append(a - b); bf.append(a); bs.append(b)
    diffs = np.array(diffs)
    q95 = [float(np.quantile(diffs, 0.025)), float(np.quantile(diffs, 0.975))]
    q90 = [float(np.quantile(diffs, 0.05)), float(np.quantile(diffs, 0.95))]
    moved = int((np.abs(Sf / Sa_f - 1) > 0.01).sum() + (np.abs(Ss / Sa_s - 1) > 0.01).sum())
    return dict(sector=sector, option=opt, n_null_stars=int(n), n_stars_alpha=int(N), u_flip=float(uf), u_scr=float(us),
                s_star_flip=s_f, s_star_scr=s_s, s_star=max(s_f, s_s), diff=s_f - s_s,
                ci95_diff=q95, ci90_diff=q90, n_boot_ok=int(len(diffs)),
                ci95_flip=[float(np.quantile(bf, .025)), float(np.quantile(bf, .975))],
                ci95_scr=[float(np.quantile(bs, .025)), float(np.quantile(bs, .975))],
                gate_pass=bool(q95[0] <= 0 <= q95[1]),
                tost_pass=bool(-TOST_MARGIN < q90[0] and q90[1] < TOST_MARGIN),
                null_peaks_moved_gt1pct=moved, vetted_flip=int(vf.sum()), vetted_scr=int(vs.sum()))


def injections(q):
    inj = interim.read_parts(interim.PROD / "injections" / "*" / "part-*.parquet")
    inj = inj[(inj.status == "ok") & inj.sector.isin(SECTORS) & ~inj.holdout.astype(bool)]
    inj = inj.merge(q[["tic", "sector", "beta_max", "status", "n_cad"]].rename(columns={"status": "q_status"}),
                    on=["tic", "sector"], how="left")
    inj = inj[inj.q_status.eq("ok") & (inj.beta_max <= quality.rule()["threshold"])].copy()
    inj["S_exact"] = inj.depth_rec * np.sqrt(inj.n_tr_rec) / inj.cdpp_d_rec
    inj["allowed"] = (inj.depth_inj < rv.MAX_DEPTH) & (inj.rp_inj < rv.MAX_RP)
    failed = inj.failed_checks.fillna("").str.split(",").map(lambda f: set(f) - {"", "low_s"})
    inj["clean"] = failed.map(len) == 0
    inj["inbin"] = inj.allowed & (inj.snr_exp >= SNR_BIN[0]) & (inj.snr_exp < SNR_BIN[1])
    return inj.drop_duplicates("key")


def completeness(inj, q, s_star):
    """-> per option: pooled vetted completeness in the bin, per-sector, and the per-row vetted indicator."""
    out, ind = {}, {}
    for opt in OPTIONS:
        v_all = np.zeros(len(inj), bool)
        per = {}
        for s in SECTORS:
            m = (inj.sector == s).values
            x = inj[m]
            pri = stat_options.prior(q[q.sector == s])
            med = float(q[q.sector == s].n_cad.median())
            vi, _ = stat_options.variants(x.S_exact.values, x.beta_d.values, x.dur_rec.values, x.n_cad.fillna(med).values, pri)
            det = x.matched.astype(bool).values & (vi[opt] >= s_star[(s, opt)])
            v_all[m] = det & x.clean.values
            b = x.inbin.values
            per[str(s)] = dict(n=int(b.sum()), vetted=float(v_all[m][b].mean()) if b.any() else None)
        b = inj.inbin.values
        out[opt] = dict(n=int(b.sum()), vetted=float(v_all[b].mean()), per_sector=per)
        ind[opt] = v_all[b]
    return out, ind


def choose(res, comp, ind):
    gate = {o: all(res[(s, o)]["gate_pass"] for s in SECTORS) for o in OPTIONS}
    tost = {o: all(res[(s, o)]["tost_pass"] for s in SECTORS) for o in OPTIONS}
    moved = {o: sum(res[(s, o)]["null_peaks_moved_gt1pct"] for s in SECTORS) for o in OPTIONS}

    def pick(passing, tie_fn):
        if not passing:
            return "b", "no option passes the gate: b (most conservative) per rule 4"
        top = max(comp[o]["vetted"] for o in passing)
        tied = [o for o in passing if tie_fn(o, top)]
        best = min(tied, key=lambda o: (moved[o], OPTIONS.index(o)))
        return best, f"tied within rule 3: {tied}; fewest null peaks moved > 1%: {best} ({moved[best]})"

    passing = [o for o in OPTIONS if gate[o]]
    chosen, why = pick(passing, lambda o, top: top - comp[o]["vetted"] <= TIE_PP + 1e-12)
    # disclosures: TOST gate and paired-bootstrap ties
    rng = np.random.default_rng(SEED)
    n = len(next(iter(ind.values())))
    boots = {o: [] for o in OPTIONS}
    for _ in range(N_BOOT):
        i = rng.integers(0, n, n)
        for o in OPTIONS:
            boots[o].append(ind[o][i].mean())
    boots = {o: np.array(v) for o, v in boots.items()}
    tpass = [o for o in OPTIONS if tost[o]]

    def paired_tie(o, top):
        lead = max(tpass, key=lambda k: comp[k]["vetted"])
        d = boots[lead] - boots[o]
        return o == lead or np.quantile(d, 0.025) <= 0 <= np.quantile(d, 0.975)
    alt, alt_why = pick(tpass, paired_tie) if tpass else ("b", "no option passes the TOST gate")
    return dict(gate=gate, tost=tost, moved=moved, passing=passing, chosen=chosen, reason=why,
                disclosure_tost_paired=dict(passing=tpass, chosen=alt, reason=alt_why))


def main(out, flip_tag="invert", scr_tag="scramble_r0", workers=8):
    t0 = time.time()
    q = quality.load()
    with ProcessPoolExecutor(max_workers=workers) as ex:
        rows = list(ex.map(sector_option, [(s, o, flip_tag, scr_tag) for s in SECTORS for o in OPTIONS]))
    res = {(r["sector"], r["option"]): r for r in rows}
    inj = injections(q)
    comp, ind = completeness(inj, q, {k: v["s_star"] for k, v in res.items()})
    ch = choose(res, comp, ind)
    rep = dict(rule="docs/preregistration/beta_treatment.md (registered 2026-11-10), applied as written",
               inputs=dict(flip=flip_tag, scrambled=scr_tag, sectors=list(SECTORS), n_boot=N_BOOT, seed=SEED,
                           injections_in_bin=int(inj.inbin.sum())),
               per_sector_option=[res[(s, o)] for s in SECTORS for o in OPTIONS], completeness_15_25=comp,
               decision=ch, runtime_s=time.time() - t0)
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    (out / "beta_selection.json").write_text(json.dumps(rep, indent=1, default=float))
    print(json.dumps(dict(decision=ch, completeness={o: comp[o]["vetted"] for o in OPTIONS},
                          gate={f"{s}{o}": [round(res[(s, o)]["s_star_flip"], 3), round(res[(s, o)]["s_star_scr"], 3),
                                            [round(x, 3) for x in res[(s, o)]["ci95_diff"]]] for s in SECTORS for o in OPTIONS}),
                     indent=1, default=float))
    return rep


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--flip", default="invert")
    ap.add_argument("--scr", default="scramble_r0")
    ap.add_argument("--workers", type=int, default=8)
    a = ap.parse_args()
    main(a.out, a.flip, a.scr, a.workers)
