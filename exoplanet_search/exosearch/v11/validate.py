"""Check (b) before any v1.1 production run: prewhitening + comb guard and the duration cap, on injections only.

Hosts and seeds: the v1.0 injection hosts of sector 70 (exosearch.inject.select_hosts, injection_v1) and their run-0
seeds, re-injected through the v1.1 code. Sets (each host star's tasks run on one host):
  main  every host, run-0 draw:            arms pw (prewhitening, comb guard on), nopw (off), noguard (guard off)
  usp   ~25% of hosts (seeded), P drawn log-uniform 0.5-1 d (USPs; harmonics above 5 c/d): same three arms
  dens  ~37.5% of hosts (seeded), "true" density always != TIC (binary / subgiant / log-uniform 0.01-100x, the
        registered wrong-density classes) while the search caps with TIC: arms cap (on) and nocap
Pass criteria (fixed here before any validation row exists; the report states each one):
  recovered = matched (P within 1%, t0 within 0.5 T14) and S >= 10 (reference; v1.1 S* is not frozen);
  rows: allowed range (depth < 5%, Rp < 25 R_E) and expected S >= 10.
  P1 recovery: in every duration, period and USP bin with n >= 30, recovery(pw) - recovery(nopw) is not below
     -2 percentage points with a paired-bootstrap 95% CI entirely below 0.
  P2 depth bias: in the same bins, median over rows matched in both arms of depth_rec/depth_inj (pw minus nopw) has
     |difference| <= 0.02, or a 95% CI that includes 0.
  P3 cap: in the dens set, for duration-error factors <= 3 (inside the 3x margin), recovery(cap) - recovery(nocap)
     is not below -2 pp with a 95% CI entirely below 0. Factors > 3 are reported, not judged.
  python -m exosearch.v11.validate build --weights h2:0.85,h3:0.15 --out DIR
  python -m exosearch.v11.validate report DIR/out_h2 DIR/out_h3 --out docs/v1.1_validation
"""
from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMBA_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[_v] = "1"          # before numpy is imported: one thread per worker (reproducible sums, no oversubscription)

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from . import run, settings

S_REF = 10.0
MIN_N = 30
LOSS_PP = 0.02
DEPTH_TOL = 0.02
N_BOOT = 2000
SEED = 20261110
DUR_BINS = [0, 1, 2, 3, 6, 100]                 # hours
PER_BINS = [0.5, 1, 3, 10, 30]                  # days
FAC_BINS = [0, 1, 2, 3, 5]                       # duration factor (rho_tic / rho_true)^(1/3)


def _u(tag, tic):
    return int.from_bytes(hashlib.sha256(f"v1.1-validation:{tag}:{int(tic)}".encode()).digest(), "big") / 2 ** 256


def build_validation(weights, out, sector=70):
    from .. import inject, store
    cfg_v1 = json.loads((settings.ROOT / "configs" / "injection_v1.json").read_text())
    hosts = set(inject.select_hosts(sector, cfg_v1, store.connect()).tic.astype(int))
    jobs = run.build([sector], weights, out, production=False, tics_filter=hosts, tasks_override=[], write=False)
    keep = []
    for j in jobs:
        rm = run.host_rm(j)
        if not rm:
            continue
        j.update(rstar=rm[0], mstar=rm[1])
        tic = j["tic"]
        seed = inject.seed_for(cfg_v1["seed"], sector, 0, tic)
        base = dict(mode="inject", run=0, seed=seed, wrong_density=False)
        tasks = [dict(base, tag="main_pw"), dict(base, tag="main_nopw", prewhiten=False),
                 dict(base, tag="main_noguard", comb_guard=False)]
        if _u("usp", tic) < 0.25:
            P = float(np.exp(np.log(0.5) + _u("usp-P", tic) * (np.log(1.0) - np.log(0.5))))
            ub = dict(base, draws=dict(P=P))
            tasks += [dict(ub, tag="usp_pw"), dict(ub, tag="usp_nopw", prewhiten=False), dict(ub, tag="usp_noguard", comb_guard=False)]
        if 0.25 <= _u("usp", tic) < 0.625:
            db = dict(base, wrong_density=True, force_density=True)
            tasks += [dict(db, tag="dens_cap"), dict(db, tag="dens_nocap", cap=False)]
        j["tasks"] = tasks
        keep.append(j)
    cfg = settings.cfg()
    cfg_inj = json.loads((settings.ROOT / "configs" / "injection_v1.1.json").read_text())
    run.write_manifests(keep, out, cfg, cfg_inj, weights, name="validation")
    return keep


# ---------------------------------------------------------------- report
def _paired(a, b, rng):
    n = len(a)
    if n == 0:
        return None, [None, None]
    d = []
    for _ in range(N_BOOT):
        i = rng.integers(0, n, n)
        d.append(a[i].mean() - b[i].mean())
    return float(a.mean() - b.mean()), [float(np.quantile(d, .025)), float(np.quantile(d, .975))]


def _paired_median(a, b, rng):
    n = len(a)
    if n < 5:
        return None, [None, None]
    d = a - b
    bs = [np.median(d[rng.integers(0, n, n)]) for _ in range(N_BOOT)]
    return float(np.median(d)), [float(np.quantile(bs, .025)), float(np.quantile(bs, .975))]


def compare_arms(df, set_, arm_a, arm_b, by, bins, label, rng, judge="recovery"):
    a = df[df.task == f"{set_}_{arm_a}"].set_index("tic")
    b = df[df.task == f"{set_}_{arm_b}"].set_index("tic")
    common = a.index.intersection(b.index)
    a, b = a.loc[common], b.loc[common]
    sel = a.allowed & (a.snr_exp >= S_REF)
    a, b = a[sel], b[sel]
    x = a[by] * (24 if by == "dur_inj" else 1)
    cat = pd.cut(x, bins, right=False)
    rows = []
    for k, idx in a.groupby(cat, observed=True).groups.items():
        ra, rb = a.loc[idx].recovered.values.astype(float), b.loc[idx].recovered.values.astype(float)
        dr, ci = _paired(ra, rb, rng)
        both = a.loc[idx].matched.fillna(False).values.astype(bool) & b.loc[idx].matched.fillna(False).values.astype(bool)
        dd, dci = _paired_median(a.loc[idx].depth_ratio.values[both], b.loc[idx].depth_ratio.values[both], rng)
        n = int(len(idx))
        r = dict(bin=f"{label} {k}", n=n, rec_a=float(ra.mean()), rec_b=float(rb.mean()), d_rec=dr, d_rec_ci95=ci,
                 n_both_matched=int(both.sum()),
                 depth_ratio_a=float(np.median(a.loc[idx].depth_ratio.values[both])) if both.any() else None,
                 depth_ratio_b=float(np.median(b.loc[idx].depth_ratio.values[both])) if both.any() else None,
                 d_depth_ratio=dd, d_depth_ratio_ci95=dci)
        judged = n >= MIN_N
        r["judged"] = judged
        r["fail_recovery"] = bool(judged and dr is not None and dr < -LOSS_PP and ci[1] < 0)
        r["fail_depth"] = bool(judge == "recovery" and judged and dd is not None and abs(dd) > DEPTH_TOL
                               and not (dci[0] <= 0 <= dci[1]))
        rows.append(r)
    return rows


def report(dirs, out):
    df = pd.concat([run.load_rows(d) for d in dirs], ignore_index=True)
    df = df[df.status == "ok"].copy()
    df["allowed"] = (df.depth_inj < 0.05) & (df.rp_inj < 25)
    df["recovered"] = df.matched.astype(bool) & (df.S >= S_REF)
    df["dur_factor"] = df.dens_factor ** (-1 / 3)
    rng = np.random.default_rng(SEED)
    res = dict(criteria=dict(S_ref=S_REF, min_n=MIN_N, loss_pp=LOSS_PP, depth_tol=DEPTH_TOL, n_boot=N_BOOT,
                             rows="allowed range and expected S >= 10", recovered="matched and S >= 10"),
               n_rows=int(len(df)), hosts=sorted(df.host_id.dropna().unique().tolist()),
               commit=sorted(df.commit.dropna().unique().tolist()), counts=df.task.value_counts().to_dict())
    tabs = {}
    tabs["prewhitening_main_by_duration"] = compare_arms(df, "main", "pw", "nopw", "dur_inj", DUR_BINS, "T14 h", rng)
    tabs["prewhitening_main_by_period"] = compare_arms(df, "main", "pw", "nopw", "P_inj", PER_BINS, "P d", rng)
    tabs["prewhitening_usp_by_duration"] = compare_arms(df, "usp", "pw", "nopw", "dur_inj", DUR_BINS, "USP T14 h", rng)
    tabs["comb_guard_main_by_period"] = compare_arms(df, "main", "pw", "noguard", "P_inj", PER_BINS, "P d", rng, judge="info")
    tabs["comb_guard_usp_by_duration"] = compare_arms(df, "usp", "pw", "noguard", "dur_inj", DUR_BINS, "USP T14 h", rng, judge="info")
    cap = compare_arms(df, "dens", "cap", "nocap", "dur_factor", FAC_BINS, "duration factor", rng, judge="cap")
    for r in cap:
        lo = float(r["bin"].split("[")[-1].split(",")[0]) if "[" in r["bin"] else 0
        inside = lo < 3
        r["judged"] = r["judged"] and inside
        r["fail_recovery"] = r["fail_recovery"] and inside
    tabs["cap_dens_by_duration_factor"] = cap
    d = df[df.task == "dens_cap"]
    tabs["cap_dens_by_class"] = [dict(cls=c, n=int(len(g)), recovered_cap=float(g.recovered.mean()),
                                      recovered_nocap=float(df[(df.task == "dens_nocap") & df.tic.isin(g.tic)].recovered.mean()))
                                 for c, g in d.groupby("dens_class")]
    for k in ("main", "usp"):
        x = df[df.task == f"{k}_pw"]
        tabs[f"prewhiten_activity_{k}"] = dict(stars=int(len(x)), any_terms=int((x.pw_prewhiten_n > 0).sum()),
                                              comb_protected=int((x.pw_prewhiten_comb > 0).sum()))
    p1 = [r["bin"] for k in ("prewhitening_main_by_duration", "prewhitening_main_by_period", "prewhitening_usp_by_duration")
          for r in tabs[k] if r["fail_recovery"]]
    p2 = [r["bin"] for k in ("prewhitening_main_by_duration", "prewhitening_main_by_period", "prewhitening_usp_by_duration")
          for r in tabs[k] if r["fail_depth"]]
    p3 = [r["bin"] for r in cap if r["fail_recovery"]]
    res.update(tables=tabs, P1_fail_bins=p1, P2_fail_bins=p2, P3_fail_bins=p3, all_pass=not (p1 or p2 or p3))
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    (out / "validation.json").write_text(json.dumps(res, indent=1, default=float))
    print(json.dumps(dict(all_pass=res["all_pass"], P1=p1, P2=p2, P3=p3, counts=res["counts"]), indent=1))
    return res


def main():
    ap = argparse.ArgumentParser()
    sp = ap.add_subparsers(dest="cmd", required=True)
    b = sp.add_parser("build"); b.add_argument("--weights", required=True); b.add_argument("--out", required=True)
    r = sp.add_parser("report"); r.add_argument("dirs", nargs="+"); r.add_argument("--out", required=True)
    a = ap.parse_args()
    if a.cmd == "build":
        build_validation([(h, float(x)) for h, x in (p.split(":") for p in a.weights.split(","))], a.out)
    else:
        report(a.dirs, a.out)


if __name__ == "__main__":
    main()
