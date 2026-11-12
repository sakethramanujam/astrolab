"""Reader/adapter: the registered beta selection (exosearch.v11.beta_select, unchanged) on v1.1 null rows.

The v1.1 amendment (2026-11-11, section 5) reruns the beta TOST as a new registration on the v1.1 sector-70 nulls
(flipped `invert` and scrambled `scramble_r0`, written by `exosearch.v11.run work`). beta_select reads the v1.0
`exosearch.nulltest` products; this module reads the v1.1 row directories instead and hands the same quantities to the
unchanged rule. Nothing here is on the search path (not imported by exosearch.v11.run or exosearch.v11.pipeline).

Reader (fails closed, RowsError):
  * merge the row directories (`<dir>/parts/part-*.parquet`), dedupe by (tic, task); duplicates must agree on the
    statistic columns (|dS| <= 1e-6, identical periods) and an ok row supersedes an error row of the same task;
  * one sector (70), only the two null tasks, one tree, commit, config hash and pipeline version (optionally equal to
    the expected tree), one production flag;
  * every task row ok, counts equal to the expected counts (12,797 flipped, 12,512 scrambled; 25,309) and, when
    manifests are given, the (tic, task) set equal to the manifests' task set (nothing missing, nothing extra);
  * no scrambled row for an excluded host, and the rows' exclusion-file sha256 equals the registered file's.
Vetting: as interim.load_null (v1.0 interim layer): common-period flag among each null's own peaks
  (vet.flag_common_periods per sector and task), then vetted = no screen flag except low_s, no common_period, no deep
  REJECT flag, checks logged.
Populations (amendment section 5): tuning (non-holdout) stars passing the star-quality cut; the 285 excluded hosts have
  no scrambled rows, so the rule's pairing removes them from the flipped side; registration_checks' same-population
  guard is run on the paired tables and stops the registration if any excluded TIC is paired. alpha = 1 / N with N =
  searched sector-70 stars (the flipped-null population, holdout included) passing the quality cut, as
  interim.search_population's `kept`.
Statistic per option: the v1.1 row's S_opt_<o> (each option searched at its own peak by exosearch.v11.search).
Rule: TOST (beta_select.tost, margin 0.5 on the 90% CI), bootstrap (2,000, seed 20261110), GPD threshold choice
  (interim_review.gpd_diagnostics) and S* (beta_select._sstar), completeness (beta_select.completeness), and the choice
  (beta_select.pick, paired_boot) are beta_select's own functions; only the sector list is S70 (amendment section 4).
"""
from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMBA_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[_v] = "1"          # as beta_select: before numpy is imported (also when run as __main__)

import glob
import gzip
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from . import beta_select as bs
from . import registration_checks
from .. import deepvet, interim, interim_review as rv, vet

TASKS = ("invert", "scramble_r0")
EXPECTED = {"invert": 12797, "scramble_r0": 12512}
SECTOR = 70
PROV_COLS = ("tree", "commit", "config_hash", "pipeline_version")
STAT_COLS = ("S",) + tuple(f"S_opt_{o}" for o in bs.OPTIONS)
EXACT_COLS = ("P_grid", "period", "flags", "deep_flags", "checks_logged") + tuple(f"P_grid_opt_{o}" for o in bs.OPTIONS)
TOL = 1e-6


class RowsError(RuntimeError):
    """The v1.1 null rows are not a complete, single-tree, all-ok set. Registration stops."""


def read_dirs(dirs):
    frames = []
    for i, d in enumerate(dirs):
        fs = sorted(glob.glob(str(Path(d) / "parts" / "part-*.parquet")))
        if not fs:
            raise RowsError(f"no part files in row directory #{i}")
        x = pd.concat([pd.read_parquet(f) for f in fs], ignore_index=True)
        x["_src"] = i
        frames.append(x)
    return pd.concat(frames, ignore_index=True)


def manifest_tasks(paths):
    """(tic, task) pairs of null manifests (exosearch.v11.run.write_manifests format)."""
    keys = set()
    for p in paths:
        with gzip.open(p, "rt") as fh:
            fh.readline()
            for line in fh:
                if line.strip():
                    j = json.loads(line)
                    keys |= {(int(j["tic"]), t["tag"]) for t in j["tasks"]}
    return keys


def _single(df, col, expect=None):
    vals = sorted(map(str, df[col].dropna().unique())) if col in df else []
    if len(vals) != 1 or df[col].isna().any():
        raise RowsError(f"{col}: expected one value, found {vals[:5]}{' and missing values' if col in df and df[col].isna().any() else ''}")
    if expect is not None and vals[0] != str(expect):
        raise RowsError(f"{col} {vals[0]} != expected {expect}")
    return vals[0]


def _agree(g):
    ok = g[g.status == "ok"]
    if len(ok) < 2:
        return True
    missing = [c for c in STAT_COLS + EXACT_COLS if c not in ok]
    if missing:
        raise RowsError(f"row columns missing: {missing}")
    a = ok.iloc[0]
    for _, b in ok.iloc[1:].iterrows():
        for c in STAT_COLS:
            if not ((pd.isna(a[c]) and pd.isna(b[c])) or abs(float(a[c]) - float(b[c])) <= TOL):
                return False
        for c in EXACT_COLS:
            if not ((pd.isna(a[c]) and pd.isna(b[c])) or a[c] == b[c]):
                return False
    return True


def merge_rows(df, expected=EXPECTED, tree=None, sector=SECTOR, manifest_keys=None, exclusion_tics=None,
               exclusion_sha256=None):
    """-> (deduplicated rows, report). Raises RowsError on any departure from a complete single-tree all-ok set."""
    rep = dict(n_rows_read=int(len(df)))
    if not len(df):
        raise RowsError("no rows")
    if set(df.sector.astype(int)) != {sector}:
        raise RowsError(f"sectors {sorted(set(df.sector.astype(int)))} != [{sector}]")
    bad = set(df.task) - set(TASKS)
    if bad:
        raise RowsError(f"unexpected tasks {sorted(bad)}")
    prov = {c: _single(df, c, tree if c == "tree" else None) for c in PROV_COLS}
    prov["production"] = _single(df, "production")
    df = df.assign(_ok=(df.status == "ok")).sort_values(["tic", "task", "_ok", "_src"], ascending=[True, True, False, True])
    dup = df.duplicated(["tic", "task"], keep=False)
    rep["n_duplicate_rows"] = int(dup.sum())
    if dup.any():
        conflicts = [k for k, g in df[dup].groupby(["tic", "task"]) if not _agree(g)]
        if conflicts:
            raise RowsError(f"{len(conflicts)} duplicated (tic, task) with different results, e.g. {conflicts[:3]}")
    x = df.drop_duplicates(["tic", "task"], keep="first").drop(columns=["_ok"])
    rep["n_after_dedupe"] = int(len(x))
    nok = x[x.status != "ok"]
    if len(nok):
        raise RowsError(f"{len(nok)} task(s) without an ok row, e.g. {nok[['tic', 'task']].head(3).values.tolist()}")
    n = {t: int((x.task == t).sum()) for t in TASKS}
    rep["n_per_task"] = n
    if expected is not None:
        if n != dict(expected) or sum(n.values()) != sum(expected.values()):
            raise RowsError(f"row counts {n} (total {sum(n.values())}) != expected {dict(expected)} "
                            f"(total {sum(expected.values())})")
    if manifest_keys is not None:
        have = set(zip(x.tic.astype(int), x.task))
        miss, extra = manifest_keys - have, have - manifest_keys
        if miss or extra:
            raise RowsError(f"rows vs manifests: {len(miss)} missing (e.g. {sorted(miss)[:3]}), "
                            f"{len(extra)} extra (e.g. {sorted(extra)[:3]})")
        rep["manifest_tasks"] = len(manifest_keys)
    if exclusion_tics is not None:
        leak = set(x[x.task == "scramble_r0"].tic.astype(int)) & set(exclusion_tics)
        if leak:
            raise RowsError(f"{len(leak)} excluded host(s) have scrambled rows, e.g. {sorted(leak)[:3]}")
    if exclusion_sha256 is not None:
        _single(x, f"scramble_exclude_sha256_{sector}", exclusion_sha256)
    rep["provenance"] = prov
    return x.drop(columns=["_src"]), rep


def with_vetting(x):
    """interim.load_null's vetted column on v1.1 rows (common period per sector and task among the null's own peaks)."""
    out = []
    for _, g in x.groupby(["sector", "task"]):
        g = vet.flag_common_periods(g.assign(status="ok"))
        g["fail_common_period"] = g["flags"].fillna("").str.contains("common_period")
        out.append(g)
    d = pd.concat(out, ignore_index=True)
    fl = d["flags"].fillna("").str.split(",").map(lambda f: set(f) - {"", "low_s"})
    dfl = d.deep_flags.fillna("").str.split(",").map(set)
    d["vetted"] = (fl.map(len) == 0) & ~d.fail_common_period.fillna(False).astype(bool) \
        & dfl.map(lambda f: not (f & deepvet.REJECT)) & d.checks_logged.fillna(False).astype(bool)
    return d


def beta_frames(d, q, exclusion, sector=SECTOR):
    """-> (flipped frame, scrambled frame, report): tuning, quality-passing stars, indexed by tic, paired."""
    d = interim.with_quality(d[d.sector == sector], q)
    flip_all = d[d.task == "invert"]
    n_alpha = int((~flip_all["drop"]).sum())              # searched stars (holdout included) passing the cut
    sel = d[~d["drop"] & ~d.holdout.astype(bool)]
    xf = sel[sel.task == "invert"].drop_duplicates("tic").set_index("tic")
    xs = sel[sel.task == "scramble_r0"].drop_duplicates("tic").set_index("tic")
    guard = registration_checks.check_registration_population(xf.reset_index(), xs.reset_index(), exclusion, sector=sector)
    common = xf.index.intersection(xs.index)
    ex = set(pd.DataFrame(exclusion).tic.astype(int))
    rep = dict(n_alpha=n_alpha, flipped_tuning_quality=int(len(xf)), scrambled_tuning_quality=int(len(xs)),
               paired=int(len(common)), excluded_removed_from_flipped=int(len(set(xf.index) & ex)),
               flipped_unpaired_not_excluded=int(len(set(xf.index) - set(common) - ex)),
               scrambled_unpaired=int(len(set(xs.index) - set(common))), guard=guard)
    return xf.loc[common], xs.loc[common], rep


def sector_gate(xf, xs, opt, N, sector=SECTOR):
    """beta_select.sector_option from the paired frames on: same bootstrap, threshold choice, S* and TOST call;
    the option statistic is the v1.1 row's S_opt_<opt>."""
    p = 1.0 / N
    Sf, Ss = xf[f"S_opt_{opt}"].values.astype(float), xs[f"S_opt_{opt}"].values.astype(float)
    Sa_f, Sa_s = xf["S_opt_a"].values.astype(float), xs["S_opt_a"].values.astype(float)
    vf, vs = xf.vetted.values.astype(bool), xs.vetted.values.astype(bool)
    _, uf, _ = rv.gpd_diagnostics(pd.DataFrame(dict(S=Sf, vetted=vf)), N, sector, plot=False)
    _, us, _ = rv.gpd_diagnostics(pd.DataFrame(dict(S=Ss, vetted=vs)), N, sector, plot=False)
    n = len(xf)
    s_f, s_s = bs._sstar(Sf, vf, n, p, uf), bs._sstar(Ss, vs, n, p, us)
    rng = np.random.default_rng(bs.SEED)
    diffs, bf, bsr = [], [], []
    for _ in range(bs.N_BOOT):
        i = rng.integers(0, n, n)
        a, b = bs._sstar(Sf[i], vf[i], n, p, uf), bs._sstar(Ss[i], vs[i], n, p, us)
        if a is not None and b is not None and np.isfinite(a) and np.isfinite(b):
            diffs.append(a - b); bf.append(a); bsr.append(b)
    diffs = np.array(diffs)
    q95 = [float(np.quantile(diffs, 0.025)), float(np.quantile(diffs, 0.975))]
    q90 = [float(np.quantile(diffs, 0.05)), float(np.quantile(diffs, 0.95))]
    moved = int((np.abs(Sf / Sa_f - 1) > 0.01).sum() + (np.abs(Ss / Sa_s - 1) > 0.01).sum())
    return dict(sector=sector, option=opt, n_null_stars=int(n), n_stars_alpha=int(N), u_flip=float(uf), u_scr=float(us),
                s_star_flip=s_f, s_star_scr=s_s, s_star=max(s_f, s_s), diff=s_f - s_s,
                ci95_diff=q95, ci90_diff=q90, n_boot_ok=int(len(diffs)),
                ci95_flip=[float(np.quantile(bf, .025)), float(np.quantile(bf, .975))],
                ci95_scr=[float(np.quantile(bsr, .025)), float(np.quantile(bsr, .975))],
                ci95_contains_0=bool(q95[0] <= 0 <= q95[1]),
                gate_pass=bs.tost(q90),
                null_peaks_moved_gt1pct=moved, vetted_flip=int(vf.sum()), vetted_scr=int(vs.sum()))


def decide(res, comp=None, ind=None, sectors=(SECTOR,)):
    """beta_select.choose over the given sectors. Completeness (rule 2) is needed only when two or more options pass:
    with none, pick() gives b (rule 4); with one, rule 2 has a single candidate. Returns None when it is needed and
    missing (registration stops: needs_completeness)."""
    gate = {o: all(res[(s, o)]["gate_pass"] for s in sectors) for o in bs.OPTIONS}
    desc = {o: all(res[(s, o)]["ci95_contains_0"] for s in sectors) for o in bs.OPTIONS}
    moved = {o: sum(res[(s, o)]["null_peaks_moved_gt1pct"] for s in sectors) for o in bs.OPTIONS}
    passing = [o for o in bs.OPTIONS if gate[o]]
    diffs = None
    if comp is not None and ind is not None:
        boots = bs.paired_boot(ind)
        diffs = {f"{a}-{b}": [float(np.quantile(boots[a] - boots[b], q)) for q in (0.025, 0.975)]
                 for a in bs.OPTIONS for b in bs.OPTIONS if a < b}
    elif len(passing) >= 2:
        return None
    else:
        comp, boots = {o: {"vetted": float("nan")} for o in passing}, {}
    chosen, why, tied = bs.pick(passing, comp, boots)
    return dict(gate_tost=gate, passing=passing, chosen=chosen, reason=why, tied=tied, gate_failed=not passing,
                completeness_diff_ci95=diffs,
                descriptive=dict(ci95_contains_0=desc, null_peaks_moved_gt1pct=moved))


def injection_completeness(q, s_star):
    """beta_select.injections and .completeness on the sector-70 interim injections (the injection input the rule
    names), each option at its v1.1 S70 S*. S106 has no v1.1 S* in this run and enters with no rows."""
    if not glob.glob(str(interim.PROD / "injections" / "*" / "part-*.parquet")):
        raise FileNotFoundError("no interim injection part files under " + str(interim.PROD / "injections"))
    inj = bs.injections(q)
    inj = inj[inj.sector == SECTOR].copy()
    if not len(inj):
        raise FileNotFoundError("no sector-70 interim injections under " + str(interim.PROD / "injections"))
    st = dict(s_star)
    st.update({(s, o): np.inf for s in bs.SECTORS if s != SECTOR for o in bs.OPTIONS})
    comp, ind = bs.completeness(inj, q, st)
    return inj, comp, ind


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main(rows, out, manifests=(), tree=None, expected=EXPECTED, completeness="interim", workers=4):
    """Reader -> vetting -> population guard -> TOST per option (S70) -> decision. Writes out/beta_selection.json
    (beta_select's report layout plus a `v11_rows` section). Writes nothing else; sets no config value."""
    import time
    from concurrent.futures import ProcessPoolExecutor
    from .. import quality
    from . import settings
    t0 = time.time()
    c = settings.cfg()
    rel = settings.exclusion_files(c).get(str(SECTOR))
    if rel is None:
        raise RowsError("no registered sector-70 exclusion file")
    ex = pd.read_csv(settings.ROOT / rel)
    rows_df = read_dirs(rows)
    keys = manifest_tasks(manifests) if manifests else None
    x, rrep = merge_rows(rows_df, expected, tree, SECTOR, keys, set(ex.tic.astype(int)), sha256(settings.ROOT / rel))
    q = quality.load()
    if not (q.sector == SECTOR).any():
        raise RowsError("no sector-70 star-quality table (products/v1.0.0/analysis_layer/star_quality_s0070.parquet)")
    xf, xs, prep = beta_frames(with_vetting(x), q, ex)
    N = prep["n_alpha"]
    with ProcessPoolExecutor(max_workers=workers) as exe:
        res_rows = list(exe.map(_gate_job, [(xf, xs, o, N) for o in bs.OPTIONS]))
    res = {(r["sector"], r["option"]): r for r in res_rows}
    s_star = {k: v["s_star"] for k, v in res.items()}
    comp = ind = None
    inj_note, n_inj = None, None
    if completeness == "interim":
        try:
            inj, comp, ind = injection_completeness(q, s_star)
            n_inj = int(inj.inbin.sum())
            inj_note = "v1.0 interim sector-70 injections (beta_select.injections), scored at the v1.1 S70 S* per option"
        except FileNotFoundError as e:
            inj_note = f"unavailable: {e}"
    else:
        inj_note = "not computed (--completeness none)"
    ch = decide(res, comp, ind)
    status = "needs_completeness" if ch is None else ("gate_failed" if ch["gate_failed"] else "gate_passed")
    rep = dict(rule="docs/preregistration/beta_treatment.md (registered 2026-11-10) with amendment 1 (2026-11-11), applied as "
                    "written; v1.1 amendment 2026-11-11 sections 4-5 (sector 70, same population for both nulls)",
               inputs=dict(flip="invert", scrambled="scramble_r0", sectors=[SECTOR], n_boot=bs.N_BOOT, seed=bs.SEED,
                           injections_in_bin=n_inj, injections=inj_note,
                           statistic="v1.1 row S_opt_<option> (each option searched at its own peak)"),
               per_sector_option=[res[(SECTOR, o)] for o in bs.OPTIONS], completeness_15_25=comp,
               decision=ch, status=status,
               v11_rows=dict(reader=rrep, population=prep, exclusion_file=rel, exclusion_sha256=sha256(settings.ROOT / rel)),
               runtime_s=time.time() - t0)
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    (out / "beta_selection.json").write_text(json.dumps(rep, indent=1, default=_js))
    print(json.dumps(dict(status=status, decision=ch, n=rrep["n_per_task"], paired=prep["paired"],
                          gate={o: [round(res[(SECTOR, o)]["s_star_flip"], 3), round(res[(SECTOR, o)]["s_star_scr"], 3),
                                    [round(v, 3) for v in res[(SECTOR, o)]["ci90_diff"]]] for o in bs.OPTIONS}),
                     indent=1, default=_js))
    return rep


def _gate_job(a):
    return sector_gate(*a)


def _js(o):
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating, float)):
        return float(o)
    if isinstance(o, (set, tuple)):
        return list(o)
    return str(o)


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(description="beta selection on v1.1 null rows (reader + guard + unchanged rule)")
    ap.add_argument("--rows", nargs="+", required=True, help="row directories (each with parts/part-*.parquet)")
    ap.add_argument("--manifest", nargs="*", default=[], help="null manifests: the (tic, task) set must match exactly")
    ap.add_argument("--tree", help="expected tree hash of every row")
    ap.add_argument("--out", required=True)
    ap.add_argument("--completeness", choices=["interim", "none"], default="interim")
    ap.add_argument("--workers", type=int, default=4)
    a = ap.parse_args()
    try:
        main(a.rows, a.out, a.manifest, a.tree, completeness=a.completeness, workers=a.workers)
    except (RowsError, registration_checks.PopulationGuardError) as e:
        raise SystemExit(f"REGISTRATION STOPPED: {type(e).__name__}: {e}")
