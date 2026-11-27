#!/usr/bin/env python
"""Read-only comparison of v1.1 (vet_amend.revet) vs v1.2 (vet_final.final_verdict) verdicts on real and null rows.

Usage:
    python repro_old_vs_new.py --rows <run>/rows_known [<run>/rows_rest ...] --out <dir>

Given --rows DIR (one or more row directories, or a run directory containing rows_known / rows_rest),
computes for every ok real and null row:
  1. the original v1.1 verdict via exosearch.v11.vet_amend.revet + vet.flag_common_periods + checks_logged,
  2. the unified v1.2 verdict via runs/v1.1_s70_x86/finish.py (rebuild_scan -> exosearch.vet_final.final_verdict),
and writes a diff report (counts identical / differing, per-rule reasons for differences, needs_lc counts
per mode including inject if present) to --out DIR (summary.json + summary.md).  Never writes into any
--rows directory.
"""
from __future__ import annotations

import argparse
import gzip
import importlib.util
import json
import sys
from collections import Counter
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))   # exoplanet_search/

from exosearch import interim, vet, vet_final                  # noqa: E402
from exosearch.v11 import settings, vet_amend                  # noqa: E402

_FINISH_PATH = Path(__file__).resolve().parent / "finish.py"
_spec = importlib.util.spec_from_file_location("s70_finish", _FINISH_PATH)
finish = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(finish)


def load_rows(dirs: List[Path]) -> pd.DataFrame:
    """Read-only concatenation of row files under the given directories."""
    frames = []
    for d in dirs:
        if d.is_file() and d.suffix in (".parquet", ".csv"):
            frames.append(pd.read_parquet(d) if d.suffix == ".parquet" else pd.read_csv(d))
            continue
        if (d / "parts").is_dir():
            df = interim.read_parts(d / "parts" / "part-*.parquet")
            if len(df):
                frames.append(df)
                continue
        sub_parts = sorted(p.parent for p in d.glob("rows_*/parts/part-*.parquet"))
        if sub_parts:
            seen = []
            for sp in sub_parts:
                if sp not in seen:
                    seen.append(sp)
            for sp in seen:
                df = interim.read_parts(sp / "part-*.parquet")
                if len(df):
                    frames.append(df)
            if frames:
                continue
        fs = sorted(list(d.glob("rows_*.parquet")) + list(d.glob("rows_*.csv")))
        if not fs:
            fs = sorted(list(d.glob("*.parquet")) + list(d.glob("*.csv")))
        if not fs:
            raise SystemExit(f"no row parquet/csv files or parts/ found under {d}")
        frames.append(pd.concat([pd.read_parquet(f) if f.suffix == ".parquet" else pd.read_csv(f) for f in fs],
                                ignore_index=True))
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def load_manifests(rows_dirs: List[Path], manifest_args: Optional[List[str]] = None) -> Tuple[Dict[int, Dict[str, Any]], Dict[str, Any]]:
    """Locate and read manifest .jsonl.gz files (read-only) for TIC stellar density rho, light-curve
    filenames, and the run cfg."""
    candidates: List[Path] = []
    if manifest_args:
        for m in manifest_args:
            p = Path(m)
            if p.is_dir():
                candidates.extend(sorted(p.glob("*.jsonl.gz")))
            elif p.is_file():
                candidates.append(p)
    if not candidates:
        search_dirs: List[Path] = []
        for d in rows_dirs:
            base = d.parent if d.is_file() else d
            search_dirs.extend([base / "manifests", base.parent / "manifests", base, base.parent])
        search_dirs.append(Path(__file__).resolve().parent)
        for sd in search_dirs:
            if not sd.is_dir():
                continue
            found = [sd / f"{name}.jsonl.gz" for name in ("known_h1", "rest_h1") if (sd / f"{name}.jsonl.gz").is_file()]
            if not found:
                found = sorted(sd.glob("*.jsonl.gz"))
            if found:
                candidates = found
                break
    jobs: Dict[int, Dict[str, Any]] = {}
    cfg: Optional[Dict[str, Any]] = None
    for mf in candidates:
        with gzip.open(mf, "rt") as fh:
            first = fh.readline()
            if not first.strip():
                continue
            head = json.loads(first).get("header", {})
            if cfg is None and "cfg" in head:
                cfg = head["cfg"]
            for line in fh:
                if line.strip():
                    j = json.loads(line)
                    jobs[int(j["tic"])] = j
    if cfg is None:
        cfg = settings.cfg()
    return jobs, cfg


def _old_failing_rules(r: Dict[str, Any]) -> List[str]:
    """Failing rule list implied by the original v1.1 vet_amend.revet + common_period + checks_logged."""
    rej = sorted(vet_amend._set(r.get("reject_v2")) - {vet_final.COMMON_PERIOD})
    if vet_final.fail_common_period_of(r, r.get("fail_common_period")):
        rej.append(vet_final.COMMON_PERIOD)
    if not vet_final.checks_logged_of(r, r.get("checks_logged")):
        rej.append(vet_final.CHECKS_NOT_LOGGED)
    return rej


def run_old_v11(ok_rn: pd.DataFrame, jobs: Dict[int, Dict[str, Any]], cfg: Dict[str, Any], workers: int = 1) -> pd.DataFrame:
    """Original v1.1 finish.py path on real and null ok rows: vet_amend.revet -> flag_common_periods -> vetted_v2."""
    if not len(ok_rn):
        return pd.DataFrame()
    v = vet_amend.revet(ok_rn, jobs, cfg, workers=workers)
    out = []
    for _, g in v.groupby(["sector", "task"]):
        g = vet.flag_common_periods(g.assign(status="ok"))
        g["fail_common_period"] = g["flags"].fillna("").str.contains("common_period")
        out.append(g)
    v = pd.concat(out, ignore_index=True)
    v["vetted_v2"] = (
        v.passed_lc_checks_v2.astype(bool)
        & ~v.fail_common_period.astype(bool)
        & v.checks_logged.fillna(False).astype(bool)
    )
    rules = [_old_failing_rules(r) for r in v.to_dict("records")]
    v["old_first_failing"] = [rs[0] if rs else None for rs in rules]
    v["old_all_failing"] = [",".join(rs) for rs in rules]
    return v


def run_new_v12(ok: pd.DataFrame, jobs: Dict[int, Dict[str, Any]], cfg: Dict[str, Any]) -> pd.DataFrame:
    """New finish.py path on ok rows: rebuild_scan for real/null secondary_scan -> vet_final.final_verdict."""
    if not len(ok):
        return pd.DataFrame()
    scans = {}
    v_recs = []
    for idx, r in enumerate(ok.to_dict("records")):
        j = (jobs.get(int(r["tic"])) or {})
        sc = finish.rebuild_scan(r, j, cfg)
        scans[idx] = sc
        d = {k: x for k, x in vet_final.final_verdict_row(
                 r, checks_logged=r.get("checks_logged"), rho=j.get("rho"), scan=sc).items()
             if k not in vet_final.VERDICT_KEYS}
        d["_idx"] = idx
        v_recs.append(d)
    v = finish.flag_common_periods_by_task(pd.DataFrame(v_recs))
    lo, hi = [], []
    for r in v.to_dict("records"):
        sc = scans[int(r["_idx"])]
        kw = dict(jobs=jobs, cfg=cfg, fail_common_period=r["fail_common_period"], scan=sc)
        a = finish.path_for(r)(r, **kw)
        b = finish.path_for(r)(r, needs_lc="exempt", **kw)
        lo.append(dict(passed_final=a.passed, first_failing_rule=a.first_failing_rule,
                       all_failing=",".join(a.all_failing), status_final=a.status))
        hi.append(dict(passed_upper=b.passed, all_failing_upper=",".join(b.all_failing)))
    v = pd.concat([v.drop(columns=["_idx"]), pd.DataFrame(lo), pd.DataFrame(hi)], axis=1)
    v["vetted_v2"] = v.passed_final.astype(bool) & (v.status_final == vet_final.STATUS_OK)
    v["vetted_v2_upper"] = v.passed_upper.astype(bool)
    return v


def compare_rows(rows: pd.DataFrame, jobs: Dict[int, Dict[str, Any]], cfg: Dict[str, Any],
                 s_star: float = 13.71, workers: int = 1) -> Dict[str, Any]:
    """Compare v1.1 (vet_amend.revet) vs v1.2 (vet_final.final_verdict) across all ok real and null rows,
    and report needs_lc counts per mode (real, null, inject)."""
    ok = rows[rows.status == "ok"].copy().reset_index(drop=True)
    ok["_row_id"] = np.arange(len(ok))
    ok_rn = ok[ok["mode"] != "inject"].copy().reset_index(drop=True)

    v_old = run_old_v11(ok_rn, jobs, cfg, workers=workers)
    v_new = run_new_v12(ok, jobs, cfg)

    v_new_rn = v_new[v_new["mode"] != "inject"].copy()
    old_by_id = v_old.set_index("_row_id") if len(v_old) else pd.DataFrame()
    new_by_id = v_new_rn.set_index("_row_id") if len(v_new_rn) else pd.DataFrame()
    common_ids = sorted(set(old_by_id.index) & set(new_by_id.index))

    flips_old_pass_new_fail: Counter = Counter()
    flips_old_fail_new_pass: Counter = Counter()
    rules_in_new_only: Counter = Counter()
    rules_in_old_only: Counter = Counter()
    rule_set_transitions: Counter = Counter()
    differing_rows: List[Dict[str, Any]] = []

    n_identical = 0
    n_differing = 0
    n_identical_rules = 0
    n_differing_rules = 0

    for rid in common_ids:
        ro = old_by_id.loc[rid]
        rn = new_by_id.loc[rid]
        po = bool(ro["vetted_v2"])
        pn = bool(rn["vetted_v2"])
        old_str = str(ro["old_all_failing"] or "")
        new_str = str(rn["all_failing"] or "")
        old_set = [x for x in old_str.split(",") if x]
        new_set = [x for x in new_str.split(",") if x]
        if po == pn:
            n_identical += 1
        else:
            n_differing += 1
            if po and not pn:
                flips_old_pass_new_fail[str(rn["first_failing_rule"] or "none")] += 1
            else:
                flips_old_fail_new_pass[str(ro["old_first_failing"] or "none")] += 1
        if old_set == new_set and po == pn:
            n_identical_rules += 1
        else:
            n_differing_rules += 1
            for rname in sorted(set(new_set) - set(old_set)):
                rules_in_new_only[rname] += 1
            for rname in sorted(set(old_set) - set(new_set)):
                rules_in_old_only[rname] += 1
            rule_set_transitions[f"{old_str or 'PASS'} -> {new_str or 'PASS'}"] += 1
            differing_rows.append(dict(
                tic=int(rn["tic"]),
                sector=int(rn["sector"]) if "sector" in rn and pd.notna(rn["sector"]) else None,
                task=str(rn.get("task", "")),
                mode=str(rn.get("mode", "")),
                S=float(rn["S"]) if "S" in rn and pd.notna(rn["S"]) else None,
                old_passed=po,
                new_passed=pn,
                new_passed_upper=bool(rn["vetted_v2_upper"]),
                old_first_failing=ro["old_first_failing"],
                new_first_failing=rn["first_failing_rule"],
                old_all_failing=old_str,
                new_all_failing=new_str,
                status_final=str(rn["status_final"]),
            ))

    def _group_stat(ids: List[int]) -> Dict[str, int]:
        if not ids:
            return dict(n=0, n_identical=0, n_differing=0, n_identical_rules=0, n_differing_rules=0,
                        old_passed=0, new_passed=0, new_passed_upper=0,
                        old_passed_ge_sstar=0, new_passed_ge_sstar=0, needs_lc=0)
        sub_o = old_by_id.loc[ids]
        sub_n = new_by_id.loc[ids]
        po = sub_o["vetted_v2"].astype(bool)
        pn = sub_n["vetted_v2"].astype(bool)
        pnu = sub_n["vetted_v2_upper"].astype(bool)
        s_ge = sub_n["S"].astype(float) >= s_star
        same_p = po == pn
        same_r = same_p & (sub_o["old_all_failing"].fillna("") == sub_n["all_failing"].fillna(""))
        return dict(
            n=int(len(ids)),
            n_identical=int(same_p.sum()),
            n_differing=int((~same_p).sum()),
            n_identical_rules=int(same_r.sum()),
            n_differing_rules=int((~same_r).sum()),
            old_passed=int(po.sum()),
            new_passed=int(pn.sum()),
            new_passed_upper=int(pnu.sum()),
            old_passed_ge_sstar=int((po & s_ge).sum()),
            new_passed_ge_sstar=int((pn & s_ge).sum()),
            needs_lc=int((sub_n["status_final"] == vet_final.STATUS_NEEDS_LC).sum()),
        )

    real_ids = [i for i in common_ids if str(new_by_id.loc[i, "mode"]) == "real"]
    null_ids = [i for i in common_ids if str(new_by_id.loc[i, "mode"]) in ("invert", "scramble")]
    inv_ids = [i for i in common_ids if str(new_by_id.loc[i, "mode"]) == "invert"]
    scr_ids = [i for i in common_ids if str(new_by_id.loc[i, "mode"]) == "scramble"]

    by_mode = dict(
        real=_group_stat(real_ids),
        null=_group_stat(null_ids),
        invert=_group_stat(inv_ids),
        scramble=_group_stat(scr_ids),
    )
    by_task = {}
    if len(new_by_id):
        for t in sorted(new_by_id["task"].dropna().astype(str).unique()):
            t_ids = [i for i in common_ids if str(new_by_id.loc[i, "task"]) == t]
            by_task[t] = _group_stat(t_ids)

    nlc_all = (v_new["status_final"] == vet_final.STATUS_NEEDS_LC) if len(v_new) else pd.Series(dtype=bool)
    modes_all = v_new["mode"].astype(str) if len(v_new) else pd.Series(dtype=str)
    needs_lc = dict(
        real=int((nlc_all & (modes_all == "real")).sum()),
        null=int((nlc_all & modes_all.isin(["invert", "scramble"])).sum()),
        inject=int((nlc_all & (modes_all == "inject")).sum()),
        invert=int((nlc_all & (modes_all == "invert")).sum()),
        scramble=int((nlc_all & (modes_all == "scramble")).sum()),
        total=int(nlc_all.sum()),
    )

    return dict(
        s_star=float(s_star),
        n_rows_total=int(len(rows)),
        n_ok_total=int(len(ok)),
        n_real_null_ok=int(len(ok_rn)),
        n_inject_ok=int((ok["mode"] == "inject").sum()),
        n_compared=int(len(common_ids)),
        n_identical=int(n_identical),
        n_differing=int(n_differing),
        n_identical_rules=int(n_identical_rules),
        n_differing_rules=int(n_differing_rules),
        all_real_null_identical=bool(n_differing == 0 and n_differing_rules == 0),
        needs_lc=needs_lc,
        by_mode=by_mode,
        by_task=by_task,
        per_rule_differences=dict(
            verdict_flips_old_pass_new_fail={k: int(v) for k, v in sorted(flips_old_pass_new_fail.items())},
            verdict_flips_old_fail_new_pass={k: int(v) for k, v in sorted(flips_old_fail_new_pass.items())},
            rules_in_new_only={k: int(v) for k, v in sorted(rules_in_new_only.items())},
            rules_in_old_only={k: int(v) for k, v in sorted(rules_in_old_only.items())},
            rule_set_transitions={k: int(v) for k, v in sorted(rule_set_transitions.items())},
        ),
        differing_rows=differing_rows,
    )


def format_markdown(summ: Dict[str, Any]) -> str:
    lines = [
        "# v1.1 (`vet_amend.revet`) vs v1.2 (`vet_final.final_verdict`) Diff Report",
        "",
        f"- Real + null ok rows compared: **{summ['n_compared']}**",
        f"- Identical pass/fail verdicts: **{summ['n_identical']}**",
        f"- Differing pass/fail verdicts: **{summ['n_differing']}**",
        f"- Identical failing-rule sets: **{summ['n_identical_rules']}**",
        f"- Differing failing-rule sets: **{summ['n_differing_rules']}**",
        f"- `needs_lc` counts by mode: real = **{summ['needs_lc']['real']}**, "
        f"null = **{summ['needs_lc']['null']}** (invert = {summ['needs_lc']['invert']}, "
        f"scramble = {summ['needs_lc']['scramble']}), inject = **{summ['needs_lc']['inject']}**",
        "",
        "## Breakdown by Mode (real and null)",
        "",
        "| mode | n | identical | differing | old passed | new passed | old >= S* | new >= S* | needs_lc |",
        "|---|---|---|---|---|---|---|---|---|",
    ]
    for m, s in summ["by_mode"].items():
        lines.append(
            f"| {m} | {s['n']} | {s['n_identical']} | {s['n_differing']} | "
            f"{s['old_passed']} | {s['new_passed']} | {s['old_passed_ge_sstar']} | "
            f"{s['new_passed_ge_sstar']} | {s['needs_lc']} |"
        )
    lines += [
        "",
        "## Per-Rule Reasons for Differences",
        "",
    ]
    prd = summ["per_rule_differences"]
    if not any(prd.values()):
        lines.append("None (0 differences across all real and null rows).")
    else:
        lines.append("```json")
        lines.append(json.dumps(prd, indent=1))
        lines.append("```")
    lines.append("")
    return "\n".join(lines)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--rows", nargs="+", required=True, help="run row directories (read-only)")
    ap.add_argument("--out", required=True, help="output directory for diff report (must be disjoint from --rows)")
    ap.add_argument("--manifests", nargs="*", default=None, help="optional manifest .jsonl.gz files or directory")
    ap.add_argument("--sstar", type=float, default=13.71)
    ap.add_argument("--workers", type=int, default=1)
    a = ap.parse_args()

    rows_dirs = [Path(d) for d in a.rows]
    out = Path(a.out)
    out_res = out.resolve()
    rows_res = [(d.parent if d.is_file() else d).resolve() for d in rows_dirs]
    if any(out_res == r or out_res in r.parents or r in out_res.parents for r in rows_res):
        raise SystemExit("--out must be disjoint from the --rows directories (rows are read-only)")

    rows = load_rows(rows_dirs)
    jobs, cfg = load_manifests(rows_dirs, a.manifests)
    summ = compare_rows(rows, jobs, cfg, s_star=float(a.sstar), workers=int(a.workers))
    summ["rows_dirs"] = [d.name for d in rows_dirs]

    out.mkdir(parents=True, exist_ok=True)
    (out / "summary.json").write_text(json.dumps(summ, indent=1, default=str) + "\n")
    (out / "summary.md").write_text(format_markdown(summ))
    if summ["differing_rows"]:
        pd.DataFrame(summ["differing_rows"]).to_csv(out / "differing_rows.csv", index=False)
    print(json.dumps(summ, indent=1, default=str))


if __name__ == "__main__":
    main()
