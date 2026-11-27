#!/usr/bin/env python
"""Read-only re-scoring of a run's injection rows with the single final verdict.

Usage:
    python rescore_inject.py --rows <run>/rows_known [--rows <run>/rows_rest ...] --sstar 13.71 --out <dir>

--rows accepts one or more directories that hold the run's rows: either rows_*.parquet / rows_*.csv
files or parts/part-*.parquet shards (the layout exosearch.v11.run.load_rows reads).  Every ok
injection row is judged by exosearch.vet_final.final_verdict — exactly what finish.py applies to real,
null and inject rows now: vet_amend.amended() + common_period + checks_logged.  No light curve is
rebuilt and no manifest is read (a/R* falls back to the duration form, as vet_amend.amended does when
no TIC density is given), so deep 'secondary_scan' rows are status 'needs_lc' and are bracketed:
    survived_lower  counts them as rejected,
    survived_upper  counts them as exempted (only the undecided scan rule is waived).

Reports: injections with expected S >= S*, recovered, survived (both bounds), completeness restricted
to injections meeting the search's minimum transit count (exosearch.vet.THRESH['min_transits'], the
'few_transits' screen), per-rule first-failing / any-failing counts, and retention by injected T14.
Writes summary.json + summary.md into --out only; the rows directories are never written to.
"""
from __future__ import annotations

import argparse
import json
import sys
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))   # exoplanet_search/

from exosearch import vet, vet_final, interim    # noqa: E402

T14_BINS = (("<3h", 0.0, 0.125), ("3-6h", 0.125, 0.25), (">6h", 0.25, np.inf))   # T14 in days
MIN_TRANSITS = vet.THRESH["min_transits"]     # 'few_transits' screen of the v1.1 search


def load_rows(dirs: list[Path]) -> pd.DataFrame:
    """Concatenate the rows_*/parts row files of the given directories (read-only)."""
    frames = []
    for d in dirs:
        if (d / "parts").is_dir():
            df = interim.read_parts(d / "parts" / "part-*.parquet")
            if len(df):
                frames.append(df)
                continue
        fs = sorted(list(d.glob("rows_*.parquet")) + list(d.glob("rows_*.csv")))
        if not fs:
            raise SystemExit(f"no rows_*.parquet / rows_*.csv / parts/part-*.parquet under {d}")
        frames.append(pd.concat([pd.read_parquet(f) if f.suffix == ".parquet" else pd.read_csv(f) for f in fs],
                                ignore_index=True))
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def verdict_table(inj: pd.DataFrame) -> pd.DataFrame:
    """Final verdict columns for the ok injection rows: amended fields, then the group common-period
    flag (the same vet.flag_common_periods call finish.py makes), then both needs_lc bounds."""
    recs = [vet_final.final_verdict_row(r, checks_logged=r.get("checks_logged"), rho=None) for r in inj.to_dict("records")]
    v = pd.DataFrame([{k: x for k, x in r.items() if k not in vet_final.VERDICT_KEYS} for r in recs])
    flagged = []
    for _, g in v.groupby(["sector", "task"]):
        g = vet.flag_common_periods(g.assign(status="ok"))
        g["fail_common_period"] = g["flags"].fillna("").str.contains("common_period")
        flagged.append(g)
    v = pd.concat(flagged, ignore_index=True) if flagged else v.assign(fail_common_period=False)
    lo, hi = [], []
    for r in v.to_dict("records"):
        kw = dict(common_period=r["fail_common_period"], checks_logged=r.get("checks_logged"), rho=None)
        a = vet_final.final_verdict(r, **kw)
        b = vet_final.final_verdict(r, needs_lc="exempt", **kw)
        lo.append(dict(passed_final=a.passed, first_failing_rule=a.first_failing_rule,
                       all_failing=",".join(a.all_failing), status_final=a.status))
        hi.append(dict(passed_upper=b.passed, all_failing_upper=",".join(b.all_failing)))
    v = pd.concat([v, pd.DataFrame(lo), pd.DataFrame(hi)], axis=1)
    v["passed"] = v.passed_final.astype(bool) & (v.status_final == vet_final.STATUS_OK)
    v["passed_upper"] = v.passed_upper.astype(bool)
    return v


def _stat(g: pd.DataFrame, s_star: float) -> dict:
    n = int(len(g))
    rec = g.matched.fillna(False).astype(bool) & (g.S >= s_star)
    out = dict(n=n, expected_s_ge_sstar=int((g.snr_exp >= s_star).sum()),
               recovered=int(rec.sum()),
               survived_lower=int((rec & g.passed).sum()), survived_upper=int((rec & g.passed_upper).sum()))
    out.update(recovered_rate=round(float(rec.mean()), 4) if n else None,
               survived_rate_lower=round(float((rec & g.passed).mean()), 4) if n else None,
               survived_rate_upper=round(float((rec & g.passed_upper).mean()), 4) if n else None)
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--rows", nargs="+", required=True, help="run row directories (rows_* parquet/csv or parts/)")
    ap.add_argument("--sstar", type=float, default=13.71)
    ap.add_argument("--out", required=True, help="directory for summary.json / summary.md (never a rows dir)")
    a = ap.parse_args()
    rows_dirs = [Path(d) for d in a.rows]
    out = Path(a.out)
    out_res, rows_res = out.resolve(), [d.resolve() for d in rows_dirs]
    if any(out_res == r or out_res in r.parents or r in out_res.parents for r in rows_res):
        raise SystemExit("--out must be disjoint from the --rows directories (rows are read-only)")
    rows = load_rows(rows_dirs)
    inj_all = rows[rows["mode"] == "inject"]
    inj = inj_all[inj_all.status == "ok"].copy()
    v = verdict_table(inj)
    S = float(a.sstar)
    rec = v.matched.fillna(False).astype(bool) & (v.S >= S)
    nlc = v.status_final == vet_final.STATUS_NEEDS_LC
    per_rule_first = Counter("none" if x is None else str(x) for x in v.first_failing_rule)
    per_rule_any = Counter(x for s in v.all_failing for x in s.split(",") if x)
    per_rule_any_up = Counter(x for s in v.all_failing_upper for x in s.split(",") if x)
    met = v.n_transits.fillna(0) >= MIN_TRANSITS          # search's 'few_transits' screen
    summ = {
        "s_star": S,
        "verdict": ("exosearch.vet_final.final_verdict (vet_amend.amended + common_period + checks_logged); "
                    "no light-curve rebuilds, so deep 'secondary_scan' rows are 'needs_lc' and bracketed"),
        "rows_dirs": [d.name for d in rows_dirs],
        "min_transits": MIN_TRANSITS,
        "min_transits_source": "exosearch.vet.THRESH['min_transits'] (the v1.1 'few_transits' screen)",
        "n_total": int(len(inj_all)),
        "n_ok": int(len(inj)),
        "n_unphysical": int((inj_all.status == "unphysical").sum()),
        "n_error": int((~inj_all.status.isin(["ok", "unphysical"])).sum()),
        "expected_s_ge_sstar": int((v.snr_exp >= S).sum()),
        "recovered": dict(n=int(rec.sum()), rate=round(float(rec.mean()), 4) if len(v) else None),
        "survived": dict(lower=int((rec & v.passed).sum()), upper=int((rec & v.passed_upper).sum()),
                         needs_lc=int(nlc.sum()), needs_lc_gap=int((rec & v.passed_upper).sum() - (rec & v.passed).sum()),
                         note="recovered (matched and S >= S*) and passing the final verdict; "
                              "lower bound counts needs_lc rows as rejected, upper bound as exempted"),
        "vetted_ge_sstar": dict(lower=int((v.passed & (v.S >= S)).sum()), upper=int((v.passed_upper & (v.S >= S)).sum()),
                                note="passing the final verdict with S >= S*, whether or not matched"),
        "completeness_min_transits": dict(min_transits=MIN_TRANSITS,
                                          **_stat(v[met], S),
                                          note=f"restricted to ok injections with n_transits >= {MIN_TRANSITS}"),
        "per_rule": dict(first_failing={k: int(n) for k, n in sorted(per_rule_first.items())},
                         any_failing={k: int(n) for k, n in sorted(per_rule_any.items())},
                         any_failing_upper={k: int(n) for k, n in sorted(per_rule_any_up.items())},
                         n_needs_lc=int(nlc.sum())),
        "retention_by_t14": {},            # injected T14 (dur_inj, days)
        "by_task": {t: _stat(g, S) for t, g in v.groupby("task")},
    }
    dur = v["dur_inj"].astype(float) if "dur_inj" in v else pd.Series(np.nan, index=v.index)
    for name, lo_d, hi_d in T14_BINS:
        g = v[(dur >= lo_d) & (dur < hi_d)]
        summ["retention_by_t14"][name] = dict(_stat(g, S), dur_inj_days=[lo_d, hi_d])
    summ["retention_by_t14"]["unknown"] = dict(_stat(v[~(((dur >= T14_BINS[0][1]) & (dur < T14_BINS[-1][2])).fillna(False))], S),
                                               dur_inj_days=None, note="missing dur_inj")
    out.mkdir(parents=True, exist_ok=True)
    (out / "summary.json").write_text(json.dumps(summ, indent=1, default=str) + "\n")
    lines = [
        "# Injection re-score (single final verdict)",
        "",
        f"S* = {S}; verdict = `exosearch.vet_final.final_verdict` (vet_amend.amended + common_period + "
        f"checks_logged), the same function finish.py applies to real, null and inject rows.",
        f"Rows: {summ['n_total']} injections, {summ['n_ok']} ok, {summ['n_unphysical']} unphysical, {summ['n_error']} error.",
        f"Search minimum transit count: {MIN_TRANSITS} ({summ['min_transits_source']}).",
        "",
        f"- Expected S >= S*: **{summ['expected_s_ge_sstar']}**",
        f"- Recovered (matched and S >= S*): **{summ['recovered']['n']}** ({summ['recovered']['rate']})",
        f"- Survived (recovered and passing the verdict): **{summ['survived']['lower']}** lower / "
        f"**{summ['survived']['upper']}** upper bound (needs_lc rows: {summ['survived']['needs_lc']})",
        f"- Completeness at n_transits >= {MIN_TRANSITS}: n = {summ['completeness_min_transits']['n']}, "
        f"expected S >= S*: {summ['completeness_min_transits']['expected_s_ge_sstar']}, "
        f"recovered {summ['completeness_min_transits']['recovered_rate']}, "
        f"survived {summ['completeness_min_transits']['survived_rate_lower']} lower / "
        f"{summ['completeness_min_transits']['survived_rate_upper']} upper",
        "",
        "## Per-rule failures (ok injections)",
        "",
        "| rule | first-failing | any-failing | any-failing (upper bound) |",
        "|---|---|---|---|",
    ]
    rules = sorted(set(per_rule_first) | set(per_rule_any) | set(per_rule_any_up))
    for r in rules:
        lines.append(f"| {r} | {per_rule_first.get(r, 0)} | {per_rule_any.get(r, 0)} | {per_rule_any_up.get(r, 0)} |")
    lines += ["", "## Retention by injected T14 (dur_inj)", "",
              "| T14 | n | expected S>=S* | recovered | recovered rate | survived lower | survived upper |",
              "|---|---|---|---|---|---|---|"]
    for name, s in summ["retention_by_t14"].items():
        lines.append(f"| {name} | {s['n']} | {s['expected_s_ge_sstar']} | {s['recovered']} | "
                     f"{s['recovered_rate']} | {s['survived_lower']} | {s['survived_upper']} |")
    lines += ["", f"needs_lc rows anywhere in the table: {summ['per_rule']['n_needs_lc']} "
              f"(lower bound counts them as rejected, upper bound as exempted).", ""]
    (out / "summary.md").write_text("\n".join(lines))
    print(json.dumps(summ, indent=1, default=str))


if __name__ == "__main__":
    main()
