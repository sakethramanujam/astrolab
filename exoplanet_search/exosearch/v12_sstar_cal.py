"""v1.2 S*_cal calibration driver over S59-S64 (exosearch.v12_sstar_cal).

Assembles inputs and invokes the registered functions at commit bc05650 per rulings.md.
No statistics, vetting, holdout, eligibility, or S* logic is reimplemented here.
Run as: python -m exosearch.v12_sstar_cal --rows-dir ... --manifest ... --out-dir ... (--dry-run | --compute)
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import re
import sys
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence, Set, Tuple, Union

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

from exosearch import config, holdout, interim, quality, sstar_grid, v12_settings, vet, vet_final
from exosearch.v11 import known_planets, registration_checks

CALIBRATION_SECTORS: Tuple[int, ...] = (59, 60, 61, 62, 63, 64)
IDENTITY_COLS: Tuple[str, ...] = ("key", "tic", "sector", "mode", "holdout", "status")
KOSTOV_RELPATH = Path("kostov2025") / "kostov2025_tics_sectors_59-64_69_70.csv"
BOOT_N = 2000
BOOT_SEED = 20261110


class ExclusionMismatch(RuntimeError):
    """Rebuilt manifest exclusion differs from the registered TOI U Kostov 2025 set (counts only)."""

    def __init__(self, report: Dict[str, Any]) -> None:
        self.report = report
        super().__init__(f"exclusion mismatch (counts only): {json.dumps(report, sort_keys=True)}")


class MissingInput(FileNotFoundError):
    """Required calibration input (e.g. star-quality table or catalogue) is missing (counts only)."""


class NeedsLightCurve(RuntimeError):
    """One or more null rows require a light-curve rebuild not performed by this driver."""

    def __init__(self, counts: Dict[str, int]) -> None:
        self.counts = counts
        super().__init__(
            f"light-curve rebuild required for null rows (counts only): {json.dumps(counts, sort_keys=True)}"
        )


def load_manifest(path: Union[str, Path]) -> Tuple[Dict[str, Any], List[Dict[str, Any]]]:
    """Read a gzipped JSONL manifest (header line + job lines); mirrors exosearch.v11.run.read_manifest."""
    with gzip.open(path, "rt", encoding="utf-8") as fh:
        head_line = fh.readline()
        if not head_line.strip():
            raise MissingInput("manifest is empty")
        header = json.loads(head_line)["header"]
        jobs = [json.loads(line) for line in fh if line.strip()]
    return header, jobs


def _jobs_list(manifest: Union[str, Path, Sequence[Mapping[str, Any]]]) -> Sequence[Mapping[str, Any]]:
    if isinstance(manifest, (str, Path)):
        return load_manifest(manifest)[1]
    return manifest


def _part_files(rows_dir: Union[str, Path]) -> List[Path]:
    parts = sorted((Path(rows_dir) / "parts").glob("*.parquet"))
    if not parts:
        raise MissingInput("no parquet part files found under rows_dir/parts")
    return parts


def load_identity(rows_dir: Union[str, Path], sectors: Sequence[int] = CALIBRATION_SECTORS) -> pd.DataFrame:
    """Identity pass: read only IDENTITY_COLS for S59-S64 and verify stored holdout == holdout.is_holdout(tic)."""
    sec_list = [int(s) for s in sectors]
    sec_set = set(sec_list)
    frames: List[pd.DataFrame] = []
    for p in _part_files(rows_dir):
        tbl = pq.read_table(p, columns=list(IDENTITY_COLS), filters=[("sector", "in", sec_list)])
        if tbl.num_rows:
            frames.append(tbl.to_pandas())
    df = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(columns=list(IDENTITY_COLS))
    if len(df):
        got_secs = set(df["sector"].astype(int))
        if not got_secs.issubset(sec_set) or 69 in got_secs:
            raise RuntimeError(f"identity pass contained unexpected sectors: {sorted(got_secs)}")
        unique_tics = df["tic"].astype(int).unique()
        expected_map = {int(t): holdout.is_holdout(int(t)) for t in unique_tics}
        expected = df["tic"].astype(int).map(expected_map)
        n_bad = int((df["holdout"].astype(bool) != expected).sum())
        if n_bad:
            raise RuntimeError(f"stored holdout disagreed with exosearch.holdout.is_holdout on {n_bad} rows")
    return df


def load_science_for_keys(
    rows_dir: Union[str, Path], keys: Sequence[str], sectors: Sequence[int] = CALIBRATION_SECTORS
) -> pd.DataFrame:
    """Load science columns strictly for the final null-sample keys and guard against S69 / holdout / non-nulls."""
    key_list = sorted(set(keys))
    if not key_list:
        raise MissingInput("null-sample key list is empty")
    frames: List[pd.DataFrame] = []
    for p in _part_files(rows_dir):
        tbl = pq.read_table(p, filters=[("key", "in", key_list)])
        if tbl.num_rows:
            frames.append(tbl.to_pandas())
    if not frames:
        raise MissingInput("no science rows matched the null-sample keys")
    sci = pd.concat(frames, ignore_index=True)
    sci = sci[sci["status"] == "ok"].drop_duplicates("key").reset_index(drop=True)
    got_secs = set(sci["sector"].astype(int))
    if 69 in got_secs or not got_secs.issubset({int(s) for s in sectors}):
        raise RuntimeError(f"science guard failed: unexpected sectors {sorted(got_secs)}")
    if sci["holdout"].astype(bool).any():
        raise RuntimeError(f"science guard failed: {int(sci['holdout'].astype(bool).sum())} holdout rows present")
    got_modes = set(sci["mode"].astype(str))
    if not got_modes.issubset({"invert", "scramble"}):
        raise RuntimeError(f"science guard failed: unexpected modes {sorted(got_modes)}")
    return sci


def rebuild_exclusion(
    manifest: Union[str, Path, Sequence[Mapping[str, Any]]], sectors: Sequence[int] = CALIBRATION_SECTORS
) -> Dict[int, Set[int]]:
    """Step 1 (rulings.md Sec. 4): per sector, TICs with an invert task minus TICs with a scramble task."""
    jobs = _jobs_list(manifest)
    sec_set = {int(s) for s in sectors}
    inv: Dict[int, Set[int]] = {s: set() for s in sec_set}
    scr: Dict[int, Set[int]] = {s: set() for s in sec_set}
    for j in jobs:
        s = int(j["sector"])
        if s not in sec_set:
            continue
        tic = int(j["tic"])
        modes = {t.get("mode") for t in j.get("tasks", ())}
        if "invert" in modes:
            inv[s].add(tic)
        if "scramble" in modes:
            scr[s].add(tic)
    return {s: inv[s] - scr[s] for s in sorted(sec_set)}


def _parse_sectors_cell(val: Any) -> Set[int]:
    if val is None or (isinstance(val, float) and np.isnan(val)):
        return set()
    out: Set[int] = set()
    for tok in re.split(r"[;,\s]+", str(val).strip()):
        if not tok:
            continue
        try:
            out.add(int(float(tok)))
        except ValueError:
            continue
    return out


def registered_exclusion(
    manifest: Union[str, Path, Sequence[Mapping[str, Any]]], sectors: Sequence[int] = CALIBRATION_SECTORS
) -> Tuple[Dict[int, Set[int]], Dict[str, Dict[str, int]]]:
    """Step 2 (rulings.md Sec. 4; known_planets.load_toi_table): TOI (any disp) U Kostov 2025 per sector."""
    jobs = _jobs_list(manifest)
    sec_set = {int(s) for s in sectors}
    inv: Dict[int, Set[int]] = {s: set() for s in sec_set}
    for j in jobs:
        s = int(j["sector"])
        if s in sec_set and any(t.get("mode") == "invert" for t in j.get("tasks", ())):
            inv[s].add(int(j["tic"]))

    toi_df, _ = known_planets.load_toi_table()
    toi_tics = set(toi_df["tic"].astype(int))

    kostov_path = config.CATALOG_DIR / KOSTOV_RELPATH
    if not kostov_path.is_file():
        raise MissingInput(f"missing Kostov 2025 exclusion table: {KOSTOV_RELPATH}")
    kdf = pd.read_csv(kostov_path)
    kostov_by_sec: Dict[int, Set[int]] = {s: set() for s in sec_set}
    for tic_val, sec_cell in zip(kdf["tic"], kdf["sectors_in_59-64_69_70"]):
        if pd.isna(tic_val):
            continue
        tic = int(tic_val)
        for s in _parse_sectors_cell(sec_cell) & sec_set:
            kostov_by_sec[s].add(tic)

    sets: Dict[int, Set[int]] = {}
    counts: Dict[str, Dict[str, int]] = {}
    for s in sorted(sec_set):
        in_toi = inv[s] & toi_tics
        in_kostov = inv[s] & kostov_by_sec[s]
        union = in_toi | in_kostov
        sets[s] = union
        counts[str(s)] = {"n_toi": len(in_toi), "n_kostov": len(in_kostov), "n_union": len(union)}
    return sets, counts


def check_exclusion(
    rebuilt: Mapping[int, Set[int]], registered: Mapping[int, Set[int]]
) -> Dict[str, Any]:
    """Step 3 (rulings.md Sec. 4): verify rebuilt == registered; return sets + sha256 or raise counts-only error."""
    per_sector: Dict[str, Dict[str, int]] = {}
    mismatch = False
    for s in sorted(set(rebuilt) | set(registered)):
        rb = set(rebuilt.get(s, ()))
        rg = set(registered.get(s, ()))
        miss = rg - rb
        extra = rb - rg
        if miss or extra:
            mismatch = True
        per_sector[str(int(s))] = {
            "rebuilt": len(rb),
            "registered": len(rg),
            "missing_from_rebuilt": len(miss),
            "extra_in_rebuilt": len(extra),
        }
    if mismatch:
        raise ExclusionMismatch({"status": "mismatch", "per_sector": per_sector})
    canonical = "".join(sorted(f"{int(s)},{int(t)}\n" for s, tics in rebuilt.items() for t in tics)).encode("utf-8")
    digest = hashlib.sha256(canonical).hexdigest()
    label = "reconstructed after the run"
    return {
        "exclusion": {int(s): set(rebuilt[s]) for s in sorted(rebuilt)},
        "sha256": digest,
        "label": label,
        "report": {"status": "ok", "label": label, "sha256": digest, "per_sector": per_sector},
    }


def quality_drop(sectors: Sequence[int] = CALIBRATION_SECTORS) -> pd.DataFrame:
    """Step 4 (exosearch.quality.load; interim.with_quality semantics): load star-quality tables for all sectors."""
    sec_list = [int(s) for s in sectors]
    qdir = quality.out_dir()
    missing = [s for s in sec_list if not (qdir / f"star_quality_s{s:04d}.parquet").is_file()]
    if missing:
        raise MissingInput(f"missing star-quality tables for {len(missing)} sector(s): {missing}")
    q = quality.load(sectors=sec_list)
    present = set(q["sector"].astype(int)) if len(q) else set()
    empty_secs = [s for s in sec_list if s not in present]
    if empty_secs:
        raise MissingInput(f"empty star-quality tables for {len(empty_secs)} sector(s): {empty_secs}")
    return q


def alpha_n(identity: pd.DataFrame, q: pd.DataFrame) -> Dict[int, Dict[str, Any]]:
    """Step 5 (rulings.md Sec. 1; beta_rows.beta_frames n_alpha): N_s = quality-passing invert rows, p_s = 1/N_s."""
    m = interim.with_quality(identity, q)
    out: Dict[int, Dict[str, Any]] = {}
    for s in sorted(int(x) for x in identity["sector"].unique()):
        sub = m[(m["sector"] == s) & (m["mode"] == "invert")].drop_duplicates("tic")  # as beta_frames n_alpha (no status filter)
        n_s = int((~sub["drop"]).sum())
        if n_s <= 0:
            raise MissingInput(f"sector {s} has zero quality-passing invert rows")
        out[s] = {"N_s": n_s, "p_s": 1.0 / n_s}
    return out


def null_sample(
    identity: pd.DataFrame, q: pd.DataFrame, exclusion: Mapping[int, Set[int]]
) -> Tuple[List[str], Dict[str, Any]]:
    """Step 6 (rulings.md Sec. 4; registration_checks.check_registration_population): paired null keys + guard."""
    m = interim.with_quality(identity, q)
    sel = m[(~m["drop"]) & (~m["holdout"].astype(bool)) & (m["status"] == "ok")]
    keys: List[str] = []
    per_sector: Dict[str, Dict[str, Any]] = {}
    for s in sorted(int(x) for x in identity["sector"].unique()):
        sub = sel[sel["sector"] == s]
        real_df = sub[sub["mode"] == "real"]
        inv_df = sub[sub["mode"] == "invert"]
        scr_df = sub[sub["mode"] == "scramble"]
        ex_df = pd.DataFrame({"tic": sorted(exclusion[s])})
        guard = registration_checks.check_registration_population(inv_df, scr_df, ex_df, sector=s)
        common = set(real_df["tic"].astype(int)) & set(inv_df["tic"].astype(int)) & set(scr_df["tic"].astype(int))
        if common & set(exclusion[s]):
            raise registration_checks.PopulationGuardError(
                f"sector {s}: {len(common & set(exclusion[s]))} excluded star(s) in paired null sample"
            )
        inv_keys = sorted(inv_df[inv_df["tic"].astype(int).isin(common)]["key"].astype(str).unique())
        scr_keys = sorted(scr_df[scr_df["tic"].astype(int).isin(common)]["key"].astype(str).unique())
        keys.extend(inv_keys)
        keys.extend(scr_keys)
        per_sector[str(s)] = {
            "n_stars": len(common),
            "n_invert_rows": len(inv_keys),
            "n_scramble_rows": len(scr_keys),
            "n_null_rows": len(inv_keys) + len(scr_keys),
            "guard": guard,
        }
    return sorted(set(keys)), {"total_null_rows": len(set(keys)), "per_sector": per_sector}


def vet_nulls(
    sci: pd.DataFrame, manifest: Union[str, Path, Sequence[Mapping[str, Any]]], cfg: Mapping[str, Any]
) -> pd.DataFrame:
    """Step 7 (rulings.md Sec. 3; beta_rows.with_vetting + vet_final.final_verdict): compute passed_lc_checks_v12."""
    jobs = _jobs_list(manifest)
    job_map = {(int(j["sector"]), int(j["tic"])): j for j in jobs}
    group_col = "task" if "task" in sci.columns else "mode"
    parts: List[pd.DataFrame] = []
    for _, g in sci.groupby(["sector", group_col], sort=True):
        gf = vet.flag_common_periods(g.assign(status="ok"))
        gf["fail_common_period"] = gf["flags"].fillna("").str.contains("common_period")
        parts.append(gf)
    df = pd.concat(parts, ignore_index=True)

    passed: List[bool] = []
    n_needs_lc = 0
    n_odd_even_needs_lc = 0
    for r in df.to_dict(orient="records"):
        job = job_map.get((int(r["sector"]), int(r["tic"])), {})
        verdict = vet_final.final_verdict(
            r,
            common_period=bool(r["fail_common_period"]),
            checks_logged=r.get("checks_logged"),
            rho=job.get("rho"),
            cfg=dict(cfg),
            rstar=job.get("rad_tic"),
            e_rad=job.get("e_rad"),
            needs_lc="reject",
        )
        passed.append(bool(verdict.passed))
        if verdict.status == vet_final.STATUS_NEEDS_LC:
            n_needs_lc += 1
        if vet_final.ODD_EVEN_NEEDS_LC in verdict.all_failing:
            n_odd_even_needs_lc += 1
    if n_needs_lc > 0 or n_odd_even_needs_lc > 0:
        raise NeedsLightCurve({"needs_lc": n_needs_lc, "odd_even_needs_lc": n_odd_even_needs_lc})
    df["passed_lc_checks_v12"] = passed
    return df


def compute_sstar(sci: pd.DataFrame, alpha: Mapping[int, Mapping[str, Any]]) -> Dict[str, Any]:
    """Step 8 (rulings.md Sec. 1-3; v12_settings.calibrate_sstar_v12; sstar_grid.sstar_grid): pooled S* + diagnostics."""
    sectors = sorted(int(s) for s in sci["sector"].unique())
    total_null = int(len(sci))
    budget = sum(float(alpha[s]["p_s"]) * int((sci["sector"] == s).sum()) for s in sectors)
    S_all = sci["S"].to_numpy(dtype=float)
    vetted = sci["passed_lc_checks_v12"].astype(bool).to_numpy()
    pooled = v12_settings.calibrate_sstar_v12(S_all[vetted], n_null=total_null, p=budget / total_null, u=None)

    per_sector: Dict[str, Any] = {}
    for s in sectors:
        sub = sci[sci["sector"] == s]
        S_s = sub["S"].to_numpy(dtype=float)
        v_s = sub["passed_lc_checks_v12"].astype(bool).to_numpy()
        n_s = int(len(sub))
        b_s = float(alpha[s]["p_s"]) * n_s
        pt = sstar_grid.sstar_grid(S_s[v_s], b_s)
        rng = np.random.default_rng(BOOT_SEED)
        boots: List[float] = []
        n_fail = 0
        for _ in range(BOOT_N):
            idx = rng.integers(0, n_s, n_s)
            try:
                boots.append(float(sstar_grid.sstar_grid(S_s[idx][v_s[idx]], b_s)["s_star"]))
            except sstar_grid.SStarGridError:
                n_fail += 1
        ci68 = [float(np.percentile(boots, 16)), float(np.percentile(boots, 84))] if boots else None
        per_sector[str(s)] = {
            **pt,
            "n_null": n_s,
            "p_s": float(alpha[s]["p_s"]),
            "boot_n": BOOT_N,
            "boot_seed": BOOT_SEED,
            "boot_fail": n_fail,
            "ci68": ci68,
            "report_only": True,
        }

    unvetted = v12_settings.calibrate_sstar_v12(S_all, n_null=total_null, p=budget / total_null, u=None)
    by_mode: Dict[str, Dict[str, Any]] = {}
    for mode in ("invert", "scramble"):
        sub_m = sci[sci["mode"] == mode]
        n_m = int(len(sub_m))
        b_m = sum(float(alpha[s]["p_s"]) * int((sub_m["sector"] == s).sum()) for s in sectors)
        S_m = sub_m.loc[sub_m["passed_lc_checks_v12"].astype(bool), "S"].to_numpy(dtype=float)
        by_mode[mode] = v12_settings.calibrate_sstar_v12(S_m, n_null=n_m, p=b_m / n_m, u=None)

    return {
        "pooled": pooled,
        "per_sector_report_only": per_sector,
        "diagnostics": {
            "unvetted": unvetted,
            "invert_only": by_mode["invert"],
            "scramble_only": by_mode["scramble"],
            "max_invert_scramble_s_star": max(float(by_mode["invert"]["s_star"]), float(by_mode["scramble"]["s_star"])),
        },
    }


def n_real_eligible(
    s_star: float,
    identity: pd.DataFrame,
    manifest: Union[str, Path, Sequence[Mapping[str, Any]]],
    out_dir: Union[str, Path],
) -> Dict[str, Any]:
    """Step 9 (rulings.md Sec. 5; known_planets.freeze_eligibility_v12): count unique eligible planets in [S*, 2S*]."""
    jobs = _jobs_list(manifest)
    out_p = Path(out_dir)
    out_p.mkdir(parents=True, exist_ok=True)
    fname_map = {(int(j["sector"]), int(j["tic"])): str(j["fname"]) for j in jobs}
    all_planets: Set[Tuple[int, str]] = set()
    per_sector: Dict[str, Dict[str, Any]] = {}

    for s in sorted(int(x) for x in identity["sector"].unique()):
        real_s = identity[
            (identity["sector"] == s)
            & (identity["mode"] == "real")
            & (~identity["holdout"].astype(bool))
            & (identity["status"] == "ok")
        ]
        searched_tics = set(real_s["tic"].astype(int))

        def _lc_lookup(tic: int, sec: int = s) -> Optional[Path]:
            fn = fname_map.get((sec, int(tic)))
            return (config.lc_dir(sec) / fn) if fn else None

        df, sidecar = known_planets.freeze_eligibility_v12(
            s_star=float(s_star),
            out_file=out_p / f"eligibility_s{s:04d}.csv",
            sector=s,
            searched_tics=searched_tics,
            lc_lookup=_lc_lookup,
        )
        type_ok = (df["type"] == "ps") | ((df["type"] == "toi") & df["disp"].isin(["CP", "KP"]))
        s_ok = df["expected_S"].notna() & (df["expected_S"] >= float(s_star)) & (df["expected_S"] <= 2.0 * float(s_star))
        ntr_ok = df["n_tr"] >= 3
        clean_ok = (df["base_ineligible_reason"].fillna("") == "") & (df["error"].fillna("") == "")
        non_holdout = ~df["tic"].astype(int).map(holdout.is_holdout)
        elig = df[type_ok & s_ok & ntr_ok & clean_ok & non_holdout]
        sec_planets = {(int(r.tic), str(r.name)) for r in elig.itertuples()}
        all_planets |= sec_planets
        per_sector[str(s)] = {
            "n_searched_real_non_holdout": len(searched_tics),
            "n_total_targets": int(len(df)),
            "n_eligible_unique_planets": len(sec_planets),
            "eligibility_csv_sha256": sidecar["sha256"],
        }
    return {"n_real_eligible": len(all_planets), "per_sector": per_sector}


def _identity_counts(identity: pd.DataFrame) -> Tuple[Dict[str, Dict[str, int]], Dict[str, Any]]:
    by_sec_mode: Dict[str, Dict[str, int]] = {}
    holdout_by_sec: Dict[str, int] = {}
    for s in sorted(int(x) for x in identity["sector"].unique()):
        sub = identity[identity["sector"] == s]
        by_sec_mode[str(s)] = {str(m): int(c) for m, c in sub["mode"].value_counts().sort_index().items()}
        holdout_by_sec[str(s)] = int(sub["holdout"].astype(bool).sum())
    return by_sec_mode, {"total": sum(holdout_by_sec.values()), "per_sector": holdout_by_sec}


def main(argv: Optional[Sequence[str]] = None) -> int:
    """Step 10: CLI entry point (--dry-run: steps 1-6 counts only; --compute: full S*_cal + sha256 sidecar)."""
    ap = argparse.ArgumentParser(description="v1.2 S*_cal driver over S59-S64 (rulings.md)")
    ap.add_argument("--rows-dir", required=True, type=Path, help="Directory containing parts/*.parquet")
    ap.add_argument("--manifest", required=True, type=Path, help="Gzipped JSONL run manifest")
    ap.add_argument("--out-dir", required=True, type=Path, help="Output directory")
    mode = ap.add_mutually_exclusive_group(required=True)
    mode.add_argument("--dry-run", action="store_true", help="Run steps 1-6 only (no science columns, vetting, or S*)")
    mode.add_argument("--compute", action="store_true", help="Run full calibration and write sstar_cal_v12.json + .sha256")
    args = ap.parse_args(argv)

    identity = load_identity(args.rows_dir, CALIBRATION_SECTORS)
    by_sec_mode, holdout_counts = _identity_counts(identity)
    _, jobs = load_manifest(args.manifest)

    rebuilt = rebuild_exclusion(jobs, CALIBRATION_SECTORS)
    reg_sets, reg_counts = registered_exclusion(jobs, CALIBRATION_SECTORS)

    base_report: Dict[str, Any] = {
        "sectors": list(CALIBRATION_SECTORS),
        "identity_rows_by_sector_mode": by_sec_mode,
        "holdout_rows_skipped": holdout_counts,
        "registered_exclusion_counts": reg_counts,
        "notes": ["rows carry scramble_exclude_file_70 only, not an S59-S64 stamp (rulings.md Sec. 4)"],
    }

    try:
        ex_res = check_exclusion(rebuilt, reg_sets)
    except ExclusionMismatch as exc:
        base_report["exclusion_check"] = exc.report
        print(json.dumps(base_report, indent=2, sort_keys=True))
        return 2
    base_report["exclusion_check"] = ex_res["report"]

    if args.dry_run:
        try:
            q = quality_drop(CALIBRATION_SECTORS)
            alpha = alpha_n(identity, q)
            _, null_rep = null_sample(identity, q, ex_res["exclusion"])
            base_report["alpha"] = {str(s): alpha[s] for s in sorted(alpha)}
            base_report["null_sample"] = null_rep
        except MissingInput as exc:
            base_report["quality_error"] = str(exc)
            base_report["alpha"] = None
            base_report["null_sample"] = None
        print(json.dumps(base_report, indent=2, sort_keys=True))
        return 0

    q = quality_drop(CALIBRATION_SECTORS)
    alpha = alpha_n(identity, q)
    keys, null_rep = null_sample(identity, q, ex_res["exclusion"])
    cfg = v12_settings.require_v12_settings()
    sci = load_science_for_keys(args.rows_dir, keys, CALIBRATION_SECTORS)
    sci = vet_nulls(sci, jobs, cfg)
    sstar_res = compute_sstar(sci, alpha)
    elig_res = n_real_eligible(float(sstar_res["pooled"]["s_star"]), identity, jobs, args.out_dir)

    result = {
        **base_report,
        "alpha": {str(s): alpha[s] for s in sorted(alpha)},
        "null_sample": null_rep,
        "sstar": sstar_res,
        "n_real_eligible": elig_res,
        "v12_settings_hash": sstar_res["pooled"]["v12_settings_hash"],
        "v12_cap_rearth": sstar_res["pooled"]["v12_cap_rearth"],
    }
    args.out_dir.mkdir(parents=True, exist_ok=True)
    out_json = args.out_dir / "sstar_cal_v12.json"
    payload = json.dumps(result, indent=2, sort_keys=True) + "\n"
    out_json.write_text(payload, encoding="utf-8")
    digest = hashlib.sha256(out_json.read_bytes()).hexdigest()
    (args.out_dir / "sstar_cal_v12.json.sha256").write_text(f"{digest}  {out_json.name}\n", encoding="utf-8")
    print(payload, end="")
    return 0


if __name__ == "__main__":
    sys.exit(main())
