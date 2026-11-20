"""v1.2 S*_cal calibration driver over S59-S64 (exosearch.v12_sstar_cal).

Assembles inputs and invokes pinned project functions per docs/registration_v1.2.md,
S*_cal rulings and Deviation 2026-11-20. No search, S statistic, vetting rules, beta,
holdout, eligibility, or S* selection logic is reimplemented here.

Run as:
  python -m exosearch.v12_sstar_cal --rows-dir ... --manifest ... --out-dir ... \
      --quality-dir ... --input-hashes ... --common-period-mode ... \
      [--workers 16] [--gate-top 300] (--dry-run | --gate | --compute)

Exit codes:  # review: S7
  0: ok
  2: population / three-way mismatch
  3: gate fail
  4: gate depth exceeded
  5: needs light curve
  6: light-curve rebuild error
  7: secondary-scan consistency failure
  8: catalogue hash mismatch
  9: missing p_grid_halfwidth table
"""
from __future__ import annotations

import os

for _THREAD_ENV in (  # review: B2
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "NUMBA_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
):
    os.environ[_THREAD_ENV] = "1"
# Note: no exosearch code at the pinned checkout (bc05650) reads EXO_OFFLINE.
# The socket monkeypatch below blocks Python-level socket creation.
os.environ["EXO_OFFLINE"] = "1"

import socket as _socket_mod

_ORIG_SOCKET = _socket_mod.socket
_ORIG_CREATE_CONNECTION = _socket_mod.create_connection
_OFFLINE_ERR = "v12_sstar_cal driver is offline and blocks Python-level socket creation"


class _BlockedSocket(_ORIG_SOCKET):
    """Replacement for socket.socket: subclassing remains valid, construction raises."""

    def __new__(cls, *args: object, **kwargs: object):  # type: ignore[override]
        raise OSError(_OFFLINE_ERR)


def _blocked_create_connection(*args: object, **kwargs: object):
    raise OSError(_OFFLINE_ERR)


_socket_mod.socket = _BlockedSocket  # type: ignore[misc,assignment]
_socket_mod.create_connection = _blocked_create_connection  # type: ignore[assignment]

import argparse
import concurrent.futures
from datetime import datetime, timezone
import gzip
import hashlib
import importlib.util
import json
import math
import multiprocessing
import re
import sys
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Set, Tuple, Union

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

import exosearch
from exosearch import (
    batch,
    config,
    holdout,
    interim,
    nulltest,
    p_grid_halfwidth,
    quality,
    sstar_grid,
    v12_settings,
    vet,
    vet_final,
    worker_guard,
)
from exosearch.v11 import kernels, known_planets, pipeline, registration_checks, scramble, settings, vet_amend

CALIBRATION_SECTORS: Tuple[int, ...] = (59, 60, 61, 62, 63, 64)
IDENTITY_COLS: Tuple[str, ...] = ("key", "tic", "sector", "task", "mode", "holdout", "status")
COMMON_PERIOD_COLS: Tuple[str, ...] = ("key", "sector", "task", "mode", "holdout", "status", "period", "flags")
COMMON_PERIOD_MODES: Tuple[str, ...] = (
    "v11_task_with_holdout",
    "v11_task_no_holdout",
    "null_sample",
)
RULED_COMMON_PERIOD_MODE = "v11_task_with_holdout"
NULL_MODES: Tuple[str, ...] = ("invert", "scramble")
ALL_SEARCH_MODES: Tuple[str, ...] = ("real", "invert", "scramble")
COMMON_PERIOD_POPULATION_DESC: Dict[str, str] = {
    "v11_task_with_holdout": (
        "ok null-mode rows (invert and scramble) per (sector, task), after EB drop, "
        "quality-dropped and excluded-star rows included; includes holdout null periods "
        "for neighbour counts only (summaries never aggregate holdout or real rows)"
    ),
    "v11_task_no_holdout": (
        "ok non-holdout null-mode rows (invert and scramble) per (sector, task), after EB drop, "
        "quality-dropped and excluded-star rows included (summaries never aggregate holdout or real rows)"
    ),
    "null_sample": "null-sample rows only per (sector, task)",
}
KOSTOV_RELPATH = Path("kostov2025") / "kostov2025_tics_sectors_59-64_69_70.csv"
FROZEN_TOI_RELPATH = Path("frozen") / "exofop_toi_3df0ecc.csv"
EXPECTED_EB_DROP: Dict[int, int] = {59: 7, 60: 5, 61: 6, 62: 6, 63: 3, 64: 2}
BOOT_N = 2000
BOOT_SEED = 20261110
GATE_RANDOM_N = 200
GATE_RANDOM_SEED = 20261120
GATE_S_TOL = 1e-6
GATE_S_COLS: Tuple[str, ...] = ("S", "S_opt_a", "S_opt_b", "S_opt_c", "S_opt_d")
GATE_EXACT_COLS: Tuple[str, ...] = ("period", "P_grid", "duration", "t0")
V12_ALLOWED_CFG_DIFF_KEYS: Set[str] = {
    "name",
    "needs_lc.rebuild",
    "sstar.fine_grid_enabled",
    "vet.odd_even_mask.enabled",
    "vet.rp_cap.enabled",
    "vet.rp_cap.include_e_rad",
    "vet.too_deep.mode",
}
_BUDGET_EPS = 1e-9
PINNED_CODE_COMMIT = "bc05650d69f739fd5fea2f37629e52c5d524b3a8"
PINNED_INPUT_HASHES_SHA256 = "20484217cf60e4c9a027309ba843ed67264c2d73b9b1245009004f581d8fc1db"
PINNED_PSCOMPPARS_SHA256 = "6dc56f25bad4516a80c1b000afbaa4f02b094bb2da2476187a677a27b1f78c3e"  # review: S2
PINNED_TOI_SHA256 = "3df0eccdf26fba9a3d4bc6bf4d0a3f53c0662e6b3dc1ef7827a74a519ea34ce0"
PINNED_KOSTOV_SHA256 = "2ce6e7cc1e539845ab723e75d7d8115eea9b38fe09c82f97f04fb57a258aeff3"
PINNED_PGRID_SHA256 = "5904f4f8019a9d3134998aa15c5187e878990ebf48f69e0b779ddb68bd9e19c3"
PINNED_QUALITY_SHA256: Dict[int, str] = {  # review: S5
    59: "7f97d1d0dd86140a76d8cc654fe9e5b996ee208bfcdc39434077f2f718433f7a",
    60: "49ab5f15739839355261e8a04e3d5ba80f0c9bd919d9a456a280f16297f93cfd",
    61: "5a4a628f3cbcb3527a3b4fe2b492558dd98da51e3d4fda0cac3d918c449303b7",
    62: "c02606bfcbc5724e78ece441c853b5af3cef78958cfe31966809af2076cc34e6",
    63: "1250055fad23b47f93d6112d5cdf41be6cdb67551321e9e88672c1ff0ae708e1",
    64: "92ef33f2204b03199c5b856c9278695efa8538487dcf1615f7c1d4745dfe4673",
}
VETTED_CSV_COLS: Tuple[str, ...] = (
    "key",
    "sector",
    "S",
    "passed",
    "first_failing_rule",
    "all_failing",
    "status",
    "too_deep_v2",
)
_WORKER_FINISH: Optional[Any] = None


class ExclusionMismatch(RuntimeError):
    """Exclusion / EB deviation check failed (counts only)."""

    def __init__(self, report: Dict[str, Any]) -> None:
        self.report = report
        super().__init__(f"exclusion mismatch (counts only): {json.dumps(report, sort_keys=True)}")


class ThreeWayMismatch(RuntimeError):
    """Structural three-way match violation on post-EB-drop identity rows (counts only)."""

    def __init__(self, report: Dict[str, Any]) -> None:
        self.report = report
        super().__init__(f"three-way mismatch (counts only): {json.dumps(report, sort_keys=True)}")


class MissingInput(FileNotFoundError):
    """Required calibration input (e.g. star-quality table or catalogue) is missing (counts only)."""


class NeedsLightCurve(RuntimeError):
    """A null row still requires a light curve during all-hard vetting (counts only)."""

    def __init__(self, report: Dict[str, Any]) -> None:
        self.report = report
        super().__init__(f"light-curve rebuild required (counts only): {json.dumps(report, sort_keys=True)}")


class GateDepthExceeded(RuntimeError):
    """Vetting walk reached a row outside the reproduced gate key set (counts only)."""

    def __init__(self, report: Dict[str, Any]) -> None:
        self.report = report
        super().__init__(f"gate depth exceeded (counts only): {json.dumps(report, sort_keys=True)}")


class RebuildError(RuntimeError):  # review: S7
    """Light-curve rebuild failed during all-hard vetting (counts only)."""

    def __init__(self, report: Dict[str, Any]) -> None:
        self.report = report
        super().__init__(f"light-curve rebuild error (counts only): {json.dumps(report, sort_keys=True)}")


class ScanConsistencyError(RuntimeError):  # review: S7
    """Rebuilt secondary scan failed consistency check during vetting (counts only)."""

    def __init__(self, report: Dict[str, Any]) -> None:
        self.report = report
        super().__init__(f"secondary-scan consistency failure (counts only): {json.dumps(report, sort_keys=True)}")


# --------------------------------------------------------------------------- I/O helpers
def _sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _git_blob_sha1(data: bytes) -> str:
    """SHA-1 of a git blob, matching `git hash-object` / `git ls-tree`."""
    return hashlib.sha1(b"blob " + str(len(data)).encode("ascii") + b"\0" + data).hexdigest()


def _driver_sha256() -> str:  # review: B4
    return _sha256_bytes(Path(__file__).resolve().read_bytes())


def _driver_git_blob_sha1() -> str:
    return _git_blob_sha1(Path(__file__).resolve().read_bytes())


def _frozen_toi_path() -> Path:
    return config.CATALOG_DIR / FROZEN_TOI_RELPATH


def _p_grid_halfwidth_path() -> Path:
    """Path of the table written and read by exosearch.p_grid_halfwidth (OUT)."""
    return Path(p_grid_halfwidth.OUT)


def check_input_hashes(
    input_hashes_path: Union[str, Path],
    quality_dir: Union[str, Path],
    sectors: Sequence[int] = CALIBRATION_SECTORS,
) -> Tuple[Dict[int, str], Dict[str, Any]]:
    """Verify the pinned --input-hashes file, star-quality parquets in quality_dir, and any listed frozen pins."""
    manifest_path = Path(input_hashes_path)
    if not manifest_path.is_file():
        raise SystemExit(f"missing required --input-hashes file: {manifest_path.name}")

    manifest_bytes = manifest_path.read_bytes()
    manifest_sha = _sha256_bytes(manifest_bytes)
    if manifest_sha != PINNED_INPUT_HASHES_SHA256:
        raise SystemExit(
            f"--input-hashes sha256 mismatch: expected {PINNED_INPUT_HASHES_SHA256}, got {manifest_sha}"
        )

    hex64 = re.compile(r"^[0-9a-fA-F]{64}$")
    sec_re = re.compile(r"star_quality_s0*(\d{2})\.parquet$", re.IGNORECASE)
    quality_expected: Dict[int, str] = {}
    n_checked = 0
    mismatches: List[str] = []
    for raw_line in manifest_bytes.decode("utf-8").splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        tokens = line.split()
        if len(tokens) < 2 or not hex64.match(tokens[0]):
            raise SystemExit(f"malformed entry in --input-hashes file: {line!r}")
        got_hash = tokens[0].lower()
        entry_name = tokens[1].lstrip("*")
        base = Path(entry_name).name.lower()
        sec_m = sec_re.search(base)
        if sec_m:
            sec = int(sec_m.group(1))
            if sec in quality_expected and quality_expected[sec] != got_hash:
                raise SystemExit(f"conflicting --input-hashes entries for sector {sec}")
            quality_expected[sec] = got_hash
            if sec in PINNED_QUALITY_SHA256:
                n_checked += 1
                if got_hash != PINNED_QUALITY_SHA256[sec].lower():
                    mismatches.append(f"quality_s{sec}")
        elif "exofop_toi" in base or base == "toi":
            n_checked += 1
            if got_hash != PINNED_TOI_SHA256.lower():
                mismatches.append("toi")
        elif "pscomppars" in base:
            n_checked += 1
            if got_hash != PINNED_PSCOMPPARS_SHA256.lower():
                mismatches.append("pscomppars")
        elif "kostov" in base:
            n_checked += 1
            if got_hash != PINNED_KOSTOV_SHA256.lower():
                mismatches.append("kostov")
        elif "p_grid" in base or "pgrid" in base:
            n_checked += 1
            if got_hash != PINNED_PGRID_SHA256.lower():
                mismatches.append("p_grid_halfwidth")

    if mismatches:
        raise SystemExit(
            f"hash mismatch between --input-hashes file and pinned constants: {sorted(mismatches)}"
        )

    sec_list = [int(s) for s in sectors]
    missing_entries = [s for s in sec_list if s not in quality_expected]
    if missing_entries:
        raise SystemExit(f"--input-hashes missing star-quality entries for sector(s): {missing_entries}")

    qdir = Path(quality_dir)
    verified_quality: Dict[int, str] = {}
    for s in sec_list:
        fpath = qdir / f"star_quality_s{s:04d}.parquet"
        if not fpath.is_file():
            raise SystemExit(f"missing star-quality parquet for sector {s}: {fpath.name}")
        got = _sha256_bytes(fpath.read_bytes())
        exp = quality_expected[s]
        if got != exp:
            raise SystemExit(
                f"star-quality sha256 mismatch for sector {s}: expected {exp} from --input-hashes, got {got}"
            )
        verified_quality[s] = got

    report = {
        "filename": manifest_path.name,
        "sha256": manifest_sha,
        "expected_sha256": PINNED_INPUT_HASHES_SHA256,
        "status": "ok",
        "n_entries_checked": n_checked,
        "n_quality_verified": len(verified_quality),
        "n_mismatches": 0,
    }
    return verified_quality, report


def check_catalogue_hashes() -> Dict[str, Any]:  # review: S2
    """Verify sha256 of pinned pscomppars, frozen TOI, and Kostov 2025 catalogues."""
    # CTOI (exofop_ctoi.csv) and Prša et al. 2022 (tess_ebs_prsa2022.csv) are not
    # read on this driver path and are left unpinned.
    _, ps_path = known_planets.load_pscomppars()
    _, toi_path = known_planets.load_toi_table(_frozen_toi_path())
    kostov_path = config.CATALOG_DIR / KOSTOV_RELPATH
    if not kostov_path.is_file():
        raise MissingInput(f"missing Kostov 2025 exclusion table: {KOSTOV_RELPATH}")
    ps_got = _sha256_bytes(Path(ps_path).read_bytes())
    toi_got = _sha256_bytes(Path(toi_path).read_bytes())
    kostov_got = _sha256_bytes(kostov_path.read_bytes())
    ps_ok = ps_got == PINNED_PSCOMPPARS_SHA256
    toi_ok = toi_got == PINNED_TOI_SHA256
    kostov_ok = kostov_got == PINNED_KOSTOV_SHA256
    mismatches: List[str] = []
    if not ps_ok:
        mismatches.append("pscomppars")
    if not toi_ok:
        mismatches.append("toi")
    if not kostov_ok:
        mismatches.append("kostov")
    return {
        "status": "ok" if not mismatches else "mismatch",
        "mismatched_catalogues": mismatches,
        "pscomppars": {
            "expected_sha256": PINNED_PSCOMPPARS_SHA256,
            "got_sha256": ps_got,
            "ok": ps_ok,
        },
        "toi": {
            "relpath": (Path("catalogs") / FROZEN_TOI_RELPATH).as_posix(),
            "expected_sha256": PINNED_TOI_SHA256,
            "got_sha256": toi_got,
            "ok": toi_ok,
        },
        "kostov": {
            "relpath": (Path("catalogs") / KOSTOV_RELPATH).as_posix(),
            "expected_sha256": PINNED_KOSTOV_SHA256,
            "got_sha256": kostov_got,
            "ok": kostov_ok,
        },
    }


def check_p_grid_halfwidth(require: bool = True) -> Dict[str, Any]:  # review: S3
    """Verify presence and pinned sha256 of p_grid_halfwidth table; raise MissingInput if absent or mismatched."""
    present = p_grid_halfwidth.load() is not None
    if not present:
        if require:
            raise MissingInput("p_grid_halfwidth table is required but p_grid_halfwidth.load() returned None")
        return {
            "present": False,
            "ok": False,
            "expected_sha256": PINNED_PGRID_SHA256,
            "got_sha256": None,
        }
    pgh_path = _p_grid_halfwidth_path()
    if pgh_path is None or not pgh_path.is_file():
        if require:
            raise MissingInput("p_grid_halfwidth table file could not be located for sha256 verification")
        return {
            "present": False,
            "ok": False,
            "expected_sha256": PINNED_PGRID_SHA256,
            "got_sha256": None,
        }
    got = _sha256_bytes(pgh_path.read_bytes())
    ok = got == PINNED_PGRID_SHA256
    if require and not ok:
        raise MissingInput(
            f"p_grid_halfwidth sha256 mismatch: expected {PINNED_PGRID_SHA256}, got {got}"
        )
    return {
        "present": True,
        "ok": ok,
        "expected_sha256": PINNED_PGRID_SHA256,
        "got_sha256": got,
    }


def _write_json(path: Path, obj: Mapping[str, Any]) -> str:
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = json.dumps(obj, indent=2, sort_keys=True) + "\n"
    path.write_text(payload, encoding="utf-8")
    return payload


def _gate_keys_sha256(keys: Sequence[str]) -> str:
    return _sha256_bytes("\n".join(sorted(set(str(k) for k in keys))).encode("utf-8"))


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


def load_null_s_for_keys(rows_dir: Union[str, Path], keys: Sequence[str]) -> pd.DataFrame:
    """Load only columns (key, S) for already-validated non-holdout null-sample keys."""
    key_list = sorted(set(str(k) for k in keys))
    if not key_list:
        raise MissingInput("null-sample key list is empty")
    frames: List[pd.DataFrame] = []
    for p in _part_files(rows_dir):
        tbl = pq.read_table(p, columns=["key", "S"], filters=[("key", "in", key_list)])
        if tbl.num_rows:
            frames.append(tbl.to_pandas())
    if not frames:
        raise MissingInput("no S rows matched the null-sample keys")
    df = pd.concat(frames, ignore_index=True).drop_duplicates("key").reset_index(drop=True)
    if len(df) != len(key_list):
        raise MissingInput(f"missing {len(key_list) - len(df)} S row(s) for null-sample keys")
    return df


def load_common_period_rows_for_keys(
    rows_dir: Union[str, Path],
    keys: Sequence[str],
    identity: pd.DataFrame,
    sectors: Sequence[int] = CALIBRATION_SECTORS,
    allowed_modes: Iterable[str] = NULL_MODES,
    allow_holdout: bool = False,
) -> pd.DataFrame:
    """Load only COMMON_PERIOD_COLS for requested keys with S69, holdout (via identity), and mode guards."""
    key_list = sorted(set(str(k) for k in keys))
    if not key_list:
        raise MissingInput("common-period key list is empty")
    sec_set = {int(s) for s in sectors}
    allowed_mode_set = {str(m) for m in allowed_modes}

    id_sub = identity[identity["key"].astype(str).isin(key_list)].drop_duplicates("key")
    if len(id_sub) != len(key_list):
        raise RuntimeError(
            f"common-period guard failed: {len(key_list) - len(id_sub)} key(s) missing from identity frame"
        )
    id_secs = set(id_sub["sector"].astype(int))
    if 69 in id_secs or not id_secs.issubset(sec_set):
        raise RuntimeError(f"common-period guard failed: unexpected identity sectors {sorted(id_secs)}")
    if (not allow_holdout) and id_sub["holdout"].astype(bool).any():
        raise RuntimeError(
            f"common-period guard failed: {int(id_sub['holdout'].astype(bool).sum())} holdout rows in identity"
        )
    id_modes = set(id_sub["mode"].astype(str))
    if not id_modes.issubset(allowed_mode_set):
        raise RuntimeError(f"common-period guard failed: unexpected identity modes {sorted(id_modes)}")

    frames: List[pd.DataFrame] = []
    for p in _part_files(rows_dir):
        tbl = pq.read_table(p, columns=list(COMMON_PERIOD_COLS), filters=[("key", "in", key_list)])
        if tbl.num_rows:
            frames.append(tbl.to_pandas())
    if not frames:
        raise MissingInput("no common-period rows matched the requested keys")
    df = pd.concat(frames, ignore_index=True)
    df = df[df["status"].astype(str) == "ok"].drop_duplicates("key").reset_index(drop=True)
    got_secs = set(df["sector"].astype(int))
    if 69 in got_secs or not got_secs.issubset(sec_set):
        raise RuntimeError(f"common-period guard failed: unexpected sectors {sorted(got_secs)}")
    got_modes = set(df["mode"].astype(str))
    if not got_modes.issubset(allowed_mode_set):
        raise RuntimeError(f"common-period guard failed: unexpected modes {sorted(got_modes)}")
    if len(df) != len(key_list):
        raise MissingInput(f"missing {len(key_list) - len(df)} ok common-period row(s) for requested keys")
    return df


def load_science_for_keys(
    rows_dir: Union[str, Path],
    keys: Sequence[str],
    sectors: Sequence[int] = CALIBRATION_SECTORS,
    allowed_modes: Iterable[str] = NULL_MODES,
    allowed_real_keys: Optional[Set[str]] = None,
) -> pd.DataFrame:
    """Load science columns strictly for requested keys with S69, holdout, and mode guards."""
    key_list = sorted(set(str(k) for k in keys))
    if not key_list:
        raise MissingInput("science key list is empty")
    allowed_mode_set = {str(m) for m in allowed_modes}
    if "real" in allowed_mode_set and allowed_real_keys is None:
        raise RuntimeError("science guard failed: real mode requested without explicit allowed_real_keys")
    frames: List[pd.DataFrame] = []
    for p in _part_files(rows_dir):
        tbl = pq.read_table(p, filters=[("key", "in", key_list)])
        if tbl.num_rows:
            frames.append(tbl.to_pandas())
    if not frames:
        raise MissingInput("no science rows matched the requested keys")
    sci = pd.concat(frames, ignore_index=True)
    sci = sci[sci["status"] == "ok"].drop_duplicates("key").reset_index(drop=True)
    got_secs = set(sci["sector"].astype(int))
    if 69 in got_secs or not got_secs.issubset({int(s) for s in sectors}):
        raise RuntimeError(f"science guard failed: unexpected sectors {sorted(got_secs)}")
    if sci["holdout"].astype(bool).any():
        raise RuntimeError(f"science guard failed: {int(sci['holdout'].astype(bool).sum())} holdout rows present")
    got_modes = set(sci["mode"].astype(str))
    if not got_modes.issubset(allowed_mode_set):
        raise RuntimeError(f"science guard failed: unexpected modes {sorted(got_modes)}")
    if "real" in got_modes:
        real_keys = set(sci.loc[sci["mode"].astype(str) == "real", "key"].astype(str))
        extra_real = real_keys - (allowed_real_keys or set())
        if extra_real:
            raise RuntimeError(f"science guard failed: {len(extra_real)} unauthorized real row(s)")
    if len(sci) != len(key_list):
        raise MissingInput(f"missing {len(key_list) - len(sci)} ok science row(s) for requested keys")
    return sci


# --------------------------------------------------------------------------- (a) Population, EB drop, three-way match
def rebuild_exclusion(
    manifest: Union[str, Path, Sequence[Mapping[str, Any]]], sectors: Sequence[int] = CALIBRATION_SECTORS
) -> Dict[int, Set[int]]:
    """Per sector, TICs with an invert task minus TICs with a scramble task."""
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
) -> Tuple[Dict[int, Set[int]], Dict[int, Set[int]], Dict[str, Dict[str, int]]]:
    """Registered exclusion per sector: TOI (any disposition among searched targets) U Kostov 2025."""
    jobs = _jobs_list(manifest)
    sec_set = {int(s) for s in sectors}
    inv: Dict[int, Set[int]] = {s: set() for s in sec_set}
    for j in jobs:
        s = int(j["sector"])
        if s in sec_set and any(t.get("mode") == "invert" for t in j.get("tasks", ())):
            inv[s].add(int(j["tic"]))

    toi_df, _ = known_planets.load_toi_table(_frozen_toi_path())
    toi_tics = set(toi_df["tic"].astype(int))

    kostov_path = config.CATALOG_DIR / KOSTOV_RELPATH
    if not kostov_path.is_file():
        raise MissingInput(f"missing Kostov 2025 exclusion table: {KOSTOV_RELPATH}")
    kostov_got = _sha256_bytes(kostov_path.read_bytes())
    if kostov_got != PINNED_KOSTOV_SHA256:
        raise MissingInput(
            f"Kostov 2025 exclusion table sha256 mismatch: expected {PINNED_KOSTOV_SHA256}, got {kostov_got}"
        )
    kdf = pd.read_csv(kostov_path)
    kostov_by_sec: Dict[int, Set[int]] = {s: set() for s in sec_set}
    for tic_val, sec_cell in zip(kdf["tic"], kdf["sectors_in_59-64_69_70"]):
        if pd.isna(tic_val):
            continue
        tic = int(tic_val)
        for s in _parse_sectors_cell(sec_cell) & sec_set:
            kostov_by_sec[s].add(tic)

    sets: Dict[int, Set[int]] = {}
    kostov_in_inv: Dict[int, Set[int]] = {}
    counts: Dict[str, Dict[str, int]] = {}
    for s in sorted(sec_set):
        in_toi = inv[s] & toi_tics
        in_kostov = inv[s] & kostov_by_sec[s]
        union = in_toi | in_kostov
        sets[s] = union
        kostov_in_inv[s] = in_kostov
        counts[str(s)] = {"toi": len(in_toi), "kostov": len(in_kostov), "registered": len(union)}
    return sets, kostov_in_inv, counts


def check_exclusion_with_eb_drop(
    manifest_excl: Mapping[int, Set[int]],
    registered: Mapping[int, Set[int]],
    kostov_sets: Mapping[int, Set[int]],
    reg_counts: Mapping[str, Mapping[str, int]],
    expected_eb_drop: Mapping[int, int] = EXPECTED_EB_DROP,
) -> Dict[str, Any]:
    """Verify EB deviation (Deviation 2026-11-20): eb_drop in Kostov, exact counts, manifest_excl <= registered."""
    per_sector: Dict[str, Dict[str, int]] = {}
    eb_drop: Dict[int, Set[int]] = {}
    effective: Dict[int, Set[int]] = {}
    mismatch = False

    for s in sorted(set(manifest_excl) | set(registered) | set(expected_eb_drop)):
        me = set(manifest_excl.get(s, ()))
        rg = set(registered.get(s, ()))
        ks = set(kostov_sets.get(s, ()))
        ed = rg - me
        extra_me = me - rg
        not_in_kostov = ed - ks
        eff = me | ed
        exp_n = expected_eb_drop.get(s)
        sec_bad = bool(extra_me or not_in_kostov or (exp_n is None) or (len(ed) != exp_n) or (eff != rg))
        if sec_bad:
            mismatch = True
        eb_drop[s] = ed
        effective[s] = eff
        rc = reg_counts.get(str(s), {})
        per_sector[str(int(s))] = {
            "manifest_excl": len(me),
            "toi": int(rc.get("toi", 0)),
            "kostov": int(rc.get("kostov", len(ks))),
            "eb_drop": len(ed),
            "expected_eb_drop": -1 if exp_n is None else int(exp_n),
            "registered": len(rg),
            "effective": len(eff),
            "extra_in_manifest_excl": len(extra_me),
            "eb_drop_not_in_kostov": len(not_in_kostov),
        }

    if mismatch:
        raise ExclusionMismatch({"status": "mismatch", "per_sector": per_sector})

    canonical = "".join(sorted(f"{int(s)},{int(t)}\n" for s, tics in effective.items() for t in tics)).encode("utf-8")
    digest = _sha256_bytes(canonical)
    return {
        "eb_drop": eb_drop,
        "effective": effective,
        "registered": {int(s): set(registered[s]) for s in sorted(registered)},
        "sha256": digest,
        "report": {"status": "ok", "effective_sha256": digest, "per_sector": per_sector},
    }


def apply_eb_drop(
    identity: pd.DataFrame, eb_drop: Mapping[int, Set[int]], out_dir: Union[str, Path]
) -> Tuple[pd.DataFrame, Dict[str, Any]]:
    """Drop eb_drop star-sectors from all modes in identity and write out-dir/eb_dropped_star_sectors.csv."""
    pairs = sorted((int(s), int(t)) for s, tics in eb_drop.items() for t in tics)
    out_p = Path(out_dir)
    out_p.mkdir(parents=True, exist_ok=True)
    csv_path = out_p / "eb_dropped_star_sectors.csv"
    eb_df = pd.DataFrame(pairs, columns=["sector", "tic"])
    eb_df.to_csv(csv_path, index=False)
    digest = _sha256_bytes(csv_path.read_bytes())

    pair_set = set(pairs)
    if len(identity) and pair_set:
        keep_mask = [
            (int(s), int(t)) not in pair_set
            for s, t in zip(identity["sector"].astype(int), identity["tic"].astype(int))
        ]
        filtered = identity.loc[keep_mask].reset_index(drop=True)
    else:
        filtered = identity.reset_index(drop=True)
    n_rows_dropped = int(len(identity) - len(filtered))
    return filtered, {
        "n_star_sectors": len(pairs),
        "n_identity_rows_dropped": n_rows_dropped,
        "csv_sha256": digest,
    }


def check_three_way_match(
    identity: pd.DataFrame, registered: Mapping[int, Set[int]], sectors: Sequence[int] = CALIBRATION_SECTORS
) -> Tuple[Dict[str, Any], Dict[str, Any]]:
    """Verify three-way structural match on post-EB-drop identity rows; count status != ok per mode."""
    per_sector: Dict[str, Any] = {}
    status_errors_per_sector: Dict[str, Dict[str, int]] = {}
    status_errors_total: Dict[str, int] = {"real": 0, "invert": 0, "scramble": 0}
    total_violations = 0

    for s in sorted(int(x) for x in sectors):
        sub = identity[identity["sector"].astype(int) == s]
        reg_s = set(registered.get(s, ()))
        modes_series = sub["mode"].astype(str)
        unexpected_mode_rows = int((~modes_series.isin(["real", "invert", "scramble"])).sum())

        sec_err: Dict[str, int] = {}
        for m in ("real", "invert", "scramble"):
            cnt = int(((modes_series == m) & (sub["status"].astype(str) != "ok")).sum())
            sec_err[m] = cnt
            status_errors_total[m] += cnt
        status_errors_per_sector[str(s)] = sec_err

        counts_by_tic: Dict[int, List[int]] = {}
        for tic_val, mode_val in zip(sub["tic"].astype(int), modes_series):
            arr = counts_by_tic.get(tic_val)
            if arr is None:
                arr = [0, 0, 0]
                counts_by_tic[tic_val] = arr
            if mode_val == "real":
                arr[0] += 1
            elif mode_val == "invert":
                arr[1] += 1
            elif mode_val == "scramble":
                arr[2] += 1

        viol = {
            "unexpected_mode_rows": unexpected_mode_rows,
            "unexcluded_missing_real": 0,
            "unexcluded_missing_invert": 0,
            "unexcluded_missing_scramble": 0,
            "unexcluded_dup_real": 0,
            "unexcluded_dup_invert": 0,
            "unexcluded_dup_scramble": 0,
            "excluded_missing_real": 0,
            "excluded_missing_invert": 0,
            "excluded_dup_real": 0,
            "excluded_dup_invert": 0,
            "excluded_has_scramble": 0,
        }
        n_unexcl = 0
        n_excl = 0
        for tic_val, (nr, ni, ns) in counts_by_tic.items():
            if tic_val in reg_s:
                n_excl += 1
                if nr == 0:
                    viol["excluded_missing_real"] += 1
                elif nr > 1:
                    viol["excluded_dup_real"] += 1
                if ni == 0:
                    viol["excluded_missing_invert"] += 1
                elif ni > 1:
                    viol["excluded_dup_invert"] += 1
                if ns > 0:
                    viol["excluded_has_scramble"] += 1
            else:
                n_unexcl += 1
                if nr == 0:
                    viol["unexcluded_missing_real"] += 1
                elif nr > 1:
                    viol["unexcluded_dup_real"] += 1
                if ni == 0:
                    viol["unexcluded_missing_invert"] += 1
                elif ni > 1:
                    viol["unexcluded_dup_invert"] += 1
                if ns == 0:
                    viol["unexcluded_missing_scramble"] += 1
                elif ns > 1:
                    viol["unexcluded_dup_scramble"] += 1

        sec_viol_sum = sum(viol.values())
        total_violations += sec_viol_sum
        per_sector[str(s)] = {
            "n_stars_unexcluded": n_unexcl,
            "n_stars_excluded": n_excl,
            "n_rows_real": int((modes_series == "real").sum()),
            "n_rows_invert": int((modes_series == "invert").sum()),
            "n_rows_scramble": int((modes_series == "scramble").sum()),
            "violations": viol,
            "total_violations": sec_viol_sum,
        }

    tw_report = {
        "status": "ok" if total_violations == 0 else "mismatch",
        "total_violations": total_violations,
        "per_sector": per_sector,
    }
    err_report = {"total": status_errors_total, "per_sector": status_errors_per_sector}
    if total_violations > 0:
        raise ThreeWayMismatch({"three_way_match": tw_report, "status_errors": err_report})
    return tw_report, err_report


def quality_drop(
    quality_dir: Union[str, Path],
    expected_sha256s: Mapping[int, str],
    sectors: Sequence[int] = CALIBRATION_SECTORS,
) -> pd.DataFrame:  # review: S5
    """Read star_quality_s00XX.parquet from quality_dir, verify manifest & sidecar sha256, and apply drop rule."""
    if not expected_sha256s:
        raise SystemExit("quality_drop requires expected_sha256s from verified --input-hashes file")
    sec_list = [int(s) for s in sectors]
    qdir = Path(quality_dir)
    files = [qdir / f"star_quality_s{s:04d}.parquet" for s in sec_list]
    missing = [s for s, f in zip(sec_list, files) if not f.is_file()]
    if missing:
        raise MissingInput(f"missing star-quality tables for {len(missing)} sector(s): {missing}")

    sha_map: Dict[int, str] = {}
    for s, f in zip(sec_list, files):
        expected = expected_sha256s.get(s)
        if expected is None:
            raise SystemExit(f"missing expected sha256 from --input-hashes for sector {s}")
        pinned = PINNED_QUALITY_SHA256.get(s)
        if pinned is not None and expected.lower() != pinned.lower():
            raise SystemExit(
                f"star-quality sha256 consistency mismatch for sector {s}: "
                f"--input-hashes {expected} != pinned {pinned}"
            )
        got = _sha256_bytes(f.read_bytes())
        if got.lower() != expected.lower():
            raise MissingInput(
                f"star-quality sha256 mismatch for sector {s}: expected {expected}, got {got}"
            )
        sidecar = f.with_name(f.name + ".sha256")
        if not sidecar.is_file():
            raise MissingInput(f"missing star-quality sha256 sidecar for sector {s}")
        sidecar_tokens = sidecar.read_text(encoding="utf-8").strip().split()
        sidecar_hash = sidecar_tokens[0] if sidecar_tokens else ""
        if sidecar_hash.lower() != got.lower():
            raise MissingInput(
                f"star-quality sidecar sha256 mismatch for sector {s}: sidecar {sidecar_hash}, got {got}"
            )
        sha_map[s] = got

    thr = quality.rule()["threshold"]
    q = pd.concat([pd.read_parquet(f) for f in files], ignore_index=True)
    q["drop"] = (q["status"] != "ok") | ~(q["beta_max"] <= thr)
    present = set(q["sector"].astype(int)) if len(q) else set()
    empty_secs = [s for s in sec_list if s not in present]
    if empty_secs:
        raise MissingInput(f"empty star-quality tables for {len(empty_secs)} sector(s): {empty_secs}")
    q.attrs["sha256s"] = sha_map
    return q


def alpha_n(identity: pd.DataFrame, q: pd.DataFrame) -> Dict[int, Dict[str, Any]]:
    """N_s = quality-passing invert rows (holdout included, no status filter, post-EB-drop), p_s = 1/N_s."""
    m = interim.with_quality(identity, q)
    out: Dict[int, Dict[str, Any]] = {}
    for s in sorted(int(x) for x in identity["sector"].unique()):
        sub = m[(m["sector"] == s) & (m["mode"] == "invert")].drop_duplicates("tic")
        n_s = int((~sub["drop"]).sum())
        if n_s <= 0:
            raise MissingInput(f"sector {s} has zero quality-passing invert rows")
        out[s] = {"N_s": n_s, "p_s": 1.0 / n_s}
    return out


def null_sample(
    identity: pd.DataFrame,
    q: pd.DataFrame,
    exclusion: Mapping[int, Set[int]],
    alpha: Mapping[int, Mapping[str, Any]],
) -> Tuple[Dict[str, List[str]], Dict[str, Dict[int, int]], Dict[str, Any]]:
    """Paired null sample per mode (non-holdout, quality-pass, status ok, ok real+invert+scramble)."""
    m = interim.with_quality(identity, q)
    sel = m[(~m["drop"]) & (~m["holdout"].astype(bool)) & (m["status"] == "ok")]
    keys_by_mode: Dict[str, List[str]] = {mode: [] for mode in NULL_MODES}
    n_by_mode_sec: Dict[str, Dict[int, int]] = {mode: {} for mode in NULL_MODES}
    per_sector: Dict[str, Dict[str, Any]] = {}

    sectors = sorted(int(x) for x in identity["sector"].unique())
    for s in sectors:
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
        keys_by_mode["invert"].extend(inv_keys)
        keys_by_mode["scramble"].extend(scr_keys)
        n_by_mode_sec["invert"][s] = len(inv_keys)
        n_by_mode_sec["scramble"][s] = len(scr_keys)
        per_sector[str(s)] = {
            "n_stars": len(common),
            "n_invert": len(inv_keys),
            "n_scramble": len(scr_keys),
            "guard": guard,
        }

    by_mode_summary: Dict[str, Dict[str, Any]] = {}
    for mode in NULL_MODES:
        keys_by_mode[mode] = sorted(set(keys_by_mode[mode]))
        n_null_m = sum(n_by_mode_sec[mode].values())
        if n_null_m <= 0:
            raise MissingInput(f"mode {mode} has zero null-sample rows")
        b_m = sum(n_by_mode_sec[mode][s] * float(alpha[s]["p_s"]) for s in sectors)
        p_m = b_m / n_null_m
        need_m = int(math.floor(b_m + _BUDGET_EPS)) + 1
        by_mode_summary[mode] = {
            "n_by_sector": {str(s): int(n_by_mode_sec[mode][s]) for s in sectors},
            "n_null": int(n_null_m),
            "B": float(b_m),
            "p": float(p_m),
            "need": int(need_m),
        }

    report = {
        "by_mode": by_mode_summary,
        "per_sector": per_sector,
        "total_null_rows_union": len(set(keys_by_mode["invert"]) | set(keys_by_mode["scramble"])),
    }
    return keys_by_mode, n_by_mode_sec, report


# --------------------------------------------------------------------------- (b) Reproduction gate (--gate)
def _flatten_dict(obj: Mapping[str, Any], prefix: str = "") -> Dict[str, Any]:
    out: Dict[str, Any] = {}
    for k, v in obj.items():
        full_k = f"{prefix}.{k}" if prefix else str(k)
        if isinstance(v, Mapping):
            out.update(_flatten_dict(v, full_k))
        else:
            out[full_k] = v
    return out


def check_checkout_and_cfg(header: Mapping[str, Any]) -> Dict[str, Any]:
    """Verify checkout commit, manifest cfg equality with settings.cfg(), and manifest cfg vs v1.2 cfg."""
    got_commit = str(settings.commit())
    want_commit = str(header["pinned_commit"])
    if got_commit != PINNED_CODE_COMMIT:
        raise SystemExit(
            f"checkout commit mismatch: settings.commit()={got_commit!r} != PINNED_CODE_COMMIT {PINNED_CODE_COMMIT!r}"
        )
    if got_commit != want_commit:
        raise SystemExit(f"checkout commit mismatch: settings.commit()={got_commit!r} != manifest {want_commit!r}")

    got_cfg_json = json.dumps(settings.cfg(), sort_keys=True)
    want_cfg_json = json.dumps(header["cfg"], sort_keys=True)
    if got_cfg_json != want_cfg_json:
        raise SystemExit("checkout settings.cfg() differs from manifest header['cfg']")

    v12_cfg = v12_settings.require_v12_settings()
    flat_hdr = _flatten_dict(header["cfg"])
    flat_v12 = _flatten_dict(v12_cfg)
    diff_keys = sorted(
        k
        for k in (set(flat_hdr) | set(flat_v12)) - V12_ALLOWED_CFG_DIFF_KEYS
        if (k not in flat_hdr) or (k not in flat_v12) or (flat_hdr[k] != flat_v12[k])
    )
    if diff_keys:
        raise SystemExit(
            f"manifest header['cfg'] and frozen v1.2 config disagree on non-exempt keys: {diff_keys}"
        )
    return {
        "status": "ok",
        "pinned_commit": got_commit,
        "allowed_v12_diff_keys": sorted(V12_ALLOWED_CFG_DIFF_KEYS),
    }


def _sort_null_frame(df: pd.DataFrame) -> pd.DataFrame:
    return df.sort_values(by=["S", "key"], ascending=[False, True], kind="mergesort", na_position="last").reset_index(
        drop=True
    )


def select_known_planet_real_keys(
    identity: pd.DataFrame, sectors: Sequence[int] = CALIBRATION_SECTORS
) -> List[str]:
    """Select all post-EB-drop, non-holdout, ok real-row keys for known-planet targets per sector."""
    ps_df, _ = known_planets.load_pscomppars()
    toi_df, _ = known_planets.load_toi_table(_frozen_toi_path())
    pgh_table = p_grid_halfwidth.load()

    kp_keys: List[str] = []
    for s in sorted(int(x) for x in sectors):
        real_s = identity[
            (identity["sector"].astype(int) == s)
            & (identity["mode"].astype(str) == "real")
            & (~identity["holdout"].astype(bool))
            & (identity["status"].astype(str) == "ok")
        ]
        searched_tics = set(real_s["tic"].astype(int))
        if not searched_tics:
            continue
        targets = known_planets.select_known_planet_targets(
            sector=s,
            searched_tics=searched_tics,
            ps_df=ps_df,
            toi_df=toi_df,
            pgh_table=pgh_table,
        )
        kp_tics = set(targets["tic"].astype(int)) if len(targets) else set()
        if kp_tics:
            kp_keys.extend(sorted(real_s[real_s["tic"].astype(int).isin(kp_tics)]["key"].astype(str).unique()))
    return sorted(set(kp_keys))


def build_gate_keys(
    identity: pd.DataFrame,
    rows_dir: Union[str, Path],
    keys_by_mode: Mapping[str, Sequence[str]],
    gate_top: int,
    sectors: Sequence[int] = CALIBRATION_SECTORS,
) -> Dict[str, Any]:
    """Deterministically build the gate key set (known-planet real + top-K per null mode + random 200 null)."""
    kp_real_keys = select_known_planet_real_keys(identity, sectors=sectors)
    union_null_keys = sorted(set(keys_by_mode["invert"]) | set(keys_by_mode["scramble"]))
    s_df = load_null_s_for_keys(rows_dir, union_null_keys)

    top_by_mode: Dict[str, List[str]] = {}
    for mode in NULL_MODES:
        mode_set = set(keys_by_mode[mode])
        sub = _sort_null_frame(s_df[s_df["key"].astype(str).isin(mode_set)])
        top_by_mode[mode] = sub["key"].astype(str).head(int(gate_top)).tolist()

    rng = np.random.default_rng(GATE_RANDOM_SEED)
    n_rand = min(GATE_RANDOM_N, len(union_null_keys))
    rand_keys = sorted(str(k) for k in rng.choice(union_null_keys, size=n_rand, replace=False).tolist())

    gate_keys = sorted(set(kp_real_keys) | set(top_by_mode["invert"]) | set(top_by_mode["scramble"]) | set(rand_keys))
    return {
        "gate_keys": gate_keys,
        "sha256": _gate_keys_sha256(gate_keys),
        "known_planet_real": kp_real_keys,
        "top_invert": top_by_mode["invert"],
        "top_scramble": top_by_mode["scramble"],
        "random_null": rand_keys,
        "n_per_group": {
            "known_planet_real": len(kp_real_keys),
            "top_invert": len(top_by_mode["invert"]),
            "top_scramble": len(top_by_mode["scramble"]),
            "random_null": len(rand_keys),
            "union_unique": len(gate_keys),
        },
    }


def _precheck_gated_light_curves(
    gate_keys: Sequence[str],
    identity: pd.DataFrame,
    jobs: Sequence[Mapping[str, Any]],
) -> int:  # review: S1
    """Verify every gated key has a local light curve via exosearch.batch.lc_path(..., download=False)."""
    job_fname: Dict[Tuple[int, int], str] = {
        (int(j["sector"]), int(j["tic"])): str(j["fname"]) for j in jobs if j.get("fname")
    }
    gate_set = set(str(k) for k in gate_keys)
    sub = identity[identity["key"].astype(str).isin(gate_set)].drop_duplicates("key")
    id_map: Dict[str, Tuple[int, int]] = {
        str(k): (int(s), int(t))
        for k, s, t in zip(sub["key"], sub["sector"], sub["tic"])
    }
    checked_pairs: Dict[Tuple[int, int], bool] = {}
    n_missing = 0
    for key in gate_keys:
        pair = id_map.get(str(key))
        if pair is None:
            n_missing += 1
            continue
        if pair not in checked_pairs:
            fname = job_fname.get(pair)
            checked_pairs[pair] = bool(fname and batch.lc_path(pair[0], fname, download=False) is not None)
        if not checked_pairs[pair]:
            n_missing += 1
    if n_missing > 0:
        raise SystemExit(f"missing local light curves for {n_missing} gated key(s) (download=False precheck)")
    return len(gate_keys)


def _get_worker_finish() -> Any:
    global _WORKER_FINISH
    if _WORKER_FINISH is None:
        _WORKER_FINISH = _load_pinned_finish()
    return _WORKER_FINISH


def _recompute_gate_job(
    payload: Tuple[str, Dict[str, Any], Dict[str, Any]]
) -> Tuple[str, Optional[Dict[str, Any]], Optional[str], Optional[str], Optional[str]]:  # review: B1
    """Top-level worker for ProcessPoolExecutor: run pipeline.run_job and (for null rows) LC rebuild check."""
    kernels.install()  # review: B2
    key, job_hdr, row_dict = payload
    try:
        rows = pipeline.run_job(job_hdr)
    except Exception as exc:
        return key, None, f"recompute:{type(exc).__name__}", None, None
    if len(rows) != 1:
        return key, None, f"unexpected_row_count:{len(rows)}", None, None
    r = rows[0]
    keep = ("status",) + GATE_S_COLS + GATE_EXACT_COLS
    slim = {c: r.get(c) for c in keep if c in r}

    mode = str(row_dict.get("mode", ""))
    rebuild_status: Optional[str] = None
    rebuild_err: Optional[str] = None
    if mode in NULL_MODES:
        try:
            finish_mod = _get_worker_finish()
            v12_cfg = v12_settings.require_v12_settings()
            job_v12 = dict(job_hdr, cfg=v12_cfg)
            flat_rb = finish_mod._rebuild_flat_dumps(row_dict, job_v12)[:3]

            job = job_hdr
            tk = job["tasks"][0]
            t, f, e, good, meta, dumps = pipeline.load(job, download=False)
            pw_on = tk.get("prewhiten", True)
            # Verbatim copy of pipeline.run_job (exosearch/v11/pipeline.py) from
            # t[good], f[good], e[good] through prepare(...), omitting row.update(scr_info)
            # and qualifying prepare as pipeline.prepare.
            tg, fg, eg = t[good], f[good], e[good]
            if tk["mode"] == "invert":
                tg, fg, eg = nulltest.transform(tk["mode"], tg, fg, eg, meta, tk["seed"])
            elif tk["mode"] == "scramble":
                tg, fg, eg, scr_info = scramble.segment_rotate(tg, fg, eg, tk["seed"], job["cfg"])
            elif tk["mode"] != "real":
                raise ValueError(tk["mode"])
            gap_trim = 0.0 if tk["mode"] == "scramble" else 0.3
            d = pipeline.prepare(tg, fg, eg, job["cfg"], pw_on, tk.get("comb_guard"), gap_trim=gap_trim)
            flat_search = d["flat"]

            bit_identical = (
                len(flat_rb) == 3
                and len(flat_search) == 3
                and all(
                    np.asarray(a).shape == np.asarray(b).shape
                    and np.array_equal(a, b, equal_nan=True)
                    for a, b in zip(flat_rb, flat_search)
                )
            )
            if bit_identical:
                rebuild_status = "flat_identical"
            else:
                rebuild_status = "rebuild_fail"
                rebuild_err = "rebuild:flat_not_identical"
        except Exception as exc:
            rebuild_status = "rebuild_fail"
            rebuild_err = f"rebuild:{type(exc).__name__}"

    return key, slim, None, rebuild_status, rebuild_err


def _is_nan_or_none(val: Any) -> bool:
    if val is None:
        return True
    try:
        return bool(pd.isna(val))
    except (TypeError, ValueError):
        return False


def _compare_gate_row(
    stored: Mapping[str, Any], recomp: Optional[Mapping[str, Any]], err: Optional[str]
) -> Tuple[bool, Optional[str], Dict[str, Optional[float]], Dict[str, bool]]:
    """Compare one recomputed gate row against the stored row."""
    diffs: Dict[str, Optional[float]] = {}
    exact: Dict[str, bool] = {}
    if err is not None or recomp is None:
        return False, err or "recompute_none", diffs, exact
    status = str(recomp.get("status"))
    if status != "ok":
        return False, f"status:{status}", diffs, exact

    failed_cols: List[str] = []
    for col in GATE_S_COLS:
        if col != "S" and (col not in stored) and (col not in recomp):
            continue
        sv = stored.get(col)
        rv = recomp.get(col)
        s_na = _is_nan_or_none(sv)
        r_na = _is_nan_or_none(rv)
        if s_na and r_na:
            diffs[col] = None
            continue
        if s_na != r_na:
            failed_cols.append(f"{col}:nan_parity")
            continue
        sf = float(sv)
        rf = float(rv)
        if not (math.isfinite(sf) and math.isfinite(rf)):
            failed_cols.append(f"{col}:non_finite")
            continue
        d = abs(sf - rf)
        diffs[col] = d
        if d > GATE_S_TOL:
            failed_cols.append(f"{col}:diff>{GATE_S_TOL}")

    for col in GATE_EXACT_COLS:
        sv = stored.get(col)
        rv = recomp.get(col)
        s_na = _is_nan_or_none(sv)
        r_na = _is_nan_or_none(rv)
        if s_na and r_na:
            exact[col] = True
        elif s_na != r_na:
            exact[col] = False
        else:
            exact[col] = float(sv) == float(rv)

    if failed_cols:
        return False, ";".join(failed_cols), diffs, exact
    return True, None, diffs, exact


def run_gate(
    args: argparse.Namespace,
    header: Mapping[str, Any],
    jobs: Sequence[Mapping[str, Any]],
    identity: pd.DataFrame,
    keys_by_mode: Mapping[str, Sequence[str]],
    base_report: Dict[str, Any],
    checkout_report: Dict[str, Any],
) -> int:
    """Execute --gate: recompute all gate rows with pinned pipeline.run_job and enforce <= 1e-6 match."""
    out_p = Path(args.out_dir)
    out_p.mkdir(parents=True, exist_ok=True)
    (out_p / "gate_pass.json").unlink(missing_ok=True)  # review: B4

    gate_info = build_gate_keys(identity, args.rows_dir, keys_by_mode, args.gate_top, CALIBRATION_SECTORS)
    gate_keys = gate_info["gate_keys"]
    _precheck_gated_light_curves(gate_keys, identity, jobs)  # review: S1

    kp_real_set = set(gate_info["known_planet_real"])
    top_inv_set = set(gate_info["top_invert"])
    top_scr_set = set(gate_info["top_scramble"])
    rand_set = set(gate_info["random_null"])
    group_sets: Dict[str, Set[str]] = {
        "known_planet_real": kp_real_set,
        "top_invert": top_inv_set,
        "top_scramble": top_scr_set,
        "random_null": rand_set,
    }

    sci = load_science_for_keys(
        args.rows_dir,
        gate_keys,
        sectors=CALIBRATION_SECTORS,
        allowed_modes=("real", "invert", "scramble"),
        allowed_real_keys=kp_real_set,
    )
    stored_by_key: Dict[str, Dict[str, Any]] = {str(r["key"]): r for r in sci.to_dict(orient="records")}
    job_map: Dict[Tuple[int, int], Mapping[str, Any]] = {(int(j["sector"]), int(j["tic"])): j for j in jobs}

    payloads: List[Tuple[str, Dict[str, Any], Dict[str, Any]]] = []
    pre_errors: Dict[str, str] = {}
    for key in gate_keys:
        row = stored_by_key[key]
        sec = int(row["sector"])
        tic = int(row["tic"])
        task_tag = str(row["task"])
        base_job = job_map.get((sec, tic))
        if base_job is None:
            pre_errors[key] = "missing_manifest_job"
            continue
        matched_tasks = [tk for tk in base_job.get("tasks", ()) if str(tk.get("tag")) == task_tag]
        if len(matched_tasks) != 1:
            pre_errors[key] = f"matched_tasks:{len(matched_tasks)}"
            continue
        job_copy = dict(base_job)
        job_copy["cfg"] = header["cfg"]
        job_copy["cfg_inj"] = header.get("cfg_inj")
        job_copy["primary"] = header["cfg"]["beta_option"]
        job_copy["prov"] = {}
        job_copy["tasks"] = matched_tasks
        payloads.append((key, job_copy, dict(row)))

    results: Dict[str, Tuple[Optional[Dict[str, Any]], Optional[str], Optional[str], Optional[str]]] = {
        k: (None, e, None, None) for k, e in pre_errors.items()
    }
    if payloads:
        with concurrent.futures.ProcessPoolExecutor(
            max_workers=int(args.workers),
            mp_context=multiprocessing.get_context("spawn"),
        ) as pool:
            for key, recomp, err, rb_status, rb_err in pool.map(_recompute_gate_job, payloads):
                results[key] = (recomp, err, rb_status, rb_err)

    max_abs_diff: Dict[str, Optional[float]] = {c: None for c in GATE_S_COLS}
    exact_equal: Dict[str, int] = {c: 0 for c in GATE_EXACT_COLS}
    exact_total: Dict[str, int] = {c: 0 for c in GATE_EXACT_COLS}
    rebuild_by_group: Dict[str, Dict[str, int]] = {  # review: B1
        g: {"n_flat_identical": 0, "n_rebuild_fail": 0}
        for g in group_sets
    }
    n_pass = 0
    n_fail = 0
    n_recompute_err = 0
    n_status_not_ok = 0
    n_s_mismatch = 0
    n_rebuild_fail = 0
    fail_rows: List[Dict[str, Any]] = []

    for key in gate_keys:
        recomp, err, rb_status, rb_err = results[key]
        ok, reason, diffs, exact = _compare_gate_row(stored_by_key[key], recomp, err)
        for col, d in diffs.items():
            if d is not None:
                prev = max_abs_diff[col]
                max_abs_diff[col] = d if prev is None else max(prev, d)
        for col, eq in exact.items():
            exact_total[col] += 1
            if eq:
                exact_equal[col] += 1

        if rb_status is not None:
            rb_field = {
                "flat_identical": "n_flat_identical",
                "rebuild_fail": "n_rebuild_fail",
            }[rb_status]
            for gname, gset in group_sets.items():
                if key in gset:
                    rebuild_by_group[gname][rb_field] += 1

        rb_failed = (rb_status == "rebuild_fail") or (rb_err is not None)
        if ok and not rb_failed:
            n_pass += 1
        else:
            n_fail += 1
            reasons: List[str] = []
            if not ok:
                rstr = reason or "unknown"
                reasons.append(rstr)
                if err is not None or recomp is None:
                    n_recompute_err += 1
                elif rstr.startswith("status:"):
                    n_status_not_ok += 1
                else:
                    n_s_mismatch += 1
            if rb_failed:
                n_rebuild_fail += 1
                reasons.append(rb_err or "rebuild_fail")
            fail_rows.append({"key": key, "reason": ";".join(reasons)})

    driver_sha = _driver_sha256()
    manifest_sha = _sha256_bytes(Path(args.manifest).read_bytes())
    gate_summary: Dict[str, Any] = {
        "passed": n_fail == 0,
        "gate_top": int(args.gate_top),
        "n_total": len(gate_keys),
        "n_pass": n_pass,
        "n_fail": n_fail,
        "n_recompute_error": n_recompute_err,
        "n_status_not_ok": n_status_not_ok,
        "n_s_mismatch": n_s_mismatch,
        "n_rebuild_fail": n_rebuild_fail,
        "n_per_group": gate_info["n_per_group"],
        "rebuild_by_group": rebuild_by_group,
        "max_abs_diff": max_abs_diff,
        "exact_equality": {
            c: {"equal": exact_equal[c], "total": exact_total[c]} for c in GATE_EXACT_COLS
        },
        "s_tol": GATE_S_TOL,
        "pinned_commit": checkout_report["pinned_commit"],
        "driver_sha256": driver_sha,
        "manifest_sha256": manifest_sha,
    }

    if n_fail > 0:
        pd.DataFrame(fail_rows, columns=["key", "reason"]).to_csv(out_p / "gate_fail_keys.csv", index=False)
        _write_json(out_p / "gate_report.json", gate_summary)
        full_report = {**base_report, "checkout_check": checkout_report, "gate": gate_summary}
        payload = _write_json(out_p / "report_gate.json", full_report)
        print(payload, end="")
        return 3

    gate_pass = {**gate_summary, "sha256": gate_info["sha256"]}
    _write_json(out_p / "gate_pass.json", gate_pass)
    gate_key_rows = [
        {
            "key": k,
            "in_known_planet_real": k in kp_real_set,
            "in_top_invert": k in top_inv_set,
            "in_top_scramble": k in top_scr_set,
            "in_random_null": k in rand_set,
        }
        for k in gate_keys
    ]
    pd.DataFrame(gate_key_rows).to_csv(out_p / "gate_keys.csv", index=False)
    full_report = {**base_report, "checkout_check": checkout_report, "gate": gate_pass}
    payload = _write_json(out_p / "report_gate.json", full_report)
    print(payload, end="")
    return 0


# --------------------------------------------------------------------------- (c) Vetting and S*_cal (--compute)
def _load_pinned_finish() -> Any:
    """Import <repo>/runs/v1.1_s70_x86/finish.py via importlib."""
    repo_root = Path(exosearch.__file__).resolve().parents[1]
    finish_path = repo_root / "runs" / "v1.1_s70_x86" / "finish.py"
    if not finish_path.is_file():
        raise MissingInput(f"pinned finish module not found at {finish_path}")
    spec = importlib.util.spec_from_file_location("exosearch_v11_s70_finish", finish_path)
    if spec is None or spec.loader is None:
        raise MissingInput(f"unable to create module spec for {finish_path}")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _flag_common_periods(sci: pd.DataFrame) -> pd.DataFrame:
    """Flag common_period per (sector, task) group of the null sample."""
    group_col = "task" if "task" in sci.columns else "mode"
    parts: List[pd.DataFrame] = []
    for _, g in sci.groupby(["sector", group_col], sort=True):
        gf = vet.flag_common_periods(g.assign(status="ok"))
        gf["fail_common_period"] = gf["flags"].fillna("").str.contains("common_period")
        parts.append(gf)
    return pd.concat(parts, ignore_index=True)


def _common_period_population_identity(
    identity: pd.DataFrame,
    common_period_mode: str,
    keys_by_mode: Optional[Mapping[str, Sequence[str]]] = None,
    sectors: Sequence[int] = CALIBRATION_SECTORS,
) -> pd.DataFrame:
    """Identity rows belonging to the selected common_period flagging population (post-EB-drop)."""
    sec_set = {int(s) for s in sectors}
    if common_period_mode == "v11_task_with_holdout":
        mask = (
            identity["sector"].astype(int).isin(sec_set)
            & identity["mode"].astype(str).isin(NULL_MODES)
            & (identity["status"].astype(str) == "ok")
        )
        return identity.loc[mask].reset_index(drop=True)
    if common_period_mode == "v11_task_no_holdout":
        mask = (
            identity["sector"].astype(int).isin(sec_set)
            & identity["mode"].astype(str).isin(NULL_MODES)
            & (~identity["holdout"].astype(bool))
            & (identity["status"].astype(str) == "ok")
        )
        return identity.loc[mask].reset_index(drop=True)
    if common_period_mode == "null_sample":
        if keys_by_mode is None:
            raise SystemExit("null_sample common_period population requires keys_by_mode")
        null_key_set = set(str(k) for k in keys_by_mode.get("invert", ())) | set(
            str(k) for k in keys_by_mode.get("scramble", ())
        )
        mask = identity["sector"].astype(int).isin(sec_set) & identity["key"].astype(str).isin(null_key_set)
        return identity.loc[mask].reset_index(drop=True)
    raise SystemExit(f"unsupported common_period_mode: {common_period_mode!r}")


def _common_period_population_counts(
    identity: pd.DataFrame,
    common_period_mode: str,
    keys_by_mode: Optional[Mapping[str, Sequence[str]]] = None,
    sectors: Sequence[int] = CALIBRATION_SECTORS,
) -> Dict[str, Any]:
    """Counts-only per (sector, task) breakdown of the selected common_period population (identity-only)."""
    pop_id = _common_period_population_identity(
        identity,
        common_period_mode=common_period_mode,
        keys_by_mode=keys_by_mode,
        sectors=sectors,
    )
    group_col = "task" if "task" in pop_id.columns else "mode"
    per_sector_task: Dict[str, Dict[str, int]] = {}
    for (sec, task), g in pop_id.groupby(["sector", group_col], sort=True):
        per_sector_task.setdefault(str(int(sec)), {})[str(task)] = int(len(g))
    n_holdout = int(pop_id["holdout"].astype(bool).sum()) if len(pop_id) else 0
    n_total = int(len(pop_id))
    return {
        "mode": common_period_mode,
        "population": COMMON_PERIOD_POPULATION_DESC[common_period_mode],
        "n_total": n_total,
        "n_holdout": n_holdout,
        "n_non_holdout": n_total - n_holdout,
        "per_sector_task": per_sector_task,
    }


def _summarize_common_period_groups(df: pd.DataFrame) -> Dict[str, Dict[str, Dict[str, int]]]:
    """Counts-only per (sector, task) group size and number flagged for common_period (non-holdout null rows only)."""
    if df.empty:
        return {}
    if "mode" not in df.columns or "holdout" not in df.columns:
        raise RuntimeError("_summarize_common_period_groups requires 'mode' and 'holdout' columns")
    sub = df[df["mode"].astype(str).isin(NULL_MODES) & (~df["holdout"].astype(bool))]
    if sub.empty:
        return {}
    group_col = "task" if "task" in sub.columns else "mode"
    per_sector_task: Dict[str, Dict[str, Dict[str, int]]] = {}
    for (sec, task), g in sub.groupby(["sector", group_col], sort=True):
        per_sector_task.setdefault(str(int(sec)), {})[str(task)] = {
            "n_rows": int(len(g)),
            "n_flagged": int(g["fail_common_period"].astype(bool).sum()),
        }
    return per_sector_task


def _apply_common_period_flags(
    sci: pd.DataFrame,
    common_period_mode: str,
    rows_dir: Union[str, Path],
    identity: pd.DataFrame,
    finish_mod: Any,
    sectors: Sequence[int] = CALIBRATION_SECTORS,
) -> Tuple[pd.DataFrame, Dict[str, Any]]:
    """Compute common_period flags on sci according to common_period_mode."""
    if common_period_mode in ("v11_task_with_holdout", "v11_task_no_holdout"):
        allow_holdout = common_period_mode == "v11_task_with_holdout"
        pop_id = _common_period_population_identity(
            identity,
            common_period_mode=common_period_mode,
            sectors=sectors,
        )
        n_holdout = int(pop_id["holdout"].astype(bool).sum()) if len(pop_id) else 0
        n_total = int(len(pop_id))
        pop_keys = sorted(pop_id["key"].astype(str).unique())
        pop_df = load_common_period_rows_for_keys(
            rows_dir,
            pop_keys,
            identity=identity,
            sectors=sectors,
            allowed_modes=NULL_MODES,
            allow_holdout=allow_holdout,
        )
        flagged_pop = finish_mod.flag_common_periods_by_task(pop_df)
        if "holdout" not in flagged_pop.columns:
            holdout_map = dict(zip(pop_id["key"].astype(str), pop_id["holdout"].astype(bool)))
            flagged_pop = flagged_pop.assign(
                holdout=flagged_pop["key"].astype(str).map(holdout_map).astype(bool)
            )
        pop_key_s = flagged_pop["key"].astype(str)
        flags_map: Dict[str, Any] = dict(zip(pop_key_s, flagged_pop["flags"]))
        fail_map: Dict[str, bool] = dict(zip(pop_key_s, flagged_pop["fail_common_period"].astype(bool)))

        sci_keys = sci["key"].astype(str)
        n_missing = int((~sci_keys.isin(fail_map)).sum())
        if n_missing:
            raise SystemExit(f"common_period {common_period_mode} mapping missing {n_missing} null-sample key(s)")
        sci = sci.copy()
        sci["flags"] = sci_keys.map(flags_map)
        sci["fail_common_period"] = sci_keys.map(fail_map).astype(bool)
        report: Dict[str, Any] = {
            "mode": common_period_mode,
            "population": COMMON_PERIOD_POPULATION_DESC[common_period_mode],
            "n_total": n_total,
            "n_holdout": n_holdout,
            "n_non_holdout": n_total - n_holdout,
            "per_sector_task": _summarize_common_period_groups(flagged_pop),
            "null_sample_flagged_by_mode": {
                m: int(sci.loc[sci["mode"].astype(str) == m, "fail_common_period"].astype(bool).sum())
                for m in NULL_MODES
            },
        }
        return sci, report

    if common_period_mode == "null_sample":
        sci = _flag_common_periods(sci)
        report = {
            "mode": "null_sample",
            "population": COMMON_PERIOD_POPULATION_DESC["null_sample"],
            "n_total": int(len(sci)),
            "n_holdout": 0,
            "n_non_holdout": int(len(sci)),
            "per_sector_task": _summarize_common_period_groups(sci),
            "null_sample_flagged_by_mode": {
                m: int(sci.loc[sci["mode"].astype(str) == m, "fail_common_period"].astype(bool).sum())
                for m in NULL_MODES
            },
        }
        return sci, report

    raise SystemExit(f"unsupported common_period_mode: {common_period_mode!r}")


def _rebuild_for_row_strict(
    row: Mapping[str, Any],
    job: Mapping[str, Any],
    cfg: Mapping[str, Any],
    finish_mod: Any,
) -> Tuple[Optional[Dict[str, Any]], Dict[str, Any], Any]:  # review: B1
    """Strict light-curve rebuild following finish.rebuild_for_row's mask-on branch with NO try/except."""
    if vet_final.odd_even_mask_enabled(cfg) is not True:
        raise SystemExit("vet_final.odd_even_mask_enabled(cfg) must be True for strict rebuild")
    t, f, e, dumps = finish_mod._rebuild_flat_dumps(row, dict(job, cfg=cfg))
    gaps = vet_final.gap_edges(t, vet_final.odd_even_mask_cfg(cfg)["gap_d"])
    lc = dict(flat=(t, f, e), dumps=dumps, gap_edges=gaps)
    scan = (
        finish_mod.scan_from_flat(row, (t, f, e))
        if "secondary_scan" in vet_amend._set(row.get("deep_flags"))
        else None
    )
    return scan, lc, dumps


def _is_scan_consistency_error(exc: BaseException) -> bool:
    """True when RuntimeError originated from finish.scan_from_flat / vet_final.check_scan_sigma."""
    if "rebuilt secondary scan" in str(exc):
        return True
    tb = exc.__traceback__
    while tb is not None:
        if tb.tb_frame.f_code.co_name in ("scan_from_flat", "check_scan_sigma"):
            return True
        tb = tb.tb_next
    return False


def _stamp_and_write_vetted_csv(
    records: Sequence[Mapping[str, Any]],
    cfg: Mapping[str, Any],
    path: Path,
) -> pd.DataFrame:  # review: B5, S8
    """Stamp vetted rows with v1.2 settings hash & cap and write CSV."""
    digest, cap = v12_settings.stamp_values(cfg)
    df = pd.DataFrame(list(records), columns=list(VETTED_CSV_COLS))
    stamped = v12_settings.stamp_v12(df, digest, cap)
    path.parent.mkdir(parents=True, exist_ok=True)
    stamped.to_csv(path, index=False)
    return stamped


def _vet_mode_top_down(
    mode: str,
    mode_df: pd.DataFrame,
    need_m: int,
    gate_key_set: Set[str],
    jobs_by_sector_tic: Mapping[int, Mapping[int, Mapping[str, Any]]],
    cfg: Mapping[str, Any],
    finish_mod: Any,
    out_dir: Path,
) -> Dict[str, Any]:  # review: S7, S8
    """Walk mode's null sample in order S desc, key asc; vet until need_m pass; enforce gate & all-hard LC."""
    ordered = _sort_null_frame(mode_df)
    records: List[Dict[str, Any]] = []
    n_vetted = 0
    n_passed = 0
    last_rank = 0
    last_s: Optional[float] = None
    all_vetted = True
    partial_path = out_dir / f"vetted_{mode}.partial.csv"

    for rank_idx, row in enumerate(ordered.to_dict(orient="records"), start=1):
        if n_passed >= need_m:
            all_vetted = False
            break
        key = str(row["key"])
        if key not in gate_key_set:
            _stamp_and_write_vetted_csv(records, cfg, partial_path)
            raise GateDepthExceeded(
                {
                    "error": "gate depth exceeded",
                    "mode": mode,
                    "n_vetted": n_vetted,
                    "n_passed": n_passed,
                    "need": int(need_m),
                    "rank_reached": rank_idx,
                }
            )
        sec = int(row["sector"])
        tic = int(row["tic"])
        jobs_sec = jobs_by_sector_tic.get(sec, {})
        job = jobs_sec.get(tic) or {}
        try:
            sc, lc, dumps = _rebuild_for_row_strict(row, job, cfg, finish_mod)
        except SystemExit:
            raise
        except RuntimeError as exc:
            _stamp_and_write_vetted_csv(records, cfg, partial_path)
            if _is_scan_consistency_error(exc):
                raise ScanConsistencyError(
                    {
                        "error": "scan_consistency",
                        "mode": mode,
                        "sector": sec,
                        "rank": rank_idx,
                        "n_vetted": n_vetted,
                        "n_passed": n_passed,
                        "exception_type": type(exc).__name__,
                    }
                ) from None
            raise RebuildError(
                {
                    "error": "rebuild_error",
                    "mode": mode,
                    "sector": sec,
                    "rank": rank_idx,
                    "n_vetted": n_vetted,
                    "n_passed": n_passed,
                    "exception_type": type(exc).__name__,
                }
            ) from None
        except Exception as exc:
            _stamp_and_write_vetted_csv(records, cfg, partial_path)
            raise RebuildError(
                {
                    "error": "rebuild_error",
                    "mode": mode,
                    "sector": sec,
                    "rank": rank_idx,
                    "n_vetted": n_vetted,
                    "n_passed": n_passed,
                    "exception_type": type(exc).__name__,
                }
            ) from None

        try:
            verdict = finish_mod.row_verdict(
                row,
                jobs=jobs_sec,
                cfg=cfg,
                fail_common_period=bool(row["fail_common_period"]),
                scan=sc,
                lc=lc,
                dumps=dumps,
                needs_lc="raise",
            )
        except ValueError as exc:
            msg = str(exc).lower()
            if "no light curve" in msg or "requires a rebuilt scan" in msg:
                _stamp_and_write_vetted_csv(records, cfg, partial_path)
                raise NeedsLightCurve(
                    {
                        "error": "needs_lc",
                        "mode": mode,
                        "sector": sec,
                        "rank": rank_idx,
                        "n_vetted": n_vetted,
                        "n_passed": n_passed,
                        "exception_type": type(exc).__name__,
                    }
                ) from None
            raise
        except RuntimeError as exc:
            if _is_scan_consistency_error(exc):
                _stamp_and_write_vetted_csv(records, cfg, partial_path)
                raise ScanConsistencyError(
                    {
                        "error": "scan_consistency",
                        "mode": mode,
                        "sector": sec,
                        "rank": rank_idx,
                        "n_vetted": n_vetted,
                        "n_passed": n_passed,
                        "exception_type": type(exc).__name__,
                    }
                ) from None
            raise
        if verdict.status == vet_final.STATUS_NEEDS_LC:
            _stamp_and_write_vetted_csv(records, cfg, partial_path)
            raise NeedsLightCurve(
                {
                    "error": "needs_lc",
                    "mode": mode,
                    "sector": sec,
                    "rank": rank_idx,
                    "n_vetted": n_vetted,
                    "n_passed": n_passed,
                    "status": str(verdict.status),
                }
            )

        td_flag = bool(vet.too_deep(row.get("depth"), cfg)["flagged"])
        s_val = float(row["S"])
        passed = bool(verdict.passed)
        records.append(
            {
                "key": key,
                "sector": sec,
                "S": s_val,
                "passed": passed,
                "first_failing_rule": verdict.first_failing_rule or "",
                "all_failing": ",".join(verdict.all_failing),
                "status": str(verdict.status),
                "too_deep_v2": td_flag,
            }
        )
        n_vetted += 1
        last_rank = rank_idx
        last_s = s_val
        if passed:
            n_passed += 1

    vetted_df = _stamp_and_write_vetted_csv(records, cfg, out_dir / f"vetted_{mode}.csv")
    return {
        "records": records,
        "vetted_df": vetted_df,
        "n_vetted": n_vetted,
        "n_passed": n_passed,
        "last_vetted_rank": last_rank,
        "last_vetted_S": last_s,
        "all_vetted": all_vetted,
    }


def _per_sector_sstar_report(
    mode_df: pd.DataFrame,
    vet_res: Mapping[str, Any],
    alpha: Mapping[int, Mapping[str, Any]],
    cfg: Mapping[str, Any],
    sectors: Sequence[int] = CALIBRATION_SECTORS,
) -> Dict[str, Any]:
    """Per-sector report-only S* (with determined rule) and 2000-resample 68% bootstrap interval."""
    last_s = vet_res["last_vetted_S"]
    all_vetted = bool(vet_res["all_vetted"])
    passed_keys_ge_cut: Set[str] = {
        str(r["key"])
        for r in vet_res["records"]
        if r["passed"] and (all_vetted or last_s is None or float(r["S"]) >= float(last_s))
    }

    per_sector: Dict[str, Any] = {}
    for s in sorted(int(x) for x in sectors):
        sub = mode_df[mode_df["sector"].astype(int) == s]
        n_ms = int(len(sub))
        p_s = float(alpha[s]["p_s"])
        b_ms = n_ms * p_s
        need_ms = int(math.floor(b_ms + _BUDGET_EPS)) + 1

        keys_arr = sub["key"].astype(str).to_numpy()
        s_arr = sub["S"].to_numpy(dtype=float)
        pass_mask = np.array([k in passed_keys_ge_cut for k in keys_arr], dtype=bool)
        s_pass = s_arr[pass_mask]
        determined = bool(all_vetted or (len(s_pass) >= need_ms))

        pt: Optional[Dict[str, Any]] = None
        pt_error: Optional[str] = None
        if determined:
            try:
                pt = v12_settings.calibrate_sstar_v12(s_pass, n_null=n_ms, p=p_s, u=None)
            except sstar_grid.SStarGridError as exc:
                pt_error = str(exc)

        rng = np.random.default_rng(BOOT_SEED)
        boots: List[float] = []
        n_undeterminable = 0
        n_grid_fail = 0
        if n_ms > 0:
            for _ in range(BOOT_N):
                idx = rng.integers(0, n_ms, n_ms)
                bm = pass_mask[idx]
                if (not all_vetted) and int(bm.sum()) < need_ms:
                    n_undeterminable += 1
                    continue
                try:
                    res = sstar_grid.select_sstar(s_arr[idx][bm], n_ms, p_s, None, cfg=cfg)
                    boots.append(float(res["s_star"]))
                except sstar_grid.SStarGridError:
                    n_grid_fail += 1
        else:
            n_undeterminable = BOOT_N

        ci68 = [float(np.percentile(boots, 16)), float(np.percentile(boots, 84))] if boots else None
        per_sector[str(s)] = {
            "n_null": n_ms,
            "p_s": p_s,
            "B": float(b_ms),
            "need": int(need_ms),
            "n_vetted_pass_above_cut": int(len(s_pass)),
            "determined": determined,
            "sstar": pt,
            "sstar_error": pt_error,
            "boot_n": BOOT_N,
            "boot_seed": BOOT_SEED,
            "boot_determinable": len(boots),
            "boot_undeterminable": n_undeterminable,
            "boot_grid_fail": n_grid_fail,
            "ci68": ci68,
            "report_only": True,
        }
    return per_sector


# review: S6 (n_real_eligible removed)
def run_compute(
    args: argparse.Namespace,
    header: Mapping[str, Any],
    jobs: Sequence[Mapping[str, Any]],
    identity: pd.DataFrame,
    alpha: Mapping[int, Mapping[str, Any]],
    keys_by_mode: Mapping[str, Sequence[str]],
    null_rep: Mapping[str, Any],
    base_report: Dict[str, Any],
    dirty_count: int = 0,
    quality_sha256s: Optional[Mapping[int, str]] = None,
) -> int:
    """Execute --compute: verify gate_pass.json, vet top-down per mode, and calibrate S*_cal."""
    kernels.install()  # review: B2
    out_p = Path(args.out_dir)
    gate_pass_path = out_p / "gate_pass.json"
    if not gate_pass_path.is_file():
        raise SystemExit(f"missing required gate pass file: {gate_pass_path}")
    gate_pass = json.loads(gate_pass_path.read_text(encoding="utf-8"))

    # review: B4
    if gate_pass.get("passed") is not True:
        raise SystemExit(f"gate_pass.json field 'passed' differs: expected True, got {gate_pass.get('passed')!r}")
    cur_commit = str(settings.commit())
    if str(gate_pass.get("pinned_commit")) != cur_commit:
        raise SystemExit(
            f"gate_pass.json field 'pinned_commit' differs: {gate_pass.get('pinned_commit')!r} != {cur_commit!r}"
        )
    got_s_tol = gate_pass.get("s_tol")
    if got_s_tol is None or float(got_s_tol) != float(GATE_S_TOL):
        raise SystemExit(f"gate_pass.json field 's_tol' differs: {got_s_tol!r} != {GATE_S_TOL!r}")
    cur_driver_sha = _driver_sha256()
    if str(gate_pass.get("driver_sha256")) != cur_driver_sha:
        raise SystemExit(
            f"gate_pass.json field 'driver_sha256' differs: {gate_pass.get('driver_sha256')!r} != {cur_driver_sha!r}"
        )
    cur_manifest_sha = _sha256_bytes(Path(args.manifest).read_bytes())
    if "manifest_sha256" not in gate_pass:
        raise SystemExit("gate_pass.json missing required field 'manifest_sha256'")
    if str(gate_pass.get("manifest_sha256")) != cur_manifest_sha:
        raise SystemExit(
            f"gate_pass.json field 'manifest_sha256' differs: "
            f"{gate_pass.get('manifest_sha256')!r} != {cur_manifest_sha!r}"
        )
    if "gate_top" not in gate_pass:
        raise SystemExit("gate_pass.json missing required field 'gate_top'")
    if "sha256" not in gate_pass:
        raise SystemExit("gate_pass.json missing required field 'sha256'")

    gate_top_used = int(gate_pass["gate_top"])
    gate_info = build_gate_keys(identity, args.rows_dir, keys_by_mode, gate_top_used, CALIBRATION_SECTORS)
    if gate_info["sha256"] != str(gate_pass["sha256"]):
        raise SystemExit(
            f"gate_pass.json field 'sha256' differs: recomputed {gate_info['sha256']} != gate_pass {gate_pass['sha256']}"
        )
    _precheck_gated_light_curves(gate_info["gate_keys"], identity, jobs)  # review: S1

    gate_key_set = set(gate_info["gate_keys"])
    gate_depth_by_mode = {
        "invert": len(gate_info["top_invert"]),
        "scramble": len(gate_info["top_scramble"]),
    }

    cfg = v12_settings.require_v12_settings()
    if cfg.get("vet", {}).get("too_deep", {}).get("mode") != "flag":
        raise SystemExit("frozen v1.2 config must have vet.too_deep.mode == 'flag'")
    if vet_final.odd_even_mask_enabled(cfg) is not True:
        raise SystemExit("frozen v1.2 config must have vet.odd_even_mask.enabled == True")
    finish_mod = _load_pinned_finish()

    cfg_inj = header.get("cfg_inj")
    jobs_by_sector_tic: Dict[int, Dict[int, Dict[str, Any]]] = {}
    for j in jobs:
        sec = int(j["sector"])
        tic = int(j["tic"])
        jc = dict(j)
        if cfg_inj is not None and "cfg_inj" not in jc:
            jc["cfg_inj"] = cfg_inj
        jobs_by_sector_tic.setdefault(sec, {})[tic] = jc

    all_null_keys = sorted(set(keys_by_mode["invert"]) | set(keys_by_mode["scramble"]))
    sci = load_science_for_keys(args.rows_dir, all_null_keys, CALIBRATION_SECTORS, allowed_modes=NULL_MODES)

    # Compute common_period flags BEFORE ranking (_sort_null_frame) and the top-down vetting walk.
    sci, common_period_report = _apply_common_period_flags(
        sci=sci,
        common_period_mode=str(args.common_period_mode),
        rows_dir=args.rows_dir,
        identity=identity,
        finish_mod=finish_mod,
        sectors=CALIBRATION_SECTORS,
    )

    per_mode_report: Dict[str, Any] = {}
    per_sector_table: Dict[str, Any] = {}
    sstar_by_mode: Dict[str, Dict[str, Any]] = {}
    vetted_dfs_by_mode: Dict[str, pd.DataFrame] = {}

    for mode in NULL_MODES:
        m_info = null_rep["by_mode"][mode]
        mode_df = sci[sci["mode"].astype(str) == mode].reset_index(drop=True)
        try:
            vet_res = _vet_mode_top_down(
                mode=mode,
                mode_df=mode_df,
                need_m=int(m_info["need"]),
                gate_key_set=gate_key_set,
                jobs_by_sector_tic=jobs_by_sector_tic,
                cfg=cfg,
                finish_mod=finish_mod,
                out_dir=out_p,
            )
        except GateDepthExceeded as exc:
            err_rep = {
                **base_report,
                "common_period": common_period_report,
                "status": "gate_depth_exceeded",
                "failure": exc.report,
            }
            payload = _write_json(out_p / "report_compute.json", err_rep)
            print(payload, end="")
            return 4
        except NeedsLightCurve as exc:
            err_rep = {
                **base_report,
                "common_period": common_period_report,
                "status": "needs_lc",
                "failure": exc.report,
            }
            payload = _write_json(out_p / "report_compute.json", err_rep)
            print(payload, end="")
            return 5
        except RebuildError as exc:  # review: S7
            err_rep = {
                **base_report,
                "common_period": common_period_report,
                "status": "rebuild_error",
                "failure": exc.report,
            }
            payload = _write_json(out_p / "report_compute.json", err_rep)
            print(payload, end="")
            return 6
        except ScanConsistencyError as exc:  # review: S7
            err_rep = {
                **base_report,
                "common_period": common_period_report,
                "status": "scan_consistency",
                "failure": exc.report,
            }
            payload = _write_json(out_p / "report_compute.json", err_rep)
            print(payload, end="")
            return 7

        vetted_df = vet_res["vetted_df"]
        vetted_dfs_by_mode[mode] = vetted_df
        s_pass = [float(r["S"]) for r in vet_res["records"] if r["passed"]]
        sstar_m = v12_settings.calibrate_sstar_v12(  # review: B5
            S_vetted=s_pass,
            n_null=int(m_info["n_null"]),
            p=float(m_info["p"]),
            u=None,
            rows=vetted_df,
        )
        sstar_by_mode[mode] = sstar_m

        try:
            unvetted_m: Dict[str, Any] = v12_settings.calibrate_sstar_v12(
                S_vetted=mode_df["S"].to_numpy(dtype=float),
                n_null=int(m_info["n_null"]),
                p=float(m_info["p"]),
                u=None,
            )
        except sstar_grid.SStarGridError as exc:
            unvetted_m = {"error": str(exc)}

        per_sector_table[mode] = _per_sector_sstar_report(mode_df, vet_res, alpha, cfg, CALIBRATION_SECTORS)
        per_mode_report[mode] = {
            "n_null": int(m_info["n_null"]),
            "B": float(m_info["B"]),
            "p": float(m_info["p"]),
            "need": int(m_info["need"]),
            "n_vetted": int(vet_res["n_vetted"]),
            "n_passed": int(vet_res["n_passed"]),
            "last_vetted_rank": int(vet_res["last_vetted_rank"]),
            "last_vetted_S": vet_res["last_vetted_S"],
            "gate_depth_available": int(gate_depth_by_mode[mode]),
            "sstar": sstar_m,
            "unvetted_diagnostic": unvetted_m,
        }

    s_inv = float(sstar_by_mode["invert"]["s_star"])
    s_scr = float(sstar_by_mode["scramble"]["s_star"])
    headline_s_star = max(s_inv, s_scr)
    headline_mode = "invert" if s_inv >= s_scr else "scramble"
    headline = {
        "s_star": headline_s_star,
        "mode": headline_mode,
        "tie": s_inv == s_scr,
        "s_star_invert": s_inv,
        "s_star_scramble": s_scr,
    }
    lock_res = v12_settings.lock_sstar_v12(  # review: B5
        rows=vetted_dfs_by_mode[headline_mode],
        s_star=headline_s_star,
    )

    cat_check = base_report.get("catalogue_check")
    if not cat_check or cat_check.get("status") != "ok":
        raise SystemExit("missing verified catalogue_check for provenance")
    if quality_sha256s is None or any(s not in quality_sha256s for s in CALIBRATION_SECTORS):
        raise SystemExit("missing verified quality_sha256s for provenance")
    q_hashes = {str(int(s)): str(quality_sha256s[s]) for s in CALIBRATION_SECTORS}
    pgh_sha = base_report.get("p_grid_halfwidth_check", {}).get("got_sha256")
    if not pgh_sha:
        raise SystemExit("missing verified p_grid_halfwidth sha256 for provenance")
    input_hashes_sha = base_report.get("input_hashes", {}).get("sha256")
    if not input_hashes_sha:
        raise SystemExit("missing verified input_hashes sha256 for provenance")
    kostov_sha = str(cat_check["kostov"]["got_sha256"])
    eb_csv_sha = str(base_report.get("exclusion_check", {}).get("eb_dropped_csv", {}).get("csv_sha256", ""))
    provenance = {  # review: S4
        "pinned_commit": cur_commit,
        "driver_git_blob_sha1": _driver_git_blob_sha1(),
        "dirty_files_count": int(dirty_count),
        "driver_sha256": cur_driver_sha,
        "v12_settings_hash": str(sstar_by_mode["invert"]["v12_settings_hash"]),
        "manifest_sha256": cur_manifest_sha,
        "input_hashes_sha256": str(input_hashes_sha),
        "rows_parts_count": len(_part_files(args.rows_dir)),
        "quality_table_sha256s": q_hashes,
        "catalogue_sha256s": {
            "pscomppars": str(cat_check["pscomppars"]["got_sha256"]),
            "toi": str(cat_check["toi"]["got_sha256"]),
            "kostov": kostov_sha,
        },
        "p_grid_halfwidth_sha256": str(pgh_sha),
        "kostov_csv_sha256": kostov_sha,
        "eb_list_csv_sha256": eb_csv_sha,
        "gate_key_set_sha256": str(gate_info["sha256"]),
        "common_period_mode": str(args.common_period_mode),
        "ruled_common_period_mode": args.ruled,
        "created_utc": datetime.now(timezone.utc).isoformat(),
    }

    result = {
        **base_report,
        "gate_pass_sha256": gate_info["sha256"],
        "gate_top": gate_top_used,
        "common_period": common_period_report,
        "by_mode": per_mode_report,
        "headline": headline,
        "lock": lock_res,
        "per_sector_report_only": per_sector_table,
        "provenance": provenance,
        "v12_settings_hash": sstar_by_mode["invert"]["v12_settings_hash"],
        "v12_cap_rearth": sstar_by_mode["invert"]["v12_cap_rearth"],
    }

    out_json = out_p / "sstar_cal.json"
    payload = _write_json(out_json, result)
    digest = _sha256_bytes(out_json.read_bytes())
    (out_p / "sstar_cal.json.sha256").write_text(f"{digest}  {out_json.name}\n", encoding="utf-8")
    _write_json(out_p / "report_compute.json", result)
    print(payload, end="")
    return 0


# --------------------------------------------------------------------------- CLI
def _parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    ap = argparse.ArgumentParser(
        description="v1.2 S*_cal driver over S59-S64 (docs/registration_v1.2.md, S*_cal rulings, Deviation 2026-11-20)"
    )
    ap.add_argument("--rows-dir", required=True, type=Path, help="Directory containing parts/*.parquet")
    ap.add_argument("--manifest", required=True, type=Path, help="Gzipped JSONL run manifest")
    ap.add_argument("--out-dir", required=True, type=Path, help="Output directory")
    ap.add_argument("--quality-dir", required=True, type=Path, help="Directory holding star_quality_s00XX.parquet")
    ap.add_argument(
        "--input-hashes",
        required=True,
        type=Path,
        help="Pinned sha256 manifest file for star-quality tables and frozen inputs",
    )
    ap.add_argument("--workers", type=int, default=16, help="Worker processes for --gate (default: 16)")
    ap.add_argument("--gate-top", type=int, default=300, help="Top-K null rows per mode in --gate (default: 300)")
    ap.add_argument(
        "--common-period-mode",
        required=True,
        choices=COMMON_PERIOD_MODES,
        help="Common-period flagging population ('v11_task_with_holdout', 'v11_task_no_holdout', or 'null_sample')",
    )
    ap.add_argument(
        "--ruled",
        choices=COMMON_PERIOD_MODES,
        default=None,
        help=(
            "Ruled common-period mode ('v11_task_with_holdout', 'v11_task_no_holdout', or 'null_sample'); "
            "--gate and --compute require both --common-period-mode and --ruled to equal 'v11_task_with_holdout'"
        ),
    )
    mode = ap.add_mutually_exclusive_group(required=True)
    mode.add_argument("--dry-run", action="store_true", help="Identity-only checks and counts (no science columns)")
    mode.add_argument("--gate", action="store_true", help="Reproduction gate on known-planet real + top/random nulls")
    mode.add_argument("--compute", action="store_true", help="Top-down vetting and S*_cal calibration")
    return ap.parse_args(argv)


def _validate_common_period_args(args: argparse.Namespace) -> None:
    cp_mode = args.common_period_mode
    ruled = args.ruled
    if args.gate or args.compute:
        if cp_mode != RULED_COMMON_PERIOD_MODE or ruled != RULED_COMMON_PERIOD_MODE:
            raise SystemExit(
                f"--gate/--compute require --common-period-mode={RULED_COMMON_PERIOD_MODE!r} "
                f"and --ruled={RULED_COMMON_PERIOD_MODE!r}; "
                f"got common_period_mode={cp_mode!r}, ruled={ruled!r}"
            )
    elif args.dry_run:
        if ruled is not None and ruled != cp_mode:
            raise SystemExit(
                f"--dry-run requires --ruled to match --common-period-mode when passed; "
                f"got common_period_mode={cp_mode!r}, ruled={ruled!r}"
            )


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = _parse_args(argv)
    _validate_common_period_args(args)
    run_mode = "dry_run" if args.dry_run else ("gate" if args.gate else "compute")
    out_p = Path(args.out_dir)
    out_p.mkdir(parents=True, exist_ok=True)

    if args.gate:
        (out_p / "gate_pass.json").unlink(missing_ok=True)  # review: B4

    dirty_count = 0
    if args.gate or args.compute:
        dirty_files_list = worker_guard.require_clean(False)  # review: S4
        dirty_count = len(dirty_files_list)

    base_report: Dict[str, Any] = {
        "mode": run_mode,
        "driver_git_blob_sha1": _driver_git_blob_sha1(),
        "sectors": list(CALIBRATION_SECTORS),
        "common_period_mode": str(args.common_period_mode),
        "ruled_common_period_mode": args.ruled,
    }

    verified_quality_shas, input_hashes_report = check_input_hashes(
        args.input_hashes, args.quality_dir, CALIBRATION_SECTORS
    )
    base_report["input_hashes"] = input_hashes_report

    cat_report = check_catalogue_hashes()  # review: S2
    base_report["catalogue_check"] = cat_report
    if cat_report["status"] != "ok":
        base_report["status"] = "catalogue_hash_mismatch"
        payload = _write_json(out_p / f"report_{run_mode}.json", base_report)
        print(payload, end="")
        return 8

    try:
        pgh_report = check_p_grid_halfwidth(require=True)  # review: S3
    except MissingInput as exc:
        base_report["p_grid_halfwidth_present"] = False
        base_report["p_grid_halfwidth_check"] = {
            "present": False,
            "ok": False,
            "expected_sha256": PINNED_PGRID_SHA256,
            "got_sha256": None,
            "error": str(exc),
        }
        base_report["status"] = "missing_p_grid_halfwidth"
        payload = _write_json(out_p / f"report_{run_mode}.json", base_report)
        print(payload, end="")
        return 9
    base_report["p_grid_halfwidth_present"] = bool(pgh_report["present"])
    base_report["p_grid_halfwidth_check"] = pgh_report

    header, jobs = load_manifest(args.manifest)
    if str(header.get("pinned_commit")) != PINNED_CODE_COMMIT:
        raise SystemExit(
            f"manifest pinned_commit mismatch: expected {PINNED_CODE_COMMIT!r}, "
            f"got {header.get('pinned_commit')!r}"
        )
    checkout_report: Optional[Dict[str, Any]] = None
    if args.gate or args.compute:          # pinned checkout and settings for any recompute or vetting
        checkout_report = check_checkout_and_cfg(header)

    identity = load_identity(args.rows_dir, CALIBRATION_SECTORS)
    manifest_excl = rebuild_exclusion(jobs, CALIBRATION_SECTORS)
    reg_sets, kostov_sets, reg_counts = registered_exclusion(jobs, CALIBRATION_SECTORS)

    try:
        ex_res = check_exclusion_with_eb_drop(manifest_excl, reg_sets, kostov_sets, reg_counts, EXPECTED_EB_DROP)
    except ExclusionMismatch as exc:
        base_report["exclusion_check"] = exc.report
        payload = _write_json(out_p / f"report_{run_mode}.json", base_report)
        print(payload, end="")
        return 2

    identity, eb_drop_info = apply_eb_drop(identity, ex_res["eb_drop"], out_p)
    base_report["exclusion_check"] = {**ex_res["report"], "eb_dropped_csv": eb_drop_info}

    try:
        tw_report, status_err_report = check_three_way_match(identity, ex_res["registered"], CALIBRATION_SECTORS)
    except ThreeWayMismatch as exc:
        base_report.update(exc.report)
        payload = _write_json(out_p / f"report_{run_mode}.json", base_report)
        print(payload, end="")
        return 2

    base_report["three_way_match"] = tw_report
    base_report["status_errors"] = status_err_report

    q = quality_drop(args.quality_dir, verified_quality_shas, CALIBRATION_SECTORS)
    alpha = alpha_n(identity, q)
    keys_by_mode, _, null_rep = null_sample(identity, q, ex_res["effective"], alpha)

    base_report["alpha"] = {str(s): alpha[s] for s in sorted(alpha)}
    base_report["null_sample"] = null_rep

    if args.dry_run:
        n_inv = int(null_rep["by_mode"]["invert"]["n_null"])
        n_scr = int(null_rep["by_mode"]["scramble"]["n_null"])
        n_union = int(null_rep["total_null_rows_union"])
        base_report["common_period_population"] = _common_period_population_counts(
            identity,
            str(args.common_period_mode),
            keys_by_mode=keys_by_mode,
            sectors=CALIBRATION_SECTORS,
        )
        base_report["gate_set_sizes"] = {
            "gate_top": int(args.gate_top),
            "top_invert": min(int(args.gate_top), n_inv),
            "top_scramble": min(int(args.gate_top), n_scr),
            "random_null": min(GATE_RANDOM_N, n_union),
            "known_planet_real": None,
            "note": "known_planet_real is computed only in --gate",
        }
        payload = _write_json(out_p / "report_dry_run.json", base_report)
        print(payload, end="")
        return 0

    if args.gate:
        assert checkout_report is not None
        return run_gate(args, header, jobs, identity, keys_by_mode, base_report, checkout_report)

    return run_compute(
        args,
        header,
        jobs,
        identity,
        alpha,
        keys_by_mode,
        null_rep,
        base_report,
        dirty_count=dirty_count,
        quality_sha256s=q.attrs.get("sha256s"),
    )


if __name__ == "__main__":
    sys.exit(main())
