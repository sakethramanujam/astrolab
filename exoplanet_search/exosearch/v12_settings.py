"""v1.2 frozen settings: load, validate and hash ``configs/v1.2.json``.

The v1.2 calibration, the S69 known-planet gate and v1.2 rescoring must not silently run with the
v1.1 defaults in ``configs/v1.1.json`` (odd-even mask off, ``too_deep`` as a cut, radius cap off,
S* fine grid off).  Those defaults stay correct for the v1.1 paths, so this module adds distinct
v1.2 entry points that always load the frozen v1.2 file first and refuse to run without it.

Frozen file
-----------
``exoplanet_search/configs/v1.2.json`` (a complete config, not an overlay on ``v1.1.json``) with
``"name": "exosearch_v1.2"`` and these explicit values:

* ``vet.odd_even_mask``: ``enabled`` true, ``min_coverage`` 0.75 (X), ``min_gap_dump_dist_h`` 9.0
  (Y), ``min_gap_dump_dist_d`` 0.375, ``gap_d`` 0.5, ``oot_baseline_t14`` 1.0,
  ``oot_one_side`` true (>= 1 T14 out-of-transit baseline on one side).
* ``vet.too_deep``: ``mode`` ``"flag"`` (recorded, never a rejection) and ``max_depth`` 0.05.
* ``vet.rp_cap``: ``enabled`` true, ``cap_rearth`` 30.0 (the default), ``cap_rearth_candidates``
  ``[22, 25, 30]`` and ``include_e_rad`` true (the Rp - 2 sigma error already includes the TIC
  radius error in quadrature).
* ``sstar``: ``fine_grid_enabled`` true, ``grid_min`` 6.0, ``grid_max`` 30.0, ``grid_step`` 0.005,
  ``interpolate`` false.
* ``needs_lc.rebuild``: true (v1.2 rescoring rebuilds the light curve as in step 5b).

All other keys are copied from ``configs/v1.1.json`` so a v1.2 run has one complete file.

Content hash
------------
``REGISTERED_HASH`` is the sha256 of the canonical JSON of the loaded object:
``json.dumps(obj, sort_keys=True, separators=(",", ":")).encode("utf-8")``.  The hash is *not*
stored inside ``v1.2.json`` (that would be circular).  It is the canonical form that is hashed, not
the pretty-printed on-disk bytes.

Entry points (all call ``require_v12_settings`` first)
------------------------------------------------------
* ``calibrate_sstar_v12`` -- S* calibration on the frozen fine grid (delegates to
  ``sstar_grid.select_sstar`` with the loaded cfg).  An optional ``cap_rearth`` override must be one
  of ``cap_rearth_candidates``.
* ``exosearch.v11.known_planets.evaluate_known_planets_v12`` and
  ``run_known_planet_recovery_v12`` -- the S69 / v1.2 known-planet gate (CLI: ``report-v12`` and
  ``freeze-eligibility-v12``).
* ``runs/v1.1_s70_x86/finish.py``: ``rescore_v12_cfg`` / ``rescore_rows_v12`` -- v1.2 rescoring with
  the step-5b light-curve rebuild path.

Output stamps (step 7)
----------------------
Every v1.2 output row / result dict carries ``v12_settings_hash`` (the canonical content hash of the
frozen file *as loaded*, before any cap override) and ``v12_cap_rearth`` (the cap actually applied:
30.0 by default, 22.0 / 25.0 with a registered override).  ``stamp_v12`` writes them without
mutating the caller, and ``require_v12_stamps`` is the guard every v1.2 aggregator / scorer calls;
it raises :class:`V12SettingsError` for missing, mismatched or mixed stamps.  Search-stage rows are
never stamped or required to be.

The v1.1 functions keep their existing defaults and never import this module.
"""
from __future__ import annotations

import copy
import hashlib
import json
import math
from pathlib import Path
from typing import Any, Dict, Mapping, Optional, Union

ROOT = Path(__file__).resolve().parents[1]
CONFIG_FILE = ROOT / "configs" / "v1.2.json"
DEFAULT_PATH = CONFIG_FILE

# sha256 of the canonical JSON of configs/v1.2.json (see module docstring).
REGISTERED_HASH = "7afb7fc178138d59c84fc509a5b47e09ebc8d855d5d1d8ef69e446ee939f8756"

CAP_REARTH_CANDIDATES = (22.0, 25.0, 30.0)

# Stamps written on every v1.2 *output* row / result (never inside the search-stage rows).
STAMP_HASH_KEY = "v12_settings_hash"
STAMP_CAP_KEY = "v12_cap_rearth"


class V12SettingsError(RuntimeError, ValueError):
    """The frozen v1.2 settings are missing, at v1.1 defaults, or not the registered content."""


def _canonical(obj: Any) -> bytes:
    return json.dumps(obj, sort_keys=True, separators=(",", ":")).encode("utf-8")


def content_hash(path_or_obj: Union[str, Path, Mapping[str, Any]]) -> str:
    """sha256 of the canonical JSON of a config path or an already-loaded config object."""
    obj = path_or_obj if isinstance(path_or_obj, Mapping) else json.loads(Path(path_or_obj).read_text())
    return hashlib.sha256(_canonical(obj)).hexdigest()


def _get(cfg: Any, *keys: str, default: Any = None) -> Any:
    cur = cfg
    for k in keys:
        if not isinstance(cur, Mapping) or k not in cur:
            return default
        cur = cur[k]
    return cur


_MISSING = object()


def _approx(a: Any, b: float, tol: float = 1e-12) -> bool:
    try:
        return math.isclose(float(a), float(b), rel_tol=0.0, abs_tol=tol)
    except (TypeError, ValueError):
        return False


def _check_bool(problems, cfg, keys, expected=True):
    val = _get(cfg, *keys, default=_MISSING)
    name = ".".join(keys)
    if val is _MISSING:
        problems.append(f"{name} is missing (frozen value: {expected})")
    elif val is not expected:
        problems.append(f"{name} is {val!r} (frozen value: {expected}; v1.1 default is {not expected})")


def _check_num(problems, cfg, keys, expected):
    val = _get(cfg, *keys, default=_MISSING)
    name = ".".join(keys)
    if val is _MISSING:
        problems.append(f"{name} is missing (frozen value: {expected})")
    elif not _approx(val, expected):
        problems.append(f"{name} is {val!r} (frozen value: {expected})")


def validate_v12_values(cfg: Mapping[str, Any]) -> Mapping[str, Any]:
    """Refuse a config whose frozen v1.2 keys are missing or still at their v1.1 defaults.

    Returns ``cfg`` on success; raises :class:`V12SettingsError` naming every bad key otherwise.
    The content hash is checked separately by :func:`require_v12_settings` so a v1.1-default file
    fails on its values first.
    """
    problems = []
    if not isinstance(cfg, Mapping):
        raise V12SettingsError(f"v1.2 settings must be a JSON object, got {type(cfg).__name__}")

    _check_bool(problems, cfg, ("vet", "odd_even_mask", "enabled"), True)
    _check_num(problems, cfg, ("vet", "odd_even_mask", "min_coverage"), 0.75)
    _check_num(problems, cfg, ("vet", "odd_even_mask", "min_gap_dump_dist_h"), 9.0)
    _check_num(problems, cfg, ("vet", "odd_even_mask", "min_gap_dump_dist_d"), 0.375)
    _check_num(problems, cfg, ("vet", "odd_even_mask", "gap_d"), 0.5)
    _check_num(problems, cfg, ("vet", "odd_even_mask", "oot_baseline_t14"), 1.0)
    _check_bool(problems, cfg, ("vet", "odd_even_mask", "oot_one_side"), True)

    mode = _get(cfg, "vet", "too_deep", "mode", default=_MISSING)
    if mode is _MISSING:
        problems.append("vet.too_deep.mode is missing (frozen value: 'flag'; v1.1 default is 'cut')")
    elif str(mode).strip().lower() != "flag":
        problems.append(f"vet.too_deep.mode is {mode!r} (frozen value: 'flag'; v1.1 default is 'cut')")
    _check_num(problems, cfg, ("vet", "too_deep", "max_depth"), 0.05)

    _check_bool(problems, cfg, ("vet", "rp_cap", "enabled"), True)
    _check_num(problems, cfg, ("vet", "rp_cap", "cap_rearth"), 30.0)
    cand = _get(cfg, "vet", "rp_cap", "cap_rearth_candidates", default=_MISSING)
    if cand is _MISSING:
        problems.append("vet.rp_cap.cap_rearth_candidates is missing (frozen value: [22, 25, 30])")
    elif not (isinstance(cand, (list, tuple)) and len(cand) == 3
              and all(_approx(x, e) for x, e in zip(cand, CAP_REARTH_CANDIDATES))):
        problems.append(f"vet.rp_cap.cap_rearth_candidates is {cand!r} (frozen value: [22, 25, 30])")
    _check_bool(problems, cfg, ("vet", "rp_cap", "include_e_rad"), True)

    _check_bool(problems, cfg, ("sstar", "fine_grid_enabled"), True)
    _check_num(problems, cfg, ("sstar", "grid_min"), 6.0)
    _check_num(problems, cfg, ("sstar", "grid_max"), 30.0)
    _check_num(problems, cfg, ("sstar", "grid_step"), 0.005)
    _check_bool(problems, cfg, ("sstar", "interpolate"), False)

    _check_bool(problems, cfg, ("needs_lc", "rebuild"), True)

    if problems:
        raise V12SettingsError("v1.2 settings have v1.1-default or wrong frozen values: " + "; ".join(problems))
    return cfg


def require_v12_settings(path: Optional[Union[str, Path]] = None) -> Dict[str, Any]:
    """Load the frozen v1.2 settings or refuse to run.

    Checks, in order: (a) the file exists, (b) every frozen key is set to its v1.2 value
    (:func:`validate_v12_values`), (c) the canonical content hash equals ``REGISTERED_HASH``.
    Returns the loaded config on success.  ``path`` defaults to the committed ``configs/v1.2.json``.
    """
    p = Path(path) if path is not None else DEFAULT_PATH
    if not p.is_file():
        raise V12SettingsError(f"frozen v1.2 settings file required but missing: {p}")
    try:
        cfg = json.loads(p.read_text())
    except (OSError, json.JSONDecodeError) as exc:
        raise V12SettingsError(f"frozen v1.2 settings file {p} is unreadable: {exc}") from exc
    validate_v12_values(cfg)
    digest = content_hash(cfg)
    if digest != REGISTERED_HASH:
        raise V12SettingsError(
            f"frozen v1.2 settings content hash mismatch for {p}: {digest} != registered {REGISTERED_HASH}")
    return cfg


def with_cap_rearth(cfg: Mapping[str, Any], cap: float) -> Dict[str, Any]:
    """Return a copy of ``cfg`` with ``vet.rp_cap.cap_rearth`` set to a registered candidate.

    Only the run-time calibration choice is overridden; the frozen file is not edited or re-hashed.
    """
    cand = _get(cfg, "vet", "rp_cap", "cap_rearth_candidates", default=CAP_REARTH_CANDIDATES)
    try:
        vals = [float(x) for x in cand]
    except (TypeError, ValueError):
        vals = list(CAP_REARTH_CANDIDATES)
    if not any(_approx(cap, x) for x in vals):
        raise V12SettingsError(f"cap_rearth override {cap!r} is not one of the registered candidates {vals}")
    out = copy.deepcopy(dict(cfg))
    out.setdefault("vet", {}).setdefault("rp_cap", {})["cap_rearth"] = float(cap)
    return out


def stamp_values(cfg: Mapping[str, Any], cap_rearth: Optional[float] = None):
    """``(loaded-file content hash, applied cap)`` for a row / result stamp.

    ``cfg`` must be the object loaded by :func:`require_v12_settings` *before* any cap override, so
    the hash always describes the frozen file rather than the run-time cap choice.
    """
    digest = content_hash(cfg)
    if cap_rearth is not None:
        with_cap_rearth(cfg, cap_rearth)        # validates the candidate list
        cap = float(cap_rearth)
    else:
        cap = float(_get(cfg, "vet", "rp_cap", "cap_rearth", default=CAP_REARTH_CANDIDATES[-1]))
    return digest, cap


def _is_dataframe(obj: Any) -> bool:
    return hasattr(obj, "columns") and hasattr(obj, "to_dict") and hasattr(obj, "copy")


def _as_rows(obj: Any):
    if obj is None:
        return []
    if isinstance(obj, Mapping):
        return [dict(obj)] if obj else []
    if isinstance(obj, (list, tuple)):
        return list(obj)
    if _is_dataframe(obj):
        if len(obj) == 0:
            return []
        return obj.to_dict(orient="records")
    raise TypeError(f"v1.2 stamps need a dict, a list of dicts or a DataFrame, got {type(obj).__name__}")


def _blank(value: Any) -> bool:
    if value is None:
        return True
    if isinstance(value, float) and math.isnan(value):
        return True
    if isinstance(value, str):
        return not value.strip()
    try:
        import pandas as pd
        return bool(pd.isna(value))
    except (TypeError, ValueError, ImportError):
        return False


def _cap_value(value: Any) -> Optional[float]:
    if _blank(value):
        return None
    try:
        cap = float(value)
    except (TypeError, ValueError):
        return None
    return cap if math.isfinite(cap) else None


def stamp_v12(obj: Any, settings_hash: str, cap_rearth: float):
    """Return a copy of ``obj`` with the two v1.2 stamp fields on every row (never mutates the caller)."""
    vals = {STAMP_HASH_KEY: str(settings_hash), STAMP_CAP_KEY: float(cap_rearth)}
    if _is_dataframe(obj):
        out = obj.copy()
        for key, value in vals.items():
            out[key] = value
        return out
    if isinstance(obj, Mapping):
        out = dict(obj)
        out.update(vals)
        return out
    if isinstance(obj, (list, tuple)):
        out = []
        for row in obj:
            if not isinstance(row, Mapping):
                raise TypeError(f"v1.2 stamp rows must be dicts, got {type(row).__name__}")
            rec = dict(row)
            rec.update(vals)
            out.append(rec)
        return out
    raise TypeError(f"cannot stamp {type(obj).__name__}; expected dict, list of dicts or DataFrame")


def require_v12_stamps(rows: Any, expected_hash: Optional[str] = None,
                       expected_cap: Optional[float] = None,
                       path: Optional[Union[str, Path]] = None,
                       cap_rearth: Optional[float] = None) -> Any:
    """Refuse v1.2 inputs whose rows lack the stamps, disagree with the loaded file, or mix caps.

    ``rows`` is a DataFrame, a list of dicts or a single dict.  ``expected_hash`` defaults to the
    canonical content hash of the config loaded by :func:`require_v12_settings` (the file as loaded,
    not the ``REGISTERED_HASH`` constant).  Empty input is accepted (vacuously).
    """
    if expected_hash is None:
        expected_hash = content_hash(require_v12_settings(path))
    wanted_caps = []
    if expected_cap is not None:
        wanted_caps.append(float(expected_cap))
    if cap_rearth is not None:
        wanted_caps.append(float(cap_rearth))

    recs = _as_rows(rows)
    if not recs:
        return rows

    hashes, caps = [], []
    for i, row in enumerate(recs):
        if not isinstance(row, Mapping):
            raise V12SettingsError(f"v1.2 output row {i} is not a mapping: {type(row).__name__}")
        h = row.get(STAMP_HASH_KEY)
        c = row.get(STAMP_CAP_KEY)
        missing = [k for k, v in ((STAMP_HASH_KEY, h), (STAMP_CAP_KEY, c)) if _blank(v)]
        if missing:
            raise V12SettingsError(
                f"v1.2 output row {i} lacks stamp field(s) {missing}; refusing to aggregate unstamped rows")
        cval = _cap_value(c)
        if cval is None:
            raise V12SettingsError(f"v1.2 output row {i} has invalid {STAMP_CAP_KEY}={c!r}")
        hashes.append(str(h))
        caps.append(cval)

    unique_hashes = set(hashes)
    if len(unique_hashes) != 1:
        raise V12SettingsError(f"v1.2 output rows mix settings hashes: {sorted(unique_hashes)}")
    if hashes[0] != expected_hash:
        raise V12SettingsError(
            f"v1.2 output rows have settings hash {hashes[0]} != expected loaded-file hash {expected_hash}")

    unique_caps = set(caps)
    if len(unique_caps) != 1:
        raise V12SettingsError(f"v1.2 output rows mix cap_rearth values: {sorted(unique_caps)}")
    for wanted in wanted_caps:
        if not _approx(caps[0], wanted):
            raise V12SettingsError(
                f"v1.2 output rows have cap_rearth={caps[0]} but the aggregator expects {wanted}")
    return rows


def calibrate_sstar_v12(S_vetted, n_null, p, u, cap_rearth: Optional[float] = None,
                        path: Optional[Union[str, Path]] = None, rows=None):
    """v1.2 S* calibration: require the frozen settings, then use its fine-grid ``select_sstar``.

    S* itself does not use the radius cap; a ``cap_rearth`` override is nevertheless checked against
    the registered candidate list here so a later calibration driver cannot select an unregistered cap.
    When ``rows`` (vetted / rescore output) is supplied, its stamps are checked before computing S*.
    The result carries the loaded-file hash and the applied cap.
    """
    cfg = require_v12_settings(path)
    digest, cap = stamp_values(cfg, cap_rearth)
    if cap_rearth is not None:
        cfg = with_cap_rearth(cfg, cap_rearth)
    if rows is not None:
        require_v12_stamps(rows, expected_hash=digest, cap_rearth=cap)
    from .sstar_grid import select_sstar
    return stamp_v12(select_sstar(S_vetted, n_null, p, u, cfg=cfg), digest, cap)


def lock_sstar_v12(rows, s_star, path: Optional[Union[str, Path]] = None,
                   cap_rearth: Optional[float] = None) -> Dict[str, Any]:
    """Lock an S* value only when the v1.2 input rows carry consistent stamps.

    This is the hook a later S* lock driver calls: it requires the frozen settings, refuses rows that
    lack / mismatch / mix stamps, and returns the locked value with both stamp fields.
    """
    cfg = require_v12_settings(path)
    digest, cap = stamp_values(cfg, cap_rearth)
    require_v12_stamps(rows, expected_hash=digest, cap_rearth=cap)
    return stamp_v12({"s_star": float(s_star)}, digest, cap)
