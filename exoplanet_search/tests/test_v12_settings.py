"""v1.2 step 6: frozen settings file, content hash, and the three refusing entry points."""
import copy
import importlib.util
import inspect
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from exosearch import interim
from exosearch import sstar_grid as sg
from exosearch import v12_settings as v12
from exosearch.v11 import known_planets as kp

ROOT = Path(__file__).resolve().parents[1]
V11 = ROOT / "configs" / "v1.1.json"
V12 = ROOT / "configs" / "v1.2.json"

FINISH = ROOT / "runs" / "v1.1_s70_x86" / "finish.py"
_spec = importlib.util.spec_from_file_location("s70_finish_v12", FINISH)
fin = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(fin)


def _sample(seed=7, n=600):
    rng = np.random.default_rng(seed)
    return 6.0 + rng.exponential(1.2, n)


def _load(path):
    return json.loads(Path(path).read_text())


def _set(cfg, keys, value):
    cur = cfg
    for k in keys[:-1]:
        cur = cur[k]
    cur[keys[-1]] = value


def _write(tmp_path, cfg, name="v1.2.json"):
    p = tmp_path / name
    p.write_text(json.dumps(cfg, indent=2) + "\n")
    return p


# ---------------------------------------------------------------- 1. refuses without the file
def test_require_refuses_missing_file(tmp_path):
    missing = tmp_path / "no-settings.json"
    with pytest.raises(v12.V12SettingsError, match="required but missing"):
        v12.require_v12_settings(missing)


def test_all_v12_entry_points_refuse_missing_file(tmp_path):
    missing = tmp_path / "no-settings.json"
    S = _sample()
    with pytest.raises(v12.V12SettingsError, match="required but missing"):
        v12.calibrate_sstar_v12(S, 600, 5.0 / 600, 6.0, path=missing)
    with pytest.raises(v12.V12SettingsError, match="required but missing"):
        kp.evaluate_known_planets_v12(pd.DataFrame(), pd.DataFrame(), 13.71, settings_path=missing)
    with pytest.raises(v12.V12SettingsError, match="required but missing"):
        kp.run_known_planet_recovery_v12(13.71, search_results_path=missing, out_dir=tmp_path / "out",
                                         frozen_eligibility_file=missing, settings_path=missing)
    with pytest.raises(v12.V12SettingsError, match="required but missing"):
        fin.rescore_v12_cfg(missing)
    with pytest.raises(v12.V12SettingsError, match="required but missing"):
        fin.rescore_rows_v12([], path=missing)


# ---------------------------------------------------------------- 2. refuses on hash mismatch
def test_content_hash_is_the_registered_hash():
    assert v12.content_hash(V12) == v12.REGISTERED_HASH
    assert v12.content_hash(_load(V12)) == v12.REGISTERED_HASH


def test_hash_mismatch_other_name_is_refused(tmp_path):
    cfg = _load(V12)
    cfg["name"] = "exosearch_v1.2_unregistered"
    p = _write(tmp_path, cfg)
    assert v12.content_hash(p) != v12.REGISTERED_HASH
    with pytest.raises(v12.V12SettingsError, match="hash"):
        v12.require_v12_settings(p)


def test_hash_mismatch_extra_key_is_refused(tmp_path):
    cfg = _load(V12)
    cfg["innocuous_extra"] = 1
    p = _write(tmp_path, cfg)
    with pytest.raises(v12.V12SettingsError, match="hash"):
        v12.require_v12_settings(p)


# ---------------------------------------------------------------- 3. refuses on v1.1 defaults
def test_v11_file_fails_on_values_not_only_hash(tmp_path):
    cfg = _load(V11)
    with pytest.raises(v12.V12SettingsError) as ei:
        v12.validate_v12_values(cfg)
    msg = str(ei.value)
    for key in ("vet.odd_even_mask.enabled", "vet.too_deep.mode", "vet.rp_cap.enabled",
                "sstar.fine_grid_enabled", "needs_lc.rebuild"):
        assert key in msg
    assert "v1.1-default" in msg or "v1.1 default" in msg

    p = _write(tmp_path, cfg, "v1.1-copied-as-v1.2.json")
    with pytest.raises(v12.V12SettingsError) as ei2:
        v12.require_v12_settings(p)
    # the values gate runs before the hash gate
    assert "vet.odd_even_mask.enabled" in str(ei2.value)
    assert "hash mismatch" not in str(ei2.value)


@pytest.mark.parametrize("keys,value,bad", [
    (("vet", "odd_even_mask", "enabled"), False, "vet.odd_even_mask.enabled"),
    (("vet", "odd_even_mask", "min_coverage"), 0.80, "vet.odd_even_mask.min_coverage"),
    (("vet", "odd_even_mask", "min_gap_dump_dist_h"), 8.0, "vet.odd_even_mask.min_gap_dump_dist_h"),
    (("vet", "odd_even_mask", "min_gap_dump_dist_d"), 0.30, "vet.odd_even_mask.min_gap_dump_dist_d"),
    (("vet", "odd_even_mask", "gap_d"), 0.4, "vet.odd_even_mask.gap_d"),
    (("vet", "odd_even_mask", "oot_baseline_t14"), 0.5, "vet.odd_even_mask.oot_baseline_t14"),
    (("vet", "odd_even_mask", "oot_one_side"), False, "vet.odd_even_mask.oot_one_side"),
    (("vet", "too_deep", "mode"), "cut", "vet.too_deep.mode"),
    (("vet", "too_deep", "max_depth"), 0.06, "vet.too_deep.max_depth"),
    (("vet", "rp_cap", "enabled"), False, "vet.rp_cap.enabled"),
    (("vet", "rp_cap", "cap_rearth"), 25.0, "vet.rp_cap.cap_rearth"),
    (("vet", "rp_cap", "cap_rearth_candidates"), [22, 25, 31], "vet.rp_cap.cap_rearth_candidates"),
    (("vet", "rp_cap", "include_e_rad"), False, "vet.rp_cap.include_e_rad"),
    (("sstar", "fine_grid_enabled"), False, "sstar.fine_grid_enabled"),
    (("sstar", "grid_step"), 0.01, "sstar.grid_step"),
    (("sstar", "interpolate"), True, "sstar.interpolate"),
    (("needs_lc", "rebuild"), False, "needs_lc.rebuild"),
])
def test_individual_v11_default_flips_are_refused(keys, value, bad):
    cfg = copy.deepcopy(_load(V12))
    _set(cfg, keys, value)
    with pytest.raises(v12.V12SettingsError, match=re.escape(bad)):
        v12.validate_v12_values(cfg)


def test_missing_needs_lc_key_is_refused():
    cfg = copy.deepcopy(_load(V12))
    cfg.pop("needs_lc")
    with pytest.raises(v12.V12SettingsError, match=re.escape("needs_lc.rebuild")):
        v12.validate_v12_values(cfg)


def test_missing_include_e_rad_is_refused():
    cfg = copy.deepcopy(_load(V12))
    cfg["vet"]["rp_cap"].pop("include_e_rad")
    with pytest.raises(v12.V12SettingsError, match=re.escape("vet.rp_cap.include_e_rad")):
        v12.validate_v12_values(cfg)


# ---------------------------------------------------------------- 4. loads correctly when valid
def test_committed_file_has_the_frozen_values():
    cfg = v12.require_v12_settings()
    assert cfg["name"] == "exosearch_v1.2"
    mask = cfg["vet"]["odd_even_mask"]
    assert mask["enabled"] is True
    assert mask["min_coverage"] == 0.75
    assert mask["min_gap_dump_dist_h"] == 9.0
    assert mask["min_gap_dump_dist_d"] == 0.375
    assert mask["gap_d"] == 0.5
    assert mask["oot_baseline_t14"] == 1.0
    assert mask["oot_one_side"] is True
    assert cfg["vet"]["too_deep"] == {"mode": "flag", "max_depth": 0.05}
    cap = cfg["vet"]["rp_cap"]
    assert cap["enabled"] is True and cap["cap_rearth"] == 30.0
    assert cap["cap_rearth_candidates"] == [22, 25, 30]
    assert cap["include_e_rad"] is True
    assert cfg["sstar"] == {"fine_grid_enabled": True, "grid_min": 6.0, "grid_max": 30.0,
                            "grid_step": 0.005, "interpolate": False}
    assert cfg["needs_lc"]["rebuild"] is True


def test_calibrate_sstar_v12_uses_the_fine_grid():
    S = _sample()
    r = v12.calibrate_sstar_v12(S, 600, 5.0 / 600, 6.0)
    assert r["rule"] == "grid_discrete"
    assert r["s_star"] == sg.sstar_discrete(S, 5.0)
    assert r["s_star"] in sg.s_grid()


def test_calibrate_sstar_v12_cap_override_only_registered_candidates():
    S = _sample()
    base = v12.calibrate_sstar_v12(S, 600, 5.0 / 600, 6.0)
    for cap in (22, 25, 30):
        r = v12.calibrate_sstar_v12(S, 600, 5.0 / 600, 6.0, cap_rearth=cap)
        assert r["s_star"] == base["s_star"]
    with pytest.raises(v12.V12SettingsError, match="cap_rearth"):
        v12.calibrate_sstar_v12(S, 600, 5.0 / 600, 6.0, cap_rearth=23)


def test_with_cap_rearth_does_not_edit_the_frozen_file_or_object():
    cfg = v12.require_v12_settings()
    out = v12.with_cap_rearth(cfg, 22)
    assert out["vet"]["rp_cap"]["cap_rearth"] == 22.0
    assert cfg["vet"]["rp_cap"]["cap_rearth"] == 30.0
    assert v12.content_hash(cfg) == v12.REGISTERED_HASH
    with pytest.raises(v12.V12SettingsError, match="candidates"):
        v12.with_cap_rearth(cfg, 31)


# ---------------------------------------------------------------- 5. v1.1 paths unchanged
def test_select_sstar_v11_defaults_still_use_gpd():
    S = _sample(seed=3)
    direct = interim.gpd_threshold(S, 5000, 1.0 / 5000, 6.0)
    assert sg.select_sstar(S, 5000, 1.0 / 5000, 6.0) == direct
    assert sg.select_sstar(S, 5000, 1.0 / 5000, 6.0, cfg=_load(V11)) == direct


def test_v11_known_planet_report_does_not_require_v12(monkeypatch, tmp_path):
    def boom(*a, **k):
        raise AssertionError("v1.1 run_known_planet_recovery must not require the v1.2 settings")

    monkeypatch.setattr(v12, "require_v12_settings", boom)
    with pytest.raises(RuntimeError, match="Frozen eligibility"):
        kp.run_known_planet_recovery(10.0, out_dir=tmp_path, frozen_eligibility_file=None)


def test_finish_main_is_not_changed_to_require_v12():
    src = inspect.getsource(fin.main)
    assert "require_v12_settings" not in src
    assert "rescore_v12_cfg" not in src
    # the v1.2 entry point exists on the side and is gated
    assert callable(fin.rescore_v12_cfg) and callable(fin.rescore_rows_v12)
