"""v1.2 steps 6–7: frozen settings file, content hash, refusing entry points, and output stamps."""
import copy
import importlib.util
import inspect
import json
import re
from pathlib import Path
from types import SimpleNamespace

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


# ---------------------------------------------------------------- 6. v1.2 output stamps (step 7)
def _loaded_hash():
    return v12.content_hash(v12.require_v12_settings())


def _toy_targets():
    rows = []
    for tic, name in ((1, "TOI 1.01"), (2, "TOI 2.01")):
        rows.append(dict(tic=tic, name=name, type="toi", period=3.0 + tic, period_err=0.001,
                         t0=3000.5, dur_h=2.0, depth_ppm=2000.0, depth_src="exofop_toi",
                         dur_src="exofop_toi", t0_src="exofop_toi"))
    return pd.DataFrame(rows)


def _freeze(tmp_path, cap_rearth=None, name="frozen.csv"):
    kw = {} if cap_rearth is None else {"cap_rearth": cap_rearth}
    return kp.freeze_eligibility_v12(10.0, tmp_path / name, targets_df=_toy_targets(),
                                     lc_lookup=lambda tic: None, contratio_map={}, commit="test",
                                     config_hash="test", **kw)


def _search_stub(tmp_path):
    p = tmp_path / "search.csv"
    pd.DataFrame([dict(tic=999, mode="real", period=1.5, S=2.0, status="ok")]).to_csv(p, index=False)
    return p


def test_require_v12_stamps_accepts_consistent_rows_and_refuses_bad_ones():
    h = _loaded_hash()
    good = pd.DataFrame([{v12.STAMP_HASH_KEY: h, v12.STAMP_CAP_KEY: 30.0},
                         {v12.STAMP_HASH_KEY: h, v12.STAMP_CAP_KEY: 30.0}])
    assert v12.require_v12_stamps(good) is good
    assert v12.require_v12_stamps(good, cap_rearth=30) is good
    # empty input is vacuously accepted
    assert v12.require_v12_stamps(pd.DataFrame()) is not None
    assert v12.require_v12_stamps([]) == []

    with pytest.raises(v12.V12SettingsError, match="lacks stamp"):
        v12.require_v12_stamps(pd.DataFrame([{v12.STAMP_HASH_KEY: h}]))
    with pytest.raises(v12.V12SettingsError, match="lacks stamp"):
        v12.require_v12_stamps({v12.STAMP_CAP_KEY: 30.0})

    flipped = h[:-1] + ("0" if h[-1] != "0" else "1")
    with pytest.raises(v12.V12SettingsError, match="hash"):
        v12.require_v12_stamps(pd.DataFrame([{v12.STAMP_HASH_KEY: flipped, v12.STAMP_CAP_KEY: 30.0}]))
    with pytest.raises(v12.V12SettingsError, match="mix settings hashes"):
        v12.require_v12_stamps(pd.DataFrame([{v12.STAMP_HASH_KEY: h, v12.STAMP_CAP_KEY: 30.0},
                                             {v12.STAMP_HASH_KEY: flipped, v12.STAMP_CAP_KEY: 30.0}]))
    with pytest.raises(v12.V12SettingsError, match="mix cap_rearth"):
        v12.require_v12_stamps(pd.DataFrame([{v12.STAMP_HASH_KEY: h, v12.STAMP_CAP_KEY: 22.0},
                                             {v12.STAMP_HASH_KEY: h, v12.STAMP_CAP_KEY: 30.0}]))
    with pytest.raises(v12.V12SettingsError, match="cap_rearth"):
        v12.require_v12_stamps(good, cap_rearth=22)


def test_stamp_v12_does_not_mutate_the_caller():
    h = _loaded_hash()
    row = {"tic": 1}
    out = v12.stamp_v12(row, h, 30)
    assert v12.STAMP_HASH_KEY not in row and out[v12.STAMP_HASH_KEY] == h
    df = pd.DataFrame([{"tic": 1}])
    stamped = v12.stamp_v12(df, h, 22.0)
    assert v12.STAMP_HASH_KEY not in df.columns
    assert stamped[v12.STAMP_CAP_KEY].tolist() == [22.0]


def test_calibrate_sstar_v12_result_is_stamped():
    S = _sample()
    h = _loaded_hash()
    r = v12.calibrate_sstar_v12(S, 600, 5.0 / 600, 6.0)
    assert r["v12_settings_hash"] == h
    assert r["v12_cap_rearth"] == 30.0
    r22 = v12.calibrate_sstar_v12(S, 600, 5.0 / 600, 6.0, cap_rearth=22)
    assert r22["v12_settings_hash"] == h and r22["v12_cap_rearth"] == 22.0
    assert r22["s_star"] == r["s_star"]


def test_calibrate_sstar_v12_refuses_bad_rows():
    S = _sample()
    h = _loaded_hash()
    good = pd.DataFrame([{v12.STAMP_HASH_KEY: h, v12.STAMP_CAP_KEY: 30.0}])
    assert v12.calibrate_sstar_v12(S, 600, 5.0 / 600, 6.0, rows=good)["s_star"] > 0
    with pytest.raises(v12.V12SettingsError, match="lacks stamp"):
        v12.calibrate_sstar_v12(S, 600, 5.0 / 600, 6.0, rows=pd.DataFrame([{"tic": 1}]))
    cap22 = pd.DataFrame([{v12.STAMP_HASH_KEY: h, v12.STAMP_CAP_KEY: 22.0}])
    with pytest.raises(v12.V12SettingsError, match="cap_rearth"):
        v12.calibrate_sstar_v12(S, 600, 5.0 / 600, 6.0, cap_rearth=30, rows=cap22)
    # default applied cap is 30; rows stamped 22 must be refused even without cap_rearth=
    with pytest.raises(v12.V12SettingsError, match="cap_rearth"):
        v12.calibrate_sstar_v12(S, 600, 5.0 / 600, 6.0, rows=cap22)


def test_lock_sstar_v12_stamps_and_refuses_bad_rows():
    h = _loaded_hash()
    rows = pd.DataFrame([{v12.STAMP_HASH_KEY: h, v12.STAMP_CAP_KEY: 30.0}])
    r = v12.lock_sstar_v12(rows, 13.71)
    assert r["s_star"] == 13.71
    assert r["v12_settings_hash"] == h and r["v12_cap_rearth"] == 30.0
    with pytest.raises(v12.V12SettingsError, match="lacks stamp"):
        v12.lock_sstar_v12(pd.DataFrame([{"tic": 1}]), 13.71)
    with pytest.raises(v12.V12SettingsError, match="cap_rearth"):
        v12.lock_sstar_v12(rows, 13.71, cap_rearth=22)
    cap22 = pd.DataFrame([{v12.STAMP_HASH_KEY: h, v12.STAMP_CAP_KEY: 22.0}])
    with pytest.raises(v12.V12SettingsError, match="cap_rearth"):
        v12.lock_sstar_v12(cap22, 13.71)


def test_rescore_rows_v12_stamps_output(monkeypatch):
    monkeypatch.setattr(fin, "rebuild_for_row", lambda row, jobs, cfg: (None, None, None))
    monkeypatch.setattr(fin, "row_verdict", lambda row, **kw: SimpleNamespace(
        passed=True, first_failing_rule="", all_failing=[], status="ok"))
    h = _loaded_hash()
    out = fin.rescore_rows_v12([dict(tic=1, mode="real")])
    assert out[v12.STAMP_HASH_KEY].tolist() == [h]
    assert out[v12.STAMP_CAP_KEY].tolist() == [30.0]
    out22 = fin.rescore_rows_v12([dict(tic=1, mode="real")], cap_rearth=22)
    assert out22[v12.STAMP_HASH_KEY].tolist() == [h]
    assert out22[v12.STAMP_CAP_KEY].tolist() == [22.0]
    assert fin.rescore_v12_cfg(cap_rearth=22)["vet"]["rp_cap"]["cap_rearth"] == 22.0


def test_freeze_eligibility_v12_stamps_frame_csv_and_sidecar(tmp_path):
    h = _loaded_hash()
    df, side = _freeze(tmp_path)
    assert df[v12.STAMP_HASH_KEY].tolist() == [h, h]
    assert df[v12.STAMP_CAP_KEY].tolist() == [30.0, 30.0]
    assert side[v12.STAMP_HASH_KEY] == h and side[v12.STAMP_CAP_KEY] == 30.0
    csv = tmp_path / "frozen.csv"
    assert side["sha256"] == kp._sha256(csv)
    on_disk = pd.read_csv(csv)
    assert on_disk[v12.STAMP_HASH_KEY].tolist() == [h, h]
    assert on_disk[v12.STAMP_CAP_KEY].tolist() == [30.0, 30.0]
    back = kp.read_frozen(csv)
    assert len(back) == 2 and set(back["error"]) == {"missing light curve"}


def test_freeze_eligibility_v12_cap_override_is_stamped(tmp_path):
    h = _loaded_hash()
    df, side = _freeze(tmp_path, cap_rearth=22)
    assert df[v12.STAMP_CAP_KEY].tolist() == [22.0, 22.0]
    assert df[v12.STAMP_HASH_KEY].tolist() == [h, h]
    assert side[v12.STAMP_CAP_KEY] == 22.0 and side[v12.STAMP_HASH_KEY] == h
    with pytest.raises(v12.V12SettingsError, match="cap_rearth"):
        _freeze(tmp_path, cap_rearth=23, name="bad.csv")


def test_evaluate_known_planets_v12_stamps_and_refuses_targets(tmp_path):
    h = _loaded_hash()
    frozen, _ = _freeze(tmp_path)
    out, summ = kp.evaluate_known_planets_v12(frozen, pd.DataFrame(), 10.0)
    assert out[v12.STAMP_HASH_KEY].tolist() == [h, h]
    assert out[v12.STAMP_CAP_KEY].tolist() == [30.0, 30.0]
    assert summ[v12.STAMP_HASH_KEY] == h and summ[v12.STAMP_CAP_KEY] == 30.0

    with pytest.raises(v12.V12SettingsError, match="lacks stamp"):
        kp.evaluate_known_planets_v12(_toy_targets(), pd.DataFrame(), 10.0)
    # default applied cap is 30; a cap-22 eligibility table is refused without cap_rearth=
    frozen22, _ = _freeze(tmp_path, cap_rearth=22, name="frozen22.csv")
    with pytest.raises(v12.V12SettingsError, match="cap_rearth"):
        kp.evaluate_known_planets_v12(frozen22, pd.DataFrame(), 10.0)
    # empty targets are vacuous (existing step-6 missing-settings tests keep their shape)
    assert kp.evaluate_known_planets_v12(pd.DataFrame(), pd.DataFrame(), 10.0)[0].empty


def test_run_known_planet_recovery_v12_stamps_and_refuses_v11_freeze(tmp_path):
    h = _loaded_hash()
    frozen, _ = _freeze(tmp_path)
    res = _search_stub(tmp_path)
    out, summ = kp.run_known_planet_recovery_v12(10.0, search_results_path=res, out_dir=tmp_path / "out",
                                                 frozen_eligibility_file=tmp_path / "frozen.csv")
    assert out[v12.STAMP_HASH_KEY].tolist() == [h, h]
    assert out[v12.STAMP_CAP_KEY].tolist() == [30.0, 30.0]
    disk = pd.read_csv(tmp_path / "out" / "known_planets_recovery.csv")
    assert disk[v12.STAMP_HASH_KEY].tolist() == [h, h]
    js = json.loads((tmp_path / "out" / "known_planets_summary.json").read_text())
    assert js[v12.STAMP_HASH_KEY] == h and js[v12.STAMP_CAP_KEY] == 30.0
    assert summ[v12.STAMP_HASH_KEY] == h and summ[v12.STAMP_CAP_KEY] == 30.0

    # a v1.1 (unstamped) frozen eligibility file is refused
    v11_file = tmp_path / "v11_frozen.csv"
    kp.freeze_eligibility(10.0, v11_file, targets_df=_toy_targets(), lc_lookup=lambda tic: None,
                          contratio_map={}, commit="test", config_hash="test")
    with pytest.raises(v12.V12SettingsError, match="lacks stamp"):
        kp.run_known_planet_recovery_v12(10.0, search_results_path=res, out_dir=tmp_path / "out11",
                                         frozen_eligibility_file=v11_file)


def test_v12_wrappers_accept_cap_rearth_but_have_no_cfg_argument(tmp_path):
    h = _loaded_hash()
    frozen, _ = _freeze(tmp_path, cap_rearth=22)
    res = _search_stub(tmp_path)
    out, summ = kp.run_known_planet_recovery_v12(10.0, search_results_path=res, out_dir=tmp_path / "out22",
                                                 frozen_eligibility_file=tmp_path / "frozen.csv",
                                                 cap_rearth=22)
    assert out[v12.STAMP_CAP_KEY].tolist() == [22.0, 22.0]
    assert out[v12.STAMP_HASH_KEY].tolist() == [h, h]
    assert summ[v12.STAMP_CAP_KEY] == 22.0 and summ[v12.STAMP_HASH_KEY] == h

    # mismatch between the stamped input cap and the requested aggregator cap
    frozen30, _ = _freeze(tmp_path, name="frozen30.csv")
    with pytest.raises(v12.V12SettingsError, match="cap_rearth"):
        kp.run_known_planet_recovery_v12(10.0, search_results_path=res, out_dir=tmp_path / "badcap",
                                         frozen_eligibility_file=tmp_path / "frozen30.csv", cap_rearth=22)
    # default applied cap is 30; a cap-22 eligibility file is refused without cap_rearth=
    with pytest.raises(v12.V12SettingsError, match="cap_rearth"):
        kp.run_known_planet_recovery_v12(10.0, search_results_path=res, out_dir=tmp_path / "badcap22",
                                         frozen_eligibility_file=tmp_path / "frozen.csv")

    # the supported override is cap_rearth= only: the wrappers expose no cfg= argument
    for fn in (kp._v12_run_cfg, kp.evaluate_known_planets_v12, kp.freeze_eligibility_v12,
               kp.run_known_planet_recovery_v12):
        assert "cfg" not in inspect.signature(fn).parameters


def test_v12_wrappers_accept_every_registered_cap_and_stamp_the_loaded_hash(tmp_path):
    h = _loaded_hash()
    for cap in (22, 25, 30):
        df, side = _freeze(tmp_path, cap_rearth=cap, name=f"frozen{cap}.csv")
        assert df[v12.STAMP_HASH_KEY].tolist() == [h, h]
        assert df[v12.STAMP_CAP_KEY].tolist() == [float(cap), float(cap)]
        assert side[v12.STAMP_HASH_KEY] == h and side[v12.STAMP_CAP_KEY] == float(cap)

    with pytest.raises(v12.V12SettingsError, match="cap_rearth"):
        _freeze(tmp_path, cap_rearth=23, name="bad.csv")

    # an otherwise-frozen cfg mutated in a non-frozen key cannot be passed in at all
    cfg = v12.require_v12_settings()
    cfg["innocuous_extra"] = 1
    with pytest.raises(TypeError):
        kp.evaluate_known_planets_v12(pd.DataFrame(), pd.DataFrame(), 10.0, cfg=cfg)
    with pytest.raises(TypeError):
        kp.freeze_eligibility_v12(10.0, tmp_path / "x.csv", targets_df=_toy_targets(),
                                  lc_lookup=lambda tic: None, contratio_map={}, cfg=cfg)
    with pytest.raises(TypeError):
        kp.run_known_planet_recovery_v12(10.0, search_results_path=_search_stub(tmp_path),
                                         out_dir=tmp_path / "x", cfg=cfg)


def test_v11_freeze_output_is_not_stamped(tmp_path):
    df, side = kp.freeze_eligibility(10.0, tmp_path / "v11.csv", targets_df=_toy_targets(),
                                     lc_lookup=lambda tic: None, contratio_map={}, commit="t", config_hash="t")
    assert v12.STAMP_HASH_KEY not in df.columns
    assert v12.STAMP_CAP_KEY not in df.columns
    assert v12.STAMP_HASH_KEY not in side
