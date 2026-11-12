from pathlib import Path
import numpy as np
import pandas as pd
import pytest

from exosearch.v11 import known_planets, grid, settings


def test_grid_halfwidth():
    """Verify plateau half-width calculation from grid spacing fallback."""
    P = 3.0
    baseline = 27.0
    hw = known_planets.calculate_grid_halfwidth(P, baseline=baseline)
    df = 0.04 / (27.0 ** 2)
    expected_hw = 0.5 * (3.0 ** 2) * df
    assert np.isclose(hw, expected_hw)


def test_alias_handling_exact_tolerance():
    """Spec A4: 2x and 0.5x aliases use the EXACT same half-width tolerance and are flagged."""
    P_cat = 4.0
    hw = known_planets.calculate_grid_halfwidth(P_cat)

    # 1. Exact match within half-width
    rec, alias = known_planets.check_recovery(P_cat, P_cat + 0.5 * hw, hw)
    assert rec and alias is None

    # 2. 2x alias within exact half-width
    rec_2x, alias_2x = known_planets.check_recovery(P_cat, 2.0 * P_cat + 0.5 * hw, hw)
    assert rec_2x and alias_2x == "2x"

    # 3. 0.5x alias within exact half-width
    rec_half, alias_half = known_planets.check_recovery(P_cat, 0.5 * P_cat + 0.5 * hw, hw)
    assert rec_half and alias_half == "0.5x"

    # 4. Alias beyond half-width fails even if within 1% of period
    rec_bad_alias, _ = known_planets.check_recovery(P_cat, 2.0 * P_cat + 2.0 * hw, hw)
    assert not rec_bad_alias


def test_s_star_required():
    """Spec B3: S* is required as argument; error if not supplied or non-positive."""
    df_tgt = pd.DataFrame([{"tic": 1, "period": 3.0, "dur_h": 2.0, "depth_ppm": 1000.0, "aperture_contamination": 0.9}])
    df_search = pd.DataFrame([{"tic": 1, "period": 3.0, "S": 15.0}])

    with pytest.raises(ValueError, match="s_star must be provided"):
        known_planets.evaluate_known_planets(df_tgt, df_search, s_star=None)

    with pytest.raises(ValueError, match="s_star must be provided"):
        known_planets.evaluate_known_planets(df_tgt, df_search, s_star=-1.0)


def test_b4_contamination_model_and_flag():
    """Spec B4: Aperture contamination primary; no aperture value -> ineligible with reason;
    flag > 20% difference vs TIC contratio."""
    s_star = 10.0
    targets = [
        # Star 1: aperture present, contratio matches well (0.95 vs 1/(1+0.05) ~ 0.952)
        {"tic": 1, "name": "P1", "period": 3.0, "t0": 1.0, "dur_h": 2.4, "depth_ppm": 2000.0, "aperture_contamination": 0.95, "contratio": 0.05, "n_tr": 3, "cdpp_d": 1e-4},
        # Star 2: aperture present, but differs by > 20% vs contratio (0.50 vs 1/(1+0.01) ~ 0.99) -> flag True
        {"tic": 2, "name": "P2", "period": 3.0, "t0": 1.0, "dur_h": 2.4, "depth_ppm": 2000.0, "aperture_contamination": 0.50, "contratio": 0.01, "n_tr": 3, "cdpp_d": 1e-4},
        # Star 3: no aperture contamination -> ineligible with reason
        {"tic": 3, "name": "P3", "period": 3.0, "t0": 1.0, "dur_h": 2.4, "depth_ppm": 2000.0, "aperture_contamination": None, "contratio": 0.05, "n_tr": 3, "cdpp_d": 1e-4},
    ]
    df_tgt = pd.DataFrame(targets)
    df_search = pd.DataFrame([
        {"tic": 1, "period": 3.0, "S": 15.0},
        {"tic": 2, "period": 3.0, "S": 15.0},
        {"tic": 3, "period": 3.0, "S": 15.0},
    ])

    out_df, summary = known_planets.evaluate_known_planets(df_tgt, df_search, s_star=s_star)

    r1 = out_df[out_df.tic == 1].iloc[0]
    assert r1["eligible"] == True
    assert not r1["contamination_flag"]

    r2 = out_df[out_df.tic == 2].iloc[0]
    assert r2["eligible"] == True
    assert r2["contamination_flag"] == True

    r3 = out_df[out_df.tic == 3].iloc[0]
    assert r3["eligible"] == False
    assert r3["ineligible_reason"] == "missing_aperture_contamination"


def test_b5_no_excused_misses_missing_lc_or_runtime_error_is_failure():
    """Spec B5: No excused misses: any eligible miss, including missing LC or runtime error,
    is a FAILURE row with diagnosed cause."""
    s_star = 10.0
    targets = [
        # Star 1: Eligible (expected S >= 10, n_tr >= 2), but search missed it
        {"tic": 1, "name": "P1", "period": 2.0, "t0": 0.5, "dur_h": 2.0, "depth_ppm": 2000.0, "aperture_contamination": 1.0, "n_tr": 4, "cdpp_d": 1e-4},
        # Star 2: Runtime error during pipeline preparation
        {"tic": 2, "name": "P2", "period": 3.0, "t0": 0.5, "dur_h": 2.0, "depth_ppm": 2000.0, "aperture_contamination": 1.0, "n_tr": 4, "cdpp_d": 1e-4, "error": "FITS corrupted"},
    ]
    df_tgt = pd.DataFrame(targets)
    df_search = pd.DataFrame([
        {"tic": 1, "period": 0.5, "S": 5.0},
        {"tic": 2, "period": 0.5, "S": 5.0},
    ])

    out_df, summary = known_planets.evaluate_known_planets(df_tgt, df_search, s_star=s_star)
    assert not summary["pass"]
    assert summary["eligible_missed"] == 2

    r1 = out_df[out_df.tic == 1].iloc[0]
    assert r1["miss_cause"] == "outcompeted"

    r2 = out_df[out_df.tic == 2].iloc[0]
    assert "runtime_error: FITS corrupted" in r2["miss_cause"]


def test_b5_report_refuses_without_frozen_eligibility_file(tmp_path):
    """Spec B5: Report refuses without frozen eligibility file."""
    with pytest.raises(RuntimeError, match="Frozen eligibility file required"):
        known_planets.run_known_planet_recovery(s_star=10.0, out_dir=tmp_path, frozen_eligibility_file=None)


def test_select_known_planet_targets_dedupe_and_pscomppars_preference(monkeypatch):
    """Spec A2, A3: PS rows preferred, dedupe within registered tolerance, no round(P, 4)."""
    from exosearch import catalogs

    # Mock catalogs
    mock_toi = pd.DataFrame([
        {"tic": 101, "toi": "101.01", "period": 3.14159, "disp": "CP", "dur_h": 2.5, "depth_ppm": 1200.0, "epoch_bjd": 2458000.1},
        {"tic": 102, "toi": "102.01", "period": 5.0, "disp": "PC", "dur_h": 3.0, "depth_ppm": 800.0, "epoch_bjd": 2458000.5},
    ])
    mock_ps = pd.DataFrame([
        # Planet on TIC 101 with period slightly different (within tolerance)
        {"tic": 101, "pl_name": "Planet-101 b", "pl_orbper": 3.14160, "pl_trandep": 0.15, "pl_trandur": 2.6, "pl_tranmid": 2458000.12},
    ])

    monkeypatch.setattr(catalogs, "load_catalogs", lambda: {"toi": mock_toi, "ps": mock_ps, "ctoi": pd.DataFrame()})
    # Avoid disk catalog loading
    monkeypatch.setattr(known_planets.config, "CATALOG_DIR", Path("/nonexistent/catalog/dir"))
    targets_df = known_planets.select_known_planet_targets(sector=70, searched_tics={101, 102})
    assert len(targets_df) == 2
    row_101 = targets_df[targets_df.tic == 101].iloc[0]
    assert row_101["type"] == "ps"
    assert row_101["name"] == "Planet-101 b"
    # Real depth and duration from pscomppars
    assert np.isclose(row_101["depth_ppm"], 1500.0)
    assert np.isclose(row_101["dur_h"], 2.6)

    row_102 = targets_df[targets_df.tic == 102].iloc[0]
    assert row_102["type"] == "toi"
    assert row_102["name"] == "TOI 102.01"


def test_recovery_requires_s_found_ge_s_star():
    """Spec B3: Peak matching period with S_found < S* must NOT count as recovered."""
    P_cat = 4.0
    hw = 0.01
    s_star = 12.0

    # Period matches, but S_found = 10.0 < s_star = 12.0
    rec, alias = known_planets.check_recovery(P_cat, P_cat, hw, s_found=10.0, s_star=s_star)
    assert not rec

    # Period matches and S_found = 13.0 >= s_star = 12.0
    rec_ok, alias_ok = known_planets.check_recovery(P_cat, P_cat, hw, s_found=13.0, s_star=s_star)
    assert rec_ok


def test_freeze_eligibility_end_to_end(tmp_path, monkeypatch):
    """Spec B5: freeze-eligibility writes dated CSV and sidecar with sha256, commit, and config_hash."""
    from exosearch import catalogs
    mock_toi = pd.DataFrame([
        {"tic": 501, "toi": "501.01", "period": 3.0, "disp": "CP", "dur_h": 2.0, "depth_ppm": 2000.0, "epoch_bjd": 2458000.1},
    ])
    monkeypatch.setattr(catalogs, "load_catalogs", lambda: {"toi": mock_toi, "ps": pd.DataFrame(columns=["tic"]), "ctoi": pd.DataFrame()})
    # Patch select_known_planet_targets to return a target
    mock_targets = pd.DataFrame([
        {"tic": 501, "name": "TOI 501.01", "period": 3.0, "t0": 1000.1, "dur_h": 2.0, "depth_ppm": 2000.0, "type": "toi"}
    ])
    monkeypatch.setattr(known_planets, "select_known_planet_targets", lambda **k: mock_targets)
    
    out_csv = tmp_path / "frozen_eligibility_s70.csv"
    eval_df, sidecar = known_planets.freeze_eligibility(s_star=10.0, out_file=out_csv, sector=70)
    assert out_csv.exists()
    sidecar_p = Path(str(out_csv) + ".meta.json")
    assert sidecar_p.exists()
    assert "timestamp_utc" in sidecar
    assert "sha256" in sidecar
    assert sidecar["sha256"] != ""
    assert "commit" in sidecar
    assert "config_hash" in sidecar
