import numpy as np
import pandas as pd
import pytest

from exosearch.v11 import known_planets, grid


def test_grid_halfwidth():
    """Verify plateau half-width calculation from grid spacing."""
    P = 3.0
    baseline = 27.0
    hw = known_planets.calculate_grid_halfwidth(P, baseline=baseline)
    # df = min(durations) / baseline^2 = 0.04 / (27^2)
    # hw = 0.5 * P^2 * df
    df = 0.04 / (27.0 ** 2)
    expected_hw = 0.5 * (3.0 ** 2) * df
    assert np.isclose(hw, expected_hw)


def test_contamination_factor_pdcsap():
    """Verify contamination factor is 1.0 for PDCSAP and crowdsap is tracked."""
    assert known_planets.CONTAMINATION_FACTOR_PDCSAP == 1.0

    targets_df = pd.DataFrame([{
        "tic": 123456,
        "name": "TOI 100.01",
        "period": 3.0,
        "duration_hours": 2.4,
        "depth_ppm": 1000.0,
        "crowdsap": 0.85,
        "n_tr": 3,
        "n_raw": 3,
        "cdpp_d": 1e-4,
    }])
    search_df = pd.DataFrame([{
        "tic": 123456,
        "period": 3.0,
        "S": 20.0
    }])
    s_star = 10.0
    out_df, summary = known_planets.evaluate_known_planets(targets_df, search_df, s_star=s_star)
    assert len(out_df) == 1
    row = out_df.iloc[0]
    assert row["contamination_factor"] == 1.0
    assert row["crowdsap"] == 0.85
    # depth_exp = depth_cat * 1.0 = 1e-3
    # expected S = 1e-3 * sqrt(3) / 1e-4 = 17.3205
    assert np.isclose(row["expected_S"], 1e-3 * np.sqrt(3) / 1e-4)


def test_alias_handling():
    """Verify exact match within halfwidth, 2x alias, 0.5x alias, and unrecovered."""
    P_cat = 4.0
    hw = known_planets.calculate_grid_halfwidth(P_cat)

    # 1. Exact match within half-width
    rec, alias = known_planets.check_recovery(P_cat, P_cat + 0.5 * hw, hw)
    assert rec and alias is None

    # 2. 2x alias
    rec_2x, alias_2x = known_planets.check_recovery(P_cat, 2.0 * P_cat, hw)
    assert rec_2x and alias_2x == "2x"

    # 3. 0.5x alias
    rec_half, alias_half = known_planets.check_recovery(P_cat, 0.5 * P_cat, hw)
    assert rec_half and alias_half == "0.5x"

    # 4. Unrecovered period (e.g. 1.33x)
    rec_bad, alias_bad = known_planets.check_recovery(P_cat, 1.33 * P_cat, hw)
    assert not rec_bad and alias_bad is None


def test_required_set_logic():
    """Target is required iff expected S >= S* and n_tr >= 2."""
    targets = [
        # TIC 1: expected S >= 10, n_tr = 3 -> Required
        {"tic": 1, "name": "T1", "period": 2.0, "dur_h": 2.0, "depth_ppm": 2000.0, "n_tr": 3, "n_raw": 3, "cdpp_d": 1e-4, "crowdsap": 0.95},
        # TIC 2: expected S < 10, n_tr = 3 -> Not required
        {"tic": 2, "name": "T2", "period": 2.0, "dur_h": 2.0, "depth_ppm": 500.0, "n_tr": 3, "n_raw": 3, "cdpp_d": 1e-4, "crowdsap": 0.95},
        # TIC 3: expected S >= 10, n_tr = 1 -> Not required (few transits)
        {"tic": 3, "name": "T3", "period": 15.0, "dur_h": 2.0, "depth_ppm": 2000.0, "n_tr": 1, "n_raw": 1, "cdpp_d": 1e-4, "crowdsap": 0.95},
    ]
    df_tgt = pd.DataFrame(targets)
    # Search recovers TIC 1 and 2, but misses TIC 3
    df_search = pd.DataFrame([
        {"tic": 1, "period": 2.0, "S": 15.0},
        {"tic": 2, "period": 2.0, "S": 8.0},
        {"tic": 3, "period": 0.8, "S": 6.0},
    ])
    s_star = 10.0
    out_df, summary = known_planets.evaluate_known_planets(df_tgt, df_search, s_star=s_star)
    
    assert out_df.loc[out_df.tic == 1, "required"].iloc[0] == True
    assert out_df.loc[out_df.tic == 2, "required"].iloc[0] == False
    assert out_df.loc[out_df.tic == 3, "required"].iloc[0] == False

    # Summary: only TIC 1 is required, and it was recovered -> Pass
    assert summary["pass"] == True
    assert summary["required_targets"] == 1
    assert summary["required_recovered"] == 1
    assert summary["required_missed"] == 0


def test_miss_causes_and_pass_fail():
    """Verify bucket classification logic and pass/fail gate: any required missed fails."""
    # Build 5 missed targets covering different causes
    targets = [
        # 1. few_transits: n_raw < 2
        {"tic": 10, "name": "T10", "period": 18.0, "duration_hours": 3.0, "depth_ppm": 2000.0, "n_raw": 1, "n_tr": 1, "cdpp_d": 1e-4},
        # 2. long_duration: duration > 0.25 d (6 hours); expected S < 10 so not required
        {"tic": 20, "name": "T20", "period": 5.0, "duration_hours": 7.0, "depth_ppm": 200.0, "n_raw": 3, "n_tr": 3, "cdpp_d": 1e-4},
        # 3. masked: n_raw >= 2, but n_tr surviving < 2
        {"tic": 30, "name": "T30", "period": 8.0, "duration_hours": 2.0, "depth_ppm": 2000.0, "n_raw": 2, "n_tr": 1, "cdpp_d": 1e-4},
        # 4. below_threshold: expected S < s_star
        {"tic": 40, "name": "T40", "period": 3.0, "duration_hours": 2.0, "depth_ppm": 400.0, "n_raw": 5, "n_tr": 5, "cdpp_d": 1e-4},
        # 5. outcompeted: required (S >= 10, n_tr >= 2), but search found completely different peak
        {"tic": 50, "name": "T50", "period": 3.0, "duration_hours": 2.0, "depth_ppm": 1500.0, "n_raw": 5, "n_tr": 5, "cdpp_d": 1e-4},
    ]
    df_tgt = pd.DataFrame(targets)
    # Search found wrong period for all
    df_search = pd.DataFrame([
        {"tic": 10, "period": 0.5, "S": 5.0},
        {"tic": 20, "period": 0.5, "S": 5.0},
        {"tic": 30, "period": 0.5, "S": 5.0},
        {"tic": 40, "period": 0.5, "S": 5.0},
        {"tic": 50, "period": 0.5, "S": 25.0},
    ])
    s_star = 10.0
    out_df, summary = known_planets.evaluate_known_planets(df_tgt, df_search, s_star=s_star)

    row10 = out_df[out_df.tic == 10].iloc[0]
    assert row10["miss_cause"] == "few_transits"

    row20 = out_df[out_df.tic == 20].iloc[0]
    assert row20["miss_cause"] == "long_duration"

    row30 = out_df[out_df.tic == 30].iloc[0]
    assert row30["miss_cause"] == "masked"

    row40 = out_df[out_df.tic == 40].iloc[0]
    assert row40["miss_cause"] == "below_threshold"

    row50 = out_df[out_df.tic == 50].iloc[0]
    assert row50["miss_cause"] == "outcompeted"
    assert row50["required"] == True
    assert row50["recovered"] == False

    # Because TIC 50 is required and missed, summary['pass'] must be False
    assert summary["pass"] == False
    assert summary["required_missed"] == 1
