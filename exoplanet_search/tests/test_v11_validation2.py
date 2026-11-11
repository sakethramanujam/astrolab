import numpy as np
import pandas as pd
import pytest

from exosearch.v11 import validate_pw2, registration_checks, settings


def test_seed_disjointness():
    """Verify new seeds: base seed 20261112, runs 12-15 are strictly disjoint
    from validate.py SEED 20261110 / run 0."""
    new_base = validate_pw2.BASE_SEED
    old_base = validate_pw2.EXCLUDED_BASE_SEED
    assert new_base != old_base
    assert new_base == 20261112
    assert old_base == 20261110

    sample_tics = [1001, 1002, 1003, 20261110, 20261112]
    old_seeds = {validate_pw2.seed_for_pw2(old_base, 70, 0, tic) for tic in sample_tics}
    new_seeds = set()
    for r in validate_pw2.INJECTION_RUNS:
        assert r in (12, 13, 14, 15)
        assert r not in validate_pw2.EXCLUDED_RUNS
        for tic in sample_tics:
            s = validate_pw2.seed_for_pw2(new_base, 70, r, tic)
            new_seeds.add(s)

    # Disjoint sets
    assert len(old_seeds.intersection(new_seeds)) == 0


def test_group_split_at_5_cpd():
    """Classify star by dominant frequency before prewhitening:
    > 5 c/d is delta Scuti, <= 5 c/d is gamma Dor."""
    t = np.linspace(3000, 3027, 2000)
    # 1. Delta Scuti dominant frequency: e.g. 12 c/d
    f_dsct = 1.0 + 0.005 * np.sin(2 * np.pi * 12.0 * t)
    grp_dsct, f_dom_dsct = validate_pw2.classify_pulsator_dominant_frequency(t, f_dsct)
    assert grp_dsct == "delta_scuti"
    assert np.isclose(f_dom_dsct, 12.0, atol=0.2)

    # 2. Gamma Dor dominant frequency: e.g. 2.5 c/d
    f_gdor = 1.0 + 0.005 * np.sin(2 * np.pi * 2.5 * t)
    grp_gdor, f_dom_gdor = validate_pw2.classify_pulsator_dominant_frequency(t, f_gdor)
    assert grp_gdor == "gamma_dor"
    assert np.isclose(f_dom_gdor, 2.5, atol=0.2)

    # 3. Exactly at boundary (<= 5.0 is gamma Dor)
    f_boundary = 1.0 + 0.005 * np.sin(2 * np.pi * 5.0 * t)
    grp_b, f_dom_b = validate_pw2.classify_pulsator_dominant_frequency(t, f_boundary)
    assert grp_b == "gamma_dor"


def test_criteria_arithmetic_all_pass():
    """Check arithmetic for all three pass criteria."""
    # (i) 23 of 25 stars clean (>90%)
    gaia_results = []
    for i in range(25):
        clean = (i < 23)
        p = 3.0 if clean else 0.1  # P=3.0 -> f=0.333; P=0.1 -> f=10.0
        puls = [10.0]
        gaia_results.append({
            "tic": 1000 + i,
            "group": "delta_scuti" if i < 15 else "gamma_dor",
            "period": p,
            "pulsation_frequencies": puls
        })

    # (ii) Injections on Gaia stars vs quiet stars
    gaia_injections = []
    quiet_injections = []
    # Generate matched expected S bins
    for s_val in [12.0, 20.0, 30.0]:
        for k in range(20):
            gaia_injections.append({"tic": 1000 + k, "expected_S": s_val, "recovered": True, "group": "delta_scuti"})
            quiet_injections.append({"tic": 2000 + k, "expected_S": s_val, "recovered": True, "prewhiten": True})
            quiet_injections.append({"tic": 2000 + k, "expected_S": s_val, "recovered": True, "prewhiten": False})

    cfg = settings.cfg()
    res = validate_pw2.evaluate_pw2_report(gaia_results, gaia_injections, quiet_injections, cfg)
    assert res["pass"]
    assert res["criterion_1_real_gaia_search"]["pass"]
    assert res["criterion_1_real_gaia_search"]["total_clean_rate"] == 23 / 25
    assert res["criterion_2_gaia_vs_quiet_injections"]["pass"]
    assert res["criterion_3_quiet_prewhiten_drop"]["pass"]
    assert res["criterion_3_quiet_prewhiten_drop"]["drop_pp"] <= 0.02


def test_criteria_arithmetic_edge_cases_failing():
    """Edge cases: failing criterion 1, 2, or 3 blocks overall pass."""
    cfg = settings.cfg()

    # Case 1: Fail criterion 1 (< 90% clean, e.g. 22/25 = 88%)
    gaia_results = []
    for i in range(25):
        clean = (i < 22)
        gaia_results.append({
            "tic": 1000 + i,
            "group": "delta_scuti",
            "period": 3.0 if clean else 0.1,
            "pulsation_frequencies": [10.0]
        })
    res1 = validate_pw2.evaluate_pw2_report(gaia_results, [], [], cfg)
    assert not res1["criterion_1_real_gaia_search"]["pass"]
    assert not res1["pass"]

    # Case 2: Fail criterion 2 (Gaia recovery drops by > 5 pp vs quiet at same S)
    gaia_inj_bad = [{"tic": 1000 + k, "expected_S": 15.0, "recovered": (k < 10), "group": "delta_scuti"} for k in range(20)] # 50%
    quiet_inj_good = [{"tic": 2000 + k, "expected_S": 15.0, "recovered": True, "prewhiten": True} for k in range(20)] # 100%
    res2 = validate_pw2.evaluate_pw2_report([], gaia_inj_bad, quiet_inj_good, cfg)
    assert not res2["criterion_2_gaia_vs_quiet_injections"]["pass"]
    assert not res2["pass"]

    # Case 3: Fail criterion 3 (quiet recovery with pw drops by > 2 pp vs without)
    quiet_inj_drop = (
        [{"tic": 2000 + k, "expected_S": 15.0, "recovered": (k < 18), "prewhiten": True} for k in range(20)] + # 90%
        [{"tic": 2000 + k, "expected_S": 15.0, "recovered": True, "prewhiten": False} for k in range(20)]       # 100% (drop = 10 pp > 2 pp)
    )
    res3 = validate_pw2.evaluate_pw2_report([], [], quiet_inj_drop, cfg)
    assert not res3["criterion_3_quiet_prewhiten_drop"]["pass"]
    assert not res3["pass"]


def test_registration_guard():
    """Population guard: asserts no excluded TIC appears in paired null set."""
    # Excluded TICs: 101, 102
    df_ex = pd.DataFrame({"tic": [101, 102]})

    # 1. Valid case: excluded TICs only in flipped (or not in both), not in paired
    df_flip = pd.DataFrame({"tic": [101, 102, 201, 202, 203], "status": "ok"})
    df_scr = pd.DataFrame({"tic": [201, 202, 203, 301], "status": "ok"})
    res = registration_checks.check_registration_population(df_flip, df_scr, df_ex)
    assert res["status"] == "ok"
    assert res["paired_count"] == 3
    assert res["excluded_in_paired"] == 0

    # 2. Invalid case: an excluded TIC leaks into scrambled -> appears in paired set
    df_scr_leaked = pd.DataFrame({"tic": [101, 201, 202], "status": "ok"})
    with pytest.raises(AssertionError, match="Population guard failed"):
        registration_checks.check_registration_population(df_flip, df_scr_leaked, df_ex)
