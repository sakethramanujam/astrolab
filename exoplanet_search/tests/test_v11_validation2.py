import numpy as np
import pandas as pd
import pytest

from exosearch import inject, holdout
from exosearch.v11 import validate_pw2, registration_checks, settings


def test_seed_disjointness():
    """Verify new seeds: base seed 20261112, runs 12-15 are strictly disjoint
    from validate.py SEED 20261110 / run 0, cfg_v1 seed 70106 runs 0-11,
    and host-selection RNG seed inject.seed_for(..., sector, -1, 0)."""
    new_base = validate_pw2.BASE_SEED
    old_base = validate_pw2.EXCLUDED_BASE_SEED
    cfg_v1_seed = 70106
    assert new_base != old_base
    assert new_base == 20261112
    assert old_base == 20261110

    sample_tics = [1001, 1002, 1003, 20261110, 20261112, 70106]

    # Host selection RNG seed check
    host_seed_new = inject.seed_for(new_base, 70, -1, 0)
    host_seed_old_base = inject.seed_for(old_base, 70, -1, 0)
    host_seed_v1 = inject.seed_for(cfg_v1_seed, 70, -1, 0)
    assert host_seed_new != host_seed_old_base
    assert host_seed_new != host_seed_v1

    # Excluded task seeds
    new_seeds = set()
    for r in validate_pw2.INJECTION_RUNS:
        assert r in (12, 13, 14, 15)
        assert r not in validate_pw2.EXCLUDED_RUNS
        for tic in sample_tics:
            s = validate_pw2.seed_for_pw2(new_base, 70, r, tic)
            new_seeds.add(s)

    for r in range(12):
        for tic in sample_tics:
            assert inject.seed_for(cfg_v1_seed, 70, r, tic) not in new_seeds
    for tic in sample_tics:
        assert validate_pw2.seed_for_pw2(old_base, 70, 0, tic) not in new_seeds


def test_blocker1_negative_control_fails_on_noop_filter():
    """Spec C4 & A1: No-op filter FAILS part one; negative control on unfiltered run
    fails if unfiltered run does not fail (test_invalid = True)."""
    cfg = settings.cfg()
    s_star = 10.0

    # 23 Gaia stars with a pulsation frequency at 15.0 c/d
    # If filter is no-op, period found is e.g. 1 / 15.0 = 0.0667 d or residual S/N is 10.0 > 4.0
    gaia_results_noop = [
        {"tic": 1000 + i, "group": "delta_scuti", "period": 1.0 / 15.0, "pulsation_frequencies": [15.0], "max_residual_snr": 12.0}
        for i in range(23)
    ]
    unfiltered_noop = [
        {"tic": 1000 + i, "period": 1.0 / 15.0, "pulsation_frequencies": [15.0], "max_residual_snr": 12.0}
        for i in range(23)
    ]
    quiet_inj = [
        {"tic": 2000 + k, "expected_S": 15.0, "recovered": True, "prewhiten": True} for k in range(20)
    ] + [
        {"tic": 2000 + k, "expected_S": 15.0, "recovered": True, "prewhiten": False} for k in range(20)
    ]
    gaia_inj = [{"tic": 1000 + k, "expected_S": 15.0, "recovered": True, "group": "delta_scuti"} for k in range(23)]

    res = validate_pw2.evaluate_pw2_report(
        gaia_results_noop, gaia_inj, quiet_inj, cfg, unfiltered_gaia_results=unfiltered_noop, s_star=s_star
    )
    # Part one must fail
    assert not res["criterion_1_real_gaia_search"]["pass"]
    assert res["criterion_1_real_gaia_search"]["total_clean"] == 0
    assert not res["pass"]

    # Now test negative control trigger: if unfiltered run passes (clean_rate >= 0.90), test_invalid must be True
    unfiltered_clean = [
        {"tic": 1000 + i, "period": 3.5, "pulsation_frequencies": [15.0], "max_residual_snr": 1.0}
        for i in range(23)
    ]
    gaia_results_clean = [
        {"tic": 1000 + i, "group": "delta_scuti", "period": 3.5, "pulsation_frequencies": [15.0], "max_residual_snr": 1.0}
        for i in range(23)
    ]
    res_inv = validate_pw2.evaluate_pw2_report(
        gaia_results_clean, gaia_inj, quiet_inj, cfg, unfiltered_gaia_results=unfiltered_clean, s_star=s_star
    )
    assert res_inv["test_invalid"] == True
    assert not res_inv["pass"]


def test_pulsation_k_and_combination_harmonics_caught():
    """Spec C3: BLS top peak within one period-grid step of k/f_i (k=1..20)
    or of k/|f_i +/- f_j| (k=1..20) must be caught."""
    f1 = 12.0
    f2 = 15.0
    frozen = [f1, f2]
    baseline = 27.0
    df = 0.04 / (baseline ** 2)

    # 1. k = 3 / f1
    target_p1 = 3.0 / f1
    caught1, detail1 = validate_pw2.check_peak_pulsation_contamination(target_p1, frozen, baseline=baseline)
    assert caught1
    assert "k=3" in detail1

    # 2. k = 5 / (f1 + f2) = 5 / 27.0
    target_p2 = 5.0 / (f1 + f2)
    caught2, detail2 = validate_pw2.check_peak_pulsation_contamination(target_p2, frozen, baseline=baseline)
    assert caught2
    assert "combo" in detail2

    # 3. Clean period not near any harmonic or combo
    clean_p = 3.71234
    caught_clean, _ = validate_pw2.check_peak_pulsation_contamination(clean_p, frozen, baseline=baseline)
    assert not caught_clean


def test_b1_denominator_23_holdout_exception():
    """Spec B1: Denominator is 23 (holdout stars dropped as dated exception).
    Any listed star not evaluated counts as NOT clean."""
    cfg = settings.cfg()
    s_star = 10.0

    # 20 clean stars out of 23 -> 20/23 = 86.9% < 90% -> fails
    gaia_results = [
        {"tic": 1000 + i, "group": "delta_scuti", "period": 3.0, "pulsation_frequencies": [15.0], "max_residual_snr": 1.0}
        for i in range(20)
    ]
    # Add 3 missing
    for i in range(3):
        gaia_results.append({
            "tic": 1020 + i,
            "period": 0.0,
            "pulsation_frequencies": [],
            "unclean_reason": "missing_search_task",
        })
    # Add 2 holdout stars (must be ignored in denominator)
    gaia_results.append({"tic": 17228263, "period": 0.0, "pulsation_frequencies": [], "unclean_reason": "holdout star"})
    gaia_results.append({"tic": 279057510, "period": 0.0, "pulsation_frequencies": [], "unclean_reason": "holdout star"})

    quiet_inj = [
        {"tic": 2000 + k, "expected_S": 15.0, "recovered": True, "prewhiten": True} for k in range(20)
    ] + [
        {"tic": 2000 + k, "expected_S": 15.0, "recovered": True, "prewhiten": False} for k in range(20)
    ]
    gaia_inj = [{"tic": 1000 + k, "expected_S": 15.0, "recovered": True, "group": "delta_scuti"} for k in range(20)]

    res = validate_pw2.evaluate_pw2_report(gaia_results, gaia_inj, quiet_inj, cfg, s_star=s_star)
    c1 = res["criterion_1_real_gaia_search"]
    assert c1["total_stars"] == 23  # Exactly 23 stars!
    assert c1["total_clean"] == 20
    assert np.isclose(c1["total_clean_rate"], 20.0 / 23.0)
    assert not c1["pass"]


def test_b2_one_sided_star_level_bins_and_group_exclusion():
    """Spec B2: Part two is one-sided; variance unit is star; bins merged until n >= 5 stars;
    groups with < 5 stars are descriptive only."""
    cfg = settings.cfg()
    s_star = 10.0

    # Gaia injections: 4 gamma_dor stars (< 5 stars -> descriptive only), 10 delta_scuti stars
    gaia_inj = []
    # 4 gamma_dor stars
    for k in range(4):
        gaia_inj.append({"tic": 100 + k, "group": "gamma_dor", "expected_S": 12.0, "recovered": False, "tag": "inj", "run": 12, "seed": 1})
    # 10 delta_scuti stars (all recovered)
    for k in range(10):
        gaia_inj.append({"tic": 200 + k, "group": "delta_scuti", "expected_S": 12.0, "recovered": True, "tag": "inj", "run": 12, "seed": 1})

    # Quiet injections: 10 quiet stars
    quiet_inj = [
        {"tic": 300 + k, "expected_S": 12.0, "recovered": True, "prewhiten": True, "tag": "inj", "run": 12, "seed": 1}
        for k in range(10)
    ] + [
        {"tic": 300 + k, "expected_S": 12.0, "recovered": True, "prewhiten": False, "tag": "inj", "run": 12, "seed": 1}
        for k in range(10)
    ]

    gaia_results = [
        {"tic": 1000 + i, "group": "delta_scuti", "period": 3.0, "pulsation_frequencies": [15.0], "max_residual_snr": 1.0}
        for i in range(23)
    ]

    res = validate_pw2.evaluate_pw2_report(gaia_results, gaia_inj, quiet_inj, cfg, s_star=s_star)
    c2 = res["criterion_2_gaia_vs_quiet_injections"]
    # gamma_dor has 4 stars -> descriptive_only is True
    assert c2["groups"]["gamma_dor"]["descriptive_only"] == True
    # delta_scuti has 10 stars -> passes
    assert c2["groups"]["delta_scuti"]["descriptive_only"] == False
    assert c2["groups"]["delta_scuti"]["overall_pass"] == True


def test_a7_row_deduplication():
    """Spec A7: Dedupe rows per (TIC, tag, run, seed)."""
    cfg = settings.cfg()
    s_star = 10.0

    # Duplicate rows for the same injection
    gaia_inj_dup = [
        {"tic": 101, "group": "delta_scuti", "expected_S": 15.0, "recovered": True, "tag": "inj", "run": 12, "seed": 1},
        {"tic": 101, "group": "delta_scuti", "expected_S": 15.0, "recovered": True, "tag": "inj", "run": 12, "seed": 1},
    ]
    quiet_inj = [
        {"tic": 201, "expected_S": 15.0, "recovered": True, "prewhiten": True, "tag": "inj", "run": 12, "seed": 1},
        {"tic": 201, "expected_S": 15.0, "recovered": True, "prewhiten": False, "tag": "inj", "run": 12, "seed": 1},
    ]
    gaia_results = [
        {"tic": 1000 + i, "group": "delta_scuti", "period": 3.0, "pulsation_frequencies": [15.0], "max_residual_snr": 1.0}
        for i in range(23)
    ]
    # Duplicate real row
    gaia_results.append({"tic": 1000, "group": "delta_scuti", "period": 3.0, "pulsation_frequencies": [15.0], "max_residual_snr": 1.0})

    res = validate_pw2.evaluate_pw2_report(gaia_results, gaia_inj_dup, quiet_inj, cfg, s_star=s_star)
    assert res["criterion_1_real_gaia_search"]["total_stars"] == 23


def test_registration_guard():
    """Population guard: asserts no excluded TIC appears in paired null set."""
    df_ex = pd.DataFrame({"tic": [101, 102]})
    df_flip = pd.DataFrame({"tic": [101, 102, 201, 202, 203], "status": "ok"})
    df_scr = pd.DataFrame({"tic": [201, 202, 203, 301], "status": "ok"})
    res = registration_checks.check_registration_population(df_flip, df_scr, df_ex)
    assert res["status"] == "ok"

    df_scr_leaked = pd.DataFrame({"tic": [101, 201, 202], "status": "ok"})
    with pytest.raises(registration_checks.PopulationGuardError, match="Population guard failed"):
        registration_checks.check_registration_population(df_flip, df_scr_leaked, df_ex)
