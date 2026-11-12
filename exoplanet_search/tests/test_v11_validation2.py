import numpy as np
import pandas as pd
import pytest

from exosearch import inject
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

    # 1. Host selection RNG seed check
    host_seed_new = inject.seed_for(new_base, 70, -1, 0)
    host_seed_old_base = inject.seed_for(old_base, 70, -1, 0)
    host_seed_v1 = inject.seed_for(cfg_v1_seed, 70, -1, 0)
    assert host_seed_new != host_seed_old_base
    assert host_seed_new != host_seed_v1

    # 2. Excluded task seeds
    # Old validate.py run-0 seeds under base seed 20261110
    old_validate_seeds = {validate_pw2.seed_for_pw2(old_base, 70, 0, tic) for tic in sample_tics}
    # v1.0 injection runs 0..11 under seed 70106
    v1_injection_seeds = {
        inject.seed_for(cfg_v1_seed, 70, r, tic)
        for r in range(12)
        for tic in sample_tics
    }
    # validate.py run-0 seeds under seed 70106
    v1_run0_seeds = {inject.seed_for(cfg_v1_seed, 70, 0, tic) for tic in sample_tics}

    # New seeds under base seed 20261112, runs 12..15
    new_seeds = set()
    for r in validate_pw2.INJECTION_RUNS:
        assert r in (12, 13, 14, 15)
        assert r not in validate_pw2.EXCLUDED_RUNS
        for tic in sample_tics:
            s = validate_pw2.seed_for_pw2(new_base, 70, r, tic)
            new_seeds.add(s)

    # Disjointness checks
    assert len(old_validate_seeds.intersection(new_seeds)) == 0
    assert len(v1_injection_seeds.intersection(new_seeds)) == 0
    assert len(v1_run0_seeds.intersection(new_seeds)) == 0


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
            "f_dom": 12.0 if i < 15 else 3.0,
            "pulsation_frequencies": puls
        })

    # (ii) Injections on Gaia stars vs quiet stars
    gaia_injections = []
    quiet_injections = []
    # Generate matched expected S bins
    for s_val in [12.0, 20.0, 30.0]:
        for k in range(20):
            grp = "delta_scuti" if k < 12 else "gamma_dor"
            gaia_injections.append({"tic": 1000 + k, "expected_S": s_val, "recovered": True, "group": grp})
            quiet_injections.append({"tic": 2000 + k, "expected_S": s_val, "recovered": True, "prewhiten": True})
            quiet_injections.append({"tic": 2000 + k, "expected_S": s_val, "recovered": True, "prewhiten": False})

    cfg = settings.cfg()
    res = validate_pw2.evaluate_pw2_report(gaia_results, gaia_injections, quiet_injections, cfg)
    assert res["pass"]
    assert res["criterion_1_real_gaia_search"]["pass"]
    assert res["criterion_1_real_gaia_search"]["total_clean_rate"] == 23 / 25
    assert len(res["criterion_1_real_gaia_search"]["real_gaia_table"]) == 25
    assert "pulsation_frequencies" in res["criterion_1_real_gaia_search"]["real_gaia_table"][0]
    assert res["criterion_2_gaia_vs_quiet_injections"]["pass"]
    assert "delta_scuti" in res["criterion_2_gaia_vs_quiet_injections"]["groups"]
    assert "gamma_dor" in res["criterion_2_gaia_vs_quiet_injections"]["groups"]
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
    quiet_pw_dummy = [{"tic": 2000 + k, "expected_S": 15.0, "recovered": True, "prewhiten": True} for k in range(20)]
    quiet_nopw_dummy = [{"tic": 2000 + k, "expected_S": 15.0, "recovered": True, "prewhiten": False} for k in range(20)]
    gaia_inj_dummy = [{"tic": 1000 + k, "expected_S": 15.0, "recovered": True, "group": "delta_scuti"} for k in range(20)]
    res1 = validate_pw2.evaluate_pw2_report(gaia_results, gaia_inj_dummy, quiet_pw_dummy + quiet_nopw_dummy, cfg)
    assert not res1["criterion_1_real_gaia_search"]["pass"]
    assert not res1["pass"]

    # Case 2: Fail criterion 2 (Gaia recovery drops by > 5 pp vs quiet at same S)
    gaia_inj_bad = [{"tic": 1000 + k, "expected_S": 15.0, "recovered": (k < 10), "group": "delta_scuti"} for k in range(20)] # 50%
    quiet_inj_good = (
        [{"tic": 2000 + k, "expected_S": 15.0, "recovered": True, "prewhiten": True} for k in range(20)] + # 100%
        [{"tic": 2000 + k, "expected_S": 15.0, "recovered": True, "prewhiten": False} for k in range(20)]
    )
    # 25 clean Gaia stars to pass criterion 1
    gaia_results_good = [{"tic": 1000 + i, "group": "delta_scuti", "period": 3.0, "pulsation_frequencies": [10.0]} for i in range(25)]
    res2 = validate_pw2.evaluate_pw2_report(gaia_results_good, gaia_inj_bad, quiet_inj_good, cfg)
    assert not res2["criterion_2_gaia_vs_quiet_injections"]["pass"]
    assert not res2["pass"]

    # Case 3: Fail criterion 3 (quiet recovery with pw drops by > 2 pp vs without)
    quiet_inj_drop = (
        [{"tic": 2000 + k, "expected_S": 15.0, "recovered": (k < 18), "prewhiten": True} for k in range(20)] + # 90%
        [{"tic": 2000 + k, "expected_S": 15.0, "recovered": True, "prewhiten": False} for k in range(20)]       # 100% (drop = 10 pp > 2 pp)
    )
    gaia_inj_good = [{"tic": 1000 + k, "expected_S": 15.0, "recovered": True, "group": "delta_scuti"} for k in range(20)]
    res3 = validate_pw2.evaluate_pw2_report(gaia_results_good, gaia_inj_good, quiet_inj_drop, cfg)
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
    with pytest.raises(registration_checks.PopulationGuardError, match="Population guard failed"):
        registration_checks.check_registration_population(df_flip, df_scr_leaked, df_ex)


def test_empty_samples_fail_every_criterion():
    """F1: Every criterion must FAIL when its sample is empty or incomplete."""
    cfg = settings.cfg()

    # Completely empty
    res_empty = validate_pw2.evaluate_pw2_report([], [], [], cfg)
    assert not res_empty["pass"]
    assert not res_empty["criterion_1_real_gaia_search"]["pass"]
    assert not res_empty["criterion_2_gaia_vs_quiet_injections"]["pass"]
    assert not res_empty["criterion_3_quiet_prewhiten_drop"]["pass"]

    # Incomplete Gaia real stars (e.g. 24 of 25)
    gaia_incomplete = [{"tic": 1000 + i, "group": "delta_scuti", "period": 3.0, "pulsation_frequencies": [10.0]} for i in range(24)]
    res_incomplete = validate_pw2.evaluate_pw2_report(gaia_incomplete, [], [], cfg)
    assert not res_incomplete["criterion_1_real_gaia_search"]["pass"]
    assert res_incomplete["criterion_1_real_gaia_search"]["total_stars"] == 24


def test_missing_or_holdout_star_counts_against_criterion1_denominator_25():
    """F2: criterion 1 denominator must be the full registered list (25 stars).
    Any listed star not evaluated (holdout star, missing light curve, missing stellar parameters,
    task error, no dominant frequency) counts as NOT clean and is listed in the report with its reason."""
    cfg = settings.cfg()

    # Construct report run scenario where only 20 stars are evaluated cleanly,
    # 2 are holdout, 1 has missing light curve, 1 has missing stellar params, 1 has no dominant frequency
    gaia_eval = []
    for i in range(20):
        gaia_eval.append({
            "tic": 1000 + i,
            "group": "delta_scuti" if i < 10 else "gamma_dor",
            "period": 3.0,
            "pulsation_frequencies": [10.0],
            "f_dom": 12.0 if i < 10 else 2.0,
        })
    # Add the remaining 5 stars with reasons
    gaia_eval.append({"tic": 1020, "group": "unclassified", "period": 0.0, "f_dom": 0.0, "pulsation_frequencies": [], "unclean_reason": "holdout star"})
    gaia_eval.append({"tic": 1021, "group": "gamma_dor", "period": 0.0, "f_dom": 1.5, "pulsation_frequencies": [], "unclean_reason": "holdout star"})
    gaia_eval.append({"tic": 1022, "group": "unclassified", "period": 0.0, "f_dom": 0.0, "pulsation_frequencies": [], "unclean_reason": "missing light curve"})
    gaia_eval.append({"tic": 1023, "group": "delta_scuti", "period": 0.0, "f_dom": 8.0, "pulsation_frequencies": [], "unclean_reason": "missing stellar parameters"})
    gaia_eval.append({"tic": 1024, "group": "unclassified", "period": 0.0, "f_dom": 0.0, "pulsation_frequencies": [], "unclean_reason": "no dominant frequency"})

    quiet_inj = (
        [{"tic": 2000 + k, "expected_S": 15.0, "recovered": True, "prewhiten": True} for k in range(20)] +
        [{"tic": 2000 + k, "expected_S": 15.0, "recovered": True, "prewhiten": False} for k in range(20)]
    )
    gaia_inj = [{"tic": 1000 + k, "expected_S": 15.0, "recovered": True, "group": "delta_scuti"} for k in range(20)]

    res = validate_pw2.evaluate_pw2_report(gaia_eval, gaia_inj, quiet_inj, cfg)
    c1 = res["criterion_1_real_gaia_search"]
    assert c1["total_stars"] == 25
    assert c1["total_clean"] == 20
    assert c1["total_clean_rate"] == 20 / 25  # 80% < 90% -> fails
    assert not c1["pass"]
    assert "unclassified" in c1["groups"]
    assert c1["groups"]["unclassified"]["n_total"] == 3
    assert c1["groups"]["unclassified"]["n_clean"] == 0

    reasons = {r["tic"]: r["unclean_reason"] for r in c1["real_gaia_table"] if not r["clean_peak"]}
    assert reasons[1020] == "holdout star"
    assert reasons[1021] == "holdout star"
    assert reasons[1022] == "missing light curve"
    assert reasons[1023] == "missing stellar parameters"
    assert reasons[1024] == "no dominant frequency"


def test_gaia_star_without_stellar_params_manifest_tasks(tmp_path):
    """F3: Gaia stars without stellar params get a real task but no injection task.
    Holdout stars stay out of every arm."""
    from exosearch import holdout
    # Create synthetic test using build_pw2_manifests
    jobs = validate_pw2.build_pw2_manifests([("local", 1.0)], tmp_path)
    job_map = {j["tic"]: j for j in jobs}

    # TIC 17228263 is a holdout star -> must NOT be in jobs
    assert holdout.is_holdout(17228263)
    assert 17228263 not in job_map

    # TIC 279057510 is also a holdout star -> must NOT be in jobs
    assert holdout.is_holdout(279057510)
    assert 279057510 not in job_map

    # TIC 284704271 is a Gaia star without stellar params (not holdout)
    # -> must get a real task but no injection tasks
    assert not holdout.is_holdout(284704271)
    assert 284704271 in job_map
    j_284 = job_map[284704271]
    tags = [t["tag"] for t in j_284["tasks"]]
    assert tags == ["gaia_real_pw"]
    assert all("inj" not in t for t in tags)
