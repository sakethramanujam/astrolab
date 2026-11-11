"""Independent analysis cross-check: sample products, and recovery of known gamma/GPD parameters."""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy import stats

from analysis.crosscheck import (
    expected_snr,
    fit_gamma_cdf,
    fit_gumbel_maxima,
    fit_tail,
    json_safe,
    mann_whitney_auc,
    reliability_from_counts,
    load_threshold,
    rescore_injections,
    run_crosscheck,
    wilson_interval,
)

ROOT = Path(__file__).resolve().parents[1]
SAMPLE_INJ = ROOT / "analysis" / "sample_injections.csv"
SAMPLE_NULL = ROOT / "analysis" / "sample_null.csv"
THRESH = ROOT / "configs" / "interim_threshold_v1.json"
# Pinned copy of the per-camera S* schema the sample-run counts were first checked against.
# The live interim file later switched to a pooled per_camera object and a different sector S*.
THRESH_PER_CAMERA = ROOT / "tests" / "fixtures" / "interim_threshold_per_camera.json"


def test_wilson_interval_known_value():
    lo, hi = wilson_interval(81, 100)
    assert lo == pytest.approx(0.7222115462093562)
    assert hi == pytest.approx(0.8748524849023126)
    assert wilson_interval(0, 0) == (None, None)


def test_expected_snr_matches_schema_formula():
    # snr_exp = depth_inj / cdpp_dur * sqrt(n_tr)
    snr = expected_snr(depth=0.001, n_tr=4.0, cdpp=0.0002)
    assert snr == pytest.approx(0.001 * math.sqrt(4.0) / 0.0002)
    assert expected_snr(depth=0.001, n_tr=0.0, cdpp=0.0002) is None
    assert expected_snr(depth=0.001, n_tr=4.0, cdpp=0.0) is None


def test_rescore_uses_supplied_s_star_not_pipeline_detected_flag():
    df = pd.DataFrame({
        "tic": [1, 2, 3],
        "sector": [70, 70, 70],
        "camera": [1, 1, 2],
        "status": ["ok", "ok", "ok"],
        "matched": [True, True, True],
        # row 0: above s_min so the pipeline says detected, but below the supplied S*
        # row 1: pipeline says not detected (S < s_min) but S clears the supplied S*
        # row 2: clears S* and fails a rejecting check, so detected but not passed
        "S": [8.0, 9.0, 20.0],
        "s_min": [7.0, 10.0, 7.0],
        "detected": [True, False, True],
        "failed_checks": ["", "", "odd_even"],
        "passed_all": [True, False, False],
    })
    out = rescore_injections(df, s_star_of=lambda sector, camera: 8.5)
    assert out["detected_star"].tolist() == [False, True, True]
    assert out["passed_star"].tolist() == [False, True, False]


def test_strict_json_turns_nan_into_null(tmp_path):
    from analysis.crosscheck import dump_json

    path = tmp_path / "x.json"
    dump_json(path, {"a": np.float64("nan"), "b": np.float64("inf"), "c": np.int64(3), "d": True})
    text = path.read_text()
    assert "NaN" not in text and "Infinity" not in text
    loaded = json.loads(text)
    assert loaded == {"a": None, "b": None, "c": 3, "d": True}
    assert json_safe(float("nan")) is None


def test_gamma_cdf_recovers_known_parameters_within_ci():
    rng = np.random.default_rng(0)
    shape, scale = 3.0, 4.0
    n_stars, reps = 600, 3
    snr = rng.uniform(0.5, 45.0, size=n_stars * reps)
    tics = np.repeat(np.arange(n_stars), reps)
    detected = rng.random(snr.size) < stats.gamma.cdf(snr, shape, scale=scale)
    fit = fit_gamma_cdf(snr, detected, tics, n_boot=120, seed=0)
    assert fit["shape"] is not None and abs(fit["shape"] - shape) < 0.8
    assert fit["scale"] is not None and abs(fit["scale"] - scale) < 1.2
    assert fit["shape_ci95"][0] <= shape <= fit["shape_ci95"][1]
    assert fit["scale_ci95"][0] <= scale <= fit["scale_ci95"][1]
    assert fit["shape_ci95"][1] - fit["shape_ci95"][0] < 3.0
    assert fit["seed"] == 0


def test_gpd_tail_recovers_known_parameters_within_ci():
    rng = np.random.default_rng(21)
    xi, sigma = 0.25, 1.5
    n = 4000
    S = stats.genpareto.rvs(xi, loc=0.0, scale=sigma, size=n, random_state=rng)
    true_s = float(stats.genpareto.isf(1.0 / n, xi, loc=0.0, scale=sigma))
    fit = fit_tail(S, np.arange(n), n_stars=n, n_boot=160, seed=21)
    assert fit["rate_in_tail"] is True
    assert abs(fit["xi"] - xi) < 0.25
    assert fit["xi_ci95"][0] <= xi <= fit["xi_ci95"][1]
    assert fit["s_star_ci95"][0] <= true_s <= fit["s_star_ci95"][1]
    assert fit["s_star_ci95"][1] - fit["s_star_ci95"][0] < 80
    assert fit["seed"] == 21
    assert fit["percentile"] == 95.0


def test_gpd_flags_extrapolation_past_the_largest_null():
    rng = np.random.default_rng(0)
    S = rng.exponential(scale=1.0, size=400) + 2.0
    fit = fit_tail(S, np.arange(400), n_stars=10_000_000, n_boot=20, seed=0)
    assert fit["s_star"] is not None
    assert fit["s_star"] > fit["max_observed"]
    assert fit["extrapolated"] is True


def test_gumbel_n_eff_on_maxima_of_exponential_trials():
    rng = np.random.default_rng(4)
    n_eff_true = 100
    maxima = (-np.log(rng.random((2500, n_eff_true)))).max(axis=1)
    fit = fit_gumbel_maxima(maxima, n_stars=2500, n_boot=40, seed=4)
    assert fit["n_eff"] == pytest.approx(n_eff_true, rel=0.15)
    assert fit["n_eff_ci95"][0] <= n_eff_true <= fit["n_eff_ci95"][1]
    assert "exponential" in fit["assumption"]


def test_reliability_thompson_formula_and_wilson_propagation():
    # E = 80/100 = 0.8, N_FP=10, N_PC=40
    # R = 1 - (10/40) * (1-0.8)/0.8 = 0.9375
    out = reliability_from_counts(n_fp=10, n_pc=40, n_rejected=80, n_sim=100)
    assert out["E"] == pytest.approx(0.8)
    assert out["R"] == pytest.approx(0.9375)
    assert out["E_ci95"][0] < 0.8 < out["E_ci95"][1]
    # (1-E)/E decreases in E, so R increases in E: lower E -> lower R
    assert out["R_ci95"][0] < out["R"] < out["R_ci95"][1]
    missing = reliability_from_counts(n_fp=None, n_pc=None, n_rejected=80, n_sim=100)
    assert missing["E"] == pytest.approx(0.8)
    assert missing["R"] is None


def test_auc_perfect_separation_and_ties():
    assert mann_whitney_auc(np.array([2.0, 3.0]), np.array([0.0, 1.0])) == pytest.approx(1.0)
    assert mann_whitney_auc(np.array([1.0, 1.0]), np.array([1.0, 1.0])) == pytest.approx(0.5)


def test_load_threshold_accepts_per_camera_records_and_pooled_metadata(tmp_path):
    table, population, raw, digest = load_threshold(THRESH_PER_CAMERA)
    assert table[(70, None)] == pytest.approx(14.241890281585157)
    assert table[(70, 1)] == pytest.approx(9.479378697517577)
    assert table[(70, 2)] == pytest.approx(14.91022300165031)
    assert (70, 3) not in table  # null s_star falls back to the sector value at resolve time
    assert table[(70, 4)] == pytest.approx(11.29349536853363)
    assert population["70"]["kept"] == 10244
    assert digest

    # Current interim schema: per_camera is pooled metadata, not camera records.
    pooled = {
        "population": {"70": {"kept": 10244}},
        "sectors": {
            "70": {
                "s_star": 11.067632388390026,
                "u": 7.0,
                "per_camera": {
                    "pooled": True,
                    "note": "cameras pooled",
                    "n_null": {"1": 1875, "2": 1887},
                },
            }
        },
    }
    path = tmp_path / "pooled.json"
    path.write_text(json.dumps(pooled))
    table, population, raw, digest = load_threshold(path)
    assert table == {(70, None): pytest.approx(11.067632388390026)}
    assert population["70"]["kept"] == 10244
    # boolean / note / count entries must not be treated as camera records
    assert (70, 1) not in table


def test_sample_run_loads_the_live_interim_config(tmp_path):
    result = run_crosscheck(
        injections=SAMPLE_INJ,
        nulls=[SAMPLE_NULL],
        out_dir=tmp_path / "live",
        threshold_config=THRESH,
        holdout_file=ROOT / "holdout_tics.txt",
        n_boot=5,
        seed=1,
        make_figures=False,
    )
    assert result["small_sample"] is True
    assert result["completeness"]["n_used"] == 31
    # Live sector 70 S* is used for every camera when per_camera is pooled.
    sources = {row["source"] for row in result["threshold"]["used"]}
    assert "sector" in sources
    assert result["completeness"]["n_detected_pipeline_flag"] == 14
    assert result["completeness"]["n_detected"] >= 1


def test_sample_run_excludes_holdout_and_writes_strict_outputs(tmp_path):
    out = tmp_path / "out"
    result = run_crosscheck(
        injections=SAMPLE_INJ,
        nulls=[SAMPLE_NULL],
        out_dir=out,
        threshold_config=THRESH_PER_CAMERA,
        holdout_file=ROOT / "holdout_tics.txt",
        n_boot=30,
        seed=20261110,
    )
    assert result["seed"] == 20261110
    assert result["pipeline_version"] == "1.0.0"
    assert result["config_hash"] == "80ffc7f29943"
    assert result["git_hash"] == "50906a55b0-dirty"
    assert result["small_sample"] is True
    assert result["holdout"]["mode"] == "tuning"
    # sample_injections.csv: 40 rows, 9 holdout; sample_null.csv: 24 rows, 3 holdout
    assert result["completeness"]["n_used"] == 31
    assert result["completeness"]["n_excluded_holdout"] == 9
    assert result["null_distribution"]["n_used"] == 21
    assert result["null_distribution"]["n_excluded_holdout"] == 3
    # Pipeline `detected` is S >= s_min. Camera-specific S* drops two of those.
    assert result["completeness"]["n_detected"] == 12
    assert result["completeness"]["n_detected_pipeline_flag"] == 14
    assert result["null_distribution"]["comparison"]["invert_vs_scramble"] is None
    assert result["reliability"]["overall"]["R"] is None
    assert result["reliability"]["overall"]["n_sim"] == 1
    assert result["reliability"]["overall"]["E"] == pytest.approx(1.0)
    for name in ("crosscheck.json", "completeness.json", "null_distribution.json",
                 "tail.json", "roc.json", "reliability.json"):
        text = (out / name).read_text()
        assert "NaN" not in text and "Infinity" not in text
        body = json.loads(text)
        for key in ("seed", "pipeline_version", "config_hash", "git_hash"):
            assert body[key] == result[key]
    figures = [
        "completeness_gamma.png", "completeness_map.png", "null_distributions.png",
        "tail_stability.png", "roc.png", "reliability_effectiveness.png",
    ]
    for name in figures:
        blob = (out / "figures" / name).read_bytes()
        assert blob[:8] == b"\x89PNG\r\n\x1a\n"
        assert len(blob) > 1000


def test_holdout_only_flag_uses_the_held_out_rows(tmp_path):
    result = run_crosscheck(
        injections=SAMPLE_INJ,
        nulls=[SAMPLE_NULL],
        out_dir=tmp_path / "hold",
        threshold_config=THRESH,
        holdout_file=ROOT / "holdout_tics.txt",
        holdout_only=True,
        n_boot=5,
        seed=3,
    )
    assert result["holdout"]["mode"] == "holdout"
    assert result["completeness"]["n_used"] == 9
    assert result["null_distribution"]["n_used"] == 3
    assert result["small_sample"] is True


def test_holdout_file_excludes_a_tic_even_if_the_column_says_false(tmp_path):
    inj = pd.read_csv(SAMPLE_INJ).head(4).copy()
    inj["holdout"] = False
    marked = int(inj.iloc[0]["tic"])
    hold = tmp_path / "holdout_tics.txt"
    hold.write_text(f"# comment\n{marked}\n")
    inj_path = tmp_path / "inj.csv"
    null = pd.read_csv(SAMPLE_NULL).head(6).copy()
    null["holdout"] = False
    inj.to_csv(inj_path, index=False)
    null_path = tmp_path / "null.csv"
    null.to_csv(null_path, index=False)
    result = run_crosscheck(
        injections=inj_path, nulls=[null_path], out_dir=tmp_path / "o",
        threshold_config=THRESH, holdout_file=hold, n_boot=3, seed=1,
    )
    assert result["completeness"]["n_excluded_holdout"] == 1
    assert result["completeness"]["n_used"] == 3


def test_products_root_parquet_and_candidates(tmp_path):
    root = tmp_path / "products" / "v1.0.0"
    inj_dir = root / "injections" / "injection_v1"
    null_dir = root / "nulltests" / "scramble_r0"
    inj_dir.mkdir(parents=True)
    null_dir.mkdir(parents=True)
    inj = pd.read_csv(SAMPLE_INJ)
    null = pd.read_csv(SAMPLE_NULL)
    null["mode"] = "scramble"
    inj.to_parquet(inj_dir / "part-0.parquet", index=False)
    null.to_parquet(null_dir / "null_scramble_r0_maxS.parquet", index=False)
    # one observed PC and one observed FP in a wide bin; R should be filled in
    cand = pd.DataFrame({
        "sector": [70, 70],
        "period": [5.0, 5.0],
        "rp_earth": [2.0, 2.0],
        "disposition": ["PC", "FP"],
        "tic": [1, 2],
        "holdout": [False, False],
    })
    cand_path = tmp_path / "candidates.csv"
    cand.to_csv(cand_path, index=False)
    out = tmp_path / "out"
    result = run_crosscheck(
        products_root=root, candidates=cand_path, out_dir=out,
        threshold_config=THRESH, holdout_file=ROOT / "holdout_tics.txt",
        n_boot=5, seed=1,
    )
    assert result["inputs"]["null_modes"] == ["scramble"]
    assert result["reliability"]["candidates_path"]
    # E is defined only where a simulated null clears S*; R stays null in bins with no such null.
    assert "overall" in result["reliability"]


def test_quality_table_drops_stars_above_the_registered_beta_cut(tmp_path):
    inj = pd.read_csv(SAMPLE_INJ).head(6).copy()
    inj["holdout"] = False
    null = pd.read_csv(SAMPLE_NULL).head(4).copy()
    null["holdout"] = False
    q = pd.DataFrame({
        "tic": list(inj.tic.astype(int).head(2)) + list(null.tic.astype(int).head(1)),
        "sector": [70, 70, 70],
        "beta_max": [1.0, 4.5, 3.1],
    })
    q_path = tmp_path / "quality.csv"
    q.to_csv(q_path, index=False)
    inj_path = tmp_path / "inj.csv"
    null_path = tmp_path / "null.csv"
    inj.to_csv(inj_path, index=False)
    null.to_csv(null_path, index=False)
    result = run_crosscheck(
        injections=inj_path, nulls=[null_path], quality=q_path, out_dir=tmp_path / "o",
        threshold_config=THRESH, holdout_file=tmp_path / "missing_holdout.txt",
        n_boot=3, seed=1,
    )
    assert result["quality_cut"]["applied"] is True
    assert result["quality_cut"]["threshold"] == 3.0
    # two injection hosts are in the table; beta_max 4.5 is dropped; the other four have no
    # measurement and are excluded (same rule as the interim layer).
    assert result["completeness"]["n_used"] == 1
    assert result["null_distribution"]["n_dropped_quality"] >= 1


def test_cli_on_sample_csvs(tmp_path):
    import subprocess
    import sys
    out = tmp_path / "cli"
    proc = subprocess.run(
        [sys.executable, "-m", "analysis",
         "--injections", str(SAMPLE_INJ),
         "--null", str(SAMPLE_NULL),
         "--out", str(out),
         "--n-boot", "4",
         "--seed", "5"],
        cwd=str(ROOT), check=False, capture_output=True, text=True,
    )
    assert proc.returncode == 0, proc.stderr
    body = json.loads((out / "crosscheck.json").read_text())
    assert body["seed"] == 5
    assert body["small_sample"] is True
    assert "too small" in body["small_sample_note"].lower() or "not a measurement" in body["small_sample_note"].lower()
