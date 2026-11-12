"""End-to-end tests of run_pw2_report and build_pw2_manifests on synthetic light curves and fixture rows (r6 B1)."""
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from exosearch import batch, targets
from exosearch.v11 import run, settings, synth, validate_pw2

from lc_fixtures import write_spoc_fits

F_PULS = 12.345        # c/d, a delta Scuti mode well above 5 c/d


def _pulsator_lc(path, tic, amp=4e-3, seed=0):
    t = synth.s70_times()
    rng = np.random.default_rng(seed)
    f = 1.0 + amp * np.sin(2 * np.pi * F_PULS * t + 0.3) + rng.normal(0, 3e-4, t.size)
    return write_spoc_fits(path, t, f, e=np.full(t.size, 3e-4), tic=tic)


def _setup(tmp_path, monkeypatch, tics, lc_tics, with_unfiltered=True):
    gaia_csv = tmp_path / "dsct.csv"
    pd.DataFrame(dict(tic=tics)).to_csv(gaia_csv, index=False)
    freq = tmp_path / "freqs.csv"
    pd.DataFrame(dict(tic=tics, f_i=[F_PULS] * len(tics), amp=1e-3, snr=50.0, group="delta_scuti")).to_csv(freq, index=False)
    Path(str(freq) + ".meta.json").write_text(json.dumps(dict(sha256=hashlib.sha256(freq.read_bytes()).hexdigest())))
    files = {tic: _pulsator_lc(tmp_path / f"lc_{tic}.fits", tic) for tic in lc_tics}
    monkeypatch.setattr(targets, "fetch_manifest",
                        lambda sector, refresh=False: pd.DataFrame(dict(tic=tics, fname=[f"lc_{t}.fits" for t in tics])))
    monkeypatch.setattr(batch, "lc_path", lambda sector, fname, download=True: files.get(int(fname[3:-5])))
    rows = []
    for tic in tics:
        rows.append(dict(tic=tic, task="gaia_real_pw", mode="real", status="ok", period=1.7, P_grid=1.7))
        if with_unfiltered:
            rows.append(dict(tic=tic, task="gaia_real_nopw", mode="real", status="ok", period=1.7, P_grid=1.7))
    d = tmp_path / "rows"
    (d / "parts").mkdir(parents=True)
    pd.DataFrame(rows).to_parquet(d / "parts" / "part-00000.parquet")
    seen = {}
    real_eval = validate_pw2.evaluate_pw2_report

    def spy(**kw):
        seen.update(kw)
        return real_eval(**kw)
    monkeypatch.setattr(validate_pw2, "evaluate_pw2_report", spy)
    return d, gaia_csv, freq, seen


def test_run_pw2_report_measures_residual_snr(tmp_path, monkeypatch):
    """measure_residual_ls_snr actually runs (the swallowed NameError left it at 0): the unfiltered arm keeps the
    pulsation (S/N > 4), prewhitening removes it; one registered star only -> negative control incomplete."""
    d, gaia_csv, freq, seen = _setup(tmp_path, monkeypatch, [11], [11])
    res = validate_pw2.run_pw2_report([d], tmp_path / "out", gaia_csv=gaia_csv, frozen_frequencies_file=freq, s_star=9.7)
    pw = {g["tic"]: g for g in seen["gaia_results"]}
    unf = {g["tic"]: g for g in seen["unfiltered_gaia_results"]}
    assert unf[11]["max_residual_snr"] > 4.0
    assert 0.0 < pw[11]["max_residual_snr"] < unf[11]["max_residual_snr"]
    assert res["test_invalid"] and res["test_invalid_reason"].startswith("negative_control_incomplete")
    assert not res["pass"]


def test_run_pw2_report_missing_light_curve_raises(tmp_path, monkeypatch):
    d, gaia_csv, freq, _ = _setup(tmp_path, monkeypatch, [11, 12], [11])
    with pytest.raises(RuntimeError, match="light curve missing"):
        validate_pw2.run_pw2_report([d], tmp_path / "out", gaia_csv=gaia_csv, frozen_frequencies_file=freq, s_star=9.7)


def test_run_pw2_report_requires_s_star(tmp_path, monkeypatch):
    d, gaia_csv, freq, _ = _setup(tmp_path, monkeypatch, [11], [11])
    with pytest.raises(ValueError, match="s_star"):
        validate_pw2.run_pw2_report([d], tmp_path / "out", gaia_csv=gaia_csv, frozen_frequencies_file=freq)


def test_run_pw2_report_missing_unfiltered_arm_invalid(tmp_path, monkeypatch):
    d, gaia_csv, freq, seen = _setup(tmp_path, monkeypatch, [11], [11], with_unfiltered=False)
    res = validate_pw2.run_pw2_report([d], tmp_path / "out", gaia_csv=gaia_csv, frozen_frequencies_file=freq, s_star=9.7)
    assert all(u.get("unclean_reason") == "missing_search_task" for u in seen["unfiltered_gaia_results"])
    assert res["test_invalid"] and not res["criterion_1_real_gaia_search"]["pass"]


def _gaia(n_clean, n=23):
    return [dict(tic=5000 + i, group="delta_scuti", period=1.7, pulsation_frequencies=[F_PULS],
                 max_residual_snr=1.0 if i < n_clean else 9.0) for i in range(n)]


def test_negative_control_states():
    cfg = settings.cfg()
    good = _gaia(23)
    # complete unfiltered arm that fails part one -> valid
    r = validate_pw2.evaluate_pw2_report(good, [], [], cfg, unfiltered_gaia_results=_gaia(0), s_star=9.7)
    assert not r["test_invalid"] and r["criterion_1_real_gaia_search"]["pass"]
    # unfiltered arm that passes part one -> invalid
    r = validate_pw2.evaluate_pw2_report(good, [], [], cfg, unfiltered_gaia_results=_gaia(23), s_star=9.7)
    assert r["test_invalid"] and r["test_invalid_reason"] == "negative_control_passed_part_one"
    # missing or incomplete unfiltered arm -> invalid
    r = validate_pw2.evaluate_pw2_report(good, [], [], cfg, unfiltered_gaia_results=None, s_star=9.7)
    assert r["test_invalid"] and r["test_invalid_reason"] == "negative_control_missing"
    r = validate_pw2.evaluate_pw2_report(good, [], [], cfg, unfiltered_gaia_results=_gaia(0, n=20), s_star=9.7)
    assert r["test_invalid"] and r["test_invalid_reason"].startswith("negative_control_incomplete")


def _manifest_env(tmp_path, monkeypatch, exclude_file):
    from exosearch import config, inject, store
    gaia_csv = tmp_path / "dsct.csv"
    pd.DataFrame(dict(tic=[11, 12])).to_csv(gaia_csv, index=False)
    monkeypatch.setattr(store, "connect", lambda *a, **k: None)
    monkeypatch.setattr(inject, "select_hosts", lambda s, c, con: pd.DataFrame(dict(tic=[21, 22, 23, 11])))
    monkeypatch.setattr(config, "CATALOG_DIR", tmp_path / "catalogs")
    jobs = lambda tics: [dict(tic=t, sector=70, host="h2", rad_tic=1.0, mass_tic=1.0) for t in sorted(tics)]
    monkeypatch.setattr(run, "build", lambda sectors, weights, out, production=False, tics_filter=None,
                        tasks_override=None, write=True: jobs(tics_filter))
    written = {}
    monkeypatch.setattr(run, "write_manifests", lambda keep, out, cfg, cfg_inj, weights, name="m": written.update(keep=keep))
    cfg = settings.cfg()
    cfg.setdefault("validation_pw2", {})["quiet_exclude_pulsators_file"] = exclude_file
    monkeypatch.setattr(settings, "cfg", lambda path=None: json.loads(json.dumps(cfg)))
    return gaia_csv, written


def test_build_pw2_manifests_arms_and_exclusion(tmp_path, monkeypatch):
    ex = tmp_path / "pulsators.csv"
    pd.DataFrame(dict(tic=[22])).to_csv(ex, index=False)
    gaia_csv, written = _manifest_env(tmp_path, monkeypatch, str(ex))
    keep = validate_pw2.build_pw2_manifests([("h2", 1.0)], tmp_path / "m", gaia_csv=gaia_csv)
    by = {j["tic"]: {t["tag"] for t in j["tasks"]} for j in keep}
    from exosearch import holdout
    assert holdout.is_holdout(12) and not holdout.is_holdout(11)
    assert set(by) == {11, 21, 23}                            # 22 excluded as a pulsator, 12 is a holdout star
    assert {"gaia_real_pw", "gaia_real_nopw"} <= by[11] and not any(t.startswith("quiet") for t in by[11])
    assert by[21] == {f"quiet_{a}_r{r}" for a in ("pw", "nopw") for r in validate_pw2.INJECTION_RUNS}
    seeds = {t["seed"] for j in keep for t in j["tasks"] if "seed" in t}
    assert seeds == {validate_pw2.seed_for_pw2(validate_pw2.BASE_SEED, 70, r, tic)
                     for tic in (11, 21, 23) for r in validate_pw2.INJECTION_RUNS}
    assert written["keep"] is keep


def test_build_pw2_manifests_bad_exclusion_file_raises(tmp_path, monkeypatch):
    gaia_csv, _ = _manifest_env(tmp_path, monkeypatch, str(tmp_path / "nope.csv"))
    with pytest.raises(FileNotFoundError, match="quiet_exclude_pulsators_file"):
        validate_pw2.build_pw2_manifests([("h2", 1.0)], tmp_path / "m", gaia_csv=gaia_csv)
