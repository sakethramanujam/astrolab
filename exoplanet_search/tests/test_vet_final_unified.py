"""Unified final verdict (exosearch.vet_final): real, null and inject rows are scored by ONE function.

The same row dict must give the same verdict through the real, null and inject paths of
runs/v1.1_s70_x86/finish.py, the amended rules (a)-(d) must hold in every path, and for secondary_scan
rows with a rebuilt light-curve scan (or mocked rebuild_flat) final_verdict must match vet_amend.revet_one
plus common_period and checks_logged exactly.
"""
import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from exosearch import vet_final as vf
from exosearch.v11 import vet_amend as va

_FINISH = Path(__file__).resolve().parents[1] / "runs" / "v1.1_s70_x86" / "finish.py"
_spec = importlib.util.spec_from_file_location("s70_finish", _FINISH)
fin = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(fin)

_REPRO = Path(__file__).resolve().parents[1] / "runs" / "v1.1_s70_x86" / "repro_old_vs_new.py"
_repro_spec = importlib.util.spec_from_file_location("s70_repro", _REPRO)
repro = importlib.util.module_from_spec(_repro_spec)
_repro_spec.loader.exec_module(repro)

BASE = dict(tic=1, sector=70, status="ok", mode="real", task="real", S=20.0, P_grid=3.0, period=3.0,
            depth=0.01, t0=0.5, duration=0.1, flags="", deep_flags="", rp_rjup=1.0, rp_earth=11.0,
            secondary_depth=0.0, sec_scan_sigma=1.2, checks_logged=True, matched=True)


def row(**kw):
    return dict(BASE, **kw)


PATHS = (fin.real_path, fin.null_path, fin.inject_path)


def _synth_lc(one_orbit_only: bool, sec_dip: float = 0.004, seed: int = 1):
    rng = np.random.default_rng(seed)
    t = np.arange(0.0, 24.0, 2.0 / 1440.0)
    P, t0, dur = 3.0, 0.5, 0.1
    f = 1.0 + rng.normal(0, 3e-4, t.size)
    ph = ((t - t0) / P) % 1.0
    f[np.abs(((t - t0 + 0.5 * P) % P) - 0.5 * P) < dur / 2] -= 0.01
    sec = np.abs(ph - 0.6) < (dur / P) / 2
    if one_orbit_only:
        sec &= (t > t0 + 3 * P) & (t < t0 + 4 * P)
    f[sec] -= sec_dip
    return t, f, np.full(t.size, 3e-4), P, t0


def test_same_verdict_through_all_three_paths():
    scan_pass = dict(sec_scan_sigma=6.5, sec_scan_phase=0.6, sec_scan_depth=1e-5, duration_fine=0.1,
                     t0_fine=0.5, loo_min_sigma=9.0, loo_min_orbit=0)
    scan_fail = dict(sec_scan_sigma=6.5, sec_scan_phase=0.6, sec_scan_depth=0.004, duration_fine=0.1,
                     t0_fine=0.5, loo_min_sigma=9.0, loo_min_orbit=0)
    cases = [
        (row(), None),                                                                  # clean pass
        (row(mode="invert", task="invert"), None),                                      # null flavour
        (row(mode="inject", task="inject_run0", matched=True), None),                   # inject flavour
        (row(deep_flags="density_mismatch"), None),                                     # (c): warn only -> PASS
        (row(rp_rjup=25.0 / 11.209, rp_earth=25.0, deep_flags="rp_gt_2rjup"), None),    # ~25 R_Earth: passes
        (row(rp_rjup=2.8, rp_earth=31.0, deep_flags="rp_gt_2rjup"), None),              # > 30 R_Earth: rejects
        (row(deep_flags="v_shape"), None),                                              # unchanged check rejects
        (row(deep_flags="secondary_scan", sec_scan_sigma=6.5), scan_pass),              # (a) albedo exempts
        (row(deep_flags="secondary_scan", sec_scan_sigma=6.5, rp_rjup=2.5), scan_fail), # scan rejects
        (row(checks_logged=False), None),                                               # checks never logged -> reject
        (row(flags="odd_even"), None),                                                  # screen flag rejects
    ]
    for r, sc in cases:
        verdicts = [p(r, scan=sc) for p in PATHS]
        assert verdicts[0] == verdicts[1] == verdicts[2], r.get("deep_flags")
        assert verdicts[0].status == vf.STATUS_OK


def test_secondary_scan_matches_revet_one_with_mocked_rebuild(monkeypatch):
    """For synthetic secondary_scan rows, final_verdict (with supplied scan or via finish.py's rebuild_scan)
    must decide rules (a) and (d) and the sigma consistency check identically to vet_amend.revet_one."""
    lc_loo_exempt = _synth_lc(one_orbit_only=True, sec_dip=0.004)     # single-orbit dip: LOO < 5 -> (d) exempts
    lc_persistent = _synth_lc(one_orbit_only=False, sec_dip=0.0008)   # multi-orbit dip: LOO >= 5
    lc_tiny_sec = _synth_lc(one_orbit_only=False, sec_dip=0.00018)    # multi-orbit dip with small depth -> (a) exempts

    for lc_tuple in (lc_loo_exempt, lc_persistent, lc_tiny_sec):
        t, f, e, P, t0 = lc_tuple
        sc_true = va.scan_depth(t, f, e, P, t0)
        sig = float(sc_true["sec_scan_sigma"])
        assert sig > 5.0
        monkeypatch.setattr(va, "rebuild_flat", lambda job, r, _t=t, _f=f, _e=e: (_t, _f, _e))

        for rp_rjup in (1.0, 2.5, 2.85):
            for rho in (None, 1.0):
                for cp in (False, True):
                    for cl in (True, False):
                        r = row(P_grid=P, period=P, t0=t0, deep_flags="secondary_scan",
                                sec_scan_sigma=sig, rp_rjup=rp_rjup, rp_earth=rp_rjup * va.R_EARTH_PER_RJUP,
                                checks_logged=cl)
                        job = {"tic": 1, "rho": rho}
                        rev = va.revet_one(r, job, {})
                        expected_pass = bool(rev["passed_lc_checks_v2"]) and (not cp) and cl

                        # 1. Explicit scan supplied to final_verdict / final_verdict_row
                        sc = fin.rebuild_scan(r, job, {})
                        v_direct = vf.final_verdict(r, scan=sc, rho=rho, common_period=cp, checks_logged=cl)
                        d_direct = vf.final_verdict_row(r, scan=sc, rho=rho, common_period=cp, checks_logged=cl)
                        assert v_direct.status == vf.STATUS_OK
                        assert v_direct.passed is expected_pass
                        for k in ("passed_lc_checks_v2", "reject_v2", "warn_v2", "rp_earth_v2", "a_rs_v2",
                                  "a_rs_src_v2", "albedo_secondary_v2", "albedo_scan_v2",
                                  "sec_scan_depth_v2", "sec_scan_phase_v2", "loo_min_sigma_v2", "loo_min_orbit_v2"):
                            val_got, val_want = d_direct[k], rev[k]
                            if isinstance(val_want, float) and np.isnan(val_want):
                                assert np.isnan(val_got)
                            else:
                                assert val_got == val_want, k

                        # 2. Through finish.py's real_path and null_path (which rebuild automatically)
                        v_real = fin.real_path(dict(r, mode="real"), jobs={1: job}, cfg={}, fail_common_period=cp)
                        v_null = fin.null_path(dict(r, mode="invert"), jobs={1: job}, cfg={}, fail_common_period=cp)
                        assert v_real == v_null == v_direct

    # Sigma consistency check: mismatch between rebuilt scan and logged sec_scan_sigma must raise RuntimeError
    t, f, e, P, t0 = lc_persistent
    monkeypatch.setattr(va, "rebuild_flat", lambda job, r: (t, f, e))
    bad_row = row(P_grid=P, period=P, t0=t0, deep_flags="secondary_scan", sec_scan_sigma=99.9)
    with pytest.raises(RuntimeError, match="rebuilt secondary scan"):
        va.revet_one(bad_row, {"tic": 1}, {})
    with pytest.raises(RuntimeError, match="rebuilt secondary scan"):
        fin.real_path(bad_row, jobs={1: {"tic": 1}}, cfg={})
    with pytest.raises(RuntimeError, match="rebuilt secondary scan"):
        vf.final_verdict(bad_row, scan=dict(sec_scan_sigma=6.0))


def test_amended_rules_hold_in_every_path():
    for p in PATHS:
        out = p(row(deep_flags="density_mismatch"))
        assert out.passed and out.all_failing == []            # (c) density mismatch is a flag
        out = p(row(rp_rjup=25.0 / 11.209, rp_earth=25.0, deep_flags="rp_gt_2rjup"))
        assert out.passed and "rp_gt_30re" not in out.all_failing    # (b) ~25 R_Earth survives
        assert "rp_gt_2rjup" not in out.all_failing
        out = p(row(deep_flags="v_shape"))
        assert not out.passed and out.first_failing_rule == "v_shape"
        out = p(row(rp_rjup=2.8, rp_earth=31.39, deep_flags="rp_gt_2rjup"))
        assert not out.passed and out.first_failing_rule == "rp_gt_30re"


def test_terms_common_period_and_checks_logged():
    r = row()
    assert vf.final_verdict(r).passed and vf.final_verdict(r).status == vf.STATUS_OK
    v = vf.final_verdict(r, common_period=True)
    assert not v.passed and v.first_failing_rule == "common_period" and v.all_failing == ["common_period"]
    v = vf.final_verdict(row(flags="v_shape,common_period"), checks_logged=False)
    assert not v.passed and v.all_failing == ["v_shape", "common_period", "checks_not_logged"]
    assert vf.final_verdict(r, checks_logged=np.nan).all_failing == ["checks_not_logged"]
    # a common_period flag in the row's flags is the group flag term -> the row fails, but it is not
    # also counted inside the amended rejection set
    v = vf.final_verdict(row(flags="common_period"))
    assert not v.passed and v.all_failing == ["common_period"] and v.first_failing_rule == "common_period"
    v = vf.final_verdict(row(flags="common_period"), checks_logged=False)
    assert v.all_failing == ["common_period", "checks_not_logged"]


def test_needs_lc_bracket():
    r = row(mode="inject", task="inject_run0", deep_flags="secondary_scan", sec_scan_sigma=6.0)
    lo, hi = vf.final_verdict(r), vf.final_verdict(r, needs_lc="exempt")
    assert lo.status == vf.STATUS_NEEDS_LC and not lo.passed and lo.all_failing == ["secondary_scan"]
    assert hi.status == vf.STATUS_NEEDS_LC and hi.passed and hi.all_failing == []
    # only the undecided scan rule is waived; other rules still reject in the upper bound
    hi2 = vf.final_verdict(row(deep_flags="secondary_scan,v_shape", sec_scan_sigma=6.0), needs_lc="exempt")
    assert not hi2.passed and hi2.all_failing == ["v_shape"]
    # an explicit scan with matching sec_scan_sigma resolves the status
    scan = dict(sec_scan_sigma=6.0, sec_scan_phase=0.6, sec_scan_depth=0.0, duration_fine=0.1,
                t0_fine=0.5, loo_min_sigma=9.0, loo_min_orbit=0)
    ok = vf.final_verdict(row(deep_flags="secondary_scan", sec_scan_sigma=6.0), scan=scan)
    assert ok.status == vf.STATUS_OK and ok.passed
    # inject_path leaves secondary_scan rows bracketed as needs_lc when scan is not supplied
    assert fin.inject_path(r).status == vf.STATUS_NEEDS_LC and not fin.inject_path(r).passed
    assert fin.inject_path(r, needs_lc="exempt").passed


def test_repro_old_vs_new_identical_on_real_and_null(monkeypatch, tmp_path):
    t, f, e, P, t0 = _synth_lc(one_orbit_only=True, sec_dip=0.004)
    sig = float(va.scan_depth(t, f, e, P, t0)["sec_scan_sigma"])
    monkeypatch.setattr(va, "rebuild_flat", lambda job, r: (t, f, e))

    rows_df = pd.DataFrame([
        row(tic=1, mode="real", task="real"),
        row(tic=2, mode="real", task="real", deep_flags="secondary_scan", P_grid=P, period=P, t0=t0, sec_scan_sigma=sig),
        row(tic=3, mode="invert", task="invert", deep_flags="secondary_scan", P_grid=P, period=P, t0=t0, sec_scan_sigma=sig),
        row(tic=4, mode="scramble", task="scramble_r0", deep_flags="density_mismatch"),
        row(tic=5, mode="inject", task="inject_run0", deep_flags="secondary_scan", sec_scan_sigma=sig),
    ])
    rows_dir = tmp_path / "rows_known"
    rows_dir.mkdir()
    rows_df.to_parquet(rows_dir / "rows_0.parquet")

    summ = repro.compare_rows(rows_df, jobs={1: {"tic": 1, "rho": 1.0}, 2: {"tic": 2, "rho": 1.0}}, cfg={})
    assert summ["n_compared"] == 4
    assert summ["n_identical"] == 4
    assert summ["n_differing"] == 0
    assert summ["n_identical_rules"] == 4
    assert summ["n_differing_rules"] == 0
    assert summ["all_real_null_identical"] is True
    assert summ["needs_lc"] == dict(real=0, null=0, inject=1, invert=0, scramble=0, total=1)


# ---------------------------------------------------------------- v1.2 step 2: odd-even transit mask in unified vetter
def _cfg_mask(enabled: bool = True, **overrides):
    from exosearch.v11 import settings
    c = settings.cfg()
    m = dict(c["vet"]["odd_even_mask"], enabled=enabled, **overrides)
    c["vet"] = dict(c["vet"], odd_even_mask=m)
    return c


def test_masked_odd_even_excludes_corrupted_partial_or_gappy_transit():
    """Plant a deep partial/dump-adjacent transit on one parity that causes unmasked odd-even to fail,
    and verify that when the mask is on it is excluded, depths match the clean kept transits, and the
    unified verdict passes across real, null and inject paths."""
    from astropy.timeseries import BoxLeastSquares
    from exosearch.v11 import vet11

    dt = 0.005
    P, t0, D = 3.0, 1.0, 0.10
    depth_true = 0.002
    t_full = np.arange(0.0, 20.0, dt)
    e_full = np.full(t_full.size, 1e-4)
    f_full = np.ones_like(t_full)

    # Epochs 0..6 at t_mid = 1, 4, 7, 10, 13, 16, 19:
    # Even epochs: 0 (1.0), 2 (7.0), 4 (13.0), 6 (19.0)
    # Odd epochs:  1 (4.0), 3 (10.0), 5 (16.0)
    ep = np.round((t_full - t0) / P).astype(int)
    in_tr = (t_full >= (t0 + ep * P - 0.5 * D - 1e-9)) & (t_full < (t0 + ep * P + 0.5 * D - 1e-9))
    f_full[in_tr] -= depth_true

    # Plant a 10x deeper corrupted transit on odd epoch 1 (t_mid = 4.0) and make it partial (50% coverage < 0.75)
    f_full[in_tr & (ep == 1)] -= 9.0 * depth_true
    drop_ep1_half = (t_full >= 3.95 - 1e-9) & (t_full < 4.00 - 1e-9)  # removes 10 of 20 cadences -> coverage = 0.50
    t = t_full[~drop_ep1_half]
    f = f_full[~drop_ep1_half]
    e = e_full[~drop_ep1_half]

    # Unmasked BoxLeastSquares.compute_stats includes the partial corrupted transit on epoch 1 -> fails odd_even
    st = BoxLeastSquares(t, f, dy=e).compute_stats(P, D, t0)
    best = dict(depth_exact=depth_true, P_grid=P, t0=t0, duration=D, n_dump_epochs=0, n_epochs=7)
    scr_off = vet11.screens(st, best, {"radius": 1.0}, _cfg_mask(False), lc=(t, f, e))
    assert "odd_even" in scr_off["flags"]

    # Masked odd-even excludes epoch 1 (coverage 0.50 < 0.75) and keeps odd epochs 3 & 5 (n_odd=2, n_even=4)
    oe_res = vf.masked_odd_even(t, f, e, P, t0, D, cfg=_cfg_mask(True))
    assert 1 not in oe_res["mask"]["kept_odd_epochs"]
    assert list(oe_res["mask"]["kept_odd_epochs"]) == [3, 5]
    assert list(oe_res["mask"]["kept_even_epochs"]) == [0, 2, 4, 6]
    assert np.isclose(oe_res["depth_odd"], depth_true, atol=1e-12)
    assert np.isclose(oe_res["depth_even"], depth_true, atol=1e-12)
    assert oe_res["odd_even_sigma"] < 1e-6
    assert oe_res["flag"] is None and oe_res["passed"] is True

    # screens() with mask on also passes
    scr_on = vet11.screens(st, best, {"radius": 1.0}, _cfg_mask(True), lc=(t, f, e))
    assert "odd_even" not in scr_on["flags"] and "odd_even_too_few" not in scr_on["flags"]

    # Unified final_verdict: even when the stored row had unmasked flags="odd_even", mask-on recomputes from lc and passes
    r_stored = row(P_grid=P, period=P, t0=t0, duration=D, flags="odd_even")
    assert not vf.final_verdict(r_stored, cfg=_cfg_mask(False)).passed
    for p in PATHS:
        v_on = p(r_stored, cfg=_cfg_mask(True), lc=(t, f, e))
        assert v_on.passed and v_on.all_failing == [] and v_on.status == vf.STATUS_OK


def test_odd_even_too_few_flag_and_verdict():
    """If < 2 odd or < 2 even transits survive the mask, `odd_even_too_few` is set, `odd_even` is NOT set,
    and the row does not pass in any mode."""
    dt = 0.005
    P, t0, D = 3.0, 1.0, 0.10
    # Baseline [0.0, 11.0): epochs 0 (1.0), 1 (4.0), 2 (7.0), 3 (10.0) -> 2 even (0, 2) and 2 odd (1, 3)
    t = np.arange(0.0, 11.0, dt)
    f = np.ones_like(t)
    e = np.full(t.size, 1e-4)
    ep = np.round((t - t0) / P).astype(int)
    in_tr = (t >= (t0 + ep * P - 0.5 * D - 1e-9)) & (t < (t0 + ep * P + 0.5 * D - 1e-9))
    f[in_tr] -= 0.002

    # Put a momentum dump at t = 4.10 (within Y = 0.375 d of odd epoch 1 at t_mid = 4.0),
    # leaving only 1 odd transit (epoch 3) and 2 even transits (epochs 0, 2) -> n_odd = 1 < 2!
    dumps = np.array([4.10])
    oe_res = vf.masked_odd_even(t, f, e, P, t0, D, cfg=_cfg_mask(True), dumps=dumps)
    assert oe_res["n_odd"] == 1 and oe_res["n_even"] == 2
    assert oe_res["too_few"] is True
    assert oe_res["odd_even_fail"] is False
    assert oe_res["flag"] == vf.ODD_EVEN_TOO_FEW
    assert oe_res["passed"] is False

    r = row(P_grid=P, period=P, t0=t0, duration=D, flags="odd_even")
    for mode, p in zip(("real", "invert", "inject"), PATHS):
        rm = dict(r, mode=mode)
        v = p(rm, cfg=_cfg_mask(True), lc=(t, f, e), dumps=dumps)
        assert not v.passed
        assert v.first_failing_rule == "odd_even_too_few"
        assert v.all_failing == ["odd_even_too_few"]
        assert "odd_even" not in v.all_failing
        assert v.status == vf.STATUS_OK


def test_mask_off_reproduces_head_verdicts_without_lc():
    """With the mask off (the default), final_verdict and final_verdict_row match HEAD logic on all rows
    without requiring a light curve."""
    cfg_default = _cfg_mask(False)
    test_rows = [
        row(),
        row(flags="odd_even"),
        row(flags="odd_even,secondary", secondary_depth=1e-6),
        row(deep_flags="density_mismatch"),
        row(deep_flags="v_shape"),
        row(rp_rjup=2.8, rp_earth=31.0, deep_flags="rp_gt_2rjup"),
        row(flags="common_period"),
        row(checks_logged=False),
    ]
    for r in test_rows:
        v_none = vf.final_verdict(r)
        v_off = vf.final_verdict(r, cfg=cfg_default)
        v_empty = vf.final_verdict(r, cfg={})
        assert v_none == v_off == v_empty
        d_none = vf.final_verdict_row(r)
        d_off = vf.final_verdict_row(r, cfg=cfg_default)
        for k in d_none:
            val1, val2 = d_none[k], d_off[k]
            if isinstance(val1, float) and np.isnan(val1):
                assert np.isnan(val2)
            else:
                assert val1 == val2, k


def test_mask_on_without_lc_does_not_silently_fallback_and_paths_agree():
    """When the mask is on and no light curve is provided, final_verdict must not fall back to the
    unmasked odd_even flag: it reports STATUS_NEEDS_LC (failing with odd_even_needs_lc under 'reject',
    or raising ValueError under 'raise'), and real/null/inject paths agree."""
    cfg_on = _cfg_mask(True)
    r_unmasked_fail = row(flags="odd_even")
    r_unmasked_pass = row(flags="")

    for r in (r_unmasked_fail, r_unmasked_pass):
        verdicts = [p(r, cfg=cfg_on) for p in PATHS]
        assert verdicts[0] == verdicts[1] == verdicts[2]
        v = verdicts[0]
        assert v.status == vf.STATUS_NEEDS_LC
        assert not v.passed
        assert "odd_even" not in v.all_failing
        assert vf.ODD_EVEN_NEEDS_LC in v.all_failing

        with pytest.raises(ValueError, match="odd_even_mask is enabled"):
            vf.final_verdict(r, cfg=cfg_on, needs_lc="raise")


def test_real_null_inject_agree_when_mask_on_for_failing_odd_even():
    """When >= 2 odd and >= 2 even transits survive the mask and their depths differ by > 3 sigma and > 10%,
    all three paths (real, null, inject) reject with 'odd_even' identically."""
    dt = 0.005
    P, t0, D = 3.0, 1.0, 0.10
    t = np.arange(0.0, 14.0, dt)  # epochs 0 (1.0), 1 (4.0), 2 (7.0), 3 (10.0), 4 (13.0) -> 3 even, 2 odd
    f = np.ones_like(t)
    e = np.full(t.size, 1e-4)
    ep = np.round((t - t0) / P).astype(int)
    in_tr = (t >= (t0 + ep * P - 0.5 * D - 1e-9)) & (t < (t0 + ep * P + 0.5 * D - 1e-9))
    f[in_tr & (ep % 2 == 0)] -= 0.0020
    f[in_tr & (ep % 2 == 1)] -= 0.0040  # persistent odd/even difference across both odd transits

    cfg_on = _cfg_mask(True)
    for mode in ("real", "invert", "scramble", "inject"):
        r = row(mode=mode, P_grid=P, period=P, t0=t0, duration=D, flags="")
        verdicts = [p(r, cfg=cfg_on, lc=(t, f, e)) for p in PATHS]
        assert verdicts[0] == verdicts[1] == verdicts[2]
        assert not verdicts[0].passed
        assert verdicts[0].all_failing == ["odd_even"]
        assert verdicts[0].status == vf.STATUS_OK
