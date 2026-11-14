"""Unified final verdict (exosearch.vet_final): real, null and inject rows are scored by ONE function.

The same row dict must give the same verdict through the real, null and inject paths of
runs/v1.1_s70_x86/finish.py, and the amended rules must hold in every path: density_mismatch alone
never rejects, and rp ~ 25 R_Earth is not rejected (the old rp_gt_2rjup cap is dropped for the
30 R_Earth cap).
"""
import importlib.util
from pathlib import Path

import numpy as np

from exosearch import vet_final as vf

_FINISH = Path(__file__).resolve().parents[1] / "runs" / "v1.1_s70_x86" / "finish.py"
_spec = importlib.util.spec_from_file_location("s70_finish", _FINISH)
fin = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(fin)

BASE = dict(tic=1, sector=70, status="ok", mode="real", task="real", S=20.0, P_grid=3.0, period=3.0,
            depth=0.01, t0=0.5, duration=0.1, flags="", deep_flags="", rp_rjup=1.0, rp_earth=11.0,
            secondary_depth=0.0, checks_logged=True, matched=True)


def row(**kw):
    return dict(BASE, **kw)


PATHS = (fin.real_path, fin.null_path, fin.inject_path)


def test_same_verdict_through_all_three_paths():
    cases = [
        row(),                                                  # clean pass
        row(mode="invert", task="invert"),                      # null flavour of the same signal
        row(mode="inject", task="inject_run0", matched=True),   # inject flavour
        row(deep_flags="density_mismatch"),                     # (c): warn only -> PASS
        row(rp_rjup=25.0 / 11.209, rp_earth=25.0, deep_flags="rp_gt_2rjup"),   # ~25 R_Earth: passes
        row(rp_rjup=2.8, rp_earth=31.0, deep_flags="rp_gt_2rjup"),            # > 30 R_Earth: rejects
        row(deep_flags="v_shape"),                              # unchanged check still rejects
        row(deep_flags="secondary_scan"),                       # needs_lc (no light curve rebuilt)
        row(checks_logged=False),                               # checks never logged -> reject
        row(flags="odd_even"),                                  # screen flag still rejects
    ]
    for r in cases:
        verdicts = [p(r) for p in PATHS]
        assert verdicts[0] == verdicts[1] == verdicts[2], r.get("deep_flags")
        assert verdicts[0].status in (vf.STATUS_OK, vf.STATUS_NEEDS_LC)


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
    r = row(deep_flags="secondary_scan")
    lo, hi = vf.final_verdict(r), vf.final_verdict(r, needs_lc="exempt")
    assert lo.status == vf.STATUS_NEEDS_LC and not lo.passed and lo.all_failing == ["secondary_scan"]
    assert hi.status == vf.STATUS_NEEDS_LC and hi.passed and hi.all_failing == []
    # only the undecided scan rule is waived; other rules still reject in the upper bound
    hi2 = vf.final_verdict(row(deep_flags="secondary_scan,v_shape"), needs_lc="exempt")
    assert not hi2.passed and hi2.all_failing == ["v_shape"]
    # an explicit scan resolves the status (the caller rebuilt the light curve)
    scan = dict(sec_scan_sigma=1.0, sec_scan_phase=0.6, sec_scan_depth=0.0, duration_fine=0.1,
                t0_fine=0.5, loo_min_sigma=9.0, loo_min_orbit=0)
    ok = vf.final_verdict(row(deep_flags="secondary_scan"), scan=scan)
    assert ok.status == vf.STATUS_OK and ok.passed
    # every path brackets the same row identically
    for p in PATHS:
        assert p(r).status == vf.STATUS_NEEDS_LC and not p(r).passed
        assert p(r, needs_lc="exempt").passed
