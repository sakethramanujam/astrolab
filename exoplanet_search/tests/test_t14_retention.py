"""v1.2 step 4: T14 retention diagnostic (reports only)."""
import copy

import numpy as np
import pytest

from exosearch import t14_retention as tr
from exosearch import vet_final as vf

P, DEPTH, SIG = 9.3, 3e-3, 5e-4


def _lc(seed=0):
    t = np.arange(0, 27.0, 2 / 1440.0) + 3000.0
    t = t[(t < 3013.0) | (t > 3014.0)]
    f = 1 + np.random.default_rng(seed).normal(0, SIG, t.size)
    return t, f, np.full(t.size, SIG)


def _case(t14_h, seed=0, t0=3004.1):
    t, f, e = _lc(seed)
    return dict(t=t, f=tr.plant_box(t, f, t0, P, t14_h, DEPTH), e=e, t0=t0, period=P, t14_h=t14_h, depth=DEPTH)


def test_bin_classification():
    short, mid, long_ = tr.bin_labels()
    assert (short, mid, long_) == ("T14 < 3 h", "3 <= T14 <= 6 h", "T14 > 6 h")
    assert tr.t14_bin(1.5) == short and tr.t14_bin(2.99) == short
    assert tr.t14_bin(3.0) == mid and tr.t14_bin(4.5) == mid and tr.t14_bin(6.0) == mid
    assert tr.t14_bin(6.01) == long_ and tr.t14_bin(12.0) == long_


def test_diagnostic_breakdown_and_qualitative_inequality():
    cases = [_case(h, seed=s) for h in (1.5, 2.5, 4.0, 12.0) for s in (0, 1)]
    out = tr.t14_retention_diagnostic(cases)
    short, mid, long_ = out["labels"]
    b = out["bins"]
    assert [b[k]["n"] for k in (short, mid, long_)] == [4, 2, 2]
    assert all(b[k]["n_valid"] == b[k]["n"] for k in b)
    # short transit on a quiet light curve: nearly all depth kept; long (12 h = 0.5 d of a 0.75 d window): much less
    assert b[short]["mean_retention"] > 0.85
    assert b[long_]["mean_retention"] < 0.6
    assert b[short]["mean_retention"] - b[long_]["mean_retention"] > 0.3
    assert b[short]["frac_kept"] == 1.0 and b[long_]["frac_kept"] == 0.0
    assert b[mid]["mean_retention"] > b[long_]["mean_retention"]
    assert out["tolerance"] == 0.5 and out["edges_h"] == [3.0, 6.0]
    assert [r[1] for r in out["retentions"]][:2] == [short, short]


def test_box_depth_estimator_exact_on_flat_noiseless():
    t = np.arange(0, 20.0, 0.01)
    f = np.ones_like(t)
    f[np.abs(((t - 5.0 + 5.0) % 10.0) - 5.0) < 0.05] -= 0.01          # T14 = 2.4 h
    d = tr.box_depth(t, f, np.full(t.size, 1e-3), 5.0, 10.0, 2.4)
    assert d == pytest.approx(0.01)
    assert np.isnan(tr.box_depth(t, f, np.full(t.size, 1e-3), 5.0, 10.0, 0.0))


def test_empty_bin_and_no_verdict_fields():
    out = tr.t14_retention_diagnostic([_case(2.0)])
    assert out["bins"][out["labels"][2]]["mean_retention"] is None
    assert out["bins"][out["labels"][2]]["n"] == 0
    assert not any(k in out for k in ("passed", "verdict", "final_verdict", "flags", "s_star", "S"))


def test_diagnostic_does_not_change_final_verdict():
    row = dict(tic=1, sector=70, status="ok", mode="real", task="real", S=20.0, P_grid=3.0, period=3.0,
               depth=0.01, t0=0.5, duration=0.1, flags="", deep_flags="", rp_rjup=1.0, rp_earth=11.0,
               secondary_depth=0.0, sec_scan_sigma=1.2, checks_logged=True, matched=True)
    before_row = copy.deepcopy(row)
    v0 = vf.final_verdict(row)
    tr.t14_retention_diagnostic([_case(2.0), _case(12.0)])
    assert row == before_row
    assert vf.final_verdict(row) == v0
    assert vf.final_verdict_row(dict(row)) == vf.final_verdict_row(dict(before_row))


def test_cfg_defaults():
    assert tr.t14_cfg(None) == dict(bin_edges_h=[3.0, 6.0], tolerance=0.5)
    assert tr.t14_cfg({"t14_retention": {"bin_edges_h": [2, 5]}})["bin_edges_h"] == [2.0, 5.0]
