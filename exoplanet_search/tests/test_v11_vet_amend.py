"""Vetter amendment 2026-11-12 (exosearch.v11.vet_amend): rules (a)-(d)."""
import numpy as np

from exosearch import deepvet
from exosearch.v11 import vet_amend as va

BASE = dict(P_grid=3.0, depth=0.01, flags="", deep_flags="", rp_rjup=1.0, rp_earth=11.0, secondary_depth=0.0,
            duration=0.1)


def row(**kw):
    return dict(BASE, **kw)


def test_unchanged_checks_still_reject():
    assert not va.amended(row(deep_flags="v_shape"))["passed_lc_checks_v2"]
    assert not va.amended(row(flags="odd_even"))["passed_lc_checks_v2"]
    assert va.amended(row())["passed_lc_checks_v2"]
    assert va.REJECT_V2 == (set(deepvet.REJECT) - {"rp_gt_2rjup", "density_mismatch"}) | {"rp_gt_30re"}


def test_a_albedo_exemption():
    # a/R* = 8 from the duration fallback; (Rp/a)^2 = 0.01 / 64; albedo 0.5 -> depth_sec = 0.5 * 0.01 / 64
    P, dur = 3.0, 3.0 / (np.pi * 8.0)
    low = dict(sec_scan_depth=0.5 * 0.01 / 64, duration_fine=dur, loo_min_sigma=9.0)
    high = dict(low, sec_scan_depth=2.0 * 0.01 / 64)
    r = row(P_grid=P, deep_flags="secondary_scan")
    out = va.amended(r, low)
    assert out["passed_lc_checks_v2"] and np.isclose(out["albedo_scan_v2"], 0.5)
    assert not va.amended(r, high)["passed_lc_checks_v2"]
    # Rp >= 30 R_Earth: no exemption
    assert not va.amended(row(P_grid=P, deep_flags="secondary_scan", rp_rjup=2.8), low)["passed_lc_checks_v2"]
    # phase-0.5 screen secondary uses the row's secondary_depth
    rs = row(P_grid=P, flags="secondary", duration=dur, secondary_depth=0.5 * 0.01 / 64)
    assert va.amended(rs)["passed_lc_checks_v2"]
    assert not va.amended(dict(rs, secondary_depth=3 * 0.01 / 64))["passed_lc_checks_v2"]


def test_a_rs_from_density():
    a, src = va.a_over_rstar(3.0, 1.0)                       # solar density (TIC rho is in solar units)
    assert src == "tic_density" and 8.6 < a < 8.9          # 215 * (3 / 365.25)^(2/3) = 8.75


def test_b_radius_cap_30_replaces_both_caps():
    assert va.amended(row(rp_rjup=2.5, deep_flags="rp_gt_2rjup", flags="rp_too_large"))["passed_lc_checks_v2"]
    out = va.amended(row(rp_rjup=2.8, deep_flags="rp_gt_2rjup"))
    assert not out["passed_lc_checks_v2"] and "rp_gt_30re" in out["reject_v2"]
    assert not va.amended(row(rp_rjup=np.nan, rp_earth=31.0, flags="rp_too_large"))["passed_lc_checks_v2"]


def test_c_density_mismatch_is_a_flag():
    out = va.amended(row(deep_flags="density_mismatch"))
    assert out["passed_lc_checks_v2"] and "density_mismatch" in out["warn_v2"]


def _lc(one_orbit_only, seed=1):
    rng = np.random.default_rng(seed)
    t = np.arange(0.0, 24.0, 2.0 / 1440.0)
    P, t0, dur = 3.0, 0.5, 0.1
    f = 1.0 + rng.normal(0, 3e-4, t.size)
    ph = ((t - t0) / P) % 1.0
    f[np.abs(((t - t0 + 0.5 * P) % P) - 0.5 * P) < dur / 2] -= 0.01
    sec = np.abs(ph - 0.6) < (dur / P) / 2
    if one_orbit_only:
        sec &= (t > t0 + 3 * P) & (t < t0 + 4 * P)
        f[sec] -= 0.004
    else:
        f[sec] -= 0.0008
    return t, f, np.full(t.size, 3e-4), P, t0


def test_d_leave_one_orbit_out():
    t, f, e, P, t0 = _lc(one_orbit_only=True)
    sc = va.scan_depth(t, f, e, P, t0)
    assert sc["sec_scan_sigma"] > 5
    loo = va.loo_scan_min(t, f, e, P, t0, sc["t0_fine"])
    assert loo["loo_min_sigma"] < 5 and loo["loo_min_orbit"] == 3
    t, f, e, P, t0 = _lc(one_orbit_only=False)
    sc2 = va.scan_depth(t, f, e, P, t0)
    loo2 = va.loo_scan_min(t, f, e, P, t0, sc2["t0_fine"])
    assert sc2["sec_scan_sigma"] > 5 and loo2["loo_min_sigma"] >= 5
    r = row(deep_flags="secondary_scan", rp_rjup=2.5)                      # large albedo: (a) does not apply
    assert va.amended(r, dict(sc, **loo, sec_scan_depth=0.004))["passed_lc_checks_v2"]
    assert not va.amended(r, dict(sc2, **loo2, sec_scan_depth=0.004))["passed_lc_checks_v2"]


def test_gate_uses_amended_verdict_when_present():
    import pandas as pd
    from exosearch.v11 import known_planets
    d = pd.DataFrame([dict(tic=1, period=3.0, S=20.0, status="ok", passed_lc_checks=0.0, passed_lc_checks_v2=1.0)])
    assert known_planets._search_index(d)[1]["vetted"] is True
    assert known_planets._search_index(d.drop(columns="passed_lc_checks_v2"))[1]["vetted"] is False
