import json

import numpy as np
import pytest

from exosearch.v11 import grid, prewhiten, run, search, settings, vet11


def _cfg():
    return settings.cfg()


# ---------------------------------------------------------------- grid and duration cap
def test_grid_ascending_same_span_and_formula():
    t = np.linspace(3000, 3027, 5000)
    g = grid.period_grid(t, search.DURATIONS)
    assert np.all(np.diff(g) > 0) and g[0] >= 0.5 and g[-1] <= 13.5 + 1e-9
    base = 27.0
    df = search.DURATIONS.min() / base ** 2
    f = (1 / 13.5) + df * np.arange(len(g), dtype=float)
    assert np.array_equal(g, (1.0 / f[f < 2.0])[::-1])
    assert grid.grid_sha256(g) == grid.grid_sha256(g.copy())


def test_q_circ_sun_and_cap():
    # Sun at 365.25 d: T14 ~ 13 h
    T = grid.q_circ(365.25, 1.0) * 365.25 * 24
    assert 12 < T < 14
    assert np.isclose(grid.d_max(2.0, None), 0.2)
    assert np.isclose(grid.d_max(2.0, 1.0), 3 * grid.q_circ(2.0, 1.0) * 2.0)
    P = np.array([0.5, 1, 5, 10])
    m = grid.allowed(P, 0.1, 50.0)           # dense M dwarf: 2.4-h boxes not allowed at short P
    assert not m[0] and m[-1]


# ---------------------------------------------------------------- beta options
def test_beta_options():
    sb = search.sigma_beta(18000, 2 / 1440, 0.1)
    assert np.isclose(sb, 1 / np.sqrt(2 * 18000 * (2 / 1440) / 0.1))
    assert search.beta_prime("a", 0.8, sb) == 0.8
    assert search.beta_prime("b", 0.8, sb) == 1.0 and search.beta_prime("b", 1.2, sb) == 1.2
    assert np.isclose(search.beta_prime("c", 0.5, sb), 1 - 2 * sb)
    d = search.beta_prime("d", 0.8, sb, (1.0, 0.04))
    assert 0.8 < d < 1.0
    with pytest.raises(ValueError):
        search.beta_prime("d", 0.8, sb, None)


def test_production_needs_beta_option():
    c = dict(_cfg(), beta_option=None)
    with pytest.raises(SystemExit):
        run.primary_option(c)
    assert run.primary_option(c, "c", dry_run=True) == "c"
    assert run.primary_option(dict(c, beta_option="d")) == "d"
    with pytest.raises(SystemExit):
        run.primary_option(dict(c, beta_option="d"), "a")


def test_config_hash_depends_on_beta_option():
    assert settings.config_hash("a") != settings.config_hash("c")


# ---------------------------------------------------------------- exact depth and momentum dumps
def _box_lc(P=3.0, t0=3001.0, D=0.1, depth=2e-3, noise=5e-4, seed=0):
    rng = np.random.default_rng(seed)
    t = np.arange(3000, 3027, 10 / 1440)
    f = 1 + rng.normal(0, noise, t.size)
    f[search.in_transit(t, P, t0, D)] -= depth
    return t, f, np.full(t.size, noise)


def test_exact_box_matches_compute_stats():
    from astropy.timeseries import BoxLeastSquares
    t, f, e = _box_lc()
    r = search.exact_box(t, f, e, 3.0, 3001.0, 0.1)
    st = BoxLeastSquares(t, f, dy=e).compute_stats(3.0, 0.1, 3001.0)
    assert np.isclose(r["depth"], st["depth"][0], rtol=1e-9)
    assert r["n_dump_epochs"] == 0 and np.isclose(r["n_tr"], r["n_tr_all"])


def test_dump_epochs_removed():
    t, f, e = _box_lc()
    r0 = search.exact_box(t, f, e, 3.0, 3001.0, 0.1)
    r = search.exact_box(t, f, e, 3.0, 3001.0, 0.1, dumps=np.array([3004.05]), window=0.1)
    assert r["n_dump_epochs"] == 1 and r["n_tr"] < r0["n_tr"]


def test_dump_times_from_quality():
    t = np.arange(0, 10, 0.01)
    q = np.zeros(t.size, int)
    q[100:103] = 32
    q[500:502] = 32 | 8
    d = vet11.dump_times(t, q)
    assert len(d) == 2 and abs(d[0] - 1.01) < 1e-9


# ---------------------------------------------------------------- plateau tie-break
def test_plateau_centre_and_lowest_disjoint():
    periods = np.linspace(1, 2, 101)
    p = np.zeros(101)
    p[40:45] = 10.0           # plateau of width 5
    p[80] = 10.0              # disjoint peak at higher period
    per_d = [dict(D=0.1, p=p)]
    plats, inb = search.plateaus_in_band(per_d, periods)
    assert len(plats) == 2 and plats[0]["width"] == 5 and plats[0]["centre"] == 42
    step = periods[1] - periods[0]
    assert np.isclose(plats[0]["hw"], (periods[44] - periods[40] + step) / 2)
    cands, *_ = search.candidates(per_d, periods, 1, 0.01)
    assert cands[0]["idx"] == 42


# ---------------------------------------------------------------- prewhitening and comb guard
def test_prewhiten_removes_pulsation():
    rng = np.random.default_rng(1)
    t = np.arange(3000, 3027, 2 / 1440)
    f = 1 + 3e-3 * np.sin(2 * np.pi * 22.97 * t) + rng.normal(0, 1e-3, t.size)
    r, log = prewhiten.prewhiten(t, f, _cfg()["prewhiten"])
    terms = json.loads(log["prewhiten_terms"])
    assert log["prewhiten_n"] >= 1 and abs(terms[0][0] - 22.97) < 0.01
    assert np.std(r) < 1.3e-3


def test_comb_guard_protects_short_period_transit():
    rng = np.random.default_rng(2)
    t = np.arange(3000, 3027, 2 / 1440)
    f = 1 + rng.normal(0, 3e-4, t.size)
    f[search.in_transit(t, 0.5, 3000.1, 1 / 24)] -= 5e-3         # USP-like box: harmonics above 5 c/d
    _, on = prewhiten.prewhiten(t, f, _cfg()["prewhiten"], comb_guard=True)
    _, off = prewhiten.prewhiten(t, f, _cfg()["prewhiten"], comb_guard=False)
    assert on["prewhiten_comb"] >= 1 and on["prewhiten_n"] < off["prewhiten_n"]


def test_box_ladder_small_duty_cycle():
    r2, r3 = prewhiten.box_ladder(1, np.array([0.01]))
    assert r2[0] > 0.99 and r3[0] > 0.99
    assert prewhiten.is_comb(6.0, 1.0, 0.95, 0.9)
    assert not prewhiten.is_comb(6.0, 1.0, 0.05, 0.02)


# ---------------------------------------------------------------- vetting rule
def _st(d_odd, d_even, e, d_sec=0.0, e_sec=1e-4, depth=1e-3):
    return dict(depth=(depth, 1e-5), depth_odd=(d_odd, e), depth_even=(d_even, e), depth_phased=(d_sec, e_sec),
                per_transit_count=np.ones(5))


def _best(depth=1e-3):
    return dict(depth_exact=depth, P_grid=3.0, duration=0.1, n_dump_epochs=0, n_epochs=5)


def test_odd_even_needs_significance_and_size():
    c = _cfg()
    sig_small = vet11.screens(_st(1.02e-3, 0.98e-3, 1e-6), _best(), {"radius": 1.0}, c)   # 28 sigma, 4%
    assert "odd_even" not in sig_small["flags"]
    both = vet11.screens(_st(1.2e-3, 0.8e-3, 1e-5), _best(), {"radius": 1.0}, c)         # 28 sigma, 40%
    assert "odd_even" in both["flags"]
    big_insig = vet11.screens(_st(1.2e-3, 0.8e-3, 2e-4), _best(), {"radius": 1.0}, c)    # 1.4 sigma
    assert "odd_even" not in big_insig["flags"]


def test_secondary_needs_significance_and_size():
    c = _cfg()
    small = vet11.screens(_st(1e-3, 1e-3, 1e-5, d_sec=5e-5, e_sec=1e-5), _best(), {"radius": 1.0}, c)   # 5 sigma, 5%
    assert "secondary" not in small["flags"]
    big = vet11.screens(_st(1e-3, 1e-3, 1e-5, d_sec=2e-4, e_sec=1e-5), _best(), {"radius": 1.0}, c)     # 20 sigma, 20%
    assert "secondary" in big["flags"]


def test_momentum_dump_reject_and_warn():
    c = _cfg()
    st = _st(1e-3, 1e-3, 1e-5)
    st["per_transit_count"] = np.ones(2)
    r = vet11.screens(st, dict(_best(), n_dump_epochs=1, n_epochs=2), {"radius": 1.0}, c)
    assert "momentum_dump" in r["flags"] and r["period_ambiguous"]
    w = vet11.screens(_st(1e-3, 1e-3, 1e-5), dict(_best(), n_dump_epochs=1, n_epochs=5), {"radius": 1.0}, c)
    assert "momentum_dump" not in w["flags"] and "momentum_dump_epoch" in w["warnings"] and not w["period_ambiguous"]


# ---------------------------------------------------------------- end to end on a synthetic light curve
def test_search_recovers_transit_with_oot_term(synthetic_lc):
    from exosearch import lightcurve
    from exosearch.v11 import kernels
    kernels.install()
    t, f, e, truth = synthetic_lc
    tf, ff, ef, *_ = lightcurve.detrend(t, f, e)
    tb, fb, eb = lightcurve.bin_lc(tf, ff, ef)
    best, sens = search.search(tb, fb, eb, tf, ff, rho=1.0, cfg=_cfg(), options=("a", "b", "c"), primary="a")
    assert abs(best["period"] / truth["period"] - 1) < 0.01 and best["S"] > 15
    assert np.isclose(best["S"], best["S_exact"] * np.sqrt(1 - best["duration"] / best["P_grid"]))
    assert sens["b"]["S"] <= sens["a"]["S"] * (1 + 1e-12) or best["beta_d"] < 1
    assert best["P_grid_halfwidth"] > 0


def test_shard_assignment_stable():
    w = [("h1", 0.2), ("h2", 0.7), ("h3", 0.1)]
    hs = [run.assign(t, w) for t in range(1000, 6000)]
    assert run.assign(1234, w) == run.assign(1234, w)
    share = {h: hs.count(h) / len(hs) for h in ("h1", "h2", "h3")}
    assert abs(share["h2"] - 0.7) < 0.03 and abs(share["h3"] - 0.1) < 0.02
