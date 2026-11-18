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


def test_comb_guard_false_comb_rate_on_pure_pulsators():
    """Comb guard false-comb rate on pure pulsators must be <= 5% (Spec A9)."""
    from exosearch.v11 import synth
    cfg = _cfg()
    rate, count = synth.false_comb_rate(cfg, range(50))
    assert rate <= 0.05, f"False comb rate {rate:.2f} ({count}/50) exceeds 5% threshold"


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


def test_candidate_order_ignores_last_bit_ties():
    periods = np.linspace(1, 10, 1001)
    p = np.zeros(1001)
    p[100] = 20.0                         # global peak
    p[700] = 10.0 * (1 + 1e-15)           # two disjoint, exactly tied (last bit) secondary peaks
    p[300] = 10.0
    c1, *_ = search.candidates([dict(D=0.1, p=p)], periods, 2, 0.01)
    p2 = p.copy(); p2[700], p2[300] = 10.0, 10.0 * (1 + 1e-15)
    c2, *_ = search.candidates([dict(D=0.1, p=p2)], periods, 2, 0.01)
    assert c1[1]["idx"] == c2[1]["idx"] == 300   # lowest period among the tie, on either host


# ---------------------------------------------------------------- T3 shared fixes (C1-C4)
import importlib.util
import subprocess
from pathlib import Path

import pandas as pd

from exosearch.v11 import registration_checks, scramble


def _exclusion_cfg(rel):
    c = settings.cfg()
    c["scramble"] = dict(c["scramble"], exclude_hosts_file={"70": rel})
    return c


def test_missing_exclusion_file_is_a_hard_error(monkeypatch):
    c = _exclusion_cfg("configs/does_not_exist.csv")
    with pytest.raises(FileNotFoundError):
        settings.exclusion_files(c)
    monkeypatch.setattr(settings, "cfg", lambda path=settings.CFG_FILE: c)
    with pytest.raises(FileNotFoundError):
        settings.provenance("h1", None, False)
    with pytest.raises(FileNotFoundError):
        settings.config_hash()
    with pytest.raises(FileNotFoundError):
        run.settings.exclusion_tics(70, c)


def test_exclusion_file_without_tic_column_is_a_hard_error(tmp_path, monkeypatch):
    (tmp_path / "ex.csv").write_text("host,x\n1,2\n")
    monkeypatch.setattr(settings, "ROOT", tmp_path)
    with pytest.raises(ValueError, match="tic"):
        settings.exclusion_files(_exclusion_cfg("ex.csv"))


def test_exclusion_content_sha_enters_config_hash(tmp_path, monkeypatch):
    (tmp_path / "ex.csv").write_text("tic\n1\n2\n")
    c = _exclusion_cfg("ex.csv")
    monkeypatch.setattr(settings, "ROOT", tmp_path)
    monkeypatch.setattr(settings, "cfg", lambda path=settings.CFG_FILE: dict(c))
    h1 = settings.config_hash()
    (tmp_path / "ex.csv").write_text("tic\n1\n2\n3\n")
    assert settings.config_hash() != h1
    assert settings.config()["scramble_exclude_sha256"]["70"]


def test_provenance_records_tree_next_to_commit():
    p = settings.provenance("h1", None, False)
    tree = subprocess.check_output(["git", "-C", str(settings.ROOT), "rev-parse", "HEAD^{tree}"], text=True).strip()
    assert p["tree"] == tree and p["commit"] != p["tree"] and "scramble_exclude_sha256_70" in p


def test_compare_reports_commit_tree_and_fails_on_tree_mismatch(monkeypatch):
    def rows(host, tree):
        return pd.DataFrame(dict(key=["a", "b"], status="ok", S=[1.0, 2.0], period=[1.0, 2.0], P_grid=[1.0, 2.0],
                                 duration=[0.1, 0.1], P_grid_halfwidth=[0.01, 0.01], host_id=host, arch="x",
                                 lib_versions="{}", commit=f"c-{host}", tree=tree, runtime_s=1.0))
    frames = {"d1": rows("h1", "T1"), "d2": rows("h2", "T1"), "d3": rows("h3", "T2")}
    monkeypatch.setattr(run, "load_rows", lambda d: frames[d])
    res = run.compare(["d1", "d2", "d3"])
    assert [h["tree"] for h in res["hosts"]] == ["T1", "T1", "T2"]
    assert [h["commit"] for h in res["hosts"]] == ["c-h1", "c-h2", "c-h3"]
    assert not res["tree_match"] and not res["all_pass"]
    res = run.compare(["d1", "d2"])
    assert res["tree_match"] and res["all_pass"]


def test_zero_rotation_offset_is_redrawn():
    t = np.concatenate([np.arange(0, 1, 0.01), np.arange(3, 4, 0.01)])
    f = np.arange(t.size, dtype=float)
    for seed in range(50):
        _, fo, _, info = scramble.segment_rotate(t, f, f.copy(), seed, dict(scramble=dict(gap_d=0.5, gap_trim_d=0.0)))
        assert all(o >= 1 for o in info["scr_offsets"])
    # two-sample segments: the only non-zero offset is 1, drawn every time
    t2 = np.array([0.0, 0.01, 5.0, 5.01, 10.0, 10.01])
    for seed in range(20):
        _, _, _, info = scramble.segment_rotate(t2, np.arange(6.0), np.arange(6.0), seed,
                                                dict(scramble=dict(gap_d=0.5, gap_trim_d=0.0)))
        assert info["scr_offsets"] == [1, 1, 1]


def test_registration_guard_is_a_dedicated_exception_and_rejects_empty_sets():
    flip = pd.DataFrame(dict(tic=[1, 2, 3], sector=70))
    scr = pd.DataFrame(dict(tic=[2, 3], sector=70))
    with pytest.raises(registration_checks.PopulationGuardError):
        registration_checks.check_registration_population(flip, scr, pd.DataFrame(dict(tic=[2])))
    with pytest.raises(registration_checks.PopulationGuardError, match="empty"):
        registration_checks.check_registration_population(flip, scr, pd.DataFrame(dict(tic=[], sector=[])))
    with pytest.raises(registration_checks.PopulationGuardError, match="empties"):
        registration_checks.check_registration_population(flip, scr, pd.DataFrame(dict(tic=[9], sector=[106])), sector=70)
    with pytest.raises(registration_checks.PopulationGuardError, match="tic"):
        registration_checks.check_registration_population(flip, scr, pd.DataFrame(dict(x=[1])))
    assert not issubclass(registration_checks.PopulationGuardError, AssertionError)


def _launcher():
    p = Path(__file__).resolve().parents[1] / "scripts" / "v11_launch.py"
    spec = importlib.util.spec_from_file_location("v11_launch_t", p)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def test_launcher_runs_guard_before_beta_select_and_stops_on_failure(tmp_path, monkeypatch):
    L = _launcher()
    monkeypatch.setattr(L, "WT", tmp_path / "wt")
    (tmp_path / "wt" / L.SUB).mkdir(parents=True)
    monkeypatch.setattr(L, "worktree", lambda *a, **k: None)
    monkeypatch.setattr(Path, "symlink_to", lambda self, target: None)
    calls = []

    def fake_sh(*cmd, **kw):
        calls.append(" ".join(cmd))
        if "registration_checks" in " ".join(cmd):
            raise RuntimeError("POPULATION GUARD FAILED")
        return ""
    monkeypatch.setattr(L, "sh", fake_sh)
    with pytest.raises(RuntimeError, match="POPULATION GUARD"):
        L.register({})
    assert any("registration_checks" in c for c in calls)
    assert not any("beta_select" in c for c in calls)          # registration stopped before the beta computation

    # main() turns the failure into the 'error' state and launches nothing
    monkeypatch.setattr(L, "STATUS", tmp_path / "status.json")
    monkeypatch.setattr(L, "OUTROOT", tmp_path)
    monkeypatch.setattr(L, "r0_done", lambda: True)
    launched = []
    monkeypatch.setattr(L, "launch", lambda *a, **k: launched.append(1))
    L.main(False)
    st = json.loads((tmp_path / "status.json").read_text())
    assert st["state"] == "error" and not launched


def test_launcher_orders_guard_then_beta_select(tmp_path, monkeypatch):
    L = _launcher()
    monkeypatch.setattr(L, "WT", tmp_path / "wt")
    (tmp_path / "wt" / L.SUB).mkdir(parents=True)
    monkeypatch.setattr(L, "worktree", lambda *a, **k: None)
    monkeypatch.setattr(Path, "symlink_to", lambda self, target: None)
    calls = []

    def fake_sh(*cmd, **kw):
        calls.append(" ".join(cmd))
        raise RuntimeError("stop")
    monkeypatch.setattr(L, "sh", fake_sh)
    with pytest.raises(RuntimeError):
        L.register({})
    assert "registration_checks" in calls[0]


# ---------------------------------------------------------------- v1.2 step 2: odd-even transit mask
def test_odd_even_mask_default_settings_in_config():
    c = _cfg()
    mcfg = vet11.odd_even_mask_cfg(c)
    assert mcfg["enabled"] is False
    assert mcfg["min_coverage"] == 0.75
    assert mcfg["min_gap_dump_dist_d"] == 0.375
    assert mcfg["min_gap_dump_dist_h"] == 9.0
    assert mcfg["gap_d"] == 0.5
    assert mcfg["oot_baseline_t14"] == 1.0
    assert mcfg["oot_one_side"] is True
    # Configs lacking the new keys behave as mask-off
    assert vet11.odd_even_mask_enabled({}) is False
    assert vet11.odd_even_mask_enabled(dict(vet=dict(odd_even_sigma=3.0, odd_even_rel=0.1))) is False


def test_odd_even_mask_coverage_boundary():
    # dt = 0.01 d, D = 0.20 d -> expected 20 cadences in [t_mid - 0.10, t_mid + 0.10)
    dt = 0.01
    P, t0, D = 3.0, 1.0, 0.20
    t_full = np.arange(0.0, 20.0, dt)
    # Epoch 0 at t_mid = 1.0 ([0.90, 1.10)): drop 6 of 20 cadences -> 14/20 = 0.70 < 0.75 (excluded)
    # Epoch 1 at t_mid = 4.0 ([3.90, 4.10)): drop 5 of 20 cadences -> 15/20 = 0.75 == X (exact boundary, kept)
    # Epoch 2 at t_mid = 7.0 ([6.90, 7.10)): drop 2 of 20 cadences -> 18/20 = 0.90 > 0.75 (kept)
    drop0 = (t_full >= 0.92 - 1e-9) & (t_full < 0.98 - 1e-9)   # 6 cadences removed (0.06 d < gap_d=0.5 d)
    drop1 = (t_full >= 3.92 - 1e-9) & (t_full < 3.97 - 1e-9)   # 5 cadences removed (15/20 = 0.75)
    drop2 = (t_full >= 6.92 - 1e-9) & (t_full < 6.94 - 1e-9)   # 2 cadences removed (18/20 = 0.90)
    t = t_full[~(drop0 | drop1 | drop2)]

    m = vet11.odd_even_transit_mask(t, P, t0, D, min_coverage=0.75)
    ep_map = {int(ep): i for i, ep in enumerate(m["epochs"])}
    assert np.isclose(m["coverage"][ep_map[0]], 0.70) and not m["pass_coverage"][ep_map[0]] and not m["keep"][ep_map[0]]
    assert np.isclose(m["coverage"][ep_map[1]], 0.75) and m["pass_coverage"][ep_map[1]] and m["keep"][ep_map[1]]
    assert np.isclose(m["coverage"][ep_map[2]], 0.90) and m["pass_coverage"][ep_map[2]] and m["keep"][ep_map[2]]


def test_odd_even_mask_gap_and_dump_distance_boundaries():
    dt = 0.005
    P, t0, D = 3.0, 1.0, 0.10
    Y = 0.375  # 9 hours
    t_full = np.arange(0.0, 22.0, dt)
    # Create a data gap > 0.5 d whose right gap edge is at t = 4.0 - Y = 3.625 (exact distance Y from epoch 1 at 4.0)
    # and left gap edge is at 3.625 - 0.600 = 3.025 (> 0.5 d gap).
    # Create another gap whose left gap edge is at 7.0 + (Y - 0.005) = 7.370 (distance 0.370 < Y from epoch 2 at 7.0).
    gap1 = (t_full > 3.025 + 1e-9) & (t_full < 3.625 - 1e-9)
    gap2 = (t_full > 7.370 + 1e-9) & (t_full < 8.000 - 1e-9)
    t = t_full[~(gap1 | gap2)]

    # Momentum dump at exact distance Y from epoch 3 (t_mid = 10.0 -> dump at 10.375)
    # and closer than Y from epoch 4 (t_mid = 13.0 -> dump at 13.370, distance 0.370 < Y)
    dumps = np.array([10.375, 13.370])

    m = vet11.odd_even_transit_mask(t, P, t0, D, dumps=dumps, min_gap_dump_dist_d=Y)
    ep_map = {int(ep): i for i, ep in enumerate(m["epochs"])}

    # Epoch 1 (t_mid = 4.0): nearest gap edge at 3.625 -> distance == 0.375 == Y -> KEPT
    assert np.isclose(m["gap_dist"][ep_map[1]], Y) and m["pass_gap"][ep_map[1]] and m["keep"][ep_map[1]]
    # Epoch 2 (t_mid = 7.0): nearest gap edge at 7.370 -> distance == 0.370 < Y -> EXCLUDED
    assert np.isclose(m["gap_dist"][ep_map[2]], 0.370) and not m["pass_gap"][ep_map[2]] and not m["keep"][ep_map[2]]
    # Epoch 3 (t_mid = 10.0): dump at 10.375 -> distance == 0.375 == Y -> KEPT
    assert np.isclose(m["dump_dist"][ep_map[3]], Y) and m["pass_dump"][ep_map[3]] and m["keep"][ep_map[3]]
    # Epoch 4 (t_mid = 13.0): dump at 13.370 -> distance == 0.370 < Y -> EXCLUDED
    assert np.isclose(m["dump_dist"][ep_map[4]], 0.370) and not m["pass_dump"][ep_map[4]] and not m["keep"][ep_map[4]]


def test_odd_even_mask_oot_baseline_one_side_and_both_sides():
    dt = 0.01
    P, D = 2.0, 0.10
    # Epoch 0 at t0 = 0.15 -> ingress = 0.10, egress = 0.20.
    # Start light curve at t = 0.00 -> exactly 1.0 * T14 = 0.10 d before ingress, and Plenty after egress.
    # If we instead start at t = 0.05 -> only 0.5 * T14 before ingress, but plenty after egress (1 side has >= 1 T14).
    t_one_side = np.arange(0.05, 10.0, dt)
    m_def = vet11.odd_even_transit_mask(t_one_side, P, 0.15, D, oot_baseline_t14=1.0, oot_one_side=True)
    assert not m_def["pass_oot_before"][0] and m_def["pass_oot_after"][0]
    assert m_def["pass_oot"][0] and m_def["keep"][0]

    # When both sides are required (oot_one_side=False), epoch 0 is excluded
    m_both = vet11.odd_even_transit_mask(t_one_side, P, 0.15, D, oot_baseline_t14=1.0, oot_one_side=False)
    assert not m_both["pass_oot"][0] and not m_both["keep"][0]

    # A transit with < 1.0 * T14 on BOTH sides (e.g. isolated segment [0.05, 0.25) around [0.10, 0.20)) is excluded
    t_neither = np.arange(0.05, 0.25, dt)
    m_neither = vet11.odd_even_transit_mask(t_neither, P, 0.15, D, oot_baseline_t14=1.0, oot_one_side=True)
    assert not m_neither["pass_oot_before"][0] and not m_neither["pass_oot_after"][0]
    assert not m_neither["pass_oot"][0] and not m_neither["keep"][0]


def test_screens_mask_on_without_lc_raises():
    # Search-time screens() must not silently fall back to the unmasked BLS odd-even flag.
    c = _cfg()
    m = dict(c["vet"]["odd_even_mask"], enabled=True)
    c["vet"] = dict(c["vet"], odd_even_mask=m)
    with pytest.raises(ValueError, match="odd_even_mask is enabled"):
        vet11.screens(_st(1.2e-3, 0.8e-3, 1e-5), _best(), {"radius": 1.0}, c)
