"""v1.2 step 4: S* on the 0.005 grid (default off; v1.1 GPD path unchanged)."""
import json
import math
from pathlib import Path

import numpy as np
import pytest

from exosearch import interim
from exosearch import sstar_grid as sg

CFG = json.loads((Path(__file__).resolve().parents[1] / "configs" / "v1.1.json").read_text())


def _sample(seed=0, n=400):
    rng = np.random.default_rng(seed)
    return 6.0 + rng.exponential(1.2, n)


def test_grid_shape():
    g = sg.s_grid()
    assert g[0] == 6.0 and g[-1] == 30.0
    assert len(g) == 4801
    assert np.allclose(np.diff(g), 0.005)
    assert np.all(np.round(g, 3) == g)
    assert np.allclose(g / 0.005, np.round(g / 0.005))


@pytest.mark.parametrize("seed,B", [(0, 3.0), (1, 1.0), (2, 7.9), (3, 0.0), (4, 25.0)])
def test_discrete_rule(seed, B):
    S = _sample(seed)
    s = sg.sstar_discrete(S, B)
    g = sg.s_grid()
    assert s in g
    cap = math.floor(B)
    assert int((S >= s).sum()) <= cap
    i = int(np.argmax(g == s))
    if i > 0:
        assert int((S >= g[i - 1]).sum()) > cap


def test_discrete_first_grid_point():
    assert sg.sstar_discrete([6.0, 6.0, 7.0], 5.0) == 6.0
    assert sg.sstar_discrete([], 0.0) == 6.0


def test_crossing_between_005_spaced_values_resolved_finely():
    # B=1 survivor allowed: the 2nd-highest value 10.012 must be excluded -> S* = first grid point above it.
    S = np.r_[np.linspace(6.0, 9.0, 50), 10.012, 10.5]
    s = sg.sstar_discrete(S, 1.0)
    assert s == pytest.approx(10.015)          # smallest 0.005 multiple > 10.012 that leaves 1 survivor
    assert s < 10.05                           # not forced up to the next 0.05 grid value (10.05)
    assert sg.sstar_discrete(S, 0.0) == pytest.approx(10.505)


def test_budget_float_noise():
    S = np.r_[np.full(3, 8.0), 9.0]
    assert sg.sstar_discrete(S, 0.1 * 30) == sg.sstar_discrete(S, 3.0)       # 0.1*30 = 3.0000000000000004


def test_no_solution_raises_and_bad_budget():
    with pytest.raises(sg.SStarGridError):
        sg.sstar_discrete([31.0, 32.0, 33.0], 1.0)
    with pytest.raises(sg.SStarGridError):
        sg.sstar_discrete([7.0], -1.0)
    with pytest.raises(sg.SStarGridError):
        sg.sstar_discrete([7.0], float("nan"))


def test_nonfinite_ignored():
    S = np.r_[_sample(0), np.nan, -np.inf]
    assert sg.sstar_discrete(S, 3.0) == sg.sstar_discrete(_sample(0), 3.0)


def test_interpolation_is_guarded_and_optional():
    for seed in range(5):
        S = _sample(seed)
        for B in (1.5, 3.0, 6.7):
            d = sg.sstar_discrete(S, B)
            r = sg.sstar_grid(S, B, interpolate=True)
            assert r["s_star"] >= d
            assert r["n_survivors"] <= math.floor(B)
            assert r["s_star_discrete"] == d
            assert r["s_star"] in sg.s_grid()
    off = sg.sstar_grid(_sample(0), 3.0)
    assert off["s_star_interp"] is None and off["s_star"] == off["s_star_discrete"] and off["rule"] == "grid_discrete"


def test_pava():
    assert np.allclose(sg.pava_nonincreasing([3, 1, 2, 0]), [3, 1.5, 1.5, 0])
    y = np.array([5.0, 4.0, 4.5, 2.0, 2.5, 1.0])
    out = sg.pava_nonincreasing(y)
    assert np.all(np.diff(out) <= 1e-12)


# --- default off / v1.1 path unchanged ------------------------------------------------------------------
def test_defaults_off_in_config_and_missing_keys():
    c = sg.sstar_cfg(CFG)
    assert c["fine_grid_enabled"] is False and c["interpolate"] is False
    assert (c["grid_min"], c["grid_max"], c["grid_step"]) == (6.0, 30.0, 0.005)
    assert sg.sstar_cfg(None)["fine_grid_enabled"] is False
    assert sg.sstar_cfg({})["fine_grid_enabled"] is False
    assert sg.sstar_cfg({"sstar": {}})["fine_grid_enabled"] is False
    assert sg.sstar_cfg({"sstar": {"fine_grid_enabled": False}})["fine_grid_enabled"] is False
    assert sg.sstar_grid_enabled({"sstar": {"fine_grid_enabled": True}}) is True
    assert CFG["t14_retention"]["bin_edges_h"] == [3.0, 6.0]


@pytest.mark.parametrize("cfg", [None, {}, {"sstar": {}}, {"sstar": {"fine_grid_enabled": False}}, CFG])
def test_off_matches_gpd(cfg):
    S = _sample(7, 600) + 0.3
    direct = interim.gpd_threshold(S, 5000, 1.0 / 5000, 6.0)
    via = sg.select_sstar(S, 5000, 1.0 / 5000, 6.0, cfg=cfg)
    assert via == direct and direct["s_star"] is not None


def test_on_uses_grid():
    S = _sample(7, 600)
    cfg = {"sstar": {"fine_grid_enabled": True}}
    r = sg.select_sstar(S, 600, 5.0 / 600, 6.0, cfg=cfg)
    assert r["rule"] == "grid_discrete" and r["s_star"] == sg.sstar_discrete(S, 5.0)
    assert r["s_star"] in sg.s_grid()


def test_default_constant_is_off():
    assert sg.DEFAULT_SSTAR["fine_grid_enabled"] is False
