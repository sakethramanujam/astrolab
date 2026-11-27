"""v1.2 step 4: S* evaluated on a fine grid (default OFF; the v1.1 GPD path is untouched).

The v1.1 S* comes from `interim.gpd_threshold` (used by `v11.beta_select`).  This module adds an *optional*
discrete alternative on the grid `S_grid = np.arange(6.0, 30.0 + 1e-12, 0.005)` (4801 points, values rounded to
3 decimals).  Nothing here is called unless `sstar.fine_grid_enabled` is true in the config; `select_sstar` with the
setting off (missing keys, false) is a plain call to `interim.gpd_threshold`.  No false-alarm budget is defined
here: the caller passes the survivor budget `B` (e.g. `p * n_null`).
"""
from __future__ import annotations

import math
from typing import Any, Dict, Optional

import numpy as np

GRID_MIN, GRID_MAX, GRID_STEP = 6.0, 30.0, 0.005
DECIMALS = 3
_BUDGET_EPS = 1e-9          # absorbs float noise in B = p * n_null (e.g. 0.1 * 30 = 3.0000000000000004)

DEFAULT_SSTAR = dict(fine_grid_enabled=False, grid_min=GRID_MIN, grid_max=GRID_MAX, grid_step=GRID_STEP,
                     interpolate=False)


class SStarGridError(ValueError):
    """No grid point satisfies the survivor budget (or the inputs are unusable)."""


def sstar_cfg(cfg: Optional[Dict[str, Any]] = None, sstar: Optional[Any] = None) -> Dict[str, Any]:
    """Resolve the S*-grid settings from `cfg["sstar"]` (and an optional bool / dict override).

    Configs without the keys resolve to `DEFAULT_SSTAR`: grid off, i.e. the v1.1 GPD S*.
    """
    out = dict(DEFAULT_SSTAR)
    block = cfg.get("sstar") if isinstance(cfg, dict) else None
    for src in (block, sstar):
        if isinstance(src, bool):
            out["fine_grid_enabled"] = src
        elif isinstance(src, dict):
            for k in ("fine_grid_enabled", "interpolate"):
                if src.get(k) is not None:
                    out[k] = bool(src[k])
            for k in ("grid_min", "grid_max", "grid_step"):
                if src.get(k) is not None:
                    out[k] = float(src[k])
    return out


def sstar_grid_enabled(cfg: Optional[Dict[str, Any]] = None, sstar: Optional[Any] = None) -> bool:
    return bool(sstar_cfg(cfg, sstar)["fine_grid_enabled"])


def s_grid(grid_min: float = GRID_MIN, grid_max: float = GRID_MAX, step: float = GRID_STEP) -> np.ndarray:
    """`np.arange(grid_min, grid_max + 1e-12, step)` rounded to 3 decimals (the serialized resolution).

    With the defaults: 4801 points, 6.000, 6.005, ..., 30.000 (30.0 included).
    """
    return np.round(np.arange(grid_min, grid_max + 1e-12, step), DECIMALS)


def _clean(S_vetted) -> np.ndarray:
    S = np.asarray(S_vetted, float).ravel()
    return S[np.isfinite(S)]


def n_survivors(S_vetted, s) -> np.ndarray:
    """n_null(S >= s) for each s (vectorised, exact comparison against the rounded grid values)."""
    S = np.sort(_clean(S_vetted))
    return len(S) - np.searchsorted(S, np.atleast_1d(np.asarray(s, float)), side="left")


def sstar_discrete(S_vetted, budget: float, grid: Optional[np.ndarray] = None) -> float:
    """Discrete S*: the smallest grid point `s` with `n_null(S >= s) <= floor(B)`.

    `S_vetted` are the vetted null S values (non-finite values are ignored); `B` is the survivor budget
    (`floor(B + 1e-9)` is used so that float noise in B = p * n_null cannot change the integer budget).
    Survivor counts are non-increasing in `s`, so the answer is the first grid point meeting the budget.
    If even the last grid point (30.0) has more than floor(B) survivors there is no valid S* on the grid and
    `SStarGridError` is raised (no silent clamp to 30.0, which would exceed the budget).  Raises as well for B < 0.
    """
    if not np.isfinite(budget) or budget < 0:
        raise SStarGridError(f"invalid survivor budget {budget!r}")
    g = s_grid() if grid is None else np.asarray(grid, float)
    cap = int(math.floor(budget + _BUDGET_EPS))
    ok = n_survivors(S_vetted, g) <= cap
    if not ok.any():
        raise SStarGridError(f"no grid point <= {g[-1]} has n_null(S >= s) <= {cap}")
    return float(g[int(np.argmax(ok))])


def pava_nonincreasing(y) -> np.ndarray:
    """Pool-adjacent-violators: least-squares non-increasing fit of `y` (equal weights)."""
    y = np.asarray(y, float)
    blocks = []                                   # [mean, count]
    for v in y:
        blocks.append([float(v), 1])
        while len(blocks) > 1 and blocks[-2][0] < blocks[-1][0]:      # increase => violates non-increasing
            m2, c2 = blocks.pop()
            m1, c1 = blocks.pop()
            blocks.append([(m1 * c1 + m2 * c2) / (c1 + c2), c1 + c2])
    return np.concatenate([np.full(c, m) for m, c in blocks]) if blocks else y


def sstar_interpolated(S_vetted, budget: float, grid: Optional[np.ndarray] = None) -> Optional[float]:
    """Monotone log-survival interpolation (optional extra; NOT the default).

    Knots: for the distinct vetted S values x (ascending) the survivor count N(x) = n(S >= x).  log N is made
    non-increasing with PAVA, then linearly interpolated in S.  The crossing of log N with log B (continuous
    B, B >= 1) is rounded UP to the next grid point.  Returns None when B < 1, B is not below the first knot's
    count, or there are fewer than 2 distinct values (then only the discrete rule applies).  Callers must use
    `sstar_grid` which applies the guard `max(interp, discrete)`.
    """
    S = np.sort(_clean(S_vetted))
    if budget < 1 or len(np.unique(S)) < 2:
        return None
    x, first = np.unique(S, return_index=True)
    logN = pava_nonincreasing(np.log(len(S) - first))
    if not (logN[-1] <= math.log(budget) <= logN[0]):
        return None
    # logN is non-increasing in x; np.interp needs ascending xp => reverse
    s = float(np.interp(math.log(budget), logN[::-1], x[::-1]))
    g = s_grid() if grid is None else np.asarray(grid, float)
    i = int(np.searchsorted(g, s - 1e-12, side="left"))
    return float(g[min(i, len(g) - 1)])


def sstar_grid(S_vetted, budget: float, interpolate: bool = False, grid: Optional[np.ndarray] = None) -> Dict[str, Any]:
    """Fine-grid S*: discrete rule, plus (only if `interpolate`) the guarded interpolation.

    Returns dict(s_star, s_star_discrete, s_star_interp, rule, budget, n_budget, n_survivors, n_vetted).
    With `interpolate`, S* = max(S*_interp, S*_discrete) so it can never admit more empirical survivors than the
    discrete rule.  Raises `SStarGridError` when the discrete rule has no solution on the grid.
    """
    g = s_grid() if grid is None else np.asarray(grid, float)
    disc = sstar_discrete(S_vetted, budget, g)
    interp = sstar_interpolated(S_vetted, budget, g) if interpolate else None
    s = disc if interp is None else max(interp, disc)
    return dict(s_star=round(float(s), DECIMALS), s_star_discrete=round(disc, DECIMALS),
                s_star_interp=None if interp is None else round(float(interp), DECIMALS),
                rule="grid_discrete_interp_guarded" if interpolate else "grid_discrete",
                budget=float(budget), n_budget=int(math.floor(budget + _BUDGET_EPS)),
                n_survivors=int(n_survivors(S_vetted, s)[0]), n_vetted=int(len(_clean(S_vetted))))


def select_sstar(S_vetted, n_null, p, u, cfg: Optional[Dict[str, Any]] = None, sstar: Optional[Any] = None):
    """S* entry point honouring the setting.

    Off (default; keys missing or false): exactly `interim.gpd_threshold(S_vetted, n_null, p, u)` (v1.1 GPD).
    On: `sstar_grid(S_vetted, p * n_null, interpolate=cfg interpolate)` on the configured grid.
    """
    c = sstar_cfg(cfg, sstar)
    if not c["fine_grid_enabled"]:
        from . import interim
        return interim.gpd_threshold(S_vetted, n_null, p, u)
    return sstar_grid(S_vetted, p * n_null, interpolate=c["interpolate"],
                      grid=s_grid(c["grid_min"], c["grid_max"], c["grid_step"]))
