"""Synthetic pure-pulsator light curves (no transit) for measuring the T0 comb-guard false-comb rate.

Realistic multi-mode delta Scuti / gamma Doradus spectra on an S70-like cadence (2 min, a mid-sector gap and two
end gaps), white plus red (AR(1)) noise. Used by tests and by the measurement helper `false_comb_rate`.
"""
from __future__ import annotations

import numpy as np


def s70_times(cadence_min: float = 2.0) -> np.ndarray:
    t = np.arange(3000.0, 3027.0, cadence_min / 1440.0)
    keep = ~(((t > 3013.0) & (t < 3014.3)) | (t > 3026.0) | (t < 3000.3))
    return t[keep]


def pulsator(seed: int, kind: str | None = None, t: np.ndarray | None = None):
    """-> (t, flux, truth) with truth = list of (f, amp). kind in {'dsct', 'gdor'} (random if None)."""
    rng = np.random.default_rng(seed)
    t = s70_times() if t is None else t
    kind = kind or ("dsct" if rng.random() < 0.5 else "gdor")
    n_modes = int(rng.integers(4, 11))
    modes = []
    if kind == "dsct":
        fs = rng.uniform(6.0, 60.0, n_modes)
        amps = 10 ** rng.uniform(-3.7, -2.3, n_modes)
    else:
        fs = rng.uniform(0.3, 3.5, n_modes)
        amps = 10 ** rng.uniform(-3.3, -2.3, n_modes)
    for f, a in zip(fs, amps):
        modes.append((float(f), float(a)))
    # combination / harmonic terms of the two strongest modes
    modes.sort(key=lambda x: -x[1])
    if n_modes >= 2 and kind == "dsct":
        modes.append((modes[0][0] + modes[1][0], 0.15 * modes[1][1]))
        modes.append((2 * modes[0][0], 0.2 * modes[0][1]))
    flux = np.ones_like(t)
    for f, a in modes:
        flux += a * np.sin(2 * np.pi * f * t + rng.uniform(0, 2 * np.pi))
    white = rng.normal(0, 3e-4, t.size)
    red = np.zeros(t.size)
    e = rng.normal(0, 2e-4, t.size)
    phi = np.exp(-(2.0 / 1440.0) / 0.3)
    for i in range(1, t.size):
        red[i] = phi * red[i - 1] + e[i]
    return t, flux + white + red, modes


def false_comb_rate(cfg: dict, seeds=range(50)) -> tuple[float, int]:
    """Fraction of pure-pulsator light curves where the comb guard fires at least once (every firing is false)."""
    from . import prewhiten
    pw_cfg = cfg.get("prewhiten", cfg)
    n = 0
    seeds_list = list(seeds)
    for s in seeds_list:
        t, f, _ = pulsator(s)
        _, log = prewhiten.prewhiten(t, f, pw_cfg)
        n += int(log["prewhiten_comb"] > 0)
    return n / len(seeds_list), n
