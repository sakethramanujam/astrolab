import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))


@pytest.fixture
def synthetic_lc():
    """27-day TESS-like 2-min light curve with variability, a mid-sector gap and a transit."""
    rng = np.random.default_rng(1)
    t = np.arange(0, 27.0, 2 / 1440.0) + 3000.0
    t = t[(t < 3013.0) | (t > 3014.0)]
    f = 1 + 0.004 * np.sin(2 * np.pi * (t - 3000) / 4.3) + rng.normal(0, 1e-3, t.size)
    P, t0, dur, depth = 3.71, 3001.3, 0.11, 2500e-6
    ph = ((t - t0 + 0.5 * P) % P) - 0.5 * P
    f[np.abs(ph) < dur / 2] -= depth
    e = np.full_like(t, 1e-3)
    return t, f, e, dict(period=P, t0=t0, duration=dur, depth=depth)
