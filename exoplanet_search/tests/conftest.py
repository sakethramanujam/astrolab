import os
import sys
from pathlib import Path

# Thread caps before any test imports numba or BLAS
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMBA_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

# Configure fallback EXO_DATA_ROOT or EXO_OFFLINE if /media default does not exist
if "EXO_DATA_ROOT" not in os.environ and not Path("/media/nfs01/astro_data/tess").exists():
    _candidates = [
        Path(__file__).resolve().parents[1] / "data" / "tess",
        Path(__file__).resolve().parents[2] / "data" / "tess",
        Path(__file__).resolve().parents[3] / "data" / "tess",
    ]
    for _cand in _candidates:
        if _cand.exists():
            os.environ["EXO_DATA_ROOT"] = str(_cand)
            break
    else:
        os.environ.setdefault("EXO_OFFLINE", "1")

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
