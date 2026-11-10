"""numpy reference of the vectorized BLS used by the cuda/mlx stubs."""
import importlib.util
import time

import numpy as np
import pytest

from exosearch import bls, lightcurve


def _prep(lc):
    t, f, e, truth = lc
    tf, ff, ef, *_ = lightcurve.detrend(t, f, e)
    return (*lightcurve.bin_lc(tf, ff, ef, 10), truth)


def test_numpy_vector_bls_matches_astropy(synthetic_lc):
    t, f, e, truth = _prep(synthetic_lc)
    t1 = time.perf_counter()
    vec = bls.search(t, f, e, backend="numpy")
    t2 = time.perf_counter()
    cpu = bls.search(t, f, e, backend="cpu")
    print(f"numpy-vector {t2-t1:.2f}s  astropy {time.perf_counter()-t2:.2f}s")
    assert abs(vec["period"] / truth["period"] - 1) < 0.005
    assert abs(vec["period"] / cpu["period"] - 1) < 0.005
    ph = ((vec["t0"] - truth["t0"]) / truth["period"]) % 1
    assert min(ph, 1 - ph) < 0.02
    assert vec["sde"] > 9
    r = np.corrcoef(vec["power"], cpu["power"])[0, 1]
    assert r > 0.8, r


@pytest.mark.skipif(importlib.util.find_spec("cupy") is None, reason="no CuPy / CUDA on this host")
def test_cuda_matches_numpy(synthetic_lc):
    t, f, e, _ = _prep(synthetic_lc)
    assert abs(bls.search(t, f, e, backend="cuda")["period"] / bls.search(t, f, e, backend="numpy")["period"] - 1) < 1e-3


@pytest.mark.skipif(importlib.util.find_spec("mlx") is None, reason="no MLX (Apple Silicon only)")
def test_mlx_matches_numpy(synthetic_lc):
    t, f, e, _ = _prep(synthetic_lc)
    assert abs(bls.search(t, f, e, backend="mlx")["period"] / bls.search(t, f, e, backend="numpy")["period"] - 1) < 1e-3
