import numpy as np

from exosearch import bls, lightcurve


def _prep(lc):
    t, f, e, truth = lc
    tf, ff, ef, *_ = lightcurve.detrend(t, f, e)
    return (*lightcurve.bin_lc(tf, ff, ef, 10), truth)


def test_cpu_bls_recovers_injection(synthetic_lc):
    t, f, e, truth = _prep(synthetic_lc)
    out = bls.search(t, f, e, backend="cpu")
    assert abs(out["period"] / truth["period"] - 1) < 0.005
    assert out["sde"] > 9
    st = bls.stats(t, f, e, out["period"], out["duration"], out["t0"])
    assert abs(st["depth"][0] / truth["depth"] - 1) < 0.35


def test_period_grid_needs_two_transits():
    t = np.linspace(0, 27, 1000)
    g = bls.period_grid(t)
    assert g.min() >= 0.5 - 1e-9 and g.max() <= 13.5 + 1e-9 and np.all(np.diff(g) > 0)


def test_noise_has_low_sde():
    rng = np.random.default_rng(3)
    t = np.arange(0, 27, 10 / 1440)
    f = 1 + rng.normal(0, 1e-3, t.size)
    out = bls.search(t, f, np.full_like(t, 1e-3), backend="cpu")
    assert out["sde"] < 12
