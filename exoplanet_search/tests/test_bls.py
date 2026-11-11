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


def test_s_statistic_white_noise_and_transit(synthetic_lc):
    from exosearch import lightcurve
    t, f, e, truth = synthetic_lc
    tf, ff, ef, *_ = lightcurve.detrend(t, f, e)
    tb, fb, eb = lightcurve.bin_lc(tf, ff, ef)
    sw, cad, prof = bls.noise_profile(tf, ff)
    betas = [b for _, b in prof.values()]
    assert all(b > 0.6 for b in betas) and betas[0] < 1.4   # ~white at 1 h; sinusoid residual adds red noise at long D
    out = bls.search_s(tb, fb, tf, ff)
    assert abs(out["period"] / truth["period"] - 1) < 0.01 and out["S"] > 15 and out["n_tr"] >= 5


def test_s_statistic_noise_only():
    from exosearch import lightcurve
    rng = np.random.default_rng(3)
    t = np.arange(0, 27, 2 / 1440)
    f = 1 + rng.normal(0, 1e-3, t.size)
    e = np.full(t.size, 1e-3)
    tf, ff, ef, *_ = lightcurve.detrend(t, f, e)
    tb, fb, eb = lightcurve.bin_lc(tf, ff, ef)
    assert bls.search_s(tb, fb, tf, ff)["S"] < 8
