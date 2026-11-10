import numpy as np

from exosearch import lightcurve


def test_detrend_removes_variability_keeps_transit(synthetic_lc):
    t, f, e, truth = synthetic_lc
    tf, ff, ef, trend, m = lightcurve.detrend(t, f, e, window=0.75)
    assert np.std(ff) < 1.5e-3                       # 4000 ppm sine removed
    ph = ((tf - truth["t0"] + 0.5 * truth["period"]) % truth["period"]) - 0.5 * truth["period"]
    intr = np.abs(ph) < truth["duration"] / 3
    depth = np.median(ff[~intr]) - np.median(ff[intr])
    assert 0.6 * truth["depth"] < depth < 1.4 * truth["depth"]
    tb, fb, eb = lightcurve.bin_lc(tf, ff, ef, minutes=10)
    assert len(tb) < len(tf) / 4 and np.all(eb > 0)
