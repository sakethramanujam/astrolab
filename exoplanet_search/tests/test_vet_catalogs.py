import os

import numpy as np
import pytest
from astropy.timeseries import BoxLeastSquares

from exosearch import catalogs, lightcurve, store, vet


def _stats(t, f, e, truth, period=None):
    P = period or truth["period"]
    return BoxLeastSquares(t, f, dy=e).compute_stats(P, truth["duration"], truth["t0"])


def test_vet_passes_clean_transit(synthetic_lc):
    t, f, e, truth = synthetic_lc
    tf, ff, ef, *_ = lightcurve.detrend(t, f, e)
    st = _stats(tf, ff, ef, truth)
    out = dict(period=truth["period"], duration=truth["duration"], sde=15.0)
    v = vet.vet(st, out, {"radius": 1.0}, tf, ff)
    assert v["passed"] == 1, v["flags"]
    assert v["snr"] > 20 and v["odd_even_sigma"] < 3 and v["n_transits"] >= 6


def test_vet_flags_eclipsing_binary(synthetic_lc):
    t, f, e, truth = synthetic_lc
    f = f.copy()
    P, t0, dur = truth["period"], truth["t0"], truth["duration"]
    ep = np.round((t - t0) / P).astype(int)
    ph = ((t - t0 + 0.5 * P) % P) - 0.5 * P
    f[(np.abs(ph) < dur / 2) & (ep % 2 == 1)] -= 0.004          # unequal odd/even
    ph2 = ((t - t0) % P) - 0.5 * P
    f[np.abs(ph2) < dur / 2] -= 0.0015                           # secondary eclipse
    tf, ff, ef, *_ = lightcurve.detrend(t, f, e)
    v = vet.vet(_stats(tf, ff, ef, truth), dict(period=P, duration=dur, sde=15.0), {"radius": 1.0}, tf, ff)
    assert "odd_even" in v["flags"] and "secondary" in v["flags"] and v["passed"] == 0


def test_duration_sanity():
    assert 0.1 < vet.expected_duration(3.0, 1.0) < 0.15     # ~2.9 h for a 3 d orbit round the Sun


def test_store_upsert_idempotent(tmp_path):
    con = store.connect(tmp_path / "r.sqlite")
    store.upsert(con, dict(tic=1, sector=70, status="ok", sde=9.0))
    store.upsert(con, dict(tic=1, sector=70, status="ok", sde=10.0))
    df = store.frame(con, 70)
    assert len(df) == 1 and df.sde[0] == 10.0 and store.done_tics(con, 70) == {1}


@pytest.mark.skipif(os.environ.get("EXO_OFFLINE") == "1", reason="offline")
def test_crossmatch_known_planet():
    cats = catalogs.load_catalogs()
    assert len(cats["toi"]) > 5000
    toi = cats["toi"].dropna(subset=["period"]).iloc[0]
    m = catalogs.crossmatch(int(toi.tic), float(toi.period) * 2, cats)      # harmonic still matches
    assert m["known_toi"] and m["toi_period_match"] == 1
    if cats["ps_ok"]:
        row = cats["ps"].dropna(subset=["pl_orbper"]).iloc[0]
        m = catalogs.crossmatch(int(row.tic), float(row.pl_orbper), cats)
        assert m["known_planet"] and m["planet_period_match"] == 1
    assert catalogs.crossmatch(1, 3.0, cats)["known_toi"] is None


def test_common_period_flag():
    import pandas as pd
    rng = np.random.default_rng(0)
    P = np.r_[rng.uniform(0.5, 13, 200), np.full(15, 6.40) * (1 + rng.normal(0, 1e-3, 15))]
    df = pd.DataFrame(dict(status="ok", period=P, flags="", passed=1))
    out = vet.flag_common_periods(df)
    assert out.tail(15)["flags"].str.contains("common_period").all() and (out.tail(15).passed == 0).all()
    assert out.head(200)["flags"].str.contains("common_period").sum() < 10
