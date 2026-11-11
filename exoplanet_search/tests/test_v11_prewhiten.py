import json
import numpy as np
import pytest

from exosearch.v11 import prewhiten, search, settings


def _cfg():
    return settings.cfg()


def test_null_light_curve_not_prewhitened():
    """Stars with pure noise (like flipped nulls) should not have terms removed."""
    t = np.arange(3000, 3027, 2 / 1440)
    for seed in range(5):
        rng = np.random.default_rng(seed)
        f = 1 + rng.normal(0, 1e-3, t.size)
        r, log = prewhiten.prewhiten(t, f, _cfg()["prewhiten"])
        assert log["prewhiten_n"] == 0
        assert np.array_equal(r, f)


def test_usp_depth_preserved_by_comb_guard():
    """USP transits with harmonics above 5 c/d should be protected without depth suppression."""
    t = np.arange(3000, 3027, 2 / 1440)
    rng = np.random.default_rng(42)
    depth = 2e-3
    P = 0.7
    dur = 1.2 / 24
    in_tr = search.in_transit(t, P, 3000.1, dur)

    f = np.ones_like(t)
    f[in_tr] -= depth
    f += rng.normal(0, 5e-4, t.size)

    r_pw, log_pw = prewhiten.prewhiten(t, f, _cfg()["prewhiten"], comb_guard=True)
    # Comb guard should recognize transit harmonics and protect them
    assert log_pw["prewhiten_comb"] >= 1
    assert log_pw["prewhiten_n"] == 0

    # Depth should be preserved (median in-transit flux difference within 2%)
    d_raw = 1.0 - np.mean(f[in_tr])
    d_pw = 1.0 - np.mean(r_pw[in_tr])
    assert np.isclose(d_pw / d_raw, 1.0, atol=0.02)


def test_transit_with_pulsation_removes_pulsation_only():
    """When a real delta Scuti pulsation and a transit coexist, prewhitening cleans the pulsation only."""
    t = np.arange(3000, 3027, 2 / 1440)
    rng = np.random.default_rng(10)
    # Transit (P = 0.8 d, T14 = 1.5 h, depth = 1.5e-3)
    in_tr = search.in_transit(t, 0.8, 3000.1, 1.5 / 24)
    f = np.ones_like(t)
    f[in_tr] -= 1.5e-3
    # Pulsation (f = 22.97 c/d, amp = 3e-3)
    f += 3e-3 * np.sin(2 * np.pi * 22.97 * t)
    f += rng.normal(0, 5e-4, t.size)

    r_pw, log_pw = prewhiten.prewhiten(t, f, _cfg()["prewhiten"], comb_guard=True)
    terms = json.loads(log_pw["prewhiten_terms"])
    combs = json.loads(log_pw["prewhiten_comb_freqs"])

    # Pulsation should be removed
    assert log_pw["prewhiten_n"] >= 1
    assert any(abs(term[0] - 22.97) < 0.05 for term in terms)
    # Transit comb should be protected
    assert log_pw["prewhiten_comb"] >= 1

    # In-transit depth should still be intact after removing the pulsation
    d_rec = 1.0 - np.mean(r_pw[in_tr])
    assert 1.3e-3 < d_rec < 1.7e-3


def test_two_halves_coherence_rejects_transient_glitch():
    """A transient signal present in only one orbit should not be removed as a pulsation."""
    t = np.arange(3000, 3027, 2 / 1440)
    rng = np.random.default_rng(7)
    t_mid = 0.5 * (t.min() + t.max())
    f = 1 + rng.normal(0, 5e-4, t.size)

    # Add a high-frequency sine burst only in the first half
    m1 = t < t_mid
    f[m1] += 2e-3 * np.sin(2 * np.pi * 18.5 * t[m1])

    r, log = prewhiten.prewhiten(t, f, _cfg()["prewhiten"])
    terms = json.loads(log["prewhiten_terms"])
    # Should reject the transient burst because it is not coherent in the second half
    assert not any(abs(term[0] - 18.5) < 0.05 for term in terms)


def test_transit_comb_matches_harmonics():
    """transit_comb_matches detects transit combs across multiple periods and harmonics."""
    t = np.arange(3000, 3027, 2 / 1440)
    fnyq = 0.5 / (2 / 1440)

    # Simulate transit at P = 0.5 d (f0 = 2 c/d, harmonic at 8 c/d)
    in_tr = search.in_transit(t, 0.5, 3000.1, 1.0 / 24)
    r = np.zeros_like(t)
    r[in_tr] -= 2e-3
    noise = 5e-4 * 2 / np.sqrt(t.size)

    a1, _ = prewhiten._fit(t, r, 8.0)
    ms = prewhiten.transit_comb_matches(t, r, 8.0, a1, noise, fnyq)
    assert len(ms) >= 1
    assert any(abs(m / 8.0 - 0.5) < 0.02 for m in ms)

    # Pure pulsation should not match
    r_puls = 3e-3 * np.sin(2 * np.pi * 22.97 * t)
    a1_p, _ = prewhiten._fit(t, r_puls, 22.97)
    ms_p = prewhiten.transit_comb_matches(t, r_puls, 22.97, a1_p, noise, fnyq)
    assert len(ms_p) == 0
