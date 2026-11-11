import json
import numpy as np
import pytest
from scipy.signal import lfilter

from exosearch.v11 import prewhiten, search, settings, validate


def _cfg():
    return settings.cfg()


def _amp_after(t, r, f):
    return max(prewhiten._fit(t, r, f + d)[0] for d in np.linspace(-0.01, 0.01, 5))


def _depth_ratio(r, clean, in_tr):
    d_clean = float(np.mean(clean[~in_tr]) - np.mean(clean[in_tr]))
    d_r = float(np.mean(r[~in_tr]) - np.mean(r[in_tr]))
    return d_r / d_clean if d_clean > 0 else 1.0


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

    clean = np.ones_like(t)
    clean[in_tr] -= depth
    f = clean + rng.normal(0, 5e-4, t.size)

    r_pw, log_pw = prewhiten.prewhiten(t, f, _cfg()["prewhiten"], comb_guard=True)
    # Comb guard should recognize transit harmonics and protect them
    assert log_pw["prewhiten_comb"] >= 1
    assert log_pw["prewhiten_n"] == 0

    # Depth should be preserved (median in-transit flux difference within 2%)
    assert np.isclose(_depth_ratio(r_pw, clean, in_tr), 1.0, atol=0.02)


def test_transit_with_pulsation_removes_pulsation_only():
    """When a real delta Scuti pulsation and a transit coexist, prewhitening cleans the pulsation only."""
    t = np.arange(3000, 3027, 2 / 1440)
    rng = np.random.default_rng(10)
    # Transit (P = 0.8 d, T14 = 1.5 h, depth = 1.5e-3)
    in_tr = search.in_transit(t, 0.8, 3000.1, 1.5 / 24)
    clean = np.ones_like(t)
    clean[in_tr] -= 1.5e-3
    # Pulsation (f = 22.97 c/d, amp = 3e-3)
    f = clean + 3e-3 * np.sin(2 * np.pi * 22.97 * t) + rng.normal(0, 5e-4, t.size)

    r_pw, log_pw = prewhiten.prewhiten(t, f, _cfg()["prewhiten"], comb_guard=True)
    terms = json.loads(log_pw["prewhiten_terms"])

    # Pulsation should be removed
    assert log_pw["prewhiten_n"] >= 1
    assert any(abs(term[0] - 22.97) < 0.05 for term in terms)
    assert _amp_after(t, r_pw, 22.97) < 0.2 * 3e-3
    # Transit comb should be protected
    assert log_pw["prewhiten_comb"] >= 1

    # In-transit depth should still be intact after removing the pulsation
    assert abs(_depth_ratio(r_pw, clean, in_tr) - 1.0) < 0.05


def test_delta_scuti_with_rotation_and_red_noise_is_removed():
    """Finding 1 regression: a peak at 5.5, 8.0, 12.0 c/d with rotation and red noise must still be removed."""
    t = np.arange(3000.0, 3027.0, 2 / 1440)
    for seed, fm in enumerate([5.5, 8.0, 12.0]):
        rng = np.random.default_rng(90 + seed)
        sigma = 4e-4
        f = 1.0 + rng.normal(0, sigma, t.size)
        f += lfilter([1.0], [1.0, -0.999], rng.normal(0, sigma * 0.03, t.size))
        fr = rng.uniform(0.3, 1.5)
        f += 2e-3 * np.sin(2 * np.pi * fr * t + rng.uniform(0, 6.28)) + 0.8e-3 * np.sin(4 * np.pi * fr * t)
        f += 3e-3 * np.sin(2 * np.pi * fm * t + rng.uniform(0, 6.28))

        r, log = prewhiten.prewhiten(t, f, _cfg()["prewhiten"])
        assert log["prewhiten_n"] >= 1
        assert _amp_after(t, r, fm) < 0.2 * 3e-3


def test_hads_harmonics_removed():
    """HADS stars with non-linear harmonics f, 2f, 3f should not be falsely identified as transit combs."""
    t = np.arange(3000, 3027, 2 / 1440)
    rng = np.random.default_rng(123)
    f = 1 + rng.normal(0, 4e-4, t.size)
    f += 5e-3 * np.sin(2 * np.pi * 6.1 * t)
    f += 2e-3 * np.sin(2 * np.pi * 12.2 * t + 1.2)
    f += 0.8e-3 * np.sin(2 * np.pi * 18.3 * t + 2.5)

    r, log = prewhiten.prewhiten(t, f, _cfg()["prewhiten"])
    assert log["prewhiten_n"] >= 3
    for fm in (6.1, 12.2, 18.3):
        assert _amp_after(t, r, fm) < 0.25 * 5e-3


def test_comb_exclusion_coverage_bounded():
    """Finding 2: comb exclusion windows must never blank the spectrum; coverage must stay <= 20%."""
    t = np.arange(3000, 3027, 2 / 1440)
    rng = np.random.default_rng(3)
    clean = np.ones_like(t)
    in_tr = search.in_transit(t, 0.5, 3000.1, 1.0 / 24)
    clean[in_tr] -= 5e-3
    f = clean + rng.normal(0, 3e-4, t.size)
    # Also add a pulsation at 17.3 c/d
    f += 3e-3 * np.sin(2 * np.pi * 17.3 * t)

    r, log = prewhiten.prewhiten(t, f, _cfg()["prewhiten"])
    # Comb should be detected for the transit and the pulsation at 17.3 should still be removed
    assert log["prewhiten_comb"] >= 1
    assert log["prewhiten_n"] >= 1
    assert _amp_after(t, r, 17.3) < 0.2 * 3e-3


def test_short_duration_transit_at_low_snr():
    """Short-duration (T14 < 1 h) transit at lower depth is protected by comb guard."""
    t = np.arange(3000, 3027, 2 / 1440)
    rng = np.random.default_rng(55)
    in_tr = search.in_transit(t, 0.45, 3000.1, 0.6 / 24)
    clean = np.ones_like(t)
    clean[in_tr] -= 1.8e-3
    f = clean + rng.normal(0, 4e-4, t.size)

    r, log = prewhiten.prewhiten(t, f, _cfg()["prewhiten"])
    assert log["prewhiten_comb"] >= 1
    assert log["prewhiten_n"] == 0
    assert abs(_depth_ratio(r, clean, in_tr) - 1.0) < 0.05


def test_two_halves_coherence_rejects_transient_glitch():
    """A transient signal present in only one orbit should not be removed as a pulsation."""
    t = np.arange(3000, 3027, 2 / 1440)
    rng = np.random.default_rng(7)
    t_mid = 0.5 * (t.min() + t.max())
    f = 1 + rng.normal(0, 5e-4, t.size)

    m1 = t < t_mid
    f[m1] += 2e-3 * np.sin(2 * np.pi * 18.5 * t[m1])

    r, log = prewhiten.prewhiten(t, f, _cfg()["prewhiten"])
    terms = json.loads(log["prewhiten_terms"])
    assert not any(abs(term[0] - 18.5) < 0.05 for term in terms)


def test_two_halves_coherence_keeps_beating_pair():
    """Two close modes that beat must not be rejected by the coherence check."""
    t = np.arange(3000, 3027, 2 / 1440)
    rng = np.random.default_rng(8)
    f = 1 + rng.normal(0, 4e-4, t.size)
    for fm in (12.00, 12.06):
        f = f + 3e-3 * np.sin(2 * np.pi * fm * t)

    r, log = prewhiten.prewhiten(t, f, _cfg()["prewhiten"])
    assert log["prewhiten_n"] >= 2
    for fm in (12.00, 12.06):
        assert _amp_after(t, r, fm) < 0.3 * 3e-3


def test_pulsator_efficacy_check_passes():
    """Verify validate.pulsator_efficacy passes on the current policy."""
    res = validate.pulsator_efficacy(_cfg(), n=8)
    assert res["pass_"]
    assert res["median_removed"] >= 0.8
    assert res["min_removed"] >= 0.5


def test_policy_removing_nothing_fails_pulsator_check(monkeypatch):
    """Verify that a dummy policy removing nothing fails the validation check."""
    def dummy_prewhiten(t, f, cfg, comb_guard=None):
        return np.asarray(f, float).copy(), dict(prewhiten_n=0, prewhiten_comb=0,
                                                 prewhiten_terms="[]", prewhiten_comb_freqs="[]")
    monkeypatch.setattr(prewhiten, "prewhiten", dummy_prewhiten)
    res = validate.pulsator_efficacy(_cfg(), n=4)
    assert not res["pass_"]


def test_deterministic():
    """Two consecutive prewhitening runs on identical inputs produce bitwise equal outputs."""
    t = np.arange(3000, 3027, 2 / 1440)
    rng = np.random.default_rng(11)
    f = 1 + rng.normal(0, 4e-4, t.size)
    in_tr = search.in_transit(t, 0.9, 3000.1, 1.5 / 24)
    f[in_tr] -= 3e-3
    f += 3e-3 * np.sin(2 * np.pi * 14.2 * t)

    a, log_a = prewhiten.prewhiten(t, f, _cfg()["prewhiten"])
    b, log_b = prewhiten.prewhiten(t, f, _cfg()["prewhiten"])
    assert np.array_equal(a, b)
    assert log_a == log_b
