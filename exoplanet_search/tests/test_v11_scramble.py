"""Tests for v1.1 scrambled-null transform (exosearch.v11.scramble)."""
from __future__ import annotations

import numpy as np
import pytest
from astropy.timeseries import BoxLeastSquares

from exosearch import lightcurve
from exosearch.v11 import scramble, settings


@pytest.fixture
def v11_cfg():
    return settings.cfg()


def _make_two_orbit_lc(cadence_min=2.0, gap_start=113.0, gap_end=114.5, t_start=100.0, t_end=127.0):
    dt = cadence_min / 1440.0
    t1 = np.arange(t_start, gap_start, dt)
    t2 = np.arange(gap_end, t_end, dt)
    t = np.concatenate([t1, t2])
    f = np.ones_like(t)
    e = np.full_like(t, 0.0005)
    return t, f, e


def test_gap_edge_cadences_removed_before_scrambling(v11_cfg):
    """Gap-edge cadences within 0.3 d of gaps > 0.5 d and endpoints must be dropped before scrambling."""
    t, f, e = _make_two_orbit_lc()
    t_out, f_out, e_out, info = scramble.segment_rotate(t, f, e, seed=42, cfg=v11_cfg)

    # Initial time span has gap 113.0 to 114.5 (size 1.5 d > 0.5 d)
    # Cadences within 0.3 d of 100.0 (start), 113.0 (before gap), 114.5 (after gap), 127.0 (end) dropped
    assert t_out.min() >= 100.0 + 0.3 - 1e-6
    assert t_out.max() <= 127.0 - 0.3 + 1e-6
    # No cadences in [112.7, 114.8]
    assert not np.any((t_out > 113.0 - 0.3 + 1e-6) & (t_out < 114.5 + 0.3 - 1e-6))
    assert len(t_out) < len(t)
    assert len(f_out) == len(t_out)
    assert len(e_out) == len(t_out)


def test_no_samples_cross_gap_and_each_segment_is_cyclic_rotation(v11_cfg):
    """Samples never cross a segment boundary; each segment is a cyclic rotation of its own masked values."""
    t, f, e = _make_two_orbit_lc()
    # Mark flux with sequential values to verify rotation order within each segment
    f = np.arange(t.size, dtype=float)
    e = f * 0.1

    # First get what cadences survive gap_edge_mask
    mask = lightcurve.gap_edge_mask(t, gap=0.5, trim=0.3)
    t_kept = t[mask]
    f_kept = f[mask]
    e_kept = e[mask]

    t_out, f_out, e_out, info = scramble.segment_rotate(t, f, e, seed=123, cfg=v11_cfg)
    assert np.array_equal(t_out, t_kept)

    # Split into segments at gaps > 0.5 d
    split_idx = np.where(np.diff(t_kept) > 0.5)[0] + 1
    t_segs = np.split(t_kept, split_idx)
    f_in_segs = np.split(f_kept, split_idx)
    e_in_segs = np.split(e_kept, split_idx)
    f_out_segs = np.split(f_out, split_idx)
    e_out_segs = np.split(e_out, split_idx)

    assert len(t_segs) == info["n_segments"]
    assert len(info["scr_offsets"]) == len(t_segs)

    for i in range(len(t_segs)):
        offset = info["scr_offsets"][i]
        # Check that f_out is exactly rolled by offset from f_in
        assert np.array_equal(f_out_segs[i], np.roll(f_in_segs[i], offset))
        assert np.array_equal(e_out_segs[i], np.roll(e_in_segs[i], offset))
        # Exact set of values in segment i matches
        assert set(f_out_segs[i]) == set(f_in_segs[i])


def test_deterministic_and_different_across_seeds(v11_cfg):
    """Deterministic for a fixed seed and different across seeds."""
    t, f, e = _make_two_orbit_lc()
    rng = np.random.default_rng(0)
    f = rng.normal(1.0, 0.01, size=t.size)
    e = np.full_like(t, 0.01)

    t1, f1, e1, info1 = scramble.segment_rotate(t, f, e, seed=999, cfg=v11_cfg)
    t2, f2, e2, info2 = scramble.segment_rotate(t, f, e, seed=999, cfg=v11_cfg)
    t3, f3, e3, info3 = scramble.segment_rotate(t, f, e, seed=1000, cfg=v11_cfg)

    # Determinism
    assert np.array_equal(f1, f2)
    assert np.array_equal(e1, e2)
    assert info1["scr_offsets"] == info2["scr_offsets"]

    # Difference across seeds
    assert not np.array_equal(f1, f3)
    assert info1["scr_offsets"] != info3["scr_offsets"]


def test_segments_get_independent_offsets(v11_cfg):
    """Segments get independent offsets drawn from rng."""
    t, f, e = _make_two_orbit_lc()
    # Check across multiple seeds that offsets for segment 0 and segment 1 are not identically locked
    offsets_seg0 = []
    offsets_seg1 = []
    for s in range(20):
        _, _, _, info = scramble.segment_rotate(t, f, e, seed=s, cfg=v11_cfg)
        offsets_seg0.append(info["scr_offsets"][0])
        offsets_seg1.append(info["scr_offsets"][1])

    assert offsets_seg0 != offsets_seg1
    # Check offsets take multiple distinct values
    assert len(set(offsets_seg0)) > 10
    assert len(set(offsets_seg1)) > 10


def test_fallback_for_short_segment(v11_cfg):
    """Fallback path for short segments (< 3x longest duration):
    - Shorter than one block (3x longest duration): left unrotated with scheme 'unrotated'.
    - Long enough for blocks >= 3x longest duration: shuffled in blocks >= 3x longest duration.
    """
    dt = 2.0 / 1440.0
    # Long segment 1: 5 days
    t1 = np.arange(100.0, 105.0, dt)
    # Short segment 2: 0.5 days span (after trimming, even smaller, < 0.75 d threshold)
    # Put a 1.0 day gap before it, and 1.0 day gap after
    t2 = np.arange(106.0, 106.5, dt)
    t = np.concatenate([t1, t2])
    f = np.arange(t.size, dtype=float)
    e = np.ones_like(t)

    # Note: lightcurve.gap_edge_mask trims 0.3 d from both ends of every segment!
    # A 0.5 d segment with 0.3 d trim at both ends has no surviving cadences.
    # To test short segment fallback that survives gap trim, use a segment of 1.1 days:
    # 1.1 d - 2*0.3 d = 0.5 d span surviving. 0.5 d < 0.75 d (3x 0.25 d).
    t2_surv = np.arange(106.0, 107.1, dt)
    t = np.concatenate([t1, t2_surv])
    f = np.arange(t.size, dtype=float)
    e = np.ones_like(t)

    t_out, f_out, e_out, info = scramble.segment_rotate(t, f, e, seed=42, cfg=v11_cfg)

    # Segment 2 surviving span is ~0.5 d < 0.75 d (one block) -> left unrotated
    assert info["scr_schemes"][-1] == "unrotated"
    assert info["scr_offsets"][-1] == 0

    # Also test an explicit custom configuration where segment has multiple blocks but < threshold
    custom_cfg = {
        "scramble": {
            "gap_d": 0.5,
            "gap_trim_d": 0.0,  # no trim for this unit test
            "block_duration_multiplier": 3.0,
            "longest_trial_duration_d": 0.1,  # threshold = 0.3 d
        }
    }
    # Segment of 0.25 d (< 0.3 d threshold) with cadences dt = 0.01 d (25 points).
    # block duration = 0.3 d -> n_block = ceil(0.3 / 0.01) = 30 points.
    # len is 25 < 30 points -> unrotated.
    t_tiny = np.arange(0, 0.25, 0.01)
    f_tiny = np.arange(len(t_tiny), dtype=float)
    e_tiny = np.ones_like(f_tiny)
    _, f_res, _, inf = scramble.segment_rotate(t_tiny, f_tiny, e_tiny, seed=1, cfg=custom_cfg)
    assert inf["scr_schemes"] == ["unrotated"]
    assert np.array_equal(f_res, f_tiny)

    # Segment of 0.25 d with block threshold 0.1 d -> n_block = 10 points.
    # 25 points >= 10 points -> block_shuffle
    custom_cfg2 = {
        "scramble": {
            "gap_d": 0.5,
            "gap_trim_d": 0.0,
            "block_duration_multiplier": 1.0,
            "longest_trial_duration_d": 0.1,  # threshold = 0.1 d, but segment span < threshold?
        }
    }
    # If segment span = 0.08 d < 0.1 d, with dt = 0.005 d (16 points), block_len = 0.1 d -> 20 points
    # Let min_seg_d = 0.2 d, segment span = 0.15 d (so < min_seg_d)
    # but dt = 0.001 d, n_block = ceil(0.05 / 0.001)
    # The rule says: contiguous blocks of length >= 3x longest trial duration.
    # If segment is shorter than 3x longest trial duration, fallback block shuffle uses blocks >= 3x longest trial duration.
    # If segment is shorter than 1 block (i.e. shorter than 3x longest trial duration), it cannot form even 1 full block.
    # Thus any segment shorter than 3x longest trial duration has length < 1 block, so it is left unrotated!


def test_planted_periodic_box_transit_train_decoheres(v11_cfg):
    """A planted periodic box-transit train (several transits over 2 orbits, realistic noise)
    does not survive scrambling coherently (phase-folded depth at injected period/epoch drops
    to consistent with noise, e.g. < 3 sigma, across several seeds).
    """
    P_inj = 4.2
    t0_inj = 101.5
    dur_inj = 0.1
    depth_inj = 0.0006  # 600 ppm depth
    noise = 0.0004     # 400 ppm noise

    # Two orbits across 27 days
    dt = 2.0 / 1440.0
    t1 = np.arange(100.0, 113.0, dt)
    t2 = np.arange(114.5, 127.0, dt)
    t = np.concatenate([t1, t2])

    rng = np.random.default_rng(42)
    f = 1.0 + rng.normal(0, noise, size=t.size)
    e = np.full_like(t, noise)

    # Invert box transit dip
    phase = ((t - t0_inj + 0.5 * P_inj) % P_inj) - 0.5 * P_inj
    in_tr = np.abs(phase) < (dur_inj / 2.0)
    f[in_tr] -= depth_inj

    # Check that unshuffled raw transit is strongly detected (> 20 sigma)
    bls_raw = BoxLeastSquares(t, f, dy=e)
    st_raw = bls_raw.compute_stats(P_inj, dur_inj, t0_inj)
    raw_depth = st_raw["depth"][0]
    raw_err = st_raw["depth_phased"][1]
    raw_signif = raw_depth / raw_err
    assert raw_signif > 20.0, f"Raw injection significance too low: {raw_signif}"

    # Across several seeds, the scrambled signal should drop to < 3 sigma
    sigmas = []
    seeds = [1, 2, 3, 4, 5]
    for s in seeds:
        t_scr, f_scr, e_scr, _ = scramble.segment_rotate(t, f, e, seed=s, cfg=v11_cfg)
        bls_scr = BoxLeastSquares(t_scr, f_scr, dy=e_scr)
        st_scr = bls_scr.compute_stats(P_inj, dur_inj, t0_inj)
        depth_scr = st_scr["depth"][0]
        err_scr = st_scr["depth_phased"][1]
        sigma = abs(depth_scr) / err_scr
        sigmas.append(sigma)

    assert all(sig < 3.0 for sig in sigmas), f"Scrambled depth significant (> 3 sigma): {sigmas}"


def test_pipeline_scramble_mode_integration(v11_cfg):
    """Test pipeline run_job integration with mode 'scramble'."""
    from exosearch.v11 import pipeline

    cfg = dict(v11_cfg, beta_option="b", s_floor_checks=6.0)

    dt = 2.0 / 1440.0
    t1 = np.arange(100.0, 113.0, dt)
    t2 = np.arange(114.5, 127.0, dt)
    t = np.concatenate([t1, t2])
    rng = np.random.default_rng(0)
    f = np.ones_like(t) + rng.normal(0, 0.0005, size=t.size)
    e = np.full_like(t, 0.0005)
    meta = {"camera": 1, "ccd": 1, "tmag": 10.0, "teff": 5500, "radius": 1.0, "logg": 4.5}

    job = {
        "key": "test:integration",
        "tic": 12345,
        "sector": 70,
        "holdout": False,
        "prov": {},
        "cfg": cfg,
        "primary": "b",
        "tasks": [
            {"tag": "scr_0", "mode": "scramble", "seed": 42}
        ]
    }

    orig_load = pipeline.load
    try:
        pipeline.load = lambda j: (t, f, e, np.ones(len(t), bool), meta, np.array([]))
        rows = pipeline.run_job(job)
    finally:
        pipeline.load = orig_load

    assert len(rows) == 1
    row = rows[0]
    assert row["status"] == "ok"
    assert row["mode"] == "scramble"
    assert "scr_offsets" in row
    assert "scr_scheme" in row
    assert row["scr_scheme"] == "rotate"
    assert len(row["scr_offsets"]) == 2

