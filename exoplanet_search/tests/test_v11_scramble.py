"""Tests for v1.1 scrambled-null transform (exosearch.v11.scramble)."""
from __future__ import annotations

from pathlib import Path

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
    assert info["scr_scheme"] == "segment_rotation"


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


def test_single_sample_segment_rotates_trivially_by_zero(v11_cfg):
    """A 1-sample segment trivially rotates by 0 and produces valid output."""
    t = np.array([100.0, 100.001, 101.5])  # gap between index 1 and 2 is 1.499 d > 0.5 d
    # Use gap_trim_d = 0 so the points survive
    cfg = {"scramble": {"gap_d": 0.5, "gap_trim_d": 0.0}}
    f = np.array([1.0, 2.0, 3.0])
    e = np.array([0.1, 0.1, 0.1])
    t_out, f_out, e_out, info = scramble.segment_rotate(t, f, e, seed=42, cfg=cfg)

    assert info["n_segments"] == 2
    # Second segment has 1 sample; offset must be 0
    assert info["scr_offsets"][1] == 0
    assert f_out[2] == 3.0
    assert info["scr_scheme"] == "segment_rotation"


def test_planted_signal_real_cadence_check(v11_cfg):
    """Planted-signal test using real sector-70 cadence and gap structure from s0070_cadence.npz.

    Injects box transit trains at P = 1, 3, 6 d (depth giving > 20 sigma unscrambled)
    with white noise plus mild red noise across at least 10 seeds each.
    Asserts every period's max sigma < 3.
    """
    fixture_path = Path(__file__).resolve().parent / "fixtures" / "s0070_cadence.npz"
    d = np.load(fixture_path)
    t = d["t"]

    noise_white = 0.0004
    noise_red_amp = 0.0001
    dur = 0.1
    t0 = float(t[0] + 0.1)
    depth = 0.0006
    n_seeds = 10
    periods = [1.0, 3.0, 6.0]

    for P in periods:
        # Check unscrambled significance > 20 sigma
        rng_u = np.random.default_rng(0)
        w_u = rng_u.normal(0, noise_white, size=len(t))
        r_u = noise_red_amp * np.sin(2 * np.pi * (t - t[0]) / 3.5)
        f_u = 1.0 + w_u + r_u
        e_u = np.full_like(t, noise_white)
        phase_u = ((t - t0 + 0.5 * P) % P) - 0.5 * P
        f_u[np.abs(phase_u) < dur / 2.0] -= depth
        bls_u = BoxLeastSquares(t, f_u, dy=e_u)
        st_u = bls_u.compute_stats(P, dur, t0)
        unscr_sigma = float(st_u["depth"][0] / st_u["depth_phased"][1])
        assert unscr_sigma > 20.0, f"Unscrambled significance {unscr_sigma:.2f} <= 20 for P={P}"

        sigmas = []
        for s in range(n_seeds):
            rng = np.random.default_rng(s)
            w = rng.normal(0, noise_white, size=len(t))
            r = noise_red_amp * np.sin(2 * np.pi * (t - t[0]) / 3.5)
            f = 1.0 + w + r
            e = np.full_like(t, noise_white)
            phase = ((t - t0 + 0.5 * P) % P) - 0.5 * P
            f[np.abs(phase) < dur / 2.0] -= depth

            t_scr, f_scr, e_scr, _ = scramble.segment_rotate(t, f, e, seed=s, cfg=v11_cfg)
            bls = BoxLeastSquares(t_scr, f_scr, dy=e_scr)
            st = bls.compute_stats(P, dur, t0)
            sig = float(abs(st["depth"][0]) / st["depth_phased"][1])
            sigmas.append(sig)

        max_sig = max(sigmas)
        assert max_sig < 3.0, f"Period {P} max sigma {max_sig:.2f} >= 3.0 across {n_seeds} seeds"


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
    assert row["scr_scheme"] == "segment_rotation"
    assert len(row["scr_offsets"]) == 2
