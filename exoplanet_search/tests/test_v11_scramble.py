"""Tests for v1.1 scrambled-null transform (exosearch.v11.scramble)."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
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
    # undo the rule-(7) join gaps (each wrap shifts every later cadence by scr_join_gap_d)
    t_out = t_out - info["scr_join_gap_d"] * np.sum(t_out[:, None] >= np.array(info["scr_wrap_t1"])[None, :], axis=1)

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
    # rule (7): time stamps move only by the inserted join gaps (cumulative), never reorder
    shift = np.zeros(t_kept.size)
    starts = np.r_[0, np.where(np.diff(t_kept) > 0.5)[0] + 1]
    for a, k in zip(starts, info["scr_offsets"]):
        if k > 0:
            shift[a + k:] += info["scr_join_gap_d"]
    assert np.allclose(t_out, t_kept + shift, rtol=0, atol=1e-9)
    assert np.all(np.diff(t_out) > 0)

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


def test_scramble_exclude_hosts_file_loading():
    """Verify that the exclude file loads, has 285 unique TICs, and expected columns."""
    cfg = settings.cfg()
    exclude_file_rel = cfg["scramble"]["exclude_hosts_file"]["70"]
    exclude_path = settings.ROOT / exclude_file_rel
    assert exclude_path.exists(), f"Exclude file not found: {exclude_path}"

    df = pd.read_csv(exclude_path)
    assert list(df.columns) == ["tic", "toi_host", "tess_eb"]
    assert len(df) == 285
    assert df["tic"].nunique() == 285
    assert (df["toi_host"] == 1).sum() == 222
    assert (df["tess_eb"] == 1).sum() == 67


def test_manifest_build_skips_scramble_for_excluded_tics(monkeypatch, tmp_path):
    """Manifest-build test that listed TICs get no scramble task and unlisted ones do."""
    from exosearch.v11 import run

    # Read real exclude list
    cfg = settings.cfg()
    exclude_file_rel = cfg["scramble"]["exclude_hosts_file"]["70"]
    df_ex = pd.read_csv(settings.ROOT / exclude_file_rel)
    excluded_tic = int(df_ex["tic"].iloc[0])
    unlisted_tic = 999999999  # definite dummy unlisted TIC

    # Mock store.connect, store.frame, targets.fetch_manifest, star_params, inject.select_hosts
    class DummyCon:
        pass

    monkeypatch.setattr("exosearch.store.connect", lambda: DummyCon())

    store_df = pd.DataFrame([
        {"tic": excluded_tic, "status": "ok", "radius": 1.0, "tmag": 10.0},
        {"tic": unlisted_tic, "status": "ok", "radius": 1.0, "tmag": 10.0},
    ])
    monkeypatch.setattr("exosearch.store.frame", lambda con, s: store_df.copy())

    targets_df = pd.DataFrame([
        {"tic": excluded_tic, "fname": "tess_excluded.fits"},
        {"tic": unlisted_tic, "fname": "tess_unlisted.fits"},
    ])
    monkeypatch.setattr("exosearch.targets.fetch_manifest", lambda s: targets_df.copy())

    monkeypatch.setattr("exosearch.v11.run.star_params", lambda s, tics, log=print: {t: {} for t in tics})
    monkeypatch.setattr("exosearch.inject.select_hosts", lambda s, cfg, con: pd.DataFrame(columns=["tic"]))

    jobs = run.build([70], [("h1", 1.0)], tmp_path, production=True, write=False)
    assert len(jobs) == 2

    job_by_tic = {j["tic"]: j for j in jobs}
    excluded_modes = [t["mode"] for t in job_by_tic[excluded_tic]["tasks"]]
    unlisted_modes = [t["mode"] for t in job_by_tic[unlisted_tic]["tasks"]]

    # Listed TIC gets real and invert, but NO scramble tasks
    assert "real" in excluded_modes
    assert "invert" in excluded_modes
    assert "scramble" not in excluded_modes

    # Unlisted TIC gets real, invert, and scramble tasks
    assert "real" in unlisted_modes
    assert "invert" in unlisted_modes
    assert "scramble" in unlisted_modes
    assert excluded_modes.count("real") == 1
    assert excluded_modes.count("invert") == 1


def test_planted_signal_real_cadence_check(v11_cfg):
    """Planted-signal test using real sector-70 cadence and gap structure from s0070_cadence.npz.

    Cyclic rotation preserves periodicity inside each segment, so scrambling does not
    erase coherent transits (max sigma stays > 3 at fixed (P, t0) across seeds).
    This is why known hosts (TOIs and EBs) are explicitly excluded from the scrambled null.

    This test asserts that scramble_planted_check is deterministic and exactly reproduces
    the committed sigma table:
      P=1.0: median 3.70, max 7.70
      P=3.0: median 1.81, max 7.26
      P=6.0: median 0.97, max 7.82
    """
    from scripts.scramble_planted_check import run_planted_check

    fixture_path = Path(__file__).resolve().parent / "fixtures" / "s0070_cadence.npz"
    results = run_planted_check(fixture_path, output_file=None, n_seeds=10, periods=(1.0, 3.0, 6.0))

    # Amendment 2026-11-12 (join gaps, rule 7): the inserted gaps move every cadence after a wrap by 0.5 d + 1 cadence,
    # which also breaks a planted signal's phase coherence across the wrap. Before the amendment (fddfb25):
    # P=1.0 3.70/7.70, P=3.0 1.81/7.26, P=6.0 0.97/7.82.
    expected = {
        1.0: {"median_sigma": 3.98, "max_sigma": 5.23},
        3.0: {"median_sigma": 1.57, "max_sigma": 3.34},
        6.0: {"median_sigma": 1.40, "max_sigma": 2.88},
    }

    assert len(results) == 3
    for r in results:
        exp = expected[r["P"]]
        assert round(r["median_sigma"], 2) == exp["median_sigma"], f"P={r['P']} median mismatch"
        assert round(r["max_sigma"], 2) == exp["max_sigma"], f"P={r['P']} max mismatch"
        assert r["n_seeds"] == 10


def test_pipeline_scramble_mode_integration(v11_cfg, monkeypatch):
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

    monkeypatch.setattr(pipeline, "load", lambda j: (t, f, e, np.ones(len(t), bool), meta, np.array([])))
    trims = []
    orig_prepare = pipeline.prepare
    monkeypatch.setattr(pipeline, "prepare", lambda *a, **k: (trims.append(k.get("gap_trim")), orig_prepare(*a, **k))[1])
    rows = pipeline.run_job(job)
    assert trims == [0.0]                         # scrambled path: the detrender's gap trim is off (no double trimming)

    assert len(rows) == 1
    row = rows[0]
    assert row["status"] == "ok", row.get("error")
    assert row["mode"] == "scramble"
    assert "scr_offsets" in row
    assert "scr_scheme" in row
    assert row["scr_scheme"] == "segment_rotation"
    assert len(row["scr_offsets"]) == 2 and all(o >= 1 for o in row["scr_offsets"])


def test_join_gap_equals_detrender_break_tolerance_plus_cadence(v11_cfg):
    """Amendment 2026-11-12: inserted offset = the detrender's break_tolerance + one cadence, exactly."""
    import inspect, re
    src = inspect.getsource(lightcurve.detrend)
    bt = float(re.search(r"break_tolerance=([0-9.]+)", src).group(1))
    assert scramble.BREAK_TOLERANCE_D == bt == 0.5
    t, f, e = _make_two_orbit_lc()
    _, _, _, info = scramble.segment_rotate(t, f, e, seed=5, cfg=v11_cfg)
    kept = t[lightcurve.gap_edge_mask(t, gap=0.5, trim=0.3)]
    cad = float(np.median(np.diff(kept)))
    assert info["scr_cadence_d"] == cad
    assert info["scr_join_gap_d"] == bt + cad
    assert abs(cad - 2.0 / 1440.0) < 1e-9


def _box_snr_at(t, f, edge_t, D, side):
    tf, fl, _, _, _ = lightcurve.detrend(t, f, np.full_like(t, 1e-4), gap_trim=0.0)
    c = edge_t - D / 2 if side < 0 else edge_t + D / 2
    m = np.abs(tf - c) < D / 2
    sig = 1.4826 * np.median(np.abs(fl - np.median(fl)))
    return float((np.median(fl) - fl[m].mean()) / (sig / np.sqrt(m.sum())))


@pytest.mark.parametrize("seed", [3, 7, 11])
def test_linear_ramp_wrap_gives_break_and_no_box_dip(seed):
    """A rotated linear ramp: the detrender breaks at the wrap and no box dip appears there (beyond the detrender's
    ordinary edge effect at a real segment start); without the inserted gap the wrap step makes a strong dip."""
    from wotan.gaps import get_gaps_indexes
    dt = 2.0 / 1440.0
    t = np.arange(100.0, 112.0, dt)
    rng = np.random.default_rng(1)
    f = 1 + 1e-3 * (t - 100.0) / 12.0 + rng.normal(0, 1e-4, t.size)
    e = np.full_like(t, 1e-4)
    cfg = {"scramble": {"gap_d": 0.5, "gap_trim_d": 0.3}}
    to, fo, eo, info = scramble.segment_rotate(t, f, e, seed, cfg)
    assert len(info["scr_wrap_t1"]) == 1
    i = int(np.searchsorted(to, info["scr_wrap_t1"][0]))
    assert i in set(get_gaps_indexes(to, scramble.BREAK_TOLERANCE_D).tolist())       # detrender breaks at the wrap
    tk = t[lightcurve.gap_edge_mask(t, gap=0.5, trim=0.3)]
    fk = f[lightcurve.gap_edge_mask(t, gap=0.5, trim=0.3)]
    for D in (0.04, 0.08, 0.25):
        fixed = max(_box_snr_at(to, fo, info["scr_wrap_t0"][0], D, -1), _box_snr_at(to, fo, info["scr_wrap_t1"][0], D, +1))
        unfixed = max(_box_snr_at(tk, fo, tk[i - 1], D, -1), _box_snr_at(tk, fo, tk[i], D, +1))
        real_edge = max(_box_snr_at(tk, fk, tk[0], D, +1), _box_snr_at(tk, fk, tk[-1], D, -1))
        assert unfixed > 10.0
        assert fixed < 3.0
        assert fixed <= real_edge + 0.5
