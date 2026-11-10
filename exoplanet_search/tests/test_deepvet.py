import numpy as np

from exosearch import deepvet


def _stamp(eclipse_on, n=2000, P=2.0, t0=0.3, dur=0.12, depth=0.3):
    rng = np.random.default_rng(0)
    t = np.linspace(0, 20, n)
    yy, xx = np.mgrid[0:11, 0:11]
    psf = lambda y, x: np.exp(-((yy - y) ** 2 + (xx - x) ** 2) / 2.0)
    tgt, nb = 1000 * psf(5, 5), 200 * psf(5, 8.5)
    ph = ((t - t0 + 0.5 * P) % P) - 0.5 * P
    intr = np.abs(ph) < dur / 2
    flux = np.empty((n, 11, 11))
    for i in range(n):
        a, b = tgt.copy(), nb.copy()
        if intr[i]:
            if eclipse_on == "target":
                a *= 1 - depth * 0.2
            else:
                b *= 1 - depth
        flux[i] = a + b + rng.normal(0, 1.0, (11, 11))
    return t, flux, P, t0, dur


def test_centroid_on_target_passes():
    t, f, P, t0, dur = _stamp("target")
    r = deepvet.centroid_core(t, f, P, t0, dur, (5, 5))
    assert r["n_epochs"] >= 8 and r["diff_peak_snr"] > 3
    assert "centroid_offset" not in deepvet.flags_for(r) and "dip_off_target" not in deepvet.flags_for(r)


def test_centroid_on_neighbour_flagged():
    t, f, P, t0, dur = _stamp("neighbour")
    r = deepvet.centroid_core(t, f, P, t0, dur, (5, 5))
    fl = deepvet.flags_for(r)
    assert "dip_off_target" in fl and "centroid_offset" in fl, r


def test_mimic_threshold():
    assert abs(deepvet.max_dm_to_mimic(0.01) - 4.99) < 0.02      # 1% dip: neighbour up to ~5 mag fainter
    fl = deepvet.flags_for(dict(diff_peak_snr=10, centroid_offset_px=0.1, centroid_offset_sigma=1,
                                diff_peak_dist_px=0.2, gaia_n_can_mimic=1, gaia_mimic_min_sep=3.0))
    assert fl == ["unresolved_neighbour"]


def test_eb_flags_and_harmonics():
    assert deepvet._harmonic(2.0, 4.0) and not deepvet._harmonic(2.0, 3.1)
    fl = deepvet.flags_for(dict(eb_match="TESS-EB TIC 1", v_ratio=0.3, v_sigma=20, sec_scan_sigma=6, rp_rjup=2.5))
    assert set(fl) >= {"known_eb", "v_shape", "secondary_scan", "rp_gt_2rjup"}
    assert deepvet.flags_for(dict(v_ratio=0.9, v_sigma=20, sec_scan_sigma=1, rp_rjup=0.3)) == []
