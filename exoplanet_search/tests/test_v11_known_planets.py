import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from exosearch.v11 import known_planets, pipeline, search, settings

from lc_fixtures import s70_like_times, transit_flux, write_spoc_fits


def test_grid_halfwidth_fallback():
    """locate's own fallback 0.5 * P^2 * df when no plateau information exists for the star."""
    hw = known_planets.calculate_grid_halfwidth(3.0, baseline=27.0)
    assert np.isclose(hw, 0.5 * 9.0 * 0.04 / 27.0 ** 2)


def test_grid_halfwidth_uses_star_plateau_not_other_tic():
    """r6 #5: the star's own plateau half-width; another TIC's row in p_grid_halfwidth_v1.csv is never used."""
    plats = [[2.990, 3.010, 21, 0.0105]]
    assert np.isclose(known_planets.get_plateau_halfwidth(3.0, plateaus_list=json.dumps(plats), grid_df=1e-5), 0.0105)
    # off-plateau period -> locate fallback with the star's own grid spacing
    assert np.isclose(known_planets.get_plateau_halfwidth(5.0, plateaus_list=plats, grid_df=1e-5), 0.5 * 25 * 1e-5)
    tab = pd.DataFrame([dict(tic=999, sector=70, period_v10=3.0, P_grid_halfwidth=0.02, P_grid_halfwidth_in_band=True),
                        dict(tic=7, sector=70, period_v10=3.0005, P_grid_halfwidth=0.004, P_grid_halfwidth_in_band=True)])
    fb = 0.5 * 9.0 * 0.04 / 27.0 ** 2
    assert np.isclose(known_planets.get_plateau_halfwidth(3.0, tic=5, sector=70, pgh_table=tab), fb)
    assert np.isclose(known_planets.get_plateau_halfwidth(3.0, tic=7, sector=70, pgh_table=tab), 0.004)
    assert np.isclose(known_planets.get_plateau_halfwidth(3.0, tic=7, sector=71, pgh_table=tab), fb)


def test_alias_reported_not_recovered():
    """2x / 0.5x use the same half-width, are reported as match types, and are not recoveries."""
    P, hw = 4.0, 0.01
    assert known_planets.check_recovery(P, P + 0.5 * hw, hw) == (True, "direct")
    assert known_planets.check_recovery(P, 2 * P + 0.5 * hw, hw) == (False, "2x")
    assert known_planets.check_recovery(P, 0.5 * P - 0.5 * hw, hw) == (False, "0.5x")
    assert known_planets.check_recovery(P, 2 * P + 2 * hw, hw) == (False, "none")
    assert known_planets.check_recovery(P, P, hw, s_found=10.0, s_star=12.0) == (False, "direct")
    assert known_planets.check_recovery(P, P, hw, s_found=13.0, s_star=12.0) == (True, "direct")


def test_s_star_required():
    df_tgt = pd.DataFrame([{"tic": 1, "period": 3.0, "dur_h": 2.0, "depth_ppm": 1000.0, "t0": 1.0,
                            "aperture_contamination": 0.9}])
    for bad in (None, -1.0, 0.0):
        with pytest.raises(ValueError, match="s_star must be provided"):
            known_planets.evaluate_known_planets(df_tgt, pd.DataFrame(), s_star=bad)


def test_contamination_flag_and_missing_aperture():
    base = dict(period=3.0, t0=1.0, dur_h=2.4, depth_ppm=2000.0, n_raw=3, n_tr=3, cdpp_d=1e-4)
    df_tgt = pd.DataFrame([dict(base, tic=1, name="P1", aperture_contamination=0.95, contratio=0.05),
                           dict(base, tic=2, name="P2", aperture_contamination=0.50, contratio=0.01),
                           dict(base, tic=3, name="P3", aperture_contamination=np.nan, contratio=0.05)])
    out, _ = known_planets.evaluate_known_planets(df_tgt, pd.DataFrame(), s_star=10.0)
    r = out.set_index("tic")
    assert r.loc[1, "eligibility"] == "eligible" and not r.loc[1, "contamination_flag"]
    assert r.loc[2, "eligibility"] == "eligible" and r.loc[2, "contamination_flag"]
    assert r.loc[3, "eligibility"] == "ineligible" and r.loc[3, "ineligible_reason"] == "missing_aperture_contamination"
    # expected S = depth * aperture * sqrt(n_tr) / CDPP_D
    assert np.isclose(r.loc[1, "expected_S"], 2000e-6 * 0.95 * np.sqrt(3) / 1e-4)


def test_threshold_applied_to_frozen_expected_s():
    base = dict(t0=1.0, dur_h=2.0, depth_ppm=1000.0, n_raw=4, n_tr=4, cdpp_d=1e-4, aperture_contamination=1.0)
    df_tgt = pd.DataFrame([dict(base, tic=1, name="a", period=2.0, depth_ppm=600.0),     # S = 12
                           dict(base, tic=2, name="b", period=3.0)])                      # S = 20
    lo, _ = known_planets.evaluate_known_planets(df_tgt, pd.DataFrame(), s_star=9.70)
    hi, _ = known_planets.evaluate_known_planets(df_tgt, pd.DataFrame(), s_star=14.52)
    assert list(lo.eligibility) == ["eligible", "eligible"]
    assert list(hi.eligibility) == ["ineligible", "eligible"]
    assert hi.iloc[0].ineligible_reason == "expected_s_below_threshold"


def test_select_targets_pscomppars_first_with_toi_fill_and_dedupe():
    toi = pd.DataFrame([
        dict(tic=101, toi="101.01", period=3.14159, disp="CP", dur_h=2.5, depth_ppm=1200.0, epoch_bjd=2460000.1),
        dict(tic=102, toi="102.01", period=5.0, disp="PC", dur_h=3.0, depth_ppm=800.0, epoch_bjd=2460000.5),
        dict(tic=103, toi="103.01", period=1.5, disp="FP", dur_h=1.0, depth_ppm=500.0, epoch_bjd=2460000.2),
        dict(tic=104, toi="104.01", period=2.2, disp="", dur_h=1.5, depth_ppm=700.0, epoch_bjd=2460000.3),
    ])
    ps = pd.DataFrame([dict(tic=101, pl_name="Planet-101 b", pl_orbper=3.14160, pl_trandep=0.15, pl_trandur=np.nan,
                            pl_tranmid=2460000.12)])
    out = known_planets.select_known_planet_targets(sector=70, searched_tics={101, 102, 103, 104}, ps_df=ps,
                                                    toi_df=toi, pgh_table=pd.DataFrame())
    r = out.set_index("tic")
    assert sorted(r.index) == [101, 102, 104]                   # FP dropped, TOI 101.01 deduped onto the PS row
    assert r.loc[101, "type"] == "ps" and r.loc[101, "name"] == "Planet-101 b"
    assert np.isclose(r.loc[101, "depth_ppm"], 1500.0) and r.loc[101, "depth_src"] == "pscomppars"
    assert np.isclose(r.loc[101, "dur_h"], 2.5) and r.loc[101, "dur_src"] == "exofop_toi 101.01"
    assert np.isclose(r.loc[101, "t0"], 3000.12)
    assert r.loc[102, "name"] == "TOI 102.01" and r.loc[102, "depth_src"] == "exofop_toi"


def test_pscomppars_cache_requires_transit_columns(tmp_path):
    p = tmp_path / "nasa_pscomppars_full_20261112.csv"
    pd.DataFrame(dict(pl_name=["x"], tic_id=["TIC 5"], pl_orbper=[1.0], disc_facility=["TESS"])).to_csv(p, index=False)
    Path(str(p) + ".meta.json").write_text(json.dumps(dict(sha256=known_planets._sha256(p))))
    with pytest.raises(ValueError, match="lacks columns"):
        known_planets.load_pscomppars(p)
    with pytest.raises(FileNotFoundError):
        known_planets.load_pscomppars(cache_dir=tmp_path / "empty")


def test_report_refuses_without_frozen_eligibility_file(tmp_path):
    with pytest.raises(RuntimeError, match="Frozen eligibility file required"):
        known_planets.run_known_planet_recovery(s_star=10.0, out_dir=tmp_path, frozen_eligibility_file=None)


def _search_row(tic, lc_path, cfg):
    """Real v1.1 search (option b) on a synthetic light curve -> one real-search row."""
    from exosearch import lightcurve
    t, f, e, q, _ = lightcurve.load_spoc_raw(lc_path)
    good = lightcurve.good_mask(f, e, q)
    med = np.nanmedian(f[good])
    d = pipeline.prepare(t[good], f[good] / med, e[good] / med, cfg, prewhiten_on=True, gap_trim=0.3)
    best, _ = search.search(*d["binned"], *d["flat"][:2], cfg=cfg, options=("b",), primary="b")
    return dict(tic=tic, mode="real", status="ok", period=best["period"], P_grid=best["P_grid"], S=best["S"],
                plateaus=best["plateaus"], grid_df=best["grid_df"], passed_lc_checks=1)


def test_freeze_then_report_end_to_end(tmp_path):
    """Frozen CSV round trip into run_known_planet_recovery on synthetic light curves and a real v1.1 search:
    eligibility and P_cat survive the round trip, a direct match is recovered, an alias is reported and fails,
    a missing light curve is a required error row (never ineligible), and expected S is never recomputed."""
    cfg = settings.cfg()
    P, t0, dur, depth = 3.71, 3001.3, 0.11, 2500e-6
    t = s70_like_times()
    lc = write_spoc_fits(tmp_path / "lc1.fits", t, transit_flux(t, P, t0, dur, depth), crowdsap=0.98, tic=1)
    row1 = _search_row(1, lc, cfg)
    # Catalogue period placed on the star's own reported plateau (this test checks the wiring and the round trip;
    # the grid's period accuracy is measured on the real targets, not assumed here).
    P = row1["period"]
    tg = pd.DataFrame([
        dict(tic=1, name="TOI 1.01", type="toi", period=P, t0=t0, dur_h=dur * 24, depth_ppm=depth * 1e6),
        dict(tic=2, name="TOI 2.01", type="toi", period=2.5, t0=3000.4, dur_h=2.0, depth_ppm=3000.0),
        dict(tic=3, name="TOI 3.01", type="toi", period=P / 2, t0=t0, dur_h=dur * 24, depth_ppm=depth * 1e6),
    ])
    for c in ("depth_src", "dur_src", "t0_src"):
        tg[c] = "exofop_toi"
    lookup = {1: lc, 2: None, 3: lc}.get
    frozen = tmp_path / "frozen.csv"
    df, meta = known_planets.freeze_eligibility(9.70, frozen, targets_df=tg, lc_lookup=lookup, contratio_map={},
                                                cfg=cfg, commit="test", config_hash="test")
    assert meta["sha256"] == known_planets._sha256(frozen)

    back = known_planets.read_frozen(frozen)
    r = back.set_index("tic")
    assert r.loc[1, "eligibility"] == "eligible" and np.isclose(r.loc[1, "period"], P)
    assert r.loc[1, "expected_S"] > 20 and r.loc[1, "n_tr"] >= 6
    assert r.loc[2, "eligibility"] == "error" and r.loc[2, "error"] == "missing light curve" and r.loc[2, "required"]
    assert r.loc[3, "eligibility"] == "eligible"

    res = pd.DataFrame([row1, dict(row1, tic=3)])
    res_p = tmp_path / "real_rows.csv"
    res.to_csv(res_p, index=False)
    for s_star in (9.70, 14.52):
        out, summ = known_planets.run_known_planet_recovery(s_star, res_p, tmp_path / f"out{s_star}", sector=70,
                                                            frozen_eligibility_file=frozen, cfg=cfg)
        o = out.set_index("tic")
        assert o.loc[1, "recovered"] and o.loc[1, "match_type"] == "direct"
        assert np.isclose(o.loc[1, "expected_S"], r.loc[1, "expected_S"])          # frozen, not recomputed
        assert o.loc[3, "match_type"] == "2x" and not o.loc[3, "recovered"] and o.loc[3, "stage"] == "outranked"
        assert o.loc[2, "required"] and not o.loc[2, "recovered"] and o.loc[2, "stage"] == "error"
        assert summ["required_targets"] == 3 and summ["eligible_recovered"] == 1 and not summ["pass"]
        assert summ["required_aliases"] == 1

    # tampering with the frozen file is refused
    frozen.write_text(frozen.read_text().replace("TOI 1.01", "TOI 9.01"))
    with pytest.raises(ValueError, match="sha256 mismatch"):
        known_planets.run_known_planet_recovery(9.70, res_p, tmp_path / "x", frozen_eligibility_file=frozen, cfg=cfg)


def test_runtime_error_during_freeze_is_error_row(tmp_path):
    bad = tmp_path / "bad.fits"
    bad.write_bytes(b"not a fits file")
    tg = pd.DataFrame([dict(tic=9, name="TOI 9.01", type="toi", period=2.0, t0=3000.5, dur_h=2.0, depth_ppm=2000.0,
                            depth_src="exofop_toi", dur_src="exofop_toi", t0_src="exofop_toi")])
    df, _ = known_planets.freeze_eligibility(10.0, tmp_path / "f.csv", targets_df=tg, lc_lookup=lambda t: bad,
                                             contratio_map={}, commit="x", config_hash="x")
    row = df.iloc[0]
    assert row.eligibility == "error" and row.required and row.error and row.ineligible_reason == ""


def test_losing_stage():
    ls = known_planets.losing_stage
    assert ls("none", 5, 10, True, error="missing light curve") == "error"
    assert ls("none", 20, 10, True, n_tr_cat=1) == "masked"
    assert ls("direct", 8, 10, True, depth_ratio_cat=0.3) == "detrending"
    assert ls("direct", 8, 10, True, depth_ratio_cat=0.9) == "below_S*"
    assert ls("none", 30, 10, True, s_cat=15) == "outranked"
    assert ls("2x", 30, 10, True, s_cat=6, depth_ratio_cat=0.9) == "below_S*"
    assert ls("direct", 20, 10, False) == "vetting"
    assert ls("direct", 20, 10, True) == ""
