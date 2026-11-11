"""Independent recovery check of a candidate in TESS FFI cutouts (TESScut), CPU only. Analysis layer; v1.0.0 untouched.

Photometry: threshold aperture on a 15x15 cutout, per-cadence background from faint out-of-aperture pixels,
pixel-level decorrelation (lightkurve PLDCorrector), cadences with high background or bad quality dropped.
Then the frozen v1.0.0 detrend (biweight 0.75 d, gap trim) and the frozen S search, run *blind* per sector and on
the stitched sectors, before looking at the candidate ephemeris. The candidate's own 2-min sector is also run
through the identical FFI photometry as a positive control.
At the candidate period: box depth from per-epoch depths (red-noise-aware error = scatter of epoch depths),
S at fixed (P, T0, D), odd/even, and the ephemeris link between the FFI sectors and the reference sector
(T0 / P uncertainties from the Delta chi^2 = 1 width of the box-S^2 scan).
  python -m exosearch.ffi_check --tic 318707987 --period 0.9595540 --t0 3210.0051 --dur 0.06 --sectors 42,43 --control 70
"""
from __future__ import annotations

import argparse
import json
import warnings
from pathlib import Path

import numpy as np

from . import bls, config, task

warnings.filterwarnings("ignore")


class NotOnSilicon(RuntimeError):
    pass


def on_silicon(tpf, tmag):
    """TESScut footprints are approximate: a cutout can be served for a sector in which the star fell off the
    science pixels (gap / collateral), giving an image of ~0 e/s. Compare the brightest central pixel with the flux
    expected for the star (~15000 e/s at Tmag 10; a single pixel holds >~ 20% of it)."""
    med = np.nanmedian(tpf.flux.value, axis=0)
    c = med.shape[0] // 2
    peak = float(np.nanmax(med[c - 2:c + 3, c - 2:c + 3]))
    expected = 15000.0 * 10 ** (-0.4 * (tmag - 10.0)) if tmag is not None else None
    ok = peak > 0.01 * expected if expected else peak > 10.0
    return ok, dict(central_peak_e_per_s=peak, expected_total_e_per_s=expected, image_median_e_per_s=float(np.nanmedian(med)))


def photometry(tic, sector, cache, size=15, tmag=None):
    import lightkurve as lk
    sr = lk.search_tesscut(f"TIC {tic}", sector=sector)
    tpf = sr.download(cutout_size=size, download_dir=str(cache), quality_bitmask="hard")
    ok_si, si = on_silicon(tpf, tmag)
    if not ok_si:
        raise NotOnSilicon(json.dumps(si))
    ap = tpf.create_threshold_mask(threshold=5, reference_pixel="center")
    flux = tpf.flux.value
    med = np.nanmedian(flux, axis=0)
    bkg_mask = (~ap) & (med < np.nanpercentile(med, 50))
    bkg = np.nanmedian(flux[:, bkg_mask], axis=1)
    tpf_b = tpf - bkg[:, None, None] * tpf.flux.unit
    ok = np.isfinite(bkg) & (bkg < np.nanpercentile(bkg, 90))          # drop the worst scattered light
    lc = tpf_b[ok].to_lightcurve(aperture_mask=ap)
    try:
        corr = lk.correctors.PLDCorrector(tpf_b[ok], aperture_mask=ap)
        lc = corr.correct(pca_components=3, aperture_mask=ap, restore_trend=True)
    except Exception as exc:                                            # PLD failure -> simple aperture photometry
        print(f"  S{sector}: PLD failed ({exc}); using SAP")
    lc = lc.remove_nans().normalize()
    t = np.asarray(lc.time.value, float)
    f = np.asarray(lc.flux.value, float)
    e = np.asarray(lc.flux_err.value, float)
    e = np.where(np.isfinite(e) & (e > 0), e, np.nanstd(f))
    return t, f, e, dict(n_ap=int(ap.sum()), cadence_min=float(np.median(np.diff(t)) * 1440), n=int(t.size), sector=sector)


def prepared(t, f, e, sector):
    meta = dict(sector=sector, camera=None, radius=None)
    cad = np.median(np.diff(t)) * 1440
    return task.prepare_arrays(t, f, e, meta, bin_minutes=10.0 if cad < 9 else 0)


def box_s(tf, ff, cdpp_d, P, t0, D):
    ph = ((tf - t0 + 0.5 * P) % P) - 0.5 * P
    inn = np.abs(ph) < D / 2
    if inn.sum() < 3:
        return 0.0, 0.0, 0.0
    depth = np.median(ff[~inn]) - np.mean(ff[inn])
    cad = np.median(np.diff(tf))
    ntr = inn.sum() * cad / D
    return float(depth * np.sqrt(ntr) / cdpp_d), float(depth), float(ntr)


def epoch_depths(tf, ff, P, t0, D):
    out = []
    n = np.round((tf - t0) / P)
    for k in np.unique(n):
        tc = t0 + k * P
        dt = np.abs(tf - tc)
        ins, oot = dt < D / 2, (dt > 0.75 * D) & (dt < 2.0 * D)
        if ins.sum() >= 3 and oot.sum() >= 6:
            out.append((int(k), tc, float(np.median(ff[oot]) - np.mean(ff[ins]))))
    return np.array(out) if out else np.zeros((0, 3))


def at_ephemeris(tf, ff, P, t0, D, sigma_cdpp):
    ep = epoch_depths(tf, ff, P, t0, D)
    d = ep[:, 2]
    res = dict(n_epochs=int(len(d)))
    if len(d) < 3:
        return res
    mean, sem = float(d.mean()), float(d.std(ddof=1) / np.sqrt(len(d)))
    odd, even = d[ep[:, 0] % 2 == 1], d[ep[:, 0] % 2 == 0]
    S, depth, ntr = box_s(tf, ff, sigma_cdpp, P, t0, D)
    res.update(depth_ppm=1e6 * mean, depth_err_ppm=1e6 * sem, snr_epoch_scatter=mean / sem if sem > 0 else None,
               S_fixed=S, box_depth_ppm=1e6 * depth, n_tr=ntr,
               odd_ppm=1e6 * float(odd.mean()) if len(odd) else None, even_ppm=1e6 * float(even.mean()) if len(even) else None,
               odd_even_sigma=float(abs(odd.mean() - even.mean()) / np.hypot(odd.std(ddof=1) / np.sqrt(len(odd)), even.std(ddof=1) / np.sqrt(len(even))))
               if len(odd) > 1 and len(even) > 1 else None)
    return res


def scan_width(tf, ff, cd, P, t0, D, what, span, n=601):
    """Best value and Delta(S^2)=1 half-width of a 1-D scan in P or T0 (box statistic)."""
    grid = np.linspace(-span, span, n)
    s2 = []
    for g in grid:
        s2.append(box_s(tf, ff, cd, P + g if what == "P" else P, t0 + g if what == "t0" else t0, D)[0] ** 2)
    s2 = np.array(s2)
    i = int(np.argmax(s2))
    inside = grid[s2 >= s2[i] - 1.0]
    return float(grid[i]), float((inside.max() - inside.min()) / 2) if inside.size > 1 else float(span / n)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tic", type=int, required=True)
    ap.add_argument("--period", type=float, required=True)
    ap.add_argument("--t0", type=float, required=True)
    ap.add_argument("--dur", type=float, required=True)
    ap.add_argument("--sectors", default="42,43")
    ap.add_argument("--control", type=int, default=None, help="reference sector run through the same FFI photometry")
    ap.add_argument("--tmag", type=float, default=None, help="TESS mag, for the on-silicon flux check")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    out_dir = Path(a.out) if a.out else config.DATA_ROOT / "products" / "v1.0.0" / "analysis_layer" / "ffi_check" / f"TIC{a.tic}"
    out_dir.mkdir(parents=True, exist_ok=True)
    cache = out_dir / "cache"
    res = dict(tic=a.tic, ref=dict(period=a.period, t0=a.t0, duration=a.dur), sectors={})
    data = {}
    try:                                                                # tess-point, if installed: authoritative pointing
        from tess_stars2px import tess_stars2px_function_entry as t2p
        from astroquery.mast import Catalogs
        row = Catalogs.query_object(f"TIC {a.tic}", catalog="TIC", radius=0.001)[0]
        r = t2p(a.tic, float(row["ra"]), float(row["dec"]))
        res["tess_point_on_silicon"] = [dict(sector=int(x), camera=int(c), ccd=int(d), col=float(cc), row=float(rr))
                                        for x, c, d, cc, rr in zip(r[3], r[4], r[5], r[6], r[7])]
    except Exception as exc:
        res["tess_point_on_silicon"] = f"unavailable ({type(exc).__name__})"
    res["not_on_silicon"] = {}
    for s in [int(x) for x in a.sectors.split(",")] + ([a.control] if a.control else []):
        try:
            t, f, e, info = photometry(a.tic, s, cache, tmag=a.tmag)
        except NotOnSilicon as exc:
            res["not_on_silicon"][str(s)] = json.loads(str(exc))
            print(f"S{s}: target not on science pixels {exc}")
            continue
        d = prepared(t, f, e, s)
        data[s] = d
        tb, fb, _ = d["binned"]
        tf, ff, _ = d["flat"]
        blind = bls.search_s(tb, fb, tf, ff)
        top = dict(period=blind["period"], t0=blind["t0"], duration=blind["duration"], S=blind["S"], sde=blind["sde"],
                   depth_ppm=1e6 * blind["depth_box"], n_tr=blind["n_tr"])
        # rank of the candidate period in the blind periodogram (any harmonic within 0.5%)
        per, pw = blind["periods"], blind["power"]
        near = np.abs(per / a.period - 1) < 0.005
        rank = int((pw > pw[near].max()).sum()) + 1 if near.any() else None
        cd = blind["cdpp_by_d"].get(f"{min(bls.DURATIONS, key=lambda x: abs(x - a.dur)):.2f}")
        info.update(blind_top=top, candidate_period_rank=rank, candidate_period_peak_S=float(pw[near].max()) if near.any() else None,
                    sigma_w_ppm=1e6 * blind["sigma_w"], cdpp_at_dur_ppm=1e6 * cd)
        res["sectors"][str(s)] = info
        print(s, json.dumps(info, default=float))
    # stitched FFI sectors: blind + at the candidate period (fine local refit of P and T0 inside the FFI data only)
    ffi = [int(x) for x in a.sectors.split(",") if int(x) in data]
    if not ffi:
        res["verdict"] = ("not testable: the target was not on science pixels in any requested FFI sector "
                          f"({', '.join(sorted(res['not_on_silicon']))}); no independent recovery is possible from these sectors")
        if a.control and a.control in data:
            tr, fr = data[a.control]["flat"][0], data[a.control]["flat"][1]
            rb = bls.search_s(*data[a.control]["binned"][:2], tr, fr)
            cdr = rb["cdpp_by_d"][f"{min(bls.DURATIONS, key=lambda x: abs(x - a.dur)):.2f}"]
            res["control_at_ref"] = at_ephemeris(tr, fr, a.period, a.t0, a.dur, cdr)
            mid = 0.5 * (tr.min() + tr.max())
            gap = tr[np.argmax(np.diff(tr)) + 1] if np.diff(tr).max() > 0.5 else mid
            res["control_orbits_at_ref"] = {name: at_ephemeris(tr[m], fr[m], a.period, a.t0, a.dur, cdr)
                                            for name, m in (("orbit1", tr < gap), ("orbit2", tr >= gap))}
        (out_dir / "ffi_check.json").write_text(json.dumps(res, indent=1, default=float))
        print(json.dumps({k: v for k, v in res.items() if k != "sectors"}, indent=1, default=float))
        return
    tf = np.concatenate([data[s]["flat"][0] for s in ffi]); ff = np.concatenate([data[s]["flat"][1] for s in ffi])
    tb = np.concatenate([data[s]["binned"][0] for s in ffi]); fb = np.concatenate([data[s]["binned"][1] for s in ffi])
    blind = bls.search_s(tb, fb, tf, ff)
    per, pw = blind["periods"], blind["power"]
    near = np.abs(per / a.period - 1) < 0.005
    res["stitched_blind"] = dict(top=dict(period=blind["period"], t0=blind["t0"], duration=blind["duration"], S=blind["S"], sde=blind["sde"],
                                          depth_ppm=1e6 * blind["depth_box"]),
                                 candidate_period_rank=int((pw > pw[near].max()).sum()) + 1 if near.any() else None,
                                 candidate_period_peak_S=float(pw[near].max()) if near.any() else None,
                                 best_period_near_candidate=float(per[near][np.argmax(pw[near])]) if near.any() else None)
    D = a.dur
    cd = blind["cdpp_by_d"][f"{min(bls.DURATIONS, key=lambda x: abs(x - D)):.2f}"]
    # FFI-only local ephemeris near the candidate period
    Pn = res["stitched_blind"]["best_period_near_candidate"] or a.period
    mid = 0.5 * (tf.min() + tf.max())
    # start T0 near the middle of the FFI data from the BLS phase at Pn
    r0 = bls.BoxLeastSquares(tb, fb).power([Pn], [D], objective="snr")
    t0n = float(r0.transit_time[0]); t0n += np.round((mid - t0n) / Pn) * Pn
    dP, sP = scan_width(tf, ff, cd, Pn, t0n, D, "P", 0.004)
    Pf = Pn + dP
    dT, sT = scan_width(tf, ff, cd, Pf, t0n, D, "t0", 0.3 * Pf)
    T0f = t0n + dT
    res["ffi_ephemeris"] = dict(period=Pf, period_err=sP, t0=T0f, t0_err=sT)
    res["ffi_at_ffi_ephemeris"] = at_ephemeris(tf, ff, Pf, T0f, D, cd)
    # reference ephemeris propagated back
    n_cyc = (a.t0 - T0f) / Pf
    res["link"] = dict(cycles_ffi_to_ref=float(n_cyc))
    if a.control:
        tr, fr = data[a.control]["flat"][0], data[a.control]["flat"][1]
        cdr = None
        rb = bls.search_s(*data[a.control]["binned"][:2], tr, fr)
        cdr = rb["cdpp_by_d"][f"{min(bls.DURATIONS, key=lambda x: abs(x - D)):.2f}"]
        dPr, sPr = scan_width(tr, fr, cdr, a.period, a.t0, D, "P", 0.004)
        dTr, sTr = scan_width(tr, fr, cdr, a.period + dPr, a.t0, D, "t0", 0.2)
        res["ref_ephemeris_from_control"] = dict(period=a.period + dPr, period_err=sPr, t0=a.t0 + dTr, t0_err=sTr)
        res["control_at_ref"] = at_ephemeris(tr, fr, a.period, a.t0, D, cdr)
        Pr, ePr, Tr, eTr = a.period + dPr, sPr, a.t0 + dTr, sTr
        # integer-cycle linkage: joint period from T0 difference must agree with both sector periods
        N = np.round((Tr - T0f) / Pr)
        Pj = (Tr - T0f) / N
        ePj = np.hypot(sT, eTr) / abs(N)
        alias = abs(Pr - Pf) / np.hypot(ePr, sP)
        # phase of the FFI events under the reference ephemeris (propagated with the reference period error)
        dphase_days = ((T0f - Tr + 0.5 * Pr) % Pr) - 0.5 * Pr
        prop_err = np.hypot(abs(N) * ePr, np.hypot(sT, eTr))
        res["link"].update(N_cycles=int(N), joint_period=float(Pj), joint_period_err=float(ePj),
                           period_agreement_sigma=float(alias),
                           ffi_minus_ref_phase_days=float(dphase_days), ref_ephemeris_propagated_err_days=float(prop_err),
                           phase_consistent_within_2sigma=bool(abs(dphase_days) < 2 * prop_err),
                           propagated_err_in_durations=float(prop_err / D))
        res["ffi_at_ref_ephemeris"] = at_ephemeris(tf, ff, a.period, a.t0, D, cd)
    (out_dir / "ffi_check.json").write_text(json.dumps(res, indent=1, default=float))
    print(json.dumps({k: v for k, v in res.items() if k != "sectors"}, indent=1, default=float))


if __name__ == "__main__":
    main()
