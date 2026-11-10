"""Deeper vetting, run only on candidates (stars that passed the per-star screens).

Eclipsing-binary checks
  * cross-match the TESS EB catalogue (Prsa et al. 2022, VizieR J/ApJS/258/16) by TIC and by position
    (an EB on a neighbour with the same period is a blend source), and the Kepler EB catalogue
    (Kirk et al. 2016, J/AJ/151/68) by position where the fields overlap;
  * transit shape: depth in the central third vs. the outer two thirds (V-shape vs flat bottom),
    on a finer duration grid at the BLS period;
  * secondary-eclipse scan over all phases 0.15-0.85, not only 0.5 (eccentric binaries);
  * radius sanity: implied companion radius from the TIC stellar radius, flag > 2 R_Jup.

Blend checks
  * difference-image centroid test on the SPOC target-pixel file: average (out-of-transit minus
    in-transit) image over all transits, locate where the light was lost, and compare with the
    out-of-transit centroid.  If the dip happens on a neighbouring pixel, it isn't our star.
  * Gaia DR3 neighbours within ~2 TESS pixels (42") that are bright enough to produce the
    observed depth even if they were completely eclipsed.

    python -m exosearch.deepvet --sector 70 [--only-review] [--redo]
"""
from __future__ import annotations

import argparse
import datetime as dt
import io
import json
import warnings

import numpy as np
import pandas as pd
import requests

from . import config, store, targets

TESS_PIX_ARCSEC = 21.0
GAIA_TAP = "https://gea.esac.esa.int/tap-server/tap/sync"
MAST_FILE = "https://mast.stsci.edu/api/v0.1/Download/file/?uri=mast:TESS/product/"

COLS = {
    "tic": "INTEGER", "sector": "INTEGER", "n_epochs": "INTEGER", "diff_peak_snr": "REAL",
    "centroid_offset_px": "REAL", "centroid_offset_sigma": "REAL", "diff_peak_dist_px": "REAL",
    "gaia_n_neighbours": "INTEGER", "gaia_n_can_mimic": "INTEGER", "gaia_mimic_min_sep": "REAL",
    "gaia_mimic_brightest_dm": "REAL", "eb_match": "TEXT", "v_ratio": "REAL", "v_sigma": "REAL",
    "duration_fine": "REAL", "sec_scan_sigma": "REAL", "sec_scan_phase": "REAL", "rp_rjup": "REAL",
    "deep_flags": "TEXT", "deep_notes": "TEXT", "updated_at": "TEXT",
}
R_JUP_RSUN = 0.10045
VIZIER = ["https://vizier.cds.unistra.fr/viz-bin/asu-tsv", "https://vizier.cfa.harvard.edu/viz-bin/asu-tsv"]
THRESH = dict(v_ratio=0.45, v_sigma=5.0, sec_scan_sigma=4.0, max_rp_rjup=2.0, eb_match_arcsec=42.0,
              centroid_px=0.5, centroid_sigma=3.0, diff_dist_px=1.5, diff_snr=3.0, gaia_radius_px=2.0,
              unresolved_arcsec=10.0)


def connect():
    con = store.connect()
    cols = ", ".join(f"{k} {v}" for k, v in COLS.items())
    con.execute(f"CREATE TABLE IF NOT EXISTS deepvet ({cols}, PRIMARY KEY (tic, sector))")
    have = {r[1] for r in con.execute("PRAGMA table_info(deepvet)")}
    for k, v in COLS.items():          # forward-compatible schema
        if k not in have:
            con.execute(f"ALTER TABLE deepvet ADD COLUMN {k} {v}")
    return con


def upsert(con, row):
    row = {k: row.get(k) for k in COLS}
    con.execute(f"INSERT OR REPLACE INTO deepvet ({', '.join(row)}) VALUES ({', '.join('?' for _ in row)})",
                list(row.values()))
    con.commit()


def frame(con, sector=None):
    q = "SELECT * FROM deepvet" + (" WHERE sector=?" if sector is not None else "")
    return pd.read_sql_query(q, con, params=(sector,) if sector is not None else None)


# ---------------------------------------------------------------- difference-image centroids
def fetch_tpf(sector, lc_fname):
    fname = lc_fname.replace("_lc.fits", "_tp.fits")
    return targets.download_one(fname, MAST_FILE + fname, config.sector_dir(sector) / "tp")


def _centroid3(img, iy, ix):
    ny, nx = img.shape
    y0, y1, x0, x1 = max(iy - 1, 0), min(iy + 2, ny), max(ix - 1, 0), min(ix + 2, nx)
    sub = np.clip(img[y0:y1, x0:x1], 0, None)
    if sub.sum() <= 0:
        return float(iy), float(ix)
    yy, xx = np.mgrid[y0:y1, x0:x1]
    return float((sub * yy).sum() / sub.sum()), float((sub * xx).sum() / sub.sum())


def centroid_core(t, flux, period, t0, duration, target_yx=None, n_boot=200, seed=0):
    """Difference-image centroid offset from cadence-by-pixel arrays (flux: [n, ny, nx])."""
    n0, n1 = np.ceil((t.min() - t0) / period), np.floor((t.max() - t0) / period)
    diffs = []
    for n in np.arange(n0, n1 + 1):
        tc = t0 + n * period
        dtc = np.abs(t - tc)
        ins = dtc < 0.4 * duration
        out = (dtc > 0.75 * duration) & (dtc < 1.75 * duration)
        if ins.sum() >= 3 and out.sum() >= 6:
            diffs.append(np.nanmean(flux[out], axis=0) - np.nanmean(flux[ins], axis=0))
    res = dict(n_epochs=len(diffs))
    if len(diffs) < 2:
        res["deep_notes"] = "centroid: <2 usable transits in TPF"
        return res
    D = np.array(diffs)
    oot = np.nanmedian(flux, axis=0)
    mean = np.nanmean(D, axis=0)
    mad = 1.4826 * np.nanmedian(np.abs(mean - np.nanmedian(mean)))
    iy, ix = np.unravel_index(np.nanargmax(np.nan_to_num(mean, nan=-np.inf)), mean.shape)
    res["diff_peak_snr"] = float(mean[iy, ix] / mad) if mad > 0 else None
    if target_yx is None or not np.all(np.isfinite(target_yx)):
        target_yx = np.unravel_index(np.nanargmax(np.nan_to_num(oot)), oot.shape)
    # reference = out-of-transit light centroid around the target (absorbs WCS errors)
    oy, ox = _centroid3(np.nan_to_num(oot), int(round(target_yx[0])), int(round(target_yx[1])))
    dy, dx = _centroid3(np.nan_to_num(mean), iy, ix)
    rng = np.random.default_rng(seed)
    boots = []
    for _ in range(n_boot):
        m = np.nanmean(D[rng.integers(0, len(D), len(D))], axis=0)
        by, bx = np.unravel_index(np.nanargmax(np.nan_to_num(m, nan=-np.inf)), m.shape)
        cy, cx = _centroid3(np.nan_to_num(m), by, bx)
        boots.append(np.hypot(cy - oy, cx - ox))
    off = float(np.hypot(dy - oy, dx - ox))
    sig = float(np.std(boots)) or 1e-3
    res.update(centroid_offset_px=off, centroid_offset_sigma=off / sig,
               diff_peak_dist_px=float(np.hypot(iy - oy, ix - ox)))
    return res


def centroid_test(tpf_path, period, t0, duration):
    import lightkurve as lk
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        tpf = lk.read(tpf_path, quality_bitmask="hard")
    t = np.asarray(tpf.time.value, float)
    flux = np.asarray(tpf.flux.value, float)
    ok = np.isfinite(t)
    try:
        px = tpf.wcs.world_to_pixel_values(tpf.ra, tpf.dec)
        tyx = (float(px[1]), float(px[0]))
    except Exception:
        tyx = None
    return centroid_core(t[ok], flux[ok], period, t0, duration, tyx)


# ---------------------------------------------------------------- Gaia neighbours
def _gaia_mast(ra, dec, radius_arcsec):
    import astropy.units as u
    from astropy.coordinates import SkyCoord
    from astroquery.mast import Catalogs
    t = Catalogs.query_region(SkyCoord(ra, dec, unit="deg"), radius=radius_arcsec * u.arcsec,
                              catalog="Gaia", version=3)
    return t.to_pandas()[["source_id", "ra", "dec", "phot_g_mean_mag"]]


def _gaia_esa(ra, dec, radius_arcsec):
    q = (f"SELECT source_id, ra, dec, phot_g_mean_mag FROM gaiadr3.gaia_source WHERE "
         f"1=CONTAINS(POINT('ICRS',ra,dec), CIRCLE('ICRS',{ra},{dec},{radius_arcsec/3600.0}))")
    r = requests.get(GAIA_TAP, params=dict(REQUEST="doQuery", LANG="ADQL", FORMAT="csv", QUERY=q), timeout=120)
    r.raise_for_status()
    return pd.read_csv(io.StringIO(r.text))


def gaia_neighbours(ra, dec, radius_arcsec):
    """Gaia DR3 cone (MAST mirror first, ESA archive as fallback), cached per position."""
    cache = config.CATALOG_DIR / "gaia_cone" / f"{ra:.5f}_{dec:.5f}_{radius_arcsec:.0f}.csv"
    if cache.exists():
        return pd.read_csv(cache)
    try:
        df = _gaia_mast(ra, dec, radius_arcsec)
    except Exception:
        df = _gaia_esa(ra, dec, radius_arcsec)
    df = df.dropna(subset=["phot_g_mean_mag"]).copy()
    cosd = np.cos(np.radians(dec))
    df["sep"] = 3600 * np.hypot((df.ra - ra) * cosd, df.dec - dec)
    cache.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(cache, index=False)
    return df


def max_dm_to_mimic(depth):
    """A neighbour this many mag fainter, if 100% eclipsed, dilutes to exactly `depth`."""
    depth = min(max(depth, 1e-7), 0.99)
    return float(-2.5 * np.log10(depth / (1 - depth)))


def gaia_test(ra, dec, tmag, depth):
    g = gaia_neighbours(ra, dec, THRESH["gaia_radius_px"] * TESS_PIX_ARCSEC)
    if g.empty:
        return dict(gaia_n_neighbours=0, gaia_n_can_mimic=0)
    # the target: brightest source within 4"
    near = g[g.sep < 4.0]
    tgt = near.sort_values("phot_g_mean_mag").iloc[0] if len(near) else g.sort_values("sep").iloc[0]
    nb = g[g.source_id != tgt.source_id].copy()
    nb["dm"] = nb.phot_g_mean_mag - tgt.phot_g_mean_mag
    mimic = nb[nb.dm <= max_dm_to_mimic(depth)]
    return dict(gaia_n_neighbours=int(len(nb)), gaia_n_can_mimic=int(len(mimic)),
                gaia_mimic_min_sep=float(mimic.sep.min()) if len(mimic) else None,
                gaia_mimic_brightest_dm=float(mimic.dm.min()) if len(mimic) else None)


# ---------------------------------------------------------------- eclipsing binaries
def _vizier(source, cols, name):
    from . import catalogs
    last = None
    for base in VIZIER:
        url = f"{base}?-source={source}&-out={','.join(cols)}&-out.max=unlimited"
        try:
            path = config.CATALOG_DIR / name
            meta = config.CATALOG_DIR / (name + ".meta.json")
            if path.exists() and meta.exists():
                return catalogs._cached_csv(name, url, max_age_days=90)
            r = requests.get(url, timeout=300)
            r.raise_for_status()
            lines = [ln for ln in r.text.splitlines() if not ln.startswith("#")]
            while lines and not lines[0].strip():
                lines.pop(0)
            hdr, body = lines[0], []
            for ln in lines[3:]:                 # first table only: stop at the blank separator
                if not ln.strip():
                    break
                if not ln.startswith("-"):
                    body.append(ln)
            df = pd.read_csv(io.StringIO("\n".join([hdr] + body)), sep="\t")
            if df.empty:
                raise IOError("empty")
            config.CATALOG_DIR.mkdir(parents=True, exist_ok=True)
            df.to_csv(path, index=False)
            meta.write_text(json.dumps(dict(fetched_at=dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
                                            rows=int(len(df)), source=base + f"?-source={source}", last_error=None)))
            return df
        except Exception as exc:
            last = exc
    raise RuntimeError(f"VizieR {source} unavailable: {last}")


_EB = {}


def eb_catalogs():
    if not _EB:
        try:
            t = _vizier("J/ApJS/258/16/tess-ebs", ["TIC", "RAJ2000", "DEJ2000", "Per", "Morph"], "tess_ebs_prsa2022.csv")
            _EB["tess"] = t.rename(columns={"RAJ2000": "ra", "DEJ2000": "dec"})
        except Exception as exc:
            warnings.warn(str(exc))
            _EB["tess"] = pd.DataFrame(columns=["TIC", "ra", "dec", "Per", "Morph"])
        try:
            k = _vizier("J/AJ/151/68", ["KIC", "Per", "_RA", "_DE", "Morph"], "kepler_ebs_kirk2016.csv")
            _EB["kepler"] = k.rename(columns={"_RA": "ra", "_DE": "dec"})
        except Exception as exc:
            warnings.warn(str(exc))
            _EB["kepler"] = pd.DataFrame(columns=["KIC", "ra", "dec", "Per", "Morph"])
    return _EB


def _harmonic(p, q, tol=0.01):
    return any(abs(p * k / q - 1) < tol for k in (1, 2, 0.5, 3, 1 / 3))


def eb_match(tic, ra, dec, period):
    out = []
    for name, cat in eb_catalogs().items():
        if cat.empty:
            continue
        idcol = "TIC" if name == "tess" else "KIC"
        if name == "tess":
            same = cat[pd.to_numeric(cat.TIC, errors="coerce") == tic]
            for r in same.itertuples():
                out.append(f"TESS-EB TIC {tic} P={r.Per}")
        if ra == ra and dec == dec:
            cosd = np.cos(np.radians(dec))
            sep = 3600 * np.hypot((pd.to_numeric(cat.ra, errors="coerce") - ra) * cosd,
                                  pd.to_numeric(cat.dec, errors="coerce") - dec)
            near = cat[(sep < THRESH["eb_match_arcsec"]) & (sep > (1.0 if name == "tess" else -1))]
            for r, sp in zip(near.itertuples(), sep[near.index]):
                per = pd.to_numeric(getattr(r, "Per"), errors="coerce")
                if name == "tess" and int(getattr(r, idcol)) == tic:
                    continue
                if per == per and _harmonic(period, per):
                    out.append(f"{'TESS' if name == 'tess' else 'Kepler'}-EB {idcol} {getattr(r, idcol)} "
                               f"{sp:.0f}\" away, P={per}")
    return "; ".join(out) or None


def shape_and_secondary(row, lc_path):
    from astropy.timeseries import BoxLeastSquares
    from .task import prepare
    d = prepare(lc_path)
    t, f, e = d["flat"]
    P, t0 = float(row.period), float(row.t0)
    durs = np.linspace(0.02, min(0.3, 0.2 * P), 40)
    bls = BoxLeastSquares(t, f, dy=e)
    pw = bls.power([P], durs, objective="likelihood")
    j = int(np.nanargmax(pw.power))
    dur, t0 = float(pw.duration[j]), float(pw.transit_time[j])
    x = (((t - t0 + 0.5 * P) % P) - 0.5 * P) / dur          # in units of the duration
    c = np.abs(x) < 1 / 6
    o = (np.abs(x) >= 1 / 6) & (np.abs(x) < 0.5)
    oot = (np.abs(x) > 1.0) & (np.abs(x) < 3.0)
    base = np.median(f[oot]) if oot.sum() > 10 else 1.0
    dc, do = base - np.mean(f[c]), base - np.mean(f[o])
    sc = np.std(f[oot]) / np.sqrt(max(c.sum(), 1)) if oot.sum() > 10 else np.nan
    res = dict(duration_fine=dur, v_ratio=float(do / dc) if dc > 0 else None,
               v_sigma=float(dc / sc) if sc and sc > 0 else None)
    # secondary scan: box of the same duration at every phase away from the primary
    ph = ((t - t0) / P) % 1.0
    sd = np.std(f[oot]) if oot.sum() > 10 else np.std(f)
    best = (0.0, None)
    w = dur / P
    for p0 in np.arange(0.15, 0.85, max(w / 3, 0.002)):
        m = np.abs(ph - p0) < w / 2
        if m.sum() < 5:
            continue
        sig = (base - np.mean(f[m])) / (sd / np.sqrt(m.sum()))
        if sig > best[0]:
            best = (float(sig), float(p0))
    res.update(sec_scan_sigma=best[0], sec_scan_phase=best[1])
    return res


# ---------------------------------------------------------------- driver
def flags_for(r):
    f = []
    if r.get("diff_peak_snr") is not None and r["diff_peak_snr"] >= THRESH["diff_snr"]:
        if (r.get("centroid_offset_px") or 0) > THRESH["centroid_px"] and \
                (r.get("centroid_offset_sigma") or 0) > THRESH["centroid_sigma"]:
            f.append("centroid_offset")
        if (r.get("diff_peak_dist_px") or 0) > THRESH["diff_dist_px"]:
            f.append("dip_off_target")
    elif r.get("n_epochs", 0) >= 2:
        f.append("centroid_inconclusive")
    if r.get("eb_match"):
        f.append("known_eb")
    if r.get("v_ratio") is not None and (r.get("v_sigma") or 0) >= THRESH["v_sigma"] and r["v_ratio"] < THRESH["v_ratio"]:
        f.append("v_shape")
    if (r.get("sec_scan_sigma") or 0) > THRESH["sec_scan_sigma"]:
        f.append("secondary_scan")
    if (r.get("rp_rjup") or 0) > THRESH["max_rp_rjup"]:
        f.append("rp_gt_2rjup")
    if r.get("gaia_n_can_mimic"):
        if (r.get("gaia_mimic_min_sep") or 99) < THRESH["unresolved_arcsec"]:
            f.append("unresolved_neighbour")      # too close for the centroid test to rule out
        elif "centroid_inconclusive" in f or r.get("diff_peak_snr") is None:
            f.append("neighbour_could_mimic")
    return f


def vet_star(row, fname):
    out = dict(tic=int(row.tic), sector=int(row.sector),
               updated_at=dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"))
    notes = []
    try:
        out.update(centroid_test(fetch_tpf(int(row.sector), fname), row.period, row.t0, row.duration))
    except Exception as exc:
        notes.append(f"centroid failed: {type(exc).__name__}: {exc}"[:160])
    try:
        if row.ra == row.ra and row.dec == row.dec:
            out.update(gaia_test(row.ra, row.dec, row.tmag, row.depth))
    except Exception as exc:
        notes.append(f"gaia failed: {type(exc).__name__}: {exc}"[:160])
    try:
        out["eb_match"] = eb_match(int(row.tic), row.ra, row.dec, float(row.period))
    except Exception as exc:
        notes.append(f"eb match failed: {exc}"[:160])
    try:
        lc = targets.download_one(fname, MAST_FILE + fname, config.lc_dir(int(row.sector)))
        out.update(shape_and_secondary(row, lc))
    except Exception as exc:
        notes.append(f"shape failed: {type(exc).__name__}: {exc}"[:160])
    if row.radius == row.radius and row.radius and row.depth > 0:
        out["rp_rjup"] = float(np.sqrt(row.depth) * row.radius / R_JUP_RSUN)
    else:
        notes.append("no TIC stellar radius")
    if out.get("deep_notes"):
        notes.append(out["deep_notes"])
    out["deep_flags"] = ",".join(flags_for(out))
    out["deep_notes"] = "; ".join(notes) or None
    return out


def candidates(df, only_review=False):
    c = df[(df.status == "ok") & (df.passed == 1)]
    if only_review:
        known = (c.toi_period_match == 1) | (c.ctoi_period_match == 1) | (c.planet_period_match == 1)
        c = c[~known]
    return c


def run(sector, only_review=False, redo=False, prune=True, log=print):
    con = connect()
    df = store.frame(con, sector)
    cand = candidates(df, only_review)
    done = set() if redo else set(frame(con, sector).tic)
    man = targets.fetch_manifest(sector).set_index("tic")
    todo = cand[~cand.tic.isin(done)]
    log(f"deep vetting sector {sector}: {len(cand)} candidates, {len(todo)} to do")
    for i, r in enumerate(todo.itertuples(), 1):
        if r.tic not in man.index:
            continue
        upsert(con, vet_star(r, man.loc[r.tic].fname))
        if i % 10 == 0:
            log(f"  {i}/{len(todo)}")
    if prune:
        for p in (config.sector_dir(sector) / "tp").glob("*.fits"):
            p.unlink()
    return frame(con, sector)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sector", type=int, required=True)
    ap.add_argument("--only-review", action="store_true", help="skip candidates that match known objects")
    ap.add_argument("--redo", action="store_true")
    ap.add_argument("--keep-tpf", action="store_true")
    a = ap.parse_args()
    d = run(a.sector, a.only_review, a.redo, prune=not a.keep_tpf)
    print(d.deep_flags.fillna("").str.split(",").explode().value_counts().to_string())


if __name__ == "__main__":
    main()
