"""Stellar inputs shared by the injection, density and reliability code.

* quadratic limb darkening for the TESS band: Claret (2017, A&A 600, A30), PHOENIX models,
  quasi-spherical r-method, solar metallicity (the only metallicity in the VizieR table)
* TIC v8.2 stellar parameters (MAST Catalogs, cached per star)
* Gaia DR3 RUWE / non_single_star for the TIC's Gaia source (ESA TAP, cached per star)
All downloads are cached under catalogs/ so reruns are offline and reproducible.
"""
from __future__ import annotations

import io
import os
import time

import numpy as np
import pandas as pd
import requests

from . import config

CLARET_SRC = "J/A+A/600/A30/tableab"
CLARET_URLS = ["https://vizier.cfa.harvard.edu/viz-bin/asu-tsv", "https://vizier.cds.unistra.fr/viz-bin/asu-tsv"]
GAIA_TAP = "https://gea.esac.esa.int/tap-server/tap/sync"
G_CGS, MSUN_G, RSUN_CM = 6.674e-8, 1.989e33, 6.957e10
RHO_SUN = MSUN_G / (4 / 3 * np.pi * RSUN_CM ** 3)          # 1.41 g/cc
TIC_COLS = ["ID", "Tmag", "Teff", "e_Teff", "logg", "e_logg", "rad", "e_rad", "mass", "e_mass",
            "rho", "e_rho", "GAIA", "lumclass", "contratio", "ra", "dec"]


# ---------------------------------------------------------------- limb darkening
_LD = None


def claret_table():
    path = config.CATALOG_DIR / "claret2017_tess_quadratic.csv"
    if not path.exists():
        last = None
        for base in CLARET_URLS:
            try:
                r = requests.get(f"{base}?-source={CLARET_SRC.replace('+', '%2B')}&-out.all&-out.max=999999", timeout=120)
                r.raise_for_status()
                lines = [l for l in r.text.splitlines() if l and not l.startswith("#")]
                hdr = next(i for i, l in enumerate(lines) if l.startswith("recno"))
                df = pd.read_csv(io.StringIO("\n".join([lines[hdr]] + lines[hdr + 3:])), sep="\t")
                df.columns = [c.strip() for c in df.columns]
                df = df[["logg", "Teff", "Z", "aLSM", "bLSM", "Type", "Mod"]]
                df["Type"], df["Mod"] = df.Type.str.strip(), df.Mod.str.strip()
                df.to_csv(path, index=False)
                break
            except Exception as exc:     # mirror down -> try the next one
                last = exc
        else:
            raise RuntimeError(f"Claret table unavailable: {last}")
    return pd.read_csv(path)


def quad_ld(teff, logg):
    """(u1, u2) at the nearest grid point in (Teff, logg); PHOENIX-COND, r-method; PHOENIX-DRIFT below its range."""
    global _LD
    if _LD is None:
        t = claret_table()
        t = t[t.Type == "r"]
        pc = t[t.Mod == "PC"]
        pd_ = t[(t.Mod == "PD") & (t.Teff < pc.Teff.min())]
        _LD = pd.concat([pc, pd_])[["Teff", "logg", "aLSM", "bLSM"]].to_numpy(float)
    teff = 5800.0 if teff is None or not np.isfinite(teff) else float(teff)
    logg = 4.4 if logg is None or not np.isfinite(logg) else float(logg)
    d = ((_LD[:, 0] - teff) / 250.0) ** 2 + ((_LD[:, 1] - logg) / 0.5) ** 2
    j = int(np.argmin(d))
    return float(_LD[j, 2]), float(_LD[j, 3])


# ---------------------------------------------------------------- TIC v8.2
def tic_params(tics, chunk=400, log=print):
    """TIC v8.2 rows for `tics` (cached in catalogs/tic82_params.parquet)."""
    from astroquery.mast import Catalogs
    path = config.CATALOG_DIR / "tic82_params.parquet"
    have = pd.read_parquet(path) if path.exists() else None
    known = set() if have is None else set(have.ID.astype(int))
    want = sorted(set(int(t) for t in tics) - known)
    new = []
    for i in range(0, len(want), chunk):
        ids = [str(t) for t in want[i:i + chunk]]
        for attempt in range(4):
            try:
                r = Catalogs.query_criteria(catalog="Tic", ID=ids).to_pandas()
                break
            except Exception as exc:
                log(f"  TIC query retry {attempt}: {exc}")
                time.sleep(10 * (attempt + 1))
        else:
            raise RuntimeError("TIC query failed")
        r = r[[c for c in TIC_COLS if c in r.columns]]
        r["ID"] = r.ID.astype(int)
        new.append(r)
    if new:
        have = pd.concat(([] if have is None else [have]) + new, ignore_index=True).drop_duplicates("ID")
        for c in ("GAIA", "lumclass"):
            have[c] = have[c].astype(str)
        have.to_parquet(path, index=False)
    return have[have.ID.isin([int(t) for t in tics])].set_index("ID")


def tic_rho(row):
    """(rho, sigma_rho) in solar units from TIC; derived from M and R when rho is missing."""
    rho, erho = row.get("rho"), row.get("e_rho")
    if rho is not None and np.isfinite(rho) and rho > 0:
        return float(rho), float(erho) if erho is not None and np.isfinite(erho) and erho > 0 else 0.3 * float(rho)
    m, r = row.get("mass"), row.get("rad")
    if m and r and np.isfinite(m) and np.isfinite(r) and r > 0:
        em = row.get("e_mass") if row.get("e_mass") and np.isfinite(row.get("e_mass")) else 0.15 * m
        er = row.get("e_rad") if row.get("e_rad") and np.isfinite(row.get("e_rad")) else 0.10 * r
        rho = m / r ** 3
        return float(rho), float(rho * np.hypot(em / m, 3 * er / r))
    return None, None


# ---------------------------------------------------------------- Gaia DR3
def gaia_dr3(dr2_ids, chunk=300, log=print, strict=False):
    """RUWE and non_single_star for each DR2 id via gaiadr3.dr2_neighbourhood (cached)."""
    path = config.CATALOG_DIR / "gaia_dr3_ruwe.parquet"
    cols = ["dr2_source_id", "source_id", "ruwe", "non_single_star", "phot_g_mean_mag", "angular_distance"]
    have = pd.read_parquet(path) if path.exists() else None
    known = set() if have is None else set(have.dr2_source_id.astype("int64"))
    ids = sorted({int(x) for x in dr2_ids if str(x).isdigit()} - known)
    if os.environ.get("EXO_GAIA_OFFLINE"):          # cache only (ESA TAP outages)
        ids = []
    new = []
    for i in range(0, len(ids), chunk):
        sub = ",".join(str(x) for x in ids[i:i + chunk])
        q = ("SELECT n.dr2_source_id, g.source_id, g.ruwe, g.non_single_star, g.phot_g_mean_mag, n.angular_distance "
             "FROM gaiadr3.dr2_neighbourhood AS n JOIN gaiadr3.gaia_source AS g ON g.source_id = n.dr3_source_id "
             f"WHERE n.dr2_source_id IN ({sub})")
        for attempt in range(4 if strict else 2):
            try:
                r = requests.post(GAIA_TAP, data=dict(REQUEST="doQuery", LANG="ADQL", FORMAT="csv", QUERY=q), timeout=90)
                r.raise_for_status()
                df = pd.read_csv(io.StringIO(r.text))
                break
            except Exception as exc:
                log(f"  Gaia query retry {attempt}: {exc}")
                time.sleep(10 * (attempt + 1))
        else:
            if strict:
                raise RuntimeError("Gaia TAP failed")
            log("  Gaia TAP unavailable: RUWE/NSS left empty for the remaining stars (warning-only checks)")
            break
        df = df.sort_values("angular_distance").drop_duplicates("dr2_source_id")   # closest DR3 match
        miss = sorted(set(ids[i:i + chunk]) - set(df.dr2_source_id.astype("int64")))
        if miss:
            df = pd.concat([df, pd.DataFrame({"dr2_source_id": np.array(miss, dtype="int64")})], ignore_index=True)
        new.append(df)
    if new:
        have = pd.concat(([] if have is None else [have]) + new, ignore_index=True).drop_duplicates("dr2_source_id")
        have["dr2_source_id"] = have.dr2_source_id.astype("int64")
        have["source_id"] = have.source_id.astype("Int64")
        have.to_parquet(path, index=False)
    if have is None:
        return pd.DataFrame(columns=cols).set_index("dr2_source_id")
    return have.set_index("dr2_source_id")
