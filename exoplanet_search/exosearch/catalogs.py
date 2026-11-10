"""Known-object cross-match: ExoFOP TOIs + CTOIs and NASA Exoplanet Archive (pscomppars)."""
from __future__ import annotations

import io
import os
import time
import warnings

import numpy as np
import pandas as pd
import requests

from . import config

TOI_URL = "https://exofop.ipac.caltech.edu/tess/download_toi.php?sort=toi&output=csv"
CTOI_URL = "https://exofop.ipac.caltech.edu/tess/download_ctoi.php?sort=ctoi&output=csv"
PS_URL = ("https://exoplanetarchive.ipac.caltech.edu/TAP/sync?query="
          "select+pl_name,tic_id,pl_orbper,disc_facility+from+pscomppars+where+tic_id+is+not+null"
          "&format=csv")


def _cached_csv(name, url, max_age_days=7):
    config.CATALOG_DIR.mkdir(parents=True, exist_ok=True)
    path = config.CATALOG_DIR / name
    fresh = path.exists() and (time.time() - path.stat().st_mtime) < max_age_days * 86400
    if not fresh:
        try:
            r = requests.get(url, timeout=180)
            r.raise_for_status()
            tmp = path.with_suffix(".part")
            tmp.write_bytes(r.content)
            os.replace(tmp, path)
        except Exception:
            if not path.exists():
                raise
    return pd.read_csv(path, comment="#", low_memory=False)


def load_catalogs():
    toi = _cached_csv("exofop_toi.csv", TOI_URL)
    toi = toi.rename(columns={"TIC ID": "tic", "Period (days)": "period",
                              "TFOPWG Disposition": "disp", "TOI": "toi"})[["tic", "toi", "period", "disp"]]
    try:
        ctoi = _cached_csv("exofop_ctoi.csv", CTOI_URL)
        ctoi = ctoi.rename(columns={"TIC ID": "tic", "Period (days)": "period", "CTOI": "ctoi"})
        ctoi = ctoi[["tic", "ctoi", "period"]]
    except Exception:
        ctoi = pd.DataFrame(columns=["tic", "ctoi", "period"])
    try:
        ps = _cached_csv("nasa_pscomppars.csv", PS_URL)
    except Exception as exc:  # archive down for maintenance etc.: TOI KP/CP dispositions still cover most
        warnings.warn(f"NASA Exoplanet Archive unavailable ({exc}); planet cross-match skipped")
        return dict(toi=toi, ctoi=ctoi, ps=pd.DataFrame(columns=["pl_name", "tic", "pl_orbper"]), ps_ok=False)
    ps["tic"] = pd.to_numeric(ps["tic_id"].astype(str).str.replace("TIC", "").str.strip(), errors="coerce")
    ps = ps.dropna(subset=["tic"]).astype({"tic": "int64"})
    return dict(toi=toi, ctoi=ctoi, ps=ps, ps_ok=True)


def _period_match(p, refs, tol=0.01):
    refs = np.asarray([r for r in refs if r == r and r > 0], float)
    if not len(refs) or not p:
        return False
    for k in (1.0, 2.0, 0.5, 3.0, 1 / 3.0):
        if np.any(np.abs(p * k / refs - 1.0) < tol):
            return True
    return False


def crossmatch(tic, period, cats):
    t = cats["toi"][cats["toi"].tic == tic]
    c = cats["ctoi"][cats["ctoi"].tic == tic]
    p = cats["ps"][cats["ps"].tic == tic]
    return dict(
        known_toi=";".join(str(x) for x in t.toi) or None,
        toi_disp=";".join(str(x) for x in t.disp.fillna("")) or None,
        toi_period_match=int(_period_match(period, t.period)),
        known_ctoi=";".join(str(x) for x in c.ctoi) or None,
        ctoi_period_match=int(_period_match(period, c.period)),
        known_planet=";".join(p.pl_name) or None,
        planet_period_match=int(_period_match(period, p.pl_orbper)),
    )
