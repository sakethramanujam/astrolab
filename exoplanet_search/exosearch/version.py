"""Frozen pipeline identity: version string + hash of every setting that can change a result.

Every reliability product (injections, null tests, labelled-set checks) stores PIPELINE_VERSION,
config_hash() and git_hash() in each row. Change any setting below -> new hash -> new version.
"""
from __future__ import annotations

import hashlib
import inspect
import json
import platform
import subprocess
from importlib import metadata
from pathlib import Path

PIPELINE_VERSION = "1.0.0"
REPO = Path(__file__).resolve().parents[2]


def _defaults(fn):
    return {k: v.default for k, v in inspect.signature(fn).parameters.items()
            if v.default is not inspect.Parameter.empty and not callable(v.default)}


def config():
    from . import bls, deepvet, lightcurve, task, vet
    return {
        "quality_mask": "QUALITY == 0 and finite TIME/PDCSAP_FLUX/PDCSAP_FLUX_ERR, ERR > 0; flux / median",
        "detrend": _defaults(lightcurve.detrend) | {"method": "biweight", "edge_cutoff": 0.0, "break_tolerance": 0.5},
        "gap_edge": _defaults(lightcurve.gap_edge_mask),
        "prepare": _defaults(task.prepare_arrays),
        "bls": {"durations_d": [float(x) for x in bls.DURATIONS], "grid": _defaults(bls.period_grid),
                "statistic": "S = depth_box sqrt(N_tr) / CDPP_D; astropy BLS 'snr' objective per D with "
                             "sigma_D = CDPP_D sqrt(D/bin); N_tr coverage-weighted; CDPP_D = 1.4826 MAD(running mean over D)",
                "method": "fast", "sde_logged_only": _defaults(bls.sde)},
        "detection_threshold": "configs/detection_threshold.json (derived from inversions; not in this hash)",
        "vet": dict(vet.THRESH),
        "deepvet": dict(deepvet.THRESH),
        "reject_flags": sorted(deepvet.REJECT),
        "warning_flags": ["centroid_shift_minor", "unresolved_neighbour", "neighbour_could_mimic",
                          "centroid_inconclusive", "ruwe_high", "gaia_nss"],
        "gaia": "DR3 via ESA TAP (cones, RUWE, non_single_star)",
        "ld": "Claret 2017 TESS quadratic, PHOENIX-COND r-method, nearest (Teff, logg) grid point",
    }


def config_hash():
    blob = json.dumps(config(), sort_keys=True, default=str).encode()
    return hashlib.sha256(blob).hexdigest()[:12]


def git_hash():
    try:
        h = subprocess.check_output(["git", "-C", str(REPO), "rev-parse", "--short=10", "HEAD"], text=True).strip()
        dirty = subprocess.run(["git", "-C", str(REPO), "diff", "--quiet", "--", "exoplanet_search/exosearch"]).returncode
        return h + ("-dirty" if dirty else "")
    except Exception:
        return "unknown"


def packages():
    out = {"python": platform.python_version()}
    for p in ("numpy", "scipy", "astropy", "wotan", "batman-package", "pandas", "pyarrow", "astroquery"):
        try:
            out[p] = metadata.version(p)
        except Exception:
            out[p] = None
    return out


def stamp():
    return dict(pipeline_version=PIPELINE_VERSION, config_hash=config_hash(), git_hash=git_hash())


if __name__ == "__main__":
    print(json.dumps(dict(stamp(), packages=packages(), config=config()), indent=1, default=str))
