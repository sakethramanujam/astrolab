"""v1.1 identity: version, configuration (configs/v1.1.json), config hash and per-row provenance.

The beta option is part of the configuration and of the hash. It is null until the registered selection rule has been
applied (docs/preregistration/beta_treatment.md); production entry points refuse to start while it is null.
Provenance carries an opaque host_id (h1, h2, h3; the mapping to machines lives outside the repo), the CPU
architecture, library versions and the full commit. Never a hostname, user name or path.
"""
from __future__ import annotations

import hashlib
import json
import os
import platform
import subprocess
from pathlib import Path

PIPELINE_VERSION = "1.1.0"
ROOT = Path(__file__).resolve().parents[2]            # exoplanet_search/
CFG_FILE = ROOT / "configs" / "v1.1.json"
HOST_IDS = ("h1", "h2", "h3")
BETA_OPTIONS = ("a", "b", "c", "d")
LIBS = ("numpy", "scipy", "astropy", "numba", "llvmlite", "wotan", "batman-package", "pandas", "pyarrow")


def cfg(path=CFG_FILE):
    return json.loads(Path(path).read_text())


def beta_prior(c=None):
    c = c or cfg()
    p = ROOT / c["beta_prior_file"]
    return json.loads(p.read_text()) if p.exists() else None


def _sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def config(beta_option=None):
    """Everything that can change a v1.1 result. beta_option overrides the file (dry runs only)."""
    from . import grid, search
    c = cfg()
    if beta_option is not None:
        c["beta_option"] = beta_option
    prior = ROOT / c["beta_prior_file"]
    return dict(pipeline_version=PIPELINE_VERSION, cfg=c,
                beta_prior_sha256=_sha(prior) if prior.exists() else None,
                durations_d=[float(x) for x in search.DURATIONS], grid=grid.GRID_SPEC,
                quality_mask="QUALITY == 0, finite TIME/PDCSAP_FLUX/PDCSAP_FLUX_ERR, ERR > 0; flux / median",
                detrend=dict(method="biweight", window=0.75, edge_cutoff=0.0, break_tolerance=0.5, clip_sigma=4.0,
                             gap_trim=0.3, numba="fastmath=False, cache=False, parallel=False"),
                bin_minutes=10.0, tie_rel=search.TIE_REL,
                statistic=("S = depth_exact * sqrt(N_tr) / CDPP_D * sqrt(1 - D/P) * beta_D / beta'_D at the chosen peak; "
                           "depth_exact = weighted in-minus-out box mean on the binned flux (momentum-dump epochs excluded); "
                           "peak = max of that S over the top-k distinct peaks of the capped, beta-treated BLS snr scan; "
                           "plateau tie-break (tie_break.md amendment 1/1a)"),
                deep_checks="exosearch.deepvet light-curve checks + exosearch.density, as v1.0.0")


def config_hash(beta_option=None):
    blob = json.dumps(config(beta_option), sort_keys=True, default=str).encode()
    return hashlib.sha256(blob).hexdigest()[:12]


def commit():
    try:
        return subprocess.check_output(["git", "-C", str(ROOT), "rev-parse", "HEAD"], text=True).strip()
    except Exception:
        return "unknown"


def lib_versions():
    import importlib.metadata as md
    v = {}
    for k in LIBS:
        try:
            v[k] = md.version(k)
        except Exception:
            v[k] = None
    v["python"] = platform.python_version()
    return v


def host_id(explicit=None):
    h = explicit or os.environ.get("EXO_HOST_ID")
    if h not in HOST_IDS:
        raise SystemExit(f"host id must be one of {HOST_IDS} (opaque; set EXO_HOST_ID or --host-id)")
    return h


def provenance(hid, beta_option, production):
    return dict(pipeline_version=PIPELINE_VERSION, config_hash=config_hash(beta_option), commit=commit(),
                host_id=hid, arch=platform.machine(), lib_versions=json.dumps(lib_versions(), sort_keys=True),
                beta_option=beta_option, production=bool(production))
