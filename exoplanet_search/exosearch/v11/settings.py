"""v1.1 identity: version, configuration (configs/v1.1.json), config hash and per-row provenance.

The beta option is part of the configuration and of the hash. It is null until the registered selection rule has been
applied (docs/preregistration/beta_treatment.md); production entry points refuse to start while it is null.
Provenance carries an opaque host_id (h1-h4; the mapping to machines lives outside the repo), the CPU
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
HOST_IDS = ("h1", "h2", "h3", "h4")
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


def exclusion_files(c=None):
    """{sector: repo-relative path} of the scrambled-null exclusion CSVs registered in the config.

    A registered file that is missing, or has no 'tic' column, is a hard error (never silently skipped).
    """
    import pandas as pd
    c = c or cfg()
    ex = c.get("scramble", {}).get("exclude_hosts_file")
    if not ex:
        return {}
    ex = ex if isinstance(ex, dict) else {"*": ex}
    out = {}
    for sec, rel in sorted(ex.items()):
        p = ROOT / rel
        if not p.is_file():
            raise FileNotFoundError(f"scramble exclusion file for sector {sec} is registered but missing: {rel}")
        if "tic" not in pd.read_csv(p, nrows=0).columns:
            raise ValueError(f"scramble exclusion file {rel} has no 'tic' column")
        out[str(sec)] = rel
    return out


def exclusion_tics(sector, c=None):
    """Excluded TIC set for a sector (empty set only when no file is registered for that sector)."""
    import pandas as pd
    files = exclusion_files(c)
    rel = files.get(str(sector), files.get("*"))
    if rel is None:
        return set()
    return set(pd.read_csv(ROOT / rel).tic.astype(int))


def config(beta_option=None):
    """Everything that can change a v1.1 result. beta_option overrides the file (dry runs only)."""
    from . import grid, search
    c = cfg()
    if beta_option is not None:
        c["beta_option"] = beta_option
    prior = ROOT / c["beta_prior_file"]
    return dict(pipeline_version=PIPELINE_VERSION, cfg=c,
                beta_prior_sha256=_sha(prior) if prior.exists() else None,
                scramble_exclude_sha256={k: _sha(ROOT / v) for k, v in exclusion_files(c).items()},
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


def tree():
    """Git tree object of HEAD (content identity of the checked-out code, independent of commit metadata)."""
    try:
        return subprocess.check_output(["git", "-C", str(ROOT), "rev-parse", "HEAD^{tree}"], text=True).strip()
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
    out = dict(pipeline_version=PIPELINE_VERSION, config_hash=config_hash(beta_option), commit=commit(), tree=tree(),
               host_id=hid, arch=platform.machine(), lib_versions=json.dumps(lib_versions(), sort_keys=True),
               beta_option=beta_option, production=bool(production))
    out["tree"] = tree()
    for sec, rel in exclusion_files(cfg()).items():
        out[f"scramble_exclude_file_{sec}"] = rel
        out[f"scramble_exclude_sha256_{sec}"] = _sha(ROOT / rel)
    return out
