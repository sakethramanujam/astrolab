"""Paths. Override the data root with EXO_DATA_ROOT on hosts without the shared disk."""
import os
from pathlib import Path

DATA_ROOT = Path(os.environ.get("EXO_DATA_ROOT", "/media/nfs01/astro_data/tess"))
CATALOG_DIR = DATA_ROOT / "catalogs"
RESULTS_DIR = DATA_ROOT / "results"


def sector_dir(sector: int) -> Path:
    return DATA_ROOT / f"s{int(sector):04d}"


def lc_dir(sector: int) -> Path:
    return sector_dir(sector) / "lc"


def plots_dir(sector: int) -> Path:
    return RESULTS_DIR / "plots" / f"s{int(sector):04d}"


def results_db() -> Path:
    return RESULTS_DIR / "exosearch.sqlite"
