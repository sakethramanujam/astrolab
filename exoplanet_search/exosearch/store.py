"""SQLite results store (one row per tic+sector; upserts make re-runs idempotent)."""
from __future__ import annotations

import sqlite3

import pandas as pd

from . import config

COLS = {
    "tic": "INTEGER", "sector": "INTEGER", "status": "TEXT", "error": "TEXT", "backend": "TEXT",
    "host": "TEXT", "runtime_s": "REAL", "n_points": "INTEGER", "tmag": "REAL", "teff": "REAL",
    "radius": "REAL", "ra": "REAL", "dec": "REAL", "period": "REAL", "t0": "REAL", "duration": "REAL",
    "sde": "REAL", "depth": "REAL", "depth_err": "REAL", "snr": "REAL", "odd_even_sigma": "REAL",
    "depth_odd": "REAL", "depth_even": "REAL", "secondary_depth": "REAL", "secondary_sigma": "REAL",
    "n_transits": "INTEGER", "dur_expected": "REAL", "dur_ratio": "REAL", "rp_earth": "REAL",
    "flags": "TEXT", "passed": "INTEGER", "known_toi": "TEXT", "toi_disp": "TEXT",
    "toi_period_match": "INTEGER", "known_ctoi": "TEXT", "ctoi_period_match": "INTEGER",
    "known_planet": "TEXT", "planet_period_match": "INTEGER", "processed_at": "TEXT",
}


def connect(path=None):
    path = path or config.results_db()
    path.parent.mkdir(parents=True, exist_ok=True)
    con = sqlite3.connect(path)
    cols = ", ".join(f"{k} {v}" for k, v in COLS.items())
    con.execute(f"CREATE TABLE IF NOT EXISTS results ({cols}, PRIMARY KEY (tic, sector))")
    return con


def upsert(con, row):
    row = {k: row.get(k) for k in COLS}
    keys = ", ".join(row)
    q = ", ".join("?" for _ in row)
    con.execute(f"INSERT OR REPLACE INTO results ({keys}) VALUES ({q})", list(row.values()))
    con.commit()


def done_tics(con, sector):
    cur = con.execute("SELECT tic FROM results WHERE sector=? AND status='ok'", (sector,))
    return {r[0] for r in cur}


def frame(con, sector=None):
    q = "SELECT * FROM results" + (" WHERE sector=?" if sector is not None else "")
    return pd.read_sql_query(q, con, params=(sector,) if sector is not None else None)
