"""Calibrate the S statistic on pure white noise with each star's real cadence and gap pattern.

For each star: keep the quality-masked time stamps, replace the flux with 1 + sigma * N(0, 1) (seeded),
run the frozen detrend + bin + S search, and record
  * the per-star max S (should follow the Gaussian max over the effective number of trials),
  * the CDPP_D estimator (1.4826 MAD of the overlapping running mean) against the truth sigma / sqrt(n_D),
    and against a non-overlapping-block estimate, for every trial duration D.
  python -m exosearch.whitenoise --sector 70 --n 300 --workers 2
Each star goes through the full chain (nulltest.run_one: search, screens, deep light-curve checks).
Output: products/v<ver>/calibration/whitenoise_s<sector>.parquet
"""
from __future__ import annotations

import argparse
import hashlib

import numpy as np
import pandas as pd

from . import batch, bls, config, holdout, lightcurve, nulltest, stellar, store, targets, task
from .version import stamp

SIGMA = 1e-3

_transform = nulltest.transform


def _white_transform(mode, t, f, e, meta, seed):
    """Full-chain hook: nulltest.run_one with mode 'white' replaces the flux by seeded Gaussian noise."""
    if mode != "white":
        return _transform(mode, t, f, e, meta, seed)
    rng = np.random.default_rng(seed)
    return t, 1.0 + SIGMA * rng.standard_normal(t.size), np.full(t.size, SIGMA)


nulltest.transform = _white_transform        # analysis-layer hook; the frozen module itself is unchanged


def _seed(sector, tic, base=20261110):
    return int.from_bytes(hashlib.sha256(f"{base}:white:{sector}:{tic}".encode()).digest()[:4], "big")


def block_cdpp(t, r, D):
    """Non-overlapping D-windows (>= half full), 1.4826 MAD and plain std of their means."""
    idx = np.floor((t - t.min()) / D).astype(int)
    n = np.bincount(idx)
    s = np.bincount(idx, r)
    cad = np.median(np.diff(t))
    ok = n >= 0.5 * D / cad
    m = s[ok] / n[ok]
    return float(1.4826 * np.median(np.abs(m - np.median(m)))), float(np.std(m)), float(np.mean(1 / np.sqrt(n[ok])))


def run_one(job):
    row = nulltest.run_one(job)              # identical search + screens + deep checks as the inverted null
    try:
        t, f, e, meta = lightcurve.load_spoc(job["path"])
        rng = np.random.default_rng(job["seed"])
        f = 1.0 + SIGMA * rng.standard_normal(t.size)
        e = np.full(t.size, SIGMA)
        d = task.prepare_arrays(t, f, e, meta)
        tb, fb, _ = d["binned"]
        tf, ff, _ = d["flat"]
        _, cad, prof = bls.noise_profile(tf, ff)      # same estimator the search used (same seed -> same flux)
        row.update(n_cad=int(tf.size), baseline=float(tf.max() - tf.min()), n_periods=int(len(bls.period_grid(tb))))
        r = ff - np.median(ff)
        row["sigma_w_true"] = SIGMA
        for D in bls.DURATIONS:
            k = f"{D:.2f}"
            # truth for a white-noise running mean with gaps: sigma / sqrt(n_D), n_D = cadences actually in the window
            j = np.searchsorted(tf, tf + D)
            n = j - np.arange(tf.size)
            ok = n >= 0.5 * D / cad
            row[f"cdpp_est_{k}"] = prof[float(D)][0]
            row[f"cdpp_true_{k}"] = float(np.median(SIGMA / np.sqrt(n[ok])))
            bm, bs, _ = block_cdpp(tf, r, D)
            row[f"cdpp_block_mad_{k}"] = bm
            row[f"cdpp_block_std_{k}"] = bs
            row[f"cdpp_runstd_{k}"] = float(np.std((np.r_[0, np.cumsum(r)][j[ok]] - np.r_[0, np.cumsum(r)][np.arange(tf.size)[ok]]) / n[ok]))
        if row.get("status") == "ok":
            kb = f"{row['duration']:.2f}"
            row["S_true"] = row["S"] * row["cdpp_d"] / row[f"cdpp_true_{kb}"]
    except Exception as exc:
        row.update(status="error", error=f"diag {type(exc).__name__}: {exc}"[:300])
    return row


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sector", type=int, default=70)
    ap.add_argument("--n", type=int, default=300)
    ap.add_argument("--workers", type=int, default=2)
    a = ap.parse_args()
    st = stamp()
    con = store.connect()
    df = store.frame(con, a.sector)
    df = df[df.status == "ok"]
    man = targets.fetch_manifest(a.sector).set_index("tic")
    tics = [int(x) for x in df.tic if int(x) in man.index and not holdout.is_holdout(int(x))]
    rng = np.random.default_rng(20261110)
    tp = stellar.tic_params(tics)
    jobs = []
    for tic in rng.permutation(tics):
        p = batch.lc_path(a.sector, man.loc[tic].fname, download=False)
        if p is not None:
            x = tp.loc[tic].to_dict() if tic in tp.index else {}
            rho, erho = stellar.tic_rho(x) if x else (None, None)
            jobs.append(dict(key=f"white:{a.sector}:{tic}", tic=int(tic), sector=a.sector, mode="white", realization=0,
                             seed=_seed(a.sector, tic), path=str(p), stamp=st, rho=rho, rho_err=erho,
                             teff=x.get("Teff"), logg=x.get("logg"), ruwe=None, nss=None))
        if len(jobs) >= a.n:
            break
    out = config.DATA_ROOT / "products" / f"v{st['pipeline_version']}" / "calibration"
    with batch.claim(a.workers):
        res = batch.pmap(run_one, jobs, out / f"whitenoise_s{a.sector:04d}_full_parts", a.workers, part_size=50)
    res = nulltest.finalize(res)
    res.to_parquet(out / f"whitenoise_s{a.sector:04d}.parquet", index=False)
    ok = res[res.status == "ok"]
    print(len(ok), "stars; max S median", ok.S.median(), "p95", ok.S.quantile(.95))


if __name__ == "__main__":
    main()
