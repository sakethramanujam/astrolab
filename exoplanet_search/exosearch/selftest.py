"""Backend self-test + benchmark on real TESS stars (used by the worker bootstrap scripts).

    python -m exosearch.selftest --backend mlx [--tics 301289516 251039147] [--sector 70]

Exits non-zero if the backend's best period disagrees with astropy BLS on a star whose
astropy SDE is >= 9 (i.e. where there is a clear answer to agree on).
"""
from __future__ import annotations

import argparse
import sys
import time

from . import bls, config, targets
from .task import prepare

DEFAULT_TICS = [301289516, 251039147, 270380593]   # TOI-4517.01, TOI-4537.01, TOI-465.01 (sector 70)


def _harm(p, q, tol=0.002):
    return any(abs(p * k / q - 1) < tol for k in (1, 2, 0.5))


def run(backend, sector=70, tics=DEFAULT_TICS, log=print):
    man = targets.fetch_manifest(sector).set_index("tic")
    ok = True
    rows = []
    for tic in tics:
        m = man.loc[tic]
        path = targets.download_one(m.fname, m.url, config.lc_dir(sector))
        tb, fb, eb = prepare(path)["binned"]
        t0 = time.perf_counter()
        ref = bls.search(tb, fb, eb, backend="cpu")
        t1 = time.perf_counter()
        got = bls.search(tb, fb, eb, backend=backend)
        t2 = time.perf_counter()
        agree = _harm(got["period"], ref["period"])
        if ref["sde"] >= 9 and not agree:
            ok = False
        rows.append((tic, ref["period"], got["period"], ref["sde"], got["sde"], t1 - t0, t2 - t1, agree))
        log(f"TIC {tic}: astropy P={ref['period']:.5f} (SDE {ref['sde']:.1f}, {t1-t0:.2f}s)  "
            f"{backend} P={got['period']:.5f} (SDE {got['sde']:.1f}, {t2-t1:.2f}s)  {'OK' if agree else 'MISMATCH'}")
    return ok, rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--backend", default="numpy", choices=["numpy", "cuda", "mlx", "cpu"])
    ap.add_argument("--sector", type=int, default=70)
    ap.add_argument("--tics", type=int, nargs="*", default=DEFAULT_TICS)
    a = ap.parse_args()
    ok, _ = run(a.backend, a.sector, a.tics)
    print("SELFTEST", "PASS" if ok else "FAIL")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
