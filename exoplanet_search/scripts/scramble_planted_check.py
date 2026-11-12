"""Planted-signal test script for v1.1 scramble rule.

Uses real sector-70 cadence and gap structure from tests/fixtures/s0070_cadence.npz,
white noise plus mild red noise, a box transit train at P = 1, 3 and 6 d
(depth giving > 20 sigma unscrambled), at least 10 seeds each.
Computes phase-folded depth significance at injected P/t0 after scrambling.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Dict, List, Tuple

# Ensure exoplanet_search package is importable when run directly as scripts/scramble_planted_check.py
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np
from astropy.timeseries import BoxLeastSquares

from exosearch.v11 import scramble, settings


def run_planted_check(
    cadence_file: Path | str,
    output_file: Path | str = "scramble_planted_check.txt",
    n_seeds: int = 10,
    periods: Tuple[float, ...] = (1.0, 3.0, 6.0),
) -> List[Dict[str, float]]:
    cfg = settings.cfg()
    d = np.load(cadence_file)
    t = d["t"]

    noise_white = 0.0004
    noise_red_amp = 0.0001
    dur = 0.1
    t0 = float(t[0] + 0.1)

    results = []

    for P in periods:
        depth = 0.0006  # gives > 20 sigma unscrambled for all P in (1, 3, 6)

        # Confirm unscrambled significance > 20 sigma
        rng_u = np.random.default_rng(0)
        w_u = rng_u.normal(0, noise_white, size=len(t))
        r_u = noise_red_amp * np.sin(2 * np.pi * (t - t[0]) / 3.5)
        f_u = 1.0 + w_u + r_u
        e_u = np.full_like(t, noise_white)
        phase_u = ((t - t0 + 0.5 * P) % P) - 0.5 * P
        f_u[np.abs(phase_u) < dur / 2.0] -= depth
        bls_u = BoxLeastSquares(t, f_u, dy=e_u)
        st_u = bls_u.compute_stats(P, dur, t0)
        unscr_sigma = float(st_u["depth"][0] / st_u["depth_phased"][1])
        assert unscr_sigma > 20.0, f"Unscrambled significance {unscr_sigma:.2f} <= 20 for P={P}"

        sigmas = []
        for s in range(n_seeds):
            rng = np.random.default_rng(s)
            w = rng.normal(0, noise_white, size=len(t))
            r = noise_red_amp * np.sin(2 * np.pi * (t - t[0]) / 3.5)
            f = 1.0 + w + r
            e = np.full_like(t, noise_white)
            phase = ((t - t0 + 0.5 * P) % P) - 0.5 * P
            f[np.abs(phase) < dur / 2.0] -= depth

            t_scr, f_scr, e_scr, _ = scramble.segment_rotate(t, f, e, seed=s, cfg=cfg)
            bls = BoxLeastSquares(t_scr, f_scr, dy=e_scr)
            st = bls.compute_stats(P, dur, t0)
            sig = float(abs(st["depth"][0]) / st["depth_phased"][1])
            sigmas.append(sig)

        med_sig = float(np.median(sigmas))
        max_sig = float(np.max(sigmas))
        results.append({
            "P": P,
            "median_sigma": med_sig,
            "max_sigma": max_sig,
            "n_seeds": n_seeds,
        })

    # Print table
    header = f"{'P':<6} {'median_sigma':<14} {'max_sigma':<12} {'n_seeds':<8}"
    lines = [header]
    for r in results:
        lines.append(f"{r['P']:<6.1f} {r['median_sigma']:<14.2f} {r['max_sigma']:<12.2f} {r['n_seeds']:<8d}")
    table_str = "\n".join(lines)
    print(table_str)

    if output_file:
        p = Path(output_file)
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(table_str + "\n")

    return results


def main() -> None:
    parser = argparse.ArgumentParser(description="Check planted transit signal significance after scrambling.")
    default_cadence = Path(__file__).resolve().parent.parent / "tests" / "fixtures" / "s0070_cadence.npz"
    parser.add_argument("--cadence", default=str(default_cadence), help="Path to cadence fixture .npz")
    parser.add_argument("--out", required=True, help="Output file path (relative or absolute)")
    parser.add_argument("--n-seeds", type=int, default=10, help="Number of seeds to test")
    args = parser.parse_args()

    run_planted_check(args.cadence, output_file=args.out, n_seeds=args.n_seeds)


if __name__ == "__main__":
    main()
