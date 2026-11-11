# Cross-host gate: dry run on the v1.0.0 search path (2026-11-10)

`python -m exosearch.host_gate` (run / compare / trace). 200 S70 stars (`configs/host_gate_tics.txt`, seed 20261110),
same commit (e0fdde34), same hash-pinned lockfile (`requirements.lock`), float64 CPU, 1 thread per worker, nice 19.
Reference h1. Gate: |ΔS| ≤ 1e-6 and |ΔS_exact| ≤ 1e-6 and identical best (P, D) for every star. Summary:
`summary_dryrun_v1.0.0.json` (host_id, arch, library versions only).

| host | arch / OS | BLAS | workers | max abs ΔS | best (P, D) differs | gate |
|---|---|---|---|---|---|---|
| h1 | x86_64 Linux | scipy-openblas 0.3.34 | 2 | 0 (reference) | 0 | pass |
| h2 | arm64 macOS | Accelerate | 12 | 1.1e-11 | 23 of 200 (17 last-bit P only; 6 different peak) | **fail** |
| h3 | x86_64 Linux | scipy-openblas 0.3.34 | 6 | 0 (bit-identical) | 0 | pass |

All 6 h2 peak flips are exact or 1e-15 ties in S between two periods with the same transits (P ≈ 6–11 d, where
two to three transits allow several periods with the same in-transit points); S agrees to 1e-11 in every star.

## Diagnosis

`host_gate trace` hashes every stage. On h2 the raw light curve and the binned time stamps are bit-identical to h1.
Two sources of difference, both on arm64 only:

* **Period grid.** `np.arange(a, b, step)` differs in the last bit on arm64 (its fill loop is compiled with
  floating-point contraction, i.e. fused multiply-add, baseline on aarch64 and not used by x86-64 baseline builds);
  `a + step * np.arange(n)` (two separate ufuncs) is bit-identical on both. This alone explains the 17 last-bit
  period differences and, through the S ties, the 6 peak flips. `np.exp` / `np.log` also differ in the last bit
  (different SIMD math libraries), though they are not on the grid path.
* **Detrending trend.** wotan's biweight is compiled by numba with `fastmath=True, cache=True`. On h2 a freshly
  JIT-compiled run gives the trend bit-identical to h1, while the code loaded from wotan's on-disk numba cache gives a
  different trend (reproducible: empty `NUMBA_CACHE_DIR` → h1 hash; second run from that cache → different hash).
  `fastmath` permits reassociation, so the result depends on how and where LLVM compiled it. Mixed fresh/cached
  workers make h2 rows inconsistent with each other as well. h1 and h3 (both x86_64) agree bit for bit.
* Package versions are identical (lockfile). BLAS differs (Accelerate vs OpenBLAS) but no BLAS call sits on the
  search path. `np.sum` agrees.

## Options for v1.1 (needs a decision; nothing changed in v1.0.0)

1. Keep arm64 out of the registered statistic (x86_64 hosts only), or
2. make the v1.1 search arm64-reproducible: period grid as `fmin + df * arange(n)` (or shipped as a file), a
   deterministic tie-break (lowest period among scan powers within 1e-9 relative of the max), and the biweight
   compiled without `fastmath` (vendored copy, or a NumPy-only biweight), with the numba cache disabled or keyed per
   host; then re-run this gate on all hosts (the `fastmath` issue is also a risk between x86_64 CPU generations), or
3. amend the gate to "same grid index for P, |ΔS| ≤ 1e-6, and ties resolved by the rule above", registered before
   any merge.

## Dry run with the adopted fixes (`--fixed`)

Fixes (interim code only, `exosearch/host_gate.py: apply_fixes`): period grid `fmin + df·arange(n)`; tie-break
|S_i − S_max| ≤ 1e-9·S_max → lowest period, then shortest duration (`docs/preregistration/tie_break.md`); wotan's 24
numba kernels recompiled with `fastmath=False, cache=False, parallel=False`. Summary: `summary_dryrun_fixed.json`.

| host | arch | max abs ΔS | max abs ΔS_exact | best (P, D) differs | stars with > 1 grid point in the tie band | gate |
|---|---|---|---|---|---|---|
| h1 | x86_64 | 0 (reference) | 0 | 0 | 25 / 200 (12.5%; max 38) | pass |
| h2 | arm64 | 7.1e-15 | 0 | 0 | 25 / 200 (same stars, same counts) | **pass** |
| h3 | x86_64 | 0 | 0 | 0 | 25 / 200 | pass |

The 6 stars that switched peaks before (38349319, 51988338, 333557762, 353927437, 365179665, 611355482) now resolve to
the same (P, D) on all hosts; they had 4–27 grid points in the tie band. The fixed grid is a different float grid
from v1.0.0, so v1.1 rows are not bit-comparable with v1.0.0 rows (S differs by ≤ 0.007 on h1 between the two paths).
Stopgap until the fixes are in the v1.1 code and the gate is re-run there: registered statistic from h1/h3 only.

Throughput on the gate (search time only, at the worker counts shown; wall time includes the MAST downloads):
h1 ~700 stars/h per worker (idle priority under the running lanes), h2 ~2,650 stars/h per worker (12 workers:
~32,000/h; 200 stars in 24 s wall), h3 ~470 stars/h per worker (6 workers: ~2,800/h).

## Dry run with the amended tie-break (plateau rule, 2026-11-10 21:30 PST)

Same 200 S70 stars and fixes, tie-break per `docs/preregistration/tie_break.md` amendment 1 (P = plateau centre,
P_grid_halfwidth = half its width, lowest period only between disjoint plateaus/peaks). Commit dcd5914a on all hosts.
Summary: `summary_dryrun_plateau.json`.

| host | arch | max abs ΔS | max abs ΔS_exact | best (P, D) differs | P_grid_halfwidth differs | gate |
|---|---|---|---|---|---|---|
| h1 | x86_64 | 0 (reference) | 0 | 0 | 0 | pass |
| h2 | arm64 | 7.1e-15 | 0 | 0 | 0 | **pass** |
| h3 | x86_64 | 0 | 0 | 0 | 0 | pass |

The 6 former switch stars resolve identically on all hosts: 38349319 and 611355482 are single plateaus (14 and 27
grid steps, half-widths 0.059 and 0.090 d); 51988338, 333557762, 353927437 and 365179665 have 6, 4, 9 and 17 disjoint
in-band plateaus/peaks (chosen plateau 1–2 steps). Still a dry run on the interim code; the registered statistic stays
h1/h3-only until the gate passes on the v1.1 code.

## Plateau statistics (`plateau_report.json`, `python -m exosearch.plateau_report`)

Fixed-grid search with the plateau rule, on the 200 S70 gate stars and on the 200 v1.0.0 S105 rows of the numba-cache
exposure check (x86_64).

| set | stars | single grid point | one plateau (width > 1) | disjoint in-band peaks | Spearman(n_tr, max width) |
|---|---|---|---|---|---|
| S70 gate | 200 | 175 | 17 (8.5%) | 8 (4.0%) | −0.35 |
| S105 v1.0 sample | 200 | 162 | 27 (13.5%) | 11 (5.5%) | −0.48 |

Of the tied stars (more than one grid point in band) 68% (S70) and 71% (S105) are plateau-only and 32% / 29% have
disjoint peaks. Plateaus wider than one step occur almost only at n_tr ≤ 2: 16 of 64 (S70) and 27 of 100 (S105) stars
with ≤ 2 transits, against 5 of 136 and 2 of 100 with ≥ 3 (max widths 36 and 35 steps at n_tr ≤ 2, ≤ 3 steps at
n_tr ≥ 3). Distance from the lowest-edge period (the original rule) to the plateau centre, over plateau stars:
median 0.0047 d (S70) / 0.0056 d (S105), p90 0.087 / 0.060 d, max 0.16 d; in units of the half-width median 0.67, max 0.97
(< 1 by construction, since the half-width includes the half-steps at the edges). Widths are reported only; no grid change is
proposed (any change follows the Ofir 2014 criterion and is registered before v1.1 production).

**v1.0.0 cached vs fresh kernels:** the 26 S105 rows that differ differ only in the last bit of P (|ΔP| ≤ 1.8e-15 d,
same grid point, identical D and T0); 26/26 stay inside one plateau and every shift is ≤ 5e-12 of the plateau
half-width. (Earlier text called these (P, D) flips; corrected.)
