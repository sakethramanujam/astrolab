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
