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

`host_gate trace` hashes every stage. On h2 the raw light curve and the binned time stamps are bit-identical to h1;
the first divergence is the detrending trend (wotan biweight, numba-compiled), and the period grid also differs.
A probe on the two hosts:

* `np.arange(a, b, step)` differs in the last bit on arm64 (its fill loop is compiled with floating-point contraction,
  i.e. fused multiply-add, which is baseline on aarch64 and not used on x86-64 baseline builds);
  `a + step * np.arange(n)` (two separate ufuncs, no fusion possible) is bit-identical on both.
* `np.exp` / `np.log` differ in the last bit (different SIMD math libraries).
* `np.sum` and BLAS-free reductions agree. numpy, scipy, numba and llvmlite versions are identical (lockfile);
  BLAS differs (Accelerate vs OpenBLAS) but no BLAS call sits on the search path.

So the arm64 failure is FMA contraction / libm, not package versions. x86_64 hosts with the same lockfile are
bit-identical.

## Options for v1.1 (needs a decision; nothing changed in v1.0.0)

1. Keep arm64 out of the registered statistic (x86_64 hosts only), or
2. make the v1.1 search arm64-reproducible: period grid as `fmin + df * arange(n)` (or shipped as a file), a
   deterministic tie-break (lowest period among scan powers within 1e-9 relative of the max), and the trend computed
   without contraction (numba `fastmath=False` is the default; check wotan's compiled loops, or replace them with a
   NumPy-only biweight), then re-run this gate, or
3. amend the gate to "same grid index for P, |ΔS| ≤ 1e-6, and ties resolved by the rule above", registered before
   any merge.
