# v1.2 registration amendment: single host h2 (native arm64), 2026-11-19

Status: AMENDMENT to `registration_v1.2.md` (tag `v1.2-registered`, ad555e1). Dated 2026-11-19, committed before the h2 run starts.

## What changes
1. **Host.** The registered v1.2 run (real, invert and scramble_r0 rows for S59-S64; invert and scramble rows for S69; S69 real still waits on the frozen eligibility list per section 7.3) runs on **host h2 only**, natively on arm64. Calibration, the S* lock and the S69 gate are computed from h2 rows only. No rows from any other host are mixed in.
2. **Partial h1 rows.** The h1 run started 2026-11-18 is stopped once the h2 run is confirmed running. Its partial rows are archived with SHA-256 checksums and are **not read** for any registered quantity.
3. **Code pin unchanged.** Commit bc05650 (tag `v1.2-step7b-cfg-guard`), search-path hash `b9e5ebea8f84b05d60491e871f2c19e3facda75bdb146f3ab46eead42dafa4f4`.

## h2 environment pin
- arm64, macOS 27.2.1, Python 3.12.14, numpy 2.5.4 (Accelerate BLAS), scipy 1.18.1, numba 0.68.1.
- Sorted `pip freeze --all` of the worker environment (89 packages): SHA-256 `9a6ba5b270b02676034af026b11127ee627f6e649bdfac674ae38e741a706f5b`.
- Single-threaded math libraries (OMP/MKL/OPENBLAS/VECLIB/NUMBA threads = 1), nice 19.

## Input data
- Light curves for S59-S64 and S69 plus catalogs: 84,680 files, every file matched against the SHA-256 list taken on h1 (0 failures).
- Manifest: the registered h1 manifest (SHA-256 `aec741e97d47353679a66dcee6ff9ffef6c817772de6bf05a613679a0b0c5f36`, 83,130 star jobs) with only the host label changed to h2 (SHA-256 `9c5a394ea7ad957b8a3126ddd4c4fec492527c1df5b58b7b37d3801562d69aeb`). Jobs, tasks, seeds and settings are otherwise identical.

## Record check (201-star S70 set, bc05650), h2 native vs h1
- period, P_grid, duration and every flag and vet column: identical.
- S, SDE and other floats: differences at about 1e-11 only. S>=6: 78 stars on both. |rho_z|>=4: 5 stars on both.
- **Known host difference: rho_z** differs by up to 0.34 on 109 of 201 stars, with 2 sign flips. rho_z is therefore not bit-reproducible across hosts. Because every registered row comes from h2, this does not affect within-run comparisons. It is logged as a v1.3 issue.
