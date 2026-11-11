# v1.1 pre-production checks

Machine-readable results that `scripts/v11_launch.py` reads from origin before it launches anything. Production
needs `all_pass: true` in both files.

## 3-host gate (`host_gate.json`)

200 stars (S70), code at 130df5dc, `run work --dry-run --beta-option c`, compared with `run compare` (|dS| <= 1e-6,
exact period, duration and P_grid_halfwidth for the primary and every option a-d).

| host | arch | rows ok | fails | max abs dS |
|---|---|---|---|---|
| h1 | x86_64 | 200 | 0 | reference |
| h2 | arm64 | 200 | 0 | 3.0e-11 |
| h3 | x86_64 | 200 | 0 | 0 (bit-identical) |

**PASS.** The first run (bb805a95) failed on one h2 star: TIC 353927437, a single transit whose disjoint
refinement peaks are exactly tied, so their order followed last-bit differences. 130df5dc orders local maxima with
the registered 1e-9 tie band (lowest period first). After that fix all three hosts agree. Same 28 prewhitened stars,
24 stars with more than one peak in the tie band, and 19 with a plateau wider than one grid step on every host.

## Prewhitening and duration-cap validation (`validation.json`)

See `docs/v1.1_validation.md`.
