# Capping the v1.0 reliability lanes (interim layer only)

Agreed 2026-11-10 with the Scientist: v1.0.0 cannot produce the frozen S* (depth audit, beta treatment,
prewhitening; see `docs/preregistration/s_star_freeze.md`), so the v1.0 nulls and injections run only as far as the
interim diagnostics need. `scripts/cap_lanes.sh` (user unit `exo-cap`) enforces it without touching a running step:
it waits for the last wanted step of each lane to log DONE, then stops the lane.

| lane | kept | dropped |
|---|---|---|
| a (nulls) | flipped null S70 + S106 (full), scrambled r0 (S70 + S106) | scrambled r1, r2 |
| b (injections) | run 0 (S70 + S106, 8,000), then run 1 for S70 only (4,000) as `exo-lane-b2` | runs 1–7 for both sectors |
| c | stellar parameters, labels | injection runs 8–14 |

Total v1.0 injections: 12,000 (S70 8,000, S106 4,000). The interim S* stays provisional and is reported with the
flipped and scrambled (r0) values side by side. Freed cores (cap 12 for batches) go to the v1.1 nulls and injections.
