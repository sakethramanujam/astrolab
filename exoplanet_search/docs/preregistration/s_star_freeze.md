# Registration: freezing S* at a fixed null size

Interim analysis layer over frozen v1.0.0. Registered 2026-11-10. Purpose: stop S* drifting as the null grows.

## Fixed null size per sector

S* is frozen once per sector, when, and only when, both of these exist for that sector:

1. **The complete flipped (inverted) null**: one realization over every star searched in the sector.
2. **Three scrambled realizations**, r0, r1, r2 (seeded ~1-day block permutations), over the same stars.

The fit uses tuning (non-holdout) stars that pass `star_quality_v1` (beta_max <= 3). Those numbers are fixed now:

| Sector | Searched | Pass quality cut | Tuning + quality (null stars per realization) | Total null star-realizations (1 flipped + 3 scrambled) |
|---|---|---|---|---|
| 70 | 12,797 | 10,244 | **8,156** | **32,624** |
| 106 | 12,996 | 10,531 | **8,377** | **33,508** |

Any other sector gets an S* only when its own complete flipped null and three scrambled realizations exist. Until then
it is "not calibrated".

## Computation (unchanged from the interim layer)

* Vetted null as defined in `analysis/SCHEMA.md`; GPD on vetted max S above u; u by the parameter-stability rule
  (lowest u >= 6 with k >= 50 whose xi is within 1 bootstrap SE of xi at every higher u with k >= 30).
* False-alarm budget alpha = 1 / N_stars, with N_stars = quality-passing stars searched in the sector (S70: 10,244).
* S*_flip from the flipped null (8,156 stars). S*_scr from the three scrambled realizations pooled (per-realization
  rate, 3 x 8,156 star-realizations). Both are reported side by side, each with a 95% star-level bootstrap CI (2,000
  resamples, seed 20261110) and the systematic spread over the stable u range.
* **Frozen S* = max(S*_flip, S*_scr).** It is written once to `configs/s_star_frozen_v1.json` with the null sizes,
  hashes of the input tables and the commit, and is not updated afterwards. Later null realizations are reported as
  checks, not refits.

## Condition to apply before freezing: red-noise floor

The S70 flipped null has a non-trivial number of peaks with beta_D < 1 (2,932 of 5,881 tuning stars; 116 of the 288
vetted peaks above S = 6; snapshot 2026-11-10, `exosearch/beta_floor.py`). Before S* is frozen, the statistic used
for the fit, for completeness and for the candidates must apply the floor CDPP_D' = max(CDPP_D, sigma_w / sqrt(N_D))
(beta_D >= 1), consistently on all populations. Post hoc on the stored peaks it moves S* from 12.07 to 11.62 on the
same sample. If the floor is implemented in the search (v1.1) rather than post hoc, the frozen S* is derived from
that run instead. Nothing is frozen by this file.

## Amendment 2026-11-10: the frozen S* comes from v1.1

v1.0.0 cannot produce the frozen S*: its statistic uses the BLS fast-method depth (inflated on some non-transit
peaks; depth audit), applies no red-noise floor and does not prewhiten pulsators. Post-hoc corrections on stored
peaks only give lower bounds. Therefore:

* **The frozen S* is derived from a v1.1 re-search** in which S uses the exact box depth, the chosen beta treatment
  (floor at 1 − 2·sigma_beta) and pulsator prewhitening are inside the statistic, and the peak is chosen with them.
* It uses **fresh v1.1 nulls at the registered sizes above** (full flipped null + three scrambled realizations per
  sector, same tuning + quality star counts), and v1.1 injections re-scored with the same statistic.
* **The v1.0 nulls feed only the interim layer** (diagnostics, provisional S*, completeness lower bounds). They are
  capped at the flipped null for S70 and S106 plus one scrambled realization (r0); see `docs/interim_lane_cap.md`.

The red-noise section above is superseded by the v1.1 statistic; the computation section is unchanged.

## Until then

The current value in `configs/interim_threshold_v1.json` is labelled provisional. It comes from a partial flipped null
and moves as that null grows.
