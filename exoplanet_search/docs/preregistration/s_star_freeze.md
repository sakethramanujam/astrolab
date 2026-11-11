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

## Until then

The current value in `configs/interim_threshold_v1.json` is labelled provisional. It comes from a partial flipped null
and moves as that null grows.
