# v1.2 catalogue cross-check: amendment, 2026-11-20

Status: AMENDMENT to `xcheck_catalogue_v1.2_2026-11-19.md` (tag `v1.2-xcheck-spec-2026-11-19`, 5ef8110). Read-only diagnostic. It does not change any registered rule, cut, S* or gate.

**This change was made after the original-spec results were seen.** It is agreed in the Astrolab room and approved by the project owner on 2026-11-20. The original-spec result stays on record below and is reported alongside every amended result.

## Original-spec result (h2 rows, non-holdout CP/KP hosts, S59-S64)
- In scope: 420 planets on 349 stars (696 planet x sector rows). Labelled before scoring: 42 outside the duration grid, 2 with incomplete catalogue values. Scored: 376 planets (600 planet x sector rows).
- Planet level (best sector row): period 118, t0 243, duration 289, depth 323; **full pass 103 / 376**.
- Planet x sector level: period 144, t0 325, duration 404, depth 463; full pass 126 / 600.
- Aliases (planet level): P/2 5, 2P 0. Period matching another planet on the same star: 18.
- Why it is amended: the period tolerance max(one grid step, 3 sigma) is far finer than a single sector can measure (the grid oversamples, and the catalogue 3 sigma is about 0.02 of a step). 106 of the 235 period-mismatch planets lie within 1% of the catalogue period, and 68 have catalogue periods outside the searched range.

## Amended rules
- **Period:** pass if |dP| <= P * T14 / (2 * T_span). P and T14 are the catalogue period and duration; T_span is that row's real data span: last minus first cadence passing the pipeline quality mask (`lightcurve.good_mask`: finite flux and error, QUALITY == 0, error > 0) in that sector's light curve.
- **t0:** pass if |dt0| <= max(T14 / 2, 3 * sigma_pred), where sigma_pred is the catalogue epoch error and period error carried forward to our epoch: sqrt(sigma_T0^2 + (n * sigma_P)^2).
- **Duration and depth:** unchanged (grid value nearest catalogue T14 or adjacent; depth ratio ours / catalogue in [0.5, 1.5]).

## Labels before scoring (none counts as a pass)
- **Outside the duration grid:** catalogue T14 > 0.25 d.
- **Outside the searchable period range:** catalogue period above the per-star maximum (baseline / 2, from that row's search grid) or below 0.5 d.
- Incomplete catalogue values (period, t0, T14 or depth missing) stay labelled as such.

## Sibling case
- If our peak matches another catalogue planet on the same star (by the amended period rule at that planet's period), the sibling is credited as recovered and the target planet is labelled **found sibling instead**. It counts as not recovered by a single-peak search. Masking the found signal and re-searching goes on the v1.3 list.
- A sibling match is classed before an alias match.

## Aliases
- Checked ratios ours / catalogue: 1:2, 2:1, 1:3, 3:1, 2:3, 3:2, each with the amended period rule applied at the aliased period. Aliases are logged, never passes.

## Detectable set for the first gate (defined now, computed later)
- A planet is **detectable** only if all hold: expected S >= S*_cal, catalogue period inside the searchable range, and catalogue T14 inside the duration grid.
- This set is computed once S*_cal exists. It is not computed by this amendment.
