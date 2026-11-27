# v1.2 catalogue cross-check spec, 2026-11-19 (committed before any row is read for it)

Purpose: test our implementation against outside truth, not against itself. Our BLS core is astropy `BoxLeastSquares`, so an astropy-vs-ours comparison would be circular and is withdrawn. Read-only diagnostic. It does not change any registered rule, cut, S* or gate.

Signed off in the Astrolab room by Scientist and RF Engineer on 2026-11-19.

## Scope
- h2 `real` rows only, sectors S59-S64, non-holdout stars (`exosearch.holdout`, SEED 20261110, FRACTION 0.20) that host confirmed or known planets (CP/KP).
- Holdout rows stay sealed. S69 rows are not read until S*_gate is frozen.

## Truth set
- The frozen local catalogues: ExoFOP TOI table and NASA Exoplanet Archive pscomppars, at the versions pinned in the data-root checksum list.

## Labelled before the run
- Planets whose catalogue T14 is over 0.25 d are labelled **outside our grid**. They are not scored as failures and they go to TLS.

## Tolerances (a planet passes only if all hold)
- **Period:** within max(one period-grid step, 3σ of the published period). P/2 or 2P is logged as an **alias**, not a pass.
- **t0:** the published ephemeris and its error are carried forward to the sector, then our t0 must be within 0.5 × catalogue T14.
- **Duration:** our chosen grid duration is the grid value nearest the catalogue T14, or the one next to it.
- **Depth:** ratio ours / catalogue within [0.5, 1.5]. A box fit to a limb-darkened transit is biased shallow, and pipelines correct for blending differently, so this is a ratio band, not a σ test.

## Follow-up
- Every mismatch is traced in our code first, before blaming the star or the catalogue.
- TLS runs only on mismatches and on the planets outside our grid, as a second search code.
