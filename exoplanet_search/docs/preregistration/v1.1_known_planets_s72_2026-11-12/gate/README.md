# v1.1 sector 72 known-planet gate (amended gate + amended vetter, S* = 13.71): **FAIL, 28 of 30**

Last gate attempt on the amended vetter (docs/preregistration/v1.1_amendment_vetter_2026-11-12.md). Per that document the
failure is reported and **the vetter is not amended again**.

* Frozen list: `../frozen_eligibility_s72_s13.71.csv` (sha256 15650be9…, pushed in 0400dbb before any S72 output).
* Search: 69 eligible known hosts of sector 72, real search only, commit 0400dbb (search-path hash fb958a4e…), host h1
  (x86_64), production rows, nice 19, 10 workers, thermal guard (pause 90 C / resume 80 C). 69 of 69 rows ok.
* Vetting: `exosearch.v11.vet_amend.revet` (rules a-d) on the rows; scored with `known_planets.run_known_planet_recovery`
  at S* = 13.71 (it uses `passed_lc_checks_v2`).

## Required (confirmed planets): 28 of 30 recovered

Both misses are at the catalogue period with S above S*, rejected by checks the amendment did not touch:

| TIC | planet | P cat (d) | P found (d) | expected S | S found | reject | detail |
|---|---|---|---|---|---|---|---|
| 468987719 | HAT-P-43 b | 3.3326 | 3.3314 | 31.32 | 37.03 | odd_even | odd 10,738 ppm vs even 12,804 ppm: 3.17 sigma and 17.6% (thresholds 3 sigma and 10%) |
| 67512645 | TOI-6894 b | 3.3708 | 3.3718 | 44.11 | 46.69 | too_deep | box depth 13.97%, above the 5% depth cap (`vet.THRESH max_depth` 0.05) (giant planet on a very low-mass M dwarf) |

All other 28 recovered directly (no aliases), vetted, S >= 13.71, including WASP-12 b, K2-34 b, GJ 3470 b and KELT-17 b.
Per-planet rows: `known_planets_recovery_s72.csv`; summary: `known_planets_summary_s72.json`.

## TOIs (scored only): 34 of 39 recovered
Not recovered: TOI 887.01 (peak at 1.333 d, S 4.96 at that peak), TOI 2782.01 (outranked by 9.558 d, which is rejected),
TOI 2819.01 (secondary), TOI 3820.01 (secondary), TOI 5427.01 (odd_even).

## Flipped-null check
The existing gate tooling does not run the flipped null on the sector 72 known hosts (real search only, as for S71); not
computed here.

## Sector 71 diagnostic at S* = 13.71 (diagnostic only, already-seen data, not a test)
S71 known hosts re-scored with all four rules at 13.71 against the S71 frozen list (frozen at 11.82; 13.71 re-applied by
the scorer): required 28 of 28 recovered (K2-260 b, expected S 13.12, drops below 13.71), including K2-34 b, WASP-12 b,
GJ 3470 b and KELT-17 b; TOIs 44 of 51. Files: `s71_diagnostic_*_s13.71.*`.
