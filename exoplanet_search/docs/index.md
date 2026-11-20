# exoplanet_search docs index

## v1.2 (current)

- [registration_v1.2.md](registration_v1.2.md) — v1.2 registration (tag `v1.2-registered`, code pinned at `bc05650`), unchanged, followed by dated sections, each verbatim from its tagged commit:
  - Amendment 2026-11-19: single host h2 (tag `v1.2-amend-h2-2026-11-19`) and Addendum A (tag `v1.2-amend-h2-2026-11-19-A`): order of remaining steps; 18 workers at normal priority, speed only.
  - Catalogue cross-check spec 2026-11-19 (tag `v1.2-xcheck-spec-2026-11-19`): read-only diagnostic against published ephemerides (103 of 376 scored planets passing all criteria).
  - Catalogue cross-check amendment 2026-11-20 (tag `v1.2-xcheck-amend-2026-11-20`, commit `e4f0f45`): 192 of 297 scored planets passing all criteria.
  - S*_cal rulings 2026-11-20 (tag `v1.2-sstar-rulings-2026-11-20`): fills registration gaps before `S*_cal` is computed.

## Open steps

In order:

1. `S*_cal` and `n_real_eligible`.
2. Injection padding (seed `20261118`, hard cap 2000 draws).
3. Vetter cut lock on S59-S64.
4. S69 gate (`S*_gate`, frozen eligibility list, zero known-planet losses).

Nothing after step 1 has started.

## v1.1 and earlier

### Files

- [interim_headlines.md](interim_headlines.md) — Interim headlines notes (see file).
- [interim_lane_cap.md](interim_lane_cap.md) — Interim lane cap notes (see file).
- [v1.0_disclosure_ties_and_numba_cache.md](v1.0_disclosure_ties_and_numba_cache.md) — v1.0 disclosure on ties and numba cache.
- [v1.1_implementation.md](v1.1_implementation.md) — v1.1 implementation notes.
- [v1.1_odd_even_secondary.md](v1.1_odd_even_secondary.md) — v1.1 odd/even and secondary check notes.
- [v1.1_plan.md](v1.1_plan.md) — v1.1 plan, followed by the 2026-11-11 v1.1 amendment and the 2026-11-11 scrambled-null amendment as dated sections.
- [v1.1_prewhitening_amendment.md](v1.1_prewhitening_amendment.md) — v1.1 prewhitening amendment.
- [v1.1_remaining.md](v1.1_remaining.md) — v1.1 remaining items.
- [v1.1_validation.md](v1.1_validation.md) — v1.1 validation notes.

### Folders

- [figures/](figures/) — Documentation figures.
- [host_gate/](host_gate/) — Host gate notes (see folder).
- [preregistration/](preregistration/) — Preregistration notes (see folder).
- [v1.1_checks/](v1.1_checks/) — v1.1 check notes (see folder).
