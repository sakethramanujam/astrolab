# v1.2 S*_cal rulings, 2026-11-20

- **Status:** RULINGS that fill gaps in `docs/registration_v1.2.md` (tag `v1.2-registered`) and the 2026-11-19 h2 amendment with Addendum A. They were found before S*_cal was computed, and are committed before S*_cal is computed. Approved by the project owner on 2026-11-20; Scientist and RF Engineer review the driver before any compute.
- **Provenance:** Each ruling comes from registered text or from code at the pinned commit `bc05650`; the source is cited under each ruling.
- **Pre-computation state:** No S*_cal, eligibility count, padding draw, cut retention, or S69 quantity has been computed when this note is committed.

## 1. Budget

- **Ruling:** Expected 1 false alarm per sector. `p_s = 1 / N_s`, where `N_s` is the number of searched stars in sector `s` passing the star-quality cut, holdout included (the flipped-null population), as `exosearch/v11/beta_rows.py` defines it (`beta_frames`: `n_alpha`; `sector_gate`: `p = 1/N`).
- **Source:** Registration Section 6.2 (`B = alpha_FAP x N_runs`) and the v1.1 convention in `exosearch/v11/beta_rows.py` / `interim.fit_sector` (`p = 1/n_stars`).

## 2. Pooling

- **Ruling:** One `S*_cal` is computed over S59-S64. The survivor budget is summed over sectors, each at its own `p_s`: `B = sum_s p_s x n_null,s` (expected to be about 12-14). Per-sector `S*` values, each with a bootstrap 68% interval, are report-only and never replace `S*_cal`.
- **Source:** Registration Section 0A.1 (one calibration set S59-S64) and Section 6.2; selector `exosearch/sstar_grid.select_sstar` via `exosearch/v12_settings.calibrate_sstar_v12`, which take the budget from the caller.

## 3. Nulls and vetting

- **Ruling:** Invert and scramble null rows are pooled and vetted with the frozen pre-calibration v1.2 vetter at `bc05650` (`exosearch/vet_final.final_verdict` under the frozen `configs/v1.2.json`, every candidate cut hard, `too_deep` as a flag). The driver computes `passed_lc_checks_v12` from it. Two diagnostics are reported and never used as `S*_cal`: the unvetted `S*`, and the v1.1-style max of separate invert and scramble values (`beta_rows.sector_gate`).
- **Source:** Registration Section 6.2 ("pooled across the paired flipped + scrambled null realizations ... after applying the same-population guard and passed_lc_checks_v12"), Section 0B.2, and the frozen settings.

## 4. Population

- **Ruling:** The S59-S64 scrambled-null exclusion set is rebuilt from the run manifest (targets with an invert task minus targets with a scramble task, per sector). It must equal the registered rule: known transiting hosts plus eclipsing binaries from the frozen TOI table and the Kostov et al. 2025 list. On any mismatch, `S*_cal` stops and the mismatch is reported. If it matches, the set is recorded with its SHA-256, labelled "reconstructed after the run". The rows carry the S70 exclusion-file stamp only (`scramble_exclude_file_70`), not an S59-S64 stamp; this is noted, not corrected. The null sample is the non-holdout stars that pass the star-quality cut and have real, invert, and scramble rows all present; the same-population guard (`exosearch/v11/registration_checks.check_registration_population`) is run on it.
- **Source:** Registration Sections 7.1 and 0B.1; `exosearch/v11/beta_rows.py` (`beta_frames`).

## 5. Eligibility (`n_real_eligible`)

- **Ruling:** As Section 0A.2.1: `S*_cal <= S_exp <= 2 S*_cal`, with `S_exp` from catalogue depth and duration and `expected_s_scale = 0.66` (the `known_planets` eligibility formula, `derive_expected`), at least 3 transits, non-holdout stars, unique planets, and CP/KP TOIs counted alongside confirmed planets. Gate disjointness is enforced when the S69 list is frozen (Step 5). Padding hosts are kept off the gate using the S69 target list from the run manifest.
- **Source:** Registration Sections 0A.1, 0A.2.1 items 1 and 5, 2.2, and 3.1.

## 6. Record

- **Ruling:** `S*_cal` is committed with SHA-256 sidecars before the first padding draw, in the order of Addendum A.
- **Source:** Addendum A item 1.

**Not changed:** No registered rule, threshold, cut, seed, or gate is changed by this note.
