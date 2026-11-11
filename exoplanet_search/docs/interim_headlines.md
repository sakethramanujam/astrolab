# Interim headline numbers, with and without the star-quality cut

Interim analysis layer over frozen v1.0.0; statistic S_exact (exact box depth, `analysis/SCHEMA.md` *Depth columns*).
Snapshot 2026-11-10 ~20:15 PST from `python -m exosearch.cut_sensitivity` (S70 flipped null still growing; S* is
provisional until the flipped and scrambled nulls are complete). "No cut" = every measured searched star, with its own
alpha = 1/N_stars.

| Sector 70 | Stars | S* (95%) | Expected post-vetting false alarms / sector | Found / vetted at S_exp 15–25 |
|---|---|---|---|---|
| quality cut (beta_max ≤ 3) | 10,244 | **11.07** (9.32–13.32; u 7.0) | 0.88 at its S* (fit on 7,148 null stars; 1 by construction at fit time) | 77% / 64% (n 195) |
| no cut | 12,797 | 11.17 (9.54–13.10; u 6.0) | 1.00 at its S*; 1.08 at the with-cut S* | 77% / 56% (n 232) |

Sector 105 (not calibrated: no own null; counted against the S70 thresholds), real search, 7,982 rows so far:
557 pass the v1.0.0 screens; 96 of those have S_exact ≥ 11.17 (no cut), and 96 have S_exact ≥ 11.07. The with-cut
count needs S105 star quality, which is queued (`exo-quality105`, 3 workers, starts when the batch cap has room).
None of these is a detection; they are not calibrated.

Completeness and S* for the beta-treatment options, with and without the cut: `docs/v1.1_odd_even_secondary.md`.
