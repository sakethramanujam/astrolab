# Registration: choosing the red-noise (beta_D) treatment for v1.1

Interim analysis layer over frozen v1.0.0. Registered 2026-11-10, before any v1.1 run.

## Disclosure

The four options below were first compared on the S70 flipped null and the S70 injections
(`python -m exosearch.stat_options`, `docs/v1.1_odd_even_secondary.md`). Option c (floor at 1 − 2·sigma_beta) was
written down as the preferred option **after** those S* effects had been seen (S* 11.13 → 10.63 on the same sample).
No candidate (real-search row) entered that comparison or enters the rule below. This file replaces that preference
with a fixed rule; if the rule picks another option, that option is used.

## Options (beta' replaces beta_D; S' = S·beta/beta'; identical on nulls, injections and the search)

| id | treatment |
|---|---|
| a | none |
| b | beta' = max(1, beta) |
| c | beta' = max(beta, 1 − 2·sigma_beta), sigma_beta = 1/sqrt(2·N_w), N_w = n_cad·cadence/D |
| d | empirical Bayes: normal prior per D (median and between-star variance of kept stars' beta_D), posterior mean |

## Selection rule (nulls and injections only)

Evaluated per sector with its own nulls (S70 and S106), tuning stars, star-quality cut applied, statistic S_exact:

1. **Null consistency (gate).** S*_flip and S*_scr (scrambled realization r0 of the v1.0 nulls; in v1.1, the
   registered scrambled set) must agree: the 95% star-level bootstrap CI of S*_flip − S*_scr contains 0, in every
   sector with both nulls. Options failing the gate in any sector are out.
2. **Completeness cost.** Among options passing the gate, choose the highest vetted completeness at expected S 15–25
   (allowed injection range), each option at its own S* = max(S*_flip, S*_scr).
3. **Tie-break.** Differences within 2 percentage points (binomial SE at these counts) count as ties; ties go to the
   option that modifies the fewest null peaks by more than 1%, i.e. the smallest intervention.
4. If no option passes the gate, use b (most conservative: the error bar is never shrunk) and report the failure.

Candidates, candidate counts and TOI recovery are not inputs to the choice. The rule is applied once, on the v1.0
interim nulls when scrambled r0 is complete for S70 and S106; the chosen treatment is then implemented inside the
v1.1 statistic. The other three options are reported as sensitivity rows next to every headline number (S*,
completeness, false-alarm rate, candidate counts), with and without the quality cut.

## Amendment 1 — 2026-11-11T06:56:32Z (UTC). Predates any scrambled r0 output.

Records the rule agreed on 2026-11-10 at 20:07 PST, which had not been written into this file. When this amendment
was committed, no scrambled r0 output existed: the v1.0 null directory held only the flipped null (`invert`, S106
still finishing) and an unrelated earlier test directory (`scramble_r9`, 16:50 PST Nov 10). There was no
`scramble_r0` directory, part file or table. The selection code (`exosearch/v11/beta_select.py`) had been run once,
as a code check with the flipped null used as both inputs; no scrambled data entered that run. The amendment
replaces rules 1 and 3 above:

1. **Gate (equivalence, TOST).** An option passes only if, in every sector with both nulls, the **90%** star-level
   bootstrap CI (paired over stars in both nulls, 2,000 resamples, seed 20261110) of dS* = S*_flip − S*_scr lies
   **entirely inside [−0.5, +0.5]**. The old check (95% CI contains 0) is no longer the gate. It is reported
   descriptively.
3. **Ties.** An option is tied with the leading option (highest vetted completeness at expected S 15–25) when the
   paired-bootstrap 95% CI (over the same injections) of their vetted-completeness difference includes 0. Among the
   leader and the options tied with it, the first in the fixed order **a, d, c, b** (smallest intervention first) is
   chosen. The "fewest null peaks modified by more than 1%" tie rule is removed. Those counts are still reported.

Rules 2 and 4 are unchanged. The rule is judged on nulls and injections only: candidates, candidate counts and TOI
recovery are not inputs. If no option passes the gate, b is recorded per rule 4 and the v1.1 production run is not
launched until the owner decides.
