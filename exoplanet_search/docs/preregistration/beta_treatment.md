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
