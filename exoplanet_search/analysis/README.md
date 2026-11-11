# Reliability cross-check

Independent statistics on the frozen exosearch v1.0.0 products. This module does not rewrite
`configs/interim_threshold_v1.json`, `configs/star_quality_v1.json`, or the interim analysis layer.
It reads the same injection and null tables and reports its own fits.

The sample files in this directory (`sample_injections.csv`, 40 rows; `sample_null.csv`, 24 inverted
rows and no scrambled rows) are a pre-freeze smoke-test extract. Every JSON this module writes on
those files has `small_sample: true`. The numbers are a check that the code runs, not a measurement
of the survey. Real runs point `--products-root` at `$EXO_DATA_ROOT` (or at `products/v1.0.0`).

```bash
cd exoplanet_search
python -m analysis \
  --injections analysis/sample_injections.csv \
  --null analysis/sample_null.csv \
  --out /tmp/exosearch-crosscheck
# real products:
python -m analysis --products-root "$EXO_DATA_ROOT" --out /tmp/exosearch-crosscheck
```

`--holdout-only` fits the held-out TICs and nothing else. The default leaves them out. A TIC is held
out when the `holdout` column is true or the TIC is listed in `holdout_tics.txt`. `--s-star` replaces
the per-sector, per-camera threshold. `--candidates` is the observed PC/FP table used for reliability;
without it, effectiveness `E` is still computed and reliability `R` is left null. `--quality` is an
optional table (`tic`, `sector`, and `beta_max` or `drop`); the cut is the registered threshold in
`configs/star_quality_v1.json` and is not re-tuned. Stars missing from that table are excluded.
`--n-boot` (default 2000) and `--seed` (default 20261110) control every bootstrap.

Outputs, each a strict JSON file (NaN and Infinity written as null) plus a PNG under `figures/`:

| file | figure | what it is |
|---|---|---|
| `completeness.json` | `completeness_gamma.png`, `completeness_map.png` | gamma CDF of detection vs expected S, and period–radius maps |
| `null_distribution.json` | `null_distributions.png` | per-star max S for each null mode, and the pooled rows |
| `tail.json` | `tail_stability.png` | GPD tail, S*, stability and mean residual life, plus a Gumbel cross-check |
| `roc.json` | `roc.png` | ROC of S, injections against null maxima |
| `reliability.json` | `reliability_effectiveness.png` | vetting effectiveness E and, when observed counts exist, reliability R |
| `crosscheck.json` |  | all of the above, with the assumption list |

Every file records `seed`, `pipeline_version`, `config_hash` and `git_hash` taken from the input rows,
and `analysis_git_hash` of the code that ran. Bootstrap streams are `seed`, `seed+1`, `seed+2` and
`seed+3` for the gamma CDF, the GPD, the Gumbel and the ROC.

## Completeness

Expected S is recomputed as `depth_inj * sqrt(n_tr) / cdpp_dur` (SCHEMA `snr_exp`). Rows with a
non-positive `n_tr` or `cdpp_dur` are left out of the fit. A detection is a period-and-epoch match
(`matched`) with `S >= S*`. An alias is not a match. The stored `detected` flag is `S >= s_min` and
is not used; `s_min` is only the v1.0.0 screening floor.

`S*` comes from `configs/interim_threshold_v1.json`: the per-camera value when that camera has one,
otherwise the sector value, otherwise the sector-70 value (the same fallback the interim layer uses
when a sector has no null of its own). A camera whose entry is null, because the interim fit had too
few vetted exceedances, uses the sector value.

The gamma CDF is the two-parameter regularized gamma distribution function, fit by an unbinned
Bernoulli likelihood. It asymptotes to 1. There is no extra plateau parameter, so a population that
never quite reaches unit completeness (aliases, edge losses) will pull the curve down. The interval
on the parameters is a star-level bootstrap: TICs are drawn with replacement and every injection on a
drawn TIC stays in the sample. The curve's band is the pointwise 2.5–97.5% range of those refits.

The period–radius map uses the fixed 8×8 logarithmic grid from `exosearch.reliability`
(0.5–13.5 d, 0.5–20 R⊕). It is not tuned on the products. Each cell carries a Wilson 95% interval
on the detection fraction, and the same cells are stored for the vetted (passed) fraction. Passing
means detected at `S*` and no rejecting check in `failed_checks` (the `low_s` token is ignored,
because the threshold was re-applied).

## Null distribution

Each row of a null product is already the maximum S of one star in one realization. Modes are the
values of the `mode` column. If that column is missing, the directory name (`invert` or `scramble`)
is used, so a file laid out as in SCHEMA.md still works.

Inverted and scrambled samples are summarized side by side, and so is the pool of all rows. The
two-sample comparison is a Kolmogorov–Smirnov test and an Anderson–Darling test. The test to quote
uses one number per TIC (the mean of that star's per-realization maxima), so several scramble
realizations of one star are not treated as independent. A pooled-row test is reported next to it
with that caveat. Scipy caps the Anderson–Darling p-value at 0.25. If either experiment is absent
the comparison is null.

## Tail and S*

For each mode, and for the pooled rows, and for each sector, excesses of per-star max S over the
sample 95th percentile are fit with a generalized Pareto distribution (scipy MLE, probability-weighted
moments if the MLE falls outside ξ ∈ (−0.9, 1.5)). This is pre-vetting. It is not the interim layer's
post-vetting fit above u = 6 on vetted inverted peaks, and the two S* values are not expected to
agree. The interim number is copied into `reference_interim_s_star` so the difference is visible.

`S*` solves `P(S > S*) = 1/N_stars` for that sector. `N_stars` is `population.kept` in the interim
config when the sector is listed there (the number of stars that passed the quality cut), not the
number of rows in a partial null extract. `--n-stars` overrides it. If the target rate is not above
the 95th percentile, the reported `S*` is the empirical quantile and `rate_in_tail` is false.

`extrapolated` is true when `S*` is larger than every observed null maximum. The threshold-stability
block refits ξ and `S*` at a grid of percentiles and records a mean-residual-life curve (mean excess
versus threshold, with a standard-error band). Stability is left unassessed when fewer than three
thresholds have at least 30 exceedances. It is flagged when `S*` changes by more than half its median
across the thresholds that do.

The Gumbel cross-check fits a Gumbel to all of the per-star maxima, not just the tail, and quotes
`N_eff = exp(location / scale)`. That is the number of i.i.d. exponential trials, with the same
scale, whose maximum has that location. It is not a count of BLS frequency samples. Its own `S*`
uses the same `1/N_stars` rate.

The bootstrap redraws TICs, recomputes the 95th percentile inside each draw, and refits. Intervals
are the 2.5 and 97.5 percentiles of the draws that produced a finite `S*`.

## ROC

Positives are injections. The score is `S` when the injection was matched and 0 otherwise, so an
unmatched injection is a miss at every positive threshold. Negatives are null maxima, pooled across
modes. AUC is the Mann–Whitney statistic (ties use midranks). The interval resamples positive TICs
and negative TICs independently. The same curve is reported overall, per sector, and per
sector×camera. Groups with an empty class are left null.

## Reliability

Per the Thompson et al. (2018, ApJS, 235, 38) form,

```
R = 1 - (N_FP,obs / N_PC,obs) * (1 - E) / E
```

`E` is the fraction of simulated nulls with `S >= S*` and checks logged that fail a rejecting vetting
check. It is computed only from the inverted and scrambled runs. `low_s` is not a vetting failure.
`ruwe_high` and `gaia_nss` are warnings in the pipeline and are not counted as failures either.
`N_FP,obs` and `N_PC,obs` come from `--candidates` when that file exists. Recognized columns are
`disposition` (or `label`), `period` (or `P`), `rp_earth` (or `rp`), and optionally `sector`, `tic`,
`holdout`. Planet-candidate labels are `PC`, `planet`, `candidate`, `KP`, `CP`; false-positive labels
are `FP`, `FA`, `EB`.

The uncertainty is a Wilson 95% interval on `E`, pushed through the formula. Observed counts are
treated as fixed. Because `(1-E)/E` decreases with `E`, the lower endpoint of `E` becomes the lower
endpoint of `R`. `R` is null when the candidates file is absent, when a bin has no planet candidates,
or when `E` is 0. The formula can come out negative; `R_floored_at_0` is also stored. Bins with no
simulated false alarm above `S*` are omitted rather than filled with a made-up `E`.

## What is deliberately not done

`analysis/make_reliability_json.py` is not installed. `exosearch.reliability` treats that filename as
the signal to replace `site/data/reliability.json`, and the Methods page reads a different schema
(binned counts, the interim threshold, per-check rates). This cross-check keeps its own files so that
hook is left to the interim summary.

## Assumptions that the schema leaves open

These are the choices above, collected in one place:

- The middle dot in the SCHEMA definition of `snr_exp` is read as multiplication.
- The tail requested here is the 95th percentile of pre-vetting per-star max S. The interim S* is a
  different estimator (post-vetting, fixed u = 6, vetted inverted peaks only).
- On a partial null extract, `N_stars` is still the searched population in the interim config, when
  that sector is present. `S*` then usually lies far beyond the largest value in the extract, and the
  extrapolation flag says so.
- Scramble realizations are separate rows. The quoted two-sample test collapses them to one number
  per star first.
- `passed_all` on a null row is false whenever `S < s_min`. Effectiveness ignores that flag and
  recomputes vetting from the failure columns.
- The candidates table is not in SCHEMA.md. The column names above are the ones this module accepts.
- A camera with a null interim `s_star` uses the sector value, as that file's note says.
