# Reliability products: per-row schema (pipeline v1.0.0)

Every row carries `pipeline_version`, `config_hash` (sha256 of every setting that can change a result,
see `exosearch/version.py`) and `git_hash` (commit of the code that produced it; `-dirty` = uncommitted).
Products live under `$EXO_DATA_ROOT/products/v<version>/` (parquet). `sample_injections.csv` and
`sample_null.csv` in this folder are small real extracts (from the pre-freeze smoke test, so their
`git_hash` ends in `-dirty`) for coding against.

Units: days for times/periods/durations (BTJD = BJD − 2457000), relative flux for depths, R⊕ for planet
radii, R☉/M☉ for stars, solar units for densities.

## Detection statistic (all products)

| column | meaning |
|---|---|
| `S` | detection statistic at the star's best peak: `depth_box · sqrt(n_tr) / cdpp_d`. Thresholded (`S >= s_min`). |
| `S_scan` | same statistic as the scan saw it (astropy BLS `snr` objective with σ_D = CDPP_D·sqrt(D/bin), fractional N_tr). Peak = argmax over (P, D, T0). |
| `depth_box` | depth reported by astropy BLS `fast` method at the S-scan peak. Enters `S` only. Can exceed the depth of any box at the reported (P, D, T0) when D/P is large: see *Depth columns* |
| `n_tr` | coverage-weighted number of transits with data: Σ_epochs min(1, n_in / n_expected) on the 10-min binned flux |
| `cdpp_d` | CDPP at the peak duration D: 1.4826·MAD of the running mean of detrended residuals over D (unbinned 2-min, windows need ≥½ the expected cadences) |
| `beta_d` | red-noise factor at D: `cdpp_d · sqrt(D/cadence) / sigma_w` (1 = white) |
| `sigma_w` | white noise per 2-min cadence: 1.4826·MAD of detrended residuals |
| `sde` | SDE of the max-over-D S periodogram. Logged only, never thresholded. |
| `snr` | white-noise BLS depth/depth_err (astropy compute_stats). Logged only. |
| `s_min` | v1.0.0 in-pipeline screening floor when the row was made (`configs/detection_threshold.json`). Not the detection threshold: re-score at S* from `configs/interim_threshold_v1.json` |
| `camera`, `ccd` | TESS camera/CCD from the light-curve header (logged; S* pools cameras) |

## Depth columns (audit 2026-11-10, interim layer)

All depths are fractional flux decrements (1e-6 = 1 ppm). Two estimators exist at the same peak (P, D, T0):

| column | where | meaning | consumers |
|---|---|---|---|
| `depth_box` | results DB (v1.0.0 rows), nulls | BLS `fast`-method depth at the S peak (10-min binned flux, phase-binned) | **`S` only** |
| `depth`, `depth_err` | results DB, nulls, labels | astropy `compute_stats` box depth at (P, D, T0): out-of-transit minus in-transit weighted mean, binned flux | `snr`; screens `too_deep` (> 0.05), `no_dip` (≤ 0), `rp_too_large` via `rp_earth`; deep checks `rp_gt_2rjup`, density fit, Gaia dilution test; site `depth_ppm`; rejected-list filter |
| `depth_odd`, `depth_even`, `secondary_depth` | results DB | `compute_stats` odd/even/phase-0.5 depths | `odd_even`, `secondary` screens |
| `rp_earth` | results DB, nulls | `sqrt(depth)·R★·109.1` | `rp_too_large` (> 25), candidate pages |
| `depth_inj` | injections | injected model depth, 1 − min(model) | `snr_exp`, completeness bins, "allowed" range (`depth_inj < 0.05`, `rp_inj < 25`) |
| `depth_rec` | injections | `compute_stats` depth at the recovered peak (same as `depth`) | injection screens and deep checks; not used for matching (P, T0 only) |
| `depth_box_implied` | reproduction export only | `S·cdpp_d_rec/sqrt(n_tr_rec)`: the `depth_box` the v1.0.0 injection rows do not store | audit |
| `depth_ppm` | site JSON, labels | `depth·1e6` (labels: TOI catalogue "Depth (ppm)") | display |

Findings (`python -m exosearch.depth_audit`): Rp, all vetoes, the site and injection matching use `depth` (correct).
`depth_box`/`depth` has median 1.01 on the S70 flipped null (p99 1.71; 0.9% of rows > 2× or `depth` ≤ 0), 1.09 on the
S105 real search (p99 2.35; 1.5%) and 1.000 on matched injections (p99 1.06; none > 2×). The excess grows with D/P and
with S: it inflates non-transit peaks, not transits. **`S_exact = depth·sqrt(n_tr)/cdpp_d`** (injections: `depth_rec`,
`n_tr_rec`, `cdpp_d_rec`) is the interim-layer statistic. The peak is not re-chosen, so per-star `S_exact` is a lower
bound on the max of that statistic. v1.0.0 products are unchanged; v1.1 computes S from the exact box depth in-search.

## Planned v1.1 provenance columns (all products; not in v1.0.0 rows)

| column | meaning |
|---|---|
| `host_id` | opaque worker id (`h1`, `h2`, `h3`); the mapping to machines is kept outside the repo |
| `arch` | CPU architecture (`x86_64`, `arm64`) |
| `lib_versions` | numpy / scipy / astropy / BLAS versions from the lockfile, as one string |
| `statistic_id` | `S_cpu_f64` for the registered statistic; any GPU / float32 BLS is a separate id with its own nulls |
| `prewhiten_n`, `prewhiten_freqs`, `prewhiten_comb` | prewhitening terms removed, their frequencies (c/d), peaks kept by the comb guard |

## Injections — `injections/injection_v1/part-*.parquet`

One row per (sector, run, TIC); config `configs/injection_v1.json`.

| column | meaning |
|---|---|
| `key` | `sector:run:tic` (unique) · `job` = batch of runs on one light curve |
| `tic`, `sector`, `sectors`, `run`, `seed` | host, sector(s) searched, run index, RNG seed of this injection |
| `holdout` | TIC in the held-out 20% (`holdout_tics.txt`) |
| `tmag`, `rstar`, `mstar`, `teff`, `logg` | TIC v8.2 host parameters used |
| `u1`, `u2` | quadratic LD (Claret 2017 TESS, nearest Teff/logg) |
| `P_inj`, `t0_inj`, `rp_inj`, `b_inj`, `k_inj`, `aR_inj` | injected period, epoch, radius (R⊕), impact parameter, Rp/R★, a/R★ |
| `depth_inj`, `dur_inj` | injected model depth (1 − min flux), T14 |
| `baseline` | span of the quality-masked data (P is drawn up to `baseline/2`) |
| `n_tr` | transits of the injected ephemeris with ≥½ their cadences after the quality mask |
| `cdpp_dur` | CDPP at D = `dur_inj` (same estimator as `cdpp_d`, injected transits masked) |
| `snr_exp` | expected S: `depth_inj / cdpp_dur · sqrt(n_tr)` |
| `P_rec`, `t0_rec`, `dur_rec`, `depth_rec` | best peak found |
| `S`, `S_scan`, `n_tr_rec`, `cdpp_d_rec`, `beta_d`, `sigma_w`, `sde`, `snr`, `s_min` | search statistics at the best peak |
| `n_transits_rec` | transits with ≥1 binned point at the recovered ephemeris (the `few_transits` screen uses it) |
| `matched` | `|P_rec/P_inj − 1| < 0.01` and `|T0_rec − T0_inj| < 0.5·dur_inj` (mod P) |
| `alias` | `2x`, `1/2x`, `3x`, `1/3x` when P_rec is that multiple of P_inj (1%); else null. Not counted as matched. |
| `detected` | `matched and S >= s_min` |
| `passed_all` | `detected` and no rejecting check failed |
| `checks_logged` | deep checks were evaluated (`matched` or `S >= 6`) |
| `fail_<check>` | per-check outcome for screens `few_transits odd_even secondary too_deep too_long variable rp_too_large no_dip common_period` and deep `v_shape secondary_scan rp_gt_2rjup density_mismatch ruwe_high gaia_nss` (null when not evaluated) |
| `failed_checks`, `warnings` | comma lists (rejecting / warning-only) |
| `v_ratio`, `sec_scan_sigma`, `rho_z`, `fit_ar`, `fit_b` | deep-check metrics |
| `status`, `error` | `ok` / `unphysical` (a/R★ ≤ 1+k) / `error` |
| `runtime_s`, `load_s`, `detrend_s`, `search_s` | profiling |

Pixel-level checks (centroid offset, dip off target, known EB, Gaia neighbours) cannot be exercised by
light-curve injections; hosts with TOIs or catalogue EBs are excluded, so these count as passed.

## Null tests — `nulltests/<invert|scramble_rN>/null_*_maxS.parquet`

One row per searched star per null realization (the per-star max-S table).

| column | meaning |
|---|---|
| `key`, `tic`, `sector`, `camera`, `ccd`, `tmag`, `holdout` | identity |
| `mode`, `realization`, `seed` | `invert` (flares masked first, flux → 2 − flux; one realization) or `scramble` (~1-day cadence blocks permuted; realizations with their own seeds) |
| `period`, `t0`, `duration`, `S`, `S_scan`, `n_tr`, `cdpp_d`, `beta_d`, `sigma_w`, `depth_box`, `sde`, `snr`, `depth`, `n_transits`, `rp_earth`, `radius`, `s_min` | best null peak |
| `flags`, `passed` | per-star screens (`low_s` here is relative to the screening floor `s_min`, not S*) |
| `checks_logged` | deep checks evaluated (`S >= 6` or passed screens) |
| `fail_<check>` | as for injections, plus `fail_common_period` from the null run's own common-period list |
| `deep_flags`, `rho_z`, `v_ratio`, `v_sigma`, `sec_scan_sigma`, `rp_rjup` | deep-check outputs |
| `passed_all` | passed screens and no rejecting deep flag |

## Labelled signals — `labelled/labelled_signals.parquet` + `labelled/lc_checks/`

Our best peak matched to a KP/CP TOI (`label=planet`) or an FP/FA TOI or TESS-EB period (`label=fp`),
with every search column from the results DB and the light-curve deep checks. Used for per-check
true/false rejection rates; tune only on `holdout == False`.

## Hook for analysis/

Put analysis code here (gamma-CDF completeness with bootstrap, GPD/Gumbel tail fits of null max S,
ROC, vetting effectiveness E, reliability). The site reads `site/data/reliability.json`; if
`analysis/` provides `make_reliability_json.py` writing that file, `scripts/publish_site.sh` picks it up
unchanged. Until then the harness writes provisional binned summaries there.

## Star quality (interim layer) — `analysis_layer/star_quality_s<sector>.parquet`

One row per searched star (real light curve, no injection or inversion), measured after the frozen detrend
with `exosearch.bls.noise_profile`. Cut: `configs/star_quality_v1.json` (frozen; registered after unblinding):
drop if `beta_max > 3.0` or `status != ok`.

| column | meaning |
|---|---|
| `key`, `tic`, `sector`, `holdout`, `tmag`, `camera`, `ccd` | identity |
| `sigma_w`, `n_cad` | white noise per cadence, cadences after the quality mask |
| `cdpp_<D>`, `beta_<D>` | CDPP and red-noise factor at each trial duration D (days) |
| `beta_max` | max over D of `beta_<D>` |
| `status` | `ok` or `error` |

## Reproduction export — `data/reproduction/` (in the repo)

`exosearch/export_repro.py` writes gzipped CSVs of the tables above (null `invert` max-S, injections, star quality)
with the profiling, error-text and file-path columns removed, plus `MANIFEST.json` (snapshot time, row counts,
sha256, pipeline/config/git hashes present). Derived interim-layer columns, listed in `MANIFEST.json`
`derived_columns`: `S_exact` (null, injections) and `depth_box_implied` (injections). Vetted null (interim layer): `status == ok`, `flags` minus `low_s`
empty, `fail_common_period` false, no `deep_flags` in the rejecting set `v_shape secondary_scan rp_gt_2rjup
density_mismatch centroid_offset dip_off_target known_eb` (`exosearch/deepvet.py` REJECT), `checks_logged` true.
