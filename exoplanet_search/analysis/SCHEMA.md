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
| `depth_box` | BLS box depth δ̂ at the peak |
| `n_tr` | coverage-weighted number of transits with data: Σ_epochs min(1, n_in / n_expected) on the 10-min binned flux |
| `cdpp_d` | CDPP at the peak duration D: 1.4826·MAD of the running mean of detrended residuals over D (unbinned 2-min, windows need ≥½ the expected cadences) |
| `beta_d` | red-noise factor at D: `cdpp_d · sqrt(D/cadence) / sigma_w` (1 = white) |
| `sigma_w` | white noise per 2-min cadence: 1.4826·MAD of detrended residuals |
| `sde` | SDE of the max-over-D S periodogram. Logged only, never thresholded. |
| `snr` | white-noise BLS depth/depth_err (astropy compute_stats). Logged only. |
| `s_min` | threshold applied when the row was made (`configs/detection_threshold.json`; `provisional` there) |
| `camera`, `ccd` | TESS camera/CCD from the light-curve header (S* is per sector and camera) |

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
| `flags`, `passed` | per-star screens (incl. `low_s` at the provisional S*) |
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
