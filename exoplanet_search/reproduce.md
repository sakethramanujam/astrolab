# Reproducing the v1.0.0 reliability results

Pipeline version `1.0.0` (git tag `exosearch-v1.0.0`). Every product row stores `pipeline_version`,
`config_hash` and `git_hash`; `python -m exosearch.version` prints the hash, the full frozen config and the
package versions. Results from code with a different `config_hash` are not comparable.

## Environment

```bash
cd exoplanet_search
python3.12 -m venv .venv && .venv/bin/pip install -r requirements.txt
export EXO_DATA_ROOT=/path/to/data          # light curves, catalogues, results, products
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
.venv/bin/python -m exosearch.version > provenance.json
.venv/bin/python -m pytest -q tests
```

Package versions used for the published products are in each product's `provenance.json`.

## Fixed seeds and splits

| what | seed / rule |
|---|---|
| held-out stars | `sha256("20261110:<tic>")/2**256 < 0.2` (`exosearch/holdout.py`, list in `holdout_tics.txt`) |
| injection hosts | `configs/injection_v1.json` seed 70106; per injection `sha256("70106:<sector>:<run>:<tic>")[:4]` |
| scramble realizations | `sha256("20261110:scramble:<r>:<sector>:<tic>")[:4]` |
| star-quality cut | `configs/star_quality_v1.json` (beta_max > 3; registered before evaluation) |
| detection threshold | `configs/interim_threshold_v1.json` (interim GPD fit, bootstrap seed 20261110, 2000 resamples); `configs/detection_threshold.json` holds only the v1.0.0 in-pipeline screening floor |
| white-noise calibration | `sha256("20261110:white:<sector>:<tic>")[:4]`, 300 tuning stars drawn with seed 20261110 |

## Steps

```bash
# 1. search the sectors (downloads SPOC 2-min light curves from MAST)
.venv/bin/python run_search.py --sector 70  --n-targets 0 --workers 12
.venv/bin/python run_search.py --sector 106 --n-targets 0 --workers 12
# keep the light curves (the continuous runner prunes FITS): products need them
mkdir -p $EXO_DATA_ROOT/keep/s0106/lc && cp -al $EXO_DATA_ROOT/s0106/lc/. $EXO_DATA_ROOT/keep/s0106/lc/

# 2. labelled signals (KP/CP vs FP/FA/EB) + light-curve deep checks
.venv/bin/python -m exosearch.labels --sectors 70,106,107 --workers 12

# 3. null tests: inversion (main) + scramble realizations (limited with one sector)
.venv/bin/python -m exosearch.nulltest --mode invert   --sectors 70,106 --workers 12
for r in 0 1 2; do .venv/bin/python -m exosearch.nulltest --mode scramble --realization $r --sectors 70,106 --workers 12; done

# 4. injection-recovery, ramped: first batch, then more runs (one injection per light curve per run)
.venv/bin/python -m exosearch.inject --config configs/injection_v1.json --runs 0   --workers 12
.venv/bin/python -m exosearch.inject --config configs/injection_v1.json --runs 1-11 --workers 12

# 5. missed known TOIs
.venv/bin/python -m exosearch.missed --sector 70

# 6. pixel-level deep vetting of candidates (TPF difference images, Gaia, EB catalogues, density, RUWE)
.venv/bin/python -m exosearch.deepvet --sector 70 --only-review
```

All batch commands are resumable (finished rows are skipped) and run fine under `nice -n 19 ionice -c3`.
On a shared machine they register worker claims in `/tmp/exosearch-claims/` so the continuous runner
(`scripts/continuous.sh`) shrinks to keep the total under `COMBINED_CAP`.

## Archive layout (for Zenodo; not uploaded)

```
products/v1.0.0/
  injections/injection_v1/{provenance.json, part-*.parquet}
  nulltests/{invert, scramble_r0, ...}/{provenance.json, parts/, null_*_maxS.parquet}
  labelled/{labelled_signals.parquet, lc_checks/}
  missed_tois/missed_tois_s0070.{csv,json}
  summary/ (site JSON + figures)
```

Plus the git tag, `analysis/SCHEMA.md`, `configs/`, `holdout_tics.txt` and this file.
Raw TESS data are not archived; they are public at MAST (SPOC 2-min light curves, sectors 70/106/107).
