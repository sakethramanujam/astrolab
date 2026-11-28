# Astrolab

TESS transit search. v1.2 searches sectors 59–64, with sector 69 nulls, using the code pinned at `bc05650`.

## Install

```bash
cd exoplanet_search
python3 -m venv .venv
.venv/bin/pip install -r requirements.txt
```

## Run

```bash
cd exoplanet_search
.venv/bin/python -m exosearch.v11.run --help
.venv/bin/python -m exosearch.v12_sstar_cal --help
```

`python -m exosearch.v11.run` is the search. `python -m exosearch.v12_sstar_cal` is the S*_cal driver. It has not been run.

## Test

```bash
cd exoplanet_search
.venv/bin/python -m pytest -q tests
```

## v1.2 result

The registered h2 run (code pin `bc05650`) finished 83,130 of 83,130 star jobs and 235,168 ok rows out of 235,171 (3 errors) for S59–S64 real, invert, and scramble, plus S69 invert and scramble. The catalogue cross-check scored 297 of 420 planets in scope, and 192 passed every amended criterion.

## Failed runs

`failed_4` and `failed_5` are run outputs on the owner's machine and are not in this repository. They stay on record. No SHA-256 for either file is stored here; point at those two outputs by checksum when the hashes are copied in.
