# Claret 2017 TESS quadratic limb-darkening table

Source: Claret, A. 2017, *Astronomy & Astrophysics*, 600, A30.
VizieR catalog `J/A+A/600/A30`, table `tableab`.

`claret2017_tess_quadratic.csv` is the cache file `exosearch.stellar.claret_table` writes
(`catalogs/claret2017_tess_quadratic.csv`): columns `logg`, `Teff`, `Z`, `aLSM`, `bLSM`,
`Type`, `Mod`. It was built from the VizieR TSV by the same steps as `claret_table`
(drop `#` lines, header on the `recno` row, skip the two following rows, keep those
columns, strip `Type` and `Mod`). `Z` is 0 only. `quad_ld` still selects the r-method
PHOENIX-COND grid and the PHOENIX-DRIFT rows below that grid's minimum Teff. A
read-back of this file matches that selection with maximum absolute difference 0.

Fetched 2026-11-27 from the CFA mirror:

`https://vizier.cfa.harvard.edu/viz-bin/asu-tsv?-source=J/A%2BA/600/A30/tableab&-out.all&-out.max=999999`

The CDS mirror `https://vizier.cds.unistra.fr/viz-bin/asu-tsv` returned HTTP 404 for the same query on that date.

SHA-256 of `claret2017_tess_quadratic.csv` (also in `claret2017_tess_quadratic.csv.sha256`):

`f39c593b6f8dd94863cef163134332940a58d242ae69e3aca03c93f65e1ca33b`

`exosearch/stellar.py` is in the v1.1 search-path hash list
(`docs/preregistration/v1.1_early_nulls_2026-11-11/searchpath_hash.py`) and is not modified.
`claret_table` already reads a local cache and calls VizieR only when that file is absent.
Tests that need coefficients monkeypatch `stellar.claret_table` to read this copy, so
`quad_ld`'s nearest-grid result is the published table.
