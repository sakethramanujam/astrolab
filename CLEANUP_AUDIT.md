# Cleanup audit (not applied, not pushed)

Status: proposal only. No files were deleted. No pull request was opened. Branch `cursor/dead-code-cleanup-7a81` exists only in this working copy.

The remote does not have `v1.2-flux-cut`, and it does not have commits `8f48f07` or `7b300d3`. Those five amendment commits are still only on the owner's machine. This audit is against a local consolidation of what *is* on the remote. Rebase onto the unpushed commits before applying anything below. Several proposals touch `v12_sstar_cal.py` and the vetter, which those commits are likely to rewrite.

## How this tree was built

Local commits on `cursor/dead-code-cleanup-7a81` (not pushed):

| Commit | What |
| --- | --- |
| `6a868bc` | Merge `origin/v1.2-vetter-unify` (`bc05650`) into `origin/exoplanet-search-v1.2` |
| `446bbe5` | Merge `origin/v1.2-amendment-h2` |
| `4cc67d3` | Add `LICENSE` from `origin/master` |
| `723a696` | Merge later `status_v12.json` snapshots from `origin/exoplanet-search` (through `3e6cd64`, 2026-11-27) |

Why this shape:

- `origin/exoplanet-search-v1.2` (`bb0cd54`) has the later S*_cal driver, the deviation and common-period ruling pasted into `exoplanet_search/docs/registration_v1.2.md`, the cross-check outputs with SHA-256s, and the site publisher. It does **not** have the grid, unified vetter, or `configs/v1.2.json`.
- `origin/v1.2-vetter-unify` (`bc05650`, tag `v1.2-step7b-cfg-guard`) is the code the driver pins (`PINNED_CODE_COMMIT`). Those modules exist on no other remote branch.
- `origin/v1.2-amendment-h2` holds the tagged pre-registration, the h2 amendment and Addendum A, the cross-check spec and amendment, the S*_cal rulings, and an **earlier** driver (`b2990d5`). The later driver was kept. The earlier driver remains in that branch's history.
- `origin/master` contributed only `LICENSE`. The teaching-lab files are already on every exoplanet branch.
- `origin/exoplanet-search` after `f1f0a20` is only `site/data/status_v12.json` refreshes.

Conflict choices that are not cleanups:

- Kept `exoplanet_search/exosearch/v12_sstar_cal.py` from `exoplanet-search-v1.2` (SHA-256 prefix `6d69becc29`). The amendment-h2 copy (`839456a638`) is the pre-deviation driver.
- Kept `exoplanet_search/docs/registration_v1.2.md` from `exoplanet-search-v1.2`. It is the registered text plus verbatim dated sections (h2 amendment, both cross-check notes, S*_cal rulings, Kostov deviation, TOI deviation, common-period ruling).
- Also kept the original separate files from `v1.2-amendment-h2`, and restored `preregistration/v1.1_amendment_2026-11-11.md` and `v1.1_scramble_amendment.md`, which `0e18dca` had folded into `v1.1_plan.md` and deleted.
- Kept the vetter-unify AI-use entries and appended the registration-draft AI-use entries. `0e18dca` had deleted the log before the vetter steps existed.

The upstream of this branch was unset so a push cannot land on `exoplanet-search-v1.2`.

## Inventory

`HEAD` has **2131** tracked files (the prior ~2111 count was one branch, before this merge added the vetter modules, amendment documents, and `LICENSE`).

| Kind | Count |
| --- | --- |
| `.json` | 953 |
| `.png` | 919 |
| `.py` | 114 |
| `.md` | 54 |
| notebooks | 7 |
| other | the rest |

Generated per-star site products:

| Set | Files | In tag `exosearch-v1.0.0` | Only on current branches | Bytes |
| --- | ---: | ---: | ---: | ---: |
| `site/plots/**/*.png` | 908 | 154 | 754 | 113,445,719 only-now; 20,855,669 kept |
| `site/data/lc/**/*.json` | 908 | 154 | 754 | 15,763,884 only-now (754 lines, one line each); 3,092,193 kept (154 lines) |

`site/data/candidates.csv` is 908 data rows and `site/data/progress.json` is the progress file. Both **paths** are in `exosearch-v1.0.0` (the contents have grown). They are not proposed for deletion.

Tag `exosearch-v1.0.0` also contains 31 plot filenames and 31 light-curve filenames that are already absent from every current branch. They remain on the tag. Nothing here deletes them again.

`failed_4` and `failed_5` do not occur as paths or as text in any commit (`git log --all -S`). There is nothing to keep and nothing to delete.

## S*_cal driver imports

Driver: `exoplanet_search/exosearch/v12_sstar_cal.py`, run as `python -m exosearch.v12_sstar_cal`. Sectors `(59, 60, 61, 62, 63, 64)`. Pinned search commit `bc05650d69f739fd5fea2f37629e52c5d524b3a8`. Frozen settings hash `7afb7fc178138d59c84fc509a5b47e09ebc8d855d5d1d8ef69e446ee939f8756`.

It also loads `exoplanet_search/runs/v1.1_s70_x86/finish.py` with `importlib` (`_load_pinned_finish`). Every name in the import block is used at least once (`batch.lc_path`, `interim.with_quality`, `nulltest.transform`, `vet_amend._set`, `worker_guard.require_clean`, and the others).

"Pinned bytes" means the blob hash equals `bc05650`. `v1.2` means `origin/exoplanet-search-v1.2`. `amend` means `origin/v1.2-amendment-h2`. `search` means `origin/exoplanet-search`. `master` has none of these files.

| Module | Where the file exists | Whose bytes |
| --- | --- | --- |
| `exosearch/__init__.py` | v1.2, vetter, amend, search, tag v1.0.0 | same on those branches |
| `exosearch/batch.py` | v1.2, vetter, amend, search, tag v1.0.0 | changed since v1.0.0; same on the four branches |
| `exosearch/config.py` | v1.2, vetter, amend, search, tag v1.0.0 | same |
| `exosearch/holdout.py` | v1.2, vetter, amend, search, tag v1.0.0 | same |
| `exosearch/interim.py` | v1.2, vetter, amend, search | same across those four; not in v1.0.0 |
| `exosearch/nulltest.py` | v1.2, vetter, amend, search, tag v1.0.0 | same |
| `exosearch/p_grid_halfwidth.py` | v1.2, vetter, amend, search | same; not in v1.0.0 |
| `exosearch/quality.py` | v1.2, vetter, amend, search | same; not in v1.0.0 |
| `exosearch/sstar_grid.py` | **vetter only** | pinned `bc05650` |
| `exosearch/v12_settings.py` | **vetter only** | pinned `bc05650` |
| `exosearch/vet.py` | all four, and tag v1.0.0 | **pinned bytes only on vetter** (`e99c1f1251`). v1.2, amend, and search still have the pre-unify file (`2c5d66fc48`) |
| `exosearch/vet_final.py` | **vetter only** | pinned `bc05650` |
| `exosearch/worker_guard.py` | v1.2, vetter, amend, search | same; not in v1.0.0 |
| `exosearch/t14_retention.py` | **vetter only** | pinned; not imported by the driver. Tests and `v12_settings` use it |
| `exosearch/v11/__init__.py` | v1.2, vetter, amend, search | same |
| `exosearch/v11/kernels.py` | v1.2, vetter, amend, search | same |
| `exosearch/v11/known_planets.py` | all four | **pinned bytes only on vetter** (`d95a489e8c`). The other three share `c86b183b0d` |
| `exosearch/v11/pipeline.py` | all four | **pinned bytes only on vetter** (`ffb724624f` vs `5239c50029`) |
| `exosearch/v11/registration_checks.py` | v1.2, vetter, amend, search | same |
| `exosearch/v11/scramble.py` | v1.2, vetter, amend, search | same |
| `exosearch/v11/settings.py` | v1.2, vetter, amend, search | same |
| `exosearch/v11/vet_amend.py` | v1.2, vetter, amend, search | same |
| `exosearch/v11/vet11.py` | all four | **pinned bytes only on vetter** (`11cfcd63a6` vs `e1c485c916`) |
| `exosearch/v11/run.py` | all four | **pinned bytes only on vetter** (`67769fe9a9` vs `f700eae7fe`) |
| `configs/v1.1.json` | all four | **pinned bytes only on vetter** (`04e9d8976d` vs `467aa91a72`) |
| `configs/v1.2.json` | **vetter only** | pinned |
| `runs/v1.1_s70_x86/finish.py` | all four | **pinned bytes only on vetter** (`09a1df3011` vs `5343f652cb`) |
| `exosearch/v12_sstar_cal.py` | v1.2 and amend, **not** on vetter and **not** at `bc05650` | this tree kept the v1.2 file. amend's file is the earlier driver |

Smoke import on this tree, with `PYTHONPATH=exoplanet_search` and the test virtualenv (`numpy`, `scipy`, `pandas`, `pyarrow`, `astropy`, `numba`, `wotan`):

```
import_ok bc05650d69f739fd5fea2f37629e52c5d524b3a8 (59, 60, 61, 62, 63, 64)
has_sstar_grid exosearch.sstar_grid
has_vet_final exosearch.vet_final
has_v12_settings 7afb7fc178138d59c84fc509a5b47e09ebc8d855d5d1d8ef69e446ee939f8756
```

No test in the tree imports the driver. A one-test smoke module should be added after the rebase, not before, so it targets the post-amendment driver.

## What was checked

- `ruff check exoplanet_search --select F401,F841,ARG,ERA`
- `vulture` at confidence 80 and 60 on `exosearch`, `scripts`, `analysis`, `runs`, and (at 80) `tests`
- Import-name use in the driver
- `git` blob comparison against `bc05650` and the four remote tips
- Search-path file list in `exoplanet_search/docs/preregistration/v1.1_early_nulls_2026-11-11/searchpath_hash.py` (the list behind hash `b9e5ebea8f84b05d60491e871f2c19e3facda75bdb146f3ab46eead42dafa4f4`)
- Byte identity with tag `exosearch-v1.0.0`
- Comment scan for lines that are commented-out statements (`# def`, `# return`, `# import`, and similar). None. The three `ERA001` hits are ordinary comments (`# review: B4`, a datetime example, and an SNR formula note).

Vulture's 60% hits are mostly false positives: `sstar_grid_enabled`, the `t14_retention` helpers, `evaluate_known_planets_v12`, `is_comb`, `false_comb_rate`, `expected_snr`, `calculate_grid_halfwidth`, and `EXCLUDED_BASE_SEED` / `EXCLUDED_RUNS` are all called from tests or from `v12_settings`. The 100% hits are the socket-blocker `kwargs`, which the signature has to keep.

No whole Python module is unused. Modules that nothing else imports are either `python -m` tools (documented in README, `reproduce.md`, or their own docstrings), preregistration scripts, or tests. They stay.

## Proposed removals

Not applied. Re-check every row that names `v12_sstar_cal.py`, `vet_final.py`, or `vet_amend.py` after the unpushed commits are fetched. If a name gains a caller, drop that row.

### 1. Generated per-star site files that are not in the v1.0.0 tag

Reason: committed output of the site publisher. Each file is one star's plot or one-line light-curve JSON. No checksum file and no registration document names these individual paths. The 154 plots and 154 light curves that are in tag `exosearch-v1.0.0` stay.

| | Files | Lines | Bytes |
| --- | ---: | ---: | ---: |
| `exoplanet_search/site/plots/**` not in the tag | 754 | n/a (binary) | 113,445,719 |
| `exoplanet_search/site/data/lc/**` not in the tag | 754 | 754 | 15,763,884 |
| **Total** | **1508** | **754** | **129,209,603** |

List them with:

```
git ls-tree -r --name-only HEAD exoplanet_search/site/plots exoplanet_search/site/data/lc \
  > /tmp/now.txt
git ls-tree -r --name-only exosearch-v1.0.0 exoplanet_search/site/plots exoplanet_search/site/data/lc \
  > /tmp/v100.txt
comm -23 <(sort /tmp/now.txt) <(sort /tmp/v100.txt)
```

Effect if applied: candidate pages for those TIC/sector pairs 404. `candidates.csv` still lists the stars. Search products, S*, and vetting are unchanged. Confirm before deleting; see "Unsure" below.

### 2. Unreachable helpers and unused names (behavior unchanged)

These are not in the search-path hash list and are not byte-identical to `exosearch-v1.0.0`.

| Item | Where | Lines | Reason |
| --- | --- | ---: | --- |
| `needs_lc_row` | `exosearch/vet_final.py` 305–318 | 14 | No caller. `test_needs_lc_bracket` calls `final_verdict` directly. |
| `bracket` | `exosearch/vet_final.py` 515–522 | 8 | No caller. Same test. |
| `ALL_SEARCH_MODES` | `exosearch/v12_sstar_cal.py` 105 | 1 | Defined, never read. Re-check after the amendment rebase. |
| `STAGES` | `exosearch/v11/known_planets.py` 77 | 1 | Never read. The stage names stay as literals in the stage function (lines 293–312). |
| `DEFAULT_PROD` | `scripts/publish_status_v12.py` 41 | 1 | Never read. |
| `SECONDARY` | `exosearch/v11/vet_amend.py` 33 | 1 | Never read. Rules (a) and (d) test the strings `"secondary"` and `"secondary_scan"` inline. |
| `any_group_evaluated` | `exosearch/v11/validate_pw2.py` 567 and 577 | 2 | Assigned, never read. `all_groups_descriptive` on the same branch is read and stays. |
| `import pandas as pd` | `exosearch/cut_sensitivity.py` 14 | 1 | F401. No other use. |
| `import pandas as pd` | `exosearch/nulldiag.py` 17 | 1 | F401. |
| `import numpy as np` | `exosearch/plateau_report.py` 14 | 1 | F401. |
| `import numpy as np` | `exosearch/rerun_rows.py` 22 | 1 | F401. |
| `import pandas as pd` | `exosearch/subsample.py` 18 | 1 | F401. |
| `import pandas as pd` | `exosearch/whitenoise.py` 18 | 1 | F401. |
| `import tempfile` | `tests/test_download.py` 2 | 1 | F401. |
| `import pytest` | `tests/test_v11_prewhiten.py` 3 | 1 | F401. |

Partial edits, one name dropped from a line that also imports names that are used:

| Edit | Line | Why the rest of the line stays |
| --- | --- | --- |
| Drop `config_hash` from `from .version import PIPELINE_VERSION, config_hash, stamp` | `exosearch/reliability.py` 22 | `PIPELINE_VERSION` and `stamp` are used. |
| Drop `grid` from `from . import grid, pipeline, run, search, settings` | `exosearch/v11/validate_pw2.py` 48 | The other four modules are used. `grid.` has no use in the file. |
| Drop `sys` from the combined import | `scripts/v11_launch.py` 18 | The other modules are used. |
| `rows, ctl = [], []` becomes `rows = []` | `exosearch/s105_diag.py` 82 | `ctl` is never appended or read. The loop still fills `rows`. |

Code subtotal if applied: **36 lines deleted**, **4 lines edited**, **0 modules deleted**. No threshold, seed, constant that is read, or numeric expression changes.

### 3. Not a removal: commented-out blocks and duplicated helpers

No commented-out statement blocks were found. Similar helpers were left in place:

- `finish.py` `_rebuild_flat_dumps` repeats the load/transform/prepare chain of `vet_amend.rebuild_flat` and also returns momentum-dump times. The comment says `rebuild_flat`'s return signature stays. Folding them would change that signature.
- SHA-256 helpers in the driver are local to the calibration record. They were not folded into `version.py`.

## Kept, and why they look unused

### Protected by the v1.0.0 tag (byte-identical)

Do not delete. The prior inventory called these teaching-lab leftovers. They are the tag.

| Path | Lines |
| --- | ---: |
| `01_fits_basics.ipynb` | 138 |
| `02_dask_distributed.ipynb` | 132 |
| `03_gpu_acceleration.ipynb` | 136 |
| `04_real_research_gaia.ipynb` | 160 |
| `05_mean_stacking.ipynb` | 107 |
| `06_multimodal_astronomy.ipynb` | 172 |
| `07_sql_and_catalogs.ipynb` | 76 |
| `check_storage.py` | 39 |
| `distributed_compute_template.py` | 32 |
| `hello_universe.py` | 41 |
| `star_finder_gpu.py` | 52 |
| `setup.sh` | 51 |
| `mkdocs.yml` | 23 |
| `docs/index.md` | 20 |

Also byte-identical to the tag, so unused imports inside them stay: `exosearch/catalogs.py` (`io`), `exosearch/learn_figs.py` (`lightcurve`), `exosearch/task.py` (`os`). `task.py` is also on the search-path hash list.

The 154 plots and 154 light curves whose paths are in the tag stay, even where the bytes have changed since the tag.

`candidates.csv` and `progress.json` stay because the paths are in the tag.

### Protected by the search-path hash

`searchpath_hash.py` hashes these, and the registered hash is `b9e5ebea8f84b05d60491e871f2c19e3facda75bdb146f3ab46eead42dafa4f4`. Editing them changes that hash even when runtime numbers would not. Leave them.

Dead-looking bits inside that set:

| Item | Where | Why it looks dead | Why it stays |
| --- | --- | --- | --- |
| `_gaia_mast` | `exosearch/deepvet.py` 156–162 | No caller. `gaia_neighbours` uses the ESA TAP path and says the MAST cones were DR2. | File is in the hash list. |
| `installed` | `exosearch/v11/kernels.py` 49–50 | No caller. | File is in the hash list. |
| `c3` from `a3, c3 = _fit(...)` | `exosearch/v11/prewhiten.py` 171 | `a3` is used in the threshold test; `c3` is not. The `_fit` call itself is live. | File is in the hash list. Replacing `c3` with `_` would still change the hash. |
| `import math` and unused `vet` | `exosearch/interim.py` 16 and 23 | F401. | File is in the hash list. |
| `import json` | `exosearch/v11/pipeline.py` 14 | F401. | File is in the hash list. |
| `import platform` | `exosearch/v11/run.py` 27 | F401. | File is in the hash list. |
| Unused parameters `t`, `f` | `exosearch/vet.py` 100 | ARG001. Callers pass them. | Signature plus hash list. |
| `durations` argument | `exosearch/v11/known_planets.py` 651 | ARG001. | Public signature. `known_planets.py` is not in the hash list, but dropping a parameter can change callers. Left alone. |

`vet.py`, `vet11.py`, `pipeline.py`, `run.py`, `known_planets.py`, `configs/v1.1.json`, and `finish.py` on this branch are the `bc05650` bytes. Do not edit them in the cleanup.

### Other kept-but-suspicious

| Item | Why it looks removable | Why it stays |
| --- | --- | --- |
| `exoplanet_search/docs/v1.1_remaining.md` | Title is "remaining work (handoff)", dated 2026-11-11 07:15 UTC, and tasks T0–T5 were later done or superseded. This is the stale handoff note. | Cited by `v1.1_prewhitening_amendment.md` and `v1.1_validation.md`. Not in the v1.0.0 tag, but it is a cited record. |
| `failed_4`, `failed_5` | Named in the keep-list. | Not in any branch, tag, or blob. Nothing to do. |
| `site/data/reliability.json` (56,593 bytes) and `site/data/review/gpd_threshold_choice_s0070.png` (60,634 bytes) | Not in the v1.0.0 tag. | `site/assets/methods.js` fetches both. Not the 908-star listing. |
| `site/data/status_v12.json` | Regenerated by the publisher. | The dashboard reads it. Latest content is the 2026-11-27 snapshot from `origin/exoplanet-search`. |
| `from . import interim_review as rv` at `reliability.py` 157 | `rv` is never used in that function. | It is the first line of a `try` that records `error` if the import raises. Deleting it changes that failure path. Line 322 uses `rv.tail_rate` and stays. |
| `find_pulsation_frequencies` | `validate_pw2.py` 127, no in-repo caller. | The source labels it a backward-compatibility alias. |
| `HASH_ONLY = {"exosearch/deepvet.py"}` | `dirty_stamps.py` 25, never read. `classify()` does not special-case it. | The comment says deepvet contributes only THRESH/REJECT to the hash. Deleting the name does not change today's output. Using it would change the disclosure classification. Left for a human. |
| `AMENDMENT` | `vet_amend.py` 24, never read. | It is the path of the vetter amendment. The rules do not consult it. |
| `docs/registration_v1.2_DRAFT.md` and `exoplanet_search/docs/registration_v1.2_DRAFT.md` | Identical 184-byte stubs. | The registration commits left stubs at both paths on purpose. |
| `docs/registration_v1.2.md` (root, 69,024 bytes) vs `exoplanet_search/docs/registration_v1.2.md` (89,788 bytes) | Two registration files. | Not duplicates. The root file is the amendment-h2 copy. The package file adds the dated sections. Both stay. |
| `sstar_refit.py` unused `sys` and `numpy` | F401, on lines that also import `gzip`, `json`, `time`, and `pandas`. | Preregistration script under `docs/preregistration/`. |
| `finish.py` names `cls` (line 78) and `j` (line 184) | F841. | `cls, rho_fac, r_fac = pipeline.draw_density(...)` advances the injection RNG. `j = (jobs or {}).get(...)` is unused in `rescore_rows_v12`, but `rebuild_for_row` on the next line reads `jobs`. Do not delete either statement. |
| Socket `_BlockedSocket.__new__` / `_blocked_create_connection` `args` and `kwargs` | Vulture 100%. | The signature is what makes the monkeypatch a valid `socket.socket` stand-in. |
| Odd-even mask, radius cap, and fine S* grid "default off" | Look like disabled flags. | Off is the v1.1 default. `configs/v1.2.json` turns them on and `require_v12_settings` refuses to run if they are still at the v1.1 default. Not dead. |
| CUDA / MLX branches in `bls_vector.py` and the worker scripts | Untested stubs. | They are the backend switch the README documents. Not deleted. |
| Every `python -m exosearch.*` module | Nothing imports many of them. | They are the v1.0 and v1.1 command-line tools. |

## Tests

Removals were not applied, so there is no "after" run. The baseline is this consolidated tree.

Smoke import: **passed** (see above).

Full suite: `cd exoplanet_search && python -m pytest tests -q --tb=no -p no:cacheprovider`

Baseline, this tree, 2026-11-27, Python 3.12, packages from `exoplanet_search/requirements.txt` except `lightkurve` and `astroquery` (not imported by the suite):

**270 passed, 5 skipped, 1 failed** in 399.57s.

The failure is `tests/test_v11_scramble.py::test_pipeline_scramble_mode_integration`. Re-run alone:

```
AssertionError: RuntimeError: Claret table unavailable: 404 Client Error: Not Found
for url: https://vizier.cds.unistra.fr/viz-bin/asu-tsv?-source=J/A%2BA/600/A30/tableab&-out.all&-out.max=999999
assert 'error' == 'ok'
```

`pipeline.run_job` calls `stellar.quad_ld`, which downloads the Claret 2017 TESS table. Both Vizier mirrors returned an error in this environment. The same `quad_ld` calls are in `pipeline.py` on `origin/exoplanet-search-v1.2` (lines 73 and 122 of that blob). The failure is the missing catalog, not a change introduced by the vetter merge. Two warnings, neither a failure: a pandas `.at` FutureWarning in `vet.py`, and a catalog-cache warning in `test_catalog_cache_falls_back`.

There is no after-cleanup run. No removal in this document was applied.

Historical context, not this run: the vetter-unify AI-use log records 274 passed, 2 skipped at `bc05650` before the S*_cal driver and the later v1.1 tests on `exoplanet-search-v1.2` existed. Do not treat that number as this tree's result.

## Unsure

1. The 1,508 post-v1.0.0 plot and light-curve files are the bulk of the old "~1,818 generated site files" estimate. The other ~310 of that estimate are the 154+154 tag members plus `candidates.csv` and `progress.json`, which stay. Deleting the 1,508 does not change search outputs, and it does remove pages the public site links by TIC. Left as a proposal until the owner confirms the site snapshot is disposable.
2. `needs_lc_row`, `bracket`, `ALL_SEARCH_MODES`, and `SECONDARY` have no caller on this tree. The unpushed driver may call them. Re-grep after the rebase before deleting.
3. `failed_4` and `failed_5` are not in the repository. If they are local run directories that were never committed, this audit cannot see them.
4. `v1.2-amendment-h2`'s S*_cal driver was not kept. It is an earlier parallel implementation (20 functions versus 62), not a strict subset of the later file. The later file is the one that matches the deviation and the common-period ruling. If the unpushed commits are based on the earlier driver, the rebase has to choose again.
5. `HASH_ONLY` looks like an unfinished classification rule in the v1.0 dirty-stamp disclosure tool. It was not "fixed", because using it would change that tool's output.
6. No second test run exists, because no removal was applied.
