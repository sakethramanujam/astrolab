# Catalogue cross-check outputs (v1.2, h2 run)

## Overview

This directory holds outputs from a read-only diagnostic comparing the search's S59-S64 real-mode peaks with published ephemerides from the NASA `pscomppars` snapshot (2026-11-13) and the TOI table.
It is not a gate and does not change any registered threshold.

Specifications:

- Original: [registration_v1.2.md](../../docs/registration_v1.2.md), section "Catalogue cross-check spec 2026-11-19" (tag `v1.2-xcheck-spec-2026-11-19`), also [xcheck_catalogue_v1.2_2026-11-19.md](../../docs/xcheck_catalogue_v1.2_2026-11-19.md)
- Amendment: [registration_v1.2.md](../../docs/registration_v1.2.md), section "Catalogue cross-check amendment 2026-11-20" (tag `v1.2-xcheck-amend-2026-11-20`, commit `e4f0f45`, written after the original result was seen; both results stay on record), also [xcheck_catalogue_v1.2_amendment_2026-11-20.md](../../docs/xcheck_catalogue_v1.2_amendment_2026-11-20.md)

## Files committed here

- `xcheck_v2.py`: amendment scoring script, committed unchanged so its SHA-256 below still matches; its docstring names the amendment file. It imports the original scoring module `xcheck.py`, which is not committed here; its SHA-256 is listed below.
- `xcheck_summary_v2.json`: aggregate counts only.

Per-planet CSVs are not committed because they contain star identifiers. They stay local and are pinned by SHA-256 below.

## Results

- Original spec: 103 of 376 scored planets pass all criteria.
- Amendment: 192 of 297 scored planets pass all criteria.
- Scope (420 planets in scope): 42 outside grid, 80 outside period range, 1 catalogue incomplete.
- Amendment per-criterion planet-level passes: period 219, t0 226, duration 243, depth 257.
- Siblings: 16 planets where a sibling on the same star was returned instead; 12 siblings credited.

## Trace of the 9 period misses (S >= 10, original spec)

Each of the 9 period misses with `S >= 10` under the original spec was attributed as follows:

| Attribution | Count | How it was attributed |
| --- | --- | --- |
| Sibling returned instead | 5 | The reported peak period and epoch matched another catalogued planet on the same star within the cross-check tolerances. |
| Eclipsing-binary signs | 2 | Folding at the reported peak showed odd/even depth difference or a secondary eclipse, consistent with an eclipsing binary rather than the catalogued planet. |
| Momentum-dump artefact | 1 | The in-transit points of the reported peak coincided with cadences flagged as reaction-wheel momentum dumps. |
| Detrending weakness | 1 | Fast stellar variability left residual structure after detrending that outranked the catalogued transit. |

No code bug was found.

Scripts used for this trace: `trace9.py`, `trace9_depth.py`, `trace9_fold.py`, and `diag_period.py` (not committed; hashes below).

## SHA-256 of all cross-check outputs

| File | SHA-256 | Committed here |
| --- | --- | --- |
| `xcheck.py` | `cf713df9668f0192f2ad6cce93956602100cc8d9ce3f94a049874ef4c0b8dfe7` | no |
| `xcheck_v2.py` | `304bf02f9732ba8a105c9900941526de4f9af2811d8af1229279620c3c3862c8` | yes |
| `xcheck_summary.json` | `2625c24acace5801d4479169d9dc094267dffc03b49b5c77cee106a42ee7eeeb` | no |
| `xcheck_summary_v2.json` | `3e43eb9851a39f406d3fba2fdb764c37ce50d274ee7ee672d35d8488c2a34665` | yes |
| `xcheck_planet.csv` | `cb0041d656d4488deea2f256d347f96b87bdfdb5ebbc58d6517e888a88c48514` | no |
| `xcheck_planet_sector.csv` | `02409bdf519f280d1121039b0bcf0fa45bdf3f000cd983366dcea9963a3ba8d0` | no |
| `xcheck_mismatches.csv` | `56860731a23c24a74cb655f7cd25d8c30646948064f685d87d5d227ade3083a8` | no |
| `xcheck_mismatches_annotated.csv` | `46c47cbbdf3dcb4717d84de140135bf1c9a9885c2d1a43ad31a868a67492512e` | no |
| `xcheck_planet_v2.csv` | `8b3a016688b7fb23ed6c337ad2179fed3fcd0c5a12f034e45895231f1ec4fac2` | no |
| `xcheck_planet_sector_v2.csv` | `4ea421fb156de45290e17ba07410260ac88e0f4a200c7f0839efd9f8e5d582c9` | no |
| `xcheck_mismatches_v2.csv` | `f2f5a88d72f525d366a776cfd86e2c30a8640ddb2cfc0a28c06e81a3dec8d780` | no |
| `diag_period.py` | `82a90fead4daa8cfb72159b7351e93f679794df0603d7ca0d105bdf36c152aad` | no |
| `diag_period_rows.csv` | `08d8d44911d56030703b691507cc030ff8576d5904f5ba7b78648328083f1f85` | no |
| `trace9.py` | `8bd9e3324feb94dda83b682280866be70cf8bd0509747ca19f241cc0b307ff74` | no |
| `trace9_depth.py` | `8da3b2aca2a5d0f6801ec37ccb7676baca6488475601b618d42f2a9f9aedf395` | no |
| `trace9_fold.py` | `f022fc6177a15d74a9508b5596022a8d8a4df426361b67ad2700c42924069a55` | no |
| `trace9_rows.csv` | `1d43b1f511ffd346456928cfd4d15adbd733e657d60d281c0c9e7ba86d4d2f34` | no |
