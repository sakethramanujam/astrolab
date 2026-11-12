# Known-TOI miss audit, sector 70, clean rerun (2026-11-11)

Reruns `python -m exosearch.missed --sector 70` on a clean committed tree
(git_hash `68452b793e`, config_hash `80ffc7f29943`; no `-dirty` suffix). The earlier product carried
git_hash `50906a55b0-dirty` and config_hash `b600e28c92ec`.

Inputs: the results database and the catalogue snapshot (exofop_toi.csv) that the original run read, copied read-only.

Result: 86 missed TOI hosts out of 222 TOI hosts in the sector-70 results table (136 period-matched, 106 of them also
passed the v1.0 screens). Buckets: below_threshold 45, few_transits 25, outcompeted 8, masked 6, long_duration 2.
Every row (tic, bucket, S at the TOI ephemeris, threshold) is identical to the dirty-tree product.

Caveats for anyone quoting these numbers:

* The threshold used for `below_threshold` is the v1.0 in-pipeline screening floor S = 7.0
  (configs/detection_threshold.json), not the interim S* = 11.07. Re-bucketed at S* = 11.07, seven of the
  eight `outcompeted` hosts (S at ephemeris 7.0 to 9.3) move to `below_threshold`: below_threshold 52,
  few_transits 25, masked 6, long_duration 2, outcompeted 1.
* "Found" (`toi_period_match`) and the competing peak come from the sector-70 rows of the results table, which
  predate the v1.0.0 S statistic (no S, pipeline_version or git_hash stored). Only S at the TOI ephemeris is
  computed with the v1.0.0 code. The v1.1 known-planet recovery test (amendment 2026-11-11) replaces this audit.
