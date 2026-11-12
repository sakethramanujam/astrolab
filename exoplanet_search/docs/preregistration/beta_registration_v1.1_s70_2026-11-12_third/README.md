# Beta-treatment registration (v1.1, sector 70), third, 2026-11-12T17:30:07Z (UTC)

Third registration under the post-hoc join-gap trim amendment (docs/preregistration/v1.1_amendment_join_gap_trim_2026-11-12.md,
commit 9bf02ac, approved by Saketh 2026-11-12; join gap: 8025edc). Pass test passed only vacuously (see diagnostics). The 2026-11-11T15:54:53Z (7326715e) and 2026-11-12 (fa3bbf8) gate failures
stay on record. configs/v1.1.json beta_option not set; production not launched.

Applied by exosearch.v11.beta_rows (reader 25c0bdc, branch v11/beta-reader; rule files identical to the null code commit)
with exosearch.v11.beta_select's rule functions: docs/preregistration/beta_treatment.md with amendment 1, as written, on
the rerun v1.1 flipped null and scrambled r0, sector 70, same population on both nulls. Per the amendment, the four
leakage stars [243015498, 2052927695, 620213822, 630275176] are removed from BOTH paired nulls after the reader's pairing (alpha's N unchanged).

* TOST gate (90% CI of dS* inside [-0.5, +0.5], sector 70): {'a': False, 'b': False, 'c': False, 'd': False}
* Passing: []
* Decision: **b**. no option passes the TOST gate: b (most conservative) per rule 4; production is not launched
* Gate failed: True

| option | S* flipped | S* scrambled | S* | 90% CI of dS* | TOST |
|---|---|---|---|---|---|
| a | 9.797 | 13.481 | 13.481 | [-11.049, 0.376] | False |
| b | 9.699 | 11.821 | 11.821 | [-6.217, 0.733] | False |
| c | 9.816 | 13.622 | 13.622 | [-11.486, 0.421] | False |
| d | 9.906 | 11.747 | 11.747 | [-5.463, 0.975] | False |

## Inputs

* Null rows: code commit `0c4f81e58ff690da685162f151593e4026bbf406` (fddfb25 + join gap + join-gap trim), tree `98f4dc6fa99ce0af68042e82432c2e0dd2f558e6`, non-production (dry-run) rows; {'tree': '98f4dc6fa99ce0af68042e82432c2e0dd2f558e6', 'commit': '0c4f81e58ff690da685162f151593e4026bbf406', 'config_hash': 'f9f10beec4f1', 'pipeline_version': '1.1.0', 'production': 'False'}
* Rows read 25309, duplicates 0, after dedupe 25309: {'invert': 12797, 'scramble_r0': 12512} (manifest tasks 25309); every task row ok.
* Same population: paired n = 7979 after the leakage drop (7983 before; flipped 8156, scrambled 7983); guard ok; alpha = 1/10244.
* Exclusion list `configs/scramble_exclude_s0070_2026-11-11.csv` sha256 `afa7c0122f3a1254177c87c85cbe7108acb9371c0b477dbf238615f03837cd97`.
* Statistic: v1.1 row S_opt_<option> (each option searched at its own peak). Completeness input: v1.0 interim sector-70 injections (beta_select.injections), scored at the v1.1 S70 S* per option.
* Search-path hash `fb958a4e3287588b376c6f57c752c0e130ffb9d3d581840e476a165c7e23ae07` at the null commit.

Full numbers: beta_selection.json. Join and leakage diagnostics and the amendment's pass test: diagnostics.md / diagnostics.json.

## Decision of record (2026-11-12, signed off by Saketh in the Astrolab room at 09:41 PT)

* S* for the v1.1 sector-70 search is **11.82** (option b). The reason of record is the rule registered in advance: beta failed the TOST gate for the third time, so S* = max(S* scrambled, S* flipped) under the most conservative option, rule 4.
* Join-gap patching of the scrambled null stops here. No further null amendments for sector 70.
* The amendment's pass test passed only vacuously. The 0.3 d edge trim removes every cadence within 0.3 d of an inserted join, and the longest trial transit duration is 0.25 d, so no trial box can reach a join and the count is 0 by construction. The trim was proposed by RF Engineer without checking it against the longest trial duration; that oversight is recorded here. The supplementary check on the vetted top 1% (inserted-gap edges 3.89% vs 2.26% chance, original edges 3.57% vs 2.22%) is the informative number.
* Dated post-hoc interpretation (2026-11-12, not a reason of record): RF Engineer's argument that the excess scrambled-null tail comes from detrender edge behaviour at joins, and that the biweight-window shoulders seen in the flipped null (6 excluded known hosts recovered above 9.70: K2-29 b, TOI-5322 b, TOI-6843.01, WASP-71 b, TOI-4608.01, TOI-5328.01) bias S* upward, i.e. conservatively.
* The merged `exoplanet-search` search code hashes to `fb958a4e…` (same as the null commit 0c4f81e), so these nulls apply to the merged code.
* Next, in order: known-planet recovery gate at S* = 11.82, then the full sector-70 v1.1 search.
