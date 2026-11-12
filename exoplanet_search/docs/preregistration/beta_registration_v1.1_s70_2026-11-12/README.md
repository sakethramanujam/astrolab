> **Decision pending diagnosis (2026-11-11, room ruling):** the S* reference null will be decided after the wrap-point and leakage checks on these same null rows, via a dated post-hoc amendment signed off by Saketh. The 2026-11-11T15:54:53Z registration (7326715e) is also kept on record.
> Reader: commit 25c0bdc (branch v11/beta-reader). Diagnostic rules, written before any diagnostic result existed: diagnostics_rules.md in this directory.
> Committed as built; the "DRAFT for review" line below describes the state at build time. No numbers changed.

# Beta-treatment registration (v1.1, sector 70), 2026-11-12T05:16:58Z (UTC)

DRAFT for review: not committed, not pushed; configs/v1.1.json beta_option not set; production not launched.

Applied by exosearch.v11.beta_rows (reader + same-population guard) with exosearch.v11.beta_select's rule functions:
docs/preregistration/beta_treatment.md with amendment 1, exactly as written, on the v1.1 flipped null and scrambled r0
(nulls and injections only), sector 70 only, same population on both nulls (v1.1 amendment 2026-11-11, section 5).
The 2026-11-11T15:54:53Z registration (7326715e) is not replaced.

* TOST gate (90% CI of dS* inside [-0.5, +0.5], sector 70): {'a': False, 'b': False, 'c': False, 'd': False}
* Passing: []
* Decision: **b**. no option passes the TOST gate: b (most conservative) per rule 4; production is not launched
* Gate failed: True

## Inputs

* Null rows: code commit `fddfb2536894ce4469864f065faca9ce772c1e32`, tree `16e4f315cf39a48c5145fa7523ff0d4fc22fdfd2`, non-production (dry-run) rows; {'tree': '16e4f315cf39a48c5145fa7523ff0d4fc22fdfd2', 'commit': 'fddfb2536894ce4469864f065faca9ce772c1e32', 'config_hash': 'f9f10beec4f1', 'pipeline_version': '1.1.0', 'production': 'False'}
* Rows read 25309, duplicates 0, after dedupe 25309: {'invert': 12797, 'scramble_r0': 12512} (manifest tasks 25309); every task row ok.
* Same population: paired tuning quality-passing stars n = 7983 (flipped 8156, scrambled 7983; 173 excluded hosts removed from the flipped side); guard ok, excluded in paired set 0; alpha = 1/10244.
* Exclusion list `configs/scramble_exclude_s0070_2026-11-11.csv` sha256 `afa7c0122f3a1254177c87c85cbe7108acb9371c0b477dbf238615f03837cd97`.
* Statistic: v1.1 row S_opt_<option> (each option searched at its own peak). Completeness input: v1.0 interim sector-70 injections (beta_select.injections), scored at the v1.1 S70 S* per option.
* Search-path hash `164a67cbe9c93b1ea72e242df2b23bc1c6a98711251458e270dc229063a618a2` at the null commit and at the reader commit `25c0bdc2fd7809d1b8e63fbc7d2cded995b09a61`.

Full numbers: beta_selection.json in this directory.
