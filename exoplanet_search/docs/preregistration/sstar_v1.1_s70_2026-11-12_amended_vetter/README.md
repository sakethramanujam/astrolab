# S* refit under the amended vetter (v1.1, sector 70), 2026-11-12

Follows the vetter amendment (`docs/preregistration/v1.1_amendment_vetter_2026-11-12.md`, e2379bd, approved by Saketh
2026-11-12 12:43 PT, "Let's go."), section "S* and the gate". Computed 2026-11-13T01:43Z (UTC).

## Result

| vetter | S* flipped (b) | S* scrambled (b) | **S* (option b)** | 90% CI of dS* | TOST |
|---|---|---|---|---|---|
| old (`beta_rows.with_vetting`) | 9.699 | 11.821 | **11.821** | [-6.217, 0.733] | False |
| amended (rules a-d) | 9.645 | 13.706 | **13.706** | [-11.980, 0.129] | False |

* **S* for the v1.1 sector-70 search under the amended vetter: 13.71** (option b, the registered option; S* = max of the
  two nulls, rule 4, unchanged). It replaces 11.82.
* **Old-vetter S* on the x86 nulls = 11.821023**, equal to the registered arm64 value (11.821023, dd79f05) to all printed
  digits: the architecture has no effect on S*; the change is entirely the vetter rule change.
* All four options, both vetters: `sstar.json` (`old`, `new`; same layout as beta_selection.json per option).

## What moved
* Verdict changes (amended vs old), fail to pass only: flipped 109 (3,514 -> 3,623 vetted), scrambled 118 (3,896 -> 4,014);
  pass to fail: 0. On the paired tuning population: flipped vetted 2,632 -> 2,652, scrambled 2,995 -> 3,019.
* The scrambled-null S* rises 11.82 -> 13.71 because the extra vetted scrambled rows sit in the tail; the GPD threshold
  choice (`interim_review.gpd_diagnostics`) moves u_scr 6.00 -> 6.25 with them. The flipped S* barely moves (9.70 -> 9.65).
* 326 null rows had `secondary_scan` and were re-vetted from the rebuilt light curve; every rebuilt scan matched its
  row (`vet_amend.revet_one` fails otherwise).
* Descriptive: option-b injection completeness (v1.0 interim S70 injections, 15-25 bin, at each S*) 0.559 at 11.82,
  0.326 at 13.71.

## Inputs and method
* Null rows: the 25,309 beta-registration tasks (12,797 flipped `invert`, 12,512 scrambled `scramble_r0`) regenerated on
  x86 from the registered Nov 11 manifests (same stars, tasks, seeds), commit e2379bd (search-path hash
  `fb958a4e3287588b376c6f57c752c0e130ffb9d3d581840e476a165c7e23ae07`), non-production rows. Hosts: h1 22,649 rows, h3
  2,660 rows (2,800 written; 140 also run on h1, all bit-identical). h3 passed a 100-task bit-identical overlap against
  h1 (130 columns, 0 differing) before its rows were used. Reader: 0 missing, 0 extra, every task ok, manifest set equal.
  Details: `null_rows_provenance.json`. No arm64 rows used.
* Old vetter: `beta_rows.with_vetting` as committed. Amended vetter: `vet_amend.revet` (rules a-d) on the same rows,
  vetted = `passed_lc_checks_v2` and no common period (`vet.flag_common_periods`, unchanged) and checks logged.
  Check: the old `vetted` equals `passed_lc_checks` & no common period & checks logged on every row (0 mismatches).
* Rule: `beta_rows.beta_frames`, `beta_rows.sector_gate` (beta_select `_sstar`, bootstrap 2,000 seed 20261110, TOST),
  `beta_rows.decide`, `beta_rows.injection_completeness`, unchanged; the third registration's four leakage stars
  (TIC 243015498, 2052927695, 620213822, 630275176) dropped from both paired nulls after pairing, alpha N = 10,244
  unchanged; paired n = 7,979. Driver: `sstar_refit.py` (adapter only; no rule reimplemented).
* No language model was used for S*, beta or vetting.

## Next
Sector 71 known-host diagnostic at 13.71 (diagnostic only), then the sector 72 known-planet list is frozen at 13.71 and
pushed before any sector 72 output; sector 72 is the last gate attempt on this vetter.
