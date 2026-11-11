# Sector-70 delta Scuti list for the prewhitening pass rule (fixed 2026-11-11)

* Query run: 2026-11-11T19:52:32Z (UTC) against the ESA Gaia archive TAP service (gaiadr3), synchronous ADQL
  (`query.adql`; `{ids}` is filled with the Gaia DR2 source_id of each sector-70 target, 2,000 per request).
* Targets: all 12,797 sector-70 SPOC 2-min targets. The TIC v8 Gaia ID (a DR2 source_id) was available for 12,763.
  DR2 to DR3 via `gaiadr3.dr2_neighbourhood`, keeping the nearest DR3 match per TIC.
* Selection: `gaiadr3.vari_classifier_result.best_class_name` containing DSCT or SXPHE. In DR3 that is the single
  combined label `DSCT|GDOR|SXPHE`, so gamma Doradus stars are included by construction.
* Result: 25 stars (`dsct_s70.csv`), all from classifier `nTransits:5+`.
* The list is fixed. It is not refreshed or edited after this commit; a change needs a dated amendment.
