# Join-gap trim rerun: pass test and diagnostics (code 0c4f81e, hash fb958a4e; amendment 9bf02ac)
* Registered pass test (top 1% vetted scrambled, n=80): k=0 of 22,213 in-transit cadences within one duration of a former
  join; p0=0; bound 0 -> PASS, but VACUOUS: the 0.3 d trim removes every cadence within 0.3 d of a former join and the
  longest trial duration is 0.25 d, so the statistic is 0 by construction. Rule-set top 1% (all rows): likewise 0/5,999.
* Supplementary (not registered), same pooled test against the new kept edges at inserted gaps: vetted 3.89% vs 2.26%
  (863 > bound 553); real (original) segment edges: 3.57% vs 2.22% (793 > 545). All rows: new edges 6.80% vs 2.74%,
  real edges 3.75% vs 2.70%. Inserted gaps now behave like real gaps; both sit above chance under this strict pooled test.
* Leakage: 0/80 vetted; 1/80 all rows (not above the shuffle q99).
* Diagnostic option-b S*: scrambled 11.82 (join gap only: 10.82; before: 14.52), flipped 9.70.
