# Registration: peak tie-break and cross-host reproducibility rules (v1.1)

Registered 2026-11-10, before any cross-host merge and before any v1.1 run. Adopted by the Scientist (gate option 2,
with option 1 as the stopgap; `docs/host_gate/README.md`).

## Tie-break

* Over all trial (P, D) of a star, let S_max be the largest scan statistic. Every grid point with
  **|S_i − S_max| ≤ 1e-9 · S_max** is a tie. Among ties the **lowest period** wins; equal periods go to the shortest
  duration. The number of grid points in that band is logged per star (`n_tie_band`).
* Stars with `n_tie_band` > 1 had their period decided by this rule rather than by the data; that fraction is
  reported with every run.

## Arithmetic

* Period grid = `fmin + df · arange(n)` in float64 (two separate operations), or a shipped grid file with its sha256.
* Detrending (biweight) is compiled without `fastmath` and without an on-disk numba cache (cache=False or an empty
  per-process cache), single-threaded.

## Hosts

* Until the gate (`exosearch/host_gate.py`: ~200 shared stars, |ΔS| ≤ 1e-6 and identical best (P, D) on every host)
  passes on all hosts with these rules, the registered statistic is computed only on x86_64 hosts (h1, h3). The
  arm64 host (h2) produces no production rows.
* GPU or float32 searches are separate statistics with their own nulls and injections, never mixed into S*.

## Amendment 1 — 2026-11-10 21:30 PST (RF Engineer). Predates every v1.1 run: no v1.1 data exists.

Only cross-host dry runs on the interim code (`exosearch/host_gate.py --fixed`, 200 S70 gate stars) and
diagnostics of stored v1.0.0 rows had been run when this was written. It replaces the "Tie-break" section above.

1. In-band grid points (|S_i − S_max| ≤ 1e-9·S_max, over all trial durations) are grouped by period index into
   **contiguous plateaus** (runs of adjacent grid indices) and **disjoint peaks** (separated runs; a single point is
   a plateau of width 1).
2. On a plateau, the reported **P is the plateau centre** (P_lo + P_hi)/2 and **sigma_P is half its width** in days
   (edges extended by half a grid step; a single point gets half a grid step). t0, D and the depth come from the
   lower-middle grid point of the plateau (`P_grid`); D is the shortest in-band duration there.
3. The **lowest-period rule chooses only between disjoint plateaus or peaks**, never within a plateau.
4. Logged per star: `n_plateaus`, `plateau_width` (grid steps, chosen plateau), `plateau_width_max`, `disjoint`
   (more than one in-band plateau/peak), `P_lowest_edge`, `sigma_P`, `n_tie_band`, and the plateau list.
5. sigma_P is carried into the candidate tables (site and exports).
6. The period grid is not changed by this amendment. Plateau widths are reported only; any future grid change follows
   the Ofir (2014, A&A 561, A138) oversampling criterion and is registered before v1.1 production.

### Amendment 1a — 2026-11-10 22:05 PST (terminology and flag; Scientist review). Still before any v1.1 run.

* The plateau half-width of rule 2 is named the **grid-degeneracy half-width**, column `P_grid_halfwidth` (days),
  everywhere (v1.1 rows, tables, exports, site). It is not a period error and is never labelled as one: it ignores the
  photometric uncertainty and P/n aliases. "sigma_P" in amendment 1 above means this quantity.
* New per-row flag **`period_ambiguous`** = number of transits ≤ 2: the true period may be P/n (n = 2, 3, …). Such
  periods are shown as approximate/ambiguous, not to full precision.
