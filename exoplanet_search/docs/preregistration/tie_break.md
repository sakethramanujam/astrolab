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
