# v1.1 S70 x86 Production Report (**UNVALIDATED**, S* = 13.71)

**Status: UNVALIDATED (S72 gate failed 28/30, `ecd2bc0`).** None of the candidates above S* is called a detection.

* Pinned search commit: `c572bee4336a2398f5a47cd3031a38582905e30d` (`c572bee4`), search-path hash `fb958a4e3287588b376c6f57c752c0e130ffb9d3d581840e476a165c7e23ae07`.
* Manifest commit: `c8987eaed2eeeb0a7e5e031b453af1148201a11e` (`c8987eae`), host `h1` (`x86_64`), 10 workers, nice 19, thermal guard (pause 90 C / resume 80 C).
* Rows: 70,024 total (`known`: 102, `rest`: 69,922); 69,969 `ok`, 55 `unphysical` (all in `inject` tasks where `a/R* <= 1 + k`), 0 errors.
* Scoring script: `../finish.py` (non-inject rows vetted with `exosearch.v11.vet_amend.revet` + `vet.flag_common_periods`; injection rows scored with the injection-recovery logic).

## Summary Counts (`summary.json`)

* **Known-planet gate at S* = 13.71 (`known/known_planets_summary.json`)**:
  * Total targets: 306; eligible: 51 (0 errors).
  * Required (confirmed planets): **19 of 19 recovered** (0 missed, 0 aliases; `pass: true`).
  * Eligible TOIs (scored only): **28 of 32 recovered**.
* **Flipped-null check on known hosts (`flipped_known_hosts.csv`)**:
  * Eligible hosts checked: 51; vanished: 50; not vanished: 1 (`K2-29 b`, invert S = 14.08 at P = 3.2586 d).
* **Null counts above S* = 13.71**:
  * `invert` (n = 12,797): **2 vetted >= S*** (221 all >= S*)
  * `scramble_r0` (n = 12,512): **3 vetted >= S*** (385 all >= S*)
  * `scramble_r1` (n = 12,512): **0 vetted >= S*** (379 all >= S*)
  * `scramble_r2` (n = 12,512): **3 vetted >= S*** (370 all >= S*)
* **Candidates above S* = 13.71 (UNVALIDATED; none called a detection)**:
  * Total vetted >= S*: **84**
  * On known planet or TOI hosts: **67**
  * Other hosts: **17**
* **Injection recovery at S* = 13.71**:
  * Total injection tasks: 6,894 (`inject_run0`: 3,447, `inject_run1`: 3,447); `ok`: 6,839 (`inject_run0`: 3,425, `inject_run1`: 3,414); `unphysical`: 55 (`inject_run0`: 22, `inject_run1`: 33); errors: 0.
  * All `ok` (n = 6,839): matched 3,107; recovered (matched & S >= S*) **2,205** (0.3224); vetted & recovered **1,165** (0.1703).
  * Allowed range (`depth_inj < 0.05` & `rp_inj < 25`, n = 6,340): matched 2,643; recovered **1,793** (0.2828); vetted & recovered **1,122** (0.1770).
  * Allowed range with `snr_exp >= 13.71` (n = 2,238): matched 2,144; recovered **1,792** (0.8007); vetted & recovered **1,122** (0.5013).
  * Allowed range with `15 <= snr_exp < 25` (n = 428): matched 392; recovered **135** (0.3154); vetted & recovered **97** (0.2266).
* **Thermal pauses (`h1`)**:
  * Pause episodes: **150** (163 paused ticks, 150 resumed events, max CPU temp 98 C).
