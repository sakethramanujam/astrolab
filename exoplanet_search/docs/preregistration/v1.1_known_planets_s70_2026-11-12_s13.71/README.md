# v1.1 S70 known-planet expected-S freeze at S* = 13.71 (2026-11-12), amended gate and vetter

Frozen and pushed before any v1.1 sector 70 production output at S* = 13.71 exists. Supersedes the 11.82 freeze
(docs/preregistration/v1.1_known_planets_2026-11-12, f27e8cf, pre known-gate amendment) for the x86 production run.
Produced by `python -m exosearch.v11.known_planets freeze-eligibility --s-star 13.71 --sector 70` at commit a9b24c3
(search-path hash fb958a4e), searched stars from the sector 70 store; same construction as the S71 and S72 freezes
(known-gate amendment 7be0dab: expected S x0.66, >= 3 transits, confirmed planets required).

**Status: unvalidated.** The sector 72 gate, the last attempt on the amended vetter, failed 28 of 30 (ecd2bc0). Saketh's
go for the sector 70 search (2026-11-12 17:40 PT) was given before that result; this run proceeds on his instruction and
its outputs are labelled unvalidated.

- Targets: 306 (confirmed 162, active TOIs 144). Eligible: 51; error rows: 0.
- **Required (confirmed planets): 19.** Eligible TOIs scored only: 32.

## Required (confirmed)

| TIC | name | disposition | P (d) | n_tr | expected S raw | expected S x0.66 | T_span (d) |
|---|---|---|---|---|---|---|---|
| 398572544 | WASP-28 b | CONFIRMED | 3.408835 | 6 | 134.29 | 88.63 | 23.66 |
| 56399553 | K2-29 b | CONFIRMED | 3.258832 | 7 | 115.68 | 76.35 | 23.17 |
| 176220787 | WASP-136 b | CONFIRMED | 5.215362 | 3 | 114.69 | 75.70 | 23.87 |
| 36592530 | WASP-75 b | CONFIRMED | 2.484193 | 7 | 114.22 | 75.38 | 23.27 |
| 267215820 | TOI-5300 b | CONFIRMED | 2.262196 | 7 | 111.53 | 73.61 | 23.39 |
| 380907135 | WASP-118 b | CONFIRMED | 4.046043 | 5 | 110.69 | 73.06 | 23.77 |
| 58858008 | TOI-6842 A b | CONFIRMED | 5.391744 | 3 | 98.05 | 64.72 | 23.54 |
| 422655579 | WASP-71 b | CONFIRMED | 2.903670 | 7 | 95.31 | 62.90 | 23.67 |
| 250111245 | TOI-5293 A b | CONFIRMED | 2.930289 | 7 | 65.53 | 43.25 | 23.53 |
| 270380593 | WASP-156 b | CONFIRMED | 3.836163 | 6 | 63.93 | 42.19 | 23.73 |
| 336691874 | TOI-5322 b | CONFIRMED | 5.423321 | 3 | 63.25 | 41.75 | 23.54 |
| 262662119 | WASP-151 b | CONFIRMED | 4.533478 | 4 | 53.98 | 35.63 | 23.36 |
| 293457754 | WASP-181 b | CONFIRMED | 4.519506 | 4 | 53.68 | 35.43 | 23.69 |
| 37718056 | WASP-158 b | CONFIRMED | 3.656333 | 5 | 43.39 | 28.64 | 23.60 |
| 28961316 | TOI-3682 b | CONFIRMED | 3.346241 | 5 | 41.92 | 27.67 | 23.32 |
| 26017005 | K2-30 b | CONFIRMED | 4.098503 | 5 | 40.84 | 26.96 | 23.33 |
| 26054627 | TOI-5349 b | CONFIRMED | 3.317921 | 6 | 33.91 | 22.38 | 23.33 |
| 58825110 | TOI-5301 b | CONFIRMED | 5.858858 | 3 | 31.02 | 20.47 | 23.48 |
| 33397739 | TOI-5292 A b | CONFIRMED | 2.021910 | 8 | 23.40 | 15.44 | 23.16 |

## Eligible TOIs (scored, not pass/fail)

| TIC | name | disposition | P (d) | n_tr | expected S raw | expected S x0.66 | T_span (d) |
|---|---|---|---|---|---|---|---|
| 422457210 | TOI 6163.01 | PC | 9.975470 | 3 | 83.74 | 55.27 | 23.47 |
| 369474541 | TOI 6843.01 | PC | 5.497787 | 4 | 70.19 | 46.32 | 23.65 |
| 4619242 | TOI 5295.01 | APC | 4.302357 | 4 | 57.27 | 37.80 | 23.32 |
| 349582831 | TOI 5100.01 | PC | 3.453214 | 6 | 51.00 | 33.66 | 23.11 |
| 408643099 | TOI 5324.01 | PC | 5.026905 | 5 | 47.83 | 31.57 | 23.69 |
| 94609570 | TOI 3655.01 | APC | 3.332565 | 6 | 47.74 | 31.51 | 23.22 |
| 327952677 | TOI 295.01 | KP | 3.396944 | 5 | 47.20 | 31.15 | 23.86 |
| 427740976 | TOI 4608.01 | PC | 2.349490 | 8 | 44.80 | 29.57 | 23.12 |
| 435336785 | TOI 4543.01 | PC | 5.774032 | 4 | 44.31 | 29.24 | 23.47 |
| 9725627 | TOI 239.01 | KP | 4.156783 | 4 | 43.99 | 29.03 | 23.73 |
| 303502327 | TOI 5327.01 | PC | 4.667118 | 4 | 43.48 | 28.69 | 23.45 |
| 4464022 | TOI 4542.01 | PC | 5.916768 | 4 | 42.32 | 27.93 | 23.17 |
| 337025682 | TOI 5302.01 | PC | 4.241769 | 4 | 40.50 | 26.73 | 23.69 |
| 34781378 | TOI 5354.01 | PC | 4.794929 | 4 | 39.75 | 26.24 | 23.38 |
| 402319411 | TOI 317.01 | PC | 3.326298 | 5 | 38.71 | 25.55 | 23.35 |
| 302202317 | TOI 5314.01 | PC | 2.535102 | 9 | 34.21 | 22.58 | 23.55 |
| 435259601 | TOI 4660.01 | PC | 1.509081 | 14 | 33.32 | 21.99 | 23.43 |
| 272146465 | TOI 5323.01 | APC | 5.804938 | 4 | 30.13 | 19.89 | 23.69 |
| 18373995 | TOI 4617.01 | PC | 3.487875 | 6 | 28.69 | 18.93 | 23.24 |
| 284173938 | TOI 3751.01 | PC | 2.915704 | 7 | 27.46 | 18.12 | 23.23 |
| 126982221 | TOI 5315.01 | PC | 2.470407 | 9 | 26.05 | 17.19 | 23.65 |
| 426448042 | TOI 5363.01 | PC | 6.750476 | 4 | 25.97 | 17.14 | 23.41 |
| 113804153 | TOI 5369.01 | PC | 4.246035 | 5 | 25.91 | 17.10 | 23.36 |
| 56457123 | TOI 5368.01 | PC | 3.612762 | 6 | 25.27 | 16.68 | 23.12 |
| 47094691 | TOI 3638.01 | PC | 4.927689 | 4 | 24.98 | 16.49 | 23.70 |
| 467684543 | TOI 5356.01 | PC | 3.194869 | 6 | 23.53 | 15.53 | 23.24 |
| 29017001 | TOI 4598.01 | PC | 8.200375 | 3 | 23.29 | 15.37 | 23.35 |
| 56655478 | TOI 5066.01 | PC | 1.345408 | 14 | 22.72 | 15.00 | 23.05 |
| 347469901 | TOI 5313.01 | PC | 3.257104 | 6 | 22.68 | 14.97 | 23.47 |
| 251039147 | TOI 4537.01 | PC | 6.662995 | 4 | 21.99 | 14.51 | 23.31 |
| 283605976 | TOI 4679.01 | PC | 6.817091 | 4 | 20.97 | 13.84 | 23.44 |
| 301035553 | TOI 5286.01 | PC | 3.298882 | 4 | 20.82 | 13.74 | 23.87 |
