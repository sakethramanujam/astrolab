# Registration: the "marginal" label

Interim analysis layer over frozen v1.0.0 (config_hash 80ffc7f29943). Registered 2026-11-10.

## Rule

A signal that passes every vetting check with S >= S* of its own sector is labelled **marginal** if the expected
number of post-vetting false alarms per sector at its S is **>= 0.1**:

    E_FA(S) = N_stars * (k / n) * GPD.sf(S - u; xi, sigma)

where the GPD (shape xi, scale sigma) is fitted to the vetted flipped-null max S above u (the same fit, u rule and
quality cut that define S*), k of n tuning null stars exceed u, and N_stars is the number of quality-passing stars
searched in that sector. Otherwise it is labelled "survives vetting". Neither label is a detection.

* Only sectors with their own flipped null have S* or E_FA; signals in other sectors are "not calibrated".
* The threshold 0.1 is fixed. It is not tuned on any candidate.
* E_FA is re-evaluated only when S* is (see `s_star_freeze.md`); once S* is frozen, E_FA uses the frozen fit.

## Disclosure

The rule was first applied in commit b5119772 (reliability survivors), before this file was committed. At that point
it labelled TIC 318707987 (S 14.8, E_FA 0.31) and TIC 242942524 (S 16.6, E_FA 0.15) marginal. The 0.1 value was chosen
then, after both signals' S values were known. From this commit on the rule is fixed and applied unchanged.

Code: `exosearch/reliability.py` (`MARGINAL_EFA = 0.1`), `exosearch/interim_review.py` (`tail_rate`).
