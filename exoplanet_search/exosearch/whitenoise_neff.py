"""Effective number of independent trials from the white-noise calibration (analysis layer; v1.0.0 untouched).

The v1.0.0 final statistic S = depth sqrt(N_tr) / CDPP_D drops the out-of-transit term of the depth error.
With it restored, sigma_delta^2 = CDPP_D^2 / N_tr * 1 / (1 - q), q = in-transit fraction, so
  S_oot = S * sqrt(1 - q),   q = D / P (uniform sampling of phase; exact in-transit cadence fraction to <1%).
The peak *selection* is by the scan statistic (astropy BLS 'snr', which already carries the out-of-transit term),
so restoring the term changes the value at the chosen peak, not which peak is chosen: re-running the scan would
reproduce the stored peaks exactly, and S_oot is computed from them.

Gumbel fit to the per-star maxima (location mu, scale beta). For the max of N independent unit Gaussians,
mu = Phi^-1(1 - 1/N), beta ~ 1 / (N phi(mu)); N_eff is read off mu (and, as a check, off beta).
Expected max = mu + 0.5772 beta. A per-trial standard deviation sd != 1 is fitted jointly by matching both.
  python -m exosearch.whitenoise_neff --sector 70
Output: products/v<ver>/calibration/whitenoise_neff_s<sector>.json and .png
"""
from __future__ import annotations

import argparse
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import optimize, stats

from . import config

EULER = 0.5772156649


def neff_from_gumbel(mu, beta):
    n_mu = 1.0 / stats.norm.sf(mu)
    n_beta = 1.0 / (beta * stats.norm.pdf(mu)) if beta > 0 else np.nan
    return float(n_mu), float(n_beta)


def joint_sd_n(mu, beta):
    """Solve for (sd, N) such that the max of N iid N(0, sd^2) has Gumbel location mu and scale beta."""
    def eq(v):
        sd, ln = v
        N = np.exp(ln)
        m = sd * stats.norm.isf(1 / N)
        b = sd / (N * stats.norm.pdf(m / sd))
        return [m - mu, b - beta]
    (sd, ln), _, ier, _ = optimize.fsolve(eq, [1.0, np.log(1e5)], full_output=True)
    if ier != 1:
        return None, None                      # no (sd, N) reproduces both Gumbel parameters (not a Gaussian max)
    return float(sd), float(np.exp(ln))


def fit(x):
    mu, beta = stats.gumbel_r.fit(x)
    n_mu, n_beta = neff_from_gumbel(mu, beta)
    sd, n_j = joint_sd_n(mu, beta)
    ks = stats.kstest(x, "gumbel_r", args=(mu, beta)).pvalue
    return dict(n=int(len(x)), median=float(np.median(x)), mean=float(np.mean(x)), gumbel_mu=float(mu), gumbel_beta=float(beta),
                expected_max=float(mu + EULER * beta), N_eff_unit_sd_from_mu=n_mu, N_eff_unit_sd_from_beta=n_beta,
                joint_trial_sd=sd, joint_N_eff=n_j, ks_p=float(ks))


def s_star_white(n_eff, n_stars, sd=1.0):
    """Gaussian white-noise threshold at alpha = 1 / N_stars over N_eff trials per star."""
    return float(sd * stats.norm.isf(1.0 / (n_eff * n_stars)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sector", type=int, default=70)
    ap.add_argument("--n-stars", type=int, default=None, help="population size for the white-noise S* (default: interim n_stars)")
    a = ap.parse_args()
    cal = config.DATA_ROOT / "products" / "v1.0.0" / "calibration"
    d = pd.read_parquet(cal / f"whitenoise_s{a.sector:04d}.parquet")
    d = d[d.status == "ok"].copy()
    d["q"] = d.duration / d.period
    d["S_oot"] = d.S * np.sqrt(1 - d.q)
    if a.n_stars is None:
        from . import interim
        a.n_stars = int(json.loads(interim.THR_FILE.read_text())["sectors"][str(a.sector)]["n_stars"])
    res = dict(sector=a.sector, n_stars_for_threshold=a.n_stars, q_median=float(d.q.median()), q_max=float(d.q.max()),
               naive_trials_per_star=float(d.n_periods.median() * 8), fits={})
    for col in ("S_scan", "S", "S_oot"):
        f = fit(d[col].values)
        f["s_star_white_unit_sd"] = s_star_white(f["N_eff_unit_sd_from_mu"], a.n_stars)
        f["s_star_white_joint"] = s_star_white(f["joint_N_eff"], a.n_stars, f["joint_trial_sd"]) if f["joint_N_eff"] else None
        res["fits"][col] = f
    fig, ax = plt.subplots(figsize=(6, 3.6), dpi=100)
    bins = np.linspace(3.5, 10.5, 36)
    for col, c in (("S_scan", "#4c78a8"), ("S", "#999999"), ("S_oot", "#e45756")):
        ax.hist(d[col], bins, histtype="step", density=True, color=c, label=col)
        f = res["fits"][col]
        xx = np.linspace(3.5, 10.5, 300)
        ax.plot(xx, stats.gumbel_r.pdf(xx, f["gumbel_mu"], f["gumbel_beta"]), color=c, lw=0.8, ls=":")
    ax.set_xlabel("per-star max S, white noise, real S70 cadence and gaps"); ax.set_ylabel("density")
    ax.set_title(f"Gumbel fits: N_eff (unit sd) S_oot = {res['fits']['S_oot']['N_eff_unit_sd_from_mu']:.2g}", fontsize=9)
    ax.legend(fontsize=8); fig.tight_layout()
    fig.savefig(cal / f"whitenoise_neff_s{a.sector:04d}.png"); plt.close(fig)
    (cal / f"whitenoise_neff_s{a.sector:04d}.json").write_text(json.dumps(res, indent=1, allow_nan=False))
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
