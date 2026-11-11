"""Period grid (fixed arithmetic) and the registered duration cap.

Grid: frequencies fmin + df * arange(n) (two separate float64 operations; np.arange's fill loop is contracted to
fused multiply-adds on arm64), same span and spacing as v1.0.0 (pmin 0.5 d, pmax = baseline / 2,
df = min(D) / baseline^2). grid_sha256() hashes the float64 bytes for the row provenance.

Duration cap (docs/v1.1_plan.md, "Registered duration-grid cap"; Hippke & Heller 2019, A&A 623, A39):
q_circ = T14 / P = (1/pi) * (3 pi / (G rho P^2))^(1/3) (b = 0, Rp << R*); trial duration D is allowed at period P
only if D <= margin * q_circ * P (margin 3). Stars without a TIC density: D <= 0.1 P.
"""
from __future__ import annotations

import hashlib

import numpy as np

G_CGS = 6.674e-8
RHO_SUN = 1.989e33 / (4 / 3 * np.pi * 6.957e10 ** 3)      # g/cm^3, as exosearch.stellar
DAY = 86400.0
GRID_SPEC = dict(pmin=0.5, min_transits=2, frequency_factor=1.0, form="fmin + df * arange(n), ascending period")


def period_grid(t, durations, pmin=0.5, min_transits=2, frequency_factor=1.0):
    """Ascending periods. Identical bytes on x86_64 and arm64."""
    baseline = float(np.max(t) - np.min(t))
    pmax = baseline / min_transits
    df = frequency_factor * float(np.min(durations)) / baseline ** 2
    fmin = 1.0 / pmax
    n = int(np.ceil((1.0 / pmin - fmin) / df))
    freqs = fmin + df * np.arange(n, dtype=float)
    freqs = freqs[freqs < 1.0 / pmin]
    return np.ascontiguousarray((1.0 / freqs)[::-1])


def grid_sha256(periods):
    return hashlib.sha256(np.ascontiguousarray(np.asarray(periods, dtype=np.float64)).tobytes()).hexdigest()[:16]


def q_circ(period_d, rho_sun):
    rho = np.asarray(rho_sun, float) * RHO_SUN
    P = np.asarray(period_d, float) * DAY
    return (1.0 / np.pi) * (3.0 * np.pi / (G_CGS * rho * P ** 2)) ** (1.0 / 3.0)


def d_max(period_d, rho_sun, margin=3.0, q_no_density=0.1):
    """Largest allowed trial duration (days) at each period."""
    P = np.asarray(period_d, float)
    if rho_sun is None or not np.isfinite(rho_sun) or rho_sun <= 0:
        return q_no_density * P
    return margin * q_circ(P, rho_sun) * P


def allowed(periods, D, rho_sun, margin=3.0, q_no_density=0.1):
    """Boolean mask over periods where trial duration D is inside the cap (and shorter than the period)."""
    P = np.asarray(periods, float)
    return (D <= d_max(P, rho_sun, margin, q_no_density)) & (D < P)
