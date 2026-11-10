"""Load SPOC light curves, mask flagged cadences, detrend with wotan biweight, bin."""
from __future__ import annotations

import numpy as np
from astropy.io import fits
from wotan import flatten


def load_spoc(path):
    with fits.open(path, memmap=False) as h:
        hdr = h[0].header
        d = h[1].data
        t = np.asarray(d["TIME"], float)
        f = np.asarray(d["PDCSAP_FLUX"], float)
        e = np.asarray(d["PDCSAP_FLUX_ERR"], float)
        q = np.asarray(d["QUALITY"], int)
    meta = {k: hdr.get(h_) for k, h_ in [("tic", "TICID"), ("sector", "SECTOR"), ("tmag", "TESSMAG"),
                                          ("teff", "TEFF"), ("radius", "RADIUS"), ("logg", "LOGG"),
                                          ("ra", "RA_OBJ"), ("dec", "DEC_OBJ")]}
    m = np.isfinite(t) & np.isfinite(f) & np.isfinite(e) & (q == 0) & (e > 0)
    t, f, e = t[m], f[m], e[m]
    if len(f) == 0:
        return t, f, e, meta
    med = np.nanmedian(f)
    return t, f / med, e / med, meta


def detrend(t, f, e, window=0.75, clip_sigma=4.0):
    """wotan biweight; drop NaN edges and upward outliers (flares). Returns t, flat, e, trend, mask."""
    flat, trend = flatten(t, f, method="biweight", window_length=window,
                          edge_cutoff=0.0, break_tolerance=0.5, return_trend=True)
    m = np.isfinite(flat)
    resid = flat[m] - np.median(flat[m])
    mad = 1.4826 * np.median(np.abs(resid))
    m &= flat < 1.0 + clip_sigma * max(mad, 1e-6)
    return t[m], flat[m], e[m], trend, m


def bin_lc(t, f, e, minutes=10.0):
    """Weighted bin to `minutes`. Keeps transits >~1 h resolvable, cuts BLS cost ~5x."""
    dt = minutes / 1440.0
    idx = np.floor((t - t.min()) / dt).astype(int)
    w = 1.0 / e ** 2
    n = idx.max() + 1
    W = np.bincount(idx, w, n)
    good = W > 0
    tb = np.bincount(idx, w * t, n)[good] / W[good]
    fb = np.bincount(idx, w * f, n)[good] / W[good]
    eb = 1.0 / np.sqrt(W[good])
    return tb, fb, eb
