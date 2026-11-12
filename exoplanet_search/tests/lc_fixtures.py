"""Small synthetic SPOC-like FITS light curves for end-to-end validation tests."""
from pathlib import Path

import numpy as np
from astropy.io import fits


def s70_like_times(cadence_min=2.0):
    t = np.arange(3000.0, 3027.0, cadence_min / 1440.0)
    return t[(t < 3013.0) | (t > 3014.0)]


def transit_flux(t, period, t0, dur, depth, noise=1e-3, seed=1):
    rng = np.random.default_rng(seed)
    f = 1.0 + rng.normal(0, noise, t.size)
    ph = ((t - t0 + 0.5 * period) % period) - 0.5 * period
    f[np.abs(ph) < dur / 2] -= depth
    return f


def write_spoc_fits(path, t, f, e=None, crowdsap=0.98, tic=1, sector=70):
    """Minimal SPOC 2-min light-curve file readable by exosearch.lightcurve.load_spoc_raw."""
    e = np.full_like(t, 1e-3) if e is None else e
    flux = 1000.0 * np.asarray(f, float)
    cols = [fits.Column(name="TIME", format="D", array=np.asarray(t, float)),
            fits.Column(name="PDCSAP_FLUX", format="E", array=flux),
            fits.Column(name="PDCSAP_FLUX_ERR", format="E", array=1000.0 * np.asarray(e, float)),
            fits.Column(name="QUALITY", format="J", array=np.zeros(t.size, np.int32))]
    h0 = fits.PrimaryHDU()
    h0.header["TICID"] = tic
    h0.header["SECTOR"] = sector
    h0.header["RADIUS"] = 1.0
    h1 = fits.BinTableHDU.from_columns(cols)
    if crowdsap is not None:
        h1.header["CROWDSAP"] = crowdsap
    fits.HDUList([h0, h1]).writeto(path, overwrite=True)
    return Path(path)
