import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.utils.data import download_file
from photutils.detection import DAOStarFinder
from astropy.stats import mad_std
import time

try:
    import cupy as cp
    GPU_AVAILABLE = True
    print("🚀 NVIDIA GPU detected via CuPy!")
except ImportError:
    GPU_AVAILABLE = False
    print("💻 GPU not found. Falling back to CPU.")

# 1. Fetch Hubble Image of a Star Field (M13 - Great Globular Cluster)
print("🌌 Fetching M13 Cluster star field...")
image_url = 'https://archive.stsci.edu/pub/vosamp/horsehead.fits' # Reusing local link for now
fits_file = download_file(image_url, cache=True)

with fits.open(fits_file) as hdul:
    data = hdul[0].data.astype(float)

# 2. Performance Comparison: Background Estimation
start_time = time.time()

if GPU_AVAILABLE:
    # Transfer data to GPU
    gpu_data = cp.array(data)
    bkg_sigma = cp.std(gpu_data) # Simplified for example
    print(f"⚡ GPU Background Sigma: {bkg_sigma:.2f}")
else:
    bkg_sigma = mad_std(data)
    print(f"🐢 CPU Background Sigma: {bkg_sigma:.2f}")

# 3. Star Detection (Using Photutils)
# This part is CPU-heavy but we can distribute it via Dask later
daofind = DAOStarFinder(fwhm=3.0, threshold=5.*bkg_sigma)
sources = daofind(data - np.median(data))

print(f"✨ Found {len(sources)} stars in the field!")
print(f"⏱️ Total processing time: {time.time() - start_time:.2f} seconds")

# 4. Plotting Results
plt.figure(figsize=(10, 8))
plt.imshow(data, cmap='viridis', origin='lower', vmax=np.percentile(data, 99))
plt.scatter(sources['xcentroid'], sources['ycentroid'], s=10, edgecolor='white', facecolor='none', alpha=0.5)
plt.title(f'Detected {len(sources)} Stars')
plt.savefig('star_detection_results.png')
print("✅ Results saved as 'star_detection_results.png'")
