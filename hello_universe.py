import os
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.utils.data import download_file
from astropy.visualization import astropy_mpl_style

plt.style.use(astropy_mpl_style)

# --- 1. Fetching Sample Space Data ---
print("🌌 Fetching a sample FITS image from Hubble/MAST...")
# Sample Hubble image (Horsehead Nebula)
image_url = 'https://archive.stsci.edu/pub/vosamp/horsehead.fits'
fits_file = download_file(image_url, cache=True)

# --- 2. Exploring the FITS Header ---
# FITS files consist of multiple Header Data Units (HDUs)
with fits.open(fits_file) as hdul:
    print("\n📦 HDU List:")
    hdul.info()
    
    # The first HDU is usually the image data
    header = hdul[0].header
    data = hdul[0].data
    
    print("\n🛰️ Satellite/Instrument:", header.get('INSTRUME', 'Unknown'))
    print("⏳ Exposure Time:", header.get('EXPTIME', 'Unknown'), "seconds")
    print("📏 Image Resolution:", data.shape)

# --- 3. Simple Image Statistics ---
print("\n📊 Image Stats:")
print(f"   Min: {data.min():.2f} | Max: {data.max():.2f}")
print(f"   Mean: {data.mean():.2f} | Std: {data.std():.2f}")

# --- 4. Basic Visualization ---
plt.figure(figsize=(10, 8))
plt.imshow(data, cmap='gray', origin='lower')
plt.colorbar()
plt.title('Hubble Space Telescope: Horsehead Nebula (FITS)')
plt.savefig('first_astro_image.png')
print("\n✅ Saved first visualization as 'first_astro_image.png'")
