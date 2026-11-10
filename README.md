# AstroCompute Lab: Networked Astronomical Data Processing

This lab focuses on learning how to process astronomical data (FITS images, lightcurves, spectroscopy) and distribute the workload across your local network devices.

## 🚀 Learning Path
1. **FITS Data Handling**: Learn to read and manipulate Flexible Image Transport System (FITS) files.
2. **Image Processing**: Basic reduction (bias, dark, flat-field correction).
3. **Photometry**: Measuring the brightness of stars.
4. **Distributed Computing**: Using Dask to parallelize star detection or pixel math across `<WORKER_HOSTNAME>`, `<MASTER_HOSTNAME>`, and other local nodes.

## 📡 Network Compute Resources
- **Master**: `<MASTER_HOSTNAME>` (<MASTER_IP>)
- **Workers**: `<WORKER_HOSTNAME>` (<WORKER_IP>), `.32` (potential node)
- **Shared Storage**: `//<MASTER_IP>/shared` (Mounted at `/media/nfs01` on Master)

### Mounting Shared Storage on Workers
To ensure all nodes see the same data at the same path:
```bash
sudo mkdir -p /mnt/astro_data
sudo mount -t cifs -o username=spiceman //<MASTER_IP>/shared /mnt/astro_data
```

## 🛠️ Getting Started
```bash
# Clone and enter the project
git clone [repo_url] astro-compute-lab
cd astro-compute-lab

# Run the unified setup script (Installs dependencies and prepares the env)
bash setup.sh

# To start practicing in Jupyter:
source venv/bin/activate
jupyter lab
```

## 🪐 Exoplanet transit search
`exoplanet_search/` is a CPU-first TESS transit-search pipeline (wotan + BLS fanned out with
Dask, vetting, TOI/planet cross-match) with CUDA and Apple MLX worker backends. See
[exoplanet_search/README.md](exoplanet_search/README.md).
