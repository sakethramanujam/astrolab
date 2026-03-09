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

## 🛠️ Getting Started
```bash
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
```
