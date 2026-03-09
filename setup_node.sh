#!/bin/bash

# --- AstroCompute Lab: Node Setup Script ---
# Run this on every device (Master and Workers) to ensure environments match.

PROJECT_DIR="$HOME/astro-compute-lab"
MASTER_IP="<MASTER_IP>" # <MASTER_HOSTNAME>

echo "🚀 Starting AstroCompute environment setup..."

# 1. Check for Python 3
if ! command -v python3 &> /dev/null; then
    echo "❌ Python 3 not found. Please install it first (sudo apt install python3-venv)."
    exit 1
fi

# 2. Create Project Directory if it doesn't exist
mkdir -p "$PROJECT_DIR"
cd "$PROJECT_DIR"

# 3. Setup Virtual Environment
if [ ! -d "venv" ]; then
    echo "📦 Creating virtual environment..."
    python3 -m venv venv
else
    echo "✅ Virtual environment already exists."
fi

# 4. Install/Update Dependencies
source venv/bin/activate
if [ -f "requirements.txt" ]; then
    echo "📥 Installing astronomical libraries (Astropy, Dask, etc.)..."
    pip install --upgrade pip
    pip install -r requirements.txt
else
    echo "⚠️ requirements.txt not found. Skipping library install."
fi

echo "--------------------------------------------------"
echo "✅ Setup Complete for $(hostname)!"
echo "--------------------------------------------------"
echo "🛰️ To start this node as a WORKER, run:"
echo "source $PROJECT_DIR/venv/bin/activate"
echo "dask-worker $MASTER_IP:8786"
echo "--------------------------------------------------"
echo "🔭 To start this node as the MASTER (Scheduler), run:"
echo "source $PROJECT_DIR/venv/bin/activate"
echo "dask-scheduler"
echo "--------------------------------------------------"
