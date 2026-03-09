#!/bin/bash

# --- AstroCompute Lab: Unified Setup Script ---
# This script handles environment creation and network configuration.

set -e # Exit on any error

PROJECT_DIR=$(pwd)
MASTER_IP="<MASTER_IP>" # The IP of this machine (<MASTER_OS_HOSTNAME>/NUC)
VENV_NAME="venv"

echo "🌌 Initializing AstroCompute Lab on $(hostname)..."

# 1. System Check
if ! command -v python3 &> /dev/null; then
    echo "❌ Error: Python 3 is not installed. Please install it with 'sudo apt install python3-venv'."
    exit 1
fi

# 2. Virtual Environment Setup
if [ ! -d "$VENV_NAME" ]; then
    echo "📦 Creating virtual environment..."
    python3 -m venv "$VENV_NAME"
else
    echo "✅ Virtual environment already exists."
fi

# 3. Dependency Installation
echo "🔌 Activating environment and installing dependencies..."
source "$VENV_NAME"/bin/activate
pip install --upgrade pip
if [ -f "requirements.txt" ]; then
    pip install -r requirements.txt
else
    echo "⚠️ Warning: requirements.txt not found. Skipping library install."
fi

# 4. Success and Instructions
echo "--------------------------------------------------"
echo "✅ Setup Complete for $(hostname)!"
echo "--------------------------------------------------"
echo "🚀 TO START THE LAB (Jupyter):"
echo "   source $VENV_NAME/bin/activate && jupyter lab"
echo ""
echo "🔭 TO USE AS THE MASTER (NUC/Main):"
echo "   source $VENV_NAME/bin/activate && dask-scheduler"
echo ""
echo "🛰️ TO USE AS A WORKER (Pi/WSL/Other):"
echo "   source $VENV_NAME/bin/activate && dask-worker $MASTER_IP:8786"
echo "--------------------------------------------------"
