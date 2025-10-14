#!/bin/bash
# Setup script for many-stats development environment

set -e  # Exit on any error

echo "🧬 Setting up many-stats development environment..."

# Check if conda is available
if ! command -v conda &> /dev/null; then
    echo "❌ Conda is not installed or not in PATH"
    echo "Please install Miniconda or Anaconda first:"
    echo "https://docs.conda.io/en/latest/miniconda.html"
    exit 1
fi

echo "✓ Conda is available"

# Create the environment
echo "📦 Creating conda environment from environment.yml..."
conda env create -f environment.yml

echo "✓ Environment created successfully"

# Activate the environment and install the package in development mode
echo "🔧 Installing many-stats in development mode..."
conda run -n many-stats pip install -e .

echo "✓ Package installed in development mode"

# Test the installation
echo "🧪 Testing the installation..."
conda run -n many-stats python -c "
import many_stats as ms
from many_stats.io.loaders import load_vcf_simple, check_bio2zarr_available
import sgkit as sg
import bio2zarr.vcf as v2z
print('✓ All imports successful!')
print(f'✓ bio2zarr available: {check_bio2zarr_available()}')
"

echo ""
echo "🎉 Setup complete! To activate the environment, run:"
echo "conda activate many-stats"
echo ""
echo "To test VCF loading, run:"
echo "conda run -n many-stats python test_vcf_api.py"


