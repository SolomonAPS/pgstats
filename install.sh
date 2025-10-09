#!/bin/bash

# Installation script for many-stats

echo "Installing many-stats..."

# Check if conda is available
if command -v conda &> /dev/null; then
    echo "Using conda to install dependencies..."
    conda env create -f environment.yml
    conda activate many-stats
    pip install -e .
else
    echo "Using pip to install dependencies..."
    pip install -r requirements.txt
    pip install -e .
fi

echo "Installation complete!"
echo ""
echo "To test the installation, run:"
echo "python -c \"import many_stats; print('Many-stats installed successfully!')\""
