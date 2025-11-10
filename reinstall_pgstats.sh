#!/bin/bash
# Script to reinstall pgstats with latest changes

echo "Reinstalling pgstats from local source..."
echo "Current directory: $(pwd)"

# Uninstall existing version
pip uninstall -y pgstats

# Install from current directory in development mode
pip install -e .

echo ""
echo "Installation complete! Verify with:"
echo "  python -c 'import pgstats; print(pgstats.__file__)'"
echo "  pgstats --version"

