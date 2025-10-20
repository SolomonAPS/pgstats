# Environment Recreation Guide

This guide shows how to recreate the `many-stats` conda environment.

## Method 1: Full Environment Recreation (Recommended)

Use the complete environment file with exact versions:

```bash
# Create environment from complete specification
conda env create -f environment.yml

# Activate the environment
conda activate many-stats

# Install the package in development mode
pip install -e .
```

## Method 2: Minimal Environment Recreation

Use the minimal environment file for a clean setup:

```bash
# Create environment from minimal specification
conda env create -f environment-minimal.yml

# Activate the environment
conda activate many-stats

# Install the package in development mode
pip install -e .
```

## Method 3: Manual Recreation

Create environment manually and install dependencies:

```bash
# Create new environment
conda create -n many-stats python=3.13

# Activate environment
conda activate many-stats

# Install core dependencies
conda install -c conda-forge -c bioconda \
    numpy pandas scipy numba xarray dask zarr \
    matplotlib seaborn jupyter pytest black flake8 mypy

# Install additional dependencies via pip
pip install sgkit bio2zarr

# Install the package in development mode
pip install -e .
```

## Method 4: Using Requirements File

```bash
# Create environment
conda create -n many-stats python=3.13

# Activate environment
conda activate many-stats

# Install dependencies
pip install -r requirements-clean.txt

# Install the package in development mode
pip install -e .
```

## Verification

After recreation, verify the installation:

```bash
# Test imports
python -c "from many_stats.core import GenomicDataset; print('✓ Import successful')"

# Test CLI
many-stats stats --help

# Run tests
pytest tests/
```

## Environment Files

- **`environment.yml`**: Complete environment with exact versions (310+ packages)
- **`environment-minimal.yml`**: Minimal environment with essential packages only
- **`requirements-clean.txt`**: Clean pip requirements without local paths
- **`requirements.txt`**: Raw pip freeze output (not recommended for recreation)

## Notes

- The complete `environment.yml` ensures exact reproducibility
- The minimal `environment-minimal.yml` is faster to install but may have version differences
- Always install the package in development mode (`pip install -e .`) for local development
- The environment includes development tools (pytest, black, flake8, mypy) for code quality
