# pgstats - Population Genetics Statistics Toolkit

A comprehensive toolkit for calculating population genetics statistics from genomic data. Built on [sgkit](https://github.com/pystatgen/sgkit).

## Installation

### Requirements

- Python >=3.8, <3.13
- Dependencies: sgkit, numba, numpy, xarray, zarr, scipy, pandas, dask, pyranges, bio2zarr

### From Source

```bash
git clone https://github.com/SolomonAPS/pgstats.git
cd pgstats

# Create conda environment with dependencies
conda env create -f environment.yml
conda activate pgstats

# Install package
pip install -e .

# Verify installation
pgstats --version
```

## Quick Start

### Command-Line Interface

```bash
# Genome-wide statistics
pgstats stats input.vcf.gz --output results.csv

# Windowed analysis (100kb windows, 50kb step)
pgstats stats input.vcf.gz --output results.csv \
    --window-size 100000 --step-size 50000

# With callable sites masking
pgstats stats input.vcf.gz --output results.csv \
    --bed non_callable.bed --window-size 50000
```

### Python API

```python
import pgstats as pg
from pgstats.core import GenomicDataset, WindowConfig

# Load VCF data
dataset = GenomicDataset("input.vcf.gz")

# Configure windowing
dataset.window_config = WindowConfig(
    window_size=100000,
    step_size=50000,
    min_variants=10
)

# Create windows and calculate statistics
dataset.create_windows()
results = dataset.calculate_windowed_stats(
    stats=['tajima_d', 'theta_pi', 'theta_w']
)

# Save results
dataset.save_results("results.csv")
```
