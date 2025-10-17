# Many-Stats

A comprehensive statistical genomics toolkit built on top of [sgkit](https://github.com/pystatgen/sgkit) with custom statistical functions implemented using Numba for high performance.

**Built with:**
- [sgkit](https://github.com/pystatgen/sgkit) - Scalable genetics toolkit for large-scale genomic data
- [sgkit documentation](https://pystatgen.github.io/sgkit/latest/) - Official sgkit documentation

## Features

- **Built on sgkit**: Leverages sgkit's efficient handling of genetic data and cloud-native formats
- **Numba-optimized**: Custom statistical functions implemented with Numba JIT compilation for high performance
- **Population genetics statistics**: Tajima's D, Fu and Li's statistics, and other neutrality tests
- **Theta estimators**: theta_pi (θπ - nucleotide diversity), theta_w (θw - Watterson's), theta_h (θh - Fay and Wu's)
- **Selection tests**: McDonald-Kreitman test, Hudson-Kreitman-Aguade test
- **Scalable**: Designed to work with large datasets using Dask for parallel processing
- **Extensible**: Easy to add custom statistics using Numba

## Installation

### From Source

```bash
git clone https://github.com/yourusername/many-stats.git
cd many-stats
pip install -e .
```

### Dependencies

- sgkit >= 0.15.0
- numba >= 0.56.0
- numpy >= 1.20.0
- dask >= 2021.0.0
- xarray >= 0.20.0
- zarr >= 2.10.0
- scipy >= 1.7.0
- pandas >= 1.3.0
- bio2zarr >= 0.3.0 (for VCF loading)

## Quick Start

### Loading VCF Data

```python
import many_stats as ms
from many_stats.io.loaders import load_vcf_simple

# Load VCF data using vcf2zarr
ds = load_vcf_simple("your_data.vcf.gz")

# Calculate population statistics
ds = ms.stats.tajima_d(ds)
ds = ms.stats.fu_li_d(ds)
ds = ms.stats.fu_li_f(ds)

# Calculate diversity statistics
ds = ms.stats.nucleotide_diversity(ds)
ds = ms.stats.theta_w(ds)

# Save results
ms.io.save_results(ds, "results.zarr")
```

### Using Simulated Data

```python
import many_stats as ms
import sgkit as sg

# Simulate genetic data
ds = sg.simulate_genotype_call_dataset(n_variant=1000, n_sample=50, missing_pct=0.1)

# Calculate population statistics
ds = ms.stats.tajima_d(ds)
ds = ms.stats.fu_li_d(ds)
ds = ms.stats.fu_li_f(ds)

# Calculate diversity statistics
ds = ms.stats.nucleotide_diversity(ds)
ds = ms.stats.theta_w(ds)

# Save results
ms.io.save_results(ds, "results.zarr")
```

## Available Statistics

### Population Statistics

- **Tajima's D**: Tests for neutrality and demographic effects
- **Fu and Li's D***: Tests for neutrality using folded spectrum
- **Fu and Li's F***: Tests for neutrality using pairwise differences

### Diversity Statistics

- **theta_pi (θπ)**: Nucleotide diversity based on average pairwise differences
- **theta_w (θw)**: Watterson's estimator based on number of segregating sites
- **theta_h (θh)**: Fay and Wu's estimator based on high-frequency derived alleles

### Selection Statistics

- **McDonald-Kreitman test**: Compares polymorphism and divergence patterns
- **Hudson-Kreitman-Aguade test**: Tests for selection across loci

### Missing Data Handling

All statistics properly handle missing data with per-site sample size adjustments:
- **Theta estimators** use actual per-site sample sizes for accurate diversity estimates
- **Neutrality tests** (Tajima's D, Fu & Li's, Zeng's E) use per-site theta values but calculate variance components with max(n) across variants, following scikit-allel and sgkit conventions
- See [Callable Sites Guide](examples/callable_sites_guide.md) for details

## Usage Examples

### Basic Analysis

```python
import many_stats as ms
import sgkit as sg

# Load data
ds = sg.load_dataset("data.zarr")

# Validate dataset
validation = ms.utils.validate_dataset(ds)
print(f"Dataset valid: {validation['valid']}")

# Check missing data
missing_stats = ms.utils.check_missing_data(ds)
print(f"Missing rate: {missing_stats['missing_rate']:.3f}")

# Calculate statistics
ds = ms.stats.tajima_d(ds)
ds = ms.stats.nucleotide_diversity(ds)

# Save results
ms.io.save_results(ds, "analysis_results.csv", format="csv")
```

### Working with Large Datasets

```python
from dask.distributed import Client

# Start Dask cluster for parallel processing
client = Client()

# Your analysis will automatically use the cluster
ds = sg.load_dataset("large_dataset.zarr")
results = ms.stats.tajima_d(ds)

# Close cluster when done
client.close()
```

### Custom Statistics

```python
import numba

@numba.njit
def custom_statistic(variant_matrix):
    """Example custom statistic implementation."""
    n_variants, n_samples = variant_matrix.shape
    result = np.zeros(n_variants)
    
    for i in range(n_variants):
        # Your custom calculation here
        result[i] = np.sum(variant_matrix[i, :]) / n_samples
    
    return result

# Use with many-stats utilities
binary_matrix = ms.utils.convert_to_variant_matrix(ds)
custom_values = custom_statistic(binary_matrix)
```

## Data Formats

Many-stats supports various genetic data formats through sgkit:

- **Zarr**: Cloud-native format for large datasets
- **PLINK**: Standard format (.bed, .bim, .fam files)
- **VCF**: Via bio2zarr conversion (vcf2zarr)

### VCF Loading with vcf2zarr

The package includes built-in support for loading VCF files using the bio2zarr package:

```python
from many_stats.io.loaders import load_vcf_simple, check_bio2zarr_available

# Check if bio2zarr is available
if check_bio2zarr_available():
    # Load VCF data
    ds = load_vcf_simple("your_data.vcf.gz")
else:
    print("Install bio2zarr: pip install bio2zarr")
```

The VCF loading process uses bio2zarr's Python API:
1. Convert VCF directly to Zarr format using `bio2zarr.vcf.convert()`
2. Load the resulting Zarr dataset with sgkit

This approach provides reliable conversion and efficient loading of VCF data with better error handling and progress reporting.

## Performance

Many-stats is optimized for performance:

- **Numba JIT compilation**: Critical functions are compiled for speed
- **Dask integration**: Automatic parallelization for large datasets
- **Memory efficient**: Works with datasets larger than RAM
- **Cloud-ready**: Supports cloud storage formats

## Contributing

Contributions are welcome! Please see our contributing guidelines for details.

### Development Setup

```bash
git clone https://github.com/yourusername/many-stats.git
cd many-stats
pip install -e ".[dev]"
pytest  # Run tests
```

## Citation

If you use many-stats in your research, please cite:

```bibtex
@software{many_stats,
  title={Many-Stats: A comprehensive statistical genomics toolkit},
  author={Your Name},
  year={2024},
  url={https://github.com/yourusername/many-stats}
}
```

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Acknowledgments

- Built on top of [sgkit](https://github.com/pystatgen/sgkit) - Special thanks to the sgkit team for creating an excellent foundation
- Statistical formulas from Wakeley (2009) *Coalescent Theory: An Introduction*
- Inspired by pylibseq and other population genetics tools

## Additional Resources

### sgkit Resources
- [sgkit GitHub Repository](https://github.com/pystatgen/sgkit)
- [sgkit Documentation](https://pystatgen.github.io/sgkit/latest/)
- [sgkit Tutorials](https://pystatgen.github.io/sgkit/latest/tutorials.html)
- [bio2zarr Documentation](https://github.com/sgkit-dev/bio2zarr) - For VCF conversion

### Related Tools
- [scikit-allel](https://scikit-allel.readthedocs.io/) - Python package for exploring and analyzing genetic variation data
- [tskit](https://tskit.dev/) - Tree sequence toolkit
- [msprime](https://tskit.dev/msprime/) - Population genetics simulator

### Learning Resources
- Wakeley J. (2009) *Coalescent Theory: An Introduction*. Roberts & Company Publishers.
- Nielsen R. & Slatkin M. (2013) *An Introduction to Population Genetics: Theory and Applications*. Sinauer Associates.
