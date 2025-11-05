# pgstats - Population Genetics Statistics Toolkit

A comprehensive toolkit for calculating population genetics statistics from genomic data. Built on [sgkit](https://github.com/pystatgen/sgkit) with high-performance implementations using Numba JIT compilation.

## Features    

- **Comprehensive statistics**: Neutrality tests, diversity measures, LD statistics, and haplotype statistics
- **Flexible windowing**: Sliding windows, fixed windows, or genome-wide analysis
- **Missing data handling**: Proper per-site sample size adjustments
- **Callable sites masking**: Account for non-callable regions when calculating theta estimators
- **Multi-region analysis**: Analyze multiple genomic regions efficiently
- **High performance**: Numba-compiled core functions for speed
- **Multiple input formats**: VCF, Zarr, or sgkit datasets

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

**Note**: Some dependencies like `pyranges` are best installed via conda rather than pip.

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

## Available Statistics

### Site Frequency Spectrum (SFS) Statistics

#### Theta Estimators (Diversity Measures)

| Statistic | Description | Formula | Reference |
|-----------|-------------|---------|-----------|
| **theta_pi** (θπ) | Nucleotide diversity based on pairwise differences | θπ = (1/C(n,2)) × Σ i(n-i)ξᵢ | Walsh & Lynch 2018, Eq. 9.4 |
| **theta_w** (θw) | Watterson's estimator based on segregating sites | θw = S/a₁ | Walsh & Lynch 2018, Eq. 9.5 |
| **theta_h** (θh) | Fay-Wu's estimator emphasizing high-frequency alleles | θh = (1/C(n,2)) × Σ i²ξᵢ | Walsh & Lynch 2018, Eq. 9.6 |
| **theta_l** (θL) | Zeng's estimator emphasizing intermediate frequencies | θL = (1/(n-1)) × Σ iξᵢ | Walsh & Lynch 2018, Eq. 9.7 |

Where:
- `S` = number of segregating sites
- `ξᵢ` = number of sites with i derived alleles
- `n` = sample size
- `a₁` = Σ(1/i) for i=1 to n-1 (harmonic number)
- `C(n,2)` = n(n-1)/2

#### Neutrality Tests

| Statistic | Description | Formula | Reference |
|-----------|-------------|---------|-----------|
| **tajima_d** | Tests for deviation from neutral evolution | D = (θπ - θw) / √Var(θπ - θw) | Walsh & Lynch 2018, Eq. 9.8 |
| **fu_li_d** | Tests using singleton frequency (folded) | D* = (S/a₁ - ((n-1)/n)η₁) / √Var(D*) | Walsh & Lynch 2018, Eq. 9.26b |
| **fu_li_f** | Combines pairwise differences and singletons | F* = (θπ - ((n-1)/n)η₁) / √Var(F*) | Walsh & Lynch 2018, Eq. 9.26e |
| **zeng_e** | Tests for selective sweeps | E = (θL - θw) / √Var(θL - θw) | Walsh & Lynch 2018, Eq. 9.28c |
| **fay_wu_h** | Tests for selection using high-frequency alleles | H = (θπ - θh) / √Var(θπ - θh) | Walsh & Lynch 2018, Eq. 9.27b |

**Unfolded versions** (require ancestral state information):
- `fu_li_d_unfolded`: Uses derived allele frequencies
- `fu_li_f_unfolded`: Uses derived allele frequencies

### Linkage Disequilibrium (LD) Statistics

| Statistic | Description | Reference |
|-----------|-------------|-----------|
| **ld_d** | Coefficient of linkage disequilibrium | Walsh & Lynch 2018, Eq. 9.1 |
| **ld_dprime** | Standardized D (D') | Walsh & Lynch 2018, Eq. 9.2 |
| **ld_r2** | Squared correlation coefficient (r²) | Walsh & Lynch 2018, Eq. 9.3 |
| **omega_statistic** | Detects recombination breakpoints | Kim & Nielsen 2004 |

### Haplotype Statistics

| Statistic | Description | Reference |
|-----------|-------------|-----------|
| **haplotype_diversity** | Probability two haplotypes differ | Walsh & Lynch 2018 |
| **garud_h1** | Haplotype homozygosity | Garud et al. 2015 |
| **garud_h12** | Combined frequency of top 2 haplotypes | Garud et al. 2015 |
| **garud_h123** | Combined frequency of top 3 haplotypes | Garud et al. 2015 |
| **garud_h2_h1** | Ratio of H2 to H1 | Garud et al. 2015 |

## Data Input Formats

pgstats supports multiple input formats:

### 1. VCF Files (.vcf, .vcf.gz)

```bash
# Direct VCF input (automatically converted to Zarr)
pgstats stats input.vcf.gz --output results.csv

# Keep Zarr file for faster re-runs
pgstats stats input.vcf.gz --output results.csv --keep-zarr

# Specify custom Zarr directory (useful on clusters)
pgstats stats input.vcf.gz --output results.csv \
    --keep-zarr --zarr-dir /scratch/zarr_cache
```

**How it works**: VCF files are automatically converted to Zarr format using `bio2zarr` for efficient processing. The Zarr file is cached by default next to your output file (or in `--zarr-dir` if specified) and reused on subsequent runs with `--keep-zarr`.

### 2. Zarr Format (.zarr)

```bash
# Use pre-converted Zarr (fastest)
pgstats stats input.zarr --output results.csv
```

**Advantages**: Zarr is a chunked, compressed array storage format optimized for parallel I/O. Pre-converting large VCF files to Zarr significantly speeds up repeated analyses.

### 3. Python API with sgkit Datasets

```python
import sgkit as sg
from pgstats.core import GenomicDataset

# Load from VCF
ds = sg.load_dataset("input.zarr")

# Or simulate data for testing
ds = sg.simulate_genotype_call_dataset(n_variant=1000, n_sample=50)

# Create GenomicDataset
genomic_ds = GenomicDataset(data_source=ds)
```

## Missing Data Handling

pgstats properly handles missing data at multiple levels:

### 1. Per-Site Sample Sizes

All theta estimators use **actual per-site sample sizes** to account for missing data:

```python
# For each variant, calculate n_i = number of non-missing samples
# Then use n_i for that site's contribution to theta
```

This ensures accurate diversity estimates even with heterogeneous missing data patterns.

### 2. Filtering Variants by Missing Data

Use `--max-missing` to exclude variants with too much missing data:

```bash
# Remove variants with >10% missing data
pgstats stats input.vcf.gz --output results.csv --max-missing 0.1
```

**What it does**: Filters out entire variants (SNPs) where the proportion of missing genotypes exceeds the threshold. This happens **before** window creation and statistics calculation.

### 3. Filtering Windows by Variant Count

Use `--min-variants` to exclude windows with too few variants:

```bash
# Only keep windows with ≥10 variants
pgstats stats input.vcf.gz --output results.csv \
    --window-size 50000 --min-variants 10
```

**What it does**: After creating windows, filters out windows that don't have enough segregating sites for reliable statistics. This happens **after** statistics are calculated.

**Key Difference**:
- `--max-missing`: Filters **variants** (rows) based on missing genotypes
- `--min-variants`: Filters **windows** based on number of variants

Both can be used together:
```bash
pgstats stats input.vcf.gz --output results.csv \
    --max-missing 0.1 \    # Remove variants with >10% missing data
    --min-variants 10 \     # Remove windows with <10 variants
    --window-size 50000
```

## Callable Sites Masking

When calculating theta estimators, the sequence length `L` must account for regions where variants **cannot** be called (e.g., low coverage, repetitive regions).

### Why This Matters

Theta estimators are normalized by sequence length:
- θπ = (sum of pairwise differences) / L
- θw = S / (a₁ × L)

If `L` includes non-callable regions, theta will be underestimated.

### Using BED Files for Masking

```bash
# Provide non-callable regions (default behavior)
pgstats stats input.vcf.gz --output results.csv \
    --bed non_callable_regions.bed \
    --window-size 50000

# Or provide callable regions only
pgstats stats input.vcf.gz --output results.csv \
    --bed callable_regions.bed \
    --bed-format callable \
    --window-size 50000
```

**BED file format** (0-based, half-open intervals):
```
chr1    1000    2000
chr1    5000    6000
chr2    0       1000
```

**How it works**:
1. For each window, pgstats calculates the overlap with BED regions using `pyranges`
2. Callable length `L` = window_size - (non-callable bases in window)
3. Theta estimators are divided by `L` instead of window_size

**Example**:
- Window: chr1:10000-20000 (10kb)
- Non-callable regions in window: 2kb
- Callable length L = 8kb
- θw = S / (a₁ × 8000) instead of S / (a₁ × 10000)

## Windowing Options

### Fixed Windows

```bash
# 100kb non-overlapping windows
pgstats stats input.vcf.gz --output results.csv \
    --window-size 100000
```

### Sliding Windows

```bash
# 100kb windows with 50kb step (50% overlap)
pgstats stats input.vcf.gz --output results.csv \
    --window-size 100000 \
    --step-size 50000
```

### Single Region Analysis

```bash
# Analyze specific region
pgstats stats input.vcf.gz --output results.csv \
    --region chr1:1000000-5000000 \
    --window-size 50000

# Whole chromosome
pgstats stats input.vcf.gz --output results.csv \
    --region chr1 \
    --window-size 100000
```

### Multi-Region Analysis

Create a BED file with regions of interest:

```
# candidate_genes.bed
chr1    1000000    1050000    gene1
chr1    2000000    2100000    gene2
chr2    500000     550000     gene3
```

Then analyze all regions:

```bash
pgstats stats input.vcf.gz --output results.csv \
    --regions-file candidate_genes.bed \
    --window-size 10000 \
    --bed non_callable.bed
```

**Output**: Results for all windows across all regions, with columns indicating which region each window belongs to.

### Window Filtering

```bash
# Only keep windows with ≥10 variants
pgstats stats input.vcf.gz --output results.csv \
    --window-size 50000 \
    --min-variants 10
```

**Tip**: Set `--min-variants 0` to keep all windows and filter later based on the `n_variants` column in the output.

## Output Format

### CSV Output (default)

```csv
contig,start,end,n_variants,callable_length,tajima_d,theta_pi,theta_w,fu_li_d,fu_li_f,zeng_e,fay_wu_h,theta_h,theta_l,ld_r2,omega_statistic,haplotype_diversity,garud_h1,garud_h12,garud_h123,garud_h2_h1
chr1,0,100000,45,98500,0.523,0.00234,0.00198,0.234,-0.123,0.456,0.123,0.00189,0.00201,0.234,12.3,0.876,0.234,0.345,0.456,0.567
chr1,50000,150000,52,99200,0.234,0.00198,0.00176,0.123,-0.234,0.345,0.234,0.00167,0.00189,0.198,15.6,0.823,0.198,0.298,0.398,0.498
```

**Columns**:
- `contig`: Chromosome/contig name
- `start`: Window start position (0-based)
- `end`: Window end position (exclusive)
- `n_variants`: Number of segregating sites in window
- `callable_length`: Effective sequence length (accounts for BED masking)
- Statistics columns: One per requested statistic

### Other Formats

```bash
# TSV output
pgstats stats input.vcf.gz --output results.tsv --format tsv

# Parquet output (efficient for large datasets)
pgstats stats input.vcf.gz --output results.parquet --format parquet
```

## Complete Examples

### Example 1: Basic Diversity Scan

```bash
pgstats stats population.vcf.gz \
    --output diversity_scan.csv \
    --stats theta_pi theta_w \
    --window-size 100000 \
    --step-size 50000 \
    --max-missing 0.1
```

### Example 2: Selection Scan with Quality Filters

```bash
pgstats stats population.vcf.gz \
    --output selection_scan.csv \
    --stats tajima_d fu_li_d fu_li_f zeng_e fay_wu_h \
    --window-size 50000 \
    --step-size 25000 \
    --bed non_callable_10x.bed \
    --max-missing 0.05 \
    --min-variants 20
```

### Example 3: Multi-Region Candidate Gene Analysis

```bash
# Create regions file
cat > candidate_genes.bed << EOF
chr2L    5000000    5100000    gene1
chr2L    8000000    8050000    gene2
chr3R    2000000    2100000    gene3
EOF

# Run analysis
pgstats stats genome.vcf.gz \
    --output candidate_genes_stats.csv \
    --regions-file candidate_genes.bed \
    --window-size 10000 \
    --bed non_callable.bed \
    --max-missing 0.1 \
    --min-variants 5
```

### Example 4: Haplotype Analysis

```bash
pgstats stats phased_data.vcf.gz \
    --output haplotype_stats.csv \
    --stats garud_h1 garud_h12 garud_h123 garud_h2_h1 haplotype_diversity \
    --window-size 20000 \
    --min-variants 10
```

## Performance Tips

1. **Use Zarr for repeated analyses**:
   ```bash
   # First run: convert and cache
   pgstats stats large.vcf.gz --output run1.csv --keep-zarr
   
   # Subsequent runs: much faster
   pgstats stats large.vcf.gz --output run2.csv --keep-zarr --region chr1
   ```

2. **Analyze chromosomes separately** for parallelization:
   ```bash
   for chr in chr1 chr2 chr3; do
       pgstats stats genome.vcf.gz \
           --output ${chr}_stats.csv \
           --region ${chr} \
           --window-size 100000 &
   done
   wait
   ```

3. **Use larger windows** for faster computation:
   ```bash
   # Faster but lower resolution
   pgstats stats data.vcf.gz --output results.csv --window-size 500000
   ```

4. **Pre-filter VCF** to reduce dataset size:
   ```bash
   bcftools view -i 'F_MISSING<0.1' input.vcf.gz | \
       pgstats stats - --output results.csv
   ```

## References

### Textbooks
- Walsh, B. & Lynch, M. (2018). *Evolution and Selection of Quantitative Traits*. Oxford University Press.
- Wakeley, J. (2009). *Coalescent Theory: An Introduction*. Roberts & Company Publishers.

### Papers
- Tajima, F. (1989). Statistical method for testing the neutral mutation hypothesis by DNA polymorphism. *Genetics*, 123(3), 585-595.
- Fu, Y. X., & Li, W. H. (1993). Statistical tests of neutrality of mutations. *Genetics*, 133(3), 693-709.
- Fay, J. C., & Wu, C. I. (2000). Hitchhiking under positive Darwinian selection. *Genetics*, 155(3), 1405-1413.
- Kim, Y., & Nielsen, R. (2004). Linkage disequilibrium as a signature of selective sweeps. *Genetics*, 167(3), 1513-1524.
- Zeng, K., Fu, Y. X., Shi, S., & Wu, C. I. (2006). Statistical tests for detecting positive selection by utilizing high-frequency variants. *Genetics*, 174(3), 1431-1439.
- Garud, N. R., Messer, P. W., Buzbas, E. O., & Petrov, D. A. (2015). Recent selective sweeps in North American Drosophila melanogaster show signatures of soft sweeps. *PLoS Genetics*, 11(2), e1005004.

### Software
- sgkit: https://github.com/pystatgen/sgkit
- bio2zarr: https://github.com/sgkit-dev/bio2zarr
- pyranges: https://github.com/biocore-ntnu/pyranges

## Citation

If you use pgstats in your research, please cite:

```bibtex
@software{pgstats,
  title={pgstats: Population Genetics Statistics Toolkit},
  author={Solomon Sloat},
  year={2024},
  url={https://github.com/SolomonAPS/pgstats}
}
```

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Support

- **Documentation**: See `docs/CLI_GUIDE.md` for detailed CLI usage
- **Issues**: https://github.com/SolomonAPS/pgstats/issues
- **Examples**: See `examples/` directory for Jupyter notebooks and guides
