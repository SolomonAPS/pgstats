# Callable Sites Guide for many-stats

This guide explains how to work with callable sites and mask non-callable regions in your population genetics analyses.

## Table of Contents

1. [Overview](#overview)
2. [BED File Format](#bed-file-format)
3. [Basic Usage](#basic-usage)
4. [Advanced Options](#advanced-options)
5. [Best Practices](#best-practices)
6. [Troubleshooting](#troubleshooting)

## Overview

Callable sites are genomic positions where you have sufficient sequencing coverage and quality to reliably call genotypes. Non-callable sites should be excluded from population genetics analyses to avoid biased statistics.

### Why Callable Sites Matter

- **Avoids bias**: Missing data at non-callable sites can bias diversity estimates
- **Improves accuracy**: Statistics are more reliable when based on high-quality data
- **Required for some analyses**: Proper callable sites handling is essential for accurate neutrality tests

### How many-stats Handles Callable Sites

many-stats uses a simple but effective approach:
1. Load callable regions from a BED file
2. Mark variants outside callable regions as missing data (-1)
3. Statistics automatically account for missing data

### Missing Data Handling in Statistics

Different statistics handle per-site missing data differently:

**Per-Site Theta Estimators** (θπ, θw, θh, θL):
- Use the actual per-site sample size for calculations
- Each variant is calculated using only non-missing samples at that site
- Provides the most accurate estimates when missingness varies across sites

**Window-Level Neutrality Tests** (Tajima's D, Fu & Li's D/F, Zeng's E):
- Per-site theta values use actual per-site sample sizes
- Harmonic numbers (a₁, a₂, b₁, b₂) and variance components use max(n) across all variants
- This approach follows scikit-allel and sgkit conventions
- Maintains theoretical consistency with the assumption of constant n in the formulas
- Balances accuracy in theta calculations with valid statistical testing

## BED File Format

### Standard BED Format

BED files use **0-based, half-open intervals**:
- **Start**: 0-based position (inclusive)
- **End**: 0-based position (exclusive)

```
chr1    0       1000
chr1    2000    3000
chr2    0       5000
```

This means:
- `chr1 0 1000` includes positions 1-1000 (1-based VCF coordinates)
- `chr1 2000 3000` includes positions 2001-3000 (1-based VCF coordinates)

### Creating a BED File

#### From Sequencing Coverage

```bash
# Using bedtools
bedtools genomecov -ibam aligned.bam -bg | \
  awk '$4 >= 10' > callable_10x.bed

# Using samtools
samtools depth -a aligned.bam | \
  awk '$3 >= 10 {print $1"\t"$2-1"\t"$2}' > callable_10x.bed
```

#### From Multiple BAM Files

```bash
# Merge coverage from multiple samples
bedtools multicov -bams *.bam -bed genome.bed | \
  awk '{sum=0; for(i=4;i<=NF;i++) sum+=$i; if(sum/NF >= 10) print $1"\t"$2"\t"$3}' \
  > callable_10x_merged.bed
```

#### Manual Creation

```bash
# Create callable sites BED file
cat > callable_sites.bed << EOF
chr1    0       1000000
chr1    2000000 3000000
chr2    0       500000
EOF
```

## Basic Usage

### Simple Callable Sites Filtering

```python
from many_stats.core.dataset import GenomicDataset, CallableSitesConfig

# Configure callable sites
callable_config = CallableSitesConfig(
    bed_file="callable_sites.bed"
)

# Load data with callable sites masking
genomic_ds = GenomicDataset(
    data_source="my_data.vcf.gz",
    callable_config=callable_config
)

# Check how many callable sites
print(f"Callable sites: {genomic_ds.callable_sites:,}")
```

### With Missing Data Filtering

Combine callable sites with missing data filtering:

```python
callable_config = CallableSitesConfig(
    bed_file="callable_sites.bed",
    max_missing=0.1  # Allow up to 10% missing data
)

genomic_ds = GenomicDataset(
    data_source="my_data.vcf.gz",
    callable_config=callable_config
)

# Filter variants with too much missing data
genomic_ds.filter_missing_data()
```

### Complete Workflow

```python
from many_stats.core.dataset import (
    GenomicDataset, 
    CallableSitesConfig, 
    WindowConfig
)

# 1. Configure callable sites and windowing
callable_config = CallableSitesConfig(
    bed_file="callable_sites.bed",
    max_missing=0.2
)

window_config = WindowConfig(
    window_size=100000,
    step_size=50000,
    min_variants=10
)

# 2. Load data
genomic_ds = GenomicDataset(
    data_source="population.vcf.gz",
    callable_config=callable_config,
    window_config=window_config
)

# 3. Filter and analyze
genomic_ds.filter_missing_data()
genomic_ds.create_windows()

# 4. Calculate statistics
window_stats = genomic_ds.calculate_windowed_stats(
    stats=['tajima_d', 'nucleotide_diversity', 'watterson_theta'],
    use_callable_sites=True
)

# 5. Save results
genomic_ds.save_results("results_with_callable_sites.csv")
```

## Advanced Options

### Strict Quality Filtering

Be more stringent about which sites to include:

```python
callable_config = CallableSitesConfig(
    bed_file="high_quality_sites.bed",
    max_missing=0.05  # Only 5% missing allowed
)
```

### Relaxed Filtering

More permissive for low-coverage data:

```python
callable_config = CallableSitesConfig(
    bed_file="callable_sites.bed",
    max_missing=0.3  # Allow up to 30% missing
)
```

### No BED File (Missing Data Only)

Use only missing data threshold without a BED file:

```python
callable_config = CallableSitesConfig(
    bed_file=None,
    max_missing=0.2
)

genomic_ds = GenomicDataset(
    data_source="my_data.vcf.gz",
    callable_config=callable_config
)

genomic_ds.filter_missing_data()
```

## Best Practices

### 1. Use Appropriate Coverage Thresholds

**High-coverage data (>30x)**:
```bash
# Stricter filtering
bedtools genomecov -ibam aligned.bam -bg | \
  awk '$4 >= 15' > callable_15x.bed
```

**Low-coverage data (5-10x)**:
```bash
# More permissive
bedtools genomecov -ibam aligned.bam -bg | \
  awk '$4 >= 5' > callable_5x.bed
```

### 2. Merge Overlapping Regions

Always merge overlapping BED intervals:

```bash
# Merge overlapping regions
bedtools merge -i callable_raw.bed > callable_merged.bed
```

### 3. Consider All Samples

For population data, ensure sites are callable across all samples:

```python
# Create BED of sites callable in ≥80% of samples
callable_config = CallableSitesConfig(
    bed_file="callable_80percent.bed",
    max_missing=0.2
)
```

### 4. Document Your Filtering

Keep track of filtering parameters:

```python
# Document in your analysis
callable_config = CallableSitesConfig(
    bed_file="callable_10x_merged.bed",  # Descriptive name
    max_missing=0.1
)

# Log the filtering
summary = genomic_ds.get_summary()
print(f"Callable sites: {summary['callable_sites']:,}")
print(f"Total variants: {summary['n_variants']:,}")
print(f"Max missing: {summary['max_missing']:.1%}")
```

### 5. Validate BED File

Check your BED file before analysis:

```bash
# Check format
head callable_sites.bed

# Check for overlaps
bedtools merge -i callable_sites.bed | wc -l
wc -l callable_sites.bed
# If these differ, you have overlaps

# Sort and merge if needed
sort -k1,1 -k2,2n callable_sites.bed | \
  bedtools merge > callable_sites_clean.bed
```

## Troubleshooting

### Issue: BED File Not Found

**Error**: `FileNotFoundError: BED file not found`

**Solution**: Check the file path

```python
from pathlib import Path

bed_path = "callable_sites.bed"
if not Path(bed_path).exists():
    print(f"BED file not found: {bed_path}")
    print(f"Current directory: {Path.cwd()}")
```

### Issue: No Callable Sites Found

**Error**: All sites marked as non-callable

**Possible Causes**:
1. Chromosome naming mismatch (chr1 vs 1)
2. Wrong coordinate system
3. BED file doesn't overlap with VCF

**Solutions**:

```python
# Check chromosome names in VCF
import sgkit as sg
ds = sg.load_dataset("my_data.zarr")
print("VCF contigs:", ds.contig_id.values)

# Check BED file contigs
import pandas as pd
bed_df = pd.read_csv("callable_sites.bed", sep='\t', 
                     names=['chrom', 'start', 'end'])
print("BED contigs:", bed_df['chrom'].unique())

# Fix chromosome naming if needed
# If VCF has "chr1" but BED has "1", update BED:
bed_df['chrom'] = 'chr' + bed_df['chrom'].astype(str)
bed_df.to_csv("callable_sites_fixed.bed", sep='\t', 
              header=False, index=False)
```

### Issue: Too Many Variants Filtered

**Problem**: Most variants are removed after filtering

**Solutions**:

1. Check missing data threshold:
```python
# Before filtering, check missing data distribution
from many_stats.utils.validation import check_missing_data

ds = sg.load_dataset("my_data.zarr")
missing_info = check_missing_data(ds)
print(f"Overall missing rate: {missing_info['missing_rate']:.3f}")
print(f"Variants with missing: {missing_info['variants_with_missing']}")

# Adjust threshold accordingly
callable_config = CallableSitesConfig(
    bed_file="callable_sites.bed",
    max_missing=0.3  # Increase if needed
)
```

2. Verify BED file coverage:
```bash
# Check how much of genome is callable
awk '{sum += $3-$2} END {print sum}' callable_sites.bed
```

### Issue: Coordinate System Confusion

**Problem**: Uncertain if BED is 0-based or 1-based

**Solution**: many-stats expects standard 0-based BED format

```python
# Standard BED (0-based, half-open)
# Position 1000 in VCF (1-based) corresponds to:
# - Start: 999 (0-based, inclusive)
# - End: 1000 (0-based, exclusive)

# Example: Include VCF position 1000
# Correct BED: chr1  999  1000
# Wrong BED:   chr1  1000 1001
```

### Issue: Performance Issues with Large BED Files

**Problem**: Slow loading with millions of BED intervals

**Solutions**:

1. Merge adjacent intervals:
```bash
bedtools merge -i large_callable.bed > merged_callable.bed
```

2. Pre-filter BED to chromosomes of interest:
```bash
# Only keep chr1
grep "^chr1\t" callable_sites.bed > callable_chr1.bed
```

## Complete Example

### Step-by-Step Analysis with Callable Sites

```python
from many_stats.core.dataset import (
    GenomicDataset,
    CallableSitesConfig,
    WindowConfig
)
import pandas as pd
import matplotlib.pyplot as plt

# Step 1: Configure analysis
callable_config = CallableSitesConfig(
    bed_file="callable_10x_merged.bed",
    max_missing=0.15  # 15% missing allowed
)

window_config = WindowConfig(
    window_size=50000,
    step_size=25000,
    min_variants=10
)

# Step 2: Load and process data
print("Loading data...")
genomic_ds = GenomicDataset(
    data_source="population.vcf.gz",
    callable_config=callable_config,
    window_config=window_config
)

# Step 3: Check initial statistics
summary = genomic_ds.get_summary()
print(f"Initial variants: {summary['n_variants']:,}")
print(f"Callable sites: {summary['callable_sites']:,}")

# Step 4: Filter missing data
print("Filtering missing data...")
genomic_ds.filter_missing_data()

# Step 5: Create windows
print("Creating windows...")
genomic_ds.create_windows()
print(f"Created {genomic_ds.get_summary()['n_windows']:,} windows")

# Step 6: Calculate statistics
print("Calculating statistics...")
window_stats = genomic_ds.calculate_windowed_stats(
    stats=['tajima_d', 'nucleotide_diversity', 'watterson_theta'],
    use_callable_sites=True
)

# Step 7: Save results
print("Saving results...")
genomic_ds.save_results("windowed_stats_callable.csv")

# Step 8: Visualize
df = pd.read_csv("windowed_stats_callable.csv")
fig, axes = plt.subplots(3, 1, figsize=(12, 8))

axes[0].plot(df['window_start'], df['tajima_d'])
axes[0].set_ylabel("Tajima's D")
axes[0].axhline(y=0, color='r', linestyle='--', alpha=0.5)

axes[1].plot(df['window_start'], df['nucleotide_diversity'])
axes[1].set_ylabel('π (Nucleotide Diversity)')

axes[2].plot(df['window_start'], df['watterson_theta'])
axes[2].set_ylabel('θw (Watterson\'s Theta)')
axes[2].set_xlabel('Genomic Position (bp)')

plt.tight_layout()
plt.savefig('population_stats_with_callable_sites.png', dpi=300)
print("Done!")
```

## Comparison: With vs Without Callable Sites

```python
# Analysis WITHOUT callable sites filtering
genomic_ds_unfiltered = GenomicDataset(
    data_source="population.vcf.gz",
    window_config=WindowConfig(window_size=50000)
)
genomic_ds_unfiltered.create_windows()
stats_unfiltered = genomic_ds_unfiltered.calculate_windowed_stats(
    stats=['tajima_d', 'nucleotide_diversity']
)

# Analysis WITH callable sites filtering
callable_config = CallableSitesConfig(bed_file="callable_sites.bed")
genomic_ds_filtered = GenomicDataset(
    data_source="population.vcf.gz",
    callable_config=callable_config,
    window_config=WindowConfig(window_size=50000)
)
genomic_ds_filtered.create_windows()
stats_filtered = genomic_ds_filtered.calculate_windowed_stats(
    stats=['tajima_d', 'nucleotide_diversity'],
    use_callable_sites=True
)

# Compare
print("Comparison:")
print(f"Unfiltered - Mean Tajima's D: {np.mean(stats_unfiltered['tajima_d'].values):.3f}")
print(f"Filtered   - Mean Tajima's D: {np.mean(stats_filtered['tajima_d'].values):.3f}")
```

## See Also

### many-stats Documentation
- [Windowing Guide](windowing_guide.md)
- [Statistics Reference](../README.md#available-statistics)
- [CLI Guide](../docs/CLI_GUIDE.md)

### File Format Resources
- [BED File Format Specification](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) - Official UCSC BED format documentation
- [VCF Format](https://samtools.github.io/hts-specs/VCFv4.3.pdf) - VCF specification
- [bedtools Documentation](https://bedtools.readthedocs.io/) - Tools for working with BED files

### sgkit Resources  
- [sgkit Documentation](https://pystatgen.github.io/sgkit/latest/) - Main sgkit documentation
- [sgkit Data Model](https://pystatgen.github.io/sgkit/latest/user_guide.html) - Understanding sgkit datasets
- [bio2zarr](https://github.com/sgkit-dev/bio2zarr) - VCF conversion tool used by many-stats

### Related Tools
- [samtools](http://www.htslib.org/) - For calculating coverage from BAM files
- [bedtools](https://bedtools.readthedocs.io/) - For manipulating BED files
- [GATK CallableLoci](https://gatk.broadinstitute.org/hc/en-us/articles/360037593151-CallableLoci-BETA-) - Identify callable regions

