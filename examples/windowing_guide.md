# Windowing Guide for many-stats

This guide explains how to use the windowing functionality in many-stats for population genetics analyses.

## Table of Contents

1. [Overview](#overview)
2. [Window Types](#window-types)
3. [Basic Usage](#basic-usage)
4. [Advanced Options](#advanced-options)
5. [Best Practices](#best-practices)
6. [Troubleshooting](#troubleshooting)

## Overview

Windowing allows you to calculate population genetics statistics across genomic regions rather than for the entire genome. This is useful for:

- Detecting regions under selection
- Identifying demographic patterns across the genome
- Reducing computational load for large datasets
- Visualizing statistic distributions along chromosomes

## Window Types

### Genome-Wide Analysis

Calculate statistics for the entire genome as a single unit:

```python
from many_stats.core.dataset import GenomicDataset, WindowConfig

# No windowing - genome-wide
window_config = WindowConfig(window_size=None)

genomic_ds = GenomicDataset(
    data_source="my_data.vcf.gz",
    window_config=window_config
)

genomic_ds.create_windows()
genome_stats = genomic_ds.calculate_genome_wide_stats(
    stats=['tajima_d', 'nucleotide_diversity', 'watterson_theta']
)
```

### Position-Based Windows

Create fixed-size windows based on genomic positions:

```python
# 50kb windows, sliding by 25kb
window_config = WindowConfig(
    window_size=50000,  # 50kb windows
    step_size=25000,     # 25kb step (overlapping)
    min_variants=10      # Minimum 10 variants per window
)

genomic_ds = GenomicDataset(
    data_source="my_data.vcf.gz",
    window_config=window_config
)

genomic_ds.create_windows()
```

### Region-Specific Analysis

Analyze specific genomic regions:

```python
# Analyze region from 1Mb to 5Mb
window_config = WindowConfig(
    window_size=100000,  # 100kb windows
    start=1000000,       # Start at 1Mb
    end=5000000,         # End at 5Mb
    min_variants=5
)
```

## Basic Usage

### Step 1: Configure Windows

```python
from many_stats.core.dataset import GenomicDataset, WindowConfig

window_config = WindowConfig(
    window_size=100000,  # 100kb windows
    step_size=50000,     # 50kb step
    min_variants=5       # Filter out windows with <5 variants
)
```

### Step 2: Create GenomicDataset

```python
genomic_ds = GenomicDataset(
    data_source="my_data.vcf.gz",
    window_config=window_config
)
```

### Step 3: Create Windows

```python
genomic_ds.create_windows()

# Check how many windows were created
summary = genomic_ds.get_summary()
print(f"Created {summary['n_windows']} windows")
```

### Step 4: Calculate Statistics

```python
# Calculate windowed statistics
window_stats = genomic_ds.calculate_windowed_stats(
    stats=['tajima_d', 'nucleotide_diversity', 'watterson_theta', 'fay_wu_theta']
)

# Access results
print(f"Windows: {len(window_stats.windows)}")
print(f"Tajima's D values: {window_stats['tajima_d'].values}")
```

### Step 5: Save Results

```python
# Save to CSV
genomic_ds.save_results("windowed_stats.csv", format='csv')

# Or save to TSV
genomic_ds.save_results("windowed_stats.tsv", format='tsv')
```

## Advanced Options

### Overlapping Windows

Create overlapping windows for smoother statistics:

```python
window_config = WindowConfig(
    window_size=100000,  # 100kb windows
    step_size=10000,     # 10kb step = 90% overlap
    min_variants=5
)
```

### Non-Overlapping Windows

Create adjacent, non-overlapping windows:

```python
window_config = WindowConfig(
    window_size=100000,   # 100kb windows
    step_size=100000,     # Same as window_size = no overlap
    min_variants=5
)
```

### Filtering Windows

Filter windows by minimum variant count:

```python
window_config = WindowConfig(
    window_size=50000,
    step_size=25000,
    min_variants=20  # Only keep windows with ≥20 variants
)
```

### Multi-Contig Datasets

Windows are automatically created per contig:

```python
# Works seamlessly with multi-chromosome VCFs
genomic_ds = GenomicDataset(
    data_source="multi_chromosome.vcf.gz",
    window_config=WindowConfig(window_size=100000)
)

genomic_ds.create_windows()
# Windows are created separately for each chromosome
```

## Best Practices

### 1. Choose Appropriate Window Size

- **Small genomes (bacteria)**: 1-10kb windows
- **Medium genomes (Drosophila)**: 10-100kb windows
- **Large genomes (human)**: 100kb-1Mb windows

```python
# For human genome
window_config = WindowConfig(
    window_size=100000,  # 100kb
    step_size=50000
)

# For bacterial genome
window_config = WindowConfig(
    window_size=5000,    # 5kb
    step_size=2500
)
```

### 2. Set Minimum Variants

Ensure statistical power by requiring minimum variants:

```python
window_config = WindowConfig(
    window_size=100000,
    min_variants=10  # At least 10 SNPs for reliable statistics
)
```

### 3. Balance Overlap and Computation

More overlap = smoother statistics but higher computation:

```python
# High resolution (slow)
WindowConfig(window_size=100000, step_size=10000)  # 90% overlap

# Medium resolution (balanced)
WindowConfig(window_size=100000, step_size=50000)  # 50% overlap

# Low resolution (fast)
WindowConfig(window_size=100000, step_size=100000) # No overlap
```

### 4. Handle Missing Data

Combine windowing with missing data filtering:

```python
from many_stats.core.dataset import CallableSitesConfig

callable_config = CallableSitesConfig(max_missing=0.2)

genomic_ds = GenomicDataset(
    data_source="my_data.vcf.gz",
    callable_config=callable_config,
    window_config=window_config
)

# Filter before windowing
genomic_ds.filter_missing_data()
genomic_ds.create_windows()
```

## Troubleshooting

### Issue: Too Few Windows Created

**Problem**: Only a few windows are created despite large dataset.

**Solutions**:
1. Reduce `min_variants` threshold
2. Increase `window_size`
3. Check that your data has sufficient SNP density

```python
# Before
WindowConfig(window_size=10000, min_variants=50)  # Too restrictive

# After
WindowConfig(window_size=50000, min_variants=10)  # More permissive
```

### Issue: Windows Have Too Few Variants

**Problem**: Many windows are filtered out due to low variant count.

**Solutions**:
1. Increase window size
2. Decrease minimum variant threshold
3. Use less stringent missing data filtering

```python
window_config = WindowConfig(
    window_size=200000,  # Larger windows
    min_variants=5       # Lower threshold
)
```

### Issue: Memory Errors with Large Datasets

**Problem**: Running out of memory with many windows.

**Solutions**:
1. Increase step size (fewer overlapping windows)
2. Process contigs separately
3. Use larger windows

```python
# Memory-efficient settings
window_config = WindowConfig(
    window_size=500000,   # Larger windows
    step_size=500000,     # No overlap
    min_variants=5
)
```

### Issue: Statistics Are NaN

**Problem**: Some windows return NaN values.

**Solutions**:
1. Check for monomorphic windows (no variation)
2. Ensure sufficient variants with `min_variants`
3. Verify data quality

```python
# Filter out low-quality windows
window_config = WindowConfig(
    window_size=100000,
    min_variants=15  # Higher threshold
)
```

## Example Workflow

Complete example for analyzing selection across the genome:

```python
from many_stats.core.dataset import GenomicDataset, WindowConfig, CallableSitesConfig
import matplotlib.pyplot as plt
import pandas as pd

# 1. Configure analysis
window_config = WindowConfig(
    window_size=100000,  # 100kb windows
    step_size=50000,     # 50kb step
    min_variants=10
)

callable_config = CallableSitesConfig(max_missing=0.1)

# 2. Load data
genomic_ds = GenomicDataset(
    data_source="population.vcf.gz",
    callable_config=callable_config,
    window_config=window_config
)

# 3. Filter and create windows
genomic_ds.filter_missing_data()
genomic_ds.create_windows()

# 4. Calculate statistics
window_stats = genomic_ds.calculate_windowed_stats(
    stats=['tajima_d', 'nucleotide_diversity', 'fay_wu_theta']
)

# 5. Save results
genomic_ds.save_results("selection_scan.csv")

# 6. Visualize (assuming single contig)
df = pd.read_csv("selection_scan.csv")
plt.figure(figsize=(12, 4))
plt.plot(df['window_start'], df['tajima_d'])
plt.xlabel('Genomic Position')
plt.ylabel("Tajima's D")
plt.title('Selection Scan Along Chromosome')
plt.savefig('tajima_d_scan.png')
```

## See Also

### many-stats Documentation
- [Callable Sites Guide](callable_sites_guide.md)
- [Statistics Reference](../README.md#available-statistics)
- [CLI Guide](../docs/CLI_GUIDE.md)

### sgkit Resources
- [sgkit Windowing Documentation](https://pystatgen.github.io/sgkit/latest/api.html#windowing) - Details on sgkit's windowing functions
- [sgkit Window API](https://pystatgen.github.io/sgkit/latest/generated/sgkit.window_by_position.html) - window_by_position documentation
- [sgkit Tutorials](https://pystatgen.github.io/sgkit/latest/tutorials.html) - sgkit usage examples

### References
- Wakeley J. (2009) *Coalescent Theory: An Introduction*. Roberts & Company Publishers.
- [Population Genetics Glossary](https://www.nature.com/scitable/topicpage/population-genetics-the-hardy-weinberg-principle-13235724/)

