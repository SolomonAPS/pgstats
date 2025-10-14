# CLI Guide for many-stats

Complete guide to using the many-stats command-line interface.

## Table of Contents

1. [Installation](#installation)
2. [Quick Start](#quick-start)
3. [Commands](#commands)
4. [Common Workflows](#common-workflows)
5. [Advanced Usage](#advanced-usage)
6. [Examples](#examples)

## Installation

After installing many-stats, the `many-stats` command should be available in your environment:

```bash
# Install from source
cd many-stats
pip install -e .

# Verify installation
many-stats --version
```

## Quick Start

### Basic genome-wide analysis

```bash
many-stats stats input.vcf.gz --output results.csv
```

### Windowed analysis

```bash
many-stats stats input.vcf.gz --output results.csv \
    --window-size 100000 \
    --step-size 50000
```

### With callable sites filtering

```bash
many-stats stats input.vcf.gz --output results.csv \
    --bed callable_sites.bed \
    --max-missing 0.1
```

## Commands

### `many-stats stats`

Calculate population genetics statistics from VCF files.

**Required Arguments:**
- `input` - Input VCF file path (.vcf, .vcf.gz, or .zarr)
- `-o, --output OUTPUT` - Output file path

**Statistics Options:**
- `-s, --stats STATS [STATS ...]` - Statistics to calculate
  - Available: `tajima_d`, `fu_li_d`, `fu_li_f`, `theta_pi`, `theta_w`, `theta_h`
  - Aliases: `pi`, `nucleotide_diversity` (theta_pi), `watterson_theta` (theta_w), `fay_wu_theta` (theta_h)
  - Default: `tajima_d theta_pi theta_w theta_h`

**Windowing Options:**
- `-w, --window-size INT` - Window size in base pairs (default: genome-wide)
- `--step-size INT` - Step size for sliding windows (default: same as window-size)
- `--min-variants INT` - Minimum variants per window (default: 5)
- `-r, --region REGION` - Genomic region (format: `chr:start-end` or `chr`)

**Filtering Options:**
- `-b, --bed FILE` - BED file with callable sites
- `--max-missing FLOAT` - Maximum missing data proportion (default: 0.0)

**Output Options:**
- `-f, --format FORMAT` - Output format: csv, tsv, parquet (default: csv)
- `--no-header` - Skip header row in output

**General Options:**
- `-v, --verbose` - Increase verbosity (-v for INFO, -vv for DEBUG)

### `many-stats info`

Display information about a VCF dataset.

**Arguments:**
- `input` - Input VCF file path
- `-b, --bed FILE` - BED file to check callable sites

**Example:**
```bash
many-stats info input.vcf.gz
many-stats info input.vcf.gz --bed callable_sites.bed
```

## Common Workflows

### 1. Genome-Wide Analysis

Calculate genome-wide statistics without windowing:

```bash
many-stats stats population.vcf.gz \
    --output genome_wide_stats.csv \
    --stats tajima_d theta_pi theta_w
```

### 2. Sliding Window Analysis

Analyze the genome in 100kb windows with 50kb overlap:

```bash
many-stats stats population.vcf.gz \
    --output windowed_stats.csv \
    --window-size 100000 \
    --step-size 50000 \
    --min-variants 10
```

### 3. Quality-Filtered Analysis

Apply callable sites masking and missing data filtering:

```bash
many-stats stats population.vcf.gz \
    --output filtered_stats.csv \
    --bed callable_10x.bed \
    --max-missing 0.2 \
    --window-size 50000
```

### 4. Region-Specific Analysis

Analyze a specific genomic region:

```bash
# Whole chromosome
many-stats stats population.vcf.gz \
    --output chr1_stats.csv \
    --region chr1 \
    --window-size 100000

# Specific region
many-stats stats population.vcf.gz \
    --output region_stats.csv \
    --region chr1:1000000-5000000 \
    --window-size 50000
```

### 5. All Statistics

Calculate all available statistics:

```bash
many-stats stats population.vcf.gz \
    --output all_stats.csv \
    --stats tajima_d fu_li_d fu_li_f \
             theta_pi theta_w theta_h \
    --window-size 100000
```

### 6. TSV Output

Output results as tab-separated values:

```bash
many-stats stats population.vcf.gz \
    --output results.tsv \
    --format tsv \
    --window-size 100000
```

## Advanced Usage

### Custom Window Parameters

```bash
# Large windows for low-density data
many-stats stats sparse_data.vcf.gz \
    --output results.csv \
    --window-size 500000 \
    --step-size 500000 \
    --min-variants 3

# Small, overlapping windows for high-resolution
many-stats stats dense_data.vcf.gz \
    --output results.csv \
    --window-size 10000 \
    --step-size 5000 \
    --min-variants 20
```

### Combining Multiple Filters

```bash
# Strict quality filters
many-stats stats population.vcf.gz \
    --output high_quality_stats.csv \
    --bed high_coverage_sites.bed \
    --max-missing 0.05 \
    --min-variants 15 \
    --window-size 100000
```

### Verbose Output

Monitor progress and debug issues:

```bash
# INFO level logging
many-stats stats -v input.vcf.gz --output results.csv

# DEBUG level logging
many-stats stats -vv input.vcf.gz --output results.csv
```

### Pipeline Integration

Use many-stats in a pipeline:

```bash
#!/bin/bash

# Process multiple VCF files
for vcf in data/*.vcf.gz; do
    basename=$(basename $vcf .vcf.gz)
    
    many-stats stats $vcf \
        --output results/${basename}_stats.csv \
        --bed callable_sites.bed \
        --window-size 100000 \
        --stats tajima_d theta_pi
done

# Combine results
cat results/*_stats.csv > combined_stats.csv
```

## Examples

### Example 1: Basic Diversity Scan

```bash
# Calculate diversity across the genome
many-stats stats sample.vcf.gz \
    --output diversity_scan.csv \
    --stats theta_pi theta_w \
    --window-size 100000 \
    --step-size 50000
```

**Output (diversity_scan.csv):**
```
window_start,window_end,n_variants,n_callable_sites,theta_pi,theta_w
0,100000,45,NULL,0.002341,0.002156
50000,150000,52,NULL,0.002789,0.002645
100000,200000,48,NULL,0.002441,0.002298
...
```

### Example 2: Selection Scan

```bash
# Look for signatures of selection
many-stats stats population.vcf.gz \
    --output selection_scan.csv \
    --stats tajima_d fu_li_d fu_li_f \
    --window-size 50000 \
    --step-size 25000 \
    --min-variants 10
```

### Example 3: Quality-Controlled Analysis

```bash
# Generate callable sites BED
bedtools genomecov -ibam aligned.bam -bg | \
    awk '$4 >= 10' | \
    bedtools merge > callable_10x.bed

# Run analysis with quality filters
many-stats stats population.vcf.gz \
    --output qc_stats.csv \
    --bed callable_10x.bed \
    --max-missing 0.15 \
    --window-size 100000 \
    --stats tajima_d nucleotide_diversity
```

### Example 4: Compare Regions

```bash
# Gene region
many-stats stats population.vcf.gz \
    --output gene_region.csv \
    --region chr2:5000000-5100000 \
    --window-size 10000

# Intergenic region
many-stats stats population.vcf.gz \
    --output intergenic_region.csv \
    --region chr2:10000000-10100000 \
    --window-size 10000
```

### Example 5: Multi-Population Analysis

```bash
# Process each population separately
for pop in popA popB popC; do
    many-stats stats ${pop}.vcf.gz \
        --output ${pop}_stats.csv \
        --bed callable_sites.bed \
        --window-size 100000 \
        --stats tajima_d theta_pi
done

# Compare results downstream
```

## Performance Tips

### For Large Files

1. **Use Zarr format** for repeated analyses:
```bash
# Convert once
vcf2zarr input.vcf.gz output.zarr

# Analyze multiple times (faster)
many-stats stats output.zarr --output results1.csv --region chr1
many-stats stats output.zarr --output results2.csv --region chr2
```

2. **Increase window size** to reduce computation:
```bash
many-stats stats large.vcf.gz \
    --output results.csv \
    --window-size 500000  # Larger windows = faster
```

3. **Process chromosomes separately**:
```bash
for chr in {1..22} X Y; do
    many-stats stats genome.vcf.gz \
        --output chr${chr}_stats.csv \
        --region chr${chr} \
        --window-size 100000 &
done
wait
```

### For Low-Coverage Data

```bash
# More permissive settings
many-stats stats lowcov.vcf.gz \
    --output results.csv \
    --max-missing 0.3 \
    --min-variants 3 \
    --window-size 200000
```

## Troubleshooting

### Error: "Input file not found"

**Problem:** File path is incorrect

**Solution:**
```bash
# Use absolute path
many-stats stats /full/path/to/input.vcf.gz --output results.csv

# Or check current directory
ls *.vcf.gz
```

### Error: "BED file not found"

**Problem:** BED file path is incorrect

**Solution:**
```bash
# Verify BED file exists
ls -l callable_sites.bed

# Use absolute path if needed
many-stats stats input.vcf.gz --output results.csv \
    --bed /full/path/to/callable_sites.bed
```

### Warning: "No callable sites found"

**Problem:** Chromosome naming mismatch between VCF and BED

**Solution:**
```bash
# Check VCF chromosomes
many-stats info input.vcf.gz

# Check BED chromosomes
cut -f1 callable_sites.bed | sort -u

# Fix BED file if needed (add/remove "chr" prefix)
sed 's/^/chr/' callable_sites.bed > callable_sites_fixed.bed
```

### Error: "Windows not created"

**Problem:** No windows meet minimum variant threshold

**Solution:**
```bash
# Reduce min-variants threshold
many-stats stats input.vcf.gz --output results.csv \
    --window-size 100000 \
    --min-variants 3  # Lower threshold

# Or increase window size
many-stats stats input.vcf.gz --output results.csv \
    --window-size 500000  # Larger windows
```

### Slow Performance

**Problem:** Analysis takes too long

**Solutions:**
```bash
# 1. Increase window size
--window-size 500000

# 2. Use non-overlapping windows
--window-size 100000 --step-size 100000

# 3. Analyze specific region
--region chr1:1000000-10000000

# 4. Convert to Zarr first
vcf2zarr input.vcf.gz output.zarr
many-stats stats output.zarr --output results.csv
```

## Output Format

### CSV Output (default)

```csv
window_start,window_end,n_variants,n_callable_sites,tajima_d,theta_pi,theta_w
0,100000,45,NULL,0.123456,0.002341,0.002156
50000,150000,52,NULL,-0.456789,0.002789,0.002645
```

### TSV Output

```
window_start	window_end	n_variants	tajima_d
0	100000	45	0.123456
50000	150000	52	-0.456789
```

### Column Descriptions

- `window_start` - Start position of window (bp)
- `window_end` - End position of window (bp)
- `n_variants` - Number of variants in window
- `n_callable_sites` - Number of callable sites (if BED provided)
- `tajima_d` - Tajima's D statistic
- `theta_pi` - Theta-pi (θπ), nucleotide diversity based on pairwise differences
- `theta_w` - Theta-w (θw), Watterson's estimator based on segregating sites
- `theta_h` - Theta-h (θh), Fay and Wu's estimator based on high-frequency derived alleles
- `fu_li_d` - Fu and Li's D* statistic
- `fu_li_f` - Fu and Li's F* statistic

## See Also

### many-stats Documentation
- [Windowing Guide](../examples/windowing_guide.md)
- [Callable Sites Guide](../examples/callable_sites_guide.md)
- [API Documentation](../README.md)
- [GitHub Repository](https://github.com/yourusername/many-stats)

### sgkit Resources
- [sgkit Documentation](https://pystatgen.github.io/sgkit/latest/) - Official sgkit documentation
- [sgkit GitHub](https://github.com/pystatgen/sgkit) - sgkit source code and issues
- [sgkit Windowing](https://pystatgen.github.io/sgkit/latest/api.html#windowing) - sgkit's windowing functions
- [bio2zarr](https://github.com/sgkit-dev/bio2zarr) - VCF to Zarr conversion tool

### File Format References
- [VCF Format Specification](https://samtools.github.io/hts-specs/VCFv4.3.pdf)
- [BED Format](https://genome.ucsc.edu/FAQ/FAQformat.html#format1)
- [Zarr Format](https://zarr.readthedocs.io/) - Cloud-native array storage

