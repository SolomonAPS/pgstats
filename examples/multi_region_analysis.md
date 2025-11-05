# Multi-Region Analysis Example

This example shows how to use the new multi-region analysis functionality in pgstats.

## Creating a Regions File

Create a BED file with the regions you want to analyze:

```bash
# candidate_genes.bed
chr1    1000000    2000000    gene1
chr1    3000000    4000000    gene2
chr2    500000     1500000    gene3
chrX    100000     800000     gene4
```

## CLI Usage

```bash
# Analyze multiple regions with callable sites masking
pgstats stats genome.vcf.gz \
    --bed callable_sites.bed \
    --regions-file candidate_genes.bed \
    --window-size 10000 \
    --output multi_region_results.csv

# Compare with single region analysis
pgstats stats genome.vcf.gz \
    --bed callable_sites.bed \
    --region chr1:1000000-2000000 \
    --window-size 10000 \
    --output single_region_results.csv
```

## Python API Usage

```python
from pgstats.core.dataset import GenomicDataset, CallableSitesConfig

# Load dataset with callable sites mask
genomic_ds = GenomicDataset(
    data_source="genome.vcf.gz",
    callable_config=CallableSitesConfig(bed_file="callable_sites.bed")
)

# Define regions to analyze
regions = [
    ("chr1", 1000000, 2000000),  # gene1
    ("chr1", 3000000, 4000000),  # gene2
    ("chr2", 500000, 1500000),   # gene3
    ("chrX", 100000, 800000)    # gene4
]

# Analyze all regions
results = genomic_ds.calculate_stats_for_regions(
    regions=regions,
    window_size=10000,
    stats=['tajima_d', 'nucleotide_diversity', 'watterson_theta'],
    min_variants=5
)

# Save results
results.to_csv("multi_region_results.csv", index=False)

# The program automatically prints mean statistics:
# Mean Statistics (multi-region (4 regions)):
# ==================================================
#   tajima_d            : 0.234567 (n=45/50 windows)
#   nucleotide_diversity: 0.001234 (n=48/50 windows)
#   watterson_theta     : 0.000987 (n=47/50 windows)
# ==================================================
# 
# Note: "n=45/50 windows" means 45 out of 50 windows had valid values.
#       Windows with insufficient variants or missing data are excluded.

# Results will have columns:
# region_contig, region_start, region_end, window_start, window_end, 
# n_variants, tajima_d, nucleotide_diversity, watterson_theta
```

## Output Format

The output includes region identifiers for each window:

| region_contig | region_start | region_end | window_start | window_end | n_variants | tajima_d | nucleotide_diversity |
|---------------|--------------|------------|--------------|------------|------------|----------|---------------------|
| chr1          | 1000000      | 2000000    | 1000000      | 1010000    | 15         | 0.23     | 0.0012              |
| chr1          | 1000000      | 2000000    | 1010000      | 1020000    | 12         | -0.45    | 0.0008              |
| chr1          | 3000000      | 4000000    | 3000000      | 3010000    | 18         | 0.67     | 0.0015              |
| ...           | ...          | ...        | ...          | ...        | ...        | ...      | ...                 |

## Key Features

- **Efficient**: Loads VCF once, analyzes multiple regions
- **Flexible**: Works with any BED file format
- **Compatible**: Works with callable sites masking
- **Comprehensive**: Includes region identifiers in output
- **Robust**: Handles errors gracefully (continues if one region fails)
- **Convenient**: Automatically prints mean statistics for quick assessment

## Mean Statistics Output

The program automatically prints mean values for all calculated statistics:

```
Mean Statistics (multi-region (3 regions)):
==================================================
  tajima_d            : 0.250000 (n=6/6 windows)
  nucleotide_diversity: 0.001233 (n=6/6 windows)
  watterson_theta     : 0.001083 (n=6/6 windows)
==================================================
```

**Understanding the output:**
- **Mean value**: Average across all valid windows
- **n=X/Y windows**: X windows had valid values out of Y total windows
- **Missing windows**: Excluded due to insufficient variants or missing data
- **Consistent format**: Same format for all analysis types (multi-region, single-region, genome-wide)

## Use Cases

- **Candidate gene analysis**: Analyze multiple genes of interest
- **Exon analysis**: Analyze all exons of a gene
- **Comparative analysis**: Compare homologous regions across chromosomes
- **Quality control**: Analyze control regions alongside targets
- **Batch processing**: Analyze hundreds of regions efficiently
