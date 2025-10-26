#!/usr/bin/env python3
"""
Test the integrated PyRanges implementation in GenomicDataset.
"""

import pandas as pd
import numpy as np
import xarray as xr
from pathlib import Path
import sys
import os

# Add the package to path
sys.path.insert(0, '/Users/sol/Documents/GitHub/many-stats')

from many_stats.core.dataset import GenomicDataset, CallableSitesConfig

# Get BED file path from previous run
bed_path = "/var/folders/v3/ybcp76dx0bx9gy8dznldx29c0000gp/T/tmp4p5mtw17.bed"

# Create a minimal mock dataset for testing
def create_mock_dataset():
    """Create a minimal xarray Dataset for testing."""
    # Create mock contig data
    contig_id = np.array(['chr1', 'chr2'])
    contig_length = np.array([1000, 1000])
    
    # Create mock variant data (we don't need real variants for L calculation)
    variant_position = np.array([100, 200, 300, 400, 500, 600, 700, 800, 900, 1000])
    variant_contig = np.array([0, 0, 0, 0, 0, 1, 1, 1, 1, 1])  # 5 variants per chr
    
    # Create mock genotype data (shape: variants, samples, ploidy)
    call_genotype = np.random.randint(0, 3, size=(10, 2, 2))  # 10 variants, 2 samples, diploid
    
    # Create xarray Dataset
    ds = xr.Dataset({
        'call_genotype': (['variants', 'samples', 'ploidy'], call_genotype),
        'variant_position': (['variants'], variant_position),
        'variant_contig': (['variants'], variant_contig),
        'contig_id': (['contigs'], contig_id),
        'contig_length': (['contigs'], contig_length),
    })
    
    return ds

# Create mock dataset
mock_dataset = create_mock_dataset()
print("Created mock dataset")

# Create GenomicDataset with BED file
callable_config = CallableSitesConfig(
    bed_file=bed_path,
    bed_format="non_callable"
)

# Create dataset instance (we'll manually set the dataset to avoid VCF loading)
genomic_ds = GenomicDataset.__new__(GenomicDataset)  # Create without calling __init__
genomic_ds.dataset = mock_dataset
genomic_ds.callable_config = callable_config
genomic_ds.callable_sites = None
genomic_ds.non_callable_mask = None

# Manually load the BED file as PyRanges
import pyranges as pr
bed_df = pd.read_csv(bed_path, sep='\t', header=None, names=['chrom', 'start', 'end'])
genomic_ds.bed_gr = pr.PyRanges(bed_df.rename(columns={'chrom': 'Chromosome', 'start': 'Start', 'end': 'End'}))

print("Created GenomicDataset instance with PyRanges BED data")

# Load expected results from corrected implementation
expected_results = pd.read_csv('corrected_implementation_results_fixed.csv')
print(f"\nLoaded {len(expected_results)} expected results from corrected implementation")

# Test integrated PyRanges implementation
print(f"\nTesting integrated PyRanges implementation:")
print("Window\t\tContig\tExpected\tActual\tMatch")
print("-" * 60)

results = []
for _, row in expected_results.iterrows():
    window_start = int(row['window_start'])
    window_stop = int(row['window_stop'])
    contig_name = row['contig_name']
    expected = int(row['expected'])
    
    actual = genomic_ds._calculate_callable_length_in_window(window_start, window_stop, contig_name)
    match = "✓" if actual == expected else "✗"
    results.append({
        'window_start': window_start,
        'window_stop': window_stop,
        'contig_name': contig_name,
        'expected': expected,
        'actual': actual,
        'match': actual == expected
    })
    print(f"[{window_start:3d}, {window_stop:3d})\t{contig_name}\t{expected:3d}bp\t\t{actual:3d}bp\t{match}")

# Summary
total_tests = len(results)
passed_tests = sum(1 for r in results if r['match'])
print(f"\nSummary: {passed_tests}/{total_tests} tests passed")

if passed_tests == total_tests:
    print("✓ All tests passed! Integrated PyRanges implementation works correctly.")
    print("Ready for production use!")
else:
    print("✗ Some tests failed. Need to debug integrated implementation.")

# Save results for comparison
results_df = pd.DataFrame(results)
results_df.to_csv('integrated_pyranges_results.csv', index=False)
print("Results saved to: integrated_pyranges_results.csv")
