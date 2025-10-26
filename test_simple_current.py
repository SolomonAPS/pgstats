#!/usr/bin/env python3
"""
Simple test of current callable length implementation.
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

print("Created GenomicDataset instance")

# Test a few simple cases manually
print("\nTesting current implementation with simple cases:")

# Test case 1: chr1 window [150, 250) - should have 100bp callable (no BED overlap)
print("\nTest 1: chr1 window [150, 250)")
print("Expected: 100bp callable (no BED regions overlap)")
try:
    result1 = genomic_ds._calculate_callable_length_in_window(150, 250)
    print(f"Actual: {result1}bp")
except Exception as e:
    print(f"Error: {e}")

# Test case 2: chr1 window [350, 450) - should have 50bp callable (BED [300,400) overlaps)
print("\nTest 2: chr1 window [350, 450)")
print("Expected: 50bp callable (BED [300,400) overlaps [350,400) = 50bp masked)")
try:
    result2 = genomic_ds._calculate_callable_length_in_window(350, 450)
    print(f"Actual: {result2}bp")
except Exception as e:
    print(f"Error: {e}")

# Test case 3: chr2 window [100, 200) - should have 100bp callable (no BED overlap)
print("\nTest 3: chr2 window [100, 200)")
print("Expected: 100bp callable (no BED regions overlap)")
try:
    result3 = genomic_ds._calculate_callable_length_in_window(100, 200)
    print(f"Actual: {result3}bp")
except Exception as e:
    print(f"Error: {e}")

# Test case 4: chr2 window [200, 300) - should have 50bp callable (BED [250,350) overlaps)
print("\nTest 4: chr2 window [200, 300)")
print("Expected: 50bp callable (BED [250,350) overlaps [250,300) = 50bp masked)")
try:
    result4 = genomic_ds._calculate_callable_length_in_window(200, 300)
    print(f"Actual: {result4}bp")
except Exception as e:
    print(f"Error: {e}")

print("\nLet's also check what BED regions are being loaded:")
try:
    bed_df = pd.read_csv(bed_path, sep='\t', header=None, names=['chrom', 'start', 'end'])
    print(bed_df)
except Exception as e:
    print(f"Error loading BED: {e}")
