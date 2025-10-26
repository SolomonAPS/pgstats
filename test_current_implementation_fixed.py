#!/usr/bin/env python3
"""
Test the current callable length implementation with proper contig handling.
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

# Load test cases
test_cases_df = pd.read_csv('test_callable_cases.csv')
print(f"Loaded {len(test_cases_df)} test cases")

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

# Test current implementation with contig specification
print("\nTesting current implementation:")
print("Window\t\tContig\tExpected\tActual\tMatch")
print("-" * 60)

# Add contig information to test cases
test_cases_with_contig = [
    # chr1 tests (contig_idx = 0)
    (50, 150, 0, 50),    # Window [50,150): BED [100,200) → overlap [100,150) = 50bp masked → 100-50=50 callable
    (150, 250, 0, 100),  # Window [150,250): No BED overlap → 100bp callable
    (250, 350, 0, 100),  # Window [250,350): No BED overlap → 100bp callable
    (350, 450, 0, 50),   # Window [350,450): BED [300,400) → overlap [350,400) = 50bp masked → 100-50=50 callable
    (450, 550, 0, 100),  # Window [450,550): No BED overlap → 100bp callable
    (550, 650, 0, 50),   # Window [550,650): BED [500,600) → overlap [550,600) = 50bp masked → 100-50=50 callable
    (650, 750, 0, 100),  # Window [650,750): No BED overlap → 100bp callable
    
    # chr2 tests (contig_idx = 1)
    (0, 100, 1, 50),     # Window [0,100): BED [50,150) → overlap [50,100) = 50bp masked → 100-50=50 callable
    (100, 200, 1, 100),  # Window [100,200): No BED overlap → 100bp callable
    (200, 300, 1, 50),   # Window [200,300): BED [250,350) → overlap [250,300) = 50bp masked → 100-50=50 callable
    (300, 400, 1, 50),   # Window [300,400): BED [250,350) → overlap [300,350) = 50bp masked → 100-50=50 callable
    (400, 500, 1, 100),  # Window [400,500): No BED overlap → 100bp callable
    
    # Edge cases
    (95, 105, 0, 5),     # Window [95,105): BED [100,200) → overlap [100,105) = 5bp masked → 10-5=5 callable
    (195, 205, 0, 5),    # Window [195,205): BED [100,200) → overlap [195,200) = 5bp masked → 10-5=5 callable
]

results = []
for window_start, window_stop, contig_idx, expected in test_cases_with_contig:
    try:
        # Call the current method
        actual = genomic_ds._calculate_callable_length_in_window(window_start, window_stop)
        match = "✓" if actual == expected else "✗"
        results.append({
            'window_start': window_start,
            'window_stop': window_stop,
            'contig_idx': contig_idx,
            'expected': expected,
            'actual': actual,
            'match': actual == expected
        })
        contig_name = mock_dataset.contig_id.values[contig_idx]
        print(f"[{window_start:3d}, {window_stop:3d})\t{contig_name}\t{expected:3d}bp\t\t{actual:3d}bp\t{match}")
    except Exception as e:
        print(f"[{window_start:3d}, {window_stop:3d})\tchr?\t{expected:3d}bp\t\tERROR: {e}")
        results.append({
            'window_start': window_start,
            'window_stop': window_stop,
            'contig_idx': contig_idx,
            'expected': expected,
            'actual': None,
            'match': False
        })

# Summary
total_tests = len(results)
passed_tests = sum(1 for r in results if r['match'])
print(f"\nSummary: {passed_tests}/{total_tests} tests passed")

if passed_tests == total_tests:
    print("✓ All tests passed! Current implementation works correctly.")
else:
    print("✗ Some tests failed. Need to debug current implementation first.")

# Save results for comparison
results_df = pd.DataFrame(results)
results_df.to_csv('current_implementation_results_fixed.csv', index=False)
print("Results saved to: current_implementation_results_fixed.csv")
