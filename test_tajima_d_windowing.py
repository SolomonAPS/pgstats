#!/usr/bin/env python3
"""
Test script to verify tajima_d works correctly with windowing.

This tests:
1. That tajima_d returns results with ['windows'] dimension
2. That the results make sense for the test data
"""

import sys
from pathlib import Path
import numpy as np
import pandas as pd

# Add the package to the path
sys.path.insert(0, str(Path(__file__).parent))

from many_stats.core.dataset import GenomicDataset, WindowConfig, CallableSitesConfig


def test_tajima_d_windowing():
    """Test tajima_d with windowed dataset."""
    print("🧬 Testing Tajima's D with Windowing")
    print("=" * 60)
    
    # Use the sample VCF files
    vcf_file = "sample_biallelic.vcf.gz"
    
    if not Path(vcf_file).exists():
        print(f"❌ Error: {vcf_file} not found")
        print("   Please ensure you're in the repository root")
        return False
    
    print(f"\n1. Loading VCF: {vcf_file}")
    ds = GenomicDataset(
        vcf_file,
        enable_profiling=False
    )
    
    print(f"   Loaded {len(ds.dataset.variants):,} variants, {len(ds.dataset.samples)} samples")
    
    # Test 1: Genome-wide (single window via sg.window_by_genome)
    print("\n2. Testing genome-wide windowing (single window)")
    ds_genome_wide = GenomicDataset(
        vcf_file,
        enable_profiling=False
    )
    ds_genome_wide.create_windows()
    
    print(f"   Created {len(ds_genome_wide.windowed_dataset.windows)} window(s)")
    print(f"   Window info:")
    if len(ds_genome_wide.windowed_dataset.windows) > 0:
        print(f"     - Dimensions: {dict(ds_genome_wide.windowed_dataset.dims)}")
        print(f"     - Has window_start_idx: {'window_start_idx' in ds_genome_wide.windowed_dataset.data_vars}")
        print(f"     - Has window_stop_idx: {'window_stop_idx' in ds_genome_wide.windowed_dataset.data_vars}")
    
    # Calculate tajima_d
    results_genome = ds_genome_wide.calculate_windowed_stats(['tajima_d'])
    
    print(f"\n   Tajima's D results:")
    print(f"     - Dimensions: {results_genome.tajima_d.dims}")
    print(f"     - Shape: {results_genome.tajima_d.shape}")
    print(f"     - Values: {results_genome.tajima_d.values}")
    
    # Verify dimensions
    if 'windows' in results_genome.tajima_d.dims:
        print("   ✅ CORRECT: tajima_d has 'windows' dimension")
    else:
        print("   ❌ ERROR: tajima_d should have 'windows' dimension")
        return False
    
    if results_genome.tajima_d.shape[0] == len(results_genome.windows):
        print("   ✅ CORRECT: Shape matches number of windows")
    else:
        print(f"   ❌ ERROR: Shape mismatch - expected {len(results_genome.windows)}, got {results_genome.tajima_d.shape[0]}")
        return False
    
    # Test 2: Multiple windows
    print("\n3. Testing multiple windows (5kb windows)")
    ds_windowed = GenomicDataset(
        vcf_file,
        window_config=WindowConfig(window_size=5000, step_size=5000),
        enable_profiling=False
    )
    ds_windowed.create_windows()
    
    print(f"   Created {len(ds_windowed.windowed_dataset.windows)} windows")
    
    # Calculate tajima_d
    results_windowed = ds_windowed.calculate_windowed_stats(['tajima_d'])
    
    print(f"\n   Tajima's D results:")
    print(f"     - Dimensions: {results_windowed.tajima_d.dims}")
    print(f"     - Shape: {results_windowed.tajima_d.shape}")
    print(f"     - First 5 values: {results_windowed.tajima_d.values[:5]}")
    print(f"     - Number of NaN values: {np.sum(np.isnan(results_windowed.tajima_d.values))}")
    
    # Verify dimensions
    if 'windows' in results_windowed.tajima_d.dims:
        print("   ✅ CORRECT: tajima_d has 'windows' dimension")
    else:
        print("   ❌ ERROR: tajima_d should have 'windows' dimension")
        return False
    
    if results_windowed.tajima_d.shape[0] == len(results_windowed.windows):
        print("   ✅ CORRECT: Shape matches number of windows")
    else:
        print(f"   ❌ ERROR: Shape mismatch - expected {len(results_windowed.windows)}, got {results_windowed.tajima_d.shape[0]}")
        return False
    
    # Test 3: Multiple statistics
    print("\n4. Testing multiple statistics together")
    results_multi = ds_windowed.calculate_windowed_stats(['tajima_d', 'theta_pi', 'theta_w'])
    
    print(f"   Results:")
    print(f"     - tajima_d dimensions: {results_multi.tajima_d.dims}")
    print(f"     - theta_pi dimensions: {results_multi.theta_pi.dims}")
    print(f"     - theta_w dimensions: {results_multi.theta_w.dims}")
    
    all_have_windows = all(
        'windows' in results_multi[stat].dims 
        for stat in ['tajima_d', 'theta_pi', 'theta_w']
        if stat in results_multi.data_vars
    )
    
    if all_have_windows:
        print("   ✅ CORRECT: All statistics have 'windows' dimension")
    else:
        print("   ❌ ERROR: Some statistics don't have 'windows' dimension")
        return False
    
    # Test 4: Convert to DataFrame
    print("\n5. Testing DataFrame conversion")
    try:
        df = results_multi.to_dataframe()
        print(f"   Data frame shape: {df.shape}")
        print(f"   Columns: {list(df.columns)}")
        print(f"\n   First few rows:")
        print(df[['window_start', 'window_stop', 'tajima_d', 'theta_pi']].head())
        
        print("   ✅ CORRECT: Successfully converted to DataFrame")
    except Exception as e:
        print(f"   ❌ ERROR: Failed to convert to DataFrame: {e}")
        return False
    
    print("\n" + "=" * 60)
    print("✅ All tests passed!")
    print("=" * 60)
    return True


if __name__ == "__main__":
    success = test_tajima_d_windowing()
    sys.exit(0 if success else 1)


