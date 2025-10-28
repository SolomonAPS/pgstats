#!/usr/bin/env python3
"""
Minimal test for tajima_d windowing using synthetic data.

This creates a small in-memory dataset to test windowing without loading VCF files.
"""

import sys
from pathlib import Path
import numpy as np
import xarray as xr
import sgkit as sg

# Add the package to the path
sys.path.insert(0, str(Path(__file__).parent))

from many_stats.stats.sfs_statistics import tajima_d


def create_minimal_dataset(n_variants=100, n_samples=10):
    """Create a minimal sgkit dataset for testing."""
    np.random.seed(42)
    
    # Create synthetic data
    positions = np.sort(np.random.randint(1, 10000, n_variants))
    genotypes = np.random.randint(0, 3, size=(n_variants, n_samples, 2))
    
    # Create dataset
    ds = xr.Dataset({
        "call_genotype": (["variants", "samples", "ploidy"], genotypes),
        "variant_position": (["variants"], positions),
        "variant_contig": (["variants"], np.zeros(n_variants, dtype=int)),
    }, coords={
        "variants": np.arange(n_variants),
        "samples": [f"s{i}" for i in range(n_samples)],
        "ploidy": [0, 1],
        "contigs": np.arange(1),
        "contig_id": (["contigs"], ["chr1"]),
        "contig_length": (["contigs"], [10000])
    })
    
    return ds


def test_tajima_d_windowing():
    """Test tajima_d with windowed dataset."""
    print("🧬 Testing Tajima's D with Windowing (Minimal Test)")
    print("=" * 60)
    
    # Test 1: Single window (genome-wide)
    print("\n1. Testing single genome-wide window")
    ds = create_minimal_dataset(n_variants=100, n_samples=10)
    
    # Create genome-wide window
    ds_windowed = sg.window_by_genome(ds)
    
    print(f"   Variants: {len(ds.variants)}")
    print(f"   Windows created: {len(ds_windowed.windows)}")
    print(f"   Window info:")
    print(f"     - window_start_idx: {ds_windowed.window_start.values}")
    print(f"     - window_stop_idx: {ds_windowed.window_stop.values}")
    print(f"     - Has window_contig: {'window_contig' in ds_windowed.data_vars or 'window_contig' in ds_windowed.coords}")
    
    # Add window indices as data variables (like our code does)
    n_windows = len(ds_windowed.windows)
    positions = ds.variant_position.values
    window_starts = ds_windowed.window_start.values
    window_stops = ds_windowed.window_stop.values
    
    window_start_positions = positions[window_starts]
    window_stop_positions = positions[window_stops - 1]
    
    # Handle window_contig
    if 'window_contig' in ds_windowed.data_vars or 'window_contig' in ds_windowed.coords:
        window_contigs = ds_windowed.window_contig.values
    else:
        window_contigs = np.zeros(n_windows, dtype=int)
    
    ds_windowed = ds_windowed.assign_coords({
        'windows': np.arange(n_windows),
        'window_start': ('windows', window_start_positions),
        'window_stop': ('windows', window_stop_positions),
        'window_contig': ('windows', window_contigs)
    })
    
    ds_windowed['window_start_idx'] = ('windows', window_starts)
    ds_windowed['window_stop_idx'] = ('windows', window_stops)
    
    print(f"\n   Dataset after windowing setup:")
    print(f"     - Dimensions: {dict(ds_windowed.dims)}")
    print(f"     - Data vars: {list(ds_windowed.data_vars)}")
    print(f"     - Coords: {list(ds_windowed.coords)}")
    
    # Calculate tajima_d
    print(f"\n   Calculating tajima_d...")
    try:
        results = tajima_d(ds_windowed)
        
        print(f"   ✅ Success!")
        print(f"   Results:")
        print(f"     - Dimensions: {results.tajima_d.dims}")
        print(f"     - Shape: {results.tajima_d.shape}")
        print(f"     - Values: {results.tajima_d.values}")
        
        # Verify dimensions
        if 'windows' in results.tajima_d.dims:
            print(f"   ✅ CORRECT: Has 'windows' dimension")
        else:
            print(f"   ❌ ERROR: Should have 'windows' dimension")
            return False
            
        if results.tajima_d.shape[0] == n_windows:
            print(f"   ✅ CORRECT: Shape matches number of windows")
        else:
            print(f"   ❌ ERROR: Shape mismatch")
            return False
    except Exception as e:
        print(f"   ❌ ERROR: {e}")
        import traceback
        traceback.print_exc()
        return False
    
    # Test 2: Multiple windows
    print("\n2. Testing multiple windows (1kb windows)")
    ds2 = create_minimal_dataset(n_variants=200, n_samples=10)
    
    # Create multiple windows
    ds_windowed2 = sg.window_by_position(ds2, size=1000, step=1000)
    
    print(f"   Variants: {len(ds2.variants)}")
    print(f"   Windows created: {len(ds_windowed2.windows)}")
    
    # Add window indices
    n_windows2 = len(ds_windowed2.windows)
    positions2 = ds2.variant_position.values
    window_starts2 = ds_windowed2.window_start.values
    window_stops2 = ds_windowed2.window_stop.values
    
    window_start_positions2 = positions2[window_starts2]
    window_stop_positions2 = positions2[window_stops2 - 1]
    
    if 'window_contig' in ds_windowed2.data_vars or 'window_contig' in ds_windowed2.coords:
        window_contigs2 = ds_windowed2.window_contig.values
    else:
        window_contigs2 = np.zeros(n_windows2, dtype=int)
    
    ds_windowed2 = ds_windowed2.assign_coords({
        'windows': np.arange(n_windows2),
        'window_start': ('windows', window_start_positions2),
        'window_stop': ('windows', window_stop_positions2),
        'window_contig': ('windows', window_contigs2)
    })
    
    ds_windowed2['window_start_idx'] = ('windows', window_starts2)
    ds_windowed2['window_stop_idx'] = ('windows', window_stops2)
    
    # Calculate tajima_d
    print(f"\n   Calculating tajima_d for {n_windows2} windows...")
    try:
        results2 = tajima_d(ds_windowed2)
        
        print(f"   ✅ Success!")
        print(f"   Results:")
        print(f"     - Dimensions: {results2.tajima_d.dims}")
        print(f"     - Shape: {results2.tajima_d.shape}")
        print(f"     - First 3 values: {results2.tajima_d.values[:3]}")
        
        # Verify dimensions
        if 'windows' in results2.tajima_d.dims:
            print(f"   ✅ CORRECT: Has 'windows' dimension")
        else:
            print(f"   ❌ ERROR: Should have 'windows' dimension")
            return False
            
        if results2.tajima_d.shape[0] == n_windows2:
            print(f"   ✅ CORRECT: Shape matches number of windows ({n_windows2})")
        else:
            print(f"   ❌ ERROR: Shape mismatch - expected {n_windows2}, got {results2.tajima_d.shape[0]}")
            return False
            
    except Exception as e:
        print(f"   ❌ ERROR: {e}")
        import traceback
        traceback.print_exc()
        return False
    
    print("\n" + "=" * 60)
    print("✅ All tests passed!")
    print("=" * 60)
    return True


if __name__ == "__main__":
    success = test_tajima_d_windowing()
    sys.exit(0 if success else 1)

