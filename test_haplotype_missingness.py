#!/usr/bin/env python3
"""
Test script for the new haplotype missing data handling features.
"""

import numpy as np
import sgkit as sg
import many_stats as ms

def test_missing_data_handling():
    """Test the new missing data handling options."""
    
    print("Testing haplotype missing data handling...")
    
    # Create a simple dataset with some missing data
    ds = sg.simulate_genotype_call_dataset(n_variant=10, n_sample=20, missing_pct=0.2)
    
    print(f"Dataset shape: {ds.sizes}")
    print(f"Missing data percentage: ~20%")
    
    # Test 1: Default behavior (sgkit approach - include missing data)
    print("\n=== Test 1: Default behavior (sgkit approach) ===")
    garud_default = ms.stats.garud_h_statistics(ds)
    print(f"H1: {garud_default.garud_h1.values[0]:.6f}")
    print(f"H12: {garud_default.garud_h12.values[0]:.6f}")
    print(f"H123: {garud_default.garud_h123.values[0]:.6f}")
    print(f"H2/H1: {garud_default.garud_h2_h1.values[0]:.6f}")
    
    # Test 2: Ignore missing data approach
    print("\n=== Test 2: Ignore missing data approach ===")
    garud_ignore = ms.stats.garud_h_statistics(ds, ignore_missing=True)
    print(f"H1: {garud_ignore.garud_h1.values[0]:.6f}")
    print(f"H12: {garud_ignore.garud_h12.values[0]:.6f}")
    print(f"H123: {garud_ignore.garud_h123.values[0]:.6f}")
    print(f"H2/H1: {garud_ignore.garud_h2_h1.values[0]:.6f}")
    
    # Test 3: Missingness threshold (percentage)
    print("\n=== Test 3: Missingness threshold (50% max missing) ===")
    garud_threshold = ms.stats.garud_h_statistics(ds, 
                                                ignore_missing=True,
                                                max_missing=0.5,
                                                missing_is_percentage=True)
    print(f"H1: {garud_threshold.garud_h1.values[0]:.6f}")
    print(f"H12: {garud_threshold.garud_h12.values[0]:.6f}")
    print(f"H123: {garud_threshold.garud_h123.values[0]:.6f}")
    print(f"H2/H1: {garud_threshold.garud_h2_h1.values[0]:.6f}")
    
    # Test 4: Missingness threshold (absolute count)
    print("\n=== Test 4: Missingness threshold (max 3 missing sites) ===")
    garud_absolute = ms.stats.garud_h_statistics(ds, 
                                               ignore_missing=True,
                                               max_missing=3,
                                               missing_is_percentage=False)
    print(f"H1: {garud_absolute.garud_h1.values[0]:.6f}")
    print(f"H12: {garud_absolute.garud_h12.values[0]:.6f}")
    print(f"H123: {garud_absolute.garud_h123.values[0]:.6f}")
    print(f"H2/H1: {garud_absolute.garud_h2_h1.values[0]:.6f}")
    
    print("\n✓ All tests completed successfully!")
    print("\nKey differences:")
    print("- ignore_missing=True: Ignores missing data in haplotype clustering")
    print("- ignore_missing=False: Includes missing data in clustering (default sgkit behavior)")
    print("- max_missing: Filters out entire haplotypes with too much missing data")
    print("- missing_is_percentage: Whether max_missing is percentage (0-1) or absolute count")
    print("\nDefaults:")
    print("- ignore_missing=False (include missing data)")
    print("- max_missing=1.0 (100% missing data allowed)")
    print("- missing_is_percentage=True (max_missing is percentage)")

if __name__ == "__main__":
    test_missing_data_handling()
