#!/usr/bin/env python3
"""
Test script to understand how sgkit handles windowing when start/end positions
extend beyond the actual variant positions.
"""

import sys
from pathlib import Path
import numpy as np
import sgkit as sg

# Add the package to the path
sys.path.insert(0, str(Path(__file__).parent))


def test_sgkit_windowing_behavior():
    """Test how sgkit handles windowing with positions beyond variant range."""
    print("🧬 Testing sgkit windowing behavior")
    print("=" * 50)
    
    # Create a simple dataset with variants at positions 100, 200, 300, 400, 500
    positions = np.array([100, 200, 300, 400, 500])
    n_variants = len(positions)
    n_samples = 10
    
    # Create genotype data
    genotypes = np.random.randint(0, 3, size=(n_variants, n_samples, 2))
    
    # Create dataset
    ds = sg.create_genotype_call_dataset(
        call_genotype=genotypes,
        variant_position=positions,
        variant_contig=np.zeros(n_variants, dtype=int),
        sample_id=[f"sample_{i}" for i in range(n_samples)],
        contig_id=["chr1"],
        contig_length=[1000]  # Set contig length to 1000
    )
    
    print(f"Created dataset with {n_variants} variants at positions: {positions}")
    print(f"Contig length: {ds.contig_length.values[0]}")
    
    # Test 1: Windows that fit within variant range
    print("\n--- Test 1: Windows within variant range ---")
    windows1 = sg.window_by_position(ds, size=150, step=100)
    print(f"Window size: 150bp, step: 100bp")
    print(f"Number of windows: {len(windows1.windows)}")
    
    for i, (start_idx, stop_idx) in enumerate(zip(windows1.window_start.values, windows1.window_stop.values)):
        start_pos = positions[start_idx] if start_idx < len(positions) else "beyond"
        stop_pos = positions[stop_idx-1] if stop_idx > 0 and stop_idx-1 < len(positions) else "beyond"
        print(f"  Window {i+1}: variants {start_idx}-{stop_idx} (positions {start_pos}-{stop_pos})")
    
    # Test 2: Windows that extend beyond variant range
    print("\n--- Test 2: Windows extending beyond variant range ---")
    windows2 = sg.window_by_position(ds, size=200, step=100)
    print(f"Window size: 200bp, step: 100bp")
    print(f"Number of windows: {len(windows2.windows)}")
    
    for i, (start_idx, stop_idx) in enumerate(zip(windows2.window_start.values, windows2.window_stop.values)):
        start_pos = positions[start_idx] if start_idx < len(positions) else "beyond"
        stop_pos = positions[stop_idx-1] if stop_idx > 0 and stop_idx-1 < len(positions) else "beyond"
        print(f"  Window {i+1}: variants {start_idx}-{stop_idx} (positions {start_pos}-{stop_pos})")
    
    # Test 3: Very large windows
    print("\n--- Test 3: Very large windows ---")
    windows3 = sg.window_by_position(ds, size=1000, step=500)
    print(f"Window size: 1000bp, step: 500bp")
    print(f"Number of windows: {len(windows3.windows)}")
    
    for i, (start_idx, stop_idx) in enumerate(zip(windows3.window_start.values, windows3.window_stop.values)):
        start_pos = positions[start_idx] if start_idx < len(positions) else "beyond"
        stop_pos = positions[stop_idx-1] if stop_idx > 0 and stop_idx-1 < len(positions) else "beyond"
        print(f"  Window {i+1}: variants {start_idx}-{stop_idx} (positions {start_pos}-{stop_pos})")
    
    # Test 4: Check what happens with offset
    print("\n--- Test 4: Windows with offset ---")
    windows4 = sg.window_by_position(ds, size=150, step=100, offset=50)
    print(f"Window size: 150bp, step: 100bp, offset: 50bp")
    print(f"Number of windows: {len(windows4.windows)}")
    
    for i, (start_idx, stop_idx) in enumerate(zip(windows4.window_start.values, windows4.window_stop.values)):
        start_pos = positions[start_idx] if start_idx < len(positions) else "beyond"
        stop_pos = positions[stop_idx-1] if stop_idx > 0 and stop_idx-1 < len(positions) else "beyond"
        print(f"  Window {i+1}: variants {start_idx}-{stop_idx} (positions {start_pos}-{stop_pos})")
    
    return ds, windows1, windows2, windows3, windows4


def test_our_region_filtering():
    """Test our region filtering behavior."""
    print("\n🧬 Testing our region filtering behavior")
    print("=" * 50)
    
    # Create dataset
    positions = np.array([100, 200, 300, 400, 500])
    n_variants = len(positions)
    n_samples = 10
    
    genotypes = np.random.randint(0, 3, size=(n_variants, n_samples, 2))
    
    ds = sg.create_genotype_call_dataset(
        call_genotype=genotypes,
        variant_position=positions,
        variant_contig=np.zeros(n_variants, dtype=int),
        sample_id=[f"sample_{i}" for i in range(n_samples)],
        contig_id=["chr1"],
        contig_length=[1000]
    )
    
    print(f"Original dataset: {n_variants} variants at positions: {positions}")
    
    # Test filtering to region that extends beyond variants
    print("\n--- Filtering to region 50-600 (extends beyond variants) ---")
    region_mask = (positions >= 50) & (positions <= 600)
    filtered_positions = positions[region_mask]
    print(f"Variants in region: {np.sum(region_mask)} at positions: {filtered_positions}")
    
    # Test filtering to region within variants
    print("\n--- Filtering to region 150-450 (within variants) ---")
    region_mask = (positions >= 150) & (positions <= 450)
    filtered_positions = positions[region_mask]
    print(f"Variants in region: {np.sum(region_mask)} at positions: {filtered_positions}")
    
    # Test filtering to region with no variants
    print("\n--- Filtering to region 50-99 (no variants) ---")
    region_mask = (positions >= 50) & (positions <= 99)
    filtered_positions = positions[region_mask]
    print(f"Variants in region: {np.sum(region_mask)} at positions: {filtered_positions}")


def main():
    """Run the tests."""
    print("🚀 sgkit Windowing Behavior Test")
    print("=" * 60)
    
    try:
        # Test sgkit's behavior
        ds, w1, w2, w3, w4 = test_sgkit_windowing_behavior()
        
        # Test our filtering
        test_our_region_filtering()
        
        print("\n✅ Tests completed!")
        print("\nKey observations:")
        print("  • sgkit creates windows based on genomic positions, not variant indices")
        print("  • Windows can extend beyond the actual variant positions")
        print("  • Empty windows (no variants) are still created")
        print("  • Our region filtering only includes variants within the specified range")
        
    except Exception as e:
        print(f"\n❌ Tests failed: {e}")
        import traceback
        traceback.print_exc()
        return False
    
    return True


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)

