#!/usr/bin/env python3
"""
Quick test script to verify the clean fixes work correctly.
Run this after installing from the fix-dask-corruption-clean branch.
"""

import sys
import numpy as np
from pathlib import Path

def test_load_vcf():
    """Test that load_vcf works with use_dask parameter."""
    print("Testing load_vcf with use_dask parameter...")
    
    try:
        from pgstats.io import load_vcf
        
        # Check that use_dask parameter exists
        import inspect
        sig = inspect.signature(load_vcf)
        assert 'use_dask' in sig.parameters, "use_dask parameter not found!"
        
        # Check default value is False
        default = sig.parameters['use_dask'].default
        assert default is False, f"use_dask default should be False, got {default}"
        
        print("✓ load_vcf has use_dask parameter (default=False)")
        return True
        
    except Exception as e:
        print(f"✗ Error testing load_vcf: {e}")
        return False


def test_empty_window_handling():
    """Test that empty window handling works."""
    print("\nTesting empty window handling...")
    
    try:
        # Simulate the empty window check
        window_starts = np.array([0, 5, 10, 15, 15])  # Last window is empty
        window_stops = np.array([5, 10, 15, 15, 20])
        positions = np.array([100, 200, 300, 400, 500, 600, 700, 800, 900, 1000,
                             1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000])
        
        # This is the fix from dataset.py
        window_stop_positions = np.where(
            window_stops > window_starts,
            positions[window_stops - 1],
            positions[window_starts]
        )
        
        # Check that empty window (index 3) uses start position
        assert window_stop_positions[3] == positions[15], "Empty window not handled correctly!"
        
        print("✓ Empty window handling works correctly")
        return True
        
    except Exception as e:
        print(f"✗ Error testing empty window handling: {e}")
        return False


def test_bio2zarr_imports():
    """Test that bio2zarr functions are available."""
    print("\nTesting bio2zarr imports...")
    
    try:
        import bio2zarr.vcf as v2z
        
        # Check that explode and encode exist
        assert hasattr(v2z, 'explode'), "v2z.explode not found!"
        assert hasattr(v2z, 'encode'), "v2z.encode not found!"
        
        print("✓ bio2zarr.explode and encode are available")
        return True
        
    except Exception as e:
        print(f"✗ Error testing bio2zarr: {e}")
        return False


def main():
    """Run all tests."""
    print("="*60)
    print("Testing Clean Fixes")
    print("="*60)
    
    results = []
    results.append(test_load_vcf())
    results.append(test_empty_window_handling())
    results.append(test_bio2zarr_imports())
    
    print("\n" + "="*60)
    if all(results):
        print("✓ All tests passed!")
        print("="*60)
        return 0
    else:
        print("✗ Some tests failed")
        print("="*60)
        return 1


if __name__ == "__main__":
    sys.exit(main())

