#!/usr/bin/env python3
"""
Simple test to verify the code changes were applied correctly.
This doesn't require any dependencies - just checks the source code.
"""

import sys
from pathlib import Path

def test_use_dask_parameter():
    """Check that use_dask parameter was added to load_vcf."""
    print("Checking use_dask parameter in loaders.py...")
    
    loaders_path = Path("pgstats/io/loaders.py")
    if not loaders_path.exists():
        print(f"✗ File not found: {loaders_path}")
        return False
    
    content = loaders_path.read_text()
    
    # Check for use_dask parameter
    if "use_dask: bool = False" not in content:
        print("✗ use_dask parameter not found or has wrong default")
        return False
    
    # Check for chunks logic
    if "chunks = 'auto' if use_dask else None" not in content:
        print("✗ chunks logic not found")
        return False
    
    print("✓ use_dask parameter added correctly")
    return True


def test_bio2zarr_two_step():
    """Check that bio2zarr two-step workflow is used."""
    print("\nChecking bio2zarr two-step workflow in loaders.py...")
    
    loaders_path = Path("pgstats/io/loaders.py")
    content = loaders_path.read_text()
    
    # Check for explode
    if "v2z.explode(" not in content:
        print("✗ v2z.explode() not found")
        return False
    
    # Check for encode
    if "v2z.encode(" not in content:
        print("✗ v2z.encode() not found")
        return False
    
    # Check for ICF cleanup
    if "shutil.rmtree(icf_path)" not in content:
        print("✗ ICF cleanup not found")
        return False
    
    # Make sure old convert() is NOT used
    if "v2z.convert(" in content:
        print("✗ Old v2z.convert() still present")
        return False
    
    print("✓ bio2zarr two-step workflow implemented correctly")
    return True


def test_empty_window_fix():
    """Check that empty window handling is fixed."""
    print("\nChecking empty window fix in dataset.py...")
    
    dataset_path = Path("pgstats/core/dataset.py")
    if not dataset_path.exists():
        print(f"✗ File not found: {dataset_path}")
        return False
    
    content = dataset_path.read_text()
    
    # Check for np.where fix
    if "np.where(" not in content:
        print("✗ np.where() not found")
        return False
    
    # Check for the specific fix
    if "window_stops > window_starts" not in content:
        print("✗ Empty window check not found")
        return False
    
    if "positions[window_starts]       # Empty window" not in content:
        print("✗ Empty window handling not found")
        return False
    
    print("✓ Empty window handling fixed correctly")
    return True


def main():
    """Run all tests."""
    print("="*60)
    print("Testing Code Changes (No Dependencies Required)")
    print("="*60)
    print()
    
    results = []
    results.append(test_use_dask_parameter())
    results.append(test_bio2zarr_two_step())
    results.append(test_empty_window_fix())
    
    print("\n" + "="*60)
    if all(results):
        print("✓ All code changes verified!")
        print("="*60)
        print("\nThe 3 critical fixes have been applied correctly:")
        print("  1. use_dask parameter (default: False)")
        print("  2. bio2zarr two-step workflow (explode → encode)")
        print("  3. Empty window handling (np.where check)")
        print("\nReady to test on the cluster!")
        return 0
    else:
        print("✗ Some code changes missing or incorrect")
        print("="*60)
        return 1


if __name__ == "__main__":
    sys.exit(main())

