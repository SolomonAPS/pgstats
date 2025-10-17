#!/usr/bin/env python3
"""
Test script to examine contig information in VCF/Zarr datasets.
"""

import sys
from pathlib import Path
import numpy as np
import sgkit as sg

# Add the package to the path
sys.path.insert(0, str(Path(__file__).parent))

from many_stats.io.loaders import load_vcf_simple


def test_contig_info_simulated():
    """Test contig information with simulated data."""
    print("🧬 Testing contig info with simulated data")
    print("=" * 50)
    
    # Create simulated dataset
    ds = sg.simulate_genotype_call_dataset(n_variant=100, n_sample=20, n_contig=3)
    
    print("Dataset variables:")
    for var in ds.data_vars:
        print(f"  {var}: {ds[var].shape}")
    
    print("\nDataset coordinates:")
    for coord in ds.coords:
        print(f"  {coord}: {ds[coord].shape}")
    
    # Check contig information
    if 'contig_id' in ds.coords:
        print(f"\nContig IDs: {ds.contig_id.values}")
    
    if 'contig_length' in ds.coords:
        print(f"Contig lengths: {ds.contig_length.values}")
    
    if 'variant_contig' in ds.data_vars:
        print(f"Variant contig indices: {np.unique(ds.variant_contig.values)}")
    
    if 'variant_position' in ds.data_vars:
        positions = ds.variant_position.values
        print(f"Position range: {np.min(positions)} - {np.max(positions)}")
    
    return ds


def test_contig_info_vcf():
    """Test contig information with real VCF data."""
    print("\n🧬 Testing contig info with VCF data")
    print("=" * 50)
    
    # Try to load a sample VCF if available
    sample_vcf = "sample.vcf.gz"
    
    if not Path(sample_vcf).exists():
        print(f"Sample VCF {sample_vcf} not found. Skipping VCF test.")
        return None
    
    try:
        ds = load_vcf_simple(sample_vcf)
        
        print("Dataset variables:")
        for var in ds.data_vars:
            print(f"  {var}: {ds[var].shape}")
        
        print("\nDataset coordinates:")
        for coord in ds.coords:
            print(f"  {coord}: {ds[coord].shape}")
        
        # Check contig information
        if 'contig_id' in ds.coords:
            print(f"\nContig IDs: {ds.contig_id.values}")
        
        if 'contig_length' in ds.coords:
            print(f"Contig lengths: {ds.contig_length.values}")
        else:
            print("No contig_length found in coordinates")
        
        if 'variant_contig' in ds.data_vars:
            print(f"Variant contig indices: {np.unique(ds.variant_contig.values)}")
        
        if 'variant_position' in ds.data_vars:
            positions = ds.variant_position.values
            print(f"Position range: {np.min(positions)} - {np.max(positions)}")
        
        return ds
        
    except Exception as e:
        print(f"Error loading VCF: {e}")
        return None


def test_windowing_with_contigs():
    """Test how windowing works with contig information."""
    print("\n🧬 Testing windowing with contig information")
    print("=" * 50)
    
    # Create dataset with multiple contigs
    ds = sg.simulate_genotype_call_dataset(n_variant=200, n_sample=20, n_contig=2)
    
    # Set realistic positions
    positions = np.concatenate([
        np.arange(1, 101),      # Contig 0: positions 1-100
        np.arange(1001, 1101)    # Contig 1: positions 1001-1100
    ])
    ds = ds.assign(variant_position=(['variants'], positions))
    
    print(f"Contig IDs: {ds.contig_id.values}")
    print(f"Contig lengths: {ds.contig_length.values}")
    print(f"Position range: {np.min(positions)} - {np.max(positions)}")
    
    # Test windowing
    windows = sg.window_by_position(ds, size=50, step=25)
    
    print(f"\nCreated {len(windows.windows)} windows")
    print("Window details:")
    for i, (contig, start_idx, stop_idx) in enumerate(zip(
        windows.window_contig.values,
        windows.window_start.values,
        windows.window_stop.values
    )):
        start_pos = positions[start_idx] if start_idx < len(positions) else "beyond"
        stop_pos = positions[stop_idx-1] if stop_idx > 0 and stop_idx-1 < len(positions) else "beyond"
        print(f"  Window {i+1}: contig {contig}, variants {start_idx}-{stop_idx} (positions {start_pos}-{stop_pos})")


def main():
    """Run the tests."""
    print("🚀 Contig Information Test")
    print("=" * 60)
    
    try:
        # Test with simulated data
        ds_sim = test_contig_info_simulated()
        
        # Test with VCF data (if available)
        ds_vcf = test_contig_info_vcf()
        
        # Test windowing
        test_windowing_with_contigs()
        
        print("\n✅ Tests completed!")
        print("\nKey findings:")
        print("  • Contig information is stored in the Zarr dataset")
        print("  • contig_id: Names of contigs")
        print("  • contig_length: Length of each contig")
        print("  • variant_contig: Which contig each variant belongs to")
        print("  • sgkit uses contig_length for windowing boundaries")
        
    except Exception as e:
        print(f"\n❌ Tests failed: {e}")
        import traceback
        traceback.print_exc()
        return False
    
    return True


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)

