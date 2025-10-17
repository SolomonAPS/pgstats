#!/usr/bin/env python3
"""
Test script for VCF loading functionality using vcf2zarr.

This script demonstrates how to load VCF data using the many-stats package.
"""

import os
import sys
from pathlib import Path

# Add the package to the path
sys.path.insert(0, str(Path(__file__).parent))

from many_stats.io.loaders import load_vcf_simple, check_bio2zarr_available


def download_sample_vcf():
    """Download sample VCF data for testing."""
    import urllib.request
    
    sample_vcf_url = "https://raw.githubusercontent.com/sgkit-dev/bio2zarr/main/tests/data/vcf/sample.vcf.gz"
    sample_vcf_path = "sample.vcf.gz"
    
    if not os.path.exists(sample_vcf_path):
        print(f"Downloading sample VCF from {sample_vcf_url}")
        urllib.request.urlretrieve(sample_vcf_url, sample_vcf_path)
        print(f"Downloaded {sample_vcf_path}")
    else:
        print(f"Sample VCF already exists: {sample_vcf_path}")
    
    return sample_vcf_path


def test_vcf_loading():
    """Test VCF loading functionality."""
    print("Testing VCF loading functionality...")
    
    # Check if bio2zarr is available
    if not check_bio2zarr_available():
        print("ERROR: bio2zarr is not available.")
        print("Please install bio2zarr package:")
        print("pip install bio2zarr")
        return False
    
    print("✓ bio2zarr is available")
    
    # Download sample VCF
    try:
        sample_vcf_path = download_sample_vcf()
        print(f"✓ Sample VCF downloaded: {sample_vcf_path}")
    except Exception as e:
        print(f"ERROR: Failed to download sample VCF: {e}")
        return False
    
    # Test loading VCF
    try:
        print("Loading VCF data...")
        dataset = load_vcf_simple(sample_vcf_path)
        
        print("✓ VCF loaded successfully!")
        print(f"Dataset shape: {dataset.sizes}")
        print(f"Available variables: {list(dataset.data_vars.keys())}")
        print(f"Available coordinates: {list(dataset.coords.keys())}")
        
        # Show some basic info about the dataset
        if 'call_genotype' in dataset:
            print(f"Genotype array shape: {dataset.call_genotype.shape}")
            print(f"Genotype array dtype: {dataset.call_genotype.dtype}")
        
        if 'variant_position' in dataset:
            print(f"Number of variants: {len(dataset.variant_position)}")
        
        if 'sample_id' in dataset:
            print(f"Number of samples: {len(dataset.sample_id)}")
            print(f"Sample IDs: {dataset.sample_id.values}")
        
        return True
        
    except Exception as e:
        print(f"ERROR: Failed to load VCF: {e}")
        import traceback
        traceback.print_exc()
        return False


if __name__ == "__main__":
    success = test_vcf_loading()
    if success:
        print("\n🎉 VCF loading test completed successfully!")
    else:
        print("\n❌ VCF loading test failed!")
        sys.exit(1)
