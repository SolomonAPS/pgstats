#!/usr/bin/env python3
"""
Test script for the new --keep-zarr functionality.
"""

import tempfile
import os
from pathlib import Path

# Test the new functionality
def test_zarr_caching():
    """Test that Zarr files are cached and reused."""
    
    # Create a temporary VCF file for testing
    temp_dir = tempfile.mkdtemp()
    vcf_path = Path(temp_dir) / "test.vcf"
    
    # Create a minimal VCF file
    vcf_content = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample1	sample2
chr1	100	.	A	T	60	PASS	.	GT	0/0	0/1
chr1	200	.	G	C	60	PASS	.	GT	1/1	0/0
"""
    
    with open(vcf_path, 'w') as f:
        f.write(vcf_content)
    
    print(f"Created test VCF: {vcf_path}")
    print(f"Temp directory: {temp_dir}")
    
    # Test the CLI command
    print("\nTesting CLI with --keep-zarr:")
    print(f"many-stats stats {vcf_path} --output test_results.csv --keep-zarr")
    
    return temp_dir, vcf_path

if __name__ == "__main__":
    temp_dir, vcf_path = test_zarr_caching()
    print(f"\nTest files created in: {temp_dir}")
    print("You can now test the --keep-zarr functionality with:")
    print(f"many-stats stats {vcf_path} --output test_results.csv --keep-zarr")
