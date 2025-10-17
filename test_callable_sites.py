#!/usr/bin/env python3
"""
Test script demonstrating the new callable sites approach.

This script shows how non-callable sites are masked by setting them to -1
in the genotype data, allowing sgkit's built-in missing data handling to
automatically account for callable sites.
"""

import sys
from pathlib import Path
import numpy as np
import pandas as pd

# Add the package to the path
sys.path.insert(0, str(Path(__file__).parent))

from many_stats.core.dataset import GenomicDataset, WindowConfig, CallableSitesConfig
import sgkit as sg


def create_test_bed_file():
    """Create a test BED file for demonstration."""
    bed_data = {
        'chrom': ['chr1', 'chr1', 'chr1'],
        'start': [1, 201, 401],  # 1-based start
        'end': [200, 400, 600]    # 0-based end (not inclusive)
    }
    
    bed_df = pd.DataFrame(bed_data)
    bed_df.to_csv('test_callable.bed', sep='\t', header=False, index=False)
    
    print("Created test BED file: test_callable.bed")
    print("Callable regions:")
    for _, row in bed_df.iterrows():
        print(f"  {row['chrom']}:{row['start']}-{row['end']}")
    
    return 'test_callable.bed'


def test_callable_sites():
    """Test the callable sites functionality."""
    print("🧬 Testing Callable Sites Implementation")
    print("=" * 50)
    
    # Create simulated dataset
    ds = sg.simulate_genotype_call_dataset(n_variant=1000, n_sample=20, missing_pct=0.05)
    
    # Set some realistic positions
    positions = np.arange(1, 1001)  # Positions 1-1000
    ds = ds.assign(variant_position=(['variants'], positions))
    
    print(f"Original dataset: {len(ds.variants)} variants, {len(ds.samples)} samples")
    
    # Check original missing data
    genotypes = ds.call_genotype.values
    original_missing = np.sum(genotypes == -1)
    print(f"Original missing genotypes: {original_missing:,}")
    
    # Create test BED file
    bed_file = create_test_bed_file()
    
    # Configure callable sites
    callable_config = CallableSitesConfig(
        bed_file=bed_file,
        max_missing=0.2
    )
    
    window_config = WindowConfig(
        window_size=100,
        step_size=50,
        min_variants=2
    )
    
    # Initialize dataset with callable sites
    genomic_ds = GenomicDataset(
        data_source=ds,
        callable_config=callable_config,
        window_config=window_config
    )
    
    print(f"\nAfter applying callable sites mask:")
    print(f"Callable sites: {genomic_ds.callable_sites:,}")
    
    # Check new missing data
    new_genotypes = genomic_ds.dataset.call_genotype.values
    new_missing = np.sum(new_genotypes == -1)
    print(f"Missing genotypes after masking: {new_missing:,}")
    print(f"Additional missing genotypes: {new_missing - original_missing:,}")
    
    # Test windowing
    print(f"\nCreating windows...")
    genomic_ds.create_windows()
    
    # Calculate statistics
    print(f"Calculating windowed statistics...")
    window_stats_ds = genomic_ds.calculate_windowed_stats(
        stats=['tajima_d', 'pi']
    )
    
    print(f"Calculated statistics for {len(window_stats_ds.windows)} windows")
    
    # Show some results
    tajima_d_values = window_stats_ds['tajima_d'].values
    pi_values = window_stats_ds['pi'].values
    
    print(f"\nWindowed statistics:")
    print(f"Tajima's D: {np.nanmin(tajima_d_values):.3f} to {np.nanmax(tajima_d_values):.3f}")
    print(f"π: {np.nanmin(pi_values):.3f} to {np.nanmax(pi_values):.3f}")
    
    # Test genome-wide statistics
    print(f"\nCalculating genome-wide statistics...")
    genome_stats = genomic_ds.calculate_genome_wide_stats(
        stats=['tajima_d', 'pi']
    )
    
    print(f"Genome-wide statistics:")
    for stat, value in genome_stats.items():
        print(f"  {stat}: {value:.4f}")
    
    # Clean up
    Path(bed_file).unlink()
    print(f"\nCleaned up test file: {bed_file}")
    
    return genomic_ds


def main():
    """Run the test."""
    print("🚀 Callable Sites Test")
    print("=" * 60)
    
    try:
        genomic_ds = test_callable_sites()
        
        print("\n✅ Test completed successfully!")
        print("\nKey features demonstrated:")
        print("  • BED file loading with 1-based start, 0-based end")
        print("  • Non-callable sites masked as -1 (missing)")
        print("  • sgkit's missing data handling automatically accounts for callable sites")
        print("  • Windowed and genome-wide statistics work seamlessly")
        print("  • No need for separate callable sites tracking in statistics")
        
    except Exception as e:
        print(f"\n❌ Test failed: {e}")
        import traceback
        traceback.print_exc()
        return False
    
    return True


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)

