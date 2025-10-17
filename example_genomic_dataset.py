#!/usr/bin/env python3
"""
Example script demonstrating the GenomicDataset framework.

This script shows how to use the new core framework for running
windowed and non-windowed population genetics statistics.
"""

import sys
from pathlib import Path
import numpy as np

# Add the package to the path
sys.path.insert(0, str(Path(__file__).parent))

from many_stats.core.dataset import GenomicDataset, WindowConfig, CallableSitesConfig
import sgkit as sg


def example_with_simulated_data():
    """Example using simulated data."""
    print("🧬 Example: Using simulated data")
    print("=" * 50)
    
    # Create simulated dataset
    ds = sg.simulate_genotype_call_dataset(n_variant=1000, n_sample=50, missing_pct=0.1)
    
    # Configure analysis
    window_config = WindowConfig(
        window_size=200,  # Smaller windows for demo
        step_size=100,
        min_variants=3
    )
    
    callable_config = CallableSitesConfig(
        max_missing=0.2  # Allow up to 20% missing data
    )
    
    # Initialize dataset
    genomic_ds = GenomicDataset(
        data_source=ds,
        callable_config=callable_config,
        window_config=window_config
    )
    
    # Get summary
    summary = genomic_ds.get_summary()
    print(f"Dataset summary: {summary}")
    
    # Filter missing data
    genomic_ds.filter_missing_data()
    
    # Create windows
    genomic_ds.create_windows()
    
    # Calculate windowed statistics
    print("\nCalculating windowed statistics...")
    window_stats_ds = genomic_ds.calculate_windowed_stats(
        stats=['tajima_d', 'pi', 'theta_w']
    )
    
    print(f"Calculated statistics for {len(window_stats_ds.windows)} windows")
    
    # Extract statistics for display
    tajima_d_values = window_stats_ds['tajima_d'].values
    pi_values = window_stats_ds['pi'].values
    
    print(f"Tajima's D range: {np.nanmin(tajima_d_values):.3f} to {np.nanmax(tajima_d_values):.3f}")
    print(f"π range: {np.nanmin(pi_values):.3f} to {np.nanmax(pi_values):.3f}")
    
    # Calculate genome-wide statistics
    print("\nCalculating genome-wide statistics...")
    genome_stats = genomic_ds.calculate_genome_wide_stats(
        stats=['tajima_d', 'pi', 'theta_w']
    )
    
    print("Genome-wide statistics:")
    for stat, value in genome_stats.items():
        print(f"  {stat}: {value:.4f}")
    
    return genomic_ds


def example_with_vcf_file():
    """Example using VCF file (if available)."""
    print("\n🧬 Example: Using VCF file")
    print("=" * 50)
    
    # This would work with a real VCF file
    vcf_path = "sample.vcf.gz"  # Replace with actual VCF file
    
    if not Path(vcf_path).exists():
        print(f"VCF file {vcf_path} not found. Skipping VCF example.")
        return None
    
    # Configure analysis with BED mask
    window_config = WindowConfig(
        window_size=1000,
        step_size=500,
        min_variants=5
    )
    
    callable_config = CallableSitesConfig(
        bed_file="callable_sites.bed",  # Replace with actual BED file
        max_missing=0.1
    )
    
    try:
        # Initialize dataset
        genomic_ds = GenomicDataset(
            data_source=vcf_path,
            callable_config=callable_config,
            window_config=window_config
        )
        
        # Process data
        genomic_ds.filter_missing_data()
        genomic_ds.create_windows()
        
        # Calculate statistics
        window_stats_ds = genomic_ds.calculate_windowed_stats(
            stats=['tajima_d', 'pi', 'theta_w', 'fu_li_d']
        )
        
        # Save results
        genomic_ds.save_results("windowed_stats.csv")
        
        return genomic_ds
        
    except Exception as e:
        print(f"Error processing VCF file: {e}")
        return None


def main():
    """Run examples."""
    print("🚀 GenomicDataset Framework Examples")
    print("=" * 60)
    
    # Run simulated data example
    genomic_ds = example_with_simulated_data()
    
    # Run VCF example (if file exists)
    example_with_vcf_file()
    
    print("\n✅ Examples completed!")
    print("\nKey features demonstrated:")
    print("  • Loading VCF data with bio2zarr")
    print("  • Configurable windowing (size, step, min variants)")
    print("  • Callable sites handling with BED files")
    print("  • Missing data filtering")
    print("  • Windowed and genome-wide statistics")
    print("  • Multiple output formats")
    
    print("\nNext steps:")
    print("  • Install bio2zarr: pip install bio2zarr")
    print("  • Try with your own VCF files")
    print("  • Add custom statistics to the framework")


if __name__ == "__main__":
    main()
