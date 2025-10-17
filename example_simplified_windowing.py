#!/usr/bin/env python3
"""
Example demonstrating the simplified windowing interface.

This script shows the three main analysis modes:
1. Genome-wide analysis (default)
2. Windowed analysis
3. Region-specific analysis
"""

import sys
from pathlib import Path
import numpy as np

# Add the package to the path
sys.path.insert(0, str(Path(__file__).parent))

from many_stats.core.dataset import GenomicDataset, WindowConfig, CallableSitesConfig
import sgkit as sg


def example_genome_wide():
    """Example: Genome-wide analysis (default when no windows specified)."""
    print("🧬 Example 1: Genome-wide Analysis")
    print("=" * 50)
    
    # Create simulated dataset
    ds = sg.simulate_genotype_call_dataset(n_variant=500, n_sample=30, missing_pct=0.1)
    
    # No window configuration = genome-wide analysis
    window_config = WindowConfig()  # All defaults
    
    genomic_ds = GenomicDataset(
        data_source=ds,
        window_config=window_config
    )
    
    # Create windows (will be genome-wide)
    genomic_ds.create_windows()
    
    # Calculate statistics
    stats_ds = genomic_ds.calculate_windowed_stats(stats=['tajima_d', 'pi'])
    
    summary = genomic_ds.get_summary()
    print(f"Analysis type: {summary['analysis_type']}")
    print(f"Number of windows: {summary['n_windows']}")
    print(f"Tajima's D: {stats_ds['tajima_d'].values[0]:.4f}")
    print(f"π: {stats_ds['pi'].values[0]:.4f}")
    
    return genomic_ds


def example_windowed():
    """Example: Windowed analysis."""
    print("\n🧬 Example 2: Windowed Analysis")
    print("=" * 50)
    
    # Create simulated dataset
    ds = sg.simulate_genotype_call_dataset(n_variant=1000, n_sample=30, missing_pct=0.1)
    
    # Configure windowed analysis
    window_config = WindowConfig(
        window_size=200,  # 200bp windows
        step_size=100,   # 100bp step (50% overlap)
        min_variants=2   # Minimum 2 variants per window
    )
    
    genomic_ds = GenomicDataset(
        data_source=ds,
        window_config=window_config
    )
    
    # Create windows
    genomic_ds.create_windows()
    
    # Calculate statistics
    stats_ds = genomic_ds.calculate_windowed_stats(stats=['tajima_d', 'pi'])
    
    summary = genomic_ds.get_summary()
    print(f"Analysis type: {summary['analysis_type']}")
    print(f"Number of windows: {summary['n_windows']}")
    print(f"Window size: {summary['window_size']}bp")
    print(f"Step size: {summary['step_size']}bp")
    
    # Show some statistics
    tajima_d_values = stats_ds['tajima_d'].values
    pi_values = stats_ds['pi'].values
    
    print(f"Tajima's D range: {np.nanmin(tajima_d_values):.3f} to {np.nanmax(tajima_d_values):.3f}")
    print(f"π range: {np.nanmin(pi_values):.3f} to {np.nanmax(pi_values):.3f}")
    
    return genomic_ds


def example_region_specific():
    """Example: Region-specific analysis."""
    print("\n🧬 Example 3: Region-specific Analysis")
    print("=" * 50)
    
    # Create simulated dataset with realistic positions
    ds = sg.simulate_genotype_call_dataset(n_variant=2000, n_sample=30, missing_pct=0.1)
    
    # Set positions to simulate a chromosome
    positions = np.arange(1, 2001)  # Positions 1-2000
    ds = ds.assign(variant_position=(['variants'], positions))
    
    # Configure region-specific analysis
    window_config = WindowConfig(
        start=500,       # Start at position 500
        end=1500,        # End at position 1500
        window_size=100, # 100bp windows within the region
        step_size=50,   # 50bp step
        min_variants=2
    )
    
    genomic_ds = GenomicDataset(
        data_source=ds,
        window_config=window_config
    )
    
    # Create windows (will be filtered to region first)
    genomic_ds.create_windows()
    
    # Calculate statistics
    stats_ds = genomic_ds.calculate_windowed_stats(stats=['tajima_d', 'pi'])
    
    summary = genomic_ds.get_summary()
    print(f"Analysis type: {summary['analysis_type']}")
    print(f"Region: {summary['start']}-{summary['end']}")
    print(f"Number of windows: {summary['n_windows']}")
    print(f"Window size: {summary['window_size']}bp")
    
    # Show some statistics
    tajima_d_values = stats_ds['tajima_d'].values
    pi_values = stats_ds['pi'].values
    
    print(f"Tajima's D range: {np.nanmin(tajima_d_values):.3f} to {np.nanmax(tajima_d_values):.3f}")
    print(f"π range: {np.nanmin(pi_values):.3f} to {np.nanmax(pi_values):.3f}")
    
    return genomic_ds


def main():
    """Run all examples."""
    print("🚀 Simplified Windowing Interface Examples")
    print("=" * 60)
    
    try:
        # Run examples
        example_genome_wide()
        example_windowed()
        example_region_specific()
        
        print("\n✅ All examples completed successfully!")
        print("\nKey features demonstrated:")
        print("  • Genome-wide analysis (default when no windows specified)")
        print("  • Windowed analysis with configurable size and step")
        print("  • Region-specific analysis with start/end positions")
        print("  • Automatic window type selection")
        print("  • Minimum variants filtering")
        
        print("\nCommand-line interface would be:")
        print("  # Genome-wide analysis")
        print("  many-stats --vcf data.vcf.gz")
        print("")
        print("  # Windowed analysis")
        print("  many-stats --vcf data.vcf.gz --window-size 1000 --step-size 500")
        print("")
        print("  # Region-specific analysis")
        print("  many-stats --vcf data.vcf.gz --start 1000000 --end 2000000 --window-size 5000")
        
    except Exception as e:
        print(f"\n❌ Examples failed: {e}")
        import traceback
        traceback.print_exc()
        return False
    
    return True


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)

