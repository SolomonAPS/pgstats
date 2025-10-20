"""
Core data structures and analysis framework for population genetics statistics.

This module provides the main classes for managing genomic data and running
windowed and non-windowed population genetics statistics while accounting for
missing data and callable sites.
"""

import numpy as np
import xarray as xr
import pandas as pd
from pathlib import Path
from typing import Optional, Union, List, Dict, Any, Tuple
import dask.array as da
from dataclasses import dataclass

from many_stats.io.loaders import load_vcf_simple
from many_stats.stats.sfs_statistics import tajima_d, fu_li_d, fu_li_f, theta_pi, theta_w, theta_h


@dataclass
class WindowConfig:
    """Configuration for windowed analysis."""
    window_size: int = 1000
    step_size: Optional[int] = None
    min_variants: int = 5
    min_callable_sites: int = 100
    
    def __post_init__(self):
        if self.step_size is None:
            self.step_size = self.window_size


@dataclass
class CallableSitesConfig:
    """Configuration for callable sites handling."""
    bed_file: Optional[str] = None
    max_missing: float = 0.0
    min_coverage: Optional[int] = None
    quality_threshold: Optional[float] = None


class GenomicDataset:
    """
    Main class for managing genomic data and running population genetics analyses.
    
    This class combines sgkit's efficient data handling with sophisticated
    windowing and callable sites functionality.
    """
    
    def __init__(self, 
                 data_source: Union[str, xr.Dataset],
                 callable_config: Optional[CallableSitesConfig] = None,
                 window_config: Optional[WindowConfig] = None):
        """
        Initialize GenomicDataset.
        
        Args:
            data_source: Path to VCF file or existing sgkit Dataset
            callable_config: Configuration for callable sites handling
            window_config: Configuration for windowed analysis
        """
        self.callable_config = callable_config or CallableSitesConfig()
        self.window_config = window_config or WindowConfig()
        
        # Load data
        if isinstance(data_source, str):
            self.dataset = load_vcf_simple(data_source)
        else:
            self.dataset = data_source
        
        # Initialize callable sites mask
        self.callable_mask = None
        self.callable_sites = None
        
        # Initialize window information
        self.windows = None
        self.window_stats = {}
        
        # Validate dataset
        self._validate_dataset()
        
        # Process callable sites if configured
        if self.callable_config.bed_file:
            self._load_callable_sites()
    
    def _validate_dataset(self):
        """Validate that the dataset has required fields."""
        required_vars = ['call_genotype', 'variant_position', 'variant_contig']
        missing_vars = [var for var in required_vars if var not in self.dataset.data_vars]
        
        if missing_vars:
            raise ValueError(f"Dataset missing required variables: {missing_vars}")
        
        # Check for missing data handling
        if 'call_genotype_mask' not in self.dataset.data_vars:
            print("Warning: No genotype mask found. Missing data will be inferred from genotype values.")
    
    def _load_callable_sites(self):
        """Load callable sites from BED file."""
        if not self.callable_config.bed_file:
            return
        
        bed_path = Path(self.callable_config.bed_file)
        if not bed_path.exists():
            raise FileNotFoundError(f"BED file not found: {bed_path}")
        
        # Load BED file
        bed_df = pd.read_csv(
            bed_path, 
            sep='\t', 
            names=['chrom', 'start', 'end'],
            usecols=[0, 1, 2]
        )
        
        # Create callable sites mask
        self.callable_mask = self._create_callable_mask(bed_df)
        self.callable_sites = np.sum(self.callable_mask)
        
        print(f"Loaded callable sites: {self.callable_sites:,} sites")
    
    def _create_callable_mask(self, bed_df: pd.DataFrame) -> np.ndarray:
        """Create boolean mask for callable sites."""
        positions = self.dataset.variant_position.values
        contigs = self.dataset.variant_contig.values
        
        # Get contig names
        contig_names = self.dataset.contig_id.values
        
        mask = np.zeros(len(positions), dtype=bool)
        
        for _, row in bed_df.iterrows():
            chrom = row['chrom']
            start = row['start']
            end = row['end']
            
            # Find matching contig
            contig_idx = None
            for i, contig_name in enumerate(contig_names):
                if str(contig_name) == str(chrom):
                    contig_idx = i
                    break
            
            if contig_idx is not None:
                # Find variants in this region
                # Convert BED (0-based) to VCF (1-based) coordinates
                # BED [start, end) -> VCF positions (start+1) to end (inclusive)
                vcf_start = start + 1
                vcf_end = end
                region_mask = (contigs == contig_idx) & (positions >= vcf_start) & (positions <= vcf_end)
                mask |= region_mask
        
        return mask
    
    def filter_missing_data(self, max_missing: Optional[float] = None):
        """
        Filter variants based on missing data threshold.
        
        Args:
            max_missing: Maximum proportion of missing data allowed
        """
        if max_missing is None:
            max_missing = self.callable_config.max_missing
        
        if max_missing <= 0:
            return
        
        # Calculate missing data proportion per variant
        genotypes = self.dataset.call_genotype.values
        n_samples = genotypes.shape[1]
        
        # Count missing genotypes (assuming -1 or missing values)
        missing_counts = np.sum(genotypes == -1, axis=(1, 2))
        missing_prop = missing_counts / (n_samples * genotypes.shape[2])
        
        # Create filter mask
        filter_mask = missing_prop <= max_missing
        
        # Apply filter
        self.dataset = self.dataset.isel(variants=filter_mask)
        
        # Update callable mask if it exists
        if self.callable_mask is not None:
            self.callable_mask = self.callable_mask[filter_mask]
            self.callable_sites = np.sum(self.callable_mask)
        
        print(f"Filtered to {np.sum(filter_mask):,} variants (removed {np.sum(~filter_mask):,} with >{max_missing:.1%} missing)")
    
    def create_windows(self, 
                      window_size: Optional[int] = None,
                      step_size: Optional[int] = None,
                      min_variants: Optional[int] = None):
        """
        Create genomic windows for analysis.
        
        Args:
            window_size: Size of each window in base pairs
            step_size: Step size for sliding windows
            min_variants: Minimum number of variants per window
        """
        if window_size is None:
            window_size = self.window_config.window_size
        if step_size is None:
            step_size = self.window_config.step_size
        if min_variants is None:
            min_variants = self.window_config.min_variants
        
        positions = self.dataset.variant_position.values
        contigs = self.dataset.variant_contig.values
        
        windows = []
        
        # Group by contig
        unique_contigs = np.unique(contigs)
        
        for contig in unique_contigs:
            contig_mask = contigs == contig
            contig_positions = positions[contig_mask]
            contig_indices = np.where(contig_mask)[0]
            
            if len(contig_positions) == 0:
                continue
            
            # Create windows for this contig
            min_pos = np.min(contig_positions)
            max_pos = np.max(contig_positions)
            
            for start in range(min_pos, max_pos, step_size):
                end = start + window_size
                
                # Find variants in this window
                window_mask = (contig_positions >= start) & (contig_positions < end)
                window_indices = contig_indices[window_mask]
                
                if len(window_indices) >= min_variants:
                    windows.append({
                        'contig': contig,
                        'start': start,
                        'end': end,
                        'variant_indices': window_indices,
                        'n_variants': len(window_indices)
                    })
        
        self.windows = windows
        print(f"Created {len(windows)} windows")
    
    def calculate_windowed_stats(self, 
                               stats: List[str] = None,
                               use_callable_sites: bool = True) -> Dict[str, List[float]]:
        """
        Calculate statistics for each window.
        
        Args:
            stats: List of statistics to calculate
            use_callable_sites: Whether to account for callable sites
            
        Returns:
            Dictionary with statistics for each window
        """
        if self.windows is None:
            raise ValueError("Windows not created. Call create_windows() first.")
        
        if stats is None:
            stats = ['tajima_d', 'pi', 'theta_w', 'theta_h']
        
        # Initialize results
        results = {stat: [] for stat in stats}
        results['window_start'] = []
        results['window_end'] = []
        results['n_variants'] = []
        results['n_callable_sites'] = []
        
        for window in self.windows:
            # Extract window data
            window_dataset = self.dataset.isel(variants=window['variant_indices'])
            
            # Calculate callable sites for this window
            n_callable = None
            if use_callable_sites and self.callable_mask is not None:
                window_callable_mask = self.callable_mask[window['variant_indices']]
                n_callable = np.sum(window_callable_mask)
            
            # Calculate statistics
            window_results = self._calculate_stats_for_window(window_dataset, stats)
            
            # Store results
            for stat in stats:
                results[stat].append(window_results[stat])
            
            results['window_start'].append(window['start'])
            results['window_end'].append(window['end'])
            results['n_variants'].append(window['n_variants'])
            results['n_callable_sites'].append(n_callable)
        
        self.window_stats = results
        return results
    
    def _calculate_stats_for_window(self, window_dataset: xr.Dataset, stats: List[str]) -> Dict[str, float]:
        """Calculate statistics for a single window."""
        results = {}
        
        for stat in stats:
            if stat == 'tajima_d':
                result_ds = tajima_d(window_dataset)
                results[stat] = np.mean(result_ds['tajima_d'].values)
            elif stat == 'pi':
                result_ds = theta_pi(window_dataset)
                results[stat] = np.mean(result_ds['pi'].values)
            elif stat == 'theta_w':
                result_ds = theta_w(window_dataset)
                results[stat] = np.mean(result_ds['theta_w'].values)
            elif stat == 'theta_h':
                result_ds = theta_h(window_dataset)
                results[stat] = np.mean(result_ds['theta_h'].values)
            elif stat == 'fu_li_d':
                result_ds = fu_li_d(window_dataset)
                results[stat] = np.mean(result_ds['fu_li_d'].values)
            elif stat == 'fu_li_f':
                result_ds = fu_li_f(window_dataset)
                results[stat] = np.mean(result_ds['fu_li_f'].values)
            else:
                results[stat] = np.nan
        
        return results
    
    def calculate_genome_wide_stats(self, stats: List[str] = None) -> Dict[str, float]:
        """
        Calculate genome-wide statistics.
        
        Args:
            stats: List of statistics to calculate
            
        Returns:
            Dictionary with genome-wide statistics
        """
        if stats is None:
            stats = ['tajima_d', 'pi', 'theta_w', 'theta_h']
        
        results = {}
        
        for stat in stats:
            if stat == 'tajima_d':
                result_ds = tajima_d(self.dataset)
                results[stat] = np.mean(result_ds['tajima_d'].values)
            elif stat == 'pi':
                result_ds = nucleotide_diversity(self.dataset)
                results[stat] = np.mean(result_ds['pi'].values)
            elif stat == 'theta_w':
                result_ds = watterson_theta(self.dataset)
                results[stat] = np.mean(result_ds['theta_w'].values)
            elif stat == 'theta_h':
                result_ds = fay_wu_theta(self.dataset)
                results[stat] = np.mean(result_ds['theta_h'].values)
            elif stat == 'fu_li_d':
                result_ds = fu_li_d(self.dataset)
                results[stat] = np.mean(result_ds['fu_li_d'].values)
            elif stat == 'fu_li_f':
                result_ds = fu_li_f(self.dataset)
                results[stat] = np.mean(result_ds['fu_li_f'].values)
            else:
                results[stat] = np.nan
        
        return results
    
    def save_results(self, output_path: str, format: str = 'csv'):
        """
        Save analysis results to file.
        
        Args:
            output_path: Path to output file
            format: Output format ('csv', 'tsv', 'parquet')
        """
        if not self.window_stats:
            raise ValueError("No window statistics calculated. Run calculate_windowed_stats() first.")
        
        df = pd.DataFrame(self.window_stats)
        
        if format == 'csv':
            df.to_csv(output_path, index=False)
        elif format == 'tsv':
            df.to_csv(output_path, sep='\t', index=False)
        elif format == 'parquet':
            df.to_parquet(output_path, index=False)
        else:
            raise ValueError(f"Unsupported format: {format}")
        
        print(f"Results saved to {output_path}")
    
    def get_summary(self) -> Dict[str, Any]:
        """Get summary information about the dataset."""
        return {
            'n_variants': len(self.dataset.variants),
            'n_samples': len(self.dataset.samples),
            'n_contigs': len(self.dataset.contigs),
            'callable_sites': self.callable_sites,
            'n_windows': len(self.windows) if self.windows else 0,
            'window_size': self.window_config.window_size,
            'step_size': self.window_config.step_size,
            'max_missing': self.callable_config.max_missing
        }

