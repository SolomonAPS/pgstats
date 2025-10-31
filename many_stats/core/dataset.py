"""
Core data structures and analysis framework for population genetics statistics.

This module provides the main classes for managing genomic data and running
windowed and non-windowed population genetics statistics while accounting for
missing data and callable sites. Built on top of sgkit's windowing capabilities.
"""

import numpy as np
import xarray as xr
import pandas as pd
import pyranges as pr
from pathlib import Path
from typing import Optional, Union, List, Dict, Any, Tuple
import dask.array as da
from dataclasses import dataclass
import sgkit as sg
import time
from functools import wraps

from many_stats.io.loaders import load_vcf_simple
from many_stats.stats.sfs_statistics import (
    tajima_d, fu_li_d, fu_li_f, fu_li_d_unfolded, fu_li_f_unfolded, zeng_e,
    theta_pi, theta_w, theta_h, theta_l, fay_wu_h
)
from many_stats.stats.haplotype_statistics import haplotype_diversity, garud_h_statistics
from many_stats.stats.ld_statistics import (
    calculate_ld_matrix, calculate_windowed_ld, omega_statistic
)


def time_operation(operation_name: str = None):
    """
    Decorator to time a method call and print the result.
    
    Args:
        operation_name: Name to display in timing output
    """
    def decorator(func):
        name = operation_name or func.__name__
        
        @wraps(func)
        def wrapper(self, *args, **kwargs):
            if hasattr(self, 'enable_profiling') and self.enable_profiling:
                start = time.time()
                result = func(self, *args, **kwargs)
                elapsed = time.time() - start
                print(f"[Timing] {name}: {elapsed:.2f}s", flush=True)
            else:
                result = func(self, *args, **kwargs)
            return result
        return wrapper
    return decorator


@dataclass
class WindowConfig:
    """Configuration for windowed analysis - minimal CLI approach."""
    window_size: Optional[int] = None  # Window size in base pairs (None = genome-wide)
    step_size: Optional[int] = None  # Step size (defaults to window_size)
    start: Optional[int] = None  # Start position for region analysis
    end: Optional[int] = None  # End position for region analysis
    min_variants: int = 5  # Minimum variants per window
    
    def __post_init__(self):
        if self.window_size is not None and self.step_size is None:
            self.step_size = self.window_size


@dataclass
class CallableSitesConfig:
    """Configuration for callable sites handling."""
    bed_file: Optional[str] = None
    bed_format: str = "non_callable"  # "non_callable" or "callable"
    max_missing: float = 0.0
    min_coverage: Optional[int] = None
    quality_threshold: Optional[float] = None


class GenomicDataset:
    """
    Main class for managing genomic data and running population genetics analyses.
    
    This class leverages sgkit's built-in windowing capabilities combined with
    callable sites functionality for comprehensive population genetics analysis.
    """
    
    def __init__(self, 
                 data_source: Union[str, xr.Dataset],
                 callable_config: Optional[CallableSitesConfig] = None,
                 window_config: Optional[WindowConfig] = None,
                 keep_zarr: bool = False,
                 zarr_dir: Optional[str] = None,
                 output_dir: Optional[str] = None,
                 enable_profiling: bool = False):
        """
        Initialize GenomicDataset.
        
        Args:
            data_source: Path to VCF file or existing sgkit Dataset
            callable_config: Configuration for callable sites handling
            window_config: Configuration for windowed analysis
            keep_zarr: Whether to keep intermediate Zarr files (VCF only)
            zarr_dir: Custom directory for Zarr files (overrides default logic)
            output_dir: Output directory for analysis (used for Zarr location with --keep-zarr)
            enable_profiling: Enable timing output for performance monitoring
        """
        self.callable_config = callable_config or CallableSitesConfig()
        self.window_config = window_config or WindowConfig()
        self.keep_zarr = keep_zarr
        self.zarr_dir = zarr_dir
        self.output_dir = output_dir
        self.enable_profiling = enable_profiling
        self.data_source_path = data_source if isinstance(data_source, str) else None
        
        # Load data
        if self.enable_profiling:
            start = time.time()
        
        if isinstance(data_source, str):
            data_path = Path(data_source)
            if data_path.suffix in ['.zarr', '.vcz'] or 'zarr' in str(data_path):
                # Load Zarr dataset
                import sgkit as sg
                self.dataset = sg.load_dataset(data_source)
            else:
                # Load VCF
                self.dataset = load_vcf_simple(data_source, keep_zarr=keep_zarr, temp_dir=zarr_dir, output_dir=output_dir)
        else:
            self.dataset = data_source
        
        if self.enable_profiling:
            elapsed = time.time() - start
            print(f"[Timing] Load VCF/Zarr: {elapsed:.2f}s", flush=True)
        
        # Initialize callable sites tracking
        self.callable_sites = None
        self.non_callable_mask = None
        
        # Initialize window information
        self.windowed_dataset = None
        self.window_stats = {}
        
        # Validate dataset
        self._validate_dataset()
        
        # Process callable sites if configured
        if self.callable_config.bed_file:
            if self.enable_profiling:
                start = time.time()
            self._load_callable_sites()
            if self.enable_profiling:
                elapsed = time.time() - start
                print(f"[Timing] Load BED and apply mask: {elapsed:.2f}s", flush=True)
    
    def _validate_dataset(self):
        """Validate that the dataset has required fields."""
        required_vars = ['call_genotype', 'variant_position', 'variant_contig']
        missing_vars = [var for var in required_vars if var not in self.dataset.data_vars]
        
        if missing_vars:
            raise ValueError(f"Dataset missing required variables: {missing_vars}")
        
        # Check for missing data handling
        if 'call_genotype_mask' not in self.dataset.data_vars:
            print("Warning: No genotype mask found. Missing data will be inferred from genotype values.")
    
    def _get_bed_cache_path(self, bed_path: Path) -> Optional[Path]:
        """
        Get the path where the BED mask cache should be stored.
        
        Returns None if caching is not enabled.
        """
        if not self.keep_zarr:
            return None
        
        # Determine cache directory
        if self.zarr_dir:
            cache_dir = Path(self.zarr_dir)
        elif self.output_dir:
            cache_dir = Path(self.output_dir)
        elif self.data_source_path:
            cache_dir = Path(self.data_source_path).parent
        else:
            return None
        
        # Generate hash of BED file content + format
        import hashlib
        hasher = hashlib.md5()
        
        # Include BED file content
        with open(bed_path, 'rb') as f:
            hasher.update(f.read())
        
        # Include BED format (callable vs non_callable matters!)
        hasher.update(self.callable_config.bed_format.encode())
        
        bed_hash = hasher.hexdigest()[:12]  # Use first 12 chars
        
        # Cache filename includes BED filename and hash
        cache_filename = f".bed_mask_{bed_path.stem}_{bed_hash}.npy"
        return cache_dir / cache_filename
    
    def _load_cached_mask(self, cache_path: Path) -> Optional[np.ndarray]:
        """Load cached BED mask if it exists and is valid."""
        if not cache_path.exists():
            return None
        
        try:
            print(f"Found cached BED mask: {cache_path.name}", flush=True)
            cached_mask = np.load(cache_path)
            
            # Validate mask dimensions match current dataset
            if len(cached_mask) != len(self.dataset.variants):
                print(f"  WARNING: Cached mask size mismatch (cache: {len(cached_mask)}, dataset: {len(self.dataset.variants)})", flush=True)
                print(f"  Ignoring cache and recomputing mask...", flush=True)
                return None
            
            print(f"  Using cached mask ({len(cached_mask):,} variants)", flush=True)
            return cached_mask
            
        except Exception as e:
            print(f"  WARNING: Error loading cached mask: {e}", flush=True)
            print(f"  Recomputing mask...", flush=True)
            return None
    
    def _save_cached_mask(self, cache_path: Path, mask: np.ndarray):
        """Save computed BED mask to cache."""
        try:
            cache_path.parent.mkdir(parents=True, exist_ok=True)
            np.save(cache_path, mask)
            print(f"  Cached mask saved: {cache_path.name}", flush=True)
        except Exception as e:
            print(f"  WARNING: Could not save mask cache: {e}", flush=True)
    
    def _load_callable_sites(self, lazy_mode: bool = False):
        """
        Load callable sites from BED file and optionally mask non-callable sites.
        
        Args:
            lazy_mode: If True, only load BED data for L calculation, don't mask yet.
                      Use this for multi-region analysis to avoid masking all variants.
        """
        if not self.callable_config.bed_file:
            return
        
        bed_path = Path(self.callable_config.bed_file)
        if not bed_path.exists():
            raise FileNotFoundError(f"BED file not found: {bed_path}")
        
        # Load BED file for PyRanges (always needed for L calculation)
        bed_df = pd.read_csv(
            bed_path, 
            sep='\t', 
            names=['chrom', 'start', 'end'],
            usecols=[0, 1, 2]
        )
        
        # Normalize chromosome names
        contig_names = self.dataset.contig_id.values
        vcf_chroms = set(contig_names.astype(str))
        bed_chroms = set(bed_df['chrom'].astype(str).unique())
        
        # Sample chromosomes to detect naming convention
        sample_vcf_chr = str(next(iter(vcf_chroms)))
        sample_bed_chr = str(next(iter(bed_chroms)))
        
        bed_df = bed_df.copy()
        if sample_bed_chr.startswith('chr') and not sample_vcf_chr.startswith('chr'):
            print(f"Normalizing: Stripping 'chr' prefix from BED chromosomes", flush=True)
            bed_df['chrom'] = bed_df['chrom'].str.replace('^chr', '', regex=True)
        elif not sample_bed_chr.startswith('chr') and sample_vcf_chr.startswith('chr'):
            print(f"Normalizing: Adding 'chr' prefix to BED chromosomes", flush=True)
            bed_df['chrom'] = 'chr' + bed_df['chrom'].astype(str)
        
        # Store BED data as PyRanges object
        self.bed_gr = pr.PyRanges(bed_df.rename(columns={'chrom': 'Chromosome', 'start': 'Start', 'end': 'End'}))
        
        # If lazy mode, just store BED data without masking
        if lazy_mode:
            print(f"Loaded BED file for L calculation (lazy mode: masking deferred)", flush=True)
            return
        
        # Non-lazy mode: apply mask to entire dataset
        # Try to load cached mask first
        cache_path = self._get_bed_cache_path(bed_path)
        if cache_path:
            cached_mask = self._load_cached_mask(cache_path)
            if cached_mask is not None:
                # Use cached mask directly
                self._apply_cached_mask(cached_mask)
                print(f"Masked non-callable variants: {self.callable_sites:,} variants remaining (from cache)", flush=True)
                return
        
        # Create non-callable sites mask and apply it (bed_df already loaded above)
        callable_mask = self._apply_callable_sites_mask(bed_df)
        
        # Save mask to cache if enabled
        if cache_path:
            self._save_cached_mask(cache_path, callable_mask)
        
        print(f"Masked non-callable variants: {self.callable_sites:,} variants remaining", flush=True)
    
    def _apply_cached_mask(self, callable_mask: np.ndarray):
        """
        Apply a pre-computed callable sites mask.
        
        Args:
            callable_mask: Boolean array where True = callable, False = non-callable
        """
        print(f"Applying cached mask to genotypes...", flush=True)
        
        # Use the cached mask directly
        non_callable_mask = ~callable_mask
        
        # Create broadcasted mask using dask to stay lazy
        mask_3d = da.broadcast_to(
            non_callable_mask[:, np.newaxis, np.newaxis],
            self.dataset.call_genotype.shape
        )
        
        # Create DataArray with same coords as call_genotype to avoid alignment issues
        mask_expanded = xr.DataArray(
            mask_3d,
            coords=self.dataset.call_genotype.coords,
            dims=self.dataset.call_genotype.dims
        )
        
        genotypes_masked = xr.where(mask_expanded, -1, self.dataset.call_genotype)
        self.dataset = self.dataset.assign(call_genotype=genotypes_masked)
        
        updated_mask = (genotypes_masked == -1)
        self.dataset = self.dataset.assign(call_genotype_mask=updated_mask)
        
        # Store information
        self.non_callable_mask = non_callable_mask
        self.callable_sites = np.sum(callable_mask)
        
        print(f"  Masked {np.sum(non_callable_mask):,} non-callable sites as missing data", flush=True)
    
    def _apply_callable_sites_mask(self, bed_df: pd.DataFrame) -> np.ndarray:
        """
        Apply callable sites mask using PyRanges for fast interval tree operations.
        
        BED file format: 0-based start (inclusive), 0-based end (exclusive)
        VCF positions are 1-based.
        
        Args:
            bed_df: DataFrame with BED regions (columns: chrom, start, end)
            
        Returns:
            Boolean array where True = callable, False = non-callable
        """
        print(f"Processing {len(bed_df):,} BED regions with PyRanges...", flush=True)
        
        positions = self.dataset.variant_position.values
        contigs = self.dataset.variant_contig.values
        contig_names = self.dataset.contig_id.values
        
        # Get actual chromosome names by indexing
        unique_contigs = np.unique(contigs)
        actual_contig_names = [str(contig_names[i]) for i in unique_contigs if i < len(contig_names)]
        
        vcf_chroms = set(actual_contig_names)
        bed_chroms = set(bed_df['chrom'].astype(str).unique())
        print(f"VCF chromosomes (first 5): {sorted(vcf_chroms)[:5]}", flush=True)
        print(f"BED chromosomes (first 5): {sorted(bed_chroms)[:5]}", flush=True)
        
        # Normalize chromosome names if needed
        # Check if BED has "chr" prefix but VCF doesn't (or vice versa)
        bed_df = bed_df.copy()
        
        # Sample chromosomes to detect naming convention
        sample_vcf_chr = str(next(iter(vcf_chroms)))
        sample_bed_chr = str(next(iter(bed_chroms)))
        
        if sample_bed_chr.startswith('chr') and not sample_vcf_chr.startswith('chr'):
            # BED has "chr" prefix, VCF doesn't - strip from BED
            print(f"Normalizing: Stripping 'chr' prefix from BED chromosomes", flush=True)
            bed_df['chrom'] = bed_df['chrom'].str.replace('^chr', '', regex=True)
        elif not sample_bed_chr.startswith('chr') and sample_vcf_chr.startswith('chr'):
            # VCF has "chr" prefix, BED doesn't - add to BED
            print(f"Normalizing: Adding 'chr' prefix to BED chromosomes", flush=True)
            bed_df['chrom'] = 'chr' + bed_df['chrom'].astype(str)
        else:
            print(f"Chromosome names match convention", flush=True)
        
        # Create PyRanges object from BED file
        # BED is 0-based, half-open [start, end)
        # But we need to convert to 1-based VCF coordinates for overlap checking
        # BED [0, 4829) should match VCF positions 1-4829, so we add 1 to start
        bed_gr = pr.PyRanges(
            chromosomes=bed_df['chrom'].astype(str),
            starts=bed_df['start'].astype(int) + 1,  # Convert to 1-based
            ends=bed_df['end'].astype(int) + 1       # Make end inclusive in 1-based
        )
        
        # Create DataFrame with variant positions and original indices
        # This is necessary because PyRanges resets indices in overlap results
        # Keep positions as 1-based since we converted BED to 1-based above
        variants_df = pd.DataFrame({
            'Chromosome': contig_names[contigs].astype(str),
            'Start': positions,       # Keep 1-based VCF position  
            'End': positions + 1,     # Make it a half-open interval [pos, pos+1)
            'variant_idx': np.arange(len(positions))  # Track original index
        })
        
        # Create PyRanges object from variants
        variants_gr = pr.PyRanges(variants_df)
        
        # Debug: Show sample data from both PyRanges
        print(f"BED PyRanges sample (first 3 regions):", flush=True)
        print(f"  {bed_gr.head(3).df[['Chromosome', 'Start', 'End']].to_string()}", flush=True)
        print(f"Variants PyRanges sample (first 3 variants):", flush=True)
        print(f"  {variants_gr.head(3).df[['Chromosome', 'Start', 'End']].to_string()}", flush=True)
        
        # Test a manual overlap check for debugging
        bed_first = bed_gr.head(1).df.iloc[0]
        var_first = variants_gr.head(1).df.iloc[0]
        print(f"\nManual overlap test:", flush=True)
        print(f"  First BED: [{bed_first['Start']}, {bed_first['End']}) on {bed_first['Chromosome']}", flush=True)
        print(f"  First variant: [{var_first['Start']}, {var_first['End']}) on {var_first['Chromosome']}", flush=True)
        overlaps = (bed_first['Start'] < var_first['End']) and (var_first['Start'] < bed_first['End']) and (bed_first['Chromosome'] == var_first['Chromosome'])
        print(f"  Should overlap? {overlaps}", flush=True)
        
        # Check chromosome overlap and data types
        bed_chrom_set = set(bed_gr.chromosomes)
        var_chrom_set = set(variants_gr.chromosomes)
        common_chroms = bed_chrom_set.intersection(var_chrom_set)
        print(f"Chromosomes in common: {sorted(common_chroms)[:5]}", flush=True)
        print(f"BED-only chromosomes: {sorted(bed_chrom_set - var_chrom_set)[:5]}", flush=True)
        print(f"VCF-only chromosomes: {sorted(var_chrom_set - bed_chrom_set)[:5]}", flush=True)
        
        # Check data types
        print(f"\nData types:", flush=True)
        print(f"  BED Chromosome dtype: {bed_gr.df['Chromosome'].dtype}", flush=True)
        print(f"  BED Start dtype: {bed_gr.df['Start'].dtype}", flush=True)
        print(f"  Variant Chromosome dtype: {variants_gr.df['Chromosome'].dtype}", flush=True)
        print(f"  Variant Start dtype: {variants_gr.df['Start'].dtype}", flush=True)
        
        print(f"Finding overlaps with interval trees...", flush=True)
        
        # Find overlapping variants using PyRanges' efficient overlap detection
        # This uses NCLS (Nested Containment Lists) internally - very fast!
        overlapping_gr = variants_gr.overlap(bed_gr)
        
        print(f"Found {len(overlapping_gr):,} overlapping variant-region pairs", flush=True)
        
        # Create boolean mask based on which variants overlap BED regions
        if self.callable_config.bed_format == "non_callable":
            # BED defines non-callable regions - start with all callable, remove overlaps
            callable_mask = np.ones(len(positions), dtype=bool)
            if len(overlapping_gr) > 0:
                # Get original indices of overlapping variants
                overlap_indices = overlapping_gr.variant_idx.values
                callable_mask[overlap_indices] = False
            else:
                print("WARNING: No overlaps found between variants and BED regions!", flush=True)
                print("         Check that chromosome names match between VCF and BED file.", flush=True)
        else:  # callable
            # BED defines callable regions - start with none callable, add overlaps
            callable_mask = np.zeros(len(positions), dtype=bool)
            if len(overlapping_gr) > 0:
                # Get original indices of overlapping variants
                overlap_indices = overlapping_gr.variant_idx.values
                callable_mask[overlap_indices] = True
            else:
                print("WARNING: No overlaps found between variants and BED regions!", flush=True)
                print("         Check that chromosome names match between VCF and BED file.", flush=True)
        
        n_non_callable = np.sum(~callable_mask)
        n_callable = np.sum(callable_mask)
        print(f"BED processing complete: {n_callable:,} callable, {n_non_callable:,} non-callable variants", flush=True)
        
        # Set non-callable sites to -1 (missing) in genotype data
        # Use xarray.where() to stay lazy with dask arrays
        non_callable_mask = ~callable_mask
        
        # Create broadcasted mask using dask to stay lazy
        mask_3d = da.broadcast_to(
            non_callable_mask[:, np.newaxis, np.newaxis],
            self.dataset.call_genotype.shape
        )
        
        # Create DataArray with same coords as call_genotype to avoid alignment issues
        mask_expanded = xr.DataArray(
            mask_3d,
            coords=self.dataset.call_genotype.coords,
            dims=self.dataset.call_genotype.dims
        )
        
        genotypes_masked = xr.where(mask_expanded, -1, self.dataset.call_genotype)
        
        # Update the dataset
        self.dataset = self.dataset.assign(call_genotype=genotypes_masked)
        
        # Update the call_genotype_mask to reflect the new missing data
        # Create mask where True = missing, False = callable (stays lazy)
        updated_mask = (genotypes_masked == -1)
        self.dataset = self.dataset.assign(call_genotype_mask=updated_mask)
        
        # Store information
        self.non_callable_mask = non_callable_mask
        self.callable_sites = np.sum(callable_mask)
        
        print(f"Masked {np.sum(non_callable_mask):,} non-callable sites as missing data", flush=True)
        
        # Return the callable mask for caching
        return callable_mask
    
    def _apply_mask_to_region(self, region_dataset: xr.Dataset, region_start: int, region_end: int, contig_name: str) -> xr.Dataset:
        """
        Apply BED mask to variants within a specific region only (optimized for multi-region analysis).
        
        Args:
            region_dataset: Filtered dataset containing only variants in the region
            region_start: Start position of the region (1-based)
            region_end: End position of the region (1-based)
            contig_name: Name of the contig/chromosome
            
        Returns:
            Dataset with BED mask applied to region variants only
        """
        if not hasattr(self, 'bed_gr') or self.bed_gr is None:
            print("No BED file loaded, skipping region masking")
            return region_dataset
        
        if self.enable_profiling:
            start = time.time()
        
        positions = region_dataset.variant_position.values
        contigs = region_dataset.variant_contig.values
        contig_names = region_dataset.contig_id.values
        
        # Create PyRanges object for variants in this region
        variants_df = pd.DataFrame({
            'Chromosome': [contig_name] * len(positions),
            'Start': positions,
            'End': positions + 1,
            'variant_idx': np.arange(len(positions))
        })
        
        variants_gr = pr.PyRanges(variants_df)
        
        # Find overlaps with BED regions
        overlapping_gr = variants_gr.overlap(self.bed_gr)
        
        # Create callable mask for this region's variants
        if self.callable_config.bed_format == "non_callable":
            callable_mask = np.ones(len(positions), dtype=bool)
            if len(overlapping_gr) > 0:
                overlap_indices = overlapping_gr.variant_idx.values
                callable_mask[overlap_indices] = False
        else:  # callable
            callable_mask = np.zeros(len(positions), dtype=bool)
            if len(overlapping_gr) > 0:
                overlap_indices = overlapping_gr.variant_idx.values
                callable_mask[overlap_indices] = True
        
        # Apply mask to genotypes
        non_callable_mask = ~callable_mask
        mask_3d = da.broadcast_to(
            non_callable_mask[:, np.newaxis, np.newaxis],
            region_dataset.call_genotype.shape
        )
        
        mask_expanded = xr.DataArray(
            mask_3d,
            coords=region_dataset.call_genotype.coords,
            dims=region_dataset.call_genotype.dims
        )
        
        genotypes_masked = xr.where(mask_expanded, -1, region_dataset.call_genotype)
        region_dataset = region_dataset.assign(call_genotype=genotypes_masked)
        
        updated_mask = (genotypes_masked == -1)
        region_dataset = region_dataset.assign(call_genotype_mask=updated_mask)
        
        n_masked = np.sum(non_callable_mask)
        n_callable = np.sum(callable_mask)
        
        if self.enable_profiling:
            elapsed = time.time() - start
            print(f"[Timing] Region BED masking: {elapsed:.2f}s ({n_masked:,} masked)", flush=True)
        else:
            print(f"  Masked {n_masked:,} non-callable sites in region ({n_callable:,} callable)", flush=True)
        
        return region_dataset
    
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
        
        # Update non-callable mask if it exists
        if self.non_callable_mask is not None:
            self.non_callable_mask = self.non_callable_mask[filter_mask]
            self.callable_sites = np.sum(~self.non_callable_mask)
        
        print(f"Filtered to {np.sum(filter_mask):,} variants (removed {np.sum(~filter_mask):,} with >{max_missing:.1%} missing)")
    
    def create_windows(self):
        """
        Create genomic windows using sgkit's windowing functions.
        
        Automatically chooses between genome-wide and windowed analysis based on configuration.
        """
        window_size = self.window_config.window_size
        step_size = self.window_config.step_size
        start = self.window_config.start
        end = self.window_config.end
        
        if window_size is None and start is None and end is None:
            # Genome-wide analysis
            print("Creating genome-wide window", flush=True)
            try:
                # window_by_genome creates a single window per contig
                # We need to ensure it handles contig boundaries correctly
                self.windowed_dataset = sg.window_by_genome(self.dataset)
            except Exception as e:
                raise ValueError(f"Failed to create genome-wide window: {e}")
        else:
            # Windowed analysis
            if start is not None and end is not None:
                print(f"Creating windows for region {start}-{end}", flush=True)
                
                # Check if dataset is already filtered (when called from calculate_stats_for_regions)
                # If the dataset variants fall within the region bounds, skip re-filtering
                positions = self.dataset.variant_position.values
                dataset_already_filtered = (
                    len(positions) > 0 and
                    positions[0] >= start and 
                    positions[-1] <= end
                )
                
                if dataset_already_filtered:
                    # Dataset is already filtered to this region, use it directly
                    region_dataset = self.dataset
                else:
                    # Filter dataset to region
                    region_dataset = self._filter_to_region(start, end)
                
                # Check if region has any variants
                if len(region_dataset.variants) == 0:
                    raise ValueError(f"No variants found in region {start}-{end}")
                
                try:
                    # Use window_by_position for base pair windows
                    # Could also use window_by_variant for variant count windows
                    self.windowed_dataset = sg.window_by_position(
                        region_dataset,
                        size=window_size or 1000,
                        step=step_size or window_size or 1000
                    )
                except Exception as e:
                    raise ValueError(f"Failed to create windows for region {start}-{end}: {e}")
            else:
                print(f"Creating position-based windows with size {window_size}", flush=True)
                try:
                    # Use window_by_position for base pair windows
                    # Could also use window_by_variant for variant count windows
                    self.windowed_dataset = sg.window_by_position(
                        self.dataset,
                        size=window_size,
                        step=step_size
                    )
                except Exception as e:
                    raise ValueError(f"Failed to create position-based windows: {e}")
        
        # Check if windowing was successful
        if self.windowed_dataset is None:
            print("Warning: Window creation failed. No windows available.", flush=True)
            return
            
        n_windows = len(self.windowed_dataset.windows)
        print(f"Created {n_windows} windows", flush=True)
        
        # Set up window coordinates properly
        
        # Get window information from sgkit output
        window_starts = self.windowed_dataset.window_start.values
        window_stops = self.windowed_dataset.window_stop.values
        
        # Handle window_contig (may not exist for single-contig datasets)
        if 'window_contig' in self.windowed_dataset.data_vars or 'window_contig' in self.windowed_dataset.coords:
            window_contigs = self.windowed_dataset.window_contig.values
        else:
            # Single contig - use contig 0 for all windows
            window_contigs = np.zeros(n_windows, dtype=int)
        
        # Get actual positions for each window
        positions = self.dataset.variant_position.values
        window_start_positions = positions[window_starts]
        window_stop_positions = positions[window_stops - 1]  # -1 because stop is exclusive
        
        # Set up proper window coordinates
        self.windowed_dataset = self.windowed_dataset.assign_coords({
            'windows': np.arange(n_windows),
            'window_start': ('windows', window_start_positions),  # Use actual positions
            'window_stop': ('windows', window_stop_positions),
            'window_contig': ('windows', window_contigs)
        })
        
        # Add window indices as data variables for LD calculation
        self.windowed_dataset['window_start_idx'] = ('windows', window_starts)
        self.windowed_dataset['window_stop_idx'] = ('windows', window_stops)
        
        
        # Check for boundary issues and warn user
        self._check_window_boundaries()
        
        # DON'T filter windows here - do it after stats calculation to avoid dimension mismatches
        # Filter windows by minimum variants if specified (store mask for later)
        # Allow min_variants=0 to keep all windows and let users filter later
        if self.window_config.min_variants >= 0:
            self._filter_windows_by_variants()
    
    def _check_window_boundaries(self):
        """Check for windows that extend beyond contig boundaries and warn user."""
        if self.windowed_dataset is None:
            return
        
        window_size = self.window_config.window_size
        if window_size is None:
            return  # No windowing, so no boundary issues
        
        # Get contig information - check if contig_length exists
        if 'contig_length' not in self.dataset.data_vars:
            print("Warning: No contig_length information available. Skipping boundary checks.")
            return
            
        contig_lengths = self.dataset.contig_length.values
        contig_ids = self.dataset.contig_id.values
        
        # Get window information (use indices, not positions)
        window_contigs = self.windowed_dataset.window_contig.values
        window_starts = self.windowed_dataset.window_start_idx.values
        window_stops = self.windowed_dataset.window_stop_idx.values
        
        # Check each window
        boundary_warnings = []
        small_window_warnings = []
        
        for i, (contig, start_idx, stop_idx) in enumerate(zip(window_contigs, window_starts, window_stops)):
            contig_length = contig_lengths[contig]
            contig_id = contig_ids[contig]
            
            # Calculate actual window size in base pairs
            if stop_idx > 0 and start_idx < len(self.dataset.variant_position):
                # Get actual genomic positions
                positions = self.dataset.variant_position.values
                start_pos = positions[start_idx] if start_idx < len(positions) else contig_length
                stop_pos = positions[stop_idx-1] if stop_idx > 0 and stop_idx-1 < len(positions) else contig_length
                actual_size = stop_pos - start_pos + 1
            else:
                actual_size = 0
            
            # Check if window extends beyond contig boundary
            if stop_idx > 0 and start_idx < len(self.dataset.variant_position):
                positions = self.dataset.variant_position.values
                if start_idx < len(positions):
                    start_pos = positions[start_idx]
                    if start_pos + window_size > contig_length:
                        boundary_warnings.append({
                            'window': i + 1,
                            'contig': contig_id,
                            'start_pos': start_pos,
                            'requested_end': start_pos + window_size,
                            'contig_length': contig_length,
                            'actual_size': actual_size
                        })
            
            # Check if window has very few variants compared to expected size
            n_variants = stop_idx - start_idx
            if n_variants > 0 and actual_size < window_size * 0.5:  # Less than 50% of expected size
                small_window_warnings.append({
                    'window': i + 1,
                    'contig': contig_id,
                    'n_variants': n_variants,
                    'actual_size': actual_size,
                    'expected_size': window_size
                })
        
        # Print warnings
        if boundary_warnings:
            print(f"\nWARNING: {len(boundary_warnings)} windows extend beyond contig boundaries:")
            for warning in boundary_warnings[:5]:  # Show first 5
                print(f"  Window {warning['window']} on {warning['contig']}: "
                      f"requested {warning['start_pos']}-{warning['requested_end']} "
                      f"(contig length: {warning['contig_length']}, actual size: {warning['actual_size']}bp)")
            if len(boundary_warnings) > 5:
                print(f"  ... and {len(boundary_warnings) - 5} more windows")
            print("  These windows will only include variants within the contig boundaries.")
        
        if small_window_warnings:
            print(f"\nWARNING: {len(small_window_warnings)} windows have fewer variants than expected:")
            for warning in small_window_warnings[:5]:  # Show first 5
                print(f"  Window {warning['window']} on {warning['contig']}: "
                      f"{warning['n_variants']} variants, {warning['actual_size']}bp "
                      f"(expected: {warning['expected_size']}bp)")
            if len(small_window_warnings) > 5:
                print(f"  ... and {len(small_window_warnings) - 5} more windows")
            print("  Consider using smaller window sizes or checking variant density.")
    
    def _filter_windows_by_variants(self):
        """
        Store information about windows that should be filtered.
        
        This method marks windows with too few variants rather than filtering
        them immediately, allowing stats calculations to work correctly.
        The actual filtering happens after stats are calculated.
        """
        if self.windowed_dataset is None:
            return
        
        # Count variants per window (use indices, not positions)
        window_starts = self.windowed_dataset.window_start_idx.values
        window_stops = self.windowed_dataset.window_stop_idx.values
        n_variants_per_window = window_stops - window_starts
        
        # Create filter mask
        self.window_filter_mask = n_variants_per_window >= self.window_config.min_variants
        
        n_filtered = np.sum(self.window_filter_mask)
        n_removed = len(self.window_filter_mask) - n_filtered
        
        if self.window_config.min_variants == 0:
            print(f"Keeping all {n_filtered} windows (min_variants=0; users can filter by n_variants column)")
        else:
            print(f"Will filter to {n_filtered} windows (removing {n_removed} with <{self.window_config.min_variants} variants after stats calculation)")
    
    def _filter_to_region(self, start: int, end: int) -> xr.Dataset:
        """Filter dataset to a specific genomic region."""
        positions = self.dataset.variant_position.values
        region_mask = (positions >= start) & (positions <= end)
        
        # Check if any variants match the region
        n_variants = np.sum(region_mask)
        if n_variants == 0:
            print(f"Filtered to region {start}-{end}: 0 variants")
            # Return empty dataset with same structure
            return self.dataset.isel(variants=slice(0, 0))
        
        filtered_dataset = self.dataset.isel(variants=region_mask)
        print(f"Filtered to region {start}-{end}: {n_variants} variants")
        
        return filtered_dataset
    
    def calculate_windowed_stats(self, 
                               stats: List[str] = None,
                               use_callable_sites: bool = True) -> xr.Dataset:
        """
        Calculate statistics for each window using sgkit's windowed dataset.
        
        Args:
            stats: List of statistics to calculate
            use_callable_sites: Whether to account for callable sites
            
        Returns:
            Dataset with windowed statistics
        """
        if self.windowed_dataset is None:
            raise ValueError("Windows not created. Call create_windows() first.")
        
        if stats is None:
            stats = [
                # SFS statistics
                'tajima_d', 'fu_li_d', 'fu_li_f', 'fu_li_d_unfolded', 'fu_li_f_unfolded',
                'zeng_e', 'fay_wu_h',
                # Theta estimators
                'theta_pi', 'theta_w', 'theta_h', 'theta_l',
                # LD statistics
                'ld_d', 'ld_dprime', 'ld_r2', 'omega_statistic',
                # Haplotype statistics
                'haplotype_diversity', 'garud_h1', 'garud_h12', 'garud_h123', 'garud_h2_h1',
                # Singleton counts
                'singletons', 'singletons_unfolded'
            ]
        
        # Filter windowed_dataset to only include windows with sufficient variants
        # This must happen BEFORE calculating statistics to avoid NaN values
        window_starts = self.windowed_dataset.window_start_idx.values
        window_stops = self.windowed_dataset.window_stop_idx.values
        n_variants_per_window = window_stops - window_starts
        
        # Store the number of windows BEFORE filtering
        n_windows_before_filtering = len(self.windowed_dataset.windows)
        
        # Apply min_variants filter if set
        if self.window_config.min_variants > 0:
            valid_windows_mask = n_variants_per_window >= self.window_config.min_variants
            # Filter the windowed dataset to only valid windows
            self.windowed_dataset = self.windowed_dataset.isel(windows=valid_windows_mask)
            # Update n_variants array to match filtered dataset
            n_variants_per_window = n_variants_per_window[valid_windows_mask]
        
        # Start with a copy of the FILTERED windowed_dataset
        result_dataset = self.windowed_dataset.copy()
        
        # Add n_variants per window
        result_dataset['n_variants'] = (['windows'], n_variants_per_window)
        
        
        # Calculate each statistic
        for stat in stats:
            
            # Time each statistic
            if self.enable_profiling:
                stat_start = time.time()
            
            if stat == 'tajima_d':
                stat_ds = tajima_d(self.windowed_dataset)
                result_dataset = result_dataset.merge(stat_ds, compat='override')
                
                # Print timing
                if self.enable_profiling:
                    elapsed = time.time() - stat_start
                    print(f"[Timing] {stat}: {elapsed:.2f}s", flush=True)
            elif stat in ['theta_pi', 'pi', 'nucleotide_diversity']:
                stat_ds = theta_pi(self.windowed_dataset)
                result_dataset = result_dataset.merge(stat_ds, compat='override')
                if self.enable_profiling:
                    elapsed = time.time() - stat_start
                    print(f"[Timing] {stat}: {elapsed:.2f}s", flush=True)
            elif stat in ['theta_w', 'watterson_theta']:
                stat_ds = theta_w(self.windowed_dataset)
                result_dataset = result_dataset.merge(stat_ds, compat='override')
            elif stat in ['theta_h', 'fay_wu_theta']:
                stat_ds = theta_h(self.windowed_dataset)
                result_dataset = result_dataset.merge(stat_ds, compat='override')
            elif stat == 'theta_l':
                stat_ds = theta_l(self.windowed_dataset)
                result_dataset = result_dataset.merge(stat_ds, compat='override')
            elif stat == 'fu_li_d':
                # Calculate both folded (D*) and unfolded (D) versions
                stat_ds_folded = fu_li_d(self.windowed_dataset, folded=True)
                stat_ds_unfolded = fu_li_d(self.windowed_dataset, folded=False)
                result_dataset = result_dataset.merge(stat_ds_folded)
                result_dataset = result_dataset.merge(stat_ds_unfolded)
            elif stat == 'fu_li_f':
                # Calculate both folded (F*) and unfolded (F) versions
                stat_ds_folded = fu_li_f(self.windowed_dataset, folded=True)
                stat_ds_unfolded = fu_li_f(self.windowed_dataset, folded=False)
                result_dataset = result_dataset.merge(stat_ds_folded)
                result_dataset = result_dataset.merge(stat_ds_unfolded)
            elif stat == 'fu_li_d_unfolded':
                # Alias for fu_li_d - already calculated above
                pass
            elif stat == 'fu_li_f_unfolded':
                # Alias for fu_li_f - already calculated above
                pass
            elif stat == 'zeng_e':
                stat_ds = zeng_e(self.windowed_dataset)
                print(f"DEBUG: zeng_e stat_ds variables: {list(stat_ds.data_vars)}", flush=True)
                print(f"DEBUG: zeng_e first 3 values: {stat_ds['zeng_e'].values[:3]}", flush=True)
                result_dataset = result_dataset.merge(stat_ds, compat='override')
                print(f"DEBUG: After merge, result_dataset has zeng_e: {'zeng_e' in result_dataset.data_vars}", flush=True)
                
                if self.enable_profiling:
                    elapsed = time.time() - stat_start
                    print(f"[Timing] {stat}: {elapsed:.2f}s", flush=True)
            elif stat == 'fay_wu_h':
                stat_ds = fay_wu_h(self.windowed_dataset)
                print(f"DEBUG: fay_wu_h stat_ds variables: {list(stat_ds.data_vars)}", flush=True)
                print(f"DEBUG: fay_wu_h first 3 values: {stat_ds['fay_wu_h'].values[:3]}", flush=True)
                result_dataset = result_dataset.merge(stat_ds, compat='override')
                print(f"DEBUG: After merge, result_dataset has fay_wu_h: {'fay_wu_h' in result_dataset.data_vars}", flush=True)
                
                if self.enable_profiling:
                    elapsed = time.time() - stat_start
                    print(f"[Timing] {stat}: {elapsed:.2f}s", flush=True)
            elif stat == 'ld_matrix':
                # Calculate LD matrix for each window
                ld_ds = calculate_windowed_ld(self.windowed_dataset)
                result_dataset = result_dataset.merge(ld_ds)
            # Handle all LD statistics at once if any are requested
            elif stat in ['ld_d', 'ld_dprime', 'ld_r2']:
                # Only calculate LD once per window
                if not hasattr(self, '_cached_ld_stats'):
                    
                    # Ensure window coordinates are properly set
                    if 'windows' not in self.windowed_dataset.dims:
                        self.windowed_dataset = self.windowed_dataset.assign_coords({
                            'windows': np.arange(len(self.windowed_dataset.window_start)),
                            'window_start': ('windows', self.windowed_dataset.window_start),
                            'window_stop': ('windows', self.windowed_dataset.window_stop)
                        })
                    
                    # Calculate all LD statistics at once
                    self._cached_ld_stats = calculate_windowed_ld(self.windowed_dataset, enable_profiling=self.enable_profiling)
                    
                    if self.enable_profiling:
                        if len(self._cached_ld_stats.windows) != len(self.windowed_dataset.windows):
                            pass
                
                # Map stat name to LD dataset variable
                ld_var_map = {
                    'ld_d': 'mean_D',
                    'ld_dprime': 'mean_D_prime',
                    'ld_r2': 'mean_r_squared'
                }
                
                # Get the actual variable name from the cached LD stats
                var_name = ld_var_map.get(stat)
                if var_name and var_name in self._cached_ld_stats:
                    # Create a new dataset with just the requested statistic
                    stat_ds = xr.Dataset({
                        var_name: self._cached_ld_stats[var_name]
                    }, coords=self._cached_ld_stats.coords)
                    
                    result_dataset = result_dataset.merge(stat_ds, compat='override')
                
                if self.enable_profiling:
                    elapsed = time.time() - stat_start
                    print(f"[Timing] {stat}: {elapsed:.2f}s", flush=True)
            elif stat == 'haplotype_diversity':
                stat_ds = haplotype_diversity(self.windowed_dataset)
                result_dataset = result_dataset.merge(stat_ds, compat='override')
            elif stat in ['garud_h1', 'garud_h12', 'garud_h123', 'garud_h2_h1']:
                # Calculate all Garud H statistics at once if any are requested
                if not hasattr(self, '_cached_garud_stats'):
                    self._cached_garud_stats = garud_h_statistics(self.windowed_dataset)
                
                # Extract the requested statistic
                stat_ds = xr.Dataset({
                    stat: self._cached_garud_stats[stat]
                }, coords=self._cached_garud_stats.coords)
                
                result_dataset = result_dataset.merge(stat_ds, compat='override')
            elif stat == 'omega_statistic':
                stat_ds = omega_statistic(self.windowed_dataset, enable_profiling=self.enable_profiling)
                print(f"DEBUG: omega_statistic stat_ds variables: {list(stat_ds.data_vars)}", flush=True)
                print(f"DEBUG: omega_statistic first 3 values: {stat_ds['omega_statistic'].values[:3]}", flush=True)
                result_dataset = result_dataset.merge(stat_ds, compat='override')
                print(f"DEBUG: After merge, result_dataset has omega_statistic: {'omega_statistic' in result_dataset.data_vars}", flush=True)
                
                if self.enable_profiling:
                    elapsed = time.time() - stat_start
                    print(f"[Timing] {stat}: {elapsed:.2f}s", flush=True)
        
        # Note: Callable sites are now handled by setting non-callable sites to -1
        # sgkit's statistics will automatically account for missing data
        
        # Normalize theta estimators by callable sites
        print(f"DEBUG: About to normalize theta by callable sites, result_dataset has {len(result_dataset.windows)} windows", flush=True)
        result_dataset_before_norm = result_dataset.copy()
        result_dataset = self._normalize_theta_by_callable_sites(result_dataset)
        print(f"DEBUG: Finished normalizing theta by callable sites", flush=True)
        
        # Now filter out windows where LD statistics failed to calculate
        # (min_variants filtering already done before stat calculation)
        
        # Check LD statistics that should have valid values if calculation succeeded
        ld_stats = {
            'mean_D': 'average D (linkage disequilibrium)',
            'mean_D_prime': "average D' (standardized D)",
            'mean_r_squared': 'average r² (squared correlation)',
            'max_r_squared': 'maximum r² in window'
        }
        valid_stats = []
        
        # Check each LD statistic
        for var, description in ld_stats.items():
            if var in result_dataset.data_vars:
                
                # Convert to float array to handle NaN checks
                values = result_dataset[var].values.astype(float)
                valid = ~np.isnan(values)  # True if LD calculation succeeded
                valid_stats.append(valid)
                
                if self.enable_profiling:
                    n_valid = np.sum(valid)
                    n_failed = len(valid) - n_valid
                    if n_failed > 0:
                        pass
        
        if valid_stats:
            # Window is valid if any LD statistic was calculated successfully
            valid_windows = np.any(valid_stats, axis=0)
            n_valid = np.sum(valid_windows)
            n_failed = len(valid_windows) - n_valid
            
            
            result_dataset = result_dataset.isel(windows=valid_windows)
            print(f"Filtered to {len(result_dataset.windows)} windows with valid LD statistics", flush=True)
        
        # Print mean statistics
        # Extract only window-related variables to avoid memory issues
        try:
            window_vars = {var: result_dataset[var].values for var in result_dataset.data_vars if var in ['n_variants'] + stats}
            stats_to_print = [stat for stat in stats if stat in window_vars]
            if stats_to_print:
                print(f"\nMean Statistics (windowed analysis):")
                print("=" * 50)
                for stat in stats_to_print:
                    values = result_dataset[stat].values
                    valid_mask = ~np.isnan(values)
                    n_valid = np.sum(valid_mask)
                    total = len(values)
                    if n_valid > 0:
                        mean_val = np.mean(values[valid_mask])
                        print(f"  {stat:20s}: {mean_val:.6f} (n={n_valid}/{total} windows)")
                    else:
                        print(f"  {stat:20s}: No valid values")
                print("=" * 50)
        except Exception as e:
            print(f"Note: Could not calculate mean statistics: {e}")
        
        self.window_stats = result_dataset
        return result_dataset
    
    
    def calculate_genome_wide_stats(self, stats: List[str] = None) -> Dict[str, float]:
        """
        Calculate genome-wide statistics.
        
        Args:
            stats: List of statistics to calculate
            
        Returns:
            Dictionary with genome-wide statistics
        """
        if stats is None:
            stats = [
                # SFS statistics
                'tajima_d', 'fu_li_d', 'fu_li_f', 'fu_li_d_unfolded', 'fu_li_f_unfolded',
                'zeng_e', 'fay_wu_h',
                # Theta estimators
                'theta_pi', 'theta_w', 'theta_h', 'theta_l',
                # LD statistics
                'ld_d', 'ld_dprime', 'ld_r2', 'omega_statistic',
                # Haplotype statistics
                'haplotype_diversity', 'garud_h1', 'garud_h12', 'garud_h123', 'garud_h2_h1',
                # Singleton counts
                'singletons', 'singletons_unfolded'
            ]
        
        results = {}
        
        for stat in stats:
            if stat == 'tajima_d':
                result_ds = tajima_d(self.dataset)
                results[stat] = np.mean(result_ds['tajima_d'].values)
            elif stat in ['theta_pi', 'pi', 'nucleotide_diversity']:
                result_ds = theta_pi(self.dataset)
                results[stat] = np.mean(result_ds['theta_pi'].values)
            elif stat in ['theta_w', 'watterson_theta']:
                result_ds = theta_w(self.dataset)
                results[stat] = np.mean(result_ds['theta_w'].values)
            elif stat in ['theta_h', 'fay_wu_theta']:
                result_ds = theta_h(self.dataset)
                results[stat] = np.mean(result_ds['theta_h'].values)
            elif stat == 'theta_l':
                result_ds = theta_l(self.dataset)
                results[stat] = np.mean(result_ds['theta_l'].values)
            elif stat == 'fu_li_d':
                result_ds = fu_li_d(self.dataset)
                results[stat] = np.mean(result_ds['fu_li_d_star'].values)
            elif stat == 'fu_li_f':
                result_ds = fu_li_f(self.dataset)
                results[stat] = np.mean(result_ds['fu_li_f_star'].values)
            elif stat == 'fu_li_d_unfolded':
                result_ds = fu_li_d_unfolded(self.dataset)
                results[stat] = np.mean(result_ds['fu_li_d'].values)
            elif stat == 'fu_li_f_unfolded':
                result_ds = fu_li_f_unfolded(self.dataset)
                results[stat] = np.mean(result_ds['fu_li_f'].values)
            elif stat == 'zeng_e':
                result_ds = zeng_e(self.dataset)
                results[stat] = np.mean(result_ds['zeng_e'].values)
            elif stat == 'fay_wu_h':
                result_ds = fay_wu_h(self.dataset)
                results[stat] = np.mean(result_ds['fay_wu_h'].values)
        
        # Print mean statistics
        print(f"\nGenome-wide Statistics:")
        print("=" * 50)
        for stat in stats:
            if stat in results and not pd.isna(results[stat]):
                print(f"  {stat:20s}: {results[stat]:.6f}")
            else:
                print(f"  {stat:20s}: No valid values")
        print("=" * 50)
        
        return results
    
    def _normalize_theta_by_callable_sites(self, window_stats: xr.Dataset) -> xr.Dataset:
        """
        Normalize theta estimators by the number of callable sites in each window.
        
        This ensures that theta estimates are per-site rather than per-window,
        accounting for the actual number of analyzable sites in each window.
        
        Args:
            window_stats: Dataset with windowed statistics (already filtered)
            
        Returns:
            Dataset with normalized theta estimators
        """
        print(f"DEBUG: _normalize_theta_by_callable_sites called with {len(window_stats.windows)} windows", flush=True)
        result = window_stats.copy()
        
        # Use vectorized batch calculation if BED file available
        if hasattr(self, 'bed_gr') and self.bed_gr is not None and self.callable_config.bed_file:
            print(f"DEBUG: Using batch callable sites calculation", flush=True)
            callable_sites_per_window = self._calculate_callable_lengths_batch(window_stats)
        else:
            # Fall back to per-window calculation (no BED file or old method)
            print(f"DEBUG: Using per-window callable sites calculation (no BED or fallback)", flush=True)
            callable_sites_per_window = self._calculate_callable_lengths_per_window(window_stats)
        
        callable_sites_per_window = np.array(callable_sites_per_window)
        
        # Normalize theta estimators by callable sites
        theta_stats = ['theta_pi', 'theta_w', 'theta_h', 'theta_l']
        
        # Check that we have the right number of callable sites
        num_windows = len(result.windows)
        num_callable_sites = len(callable_sites_per_window)
        
        if num_windows != num_callable_sites:
            print(f"Warning: Window count mismatch: {num_windows} windows vs {num_callable_sites} callable sites. Skipping normalization.", flush=True)
            return result
        
        for stat in theta_stats:
            if stat in result.data_vars:
                # Debug: Check values before normalization
                if stat == 'theta_pi':
                    print(f"DEBUG theta normalization: First 3 windows before normalization:", flush=True)
                    print(f"  theta_pi values: {result[stat].values[:3]}", flush=True)
                    print(f"  callable_sites: {callable_sites_per_window[:3]}", flush=True)
                
                # Normalize each window's values by its callable sites
                for window_idx in range(num_windows):
                    if callable_sites_per_window[window_idx] > 0:
                        result[stat].values[window_idx] /= callable_sites_per_window[window_idx]
                    else:
                        # If 0 callable sites, set to NaN (window is entirely non-callable)
                        result[stat].values[window_idx] = np.nan
                
                # Debug: Check values after normalization
                if stat == 'theta_pi':
                    print(f"  After normalization: {result[stat].values[:3]}", flush=True)
        
        return result
    
    def _calculate_callable_lengths_batch(self, window_stats: xr.Dataset) -> np.ndarray:
        """
        Vectorized calculation of callable sites for all windows at once.
        
        This is much faster than calculating callable length for each window individually.
        
        Args:
            window_stats: Dataset with windowed statistics
            
        Returns:
            Array of callable site counts per window
        """
        if self.enable_profiling:
            start = time.time()
        
        num_windows = len(window_stats.windows)
        positions = self.dataset.variant_position.values
        contigs = self.dataset.variant_contig.values
        contig_names = self.dataset.contig_id.values
        
        # Convert window indices to genomic positions
        window_start_positions = []
        window_end_positions = []
        window_contigs = []
        
        for window_idx in range(num_windows):
            # Use window_start_idx and window_stop_idx (variant array indices)
            # NOT window_start/window_stop (which are genomic positions after reassignment)
            window_start_idx = window_stats.window_start_idx.values[window_idx]
            window_stop_idx = window_stats.window_stop_idx.values[window_idx]
            
            if window_start_idx < len(positions) and window_stop_idx > 0:
                window_start_pos = positions[window_start_idx]
                window_end_pos = positions[window_stop_idx-1] if window_stop_idx-1 < len(positions) else positions[-1]
                window_length = window_end_pos - window_start_pos + 1
                
                # Determine contig for this window
                contig_idx = contigs[window_start_idx]
                contig_name = contig_names[contig_idx]
                
                window_start_positions.append(window_start_pos)
                window_end_positions.append(window_end_pos)
                window_contigs.append(contig_name)
            else:
                # Empty window
                window_start_positions.append(0)
                window_end_positions.append(0)
                window_contigs.append(None)
        
        # Create PyRanges for all windows at once
        windows_df = pd.DataFrame({
            'Chromosome': window_contigs,
            'Start': window_start_positions,
            'End': window_end_positions,
            'window_idx': range(num_windows)
        })
        
        # Remove windows with None contig
        valid_windows = windows_df['Chromosome'].notna()
        windows_gr = pr.PyRanges(windows_df[valid_windows])
        
        if len(windows_gr) == 0:
            return np.zeros(num_windows)
        
        # Get BED regions
        bed_df = self.bed_gr.df.copy()
        
        if self.callable_config.bed_format == "non_callable":
            # Calculate masked length per window
            callable_lengths = []
            
            for window_idx in range(num_windows):
                if window_contigs[window_idx] is None:
                    callable_lengths.append(0)
                    continue
                
                window_start = window_start_positions[window_idx]
                window_end = window_end_positions[window_idx]
                window_length = window_end - window_start
                
                # Filter BED regions for this contig
                contig_bed = bed_df[bed_df['Chromosome'] == window_contigs[window_idx]]
                
                # Calculate overlap
                masked_length = 0
                for _, row in contig_bed.iterrows():
                    overlap_start = max(window_start, row['Start'])
                    overlap_end = min(window_end, row['End'])
                    if overlap_start < overlap_end:
                        masked_length += overlap_end - overlap_start
                
                callable_length = window_length - masked_length
                callable_lengths.append(max(1, callable_length))
            
            if self.enable_profiling:
                elapsed = time.time() - start
                print(f"[Timing] Vectorized callable length calculation: {elapsed:.2f}s ({num_windows} windows)", flush=True)
            
            return np.array(callable_lengths)
            
        elif self.callable_config.bed_format == "callable":
            # Calculate callable length per window
            callable_lengths = []
            
            for window_idx in range(num_windows):
                if window_contigs[window_idx] is None:
                    callable_lengths.append(0)
                    continue
                
                window_start = window_start_positions[window_idx]
                window_end = window_end_positions[window_idx]
                
                # Filter BED regions for this contig
                contig_bed = bed_df[bed_df['Chromosome'] == window_contigs[window_idx]]
                
                # Calculate overlap
                callable_length = 0
                for _, row in contig_bed.iterrows():
                    overlap_start = max(window_start, row['Start'])
                    overlap_end = min(window_end, row['End'])
                    if overlap_start < overlap_end:
                        callable_length += overlap_end - overlap_start
                
                callable_lengths.append(max(1, callable_length))
            
            if self.enable_profiling:
                elapsed = time.time() - start
                print(f"[Timing] Vectorized callable length calculation: {elapsed:.2f}s ({num_windows} windows)", flush=True)
            
            return np.array(callable_lengths)
        
        else:
            raise ValueError(f"Unknown bed_format: {self.callable_config.bed_format}")
    
    def _calculate_callable_lengths_per_window(self, window_stats: xr.Dataset) -> List[int]:
        """
        Calculate callable sites per window using the old per-window method (for fallback).
        
        Args:
            window_stats: Dataset with windowed statistics
            
        Returns:
            List of callable site counts per window
        """
        if self.enable_profiling:
            start = time.time()
        
        callable_sites_per_window = []
        
        positions = self.dataset.variant_position.values
        contigs = self.dataset.variant_contig.values
        contig_names = self.dataset.contig_id.values
        
        for window_idx in range(len(window_stats.windows)):
            window_start_idx = window_stats.window_start.values[window_idx]
            window_stop_idx = window_stats.window_stop.values[window_idx]
            
            if window_start_idx < len(positions) and window_stop_idx > 0:
                window_start_pos = positions[window_start_idx]
                window_stop_pos = positions[window_stop_idx-1] if window_stop_idx-1 < len(positions) else positions[-1]
                
                # Determine contig for this window
                contig_idx = contigs[window_start_idx]
                contig_name = contig_names[contig_idx]
            else:
                window_start_pos = 0
                window_stop_pos = 0
                contig_name = None
            
            # Calculate callable length
            callable_length = self._calculate_callable_length_in_window(window_start_pos, window_stop_pos, contig_name)
            callable_sites_per_window.append(callable_length)
        
        if self.enable_profiling:
            elapsed = time.time() - start
            print(f"[Timing] Per-window callable length calculation: {elapsed:.2f}s ({len(window_stats.windows)} windows)", flush=True)
        
        return callable_sites_per_window
    
    def _calculate_callable_length_in_window(self, window_start: int, window_stop: int, contig_name: str = None) -> int:
        """
        Calculate the number of callable sites in a window based on BED overlap using PyRanges.
        
        Args:
            window_start: Start position of the window
            window_stop: End position of the window
            contig_name: Name of the contig (chromosome) for this window
            
        Returns:
            Number of callable sites in the window
        """
        if not hasattr(self, 'callable_config') or self.callable_config is None:
            # No BED file provided, entire window is callable
            return window_stop - window_start
        
        window_length = window_stop - window_start
        
        if self.callable_config.bed_file is None:
            return window_length
        
        # Check if PyRanges BED data is available
        if not hasattr(self, 'bed_gr') or self.bed_gr is None:
            # Fallback to original method if PyRanges data not available
            return self._calculate_callable_length_in_window_fallback(window_start, window_stop, contig_name)
        
        if contig_name is None:
            # If no contig specified, assume entire window is callable
            return window_length
        
        if self.callable_config.bed_format == "non_callable":
            # BED defines non-callable regions
            # Calculate masked length by finding intersections with BED regions on this contig
            masked_length = 0
            
            # Get BED regions on this contig
            bed_df = self.bed_gr.df
            overlapping_bed = bed_df[bed_df['Chromosome'] == contig_name]
            
            for _, bed_row in overlapping_bed.iterrows():
                bed_start = bed_row['Start']
                bed_end = bed_row['End']
                
                # Calculate intersection
                intersection_start = max(window_start, bed_start)
                intersection_end = min(window_stop, bed_end)
                
                if intersection_start < intersection_end:
                    masked_length += intersection_end - intersection_start
            
            callable_length = window_length - masked_length
            
        elif self.callable_config.bed_format == "callable":
            # BED defines callable regions
            # Calculate callable length by finding intersections with BED regions on this contig
            callable_length = 0
            
            # Get BED regions on this contig
            bed_df = self.bed_gr.df
            overlapping_bed = bed_df[bed_df['Chromosome'] == contig_name]
            
            for _, bed_row in overlapping_bed.iterrows():
                bed_start = bed_row['Start']
                bed_end = bed_row['End']
                
                # Calculate intersection
                intersection_start = max(window_start, bed_start)
                intersection_end = min(window_stop, bed_end)
                
                if intersection_start < intersection_end:
                    callable_length += intersection_end - intersection_start
            
        else:
            raise ValueError(f"Unknown bed_format: {self.callable_config.bed_format}")
        
        # Ensure at least 1 callable site to avoid division by zero
        return max(1, callable_length)
    
    def _calculate_callable_length_in_window_fallback(self, window_start: int, window_stop: int, contig_name: str = None) -> int:
        """
        Fallback method for calculating callable length when PyRanges is not available.
        This is the original implementation for backward compatibility.
        """
        window_length = window_stop - window_start
        
        # Load BED file
        bed_df = pd.read_csv(
            self.callable_config.bed_file,
            sep='\t',
            header=None,
            names=['chrom', 'start', 'end'],
            usecols=[0, 1, 2]
        )
        
        if self.callable_config.bed_format == "non_callable":
            # BED defines non-callable regions
            masked_length = 0
            
            for _, row in bed_df.iterrows():
                chrom = row['chrom']
                bed_start = row['start']
                bed_end = row['end']
                
                # Only process BED regions on the same contig
                if contig_name is not None and str(chrom) == str(contig_name):
                    # Calculate overlap between BED region and window
                    overlap_start = max(window_start, bed_start)
                    overlap_end = min(window_stop, bed_end)
                    
                    if overlap_start < overlap_end:
                        masked_length += overlap_end - overlap_start
            
            callable_length = window_length - masked_length
            
        elif self.callable_config.bed_format == "callable":
            # BED defines callable regions
            callable_length = 0
            
            for _, row in bed_df.iterrows():
                chrom = row['chrom']
                bed_start = row['start']
                bed_end = row['end']
                
                # Only process BED regions on the same contig
                if contig_name is not None and str(chrom) == str(contig_name):
                    # Calculate overlap between BED region and window
                    overlap_start = max(window_start, bed_start)
                    overlap_end = min(window_stop, bed_end)
                    
                    if overlap_start < overlap_end:
                        callable_length += overlap_end - overlap_start
            
        else:
            raise ValueError(f"Unknown bed_format: {self.callable_config.bed_format}")
        
        # Ensure at least 1 callable site to avoid division by zero
        return max(1, callable_length)
    
    def save_results(self, output_path: str, format: str = 'csv'):
        """
        Save analysis results to file.
        
        Args:
            output_path: Path to output file
            format: Output format ('csv', 'tsv', 'parquet', 'zarr')
        """
        if not self.window_stats:
            raise ValueError("No window statistics calculated. Run calculate_windowed_stats() first.")
        
        if format == 'zarr':
            # Save as Zarr dataset
            self.window_stats.to_zarr(output_path)
        else:
            # Extract only statistics and window metadata
            data = {}
            
            # Add window information
            if 'window_contig' in self.window_stats.coords:
                data['window_contig'] = self.window_stats.window_contig.values
            if 'window_start' in self.window_stats.coords:
                data['window_start'] = self.window_stats.window_start.values
            if 'window_stop' in self.window_stats.coords:
                data['window_stop'] = self.window_stats.window_stop.values
            
            # Add statistics (only variables that look like statistics)
            stat_names = [
                'tajima_d', 'fu_li_d', 'fu_li_f', 'zeng_e', 'fay_wu_h',
                'theta_pi', 'theta_w', 'theta_h', 'theta_l',
                # Backward compatibility
                'nucleotide_diversity', 'watterson_theta', 'fay_wu_theta', 'pi'
            ]
            
            for var_name in self.window_stats.data_vars:
                if any(stat in str(var_name).lower() for stat in stat_names):
                    data[var_name] = self.window_stats[var_name].values
            
            # Create DataFrame
            df = pd.DataFrame(data)
            
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
        summary = {
            'n_variants': len(self.dataset.variants),
            'n_samples': len(self.dataset.samples),
            'n_contigs': len(self.dataset.contigs),
            'callable_sites': self.callable_sites,
            'window_size': self.window_config.window_size,
            'step_size': self.window_config.step_size,
            'start': self.window_config.start,
            'end': self.window_config.end,
            'max_missing': self.callable_config.max_missing
        }
        
        if self.windowed_dataset is not None:
            summary['n_windows'] = len(self.windowed_dataset.windows)
            summary['analysis_type'] = 'genome-wide' if len(self.windowed_dataset.windows) == 1 else 'windowed'
        else:
            summary['n_windows'] = 0
            summary['analysis_type'] = 'not_analyzed'
        
        return summary
    
    def _print_statistic_means(self, result_dataset: xr.Dataset, stats: List[str], analysis_type: str = "analysis"):
        """
        Calculate and print mean values for all statistics.
        
        Args:
            result_dataset: Dataset with calculated statistics
            stats: List of statistic names to calculate means for
            analysis_type: Description of the analysis type for output
        """
        print(f"\nMean Statistics ({analysis_type}):")
        print("=" * 50)
        
        # Map stat names to dataset variables
        # Note: fu_li_d and fu_li_f now calculate BOTH folded and unfolded versions
        stat_map = {
            'tajima_d': 'tajima_d',
            'fu_li_d': ['fu_li_d_star', 'fu_li_d'],  # Both folded (D*) and unfolded (D)
            'fu_li_f': ['fu_li_f_star', 'fu_li_f'],  # Both folded (F*) and unfolded (F)
            'fu_li_d_unfolded': 'fu_li_d',
            'fu_li_f_unfolded': 'fu_li_f',
            'zeng_e': 'zeng_e',
            'fay_wu_h': 'fay_wu_h',
            'theta_pi': 'theta_pi',
            'theta_w': 'theta_w',
            'theta_h': 'theta_h',
            'theta_l': 'theta_l',
            'ld_d': 'mean_D',
            'ld_dprime': 'mean_D_prime',
            'ld_r2': 'mean_r_squared',
            'omega_statistic': 'omega_statistic',
            'haplotype_diversity': 'haplotype_diversity',
            'garud_h1': 'garud_h1',
            'garud_h12': 'garud_h12',
            'garud_h123': 'garud_h123',
            'garud_h2_h1': 'garud_h2_h1'
        }
        
        for stat in stats:
            var_names = stat_map.get(stat)
            if var_names is None:
                print(f"  {stat:20s}: Not calculated")
                continue
            
            # Handle both single variable names and lists (for fu_li_d/f)
            if isinstance(var_names, list):
                # Print both folded and unfolded versions
                for var_name in var_names:
                    if var_name in result_dataset:
                        values = result_dataset[var_name].values
                        valid_mask = ~np.isnan(values)
                        n_windows = np.sum(valid_mask)
                        total_windows = len(values)
                        
                        if n_windows > 0:
                            mean_val = np.mean(values[valid_mask])
                            # Add suffix to distinguish folded vs unfolded
                            suffix = " (folded)" if "_star" in var_name else " (unfolded)"
                            print(f"  {stat + suffix:20s}: {mean_val:.6f} (n={n_windows}/{total_windows} windows)")
                        else:
                            suffix = " (folded)" if "_star" in var_name else " (unfolded)"
                            print(f"  {stat + suffix:20s}: No valid values")
            else:
                # Single variable name
                if var_names in result_dataset:
                    # Calculate mean, excluding NaN values
                    values = result_dataset[var_names].values
                    valid_mask = ~np.isnan(values)
                    n_windows = np.sum(valid_mask)
                    total_windows = len(values)
                    
                    if n_windows > 0:
                        mean_val = np.mean(values[valid_mask])
                        print(f"  {stat:20s}: {mean_val:.6f} (n={n_windows}/{total_windows} windows)")
                    else:
                        print(f"  {stat:20s}: No valid values")
                else:
                    print(f"  {stat:20s}: Not calculated")
        
        print("=" * 50)
    
    def calculate_stats_for_regions(self,
                                    regions: List[Tuple[str, int, int]],
                                    window_size: int,
                                    step_size: Optional[int] = None,
                                    stats: List[str] = None,
                                    min_variants: int = 1,
                                    use_callable_sites: bool = True) -> pd.DataFrame:
        """
        Calculate statistics for multiple genomic regions.
        
        This method efficiently analyzes multiple regions by loading the dataset
        once and applying region-specific filtering for each analysis. It's
        particularly useful for analyzing candidate genes, exons, or homologous
        regions across chromosomes.
        
        Args:
            regions: List of (contig, start, end) tuples (1-based, inclusive)
                    Example: [("chr1", 1000000, 2000000), ("chr5", 3000000, 4000000)]
            window_size: Window size in bp for sliding window analysis
            step_size: Step size for sliding windows (default: window_size)
            stats: List of statistics to calculate (default: all available)
            min_variants: Minimum variants per window (default: 1)
            use_callable_sites: Whether to use callable sites mask (default: True)
            
        Returns:
            DataFrame with columns:
            - region_contig, region_start, region_end: Region identifiers
            - window_start, window_end: Window boundaries
            - n_variants: Number of variants in window
            - <stat_name>: Calculated statistics
            
        Example:
            >>> regions = [("chr1", 1000000, 2000000), ("chr5", 3000000, 4000000)]
            >>> results = genomic_ds.calculate_stats_for_regions(
            ...     regions=regions,
            ...     window_size=50000,
            ...     stats=['tajima_d', 'nucleotide_diversity']
            ... )
            >>> print(results.columns)
            ['region_contig', 'region_start', 'region_end', 'window_start', 
             'window_end', 'n_variants', 'tajima_d', 'nucleotide_diversity']
        """
        import pandas as pd
        
        if not regions:
            raise ValueError("regions list cannot be empty")
        
        if stats is None:
            stats = [
                # SFS statistics
                'tajima_d', 'fu_li_d', 'fu_li_f', 'fu_li_d_unfolded', 'fu_li_f_unfolded',
                'zeng_e', 'fay_wu_h',
                # Theta estimators
                'theta_pi', 'theta_w', 'theta_h', 'theta_l',
                # LD statistics
                'ld_d', 'ld_dprime', 'ld_r2', 'omega_statistic',
                # Haplotype statistics
                'haplotype_diversity', 'garud_h1', 'garud_h12', 'garud_h123', 'garud_h2_h1',
                # Singleton counts
                'singletons', 'singletons_unfolded'
            ]
        
        if step_size is None:
            step_size = window_size
        
        all_results = []
        
        print(f"\nAnalyzing {len(regions)} regions with {window_size:,}bp windows...")
        
        for i, (contig, region_start, region_end) in enumerate(regions, 1):
            if self.enable_profiling:
                region_start_time = time.time()
            
            print(f"\n[{i}/{len(regions)}] Analyzing region {contig}:{region_start:,}-{region_end:,}")
            
            # Update window config for this region
            self.window_config.start = region_start
            self.window_config.end = region_end
            self.window_config.window_size = window_size
            self.window_config.step_size = step_size
            self.window_config.min_variants = min_variants
            
            # Clear cached statistics from previous region
            if hasattr(self, '_cached_ld_stats'):
                delattr(self, '_cached_ld_stats')
            if hasattr(self, '_cached_garud_stats'):
                delattr(self, '_cached_garud_stats')
            
            try:
                # Filter to region first
                if self.enable_profiling:
                    start = time.time()
                region_dataset = self._filter_to_region(region_start, region_end)
                if self.enable_profiling:
                    elapsed = time.time() - start
                    print(f"[Timing] Filter to region: {elapsed:.2f}s", flush=True)
                
                # Apply BED mask to region variants only (lazy masking optimization)
                if use_callable_sites and hasattr(self, 'bed_gr') and self.bed_gr is not None:
                    region_dataset = self._apply_mask_to_region(region_dataset, region_start, region_end, contig)
                
                # Check if any variants remain in the region
                n_variants = len(region_dataset.variants)
                if n_variants == 0:
                    print(f"  WARNING: No variants in region {contig}:{region_start:,}-{region_end:,}, skipping", flush=True)
                    continue
                
                # Temporarily swap dataset to use region_dataset
                original_dataset = self.dataset
                self.dataset = region_dataset
                
                # Create windows for this region
                self.create_windows()
                
                # Calculate statistics
                region_stats = self.calculate_windowed_stats(
                    stats=stats,
                    use_callable_sites=use_callable_sites
                )
                
                # Convert to DataFrame - only extract window-related variables
                # Build dict of only the variables we need (avoid loading all variant/sample data)
                region_data = {}
                
                # Add window coordinates
                if 'window_start' in region_stats.coords:
                    region_data['window_start'] = region_stats.window_start.values
                if 'window_stop' in region_stats.coords:
                    region_data['window_stop'] = region_stats.window_stop.values
                if 'window_contig' in region_stats.coords:
                    region_data['window_contig'] = region_stats.window_contig.values
                
                # Add window statistics only
                # Map stat names to actual variable names in dataset
                stat_to_var = {
                    'tajima_d': ['tajima_d'],
                    'fu_li_d': ['fu_li_d_star', 'fu_li_d'],  # Both folded and unfolded
                    'fu_li_f': ['fu_li_f_star', 'fu_li_f'],  # Both folded and unfolded
                    'fu_li_d_unfolded': ['fu_li_d'],
                    'fu_li_f_unfolded': ['fu_li_f'],
                    'zeng_e': ['zeng_e'],
                    'fay_wu_h': ['fay_wu_h'],
                    'theta_pi': ['theta_pi'],
                    'theta_w': ['theta_w'],
                    'theta_h': ['theta_h'],
                    'theta_l': ['theta_l'],
                    'ld_d': ['mean_D'],
                    'ld_dprime': ['mean_D_prime'],
                    'ld_r2': ['mean_r_squared'],
                    'omega_statistic': ['omega_statistic'],
                    'haplotype_diversity': ['haplotype_diversity'],
                    'garud_h1': ['garud_h1'],
                    'garud_h12': ['garud_h12'],
                    'garud_h123': ['garud_h123'],
                    'garud_h2_h1': ['garud_h2_h1']
                }
                
                for stat in stats:
                    var_names = stat_to_var.get(stat, [stat])
                    for var_name in var_names:
                        if var_name in region_stats.data_vars:
                            region_data[var_name] = region_stats[var_name].values
                
                # Add n_variants if present
                if 'n_variants' in region_stats.data_vars:
                    region_data['n_variants'] = region_stats.n_variants.values
                
                # Create DataFrame
                region_df = pd.DataFrame(region_data)
                
                # Add region identifiers
                region_df.insert(0, 'region_contig', contig)
                region_df.insert(1, 'region_start', region_start)
                region_df.insert(2, 'region_end', region_end)
                
                all_results.append(region_df)
                
                print(f"  Analyzed {len(region_df)} windows")
                
                # Restore original dataset for next region
                self.dataset = original_dataset
                self.windowed_dataset = None
                
                if self.enable_profiling:
                    elapsed = time.time() - region_start_time
                    print(f"[Timing] Total region time: {elapsed:.2f}s", flush=True)
                
            except Exception as e:
                import traceback
                print(f"  ERROR: analyzing region {contig}:{region_start}-{region_end}: {e}")
                if self.enable_profiling:
                    print(f"  Full traceback:")
                    traceback.print_exc()
                # Restore original dataset even on error
                if 'original_dataset' in locals():
                    self.dataset = original_dataset
                continue
        
        if not all_results:
            raise RuntimeError("No regions were successfully analyzed")
        
        # Combine all results
        combined = pd.concat(all_results, ignore_index=True)
        
        print(f"\nMulti-region analysis complete: {len(combined)} total windows")
        print(f"  Regions analyzed: {len(all_results)}")
        print(f"  Average windows per region: {len(combined) / len(all_results):.1f}")
        
        # Print mean statistics
        self._print_statistic_means(combined, stats, f"multi-region ({len(all_results)} regions)")
        
        return combined
