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

from many_stats.io.loaders import load_vcf_simple
from many_stats.stats.sfs_statistics import (
    tajima_d, fu_li_d, fu_li_f, fu_li_d_unfolded, fu_li_f_unfolded, zeng_e,
    theta_pi, theta_w, theta_h, theta_l, fay_wu_h
)
from many_stats.stats.ld_statistics import (
    calculate_ld_matrix, calculate_windowed_ld
)


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
                 output_dir: Optional[str] = None):
        """
        Initialize GenomicDataset.
        
        Args:
            data_source: Path to VCF file or existing sgkit Dataset
            callable_config: Configuration for callable sites handling
            window_config: Configuration for windowed analysis
            keep_zarr: Whether to keep intermediate Zarr files (VCF only)
            zarr_dir: Custom directory for Zarr files (overrides default logic)
            output_dir: Output directory for analysis (used for Zarr location with --keep-zarr)
        """
        self.callable_config = callable_config or CallableSitesConfig()
        self.window_config = window_config or WindowConfig()
        self.keep_zarr = keep_zarr
        self.zarr_dir = zarr_dir
        self.output_dir = output_dir
        self.data_source_path = data_source if isinstance(data_source, str) else None
        
        # Load data
        if isinstance(data_source, str):
            data_path = Path(data_source)
            if data_path.suffix in ['.zarr'] or 'zarr' in str(data_path):
                # Load Zarr dataset
                import sgkit as sg
                self.dataset = sg.load_dataset(data_source)
            else:
                # Load VCF
                self.dataset = load_vcf_simple(data_source, keep_zarr=keep_zarr, temp_dir=zarr_dir, output_dir=output_dir)
        else:
            self.dataset = data_source
        
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
    
    def _load_callable_sites(self):
        """Load callable sites from BED file and mask non-callable sites."""
        if not self.callable_config.bed_file:
            return
        
        bed_path = Path(self.callable_config.bed_file)
        if not bed_path.exists():
            raise FileNotFoundError(f"BED file not found: {bed_path}")
        
        # Try to load cached mask first
        cache_path = self._get_bed_cache_path(bed_path)
        if cache_path:
            cached_mask = self._load_cached_mask(cache_path)
            if cached_mask is not None:
                # Use cached mask directly
                self._apply_cached_mask(cached_mask)
                print(f"Applied callable sites mask: {self.callable_sites:,} callable sites (from cache)", flush=True)
                return
        
        # Load BED file
        bed_df = pd.read_csv(
            bed_path, 
            sep='\t', 
            names=['chrom', 'start', 'end'],
            usecols=[0, 1, 2]
        )
        
        # Create non-callable sites mask and apply it
        callable_mask = self._apply_callable_sites_mask(bed_df)
        
        # Save mask to cache if enabled
        if cache_path:
            self._save_cached_mask(cache_path, callable_mask)
        
        print(f"Applied callable sites mask: {self.callable_sites:,} callable sites", flush=True)
    
    def _apply_cached_mask(self, callable_mask: np.ndarray):
        """
        Apply a pre-computed callable sites mask.
        
        Args:
            callable_mask: Boolean array where True = callable, False = non-callable
        """
        print(f"Applying cached mask to genotypes...", flush=True)
        
        # Use the cached mask directly
        non_callable_mask = ~callable_mask
        
        # Expand mask and apply to genotypes (lazy operation)
        mask_expanded = xr.DataArray(
            non_callable_mask[:, np.newaxis, np.newaxis],
            dims=['variants', 'samples', 'ploidy'],
            coords={
                'variants': self.dataset.variants,
                'samples': self.dataset.samples,
                'ploidy': self.dataset.ploidy
            }
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
        
        # Debug: Show chromosome names
        vcf_chroms = set(contig_names.astype(str))
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
        bed_gr = pr.PyRanges(
            chromosomes=bed_df['chrom'].astype(str),
            starts=bed_df['start'].astype(int),
            ends=bed_df['end'].astype(int)
        )
        
        # Create DataFrame with variant positions and original indices
        # This is necessary because PyRanges resets indices in overlap results
        variants_df = pd.DataFrame({
            'Chromosome': contig_names[contigs].astype(str),
            'Start': positions - 1,  # Convert 1-based to 0-based
            'End': positions,         # VCF position becomes end (exclusive)
            'variant_idx': np.arange(len(positions))  # Track original index
        })
        
        # Create PyRanges object from variants
        variants_gr = pr.PyRanges(variants_df)
        
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
        
        # Expand mask to match genotype dimensions (variants, samples, ploidy)
        # non_callable_mask is shape (variants,), need to broadcast to (variants, samples, ploidy)
        mask_expanded = xr.DataArray(
            non_callable_mask[:, np.newaxis, np.newaxis],
            dims=['variants', 'samples', 'ploidy'],
            coords={
                'variants': self.dataset.variants,
                'samples': self.dataset.samples,
                'ploidy': self.dataset.ploidy
            }
        )
        
        # Set all genotypes at non-callable sites to -1 (stays lazy with dask)
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
                self.windowed_dataset = sg.window_by_genome(self.dataset)
            except Exception as e:
                raise ValueError(f"Failed to create genome-wide window: {e}")
        else:
            # Windowed analysis
            if start is not None and end is not None:
                print(f"Creating windows for region {start}-{end}", flush=True)
                # Filter dataset to region first
                region_dataset = self._filter_to_region(start, end)
                
                # Check if region has any variants
                if len(region_dataset.variants) == 0:
                    raise ValueError(f"No variants found in region {start}-{end}")
                
                try:
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
        
        # Check for boundary issues and warn user
        self._check_window_boundaries()
        
        # Filter windows by minimum variants if specified
        if self.window_config.min_variants > 1:
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
        
        # Get window information
        window_contigs = self.windowed_dataset.window_contig.values
        window_starts = self.windowed_dataset.window_start.values
        window_stops = self.windowed_dataset.window_stop.values
        
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
        """Filter windows to only include those with minimum number of variants."""
        if self.windowed_dataset is None:
            return
        
        # Count variants per window
        window_starts = self.windowed_dataset.window_start.values
        window_stops = self.windowed_dataset.window_stop.values
        n_variants_per_window = window_stops - window_starts
        
        # Create filter mask
        filter_mask = n_variants_per_window >= self.window_config.min_variants
        
        # Apply filter
        self.windowed_dataset = self.windowed_dataset.isel(windows=filter_mask)
        
        n_filtered = np.sum(filter_mask)
        n_removed = len(filter_mask) - n_filtered
        print(f"Filtered to {n_filtered} windows (removed {n_removed} with <{self.window_config.min_variants} variants)")
    
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
            stats = ['tajima_d', 'theta_pi', 'theta_w', 'theta_h']
        
        # Start with the windowed dataset
        result_dataset = self.windowed_dataset.copy()
        
        # Calculate each statistic
        for stat in stats:
            if stat == 'tajima_d':
                stat_ds = tajima_d(self.windowed_dataset)
                result_dataset = result_dataset.merge(stat_ds)
            elif stat in ['theta_pi', 'pi', 'nucleotide_diversity']:
                stat_ds = theta_pi(self.windowed_dataset)
                result_dataset = result_dataset.merge(stat_ds)
            elif stat in ['theta_w', 'watterson_theta']:
                stat_ds = theta_w(self.windowed_dataset)
                result_dataset = result_dataset.merge(stat_ds)
            elif stat in ['theta_h', 'fay_wu_theta']:
                stat_ds = theta_h(self.windowed_dataset)
                result_dataset = result_dataset.merge(stat_ds)
            elif stat == 'theta_l':
                stat_ds = theta_l(self.windowed_dataset)
                result_dataset = result_dataset.merge(stat_ds)
            elif stat == 'fu_li_d':
                stat_ds = fu_li_d(self.windowed_dataset)
                result_dataset = result_dataset.merge(stat_ds)
            elif stat == 'fu_li_f':
                stat_ds = fu_li_f(self.windowed_dataset)
                result_dataset = result_dataset.merge(stat_ds)
            elif stat == 'fu_li_d_unfolded':
                stat_ds = fu_li_d_unfolded(self.windowed_dataset)
                result_dataset = result_dataset.merge(stat_ds)
            elif stat == 'fu_li_f_unfolded':
                stat_ds = fu_li_f_unfolded(self.windowed_dataset)
                result_dataset = result_dataset.merge(stat_ds)
            elif stat == 'zeng_e':
                stat_ds = zeng_e(self.windowed_dataset)
                result_dataset = result_dataset.merge(stat_ds)
            elif stat == 'fay_wu_h':
                stat_ds = fay_wu_h(self.windowed_dataset)
                result_dataset = result_dataset.merge(stat_ds)
            elif stat == 'ld_matrix':
                # Calculate LD matrix for each window
                ld_ds = calculate_windowed_ld(self.windowed_dataset)
                result_dataset = result_dataset.merge(ld_ds)
            elif stat == 'ld_d':
                # Calculate D values for each window
                ld_ds = calculate_windowed_ld(self.windowed_dataset)
                result_dataset = result_dataset.merge(ld_ds)
            elif stat == 'ld_d_prime':
                # Calculate D' values for each window
                ld_ds = calculate_windowed_ld(self.windowed_dataset)
                result_dataset = result_dataset.merge(ld_ds)
            elif stat == 'ld_r_squared':
                # Calculate r² values for each window
                ld_ds = calculate_windowed_ld(self.windowed_dataset)
                result_dataset = result_dataset.merge(ld_ds)
        
        # Note: Callable sites are now handled by setting non-callable sites to -1
        # sgkit's statistics will automatically account for missing data
        
        # Normalize theta estimators by callable sites
        result_dataset = self._normalize_theta_by_callable_sites(result_dataset)
        
        # Print mean statistics
        try:
            results_df = result_dataset.to_dataframe()
            self._print_statistic_means(results_df, stats, "windowed analysis")
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
            stats = ['tajima_d', 'theta_pi', 'theta_w', 'theta_h']
        
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
            window_stats: Dataset with windowed statistics
            
        Returns:
            Dataset with normalized theta estimators
        """
        result = window_stats.copy()
        
        # Calculate callable sites per window based on BED overlap
        callable_sites_per_window = []
        
        for window_idx in range(len(window_stats.windows)):
            window_start_idx = window_stats.window_start.values[window_idx]
            window_stop_idx = window_stats.window_stop.values[window_idx]
            
            # Convert variant indices to genomic positions
            positions = self.dataset.variant_position.values
            if window_start_idx < len(positions) and window_stop_idx > 0:
                window_start_pos = positions[window_start_idx]
                window_stop_pos = positions[window_stop_idx-1] if window_stop_idx-1 < len(positions) else positions[-1]
                window_length = window_stop_pos - window_start_pos + 1
            else:
                # Empty window
                window_start_pos = 0
                window_stop_pos = 0
                window_length = 0
            
            # Calculate callable length in this window using genomic positions
            callable_length = self._calculate_callable_length_in_window(window_start_pos, window_stop_pos)
            callable_sites_per_window.append(callable_length)
        
        callable_sites_per_window = np.array(callable_sites_per_window)
        
        # Normalize theta estimators by callable sites
        theta_stats = ['theta_pi', 'theta_w', 'theta_h', 'theta_l']
        
        for stat in theta_stats:
            if stat in result.data_vars:
                # Normalize each window's values by its callable sites
                for window_idx in range(len(result.windows)):
                    if callable_sites_per_window[window_idx] > 0:
                        result[stat].values[window_idx] /= callable_sites_per_window[window_idx]
        
        return result
    
    def _calculate_callable_length_in_window(self, window_start: int, window_stop: int) -> int:
        """
        Calculate the number of callable sites in a window based on BED overlap.
        
        Args:
            window_start: Start position of the window
            window_stop: End position of the window
            
        Returns:
            Number of callable sites in the window
        """
        if not hasattr(self, 'callable_config') or self.callable_config is None:
            # No BED file provided, entire window is callable
            return window_stop - window_start
        
        window_length = window_stop - window_start
        
        if self.callable_config.bed_file is None:
            return window_length
        
        # Load BED file
        bed_df = pd.read_csv(
            self.callable_config.bed_file,
            sep='\t',
            header=None,
            names=['chrom', 'start', 'end'],
            usecols=[0, 1, 2]
        )
        
        # Get contig names
        contig_names = self.dataset.contig_id.values
        
        if self.callable_config.bed_format == "non_callable":
            # BED defines non-callable regions (current behavior)
            # Calculate masked length in this window
            masked_length = 0
            
            for _, row in bed_df.iterrows():
                chrom = row['chrom']
                bed_start = row['start']
                bed_end = row['end']
                
                # Find matching contig
                contig_idx = None
                for i, contig_name in enumerate(contig_names):
                    if str(contig_name) == str(chrom):
                        contig_idx = i
                        break
                
                if contig_idx is not None:
                    # Calculate overlap between BED region and window
                    overlap_start = max(window_start, bed_start)
                    overlap_end = min(window_stop, bed_end)
                    
                    if overlap_start < overlap_end:
                        masked_length += overlap_end - overlap_start
            
            # Calculate callable length
            callable_length = window_length - masked_length
            
        elif self.callable_config.bed_format == "callable":
            # BED defines callable regions (new behavior)
            # Calculate callable length in this window
            callable_length = 0
            
            for _, row in bed_df.iterrows():
                chrom = row['chrom']
                bed_start = row['start']
                bed_end = row['end']
                
                # Find matching contig
                contig_idx = None
                for i, contig_name in enumerate(contig_names):
                    if str(contig_name) == str(chrom):
                        contig_idx = i
                        break
                
                if contig_idx is not None:
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
    
    def _print_statistic_means(self, results_df: pd.DataFrame, stats: List[str], analysis_type: str = "analysis"):
        """
        Calculate and print mean values for all statistics.
        
        Args:
            results_df: DataFrame with calculated statistics
            stats: List of statistic names to calculate means for
            analysis_type: Description of the analysis type for output
        """
        print(f"\nMean Statistics ({analysis_type}):")
        print("=" * 50)
        
        for stat in stats:
            if stat in results_df.columns:
                # Calculate mean, excluding NaN values
                mean_val = results_df[stat].mean()
                n_windows = len(results_df[stat].dropna())
                total_windows = len(results_df)
                
                if pd.isna(mean_val):
                    print(f"  {stat:20s}: No valid values")
                else:
                    print(f"  {stat:20s}: {mean_val:.6f} (n={n_windows}/{total_windows} windows)")
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
            stats = ['tajima_d', 'nucleotide_diversity', 'watterson_theta']
        
        if step_size is None:
            step_size = window_size
        
        all_results = []
        
        print(f"\nAnalyzing {len(regions)} regions with {window_size:,}bp windows...")
        
        for i, (contig, region_start, region_end) in enumerate(regions, 1):
            print(f"\n[{i}/{len(regions)}] Analyzing region {contig}:{region_start:,}-{region_end:,}")
            
            # Update window config for this region
            self.window_config.start = region_start
            self.window_config.end = region_end
            self.window_config.window_size = window_size
            self.window_config.step_size = step_size
            self.window_config.min_variants = min_variants
            
            try:
                # Create windows for this region
                self.create_windows()
                
                # Calculate statistics
                region_stats = self.calculate_windowed_stats(
                    stats=stats,
                    use_callable_sites=use_callable_sites
                )
                
                # Convert to DataFrame
                region_df = region_stats.to_dataframe()
                
                # Add region identifiers
                region_df.insert(0, 'region_contig', contig)
                region_df.insert(1, 'region_start', region_start)
                region_df.insert(2, 'region_end', region_end)
                
                all_results.append(region_df)
                
                print(f"  Analyzed {len(region_df)} windows")
                
            except Exception as e:
                print(f"  ERROR: analyzing region {contig}:{region_start}-{region_end}: {e}")
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
