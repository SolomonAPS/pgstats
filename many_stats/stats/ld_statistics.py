"""
Linkage Disequilibrium Statistics

This module implements classic linkage disequilibrium statistics including:
- D (coefficient of linkage disequilibrium)
- D' (standardized D)
- r² (squared correlation coefficient)

Based on Walsh and Lynch (2018) and following computational best practices
from sgkit for handling missing data and performance.

References:
- Walsh and Lynch (2018) Evolution and Selection of Quantitative Traits
- Equations: D (9.1), D' (9.2), r² (9.3)
"""

import numpy as np
import numba
import xarray as xr
from typing import Tuple, Optional


# =============================================================================
# LOW-LEVEL LD CALCULATIONS (NUMBA-COMPILED)
# =============================================================================

@numba.njit(nogil=True, fastmath=False)
def calculate_allele_frequencies(genotypes: np.ndarray) -> Tuple[float, float, int]:
    """
    Calculate allele frequencies for two loci from genotype data.
    
    Args:
        genotypes: Array of shape (n_samples, 2) where each row is [locus1_genotype, locus2_genotype]
                  Values: 0 (homozygous reference), 1 (heterozygous), 2 (homozygous alternate)
                  -1 indicates missing data
    
    Returns:
        Tuple of (p_A, p_B, n_valid) where:
        - p_A: frequency of alternate allele at locus A
        - p_B: frequency of alternate allele at locus B  
        - n_valid: number of valid samples
    """
    n_samples = genotypes.shape[0]
    
    # Count alleles
    n_A = 0.0  # alternate alleles at locus A
    n_B = 0.0  # alternate alleles at locus B
    n_valid = 0
    
    for i in range(n_samples):
        g1, g2 = genotypes[i, 0], genotypes[i, 1]
        
        # Skip if either genotype is missing
        if g1 < 0 or g2 < 0:
            continue
            
        n_valid += 1
        
        # Count alleles (0=0 alleles, 1=1 allele, 2=2 alleles)
        n_A += g1
        n_B += g2
    
    if n_valid == 0:
        return np.nan, np.nan, 0
    
    # Convert counts to frequencies
    total_alleles = 2 * n_valid
    p_A = n_A / total_alleles
    p_B = n_B / total_alleles
    
    return p_A, p_B, n_valid


@numba.njit(nogil=True, fastmath=False)
def calculate_ld_d(genotypes: np.ndarray) -> float:
    """
    Calculate D (coefficient of linkage disequilibrium).
    
    For unphased genotype data, we use the correlation-based approach:
    D = r * sqrt(p_A * (1-p_A) * p_B * (1-p_B))
    
    where r is the correlation coefficient between the two loci.
    
    Reference: Walsh and Lynch (2018) Equation 9.1
    
    Args:
        genotypes: Array of shape (n_samples, 2) with genotype data
        
    Returns:
        D value
    """
    p_A, p_B, n_valid = calculate_allele_frequencies(genotypes)
    
    if n_valid < 2 or np.isnan(p_A) or np.isnan(p_B):
        return np.nan
    
    # Calculate correlation coefficient
    r = calculate_correlation_coefficient(genotypes)
    
    if np.isnan(r):
        return np.nan
    
    # D = r * sqrt(p_A * (1-p_A) * p_B * (1-p_B))
    D = r * np.sqrt(p_A * (1 - p_A) * p_B * (1 - p_B))
    
    return D


@numba.njit(nogil=True, fastmath=False)
def calculate_correlation_coefficient(genotypes: np.ndarray) -> float:
    """
    Calculate correlation coefficient between two loci.
    
    This is the classic Pearson correlation coefficient adapted for genotype data.
    
    Args:
        genotypes: Array of shape (n_samples, 2) with genotype data
        
    Returns:
        Correlation coefficient r
    """
    n_samples = genotypes.shape[0]
    
    # Calculate means and variances
    sum_x = sum_y = sum_x2 = sum_y2 = sum_xy = 0.0
    n_valid = 0
    
    for i in range(n_samples):
        x, y = genotypes[i, 0], genotypes[i, 1]
        
        # Skip missing data
        if x < 0 or y < 0:
            continue
            
        n_valid += 1
        sum_x += x
        sum_y += y
        sum_x2 += x * x
        sum_y2 += y * y
        sum_xy += x * y
    
    if n_valid < 2:
        return np.nan
    
    # Calculate means
    mean_x = sum_x / n_valid
    mean_y = sum_y / n_valid
    
    # Calculate covariance and variances
    cov_xy = (sum_xy / n_valid) - (mean_x * mean_y)
    var_x = (sum_x2 / n_valid) - (mean_x * mean_x)
    var_y = (sum_y2 / n_valid) - (mean_y * mean_y)
    
    # Calculate correlation coefficient
    if var_x <= 0 or var_y <= 0:
        return np.nan
    
    r = cov_xy / np.sqrt(var_x * var_y)
    
    return r


@numba.njit(nogil=True, fastmath=False)
def calculate_ld_d_prime(genotypes: np.ndarray) -> float:
    """
    Calculate D' (standardized D).
    
    D' = D / D_max where D_max is the maximum possible D given allele frequencies.
    
    Reference: Walsh and Lynch (2018) Equation 9.2
    
    Args:
        genotypes: Array of shape (n_samples, 2) with genotype data
        
    Returns:
        D' value
    """
    p_A, p_B, n_valid = calculate_allele_frequencies(genotypes)
    
    if n_valid < 2 or np.isnan(p_A) or np.isnan(p_B):
        return np.nan
    
    D = calculate_ld_d(genotypes)
    
    if np.isnan(D):
        return np.nan
    
    # Calculate D_max
    if D >= 0:
        D_max = min(p_A * (1 - p_B), (1 - p_A) * p_B)
    else:
        D_max = min(p_A * p_B, (1 - p_A) * (1 - p_B))
    
    if D_max <= 0:
        return np.nan
    
    D_prime = D / D_max
    
    return D_prime


@numba.njit(nogil=True, fastmath=False)
def calculate_ld_r_squared(genotypes: np.ndarray) -> float:
    """
    Calculate r² (squared correlation coefficient).
    
    This is the classic r² measure of linkage disequilibrium.
    
    Reference: Walsh and Lynch (2018) Equation 9.3
    
    Args:
        genotypes: Array of shape (n_samples, 2) with genotype data
        
    Returns:
        r² value
    """
    r = calculate_correlation_coefficient(genotypes)
    
    if np.isnan(r):
        return np.nan
    
    return r * r


@numba.njit(nogil=True, fastmath=False)
def calculate_omega_statistic(genotypes: np.ndarray, window_size: int = 10) -> float:
    """
    Calculate Kim and Nielsen's omega (ω) statistic for detecting hard sweeps.
    
    The omega statistic contrasts linkage disequilibrium within regions
    versus between regions to detect the characteristic LD pattern of
    a completed hard sweep.
    
    Reference: Walsh and Lynch (2018) Equation 9.37
    
    Args:
        genotypes: Array of shape (n_variants, n_samples, ploidy) with genotype calls
        window_size: Size of the sliding window for analysis
        
    Returns:
        Omega statistic value
    """
    n_variants, n_samples, ploidy = genotypes.shape
    
    if n_variants < 3:  # Need at least 3 variants for meaningful analysis
        return np.nan
    
    # Calculate all pairwise r² values
    r_squared_matrix = np.zeros((n_variants, n_variants))
    
    for i in range(n_variants):
        for j in range(i + 1, n_variants):
            # Get genotypes for this pair
            gt_i = genotypes[i, :, 0]  # Use first chromosome
            gt_j = genotypes[j, :, 0]  # Use first chromosome
            
            # Combine into pair format for calculate_ld_r_squared
            pair_genotypes = np.column_stack((gt_i, gt_j))
            
            # Calculate r²
            r_squared = calculate_ld_r_squared(pair_genotypes)
            r_squared_matrix[i, j] = r_squared
            r_squared_matrix[j, i] = r_squared
    
    # Calculate omega for different window positions
    max_omega = -np.inf
    
    for start_pos in range(n_variants - window_size + 1):
        end_pos = start_pos + window_size
        
        # Define regions L and R
        l_size = window_size // 2
        r_size = window_size - l_size
        
        if l_size < 1 or r_size < 1:
            continue
        
        # Calculate Cs,l
        s = window_size
        l = l_size
        c_s_l = (l * (s - l)) / (l * (l - 1) / 2 + (s - l) * (s - l - 1) / 2)
        
        # Calculate within-L LD
        within_l = 0.0
        for i in range(start_pos, start_pos + l):
            for j in range(i + 1, start_pos + l):
                within_l += r_squared_matrix[i, j]
        
        # Calculate within-R LD
        within_r = 0.0
        for i in range(start_pos + l, end_pos):
            for j in range(i + 1, end_pos):
                within_r += r_squared_matrix[i, j]
        
        # Calculate between L and R LD
        between_lr = 0.0
        for i in range(start_pos, start_pos + l):
            for j in range(start_pos + l, end_pos):
                between_lr += r_squared_matrix[i, j]
        
        # Calculate omega
        if between_lr > 0:
            omega = c_s_l * (within_l + within_r) / between_lr
            max_omega = max(max_omega, omega)
    
    return max_omega if max_omega != -np.inf else np.nan


# =============================================================================
# HIGH-LEVEL LD FUNCTIONS (USER-FACING)
# =============================================================================

def calculate_ld_matrix(ds: xr.Dataset, 
                       call_genotype: str = "call_genotype",
                       max_distance: Optional[int] = None) -> xr.Dataset:
    """
    Calculate LD matrix for all variant pairs in the dataset.
    
    Args:
        ds: sgkit Dataset containing genotype calls
        call_genotype: Name of the genotype variable
        max_distance: Maximum distance (in bp) to calculate LD for variant pairs
        
    Returns:
        Dataset with LD statistics for all variant pairs
    """
    genotypes = ds[call_genotype].values
    n_variants, n_samples, ploidy = genotypes.shape
    
    # Convert to dosage format (sum of alleles)
    dosage = genotypes.sum(axis=2)  # Shape: (n_variants, n_samples)
    
    # Get variant positions
    positions = ds.variant_position.values
    
    # Calculate LD for all pairs
    ld_results = []
    
    for i in range(n_variants):
        for j in range(i + 1, n_variants):
            # Check distance constraint
            if max_distance is not None:
                distance = abs(positions[j] - positions[i])
                if distance > max_distance:
                    continue
            
            # Extract genotype pair
            pair_genotypes = np.column_stack([dosage[i, :], dosage[j, :]])
            
            # Calculate LD statistics
            D = calculate_ld_d(pair_genotypes)
            D_prime = calculate_ld_d_prime(pair_genotypes)
            r_squared = calculate_ld_r_squared(pair_genotypes)
            
            # Only include if we have valid results
            if not (np.isnan(D) or np.isnan(D_prime) or np.isnan(r_squared)):
                ld_results.append({
                    'variant_i': i,
                    'variant_j': j,
                    'position_i': positions[i],
                    'position_j': positions[j],
                    'distance': abs(positions[j] - positions[i]),
                    'D': D,
                    'D_prime': D_prime,
                    'r_squared': r_squared
                })
    
    # Convert to DataFrame and then to Dataset
    if not ld_results:
        # Return empty dataset
        return xr.Dataset({
            'variant_i': (['pairs'], []),
            'variant_j': (['pairs'], []),
            'position_i': (['pairs'], []),
            'position_j': (['pairs'], []),
            'distance': (['pairs'], []),
            'D': (['pairs'], []),
            'D_prime': (['pairs'], []),
            'r_squared': (['pairs'], [])
        })
    
    import pandas as pd
    df = pd.DataFrame(ld_results)
    
    # Convert to xarray Dataset
    result = xr.Dataset({
        'variant_i': (['pairs'], df['variant_i'].values),
        'variant_j': (['pairs'], df['variant_j'].values),
        'position_i': (['pairs'], df['position_i'].values),
        'position_j': (['pairs'], df['position_j'].values),
        'distance': (['pairs'], df['distance'].values),
        'D': (['pairs'], df['D'].values),
        'D_prime': (['pairs'], df['D_prime'].values),
        'r_squared': (['pairs'], df['r_squared'].values)
    })
    
    return result


def calculate_windowed_ld(ds: xr.Dataset,
                         call_genotype: str = "call_genotype",
                         enable_profiling: bool = False) -> xr.Dataset:
    """
    Calculate LD statistics within existing windows.
    
    Args:
        ds: sgkit Dataset containing genotype calls and windows
        call_genotype: Name of the genotype variable
        enable_profiling: Whether to print timing and debug info
        
    Returns:
        Dataset with LD statistics for each window
    """
    genotypes = ds[call_genotype].values
    n_variants, n_samples, ploidy = genotypes.shape
    
    # Convert to dosage format
    dosage = genotypes.sum(axis=2)
    
    # Get variant positions and window information
    positions = ds.variant_position.values
    
    if enable_profiling:
        print(f"[DEBUG] Checking dataset...", flush=True)
        print(f"[DEBUG] Dataset dimensions: {dict(ds.sizes)}", flush=True)
        print(f"[DEBUG] Dataset coordinates: {list(ds.coords)}", flush=True)
        print(f"[DEBUG] Dataset data vars: {list(ds.data_vars)}", flush=True)
        print(f"[DEBUG] Number of variants: {len(positions)}", flush=True)
        print(f"[DEBUG] Position range: {np.min(positions)}-{np.max(positions)}", flush=True)
    
    # Verify windows exist
    required_vars = ['windows', 'window_start_idx', 'window_stop_idx']
    missing_vars = [var for var in required_vars if var not in ds.dims and var not in ds.data_vars]
    if missing_vars:
        raise ValueError(f"Input dataset missing required window variables: {missing_vars}. Call create_windows() first.")
    
    # Get window indices
    window_starts = ds.window_start_idx.values
    window_stops = ds.window_stop_idx.values
    n_windows = len(ds.windows)
    
    if enable_profiling:
        print(f"[DEBUG] Using {n_windows} windows", flush=True)
        print(f"[DEBUG] Window start indices: {window_starts[:5]}... (first 5)", flush=True)
        print(f"[DEBUG] Window stop indices: {window_stops[:5]}... (first 5)", flush=True)
        print(f"[DEBUG] Window positions:", flush=True)
        for i in range(min(5, n_windows)):
            start_pos = positions[window_starts[i]]
            stop_pos = positions[window_stops[i]-1]  # -1 because stop is exclusive
            print(f"[DEBUG]   Window {i}: {start_pos}-{stop_pos} (indices {window_starts[i]}-{window_stops[i]})", flush=True)
    
    # Initialize arrays for all windows
    n_variants_per_window = np.zeros(n_windows, dtype=np.int32)
    n_pairs_per_window = np.zeros(n_windows, dtype=np.int32)
    mean_D = np.full(n_windows, np.nan)
    mean_D_prime = np.full(n_windows, np.nan)
    mean_r_squared = np.full(n_windows, np.nan)
    max_r_squared = np.full(n_windows, np.nan)
    mean_distance = np.full(n_windows, np.nan)
    
    if enable_profiling:
        print(f"[DEBUG] Starting LD calculation for {n_windows} windows", flush=True)
        print(f"[DEBUG] Input dataset dimensions: {ds.dims}", flush=True)
        import time
        start_time = time.time()
    
    for w_idx in range(n_windows):
        window_start_pos = window_starts[w_idx]
        window_end_pos = window_stops[w_idx]
        
        if enable_profiling:
            print(f"[DEBUG] Window {w_idx}: start_pos={window_start_pos}, end_pos={window_end_pos}", flush=True)
        
        # Get variants in this window using indices
        window_variants = np.arange(window_starts[w_idx], window_stops[w_idx])
        
        if enable_profiling:
            print(f"[DEBUG] Window {w_idx}: found {len(window_variants)} variants", flush=True)
            if len(window_variants) > 0:
                print(f"[DEBUG] Window {w_idx}: variant indices {window_variants[0]}-{window_variants[-1]}", flush=True)
                print(f"[DEBUG] Window {w_idx}: variant positions {positions[window_variants[0]]}-{positions[window_variants[-1]]}", flush=True)
        
        # Store number of variants even if < 2
        n_variants_per_window[w_idx] = len(window_variants)
        
        if len(window_variants) < 2:
            continue
        
        # Calculate LD for all pairs in this window
        window_ld = []
        
        for i in range(len(window_variants)):
            for j in range(i + 1, len(window_variants)):
                var_i, var_j = window_variants[i], window_variants[j]
                
                # Extract genotype pair
                pair_genotypes = np.column_stack([dosage[var_i, :], dosage[var_j, :]])
                
                # Calculate LD statistics with error handling
                try:
                    D = calculate_ld_d(pair_genotypes)
                    D_prime = calculate_ld_d_prime(pair_genotypes)
                    r_squared = calculate_ld_r_squared(pair_genotypes)
                    
                    if not (np.isnan(D) or np.isnan(D_prime) or np.isnan(r_squared)):
                        window_ld.append({
                            'D': D,
                            'D_prime': D_prime,
                            'r_squared': r_squared,
                            'distance': abs(positions[var_j] - positions[var_i])
                        })
                except Exception as e:
                    if enable_profiling:
                        print(f"[WARNING] LD calculation failed for variants {var_i}-{var_j}: {str(e)}", flush=True)
                    continue
        
        if window_ld:
            # Calculate summary statistics for this window
            import pandas as pd
            window_df = pd.DataFrame(window_ld)
            
            # Store statistics for this window
            n_pairs_per_window[w_idx] = len(window_ld)
            mean_D[w_idx] = window_df['D'].mean()
            mean_D_prime[w_idx] = window_df['D_prime'].mean()
            mean_r_squared[w_idx] = window_df['r_squared'].mean()
            max_r_squared[w_idx] = window_df['r_squared'].max()
            mean_distance[w_idx] = window_df['distance'].mean()
    
    if enable_profiling:
        elapsed = time.time() - start_time
        print(f"[Timing] LD calculation completed in {elapsed:.2f}s", flush=True)
        n_valid = np.sum(~np.isnan(mean_D))
        print(f"[DEBUG] Processed {n_valid} windows with valid LD statistics", flush=True)
    
    if enable_profiling:
        print(f"[DEBUG] Creating LD result dataset...", flush=True)
        print(f"[DEBUG] Input dataset windows dimension: {ds.sizes.get('windows', 'not found')}", flush=True)
        print(f"[DEBUG] Number of windows to create: {n_windows}", flush=True)
        print(f"[DEBUG] Number of windows with valid LD stats: {np.sum(~np.isnan(mean_D))}", flush=True)
        print(f"[DEBUG] Window starts shape: {window_starts.shape}", flush=True)
        print(f"[DEBUG] Window stops shape: {window_stops.shape}", flush=True)
        print(f"[DEBUG] n_variants_per_window shape: {n_variants_per_window.shape}", flush=True)
        print(f"[DEBUG] mean_D shape: {mean_D.shape}", flush=True)
    
    # Create result dataset with same number of windows as input
    result = xr.Dataset()
    
    if enable_profiling:
        print(f"[DEBUG] Creating LD result dataset...", flush=True)
        print(f"[DEBUG] Input dataset:", flush=True)
        print(f"[DEBUG]   Dimensions: {dict(ds.sizes)}", flush=True)
        print(f"[DEBUG]   Coordinates: {list(ds.coords)}", flush=True)
        print(f"[DEBUG]   Data vars: {list(ds.data_vars)}", flush=True)
        print(f"[DEBUG] LD statistics:", flush=True)
        print(f"[DEBUG]   n_variants shape: {n_variants_per_window.shape}", flush=True)
        print(f"[DEBUG]   mean_D shape: {mean_D.shape}", flush=True)
        print(f"[DEBUG]   mean_D_prime shape: {mean_D_prime.shape}", flush=True)
        print(f"[DEBUG]   mean_r_squared shape: {mean_r_squared.shape}", flush=True)
        print(f"[DEBUG]   max_r_squared shape: {max_r_squared.shape}", flush=True)
        print(f"[DEBUG]   mean_distance shape: {mean_distance.shape}", flush=True)
    
    # Copy variant coordinates from input
    if enable_profiling:
        print(f"[DEBUG] Using variant coordinates from input dataset", flush=True)
        print(f"[DEBUG]   variants shape: {ds.variants.shape}", flush=True)
        print(f"[DEBUG]   variant_position shape: {ds.variant_position.shape}", flush=True)
    
    result = result.assign_coords({
        'variants': ds.variants,
        'variant_position': ds.variant_position
    })
    
    if enable_profiling:
        print(f"[DEBUG] Result dataset after coords:", flush=True)
        print(f"[DEBUG]   Dimensions: {dict(result.sizes)}", flush=True)
        print(f"[DEBUG]   Coordinates: {list(result.coords)}", flush=True)
    
    # Add data variables
    result = result.assign({
        'n_variants': (['variants'], n_variants_per_window),
        'n_pairs': (['variants'], n_pairs_per_window),
        'mean_D': (['variants'], mean_D),
        'mean_D_prime': (['variants'], mean_D_prime),
        'mean_r_squared': (['variants'], mean_r_squared),
        'max_r_squared': (['variants'], max_r_squared),
        'mean_distance': (['variants'], mean_distance)
    })
    
    if enable_profiling:
        print(f"[DEBUG] Final result dataset:", flush=True)
        print(f"[DEBUG]   Dimensions: {dict(result.sizes)}", flush=True)
        print(f"[DEBUG]   Coordinates: {list(result.coords)}", flush=True)
        print(f"[DEBUG]   Data vars: {list(result.data_vars)}", flush=True)
        print(f"[DEBUG]   n_variants range: {result.n_variants.values.min()}-{result.n_variants.values.max()}", flush=True)
        print(f"[DEBUG]   n_pairs range: {result.n_pairs.values.min()}-{result.n_pairs.values.max()}", flush=True)
        print(f"[DEBUG]   mean_D range: {np.nanmin(result.mean_D.values)}-{np.nanmax(result.mean_D.values)}", flush=True)
        print(f"[DEBUG]   mean_D_prime range: {np.nanmin(result.mean_D_prime.values)}-{np.nanmax(result.mean_D_prime.values)}", flush=True)
        print(f"[DEBUG]   mean_r_squared range: {np.nanmin(result.mean_r_squared.values)}-{np.nanmax(result.mean_r_squared.values)}", flush=True)
        print(f"[DEBUG]   max_r_squared range: {np.nanmin(result.max_r_squared.values)}-{np.nanmax(result.max_r_squared.values)}", flush=True)
        print(f"[DEBUG]   mean_distance range: {np.nanmin(result.mean_distance.values)}-{np.nanmax(result.mean_distance.values)}", flush=True)
    
    return result


# =============================================================================
# CONVENIENCE FUNCTIONS
# =============================================================================

def ld_d(ds: xr.Dataset, call_genotype: str = "call_genotype", enable_profiling: bool = False) -> xr.Dataset:
    """
    Calculate D (coefficient of linkage disequilibrium) for all variant pairs.
    
    Args:
        ds: sgkit Dataset containing genotype calls
        call_genotype: Name of the genotype variable
        enable_profiling: Whether to print timing and debug info
        
    Returns:
        Dataset with D values
    """
    if enable_profiling:
        print(f"[DEBUG] Computing LD D statistic...", flush=True)
        import time
        start_time = time.time()
    
    result = calculate_ld_matrix(ds, call_genotype)
    
    if enable_profiling:
        elapsed = time.time() - start_time
        print(f"[Timing] LD D calculation completed in {elapsed:.2f}s", flush=True)
    
    return result


def ld_d_prime(ds: xr.Dataset, call_genotype: str = "call_genotype", enable_profiling: bool = False) -> xr.Dataset:
    """
    Calculate D' (standardized D) for all variant pairs.
    
    Args:
        ds: sgkit Dataset containing genotype calls
        call_genotype: Name of the genotype variable
        enable_profiling: Whether to print timing and debug info
        
    Returns:
        Dataset with D' values
    """
    if enable_profiling:
        print(f"[DEBUG] Computing LD D' statistic...", flush=True)
        import time
        start_time = time.time()
    
    result = calculate_ld_matrix(ds, call_genotype)
    
    if enable_profiling:
        elapsed = time.time() - start_time
        print(f"[Timing] LD D' calculation completed in {elapsed:.2f}s", flush=True)
    
    return result


def ld_r_squared(ds: xr.Dataset, call_genotype: str = "call_genotype", enable_profiling: bool = False) -> xr.Dataset:
    """
    Calculate r² (squared correlation coefficient) for all variant pairs.
    
    Args:
        ds: sgkit Dataset containing genotype calls
        call_genotype: Name of the genotype variable
        enable_profiling: Whether to print timing and debug info
        
    Returns:
        Dataset with r² values
    """
    if enable_profiling:
        print(f"[DEBUG] Computing LD r² statistic...", flush=True)
        import time
        start_time = time.time()
    
    result = calculate_ld_matrix(ds, call_genotype)
    
    if enable_profiling:
        elapsed = time.time() - start_time
        print(f"[Timing] LD r² calculation completed in {elapsed:.2f}s", flush=True)
    
    return result


def omega_statistic(ds: xr.Dataset, 
                   call_genotype: str = "call_genotype",
                   window_size: int = 10,
                   enable_profiling: bool = False) -> xr.Dataset:
    """
    Calculate Kim and Nielsen's omega (ω) statistic for detecting hard sweeps.
    
    The omega statistic contrasts linkage disequilibrium within regions
    versus between regions to detect the characteristic LD pattern of
    a completed hard sweep.
    
    Reference: Walsh and Lynch (2018) Equation 9.37
    
    Args:
        ds: sgkit Dataset containing genotype calls
        call_genotype: Name of the genotype variable
        window_size: Size of the sliding window for analysis
        
    Returns:
        Dataset with omega statistic value
    """
    if enable_profiling:
        print(f"[DEBUG] Computing omega statistic...", flush=True)
        import time
        start_time = time.time()
    
    genotypes = ds[call_genotype].values
    
    # Calculate omega statistic
    omega_value = calculate_omega_statistic(genotypes, window_size)
    
    # Create result dataset
    result = ds.copy()
    result["omega_statistic"] = (["statistics"], np.array([omega_value]))
    
    if enable_profiling:
        elapsed = time.time() - start_time
        print(f"[Timing] Omega statistic calculation completed in {elapsed:.2f}s", flush=True)
    
    return result
