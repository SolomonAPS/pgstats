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
    Calculate D (coefficient of linkage disequilibrium) from phased data.
    
    For phased genotype data, we calculate D directly from gamete frequencies:
    D = P_AB - p_A * p_B
    
    where P_AB is the observed frequency of the AB gamete, and p_A, p_B are
    the allele frequencies at the two loci.
    
    Reference: Walsh and Lynch (2018) Equation 2.18
    
    Args:
        genotypes: Array of shape (n_samples, 2, ploidy) with phased genotype data
                  First dimension: samples
                  Second dimension: [locus_A, locus_B]
                  Third dimension: chromosomes (ploidy)
                  Values: 0 (reference allele), 1 (alternate allele), -1 (missing)
        
    Returns:
        D value
    """
    n_samples, n_loci, ploidy = genotypes.shape
    
    if n_loci != 2:
        return np.nan
    
    # Count gamete types: AB, Ab, aB, ab
    # Also count alleles for frequency calculation
    count_AB = 0.0  # Both alternate
    count_Ab = 0.0  # A alternate, B reference
    count_aB = 0.0  # A reference, B alternate
    count_ab = 0.0  # Both reference
    n_valid_gametes = 0
    n_A = 0.0  # Total alternate alleles at locus A
    n_B = 0.0  # Total alternate alleles at locus B
    n_valid_alleles = 0
    
    for sample_idx in range(n_samples):
        for chrom_idx in range(ploidy):
            allele_A = genotypes[sample_idx, 0, chrom_idx]
            allele_B = genotypes[sample_idx, 1, chrom_idx]
        
            # Skip if either allele is missing
            if allele_A < 0 or allele_B < 0:
                continue
            
            n_valid_gametes += 1
            
            # Count gamete types
            if allele_A == 1 and allele_B == 1:
                count_AB += 1
            elif allele_A == 1 and allele_B == 0:
                count_Ab += 1
            elif allele_A == 0 and allele_B == 1:
                count_aB += 1
            else:  # allele_A == 0 and allele_B == 0
                count_ab += 1
            
            # Count alleles for frequency calculation
            if allele_A == 1:
                n_A += 1
            if allele_B == 1:
                n_B += 1
            n_valid_alleles += 2
    
    if n_valid_gametes == 0:
        return np.nan
    
    # Calculate gamete frequencies
    P_AB = count_AB / n_valid_gametes
    
    # Calculate allele frequencies
    if n_valid_alleles == 0:
        return np.nan
    p_A = n_A / n_valid_alleles
    p_B = n_B / n_valid_alleles
    
    # Calculate D = P_AB - p_A * p_B
    D = P_AB - p_A * p_B
    
    return D


@numba.njit(nogil=True, fastmath=False)
def calculate_ld_d_prime(genotypes: np.ndarray) -> float:
    """
    Calculate D' (standardized D) from phased data.
    
    D' = D / D_max where D_max is the maximum possible D given allele frequencies.
    D_max depends on the sign of D:
    - If D >= 0: D_max = min(p_A * (1-p_B), (1-p_A) * p_B)
    - If D < 0: D_max = min(p_A * p_B, (1-p_A) * (1-p_B))
    
    This normalization bounds D' between -1 and 1.
    
    Args:
        genotypes: Array of shape (n_samples, 2, ploidy) with phased genotype data
        
    Returns:
        D' value
    """
    D = calculate_ld_d(genotypes)
    
    if np.isnan(D):
        return np.nan
    
    # Calculate allele frequencies from phased data
    n_samples, n_loci, ploidy = genotypes.shape
    n_A = 0.0
    n_B = 0.0
    n_valid_alleles = 0
    
    for sample_idx in range(n_samples):
        for chrom_idx in range(ploidy):
            allele_A = genotypes[sample_idx, 0, chrom_idx]
            allele_B = genotypes[sample_idx, 1, chrom_idx]
            
            if allele_A < 0 or allele_B < 0:
                continue
            
            if allele_A == 1:
                n_A += 1
            if allele_B == 1:
                n_B += 1
            n_valid_alleles += 2
    
    if n_valid_alleles == 0:
        return np.nan
    
    p_A = n_A / n_valid_alleles
    p_B = n_B / n_valid_alleles
    
    # Calculate D_max based on sign of D
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
    Calculate r² (squared correlation coefficient) from phased data.
    
    For phased data, we calculate r² directly from D and allele frequencies:
    r² = D² / (p_A(1-p_A) * p_B(1-p_B))
    
    This is the standardized measure of linkage disequilibrium, defined as the
    squared within-gamete correlation of allele frequencies at two loci.
    
    Reference: Walsh and Lynch (2018) Equation 2.22
    
    Args:
        genotypes: Array of shape (n_samples, 2, ploidy) with phased genotype data
                  First dimension: samples
                  Second dimension: [locus_A, locus_B]
                  Third dimension: chromosomes (ploidy)
                  Values: 0 (reference allele), 1 (alternate allele), -1 (missing)
        
    Returns:
        r² value
    """
    # Calculate D from phased data
    D = calculate_ld_d(genotypes)
    
    if np.isnan(D):
        return np.nan
    
    # Calculate allele frequencies from phased data
    n_samples, n_loci, ploidy = genotypes.shape
    n_A = 0.0
    n_B = 0.0
    n_valid_alleles = 0
    
    for sample_idx in range(n_samples):
        for chrom_idx in range(ploidy):
            allele_A = genotypes[sample_idx, 0, chrom_idx]
            allele_B = genotypes[sample_idx, 1, chrom_idx]
    
            if allele_A < 0 or allele_B < 0:
                continue
            
            if allele_A == 1:
                n_A += 1
            if allele_B == 1:
                n_B += 1
            n_valid_alleles += 2
    
    if n_valid_alleles == 0:
        return np.nan
    
    p_A = n_A / n_valid_alleles
    p_B = n_B / n_valid_alleles
    
    # Calculate denominator: p_A(1-p_A) * p_B(1-p_B)
    denominator = p_A * (1 - p_A) * p_B * (1 - p_B)
    
    if denominator <= 0:
        return np.nan
    
    # Calculate r² = D² / (p_A(1-p_A) * p_B(1-p_B))
    r_squared = (D * D) / denominator
    
    return r_squared


@numba.njit(nogil=True, fastmath=False)
def calculate_omega_statistic(genotypes: np.ndarray, window_size: int = 10) -> float:
    """
    Calculate Kim and Nielsen's omega (ω) statistic for detecting hard sweeps.
    
    The omega statistic contrasts linkage disequilibrium within regions
    versus between regions to detect the characteristic LD pattern of
    a completed hard sweep.
    
    Following pylibseq implementation: iterate through each SNP as a potential
    breakpoint, calculate LD within left region, within right region, and between
    regions, then find the maximum omega across all breakpoints.
    
    Reference: Walsh and Lynch (2018) Equation 9.37
    Kim & Nielsen (2004) Genetics 167:1513-1524
    
    Args:
        genotypes: Array of shape (n_variants, n_samples, ploidy) with genotype calls
        window_size: Minimum number of variants (unused, kept for compatibility)
        
    Returns:
        Maximum omega statistic value across all potential breakpoints
    """
    n_variants, n_samples, ploidy = genotypes.shape
    
    if n_variants < 3:  # Need at least 3 variants for meaningful analysis
        return np.nan
    
    # Convert to dosage (sum across ploidy)
    dosage = np.zeros((n_variants, n_samples))
    for i in range(n_variants):
        for j in range(n_samples):
            dosage[i, j] = genotypes[i, j, :].sum()
    
    # Calculate minor allele frequencies
    # Following Kim & Nielsen (2004) and pylibseq: only consider SNPs with MAF > 1
    mafs = np.zeros(n_variants)
    for i in range(n_variants):
        allele_count = dosage[i, :].sum()
        total_alleles = n_samples * ploidy
        maf = min(allele_count, total_alleles - allele_count)
        mafs[i] = maf
    
    # Calculate all pairwise r² values
    r_squared_matrix = np.full((n_variants, n_variants), np.nan)
    
    for i in range(n_variants):
        for j in range(i + 1, n_variants):
            # Extract phased genotype pair: shape (n_samples, 2, ploidy)
            # First dimension: samples
            # Second dimension: [variant_i, variant_j]
            # Third dimension: chromosomes (ploidy)
            pair_genotypes = np.zeros((n_samples, 2, ploidy), dtype=genotypes.dtype)
            pair_genotypes[:, 0, :] = genotypes[i, :, :]  # variant i
            pair_genotypes[:, 1, :] = genotypes[j, :, :]  # variant j
            
            # Calculate r²
            r_squared = calculate_ld_r_squared(pair_genotypes)
            r_squared_matrix[i, j] = r_squared
            r_squared_matrix[j, i] = r_squared
    
    # Iterate through each SNP as a potential breakpoint (excluding first and last)
    # Following Kim & Nielsen (2004) and pylibseq: only consider SNPs with MAF > 1
    max_omega = -np.inf
    
    for breakpoint in range(1, n_variants - 1):  # positions 1 to S-2
        if mafs[breakpoint] <= 1:
            continue
        
        # L = variants at positions [0, breakpoint] (inclusive)
        # R = variants at positions [breakpoint+1, n_variants-1] (inclusive)
        l = breakpoint + 1  # Number of variants in L (0-indexed, so add 1)
        s = n_variants  # Total number of variants
        
        if l < 1 or (s - l) < 1:
            continue
        
        # Calculate sum of r² within L
        sum_rsq_L = 0.0
        for i in range(0, breakpoint + 1):
            for j in range(i + 1, breakpoint + 1):
                if not np.isnan(r_squared_matrix[i, j]):
                    sum_rsq_L += r_squared_matrix[i, j]
        
        # Calculate sum of r² within R
        sum_rsq_R = 0.0
        for i in range(breakpoint + 1, n_variants):
            for j in range(i + 1, n_variants):
                if not np.isnan(r_squared_matrix[i, j]):
                    sum_rsq_R += r_squared_matrix[i, j]
        
        # Calculate sum of r² between L and R
        sum_rsq_LR = 0.0
        for i in range(0, breakpoint + 1):
            for j in range(breakpoint + 1, n_variants):
                if not np.isnan(r_squared_matrix[i, j]):
                    sum_rsq_LR += r_squared_matrix[i, j]
        
        # Calculate C_{S,ℓ} from Equation 9.37
        # C_{S,ℓ} = ℓ(S-ℓ) / [binom(ℓ,2) + binom(S-ℓ,2)]
        numerator_c = l * (s - l)
        denominator_c = (l * (l - 1) / 2.0) + ((s - l) * (s - l - 1) / 2.0)
        
        if denominator_c == 0:
            continue
        
        c_s_l = numerator_c / denominator_c
        
        # Calculate omega
        numerator_omega = c_s_l * (sum_rsq_L + sum_rsq_R)
        
        if sum_rsq_LR > 0:
            omega = numerator_omega / sum_rsq_LR
            if np.isfinite(omega):
                max_omega = max(max_omega, omega)
    
    return max_omega if max_omega != -np.inf else np.nan


@numba.njit(nogil=True, fastmath=False)
def calculate_kelly_zns(genotypes: np.ndarray) -> float:
    """
    Calculate Kelly's Z_nS statistic.
    
    Z_nS is the average of all pairwise r² values within a region:
    Z_nS = (2 / (S(S - 1))) * Σ (from i=1 to S-1) Σ (from j=i+1 to S) r²_ij
    
    where S is the number of segregating sites and r²_ij is the squared
    correlation coefficient between sites i and j.
    
    Reference: Walsh and Lynch (2018) Equation 9.36b
    Kelly (1997) Genetics 145:833-846
    
    Args:
        genotypes: Array of shape (n_variants, n_samples, ploidy) with genotype calls
        
    Returns:
        Z_nS value (average pairwise r²)
    """
    n_variants, n_samples, ploidy = genotypes.shape
    
    if n_variants < 2:  # Need at least 2 variants for pairwise comparison
        return np.nan
    
    # Calculate all pairwise r² values and compute mean directly
    sum_rsq = 0.0
    n_valid_pairs = 0
    
    for i in range(n_variants):
        for j in range(i + 1, n_variants):
            # Extract phased genotype pair: shape (n_samples, 2, ploidy)
            # First dimension: samples
            # Second dimension: [variant_i, variant_j]
            # Third dimension: chromosomes (ploidy)
            pair_genotypes = np.zeros((n_samples, 2, ploidy), dtype=genotypes.dtype)
            pair_genotypes[:, 0, :] = genotypes[i, :, :]  # variant i
            pair_genotypes[:, 1, :] = genotypes[j, :, :]  # variant j
            
            # Calculate r²
            r_squared = calculate_ld_r_squared(pair_genotypes)
            
            if not np.isnan(r_squared):
                sum_rsq += r_squared
                n_valid_pairs += 1
    
    if n_valid_pairs == 0:
        return np.nan
    
    # Z_nS is the mean of all pairwise r² values
    # This is equivalent to: (2 / (S(S-1))) * sum(r²_ij)
    # where the normalization factor accounts for the number of pairs
    zns = sum_rsq / n_valid_pairs
    
    return zns


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
            
            # Extract phased genotype pair: shape (n_samples, 2, ploidy)
            # First dimension: samples
            # Second dimension: [variant_i, variant_j]
            # Third dimension: chromosomes (ploidy)
            pair_genotypes = np.zeros((n_samples, 2, ploidy), dtype=genotypes.dtype)
            pair_genotypes[:, 0, :] = genotypes[i, :, :]  # variant i
            pair_genotypes[:, 1, :] = genotypes[j, :, :]  # variant j
            
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
    
    # Get variant positions and window information
    positions = ds.variant_position.values
    
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
        for i in range(min(5, n_windows)):
            start_pos = positions[window_starts[i]]
            stop_pos = positions[window_stops[i]-1]  # -1 because stop is exclusive
    
    # Initialize arrays for all windows
    n_variants_per_window = np.zeros(n_windows, dtype=np.int32)
    n_pairs_per_window = np.zeros(n_windows, dtype=np.int32)
    mean_D = np.full(n_windows, np.nan)
    mean_D_prime = np.full(n_windows, np.nan)
    mean_r_squared = np.full(n_windows, np.nan)
    max_r_squared = np.full(n_windows, np.nan)
    mean_distance = np.full(n_windows, np.nan)
    kelly_zns = np.full(n_windows, np.nan)
    
    if enable_profiling:
        import time
        start_time = time.time()
    
    for w_idx in range(n_windows):
        window_start_pos = window_starts[w_idx]
        window_end_pos = window_stops[w_idx]
        
        
        # Get variants in this window using indices
        window_variants = np.arange(window_starts[w_idx], window_stops[w_idx])
        
        # Store number of variants even if < 2
        n_variants_per_window[w_idx] = len(window_variants)
        
        if len(window_variants) < 2:
            continue
        
        # Calculate LD for all pairs in this window
        window_ld = []
        
        for i in range(len(window_variants)):
            for j in range(i + 1, len(window_variants)):
                var_i, var_j = window_variants[i], window_variants[j]
                
                # Extract phased genotype pair: shape (n_samples, 2, ploidy)
                # First dimension: samples
                # Second dimension: [variant_i, variant_j]
                # Third dimension: chromosomes (ploidy)
                pair_genotypes = np.zeros((n_samples, 2, ploidy), dtype=genotypes.dtype)
                pair_genotypes[:, 0, :] = genotypes[var_i, :, :]  # variant i
                pair_genotypes[:, 1, :] = genotypes[var_j, :, :]  # variant j
                
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
            kelly_zns[w_idx] = window_df['r_squared'].mean()
    
    # Create result dataset with same number of windows as input
    result = xr.Dataset()
    
    
    # Copy window coordinates from input
    
    result = result.assign_coords({
        'windows': ds.windows
    })
    
    # Add window information as data variables if they exist
    if 'window_contig' in ds.data_vars or 'window_contig' in ds.coords:
        result = result.assign({'window_contig': ds.window_contig})
    if 'window_start' in ds.data_vars or 'window_start' in ds.coords:
        result = result.assign({'window_start': ds.window_start})
    if 'window_stop' in ds.data_vars or 'window_stop' in ds.coords:
        result = result.assign({'window_stop': ds.window_stop})
    
    
    # Add data variables with WINDOWS dimension
    result = result.assign({
        'n_variants': (['windows'], n_variants_per_window),
        'n_pairs': (['windows'], n_pairs_per_window),
        'mean_D': (['windows'], mean_D),
        'mean_D_prime': (['windows'], mean_D_prime),
        'mean_r_squared': (['windows'], mean_r_squared),
        'max_r_squared': (['windows'], max_r_squared),
        'mean_distance': (['windows'], mean_distance),
        'kelly_zns': (['windows'], kelly_zns)
    })
    
    
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
        import time
        start_time = time.time()
    
    result = calculate_ld_matrix(ds, call_genotype)
    
    if enable_profiling:
        elapsed = time.time() - start_time
        print(f"[Timing] LD r² calculation completed in {elapsed:.2f}s", flush=True)
    
    return result


def omega_statistic(ds: xr.Dataset, 
                   call_genotype: str = "call_genotype",
                   omega_window_size: int = 10,
                   enable_profiling: bool = False) -> xr.Dataset:
    """
    Calculate Kim and Nielsen's omega (ω) statistic for detecting hard sweeps.
    
    The omega statistic contrasts linkage disequilibrium within regions
    versus between regions to detect the characteristic LD pattern of
    a completed hard sweep. This calculates omega for each window in the dataset.
    
    Reference: Walsh and Lynch (2018) Equation 9.37
    
    Args:
        ds: sgkit Dataset containing genotype calls and windows
        call_genotype: Name of the genotype variable
        omega_window_size: Size of the sliding window for omega analysis (number of variants)
        enable_profiling: Whether to print timing info
        
    Returns:
        Dataset with omega statistic values per window
    """
    if enable_profiling:
        import time
        start_time = time.time()
    
    genotypes = ds[call_genotype].values
    n_variants, n_samples, ploidy = genotypes.shape
    
    # Get window information (windows are always defined)
    n_windows = len(ds.windows)
    window_starts = ds.window_start_idx.values
    window_stops = ds.window_stop_idx.values
    
    # Calculate omega for each window
    omega_values = np.full(n_windows, np.nan)
    
    for w_idx in range(n_windows):
        window_start = window_starts[w_idx]
        window_stop = window_stops[w_idx]
        window_genotypes = genotypes[window_start:window_stop, :, :]
        
        # Calculate omega for this window
        if len(window_genotypes) >= 3:  # Need at least 3 variants
            omega_values[w_idx] = calculate_omega_statistic(window_genotypes, omega_window_size)
    
    # Create result dataset
    result = ds.copy()
    result["omega_statistic"] = (["windows"], omega_values)
    
    if enable_profiling:
        elapsed = time.time() - start_time
        print(f"[Timing] Omega statistic calculation completed in {elapsed:.2f}s", flush=True)
    
    return result


def kelly_zns(ds: xr.Dataset,
              call_genotype: str = "call_genotype",
              enable_profiling: bool = False) -> xr.Dataset:
    """
    Calculate Kelly's Z_nS statistic for detecting selective sweeps.
    
    Z_nS is the average of all pairwise r² values within a region. It is
    computed as:
    Z_nS = (2 / (S(S - 1))) * Σ (from i=1 to S-1) Σ (from j=i+1 to S) r²_ij
    
    where S is the number of segregating sites and r²_ij is the squared
    correlation coefficient between sites i and j.
    
    Kelly (1997) showed that Z_nS values are largely determined by the final
    coalescent time in the sample. A small Z_nS value is consistent with a hard
    sweep or extreme bottleneck, while a partial or soft sweep increases Z_nS.
    
    This function calculates Z_nS for each window in the dataset.
    
    Reference: Walsh and Lynch (2018) Equation 9.36b
    Kelly (1997) Genetics 145:833-846
    
    Args:
        ds: sgkit Dataset containing genotype calls and windows
        call_genotype: Name of the genotype variable
        enable_profiling: Whether to print timing info
        
    Returns:
        Dataset with Z_nS statistic values per window
    """
    if enable_profiling:
        import time
        start_time = time.time()
    
    genotypes = ds[call_genotype].values
    n_variants, n_samples, ploidy = genotypes.shape
    
    # Get window information (windows are always defined)
    n_windows = len(ds.windows)
    window_starts = ds.window_start_idx.values
    window_stops = ds.window_stop_idx.values
    
    # Calculate Z_nS for each window
    zns_values = np.full(n_windows, np.nan)
    
    for w_idx in range(n_windows):
        window_start = window_starts[w_idx]
        window_stop = window_stops[w_idx]
        window_genotypes = genotypes[window_start:window_stop, :, :]
        
        # Calculate Z_nS for this window
        if len(window_genotypes) >= 2:  # Need at least 2 variants for pairwise comparison
            zns_values[w_idx] = calculate_kelly_zns(window_genotypes)
    
    # Create result dataset
    result = ds.copy()
    result["kelly_zns"] = (["windows"], zns_values)
    
    if enable_profiling:
        elapsed = time.time() - start_time
        print(f"[Timing] Kelly Z_nS calculation completed in {elapsed:.2f}s", flush=True)
    
    return result


# =============================================================================
# LD DECAY
# =============================================================================

@numba.njit(nogil=True, fastmath=False)
def _compute_ld_decay_bins(genotypes: np.ndarray,
                           positions: np.ndarray,
                           max_distance: int,
                           bin_size: int) -> Tuple[np.ndarray, np.ndarray]:
    """
    Accumulate pairwise r-squared into distance bins for LD decay.

    Positions must be sorted.  The inner loop breaks early once distance
    exceeds max_distance.

    Bins are 1-indexed by bp: [1, bin_size], [bin_size+1, 2*bin_size], ...

    Returns:
        (sum_r2, counts) arrays of length n_bins.
    """
    n_variants = genotypes.shape[0]
    n_samples = genotypes.shape[1]
    ploidy = genotypes.shape[2]

    n_bins = (max_distance + bin_size - 1) // bin_size
    sum_r2 = np.zeros(n_bins, dtype=np.float64)
    counts = np.zeros(n_bins, dtype=np.int64)

    for i in range(n_variants):
        for j in range(i + 1, n_variants):
            dist = positions[j] - positions[i]
            if dist > max_distance:
                break
            if dist <= 0:
                continue

            bin_idx = (dist - 1) // bin_size
            if bin_idx >= n_bins:
                continue

            pair = np.zeros((n_samples, 2, ploidy), dtype=genotypes.dtype)
            pair[:, 0, :] = genotypes[i, :, :]
            pair[:, 1, :] = genotypes[j, :, :]

            r2 = calculate_ld_r_squared(pair)
            if not np.isnan(r2):
                sum_r2[bin_idx] += r2
                counts[bin_idx] += 1

    return sum_r2, counts


def ld_decay(ds: xr.Dataset,
             max_distance: int,
             bin_size: int,
             min_maf: float = 0.0,
             call_genotype: str = "call_genotype") -> xr.Dataset:
    """
    Compute LD decay: mean r-squared as a function of physical distance.

    All variant pairs with inter-SNP distance in (0, max_distance] are
    considered.  Pairs are placed in distance bins of width ``bin_size`` bp
    and mean r-squared is reported per bin.

    Args:
        ds: sgkit Dataset with call_genotype and variant_position.
        max_distance: maximum pairwise distance in bp.
        bin_size: distance bin width in bp.
        min_maf: minimum minor allele frequency (0-0.5).  Variants below
                 this threshold are excluded before computing pairs.
        call_genotype: name of the genotype variable.

    Returns:
        xr.Dataset with dimension ``distance_bins`` and variables:
        bin_start, bin_end, bin_midpoint, mean_r_squared, n_pairs.
    """
    genotypes = ds[call_genotype].values
    positions = ds["variant_position"].values
    n_variants, n_samples, ploidy = genotypes.shape

    if min_maf > 0:
        keep = np.ones(n_variants, dtype=np.bool_)
        for i in range(n_variants):
            valid_alleles = 0
            alt_count = 0
            for s in range(n_samples):
                if np.any(genotypes[i, s, :] == -1):
                    continue
                for p in range(ploidy):
                    alt_count += genotypes[i, s, p]
                valid_alleles += ploidy
            if valid_alleles == 0:
                keep[i] = False
            else:
                freq = alt_count / valid_alleles
                maf = min(freq, 1.0 - freq)
                if maf < min_maf:
                    keep[i] = False
        genotypes = genotypes[keep]
        positions = positions[keep]

    sum_r2, counts = _compute_ld_decay_bins(genotypes, positions,
                                            max_distance, bin_size)

    n_bins = len(sum_r2)
    bin_starts = np.arange(n_bins, dtype=np.int64) * bin_size + 1
    bin_ends = bin_starts + bin_size - 1
    bin_midpoints = (bin_starts + bin_ends) / 2.0
    with np.errstate(invalid='ignore'):
        mean_r2 = np.where(counts > 0, sum_r2 / counts, np.nan)

    return xr.Dataset(
        {
            "bin_start": (["distance_bins"], bin_starts),
            "bin_end": (["distance_bins"], bin_ends),
            "bin_midpoint": (["distance_bins"], bin_midpoints),
            "mean_r_squared": (["distance_bins"], mean_r2),
            "n_pairs": (["distance_bins"], counts),
        },
        coords={"distance_bins": np.arange(n_bins)},
    )
