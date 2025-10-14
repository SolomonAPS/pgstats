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
    # D_max = min(p_A * (1-p_B), (1-p_A) * p_B) if D > 0
    # D_max = min(p_A * p_B, (1-p_A) * (1-p_B)) if D < 0
    
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
                         window_size: int = 1000) -> xr.Dataset:
    """
    Calculate LD statistics within windows.
    
    Args:
        ds: sgkit Dataset containing genotype calls
        call_genotype: Name of the genotype variable
        window_size: Size of windows for LD calculation
        
    Returns:
        Dataset with windowed LD statistics
    """
    genotypes = ds[call_genotype].values
    n_variants, n_samples, ploidy = genotypes.shape
    
    # Convert to dosage format
    dosage = genotypes.sum(axis=2)
    
    # Get variant positions
    positions = ds.variant_position.values
    
    # Create windows
    min_pos = np.min(positions)
    max_pos = np.max(positions)
    
    window_results = []
    
    for window_start in range(min_pos, max_pos, window_size):
        window_end = window_start + window_size
        
        # Find variants in this window
        in_window = (positions >= window_start) & (positions < window_end)
        window_variants = np.where(in_window)[0]
        
        if len(window_variants) < 2:
            continue
        
        # Calculate LD for all pairs in this window
        window_ld = []
        
        for i in range(len(window_variants)):
            for j in range(i + 1, len(window_variants)):
                var_i, var_j = window_variants[i], window_variants[j]
                
                # Extract genotype pair
                pair_genotypes = np.column_stack([dosage[var_i, :], dosage[var_j, :]])
                
                # Calculate LD statistics
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
        
        if window_ld:
            # Calculate summary statistics for this window
            import pandas as pd
            window_df = pd.DataFrame(window_ld)
            
            window_results.append({
                'window_start': window_start,
                'window_end': window_end,
                'n_variants': len(window_variants),
                'n_pairs': len(window_ld),
                'mean_D': window_df['D'].mean(),
                'mean_D_prime': window_df['D_prime'].mean(),
                'mean_r_squared': window_df['r_squared'].mean(),
                'max_r_squared': window_df['r_squared'].max(),
                'mean_distance': window_df['distance'].mean()
            })
    
    if not window_results:
        # Return empty dataset
        return xr.Dataset({
            'window_start': (['windows'], []),
            'window_end': (['windows'], []),
            'n_variants': (['windows'], []),
            'n_pairs': (['windows'], []),
            'mean_D': (['windows'], []),
            'mean_D_prime': (['windows'], []),
            'mean_r_squared': (['windows'], []),
            'max_r_squared': (['windows'], []),
            'mean_distance': (['windows'], [])
        })
    
    import pandas as pd
    df = pd.DataFrame(window_results)
    
    # Convert to xarray Dataset
    result = xr.Dataset({
        'window_start': (['windows'], df['window_start'].values),
        'window_end': (['windows'], df['window_end'].values),
        'n_variants': (['windows'], df['n_variants'].values),
        'n_pairs': (['windows'], df['n_pairs'].values),
        'mean_D': (['windows'], df['mean_D'].values),
        'mean_D_prime': (['windows'], df['mean_D_prime'].values),
        'mean_r_squared': (['windows'], df['mean_r_squared'].values),
        'max_r_squared': (['windows'], df['max_r_squared'].values),
        'mean_distance': (['windows'], df['mean_distance'].values)
    })
    
    return result


# =============================================================================
# CONVENIENCE FUNCTIONS
# =============================================================================

def ld_d(ds: xr.Dataset, call_genotype: str = "call_genotype") -> xr.Dataset:
    """
    Calculate D (coefficient of linkage disequilibrium) for all variant pairs.
    
    Args:
        ds: sgkit Dataset containing genotype calls
        call_genotype: Name of the genotype variable
        
    Returns:
        Dataset with D values
    """
    return calculate_ld_matrix(ds, call_genotype)


def ld_d_prime(ds: xr.Dataset, call_genotype: str = "call_genotype") -> xr.Dataset:
    """
    Calculate D' (standardized D) for all variant pairs.
    
    Args:
        ds: sgkit Dataset containing genotype calls
        call_genotype: Name of the genotype variable
        
    Returns:
        Dataset with D' values
    """
    return calculate_ld_matrix(ds, call_genotype)


def ld_r_squared(ds: xr.Dataset, call_genotype: str = "call_genotype") -> xr.Dataset:
    """
    Calculate r² (squared correlation coefficient) for all variant pairs.
    
    Args:
        ds: sgkit Dataset containing genotype calls
        call_genotype: Name of the genotype variable
        
    Returns:
        Dataset with r² values
    """
    return calculate_ld_matrix(ds, call_genotype)
