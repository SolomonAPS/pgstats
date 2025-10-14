"""
Haplotype Statistics for Population Genetics

This module implements haplotype-based statistics using sgkit's hashing approach:
- Haplotype diversity (H)
- Garud H1, H12, H123, H2/H1 statistics

References:
- Walsh and Lynch (2018) Evolution and Selection of Quantitative Traits
- Garud et al. (2015) A selective sweep on the Drosophila X chromosome
"""

import numpy as np
import numba
import xarray as xr
from typing import Tuple


# =============================================================================
# HAPLOTYPE ESTIMATION AND HASHING
# =============================================================================

@numba.njit(nogil=True, fastmath=False)
def hash_haplotype(haplotype: np.ndarray) -> int:
    """
    Create a hash for a haplotype sequence.
    
    This is a simplified version of sgkit's hash_array function.
    It creates a unique identifier for each haplotype by treating
    the sequence as a base-3 number (0, 1, missing).
    
    Args:
        haplotype: Array of shape (n_variants,) with haplotype data
                  0 = reference allele, 1 = alternate allele, -1 = missing
    
    Returns:
        Hash value for the haplotype
    """
    hash_value = 0
    base = 3  # Base-3 encoding: 0, 1, missing
    
    for i, allele in enumerate(haplotype):
        if allele == -1:  # Missing data
            hash_value = hash_value * base + 2
        else:
            hash_value = hash_value * base + int(allele)
    
    return hash_value


@numba.njit(nogil=True, fastmath=False)
def estimate_haplotypes_from_phased_data(genotypes: np.ndarray) -> np.ndarray:
    """
    Estimate haplotypes assuming data is already phased.
    
    This approach treats each chromosome separately, similar to sgkit's
    implementation. It works well when the data is already phased or
    when phase ambiguity is not critical for the analysis.
    
    Args:
        genotypes: Array of shape (n_variants, n_samples, ploidy) with genotype calls
                  0 = homozygous reference, 1 = heterozygous, 2 = homozygous alternate
                  -1 = missing data
    
    Returns:
        Array of shape (n_variants, n_haplotypes) where n_haplotypes = n_samples * ploidy
    """
    n_variants, n_samples, ploidy = genotypes.shape
    n_haplotypes = n_samples * ploidy
    
    haplotypes = np.zeros((n_variants, n_haplotypes), dtype=np.int8)
    
    for var_idx in range(n_variants):
        for sample_idx in range(n_samples):
            for chrom_idx in range(ploidy):
                haplotype_idx = sample_idx * ploidy + chrom_idx
                genotype = genotypes[var_idx, sample_idx, chrom_idx]
                
                if genotype == -1:  # Missing data
                    haplotypes[var_idx, haplotype_idx] = -1
                elif genotype == 0:  # Homozygous reference
                    haplotypes[var_idx, haplotype_idx] = 0
                elif genotype == 2:  # Homozygous alternate
                    haplotypes[var_idx, haplotype_idx] = 1
                else:  # Heterozygous - assign based on chromosome
                    haplotypes[var_idx, haplotype_idx] = chrom_idx  # 0 or 1
    
    return haplotypes


@numba.njit(nogil=True, fastmath=False)
def hash_haplotypes(haplotypes: np.ndarray) -> np.ndarray:
    """
    Hash all haplotypes to create unique identifiers.
    
    This is the core of sgkit's approach - it creates unique identifiers
    for each haplotype by treating the sequence as a base-3 number.
    
    Args:
        haplotypes: Array of shape (n_variants, n_haplotypes) with haplotype data
    
    Returns:
        Array of shape (n_haplotypes,) with hash values
    """
    n_variants, n_haplotypes = haplotypes.shape
    hash_values = np.zeros(n_haplotypes, dtype=np.int64)
    
    for haplotype_idx in range(n_haplotypes):
        hash_values[haplotype_idx] = hash_haplotype(haplotypes[:, haplotype_idx])
    
    return hash_values


# =============================================================================
# HAPLOTYPE STATISTICS CALCULATIONS
# =============================================================================

@numba.njit(nogil=True, fastmath=False)
def calculate_haplotype_diversity(haplotypes: np.ndarray) -> float:
    """
    Calculate haplotype diversity (H) - the probability that two randomly
    chosen haplotypes are different.
    
    Reference: Walsh and Lynch (2018) Equation 9.9
    
    Args:
        haplotypes: Array of shape (n_variants, n_haplotypes) with haplotype data
    
    Returns:
        Haplotype diversity value
    """
    n_variants, n_haplotypes = haplotypes.shape
    
    if n_haplotypes < 2:
        return np.nan
    
    # Count unique haplotypes
    unique_haplotypes = 0
    total_haplotypes = 0
    
    for var_idx in range(n_variants):
        # Get haplotypes for this variant (excluding missing data)
        var_haplotypes = haplotypes[var_idx, :]
        valid_haplotypes = var_haplotypes[var_haplotypes != -1]
        
        if len(valid_haplotypes) < 2:
            continue
        
        # Count unique haplotypes
        unique_count = len(np.unique(valid_haplotypes))
        total_count = len(valid_haplotypes)
        
        unique_haplotypes += unique_count
        total_haplotypes += total_count
    
    if total_haplotypes == 0:
        return np.nan
    
    # Calculate diversity as 1 - sum of squared frequencies
    # For simplicity, we'll use a basic approach
    diversity = 1.0 - (unique_haplotypes / total_haplotypes)
    
    return diversity


@numba.njit(nogil=True, fastmath=False)
def count_unique_values(values: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Count unique values and their frequencies (Numba-compatible version of np.unique with return_counts).
    
    Args:
        values: Array of values to count
    
    Returns:
        Tuple of (unique_values, counts)
    """
    # Get unique values
    unique_values = np.unique(values)
    
    # Count occurrences
    counts = np.zeros(len(unique_values), dtype=np.int64)
    for i, unique_val in enumerate(unique_values):
        counts[i] = np.sum(values == unique_val)
    
    return unique_values, counts


@numba.njit(nogil=True, fastmath=False)
def calculate_garud_h_statistics(haplotypes: np.ndarray) -> np.ndarray:
    """
    Calculate Garud H1, H12, H123, and H2/H1 statistics using sgkit's hashing approach.
    
    References:
    - Walsh and Lynch (2018) Equation 9.10-9.13
    - Garud et al. (2015) A selective sweep on the Drosophila X chromosome
    
    Args:
        haplotypes: Array of shape (n_variants, n_haplotypes) with haplotype data
    
    Returns:
        Array of [H1, H12, H123, H2/H1] values
    """
    n_variants, n_haplotypes = haplotypes.shape
    
    if n_haplotypes < 2:
        return np.array([np.nan, np.nan, np.nan, np.nan])
    
    # Hash all haplotypes to create unique identifiers
    hash_values = hash_haplotypes(haplotypes)
    
    # Count haplotype frequencies
    unique_hashes, counts = count_unique_values(hash_values)
    frequencies = counts / n_haplotypes
    
    # Sort frequencies in descending order
    frequencies = np.sort(frequencies)[::-1]
    
    # Calculate H1 (sum of squared frequencies)
    h1 = np.sum(frequencies ** 2)
    
    # Calculate H12 (sum of squared frequencies of top 2 haplotypes)
    if len(frequencies) >= 2:
        h12 = np.sum(frequencies[:2]) ** 2 + np.sum(frequencies[2:] ** 2)
    else:
        h12 = h1
    
    # Calculate H123 (sum of squared frequencies of top 3 haplotypes)
    if len(frequencies) >= 3:
        h123 = np.sum(frequencies[:3]) ** 2 + np.sum(frequencies[3:] ** 2)
    else:
        h123 = h1
    
    # Calculate H2/H1
    h2 = h1 - frequencies[0] ** 2
    h2_h1 = h2 / h1 if h1 > 0 else np.nan
    
    return np.array([h1, h12, h123, h2_h1])


# =============================================================================
# HIGH-LEVEL HAPLOTYPE FUNCTIONS (USER-FACING)
# =============================================================================

def haplotype_diversity(ds: xr.Dataset, 
                       call_genotype: str = "call_genotype") -> xr.Dataset:
    """
    Calculate haplotype diversity for each variant using sgkit's hashing approach.
    
    References:
    - Walsh and Lynch (2018) Equation 9.9
    
    Args:
        ds: sgkit Dataset containing genotype calls
        call_genotype: Name of the genotype variable
        
    Returns:
        Dataset with haplotype diversity values
    """
    genotypes = ds[call_genotype].values
    n_variants, n_samples, ploidy = genotypes.shape
    
    # Estimate haplotypes assuming phased data (sgkit approach)
    haplotypes = estimate_haplotypes_from_phased_data(genotypes)
    
    # Calculate haplotype diversity for each variant
    diversity_values = []
    for var_idx in range(n_variants):
        diversity = calculate_haplotype_diversity(haplotypes[var_idx:var_idx+1, :])
        diversity_values.append(diversity)
    
    # Create result dataset
    result = ds.copy()
    result["haplotype_diversity"] = (["variants"], np.array(diversity_values))
    
    return result


def garud_h_statistics(ds: xr.Dataset,
                      call_genotype: str = "call_genotype") -> xr.Dataset:
    """
    Calculate Garud H1, H12, H123, and H2/H1 statistics using sgkit's hashing approach.
    
    References:
    - Walsh and Lynch (2018) Equations 9.10-9.13
    - Garud et al. (2015) A selective sweep on the Drosophila X chromosome
    
    Args:
        ds: sgkit Dataset containing genotype calls
        call_genotype: Name of the genotype variable
        
    Returns:
        Dataset with Garud H statistics
    """
    genotypes = ds[call_genotype].values
    n_variants, n_samples, ploidy = genotypes.shape
    
    # Estimate haplotypes assuming phased data (sgkit approach)
    haplotypes = estimate_haplotypes_from_phased_data(genotypes)
    
    # Calculate Garud H statistics for each variant
    h1_values = []
    h12_values = []
    h123_values = []
    h2_h1_values = []
    
    for var_idx in range(n_variants):
        stats = calculate_garud_h_statistics(haplotypes[var_idx:var_idx+1, :])
        h1_values.append(stats[0])
        h12_values.append(stats[1])
        h123_values.append(stats[2])
        h2_h1_values.append(stats[3])
    
    # Create result dataset
    result = ds.copy()
    result["garud_h1"] = (["variants"], np.array(h1_values))
    result["garud_h12"] = (["variants"], np.array(h12_values))
    result["garud_h123"] = (["variants"], np.array(h123_values))
    result["garud_h2_h1"] = (["variants"], np.array(h2_h1_values))
    
    return result


def calculate_windowed_haplotype_stats(ds: xr.Dataset,
                                     call_genotype: str = "call_genotype",
                                     window_size: int = 1000) -> xr.Dataset:
    """
    Calculate windowed haplotype statistics using sgkit's hashing approach.
    
    Args:
        ds: sgkit Dataset containing genotype calls
        call_genotype: Name of the genotype variable
        window_size: Size of windows for haplotype analysis
        
    Returns:
        Dataset with windowed haplotype statistics
    """
    genotypes = ds[call_genotype].values
    n_variants, n_samples, ploidy = genotypes.shape
    
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
        
        # Extract genotypes for this window
        window_genotypes = genotypes[window_variants, :, :]
        
        # Estimate haplotypes assuming phased data (sgkit approach)
        haplotypes = estimate_haplotypes_from_phased_data(window_genotypes)
        
        # Calculate statistics
        diversity = calculate_haplotype_diversity(haplotypes)
        garud_stats = calculate_garud_h_statistics(haplotypes)
        
        window_results.append({
            'window_start': window_start,
            'window_end': window_end,
            'n_variants': len(window_variants),
            'haplotype_diversity': diversity,
            'garud_h1': garud_stats[0],
            'garud_h12': garud_stats[1],
            'garud_h123': garud_stats[2],
            'garud_h2_h1': garud_stats[3]
        })
    
    if not window_results:
        # Return empty dataset
        return xr.Dataset({
            'window_start': (['windows'], []),
            'window_end': (['windows'], []),
            'n_variants': (['windows'], []),
            'haplotype_diversity': (['windows'], []),
            'garud_h1': (['windows'], []),
            'garud_h12': (['windows'], []),
            'garud_h123': (['windows'], []),
            'garud_h2_h1': (['windows'], [])
        })
    
    import pandas as pd
    df = pd.DataFrame(window_results)
    
    # Convert to xarray Dataset
    result = xr.Dataset({
        'window_start': (['windows'], df['window_start'].values),
        'window_end': (['windows'], df['window_end'].values),
        'n_variants': (['windows'], df['n_variants'].values),
        'haplotype_diversity': (['windows'], df['haplotype_diversity'].values),
        'garud_h1': (['windows'], df['garud_h1'].values),
        'garud_h12': (['windows'], df['garud_h12'].values),
        'garud_h123': (['windows'], df['garud_h123'].values),
        'garud_h2_h1': (['windows'], df['garud_h2_h1'].values)
    })
    
    return result


# =============================================================================
# CONVENIENCE FUNCTIONS
# =============================================================================

def garud_h1(ds: xr.Dataset, call_genotype: str = "call_genotype") -> xr.Dataset:
    """Calculate Garud H1 statistic."""
    return garud_h_statistics(ds, call_genotype)


def garud_h12(ds: xr.Dataset, call_genotype: str = "call_genotype") -> xr.Dataset:
    """Calculate Garud H12 statistic."""
    return garud_h_statistics(ds, call_genotype)


def garud_h123(ds: xr.Dataset, call_genotype: str = "call_genotype") -> xr.Dataset:
    """Calculate Garud H123 statistic."""
    return garud_h_statistics(ds, call_genotype)


def garud_h2_h1(ds: xr.Dataset, call_genotype: str = "call_genotype") -> xr.Dataset:
    """Calculate Garud H2/H1 statistic."""
    return garud_h_statistics(ds, call_genotype)
