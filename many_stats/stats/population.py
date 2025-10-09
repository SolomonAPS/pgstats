"""
Population genetics statistics implemented with Numba for high performance.

This module contains implementations of common population genetics statistics
such as Tajima's D, Fu and Li's statistics, and other neutrality tests.
"""

import numpy as np
import numba
import xarray as xr
import sgkit as sg
from typing import Tuple, Optional


@numba.njit
def calculate_a1(n: int) -> float:
    """
    Calculate a₁, the (n-1)th harmonic number: sum(1/i) for i from 1 to n-1.
    From Wakeley (2009) Coalescent Theory, used in equation 4.35 for Tajima's D
    and throughout chapter 4 for various neutrality statistics.
    
    Args:
        n: Sample size
        
    Returns:
        float: a₁ value
    """
    result = 0.0
    for i in range(1, n):
        result += 1.0 / i
    return result


@numba.njit
def calculate_a2(n: int) -> float:
    """
    Calculate a₂, sum of 1/i² for i from 1 to n-1.
    From Wakeley (2009) Coalescent Theory, used in variance calculations
    for neutrality statistics (equation 4.35 and related equations).
    
    Args:
        n: Sample size
        
    Returns:
        float: a₂ value
    """
    result = 0.0
    for i in range(1, n):
        result += 1.0 / (i * i)
    return result


@numba.njit
def calculate_b1(n: int) -> float:
    """
    Calculate b₁ = (n+1)/(3(n-1)).
    Used in variance calculations for neutrality statistics.
    
    Args:
        n: Sample size
        
    Returns:
        float: b₁ value
    """
    return (n + 1.0) / (3.0 * (n - 1.0))


@numba.njit
def calculate_b2(n: int) -> float:
    """
    Calculate b₂ = 2(n² + n + 3)/(9n(n-1)).
    Used in variance calculations for neutrality statistics.
    
    Args:
        n: Sample size
        
    Returns:
        float: b₂ value
    """
    return (2.0 * (n * n + n + 3.0)) / (9.0 * n * (n - 1.0))


@numba.njit
def calculate_c1(n: int, a1: float) -> float:
    """
    Calculate c₁ = b₁ - 1/a₁.
    Used in variance calculations for Tajima's D.
    
    Args:
        n: Sample size
        a1: a₁ value
        
    Returns:
        float: c₁ value
    """
    return calculate_b1(n) - 1.0 / a1


@numba.njit
def calculate_c2(n: int, a1: float, a2: float) -> float:
    """
    Calculate c₂ = b₂ - (n+2)/(a₁n) + a₂/a₁².
    Used in variance calculations for Tajima's D.
    
    Args:
        n: Sample size
        a1: a₁ value
        a2: a₂ value
        
    Returns:
        float: c₂ value
    """
    return calculate_b2(n) - (n + 2.0) / (a1 * n) + a2 / (a1 * a1)


@numba.njit
def get_unfolded_sfs(variant_matrix: np.ndarray) -> Tuple[np.ndarray, int]:
    """
    Calculate unfolded site frequency spectrum from variant matrix.
    From Wakeley (2009) Coalescent Theory, section 4.3.1:
    ξᵢ is the number of sites with i derived alleles.
    
    Args:
        variant_matrix: numpy array where rows are positions and columns are samples
                       0 = ancestral, 1 = derived
        
    Returns:
        tuple: (sfs, n) where:
            sfs: Array of counts [ξ₁, ξ₂, ..., ξₙ₋₁]
            n: Sample size
    """
    n_variants, n_samples = variant_matrix.shape
    sfs = np.zeros(n_samples - 1, dtype=np.int64)
    
    for i in range(n_variants):
        derived_count = 0
        for j in range(n_samples):
            if variant_matrix[i, j] == 1:
                derived_count += 1
        
        if 1 <= derived_count <= n_samples - 1:
            sfs[derived_count - 1] += 1
    
    return sfs, n_samples


@numba.njit
def get_folded_sfs(variant_matrix: np.ndarray) -> Tuple[np.ndarray, int]:
    """
    Calculate folded site frequency spectrum.
    From Wakeley (2009) Coalescent Theory, section 4.3.1:
    When ancestral state is unknown, we use η₁ which counts the minor allele frequency.
    
    Args:
        variant_matrix: numpy array where rows are positions and columns are samples
                       0 = ancestral, 1 = derived
        
    Returns:
        tuple: (folded_sfs, n) where:
            folded_sfs: Array of counts [η₁, η₂, ..., η_{n/2}]
            n: Sample size
    """
    n_variants, n_samples = variant_matrix.shape
    folded_sfs = np.zeros(n_samples // 2, dtype=np.int64)
    
    for i in range(n_variants):
        derived_count = 0
        for j in range(n_samples):
            if variant_matrix[i, j] == 1:
                derived_count += 1
        
        minor_count = min(derived_count, n_samples - derived_count)
        if 1 <= minor_count <= n_samples // 2:
            folded_sfs[minor_count - 1] += 1
    
    return folded_sfs, n_samples


@numba.njit
def calculate_pi(variant_matrix: np.ndarray) -> float:
    """
    Calculate π, average number of pairwise differences.
    From Wakeley (2009) Coalescent Theory, equation 4.39:
    π = 1/(n choose 2) * sum(i(n-i)ξᵢ) from i=1 to n-1
    
    Args:
        variant_matrix: numpy array where rows are positions and columns are samples
        
    Returns:
        float: π value
    """
    sfs, n = get_unfolded_sfs(variant_matrix)
    n_pairs = (n * (n - 1)) / 2.0
    
    pi = 0.0
    for i in range(1, n):
        pi += i * (n - i) * sfs[i - 1]
    
    return pi / n_pairs


@numba.njit
def calculate_S(sfs: np.ndarray) -> int:
    """
    Calculate S, the number of segregating sites.
    From Wakeley (2009) Coalescent Theory, equation 4.38:
    S = sum(ξᵢ) from i=1 to n-1
    
    Args:
        sfs: Site frequency spectrum (array of counts)
        
    Returns:
        int: Number of segregating sites
    """
    S = 0
    for i in range(len(sfs)):
        S += sfs[i]
    return S


def tajima_d(ds: xr.Dataset, call_genotype: str = "call_genotype") -> xr.Dataset:
    """
    Calculate Tajima's D statistic for each variant.
    
    Args:
        ds: sgkit Dataset containing genotype calls
        call_genotype: Name of the genotype variable
        
    Returns:
        Dataset with Tajima's D values
    """
    # Convert genotypes to variant matrix format
    genotypes = ds[call_genotype].values
    n_variants, n_samples, ploidy = genotypes.shape
    
    # Convert to binary matrix (0=ancestral, 1=derived)
    # For now, assume major allele is ancestral (this should be improved)
    variant_matrix = np.zeros((n_variants, n_samples), dtype=np.int8)
    
    for i in range(n_variants):
        for j in range(n_samples):
            # Sum alleles across ploidy
            allele_sum = np.sum(genotypes[i, j, :])
            # Convert to binary (0 or 1)
            variant_matrix[i, j] = int(allele_sum > ploidy // 2)
    
    # Calculate Tajima's D for each variant
    tajima_d_values = np.zeros(n_variants)
    
    for i in range(n_variants):
        # Get SFS for this variant
        sfs, n = get_unfolded_sfs(variant_matrix[i:i+1, :])
        S = calculate_S(sfs)
        
        if S == 0:
            tajima_d_values[i] = 0.0
            continue
            
        pi = calculate_pi(variant_matrix[i:i+1, :])
        a1 = calculate_a1(n)
        a2 = calculate_a2(n)
        
        # Calculate variance components
        c1 = calculate_c1(n, a1)
        c2 = calculate_c2(n, a1, a2)
        
        # Calculate variance
        var = c1 * S + c2 * S * (S - 1)
        
        if var <= 0:
            tajima_d_values[i] = 0.0
        else:
            tajima_d_values[i] = (pi - S / a1) / np.sqrt(var)
    
    # Create output dataset
    result = ds.copy()
    result["tajima_d"] = (["variants"], tajima_d_values)
    
    return result


def fu_li_d(ds: xr.Dataset, call_genotype: str = "call_genotype") -> xr.Dataset:
    """
    Calculate Fu and Li's D* statistic for each variant.
    
    Args:
        ds: sgkit Dataset containing genotype calls
        call_genotype: Name of the genotype variable
        
    Returns:
        Dataset with Fu and Li's D* values
    """
    # Convert genotypes to variant matrix format
    genotypes = ds[call_genotype].values
    n_variants, n_samples, ploidy = genotypes.shape
    
    # Convert to binary matrix
    variant_matrix = np.zeros((n_variants, n_samples), dtype=np.int8)
    
    for i in range(n_variants):
        for j in range(n_samples):
            allele_sum = np.sum(genotypes[i, j, :])
            variant_matrix[i, j] = int(allele_sum > ploidy // 2)
    
    # Calculate Fu and Li's D* for each variant
    fu_li_d_values = np.zeros(n_variants)
    
    for i in range(n_variants):
        # Get folded SFS for this variant
        folded_sfs, n = get_folded_sfs(variant_matrix[i:i+1, :])
        S = calculate_S(folded_sfs)
        
        if S == 0:
            fu_li_d_values[i] = 0.0
            continue
            
        eta_1 = folded_sfs[0] if len(folded_sfs) > 0 else 0
        a1 = calculate_a1(n)
        a2 = calculate_a2(n)
        
        # Calculate numerator
        numerator = S / a1 - ((n - 1) / n) * eta_1
        
        # Calculate variance components (simplified)
        u_d_star = a1 - 1 - (a1 * a1 + a2) / (a1 * a1 + a2)
        v_d_star = a2 - (a1 * a1 + a2) / (a1 * a1 + a2)
        
        # Calculate variance
        var = u_d_star * S + v_d_star * S * (S - 1)
        
        if var <= 0:
            if S == 1:
                fu_li_d_values[i] = numerator / np.sqrt(abs(u_d_star))
            else:
                fu_li_d_values[i] = 0.0
        else:
            fu_li_d_values[i] = numerator / np.sqrt(var)
    
    # Create output dataset
    result = ds.copy()
    result["fu_li_d"] = (["variants"], fu_li_d_values)
    
    return result


def fu_li_f(ds: xr.Dataset, call_genotype: str = "call_genotype") -> xr.Dataset:
    """
    Calculate Fu and Li's F* statistic for each variant.
    
    Args:
        ds: sgkit Dataset containing genotype calls
        call_genotype: Name of the genotype variable
        
    Returns:
        Dataset with Fu and Li's F* values
    """
    # Convert genotypes to variant matrix format
    genotypes = ds[call_genotype].values
    n_variants, n_samples, ploidy = genotypes.shape
    
    # Convert to binary matrix
    variant_matrix = np.zeros((n_variants, n_samples), dtype=np.int8)
    
    for i in range(n_variants):
        for j in range(n_samples):
            allele_sum = np.sum(genotypes[i, j, :])
            variant_matrix[i, j] = int(allele_sum > ploidy // 2)
    
    # Calculate Fu and Li's F* for each variant
    fu_li_f_values = np.zeros(n_variants)
    
    for i in range(n_variants):
        # Get folded SFS for this variant
        folded_sfs, n = get_folded_sfs(variant_matrix[i:i+1, :])
        S = calculate_S(folded_sfs)
        
        if S == 0:
            fu_li_f_values[i] = 0.0
            continue
            
        eta_1 = folded_sfs[0] if len(folded_sfs) > 0 else 0
        pi = calculate_pi(variant_matrix[i:i+1, :])
        a1 = calculate_a1(n)
        a2 = calculate_a2(n)
        
        # Calculate numerator
        numerator = pi - ((n - 1) / n) * eta_1
        
        # Calculate variance components (simplified)
        u_f_star = a1 + a2 - (a1 * a1 + a2) / (a1 * a1 + a2)
        v_f_star = a2 - (a1 * a1 + a2) / (a1 * a1 + a2)
        
        # Calculate variance
        var = u_f_star * S + v_f_star * S * (S - 1)
        
        if var <= 0:
            if S == 1:
                fu_li_f_values[i] = numerator / np.sqrt(abs(u_f_star))
            else:
                fu_li_f_values[i] = 0.0
        else:
            fu_li_f_values[i] = numerator / np.sqrt(var)
    
    # Create output dataset
    result = ds.copy()
    result["fu_li_f"] = (["variants"], fu_li_f_values)
    
    return result
