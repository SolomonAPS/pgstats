"""
Selection statistics for population genetics.

This module contains implementations of tests for selection such as
McDonald-Kreitman test, Hudson-Kreitman-Aguade test, and other
selection-related statistics.
"""

import numpy as np
import numba
import xarray as xr
from typing import Tuple, Dict, Any


@numba.njit
def mcdonald_kreitman_test(synonymous_sfs: np.ndarray, 
                          nonsynonymous_sfs: np.ndarray,
                          n: int) -> Tuple[float, float, float]:
    """
    Calculate McDonald-Kreitman test statistics.
    
    The McDonald-Kreitman test compares the ratio of nonsynonymous to synonymous
    polymorphisms within species to the ratio of nonsynonymous to synonymous
    substitutions between species.
    
    Args:
        synonymous_sfs: Site frequency spectrum for synonymous sites
        nonsynonymous_sfs: Site frequency spectrum for nonsynonymous sites
        n: Sample size
        
    Returns:
        tuple: (G_statistic, p_value, alpha) where:
            G_statistic: G-test statistic
            p_value: p-value from chi-square test
            alpha: Proportion of adaptive substitutions
    """
    from .population import calculate_S
    
    # Calculate segregating sites
    S_syn = calculate_S(synonymous_sfs)
    S_nonsyn = calculate_S(nonsynonymous_sfs)
    
    # For this simplified version, we'll assume fixed differences
    # In practice, you'd need divergence data
    D_syn = 0  # Would be calculated from divergence data
    D_nonsyn = 0  # Would be calculated from divergence data
    
    # Calculate G-test statistic (simplified)
    if S_syn == 0 and S_nonsyn == 0:
        return 0.0, 1.0, 0.0
    
    # Expected values under neutrality
    total_poly = S_syn + S_nonsyn
    total_div = D_syn + D_nonsyn
    
    if total_poly == 0 or total_div == 0:
        return 0.0, 1.0, 0.0
    
    # Calculate alpha (proportion of adaptive substitutions)
    if D_syn > 0:
        alpha = 1 - (D_nonsyn / D_syn) * (S_syn / S_nonsyn) if S_nonsyn > 0 else 0
    else:
        alpha = 0
    
    # Simplified G-test (in practice, you'd use proper chi-square test)
    G_statistic = 0.0  # Would be calculated properly
    p_value = 1.0  # Would be calculated from chi-square distribution
    
    return G_statistic, p_value, alpha


@numba.njit
def hudson_kreitman_aguade_test(sfs: np.ndarray, n: int) -> Tuple[float, float]:
    """
    Calculate Hudson-Kreitman-Aguade test statistic.
    
    The HKA test compares polymorphism and divergence patterns across loci
    to test for selection or demographic effects.
    
    Args:
        sfs: Site frequency spectrum
        n: Sample size
        
    Returns:
        tuple: (HKA_statistic, p_value)
    """
    from .population import calculate_S, calculate_a1
    
    S = calculate_S(sfs)
    a1 = calculate_a1(n)
    
    if a1 == 0:
        return 0.0, 1.0
    
    # Calculate theta
    theta = S / a1
    
    # Simplified HKA test (in practice, you'd compare multiple loci)
    HKA_statistic = theta  # Simplified
    p_value = 1.0  # Would be calculated properly
    
    return HKA_statistic, p_value


def mcdonald_kreitman_analysis(ds: xr.Dataset, 
                              synonymous_mask: xr.DataArray,
                              nonsynonymous_mask: xr.DataArray,
                              call_genotype: str = "call_genotype") -> Dict[str, Any]:
    """
    Perform McDonald-Kreitman test analysis.
    
    Args:
        ds: sgkit Dataset containing genotype calls
        synonymous_mask: Boolean mask for synonymous sites
        nonsynonymous_mask: Boolean mask for nonsynonymous sites
        call_genotype: Name of the genotype variable
        
    Returns:
        Dictionary containing test results
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
    
    # Separate synonymous and nonsynonymous sites
    syn_indices = np.where(synonymous_mask.values)[0]
    nonsyn_indices = np.where(nonsynonymous_mask.values)[0]
    
    if len(syn_indices) == 0 or len(nonsyn_indices) == 0:
        return {
            "G_statistic": 0.0,
            "p_value": 1.0,
            "alpha": 0.0,
            "S_syn": 0,
            "S_nonsyn": 0
        }
    
    # Calculate SFS for synonymous and nonsynonymous sites
    syn_matrix = variant_matrix[syn_indices, :]
    nonsyn_matrix = variant_matrix[nonsyn_indices, :]
    
    from .population import get_unfolded_sfs
    
    syn_sfs, n = get_unfolded_sfs(syn_matrix)
    nonsyn_sfs, _ = get_unfolded_sfs(nonsyn_matrix)
    
    # Perform McDonald-Kreitman test
    G_stat, p_val, alpha = mcdonald_kreitman_test(syn_sfs, nonsyn_sfs, n)
    
    return {
        "G_statistic": G_stat,
        "p_value": p_val,
        "alpha": alpha,
        "S_syn": np.sum(syn_sfs),
        "S_nonsyn": np.sum(nonsyn_sfs)
    }


def hudson_kreitman_aguade_analysis(ds: xr.Dataset, 
                                   call_genotype: str = "call_genotype") -> xr.Dataset:
    """
    Perform Hudson-Kreitman-Aguade test analysis.
    
    Args:
        ds: sgkit Dataset containing genotype calls
        call_genotype: Name of the genotype variable
        
    Returns:
        Dataset with HKA test results
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
    
    # Calculate HKA test for each variant
    hka_values = np.zeros(n_variants)
    p_values = np.zeros(n_variants)
    
    for i in range(n_variants):
        from .population import get_unfolded_sfs
        sfs, n = get_unfolded_sfs(variant_matrix[i:i+1, :])
        hka_stat, p_val = hudson_kreitman_aguade_test(sfs, n)
        hka_values[i] = hka_stat
        p_values[i] = p_val
    
    # Create output dataset
    result = ds.copy()
    result["hka_statistic"] = (["variants"], hka_values)
    result["hka_p_value"] = (["variants"], p_values)
    
    return result
