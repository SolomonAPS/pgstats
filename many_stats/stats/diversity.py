"""
Diversity statistics for population genetics.

This module contains implementations of diversity measures such as
nucleotide diversity (π), Watterson's theta, and other diversity indices.
"""

import numpy as np
import numba
import xarray as xr
from typing import Tuple


@numba.njit
def calculate_theta_w(sfs: np.ndarray, n: int) -> float:
    """
    Calculate Watterson's theta (θw).
    From Wakeley (2009) Coalescent Theory, equation 4.40:
    θw = S/a₁ where S is the number of segregating sites and a₁ is the harmonic number.
    
    Args:
        sfs: Site frequency spectrum
        n: Sample size
        
    Returns:
        float: Watterson's theta
    """
    from .population import calculate_S, calculate_a1
    
    S = calculate_S(sfs)
    a1 = calculate_a1(n)
    
    if a1 == 0:
        return 0.0
    
    return S / a1


@numba.njit
def calculate_theta_h(sfs: np.ndarray, n: int) -> float:
    """
    Calculate Fay and Wu's theta (θh).
    From Fay and Wu (2000), equation 2:
    θh = sum(i²ξᵢ) / (n(n-1)) from i=1 to n-1
    
    Args:
        sfs: Site frequency spectrum
        n: Sample size
        
    Returns:
        float: Fay and Wu's theta
    """
    theta_h = 0.0
    for i in range(1, n):
        if i - 1 < len(sfs):
            theta_h += i * i * sfs[i - 1]
    
    return theta_h / (n * (n - 1))


@numba.njit
def calculate_theta_pi(variant_matrix: np.ndarray) -> float:
    """
    Calculate theta based on pairwise differences (π).
    This is the same as calculate_pi but renamed for clarity.
    
    Args:
        variant_matrix: Binary variant matrix
        
    Returns:
        float: θπ value
    """
    from .population import calculate_pi
    return calculate_pi(variant_matrix)


def nucleotide_diversity(ds: xr.Dataset, call_genotype: str = "call_genotype") -> xr.Dataset:
    """
    Calculate nucleotide diversity (π) for each variant.
    
    Args:
        ds: sgkit Dataset containing genotype calls
        call_genotype: Name of the genotype variable
        
    Returns:
        Dataset with nucleotide diversity values
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
    
    # Calculate π for each variant
    pi_values = np.zeros(n_variants)
    
    for i in range(n_variants):
        pi_values[i] = calculate_theta_pi(variant_matrix[i:i+1, :])
    
    # Create output dataset
    result = ds.copy()
    result["nucleotide_diversity"] = (["variants"], pi_values)
    
    return result


def watterson_theta(ds: xr.Dataset, call_genotype: str = "call_genotype") -> xr.Dataset:
    """
    Calculate Watterson's theta (θw) for each variant.
    
    Args:
        ds: sgkit Dataset containing genotype calls
        call_genotype: Name of the genotype variable
        
    Returns:
        Dataset with Watterson's theta values
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
    
    # Calculate θw for each variant
    theta_w_values = np.zeros(n_variants)
    
    for i in range(n_variants):
        sfs, n = get_unfolded_sfs(variant_matrix[i:i+1, :])
        theta_w_values[i] = calculate_theta_w(sfs, n)
    
    # Create output dataset
    result = ds.copy()
    result["watterson_theta"] = (["variants"], theta_w_values)
    
    return result


def fay_wu_theta(ds: xr.Dataset, call_genotype: str = "call_genotype") -> xr.Dataset:
    """
    Calculate Fay and Wu's theta (θh) for each variant.
    
    Args:
        ds: sgkit Dataset containing genotype calls
        call_genotype: Name of the genotype variable
        
    Returns:
        Dataset with Fay and Wu's theta values
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
    
    # Calculate θh for each variant
    theta_h_values = np.zeros(n_variants)
    
    for i in range(n_variants):
        sfs, n = get_unfolded_sfs(variant_matrix[i:i+1, :])
        theta_h_values[i] = calculate_theta_h(sfs, n)
    
    # Create output dataset
    result = ds.copy()
    result["fay_wu_theta"] = (["variants"], theta_h_values)
    
    return result


# Import the helper function
from .population import get_unfolded_sfs
