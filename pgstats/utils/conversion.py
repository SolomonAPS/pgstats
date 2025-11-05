"""
Data conversion utilities for genetic datasets.
"""

import numpy as np
import numba
import xarray as xr
from typing import Tuple


@numba.njit
def convert_genotypes_to_binary(genotypes: np.ndarray) -> np.ndarray:
    """
    Convert genotype calls to binary variant matrix.
    
    Args:
        genotypes: Array of shape (variants, samples, ploidy) with genotype calls
        
    Returns:
        Binary matrix of shape (variants, samples) where 0=ancestral, 1=derived
    """
    n_variants, n_samples, ploidy = genotypes.shape
    binary_matrix = np.zeros((n_variants, n_samples), dtype=np.int8)
    
    for i in range(n_variants):
        for j in range(n_samples):
            # Sum alleles across ploidy
            allele_sum = 0
            for k in range(ploidy):
                if genotypes[i, j, k] >= 0:  # Skip missing data
                    allele_sum += genotypes[i, j, k]
            
            # Convert to binary (0 or 1)
            binary_matrix[i, j] = int(allele_sum > ploidy // 2)
    
    return binary_matrix


@numba.njit
def convert_genotypes_to_allele_counts(genotypes: np.ndarray) -> np.ndarray:
    """
    Convert genotype calls to allele count matrix.
    
    Args:
        genotypes: Array of shape (variants, samples, ploidy) with genotype calls
        
    Returns:
        Allele count matrix of shape (variants, samples) with counts 0 to ploidy
    """
    n_variants, n_samples, ploidy = genotypes.shape
    count_matrix = np.zeros((n_variants, n_samples), dtype=np.int8)
    
    for i in range(n_variants):
        for j in range(n_samples):
            # Count non-reference alleles
            count = 0
            for k in range(ploidy):
                if genotypes[i, j, k] > 0:  # Count derived alleles
                    count += genotypes[i, j, k]
            
            count_matrix[i, j] = count
    
    return count_matrix


def convert_call_to_index(ds: xr.Dataset, call_genotype: str = "call_genotype") -> xr.Dataset:
    """
    Convert genotype calls to index format (similar to sgkit's convert_call_to_index).
    
    Args:
        ds: sgkit Dataset containing genotype calls
        call_genotype: Name of the genotype variable
        
    Returns:
        Dataset with converted genotype index
    """
    genotypes = ds[call_genotype].values
    n_variants, n_samples, ploidy = genotypes.shape
    
    # Convert to allele counts
    count_matrix = convert_genotypes_to_allele_counts(genotypes)
    
    # Create output dataset
    result = ds.copy()
    result["call_genotype_index"] = (["variants", "samples"], count_matrix)
    
    return result


def convert_to_variant_matrix(ds: xr.Dataset, call_genotype: str = "call_genotype") -> np.ndarray:
    """
    Convert sgkit dataset to binary variant matrix format.
    
    Args:
        ds: sgkit Dataset containing genotype calls
        call_genotype: Name of the genotype variable
        
    Returns:
        Binary variant matrix of shape (variants, samples)
    """
    genotypes = ds[call_genotype].values
    return convert_genotypes_to_binary(genotypes)
