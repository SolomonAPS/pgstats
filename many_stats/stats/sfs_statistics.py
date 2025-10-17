"""
Site Frequency Spectrum (SFS) based statistics for population genetics.

This module contains implementations of population genetics statistics that are
calculated from the site frequency spectrum (SFS), including:
- Theta estimators: theta_pi (π), theta_w (Watterson's), theta_h (Fay and Wu's), theta_l (Zeng's)
- Neutrality tests: Tajima's D, Fu and Li's D/D*, Fu and Li's F/F*, Zeng's E

All statistics properly handle missing data by calculating the SFS first.

References:
- Walsh and Lynch (2018) Evolution and Selection of Quantitative Traits
- Wakeley (2009) Coalescent Theory: An Introduction
- Fay and Wu (2000) Hitchhiking under positive Darwinian selection
- Tajima (1989) Statistical method for testing the neutral mutation hypothesis
- Fu and Li (1993) Statistical tests of neutrality of mutations
- Zeng et al. (2006) Statistical tests for detecting positive selection
"""

import numpy as np
import numba
import xarray as xr
import sgkit as sg
from typing import Tuple, Optional


# =============================================================================
# SITE FREQUENCY SPECTRUM CALCULATION
# =============================================================================

@numba.njit
def get_unfolded_sfs(variant_matrix: np.ndarray) -> Tuple[np.ndarray, int]:
    """
    Calculate unfolded site frequency spectrum from variant matrix.
    From Wakeley (2009) Coalescent Theory, section 4.3.1:
    ξᵢ is the number of sites with i derived alleles.
    
    Args:
        variant_matrix: numpy array where rows are positions and columns are samples
                       0 = ancestral, 1 = derived, -1 = missing
        
    Returns:
        tuple: (sfs, n) where:
            sfs: Array of counts [ξ₁, ξ₂, ..., ξₙ₋₁]
            n: Sample size (maximum possible, actual per-site n varies)
    """
    n_variants, n_samples = variant_matrix.shape
    sfs = np.zeros(n_samples - 1, dtype=np.int64)
    
    for i in range(n_variants):
        # Count non-missing samples and derived alleles for this site
        non_missing_count = 0
        derived_count = 0
        for j in range(n_samples):
            if variant_matrix[i, j] != -1:  # Not missing
                non_missing_count += 1
                if variant_matrix[i, j] == 1:
                    derived_count += 1
        
        # Only include if we have data and it's polymorphic
        if non_missing_count > 1 and 1 <= derived_count <= non_missing_count - 1:
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
                       0 = ancestral, 1 = derived, -1 = missing
        
    Returns:
        tuple: (folded_sfs, n) where:
            folded_sfs: Array of counts [η₁, η₂, ..., η_{n/2}]
            n: Sample size (maximum possible, actual per-site n varies)
    """
    n_variants, n_samples = variant_matrix.shape
    folded_sfs = np.zeros(n_samples // 2, dtype=np.int64)
    
    for i in range(n_variants):
        # Count non-missing samples and derived alleles for this site
        non_missing_count = 0
        derived_count = 0
        for j in range(n_samples):
            if variant_matrix[i, j] != -1:  # Not missing
                non_missing_count += 1
                if variant_matrix[i, j] == 1:
                    derived_count += 1
        
        # Calculate minor allele count
        minor_count = min(derived_count, non_missing_count - derived_count)
        if non_missing_count > 1 and 1 <= minor_count <= non_missing_count // 2:
            folded_sfs[minor_count - 1] += 1
    
    return folded_sfs, n_samples


# =============================================================================
# LOW-LEVEL HELPER FUNCTIONS (NUMBA-COMPILED)
# =============================================================================

@numba.njit
def get_per_site_sample_sizes(variant_matrix: np.ndarray) -> np.ndarray:
    """
    Get the actual sample size (non-missing count) for each site.
    
    Args:
        variant_matrix: numpy array where rows are positions and columns are samples
                       0 = ancestral, 1 = derived, -1 = missing
    
    Returns:
        Array of sample sizes per site
    """
    n_variants, n_samples = variant_matrix.shape
    site_n = np.zeros(n_variants, dtype=np.int64)
    
    for i in range(n_variants):
        count = 0
        for j in range(n_samples):
            if variant_matrix[i, j] != -1:
                count += 1
        site_n[i] = count
    
    return site_n

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
def calculate_c_n(n: int) -> float:
    """
    Calculate c_n = (n+1)/n - 1/a_n.
    Used in Fu and Li's D* variance calculation.
    From Fu and Li (1993), equation 9.26b.
    
    Args:
        n: Sample size
        
    Returns:
        float: c_n value
    """
    a1 = calculate_a1(n)
    return (n + 1.0) / n - 1.0 / a1


@numba.njit
def calculate_d_n(n: int) -> float:
    """
    Calculate d_n = b_n/a_n² - (2/n)(1 + 1/a_n + a_n/n) - 1/n².
    Used in Fu and Li's D* variance calculation.
    From Fu and Li (1993), equation 9.26c.
    
    Args:
        n: Sample size
        
    Returns:
        float: d_n value
    """
    a1 = calculate_a1(n)
    a2 = calculate_a2(n)
    b2 = calculate_b2(n)
    
    return b2 / (a1 * a1) - (2.0 / n) * (1.0 + 1.0 / a1 + a1 / n) - 1.0 / (n * n)


@numba.njit
def calculate_v_d_star(n: int, a1: float, a2: float) -> tuple[float, float]:
    """
    Calculate exact variance components for Fu and Li's D*.
    From Fu and Li (1993), equations 9.26b-c.
    
    Args:
        n: Sample size
        a1: a₁ value
        a2: a₂ value
        
    Returns:
        tuple: (α*, β*) variance components
    """
    b2 = calculate_b2(n)
    
    # β* = (1/(a_n² + b_n)) * [b_n/a_n² - (2/n)(1 + 1/a_n + a_n/n) - 1/n²]
    beta_star = (1.0 / (a1 * a1 + b2)) * calculate_d_n(n)
    
    # α* = (1/a_n) * [(n+1)/n - 1/a_n] - β*
    alpha_star = (1.0 / a1) * calculate_c_n(n) - beta_star
    
    return alpha_star, beta_star


@numba.njit
def calculate_v_f_star(n: int, a1: float, a2: float) -> tuple[float, float]:
    """
    Calculate exact variance components for Fu and Li's F*.
    From Fu and Li (1993), equations 9.26e-f.
    
    Args:
        n: Sample size
        a1: a₁ value
        a2: a₂ value
        
    Returns:
        tuple: (α_F, β_F) variance components
    """
    b2 = calculate_b2(n)
    
    # Calculate a_{n+1} for the formula
    a_n_plus_1 = calculate_a1(n + 1)
    
    # β_F = (1/(a_n² + b_n)) * [(2n³ + 110n² - 255n + 153)/(9n²(n-1)) + (2(n-1)a_n)/n² - (8b_n)/n]
    term1 = (2.0 * n * n * n + 110.0 * n * n - 255.0 * n + 153.0) / (9.0 * n * n * (n - 1.0))
    term2 = (2.0 * (n - 1.0) * a1) / (n * n)
    term3 = (8.0 * b2) / n
    
    beta_f = (1.0 / (a1 * a1 + b2)) * (term1 + term2 - term3)
    
    # α_F = (1/a_n) * [(4n² + 19n + 3 - 12(n+1)a_{n+1})/(3n(n-1))] - β_F
    numerator = 4.0 * n * n + 19.0 * n + 3.0 - 12.0 * (n + 1.0) * a_n_plus_1
    denominator = 3.0 * n * (n - 1.0)
    
    alpha_f = (1.0 / a1) * (numerator / denominator) - beta_f
    
    return alpha_f, beta_f


@numba.njit
def calculate_callable_sites_per_window(variant_matrix: np.ndarray) -> int:
    """
    Calculate the number of callable sites in a window.
    
    Callable sites are those that are not missing (not -1 in the genotype matrix).
    This is used to normalize theta estimators by the actual number of sites
    that can be analyzed, rather than the total window size.
    
    Args:
        variant_matrix: numpy array where rows are positions and columns are samples
        
    Returns:
        int: Number of callable sites
    """
    n_variants = variant_matrix.shape[0]
    callable_sites = 0
    
    for i in range(n_variants):
        # Check if this variant has any non-missing data
        has_data = False
        for j in range(variant_matrix.shape[1]):
            if variant_matrix[i, j] != -1:  # -1 indicates missing data
                has_data = True
                break
        
        if has_data:
            callable_sites += 1
    
    return callable_sites


# =============================================================================
# LOW-LEVEL THETA ESTIMATORS (NUMBA-COMPILED)
# =============================================================================

@numba.njit
def calculate_pi(variant_matrix: np.ndarray) -> float:
    """
    Calculate π (theta_pi), average number of pairwise differences.
    From Wakeley (2009) Coalescent Theory, equation 4.39:
    π = 1/(n choose 2) * sum(i(n-i)ξᵢ) from i=1 to n-1
    
    Args:
        variant_matrix: numpy array where rows are positions and columns are samples
                       0 = ancestral, 1 = derived, -1 = missing
        
    Returns:
        float: π value
    """
    n_variants, n_samples = variant_matrix.shape
    total_pi = 0.0
    
    for i in range(n_variants):
        # Count non-missing samples and derived alleles for this site
        non_missing = 0
        derived_count = 0
        for j in range(n_samples):
            if variant_matrix[i, j] != -1:
                non_missing += 1
                if variant_matrix[i, j] == 1:
                    derived_count += 1
        
        if non_missing > 1:
            n_pairs = (non_missing * (non_missing - 1)) / 2.0
            pi_site = (derived_count * (non_missing - derived_count)) / n_pairs
            total_pi += pi_site
    
    return total_pi


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
    # Calculate S (number of segregating sites)
    S = 0
    for i in range(len(sfs)):
        S += sfs[i]
    
    # Calculate a1 (harmonic number)
    a1 = 0.0
    for i in range(1, n):
        a1 += 1.0 / i
    
    if a1 == 0:
        return 0.0
    
    return S / a1


@numba.njit
def calculate_theta_w_per_site(variant_matrix: np.ndarray) -> float:
    """
    Calculate Watterson's theta (θw) accounting for per-site missing data.
    From Wakeley (2009) Coalescent Theory, equation 4.40:
    θw = S/a₁ where S is the number of segregating sites and a₁ is the harmonic number.
    
    Args:
        variant_matrix: numpy array where rows are positions and columns are samples
                       0 = ancestral, 1 = derived, -1 = missing
        
    Returns:
        float: Watterson's theta
    """
    n_variants, n_samples = variant_matrix.shape
    total_theta_w = 0.0
    
    for i in range(n_variants):
        # Count non-missing samples and derived alleles for this site
        non_missing = 0
        derived_count = 0
        for j in range(n_samples):
            if variant_matrix[i, j] != -1:
                non_missing += 1
                if variant_matrix[i, j] == 1:
                    derived_count += 1
        
        # Only count as segregating if polymorphic
        if non_missing > 1 and 1 <= derived_count <= non_missing - 1:
            a1 = calculate_a1(non_missing)
            theta_w_site = 1.0 / a1
            total_theta_w += theta_w_site
    
    return total_theta_w


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
def calculate_theta_h_per_site(variant_matrix: np.ndarray) -> float:
    """
    Calculate Fay and Wu's theta (θh) accounting for per-site missing data.
    From Fay and Wu (2000), equation 2:
    θh = sum(i²ξᵢ) / (n(n-1)) from i=1 to n-1
    
    Args:
        variant_matrix: numpy array where rows are positions and columns are samples
                       0 = ancestral, 1 = derived, -1 = missing
        
    Returns:
        float: Fay and Wu's theta
    """
    n_variants, n_samples = variant_matrix.shape
    total_theta_h = 0.0
    
    for i in range(n_variants):
        # Count non-missing samples and derived alleles for this site
        non_missing = 0
        derived_count = 0
        for j in range(n_samples):
            if variant_matrix[i, j] != -1:
                non_missing += 1
                if variant_matrix[i, j] == 1:
                    derived_count += 1
        
        if non_missing > 1 and 1 <= derived_count <= non_missing - 1:
            theta_h_site = (derived_count * derived_count) / (non_missing * (non_missing - 1))
            total_theta_h += theta_h_site
    
    return total_theta_h


@numba.njit
def calculate_theta_l(sfs: np.ndarray, n: int) -> float:
    """
    Calculate Zeng et al.'s theta_L (θL).
    From Zeng et al. (2006), equation 9.28a:
    θL = 1/(n-1) * sum(i·ξᵢ) from i=1 to n-1
    
    This estimator places more weight on high-frequency sites compared to θw,
    making it useful for detecting selective sweeps.
    
    Args:
        sfs: Site frequency spectrum
        n: Sample size
        
    Returns:
        float: Zeng et al.'s theta_L
    """
    theta_l = 0.0
    for i in range(1, n):
        if i - 1 < len(sfs):
            theta_l += i * sfs[i - 1]
    
    if n <= 1:
        return 0.0
    
    return theta_l / (n - 1.0)


@numba.njit
def calculate_theta_l_per_site(variant_matrix: np.ndarray) -> float:
    """
    Calculate Zeng et al.'s theta_L (θL) accounting for per-site missing data.
    From Zeng et al. (2006), equation 9.28a:
    θL = 1/(n-1) * sum(i·ξᵢ) from i=1 to n-1
    
    This estimator places more weight on high-frequency sites compared to θw,
    making it useful for detecting selective sweeps.
    
    Args:
        variant_matrix: numpy array where rows are positions and columns are samples
                       0 = ancestral, 1 = derived, -1 = missing
        
    Returns:
        float: Zeng et al.'s theta_L
    """
    n_variants, n_samples = variant_matrix.shape
    total_theta_l = 0.0
    
    for i in range(n_variants):
        # Count non-missing samples and derived alleles for this site
        non_missing = 0
        derived_count = 0
        for j in range(n_samples):
            if variant_matrix[i, j] != -1:
                non_missing += 1
                if variant_matrix[i, j] == 1:
                    derived_count += 1
        
        if non_missing > 1 and 1 <= derived_count <= non_missing - 1:
            theta_l_site = derived_count / (non_missing - 1.0)
            total_theta_l += theta_l_site
    
    return total_theta_l


# =============================================================================
# HIGH-LEVEL THETA ESTIMATORS (USER-FACING)
# =============================================================================

def theta_pi(ds: xr.Dataset, call_genotype: str = "call_genotype") -> xr.Dataset:
    """
    Calculate theta_pi (θπ) - nucleotide diversity based on pairwise differences.
    
    References:
    - Walsh and Lynch (2018) Equation 9.4
    - Wakeley (2009) Coalescent Theory, equation 4.39
    
    Formula: θπ = 1/(n choose 2) * sum(i(n-i)ξᵢ) from i=1 to n-1
    
    Args:
        ds: sgkit Dataset containing genotype calls
        call_genotype: Name of the genotype variable
        
    Returns:
        Dataset with theta_pi values
    """
    # Convert genotypes to variant matrix format
    genotypes = ds[call_genotype].values
    n_variants, n_samples, ploidy = genotypes.shape
    
    # Convert to binary matrix
    variant_matrix = np.zeros((n_variants, n_samples), dtype=np.int8)
    
    for i in range(n_variants):
        for j in range(n_samples):
            # Check for missing data first
            if np.any(genotypes[i, j, :] == -1):
                variant_matrix[i, j] = -1  # Mark as missing
            else:
                allele_sum = np.sum(genotypes[i, j, :])
                variant_matrix[i, j] = int(allele_sum > ploidy // 2)
    
    # Calculate θπ for each variant using calculate_pi
    pi_values = np.zeros(n_variants)
    
    for i in range(n_variants):
        pi_values[i] = calculate_pi(variant_matrix[i:i+1, :])
    
    # Create output dataset
    result = ds.copy()
    result["theta_pi"] = (["variants"], pi_values)
    
    return result


def theta_w(ds: xr.Dataset, call_genotype: str = "call_genotype") -> xr.Dataset:
    """
    Calculate theta_w (θw) - Watterson's estimator based on segregating sites.
    
    References:
    - Walsh and Lynch (2018) Equation 9.5
    - Wakeley (2009) Coalescent Theory, equation 4.40
    
    Formula: θw = S/a₁ where S is the number of segregating sites
    
    Args:
        ds: sgkit Dataset containing genotype calls
        call_genotype: Name of the genotype variable
        
    Returns:
        Dataset with theta_w values
    """
    # Convert genotypes to variant matrix format
    genotypes = ds[call_genotype].values
    n_variants, n_samples, ploidy = genotypes.shape
    
    # Convert to binary matrix
    variant_matrix = np.zeros((n_variants, n_samples), dtype=np.int8)
    
    for i in range(n_variants):
        for j in range(n_samples):
            # Check for missing data first
            if np.any(genotypes[i, j, :] == -1):
                variant_matrix[i, j] = -1  # Mark as missing
            else:
                allele_sum = np.sum(genotypes[i, j, :])
                variant_matrix[i, j] = int(allele_sum > ploidy // 2)
    
    # Calculate θw for each variant
    theta_w_values = np.zeros(n_variants)
    
    for i in range(n_variants):
        theta_w_values[i] = calculate_theta_w_per_site(variant_matrix[i:i+1, :])
    
    # Create output dataset
    result = ds.copy()
    result["theta_w"] = (["variants"], theta_w_values)
    
    return result


def theta_h(ds: xr.Dataset, call_genotype: str = "call_genotype") -> xr.Dataset:
    """
    Calculate theta_h (θh) - Fay and Wu's estimator based on high-frequency derived alleles.
    
    References:
    - Walsh and Lynch (2018) Equation 9.6
    - Fay and Wu (2000), equation 2
    
    Formula: θh = sum(i²ξᵢ) / (n(n-1))
    
    Args:
        ds: sgkit Dataset containing genotype calls
        call_genotype: Name of the genotype variable
        
    Returns:
        Dataset with theta_h values
    """
    # Convert genotypes to variant matrix format
    genotypes = ds[call_genotype].values
    n_variants, n_samples, ploidy = genotypes.shape
    
    # Convert to binary matrix
    variant_matrix = np.zeros((n_variants, n_samples), dtype=np.int8)
    
    for i in range(n_variants):
        for j in range(n_samples):
            # Check for missing data first
            if np.any(genotypes[i, j, :] == -1):
                variant_matrix[i, j] = -1  # Mark as missing
            else:
                allele_sum = np.sum(genotypes[i, j, :])
                variant_matrix[i, j] = int(allele_sum > ploidy // 2)
    
    # Calculate θh for each variant
    theta_h_values = np.zeros(n_variants)
    
    for i in range(n_variants):
        theta_h_values[i] = calculate_theta_h_per_site(variant_matrix[i:i+1, :])
    
    # Create output dataset
    result = ds.copy()
    result["theta_h"] = (["variants"], theta_h_values)
    
    return result


def theta_l(ds: xr.Dataset, call_genotype: str = "call_genotype") -> xr.Dataset:
    """
    Calculate theta_L (θL) - Zeng et al.'s estimator emphasizing high-frequency sites.
    
    References:
    - Walsh and Lynch (2018) Equation 9.7
    - Zeng et al. (2006), equation 9.28a
    
    Formula: θL = 1/(n-1) * sum(i·ξᵢ) from i=1 to n-1
    
    This estimator is useful for detecting selective sweeps as it places more
    weight on high-frequency alleles compared to Watterson's theta.
    
    Args:
        ds: sgkit Dataset containing genotype calls
        call_genotype: Name of the genotype variable
        
    Returns:
        Dataset with theta_L values
    """
    # Convert genotypes to variant matrix format
    genotypes = ds[call_genotype].values
    n_variants, n_samples, ploidy = genotypes.shape
    
    # Convert to binary matrix
    variant_matrix = np.zeros((n_variants, n_samples), dtype=np.int8)
    
    for i in range(n_variants):
        for j in range(n_samples):
            # Check for missing data first
            if np.any(genotypes[i, j, :] == -1):
                variant_matrix[i, j] = -1  # Mark as missing
            else:
                allele_sum = np.sum(genotypes[i, j, :])
                variant_matrix[i, j] = int(allele_sum > ploidy // 2)
    
    # Calculate θL for each variant
    theta_l_values = np.zeros(n_variants)
    
    for i in range(n_variants):
        theta_l_values[i] = calculate_theta_l_per_site(variant_matrix[i:i+1, :])
    
    # Create output dataset
    result = ds.copy()
    result["theta_l"] = (["variants"], theta_l_values)
    
    return result


# =============================================================================
# HIGH-LEVEL NEUTRALITY TESTS (USER-FACING)
# =============================================================================

def tajima_d(ds: xr.Dataset, call_genotype: str = "call_genotype") -> xr.Dataset:
    """
    Calculate Tajima's D statistic for each variant.
    
    References:
    - Walsh and Lynch (2018) Equation 9.8
    - Tajima (1989) and Wakeley (2009) equation 4.35
    
    Tajima's D tests for neutrality by comparing theta_pi and theta_w.
    
    Note: Following scikit-allel and sgkit conventions, harmonic numbers (a1, a2, etc.)
    are calculated using the maximum observed sample size across all variants.
    This maintains theoretical consistency with the assumption of constant n.
    
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
    variant_matrix = np.zeros((n_variants, n_samples), dtype=np.int8)
    
    for i in range(n_variants):
        for j in range(n_samples):
            # Check for missing data first
            if np.any(genotypes[i, j, :] == -1):
                variant_matrix[i, j] = -1  # Mark as missing
            else:
                # Sum alleles across ploidy
                allele_sum = np.sum(genotypes[i, j, :])
                # Convert to binary (0 or 1)
                variant_matrix[i, j] = int(allele_sum > ploidy // 2)
    
    # Calculate max sample size across all variants (following scikit-allel/sgkit)
    max_n = 0
    for i in range(n_variants):
        site_n = 0
        for j in range(n_samples):
            if variant_matrix[i, j] != -1:
                site_n += 1
        if site_n > max_n:
            max_n = site_n
    
    if max_n <= 1:
        # No valid data, return zeros
        result = ds.copy()
        result["tajima_d"] = (["variants"], np.zeros(n_variants))
        return result
    
    # Calculate harmonic numbers once using max_n
    a1 = calculate_a1(max_n)
    a2 = calculate_a2(max_n)
    c1 = calculate_c1(max_n, a1)
    c2 = calculate_c2(max_n, a1, a2)
    
    # Calculate Tajima's D for each variant
    tajima_d_values = np.zeros(n_variants)
    
    for i in range(n_variants):
        # Get SFS for this variant
        sfs, n_max = get_unfolded_sfs(variant_matrix[i:i+1, :])
        S = calculate_S(sfs)
        
        if S == 0:
            tajima_d_values[i] = 0.0
            continue
            
        pi = calculate_pi(variant_matrix[i:i+1, :])
        
        # Calculate variance using max_n harmonic numbers
        var = c1 * S + c2 * S * (S - 1)
        
        if var <= 0:
            tajima_d_values[i] = 0.0
        else:
            tajima_d_values[i] = (pi - S / a1) / np.sqrt(var)
    
    # Create output dataset
    result = ds.copy()
    result["tajima_d"] = (["variants"], tajima_d_values)
    
    return result


def fu_li_d(ds: xr.Dataset, call_genotype: str = "call_genotype", folded: bool = True) -> xr.Dataset:
    """
    Calculate Fu and Li's D or D* statistic for each variant.
    
    References:
    - Walsh and Lynch (2018) Equation 9.26b (D*) and 9.26c (D)
    - Fu and Li (1993) Statistical tests of neutrality of mutations
    
    Fu and Li's D uses the unfolded SFS (ζᵢ) while D* uses the folded SFS (ηᵢ).
    Both are sensitive to population size changes.
    
    Note: Following scikit-allel and sgkit conventions, harmonic numbers (a1, a2, etc.)
    are calculated using the maximum observed sample size across all variants.
    This maintains theoretical consistency with the assumption of constant n.
    
    Args:
        ds: sgkit Dataset containing genotype calls
        call_genotype: Name of the genotype variable
        folded: If True, use folded SFS (D*); if False, use unfolded SFS (D)
        
    Returns:
        Dataset with Fu and Li's D or D* values
    """
    # Convert genotypes to variant matrix format
    genotypes = ds[call_genotype].values
    n_variants, n_samples, ploidy = genotypes.shape
    
    # Convert to binary matrix
    variant_matrix = np.zeros((n_variants, n_samples), dtype=np.int8)
    
    for i in range(n_variants):
        for j in range(n_samples):
            # Check for missing data first
            if np.any(genotypes[i, j, :] == -1):
                variant_matrix[i, j] = -1  # Mark as missing
            else:
                allele_sum = np.sum(genotypes[i, j, :])
                variant_matrix[i, j] = int(allele_sum > ploidy // 2)
    
    # Calculate max sample size across all variants (following scikit-allel/sgkit)
    max_n = 0
    for i in range(n_variants):
        site_n = 0
        for j in range(n_samples):
            if variant_matrix[i, j] != -1:
                site_n += 1
        if site_n > max_n:
            max_n = site_n
    
    if max_n <= 1:
        # No valid data, return zeros
        result = ds.copy()
        stat_name = "fu_li_d" if not folded else "fu_li_d_star"
        result[stat_name] = (["variants"], np.zeros(n_variants))
        return result
    
    # Calculate harmonic numbers once using max_n
    a1 = calculate_a1(max_n)
    a2 = calculate_a2(max_n)
    u_d_star, v_d_star = calculate_v_d_star(max_n, a1, a2)
    
    # Calculate Fu and Li's D or D* for each variant
    fu_li_d_values = np.zeros(n_variants)
    
    for i in range(n_variants):
        if folded:
            # Use folded SFS (D*)
            sfs, n_max = get_folded_sfs(variant_matrix[i:i+1, :])
            singleton_count = sfs[0] if len(sfs) > 0 else 0
        else:
            # Use unfolded SFS (D)
            sfs, n_max = get_unfolded_sfs(variant_matrix[i:i+1, :])
            singleton_count = sfs[0] if len(sfs) > 0 else 0
            
        S = calculate_S(sfs)
        
        if S == 0:
            fu_li_d_values[i] = 0.0
            continue
        
        # Calculate numerator: S/a₁ - singleton estimator
        # Note: Using max_n for scaling factor to maintain consistency
        if folded:
            # D*: S/a₁ - ((n-1)/n)η₁
            numerator = S / a1 - ((max_n - 1) / max_n) * singleton_count
        else:
            # D: S/a₁ - ζ₁
            numerator = S / a1 - singleton_count
        
        # Calculate variance using max_n harmonic numbers
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
    stat_name = "fu_li_d" if not folded else "fu_li_d_star"
    result[stat_name] = (["variants"], fu_li_d_values)
    
    return result


def fu_li_f(ds: xr.Dataset, call_genotype: str = "call_genotype", folded: bool = True) -> xr.Dataset:
    """
    Calculate Fu and Li's F or F* statistic for each variant.
    
    References:
    - Walsh and Lynch (2018) Equation 9.26e (F*) and 9.26f (F)
    - Fu and Li (1993) Statistical tests of neutrality of mutations
    
    Fu and Li's F uses the unfolded SFS (ζᵢ) while F* uses the folded SFS (ηᵢ).
    Both combine pairwise differences with singleton information.
    
    Note: Following scikit-allel and sgkit conventions, harmonic numbers (a1, a2, etc.)
    are calculated using the maximum observed sample size across all variants.
    This maintains theoretical consistency with the assumption of constant n.
    
    Args:
        ds: sgkit Dataset containing genotype calls
        call_genotype: Name of the genotype variable
        folded: If True, use folded SFS (F*); if False, use unfolded SFS (F)
        
    Returns:
        Dataset with Fu and Li's F or F* values
    """
    # Convert genotypes to variant matrix format
    genotypes = ds[call_genotype].values
    n_variants, n_samples, ploidy = genotypes.shape
    
    # Convert to binary matrix
    variant_matrix = np.zeros((n_variants, n_samples), dtype=np.int8)
    
    for i in range(n_variants):
        for j in range(n_samples):
            # Check for missing data first
            if np.any(genotypes[i, j, :] == -1):
                variant_matrix[i, j] = -1  # Mark as missing
            else:
                allele_sum = np.sum(genotypes[i, j, :])
                variant_matrix[i, j] = int(allele_sum > ploidy // 2)
    
    # Calculate max sample size across all variants (following scikit-allel/sgkit)
    max_n = 0
    for i in range(n_variants):
        site_n = 0
        for j in range(n_samples):
            if variant_matrix[i, j] != -1:
                site_n += 1
        if site_n > max_n:
            max_n = site_n
    
    if max_n <= 1:
        # No valid data, return zeros
        result = ds.copy()
        stat_name = "fu_li_f" if not folded else "fu_li_f_star"
        result[stat_name] = (["variants"], np.zeros(n_variants))
        return result
    
    # Calculate harmonic numbers once using max_n
    a1 = calculate_a1(max_n)
    a2 = calculate_a2(max_n)
    u_f_star, v_f_star = calculate_v_f_star(max_n, a1, a2)
    
    # Calculate Fu and Li's F or F* for each variant
    fu_li_f_values = np.zeros(n_variants)
    
    for i in range(n_variants):
        if folded:
            # Use folded SFS (F*)
            sfs, n_max = get_folded_sfs(variant_matrix[i:i+1, :])
            singleton_count = sfs[0] if len(sfs) > 0 else 0
        else:
            # Use unfolded SFS (F)
            sfs, n_max = get_unfolded_sfs(variant_matrix[i:i+1, :])
            singleton_count = sfs[0] if len(sfs) > 0 else 0
            
        S = calculate_S(sfs)
        
        if S == 0:
            fu_li_f_values[i] = 0.0
            continue
            
        pi = calculate_pi(variant_matrix[i:i+1, :])
        
        # Calculate numerator: π - singleton estimator
        # Note: Using max_n for scaling factor to maintain consistency
        if folded:
            # F*: π - ((n-1)/n)η₁
            numerator = pi - ((max_n - 1) / max_n) * singleton_count
        else:
            # F: π - ζ₁
            numerator = pi - singleton_count
        
        # Calculate variance using max_n harmonic numbers
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
    stat_name = "fu_li_f" if not folded else "fu_li_f_star"
    result[stat_name] = (["variants"], fu_li_f_values)
    
    return result


# Convenience functions for unfolded versions
def fu_li_d_unfolded(ds: xr.Dataset, call_genotype: str = "call_genotype") -> xr.Dataset:
    """
    Calculate Fu and Li's D statistic using unfolded SFS.
    
    This is a convenience function that calls fu_li_d with folded=False.
    
    Args:
        ds: sgkit Dataset containing genotype calls
        call_genotype: Name of the genotype variable
        
    Returns:
        Dataset with Fu and Li's D values
    """
    return fu_li_d(ds, call_genotype, folded=False)


def fu_li_f_unfolded(ds: xr.Dataset, call_genotype: str = "call_genotype") -> xr.Dataset:
    """
    Calculate Fu and Li's F statistic using unfolded SFS.
    
    This is a convenience function that calls fu_li_f with folded=False.
    
    Args:
        ds: sgkit Dataset containing genotype calls
        call_genotype: Name of the genotype variable
        
    Returns:
        Dataset with Fu and Li's F values
    """
    return fu_li_f(ds, call_genotype, folded=False)


def zeng_e(ds: xr.Dataset, call_genotype: str = "call_genotype") -> xr.Dataset:
    """
    Calculate Zeng et al.'s E statistic for detecting selective sweeps.
    
    References:
    - Walsh and Lynch (2018) Equation 9.28c
    - Zeng et al. (2006), equation 9.28b
    
    Formula: E = (θ_L - θ_S) / σ(E)
    
    Where θ_S is Watterson's theta (theta_w) and θ_L emphasizes high-frequency sites.
    The E test is powerful for detecting selective sweeps and can persist longer
    after a sweep (up to 2N generations) compared to other tests.
    
    A negative E indicates an excess of low-frequency sites, which occurs
    immediately after a selective sweep.
    
    Note: Following scikit-allel and sgkit conventions, harmonic numbers (a1, a2, etc.)
    are calculated using the maximum observed sample size across all variants.
    This maintains theoretical consistency with the assumption of constant n.
    
    Args:
        ds: sgkit Dataset containing genotype calls
        call_genotype: Name of the genotype variable
        
    Returns:
        Dataset with Zeng's E values
    """
    # Convert genotypes to variant matrix format
    genotypes = ds[call_genotype].values
    n_variants, n_samples, ploidy = genotypes.shape
    
    # Convert to binary matrix
    variant_matrix = np.zeros((n_variants, n_samples), dtype=np.int8)
    
    for i in range(n_variants):
        for j in range(n_samples):
            # Check for missing data first
            if np.any(genotypes[i, j, :] == -1):
                variant_matrix[i, j] = -1  # Mark as missing
            else:
                allele_sum = np.sum(genotypes[i, j, :])
                variant_matrix[i, j] = int(allele_sum > ploidy // 2)
    
    # Calculate max sample size across all variants (following scikit-allel/sgkit)
    max_n = 0
    for i in range(n_variants):
        site_n = 0
        for j in range(n_samples):
            if variant_matrix[i, j] != -1:
                site_n += 1
        if site_n > max_n:
            max_n = site_n
    
    if max_n <= 1:
        # No valid data, return zeros
        result = ds.copy()
        result["zeng_e"] = (["variants"], np.zeros(n_variants))
        return result
    
    # Calculate harmonic numbers once using max_n
    a1 = calculate_a1(max_n)
    a2 = calculate_a2(max_n)
    b1 = calculate_b1(max_n)
    b2 = calculate_b2(max_n)
    
    # Pre-calculate variance components using max_n
    term1 = (max_n / (2.0 * (max_n - 1.0))) - (1.0 / a1)
    term2 = b2 + 2.0 * (max_n / (max_n - 1.0)) ** 2 * b2
    term2 -= 2.0 * (max_n * b2 - max_n + 1.0) / ((max_n - 1.0) * a1)
    term2 -= (3.0 * max_n + 1.0) / (max_n - 1.0)
    
    # Calculate Zeng's E for each variant
    zeng_e_values = np.zeros(n_variants)
    
    for i in range(n_variants):
        # Get unfolded SFS for this variant
        sfs, n_max = get_unfolded_sfs(variant_matrix[i:i+1, :])
        S = calculate_S(sfs)
        
        if S == 0:
            zeng_e_values[i] = 0.0
            continue
        
        # Calculate theta_L and theta_w (theta_S)
        theta_l_val = calculate_theta_l_per_site(variant_matrix[i:i+1, :])
        theta_w_val = calculate_theta_w_per_site(variant_matrix[i:i+1, :])
        
        # Calculate variance of E according to equation 9.28c
        # σ²(E) = [n/(2(n-1)) - 1/a_n]θ + [b_n + 2(n/(n-1))² * b_n - 2(nbn-n+1)/((n-1)a_n) - 3n+1/(n-1)]θ²
        # Using max_n for variance components and theta_w as estimate
        
        var_e = term1 * theta_w_val + term2 * theta_w_val * theta_w_val
        
        if var_e <= 0:
            if S == 1:
                zeng_e_values[i] = (theta_l_val - theta_w_val) / np.sqrt(abs(term1))
            else:
                zeng_e_values[i] = 0.0
        else:
            zeng_e_values[i] = (theta_l_val - theta_w_val) / np.sqrt(var_e)
    
    # Create output dataset
    result = ds.copy()
    result["zeng_e"] = (["variants"], zeng_e_values)
    
    return result


def singletons(ds: xr.Dataset, call_genotype: str = "call_genotype", folded: bool = True) -> xr.Dataset:
    """
    Calculate singleton count for each variant.
    
    References:
    - Walsh and Lynch (2018) Equation 9.26b (folded) and 9.26c (unfolded)
    - Fu and Li (1993) Statistical tests of neutrality of mutations
    
    Singletons are sites where an allele appears only once in the sample.
    
    Args:
        ds: sgkit Dataset containing genotype calls
        call_genotype: Name of the genotype variable
        folded: If True, count minor allele singletons (η₁); 
                if False, count derived allele singletons (ζ₁)
        
    Returns:
        Dataset with singleton counts
    """
    # Convert genotypes to variant matrix format
    genotypes = ds[call_genotype].values
    n_variants, n_samples, ploidy = genotypes.shape
    
    # Convert to binary matrix
    variant_matrix = np.zeros((n_variants, n_samples), dtype=np.int8)
    
    for i in range(n_variants):
        for j in range(n_samples):
            # Check for missing data first
            if np.any(genotypes[i, j, :] == -1):
                variant_matrix[i, j] = -1  # Mark as missing
            else:
                allele_sum = np.sum(genotypes[i, j, :])
                variant_matrix[i, j] = int(allele_sum > ploidy // 2)

    # Calculate singletons for each variant
    singleton_values = np.zeros(n_variants)
    
    for i in range(n_variants):
        if folded:
            # Use folded SFS (η₁)
            sfs, n = get_folded_sfs(variant_matrix[i:i+1, :])
            singleton_count = sfs[0] if len(sfs) > 0 else 0
        else:
            # Use unfolded SFS (ζ₁)
            sfs, n = get_unfolded_sfs(variant_matrix[i:i+1, :])
            singleton_count = sfs[0] if len(sfs) > 0 else 0
            
        singleton_values[i] = singleton_count
    
    # Create output dataset
    result = ds.copy()
    if folded:
        result["singletons_folded"] = (["variants"], singleton_values)
    else:
        result["singletons_unfolded"] = (["variants"], singleton_values)
    
    return result

