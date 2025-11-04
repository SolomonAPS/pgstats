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
def calculate_v_h(n: int, a1: float, a2: float) -> tuple[float, float]:
    """
    Calculate variance components for Fay and Wu's H test.
    From Walsh and Lynch (2018) Equation 9.27c.
    
    The variance of (θ_π - θ_H) has the form: Var(θ_π - θ_H) = u_H * θ + v_H * θ²
    
    Exact textbook formula:
    [(n-2) / 6(n-1)]θ + [(18n²)(3n+2)b_{n+1} - (88n³ + 9n² - 13n + 6) / 9n(n-1)²]θ²
    
    Args:
        n: Sample size
        a1: a₁ value
        a2: a₂ value
        
    Returns:
        tuple: (u_H, v_H) variance components
    """
    # Calculate b_{n+1} for the formula
    b_n_plus_1 = calculate_b2(n + 1)
    
    # u_H = (n-2) / 6(n-1)
    u_h = (n - 2.0) / (6.0 * (n - 1.0))
    
    # v_H = [(18n²)(3n+2)b_{n+1} - (88n³ + 9n² - 13n + 6)] / 9n(n-1)²
    numerator = (18.0 * n * n) * (3.0 * n + 2.0) * b_n_plus_1 - (88.0 * n * n * n + 9.0 * n * n - 13.0 * n + 6.0)
    denominator = 9.0 * n * (n - 1.0) * (n - 1.0)
    
    v_h = numerator / denominator
    
    return u_h, v_h


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
    
    Textbook formula (Walsh & Lynch 2018, Equation 9.4):
    θπ = (1/C(n,2)) * Σ[i=1 to n-1] i(n-i)ξᵢ
    
    Where:
    - C(n,2) = n(n-1)/2 (number of pairs)
    - ξᵢ = number of sites with i derived alleles
    - n = sample size
    
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
    
    # Get window information (windows are always defined)
    n_windows = len(ds.windows)
    window_starts = ds.window_start_idx.values
    window_stops = ds.window_stop_idx.values
    
    # Calculate θπ for each window
    pi_values = np.zeros(n_windows)
    
    for w_idx in range(n_windows):
        # Extract variants in this window
        window_start = window_starts[w_idx]
        window_stop = window_stops[w_idx]
        window_variant_matrix = variant_matrix[window_start:window_stop, :]
        
        # Calculate θπ for the entire window
        pi_values[w_idx] = calculate_pi(window_variant_matrix)
    
    # Create output dataset
    result = ds.copy()
    result["theta_pi"] = (["windows"], pi_values)
    
    return result


def theta_w(ds: xr.Dataset, call_genotype: str = "call_genotype") -> xr.Dataset:
    """
    Calculate theta_w (θw) - Watterson's estimator based on segregating sites.
    
    References:
    - Walsh and Lynch (2018) Equation 9.5
    - Wakeley (2009) Coalescent Theory, equation 4.40
    
    Textbook formula (Walsh & Lynch 2018, Equation 9.5):
    θw = S/a₁
    
    Where:
    - S = number of segregating sites
    - a₁ = Σ[i=1 to n-1] 1/i (harmonic number)
    - n = sample size
    
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
    
    # Get window information (windows are always defined)
    n_windows = len(ds.windows)
    window_starts = ds.window_start_idx.values
    window_stops = ds.window_stop_idx.values
    
    # Calculate θw for each window
    theta_w_values = np.zeros(n_windows)
    
    for w_idx in range(n_windows):
        # Extract variants in this window
        window_start = window_starts[w_idx]
        window_stop = window_stops[w_idx]
        window_variant_matrix = variant_matrix[window_start:window_stop, :]
        
        # Calculate θw for the entire window
        theta_w_values[w_idx] = calculate_theta_w_per_site(window_variant_matrix)
    
    # Create output dataset
    result = ds.copy()
    result["theta_w"] = (["windows"], theta_w_values)
    
    return result


def theta_h(ds: xr.Dataset, call_genotype: str = "call_genotype") -> xr.Dataset:
    """
    Calculate theta_h (θh) - Fay and Wu's estimator based on high-frequency derived alleles.
    
    References:
    - Walsh and Lynch (2018) Equation 9.6
    - Fay and Wu (2000), equation 2
    
    Textbook formula (Walsh & Lynch 2018, Equation 9.6):
    θh = (1/C(n,2)) * Σ[i=1 to n-1] i²ξᵢ
    
    Where:
    - C(n,2) = n(n-1)/2 (number of pairs)
    - ξᵢ = number of sites with i derived alleles
    - n = sample size
    
    This estimator emphasizes high-frequency derived alleles.
    
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
    
    # Get window information (windows are always defined)
    n_windows = len(ds.windows)
    window_starts = ds.window_start_idx.values
    window_stops = ds.window_stop_idx.values
    
    # Calculate θh for each window
    theta_h_values = np.zeros(n_windows)
    
    for w_idx in range(n_windows):
        # Extract variants in this window
        window_start = window_starts[w_idx]
        window_stop = window_stops[w_idx]
        window_variant_matrix = variant_matrix[window_start:window_stop, :]
        
        # Calculate θh for the entire window
        theta_h_values[w_idx] = calculate_theta_h_per_site(window_variant_matrix)
    
    # Create output dataset
    result = ds.copy()
    result["theta_h"] = (["windows"], theta_h_values)
    
    return result


def theta_l(ds: xr.Dataset, call_genotype: str = "call_genotype") -> xr.Dataset:
    """
    Calculate theta_L (θL) - Zeng et al.'s estimator emphasizing high-frequency sites.
    
    References:
    - Walsh and Lynch (2018) Equation 9.7
    - Zeng et al. (2006), equation 9.28a
    
    Textbook formula (Walsh & Lynch 2018, Equation 9.7):
    θL = (1/(n-1)) * Σ[i=1 to n-1] i·ξᵢ
    
    Where:
    - ξᵢ = number of sites with i derived alleles
    - n = sample size
    
    This estimator places more weight on high-frequency sites compared to θw,
    making it useful for detecting selective sweeps.
    
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
    
    # Get window information (windows are always defined)
    n_windows = len(ds.windows)
    window_starts = ds.window_start_idx.values
    window_stops = ds.window_stop_idx.values
    
    # Calculate θL for each window
    theta_l_values = np.zeros(n_windows)
    
    for w_idx in range(n_windows):
        # Extract variants in this window
        window_start = window_starts[w_idx]
        window_stop = window_stops[w_idx]
        window_variant_matrix = variant_matrix[window_start:window_stop, :]
        
        # Calculate θL for the entire window
        theta_l_values[w_idx] = calculate_theta_l_per_site(window_variant_matrix)
    
    # Create output dataset
    result = ds.copy()
    result["theta_l"] = (["windows"], theta_l_values)
    
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
    
    Textbook formula (Walsh & Lynch 2018, Equation 9.8):
    D = (θπ - θw) / sqrt(Var(θπ - θw))
    
    Where:
    - θπ = nucleotide diversity (pairwise differences)
    - θw = Watterson's theta (segregating sites)
    - Var(θπ - θw) = c₁S + c₂S(S-1)
    - c₁ = b₁ - 1/a₁
    - c₂ = b₂ - (n+2)/(a₁n) + a₂/a₁²
    - b₁ = (n+1)/(3(n-1))
    - b₂ = 2(n²+n+3)/(9n(n-1))
    - a₁ = Σ[i=1 to n-1] 1/i, a₂ = Σ[i=1 to n-1] 1/i²
    - S = number of segregating sites
    
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
    
    # Get window information (windows are always defined)
    n_windows = len(ds.windows)
    window_starts = ds.window_start_idx.values
    window_stops = ds.window_stop_idx.values
    
    # Calculate Tajima's D for each window
    tajima_d_values = np.zeros(n_windows)
    
    for w_idx in range(n_windows):
        # Extract variants in this window
        window_start = window_starts[w_idx]
        window_stop = window_stops[w_idx]
        window_variant_matrix = variant_matrix[window_start:window_stop, :]
        
        # Calculate max sample size for this window
        window_max_n = 0
        for i in range(len(window_variant_matrix)):
            site_n = 0
            for j in range(n_samples):
                if window_variant_matrix[i, j] != -1:
                    site_n += 1
            if site_n > window_max_n:
                window_max_n = site_n
        
        if window_max_n <= 1:
            tajima_d_values[w_idx] = np.nan
            continue
        
        # Calculate harmonic numbers for this window
        a1 = calculate_a1(window_max_n)
        a2 = calculate_a2(window_max_n)
        c1 = calculate_c1(window_max_n, a1)
        c2 = calculate_c2(window_max_n, a1, a2)
        
        # Get SFS for entire window
        sfs, n_max = get_unfolded_sfs(window_variant_matrix)
        S = calculate_S(sfs)
        
        if S == 0:
            tajima_d_values[w_idx] = np.nan
            continue
        
        pi = calculate_pi(window_variant_matrix)
        
        # Calculate variance using window_max_n
        # From Walsh & Lynch 2018, Equation 9.24a (page 305):
        # Var(D) = √(α_D*S + β_D*S²) where α_D = c₁, β_D = c₂
        var = c1 * S + c2 * S * (S - 1)
        
        if var <= 0:
            tajima_d_values[w_idx] = np.nan
        else:
            tajima_d_values[w_idx] = (pi - S / a1) / np.sqrt(var)
    
    # Create output dataset
    result = ds.copy()
    result["tajima_d"] = (["windows"], tajima_d_values)
    
    return result


def fu_li_d(ds: xr.Dataset, call_genotype: str = "call_genotype", folded: bool = True) -> xr.Dataset:
    """
    Calculate Fu and Li's D or D* statistic for each variant.
    
    References:
    - Walsh and Lynch (2018) Equation 9.26b (D*) and 9.26c (D)
    - Fu and Li (1993) Statistical tests of neutrality of mutations
    
    Textbook formulas (Walsh & Lynch 2018):
    
    D* (folded, Equation 9.26b):
    D* = (S/a₁ - ((n-1)/n)η₁) / sqrt(Var(D*))
    
    D (unfolded, Equation 9.26c):
    D = (S/a₁ - ζ₁) / sqrt(Var(D))
    
    Where:
    - S = number of segregating sites
    - a₁ = Σ[i=1 to n-1] 1/i (harmonic number)
    - η₁ = number of sites with minor allele count = 1 (folded SFS)
    - ζ₁ = number of sites with derived allele count = 1 (unfolded SFS)
    - n = sample size
    
    Variance components (Equations 9.26b-c):
    Var(D*) = α*S + β*S(S-1)
    Var(D) = αS + βS(S-1)
    
    Where α*, β*, α, β are functions of harmonic numbers a₁, a₂, b₁, b₂.
    
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
    
    # Get window information (windows are always defined)
    n_windows = len(ds.windows)
    window_starts = ds.window_start_idx.values
    window_stops = ds.window_stop_idx.values
    
    # Calculate Fu and Li's D or D* for each window
    fu_li_d_values = np.zeros(n_windows)
    
    for w_idx in range(n_windows):
        # Extract variants in this window
        window_start = window_starts[w_idx]
        window_stop = window_stops[w_idx]
        window_variant_matrix = variant_matrix[window_start:window_stop, :]
        
        # Calculate max sample size for this window
        window_max_n = 0
        for i in range(len(window_variant_matrix)):
            site_n = 0
            for j in range(n_samples):
                if window_variant_matrix[i, j] != -1:
                    site_n += 1
            if site_n > window_max_n:
                window_max_n = site_n
        
        if window_max_n <= 1:
            fu_li_d_values[w_idx] = np.nan
            continue
        
        # Calculate harmonic numbers for this window
        a1 = calculate_a1(window_max_n)
        a2 = calculate_a2(window_max_n)
        u_d_star, v_d_star = calculate_v_d_star(window_max_n, a1, a2)
        
        # Get SFS for entire window
        if folded:
            # Use folded SFS (D*)
            sfs, n_max = get_folded_sfs(window_variant_matrix)
            singleton_count = sfs[0] if len(sfs) > 0 else 0
        else:
            # Use unfolded SFS (D)
            sfs, n_max = get_unfolded_sfs(window_variant_matrix)
            singleton_count = sfs[0] if len(sfs) > 0 else 0
            
        S = calculate_S(sfs)
        
        if S == 0:
            fu_li_d_values[w_idx] = np.nan
            continue
        
        # Calculate numerator: S/a₁ - singleton estimator
        if folded:
            # D*: S/a₁ - ((n-1)/n)η₁
            numerator = S / a1 - ((window_max_n - 1) / window_max_n) * singleton_count
        else:
            # D: S/a₁ - ζ₁
            numerator = S / a1 - singleton_count
        
        # Calculate variance using window_max_n
        # From Walsh & Lynch 2018, Equation 9.26a-b (page 305):
        # Var(D*) = √(α*S + β*S(S-1))
        var = u_d_star * S + v_d_star * S * (S - 1)
        
        if var <= 0:
            if S == 1:
                fu_li_d_values[w_idx] = numerator / np.sqrt(abs(u_d_star))
            else:
                fu_li_d_values[w_idx] = np.nan
        else:
            fu_li_d_values[w_idx] = numerator / np.sqrt(var)
    
    # Create output dataset
    result = ds.copy()
    stat_name = "fu_li_d" if not folded else "fu_li_d_star"
    result[stat_name] = (["windows"], fu_li_d_values)
    
    return result


def fu_li_f(ds: xr.Dataset, call_genotype: str = "call_genotype", folded: bool = True) -> xr.Dataset:
    """
    Calculate Fu and Li's F or F* statistic for each variant.
    
    References:
    - Walsh and Lynch (2018) Equation 9.26e (F*) and 9.26f (F)
    - Fu and Li (1993) Statistical tests of neutrality of mutations
    
    Textbook formulas (Walsh & Lynch 2018):
    
    F* (folded, Equation 9.26e):
    F* = (θπ - ((n-1)/n)η₁) / sqrt(Var(F*))
    
    F (unfolded, Equation 9.26f):
    F = (θπ - ζ₁) / sqrt(Var(F))
    
    Where:
    - θπ = nucleotide diversity (pairwise differences)
    - η₁ = number of sites with minor allele count = 1 (folded SFS)
    - ζ₁ = number of sites with derived allele count = 1 (unfolded SFS)
    - n = sample size
    
    Variance components (Equations 9.26e-f):
    Var(F*) = α_F*S + β_F*S(S-1)
    Var(F) = α_FS + β_FS(S-1)
    
    Where α_F*, β_F*, α_F, β_F are functions of harmonic numbers a₁, a₂, b₁, b₂.
    
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
    
    # Get window information (windows are always defined)
    n_windows = len(ds.windows)
    window_starts = ds.window_start_idx.values
    window_stops = ds.window_stop_idx.values
    
    # Calculate Fu and Li's F or F* for each window
    fu_li_f_values = np.zeros(n_windows)
    
    for w_idx in range(n_windows):
        # Extract variants in this window
        window_start = window_starts[w_idx]
        window_stop = window_stops[w_idx]
        window_variant_matrix = variant_matrix[window_start:window_stop, :]
        
        # Calculate max sample size for this window
        window_max_n = 0
        for i in range(len(window_variant_matrix)):
            site_n = 0
            for j in range(n_samples):
                if window_variant_matrix[i, j] != -1:
                    site_n += 1
            if site_n > window_max_n:
                window_max_n = site_n
        
        if window_max_n <= 1:
            fu_li_f_values[w_idx] = np.nan
            continue
        
        # Calculate harmonic numbers for this window
        a1 = calculate_a1(window_max_n)
        a2 = calculate_a2(window_max_n)
        u_f_star, v_f_star = calculate_v_f_star(window_max_n, a1, a2)
        
        # Get SFS for entire window
        if folded:
            # Use folded SFS (F*)
            sfs, n_max = get_folded_sfs(window_variant_matrix)
            singleton_count = sfs[0] if len(sfs) > 0 else 0
        else:
            # Use unfolded SFS (F)
            sfs, n_max = get_unfolded_sfs(window_variant_matrix)
            singleton_count = sfs[0] if len(sfs) > 0 else 0
            
        S = calculate_S(sfs)
        
        if S == 0:
            fu_li_f_values[w_idx] = np.nan
            continue
            
        # Calculate θπ for entire window
        pi = calculate_pi(window_variant_matrix)
        
        # Calculate numerator: π - singleton estimator
        if folded:
            # F*: π - ((n-1)/n)η₁
            numerator = pi - ((window_max_n - 1) / window_max_n) * singleton_count
        else:
            # F: π - ζ₁
            numerator = pi - singleton_count
        
        # Calculate variance using window harmonic numbers
        # From Walsh & Lynch 2018, Equation 9.26d-e (page 305):
        # Var(F*) = √(α_F*S + β_F*S(S-1))
        var = u_f_star * S + v_f_star * S * (S - 1)
        
        if var <= 0:
            if S == 1:
                fu_li_f_values[w_idx] = numerator / np.sqrt(abs(u_f_star))
            else:
                fu_li_f_values[w_idx] = np.nan
        else:
            fu_li_f_values[w_idx] = numerator / np.sqrt(var)
    
    # Create output dataset
    result = ds.copy()
    stat_name = "fu_li_f" if not folded else "fu_li_f_star"
    result[stat_name] = (["windows"], fu_li_f_values)
    
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
    
    Textbook formula (Walsh & Lynch 2018, Equation 9.28c):
    E = (θL - θw) / sqrt(Var(θL - θw))
    
    Where:
    - θL = Zeng et al.'s theta emphasizing high-frequency sites
    - θw = Watterson's theta (segregating sites)
    - Var(θL - θw) = [n/(2(n-1)) - 1/a₁]θ + [b₂ + 2(n/(n-1))²b₂ - 2(nb₂-n+1)/((n-1)a₁) - (3n+1)/(n-1)]θ²
    - a₁ = Σ[i=1 to n-1] 1/i, b₂ = 2(n²+n+3)/(9n(n-1))
    - n = sample size
    
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
    
    # Get window information (windows are always defined)
    n_windows = len(ds.windows)
    window_starts = ds.window_start_idx.values
    window_stops = ds.window_stop_idx.values
    
    # Calculate Zeng's E for each window
    zeng_e_values = np.full(n_windows, np.nan)
    
    if n_windows > 0:
        print(f"DEBUG zeng_e: n_windows={n_windows}, first window variants: {window_stops[0] - window_starts[0]}", flush=True)
    
    for w_idx in range(n_windows):
        # Extract variants in this window
        window_start = window_starts[w_idx]
        window_stop = window_stops[w_idx]
        window_variant_matrix = variant_matrix[window_start:window_stop, :]
        
        # Calculate max sample size for this window
        window_max_n = 0
        for i in range(len(window_variant_matrix)):
            site_n = 0
            for j in range(n_samples):
                if window_variant_matrix[i, j] != -1:
                    site_n += 1
            if site_n > window_max_n:
                window_max_n = site_n
        
        if window_max_n <= 1:
            zeng_e_values[w_idx] = np.nan
            continue
        
        # Calculate harmonic numbers for this window
        a1 = calculate_a1(window_max_n)
        b2 = calculate_b2(window_max_n)
        
        # Pre-calculate variance components using window_max_n
        # From Walsh & Lynch 2018, Equation 9.28c (page 307)
        term1 = (window_max_n / (2.0 * (window_max_n - 1.0))) - (1.0 / a1)
        term2 = b2 + 2.0 * (window_max_n / (window_max_n - 1.0)) ** 2 * b2
        term2 -= 2.0 * (window_max_n * b2 - window_max_n + 1.0) / ((window_max_n - 1.0) * a1)
        term2 -= (3.0 * window_max_n + 1.0) / (window_max_n - 1.0)
        
        # DEBUG: Print variance components for first window
        if w_idx == 0:
            print(f"DEBUG zeng_e variance components: n={window_max_n}, a1={a1:.4f}, b2={b2:.4f}, term1={term1:.4f}, term2={term2:.4f}", flush=True)
        
        # Get unfolded SFS for entire window
        sfs, n_max = get_unfolded_sfs(window_variant_matrix)
        S = calculate_S(sfs)
        
        if S == 0:
            zeng_e_values[w_idx] = np.nan
            continue
        
        # Calculate theta_L and theta_w from SFS (not per-site)
        theta_l_val = calculate_theta_l(sfs, window_max_n)
        theta_w_val = calculate_theta_w(sfs, window_max_n)
        
        # DEBUG: Check if theta values are zero
        if w_idx < 3:  # Only print for first few windows
            b1_temp = calculate_b1(window_max_n)
            theta_temp = S / a1
            theta_sq_temp = S * (S - 1) / (a1 * a1 + b1_temp)
            var_e_temp = term1 * theta_temp + term2 * theta_sq_temp
            print(f"DEBUG zeng_e window {w_idx}: S={S}, theta_l={theta_l_val}, theta_w={theta_w_val}, var_e={var_e_temp}", flush=True)
        
        # Calculate variance of E
        # From Walsh & Lynch 2018, Equation 9.28c (page 307):
        # Var(θ_L - θ_W) uses scaled theta from Equation 9.21b (page 301):
        # θ → S/a₁ and θ² → S(S-1)/(a₁² + b₁)
        # This is DIFFERENT from Tajima's D which uses S directly
        b1 = calculate_b1(window_max_n)
        theta_for_var = S / a1
        theta_sq_for_var = S * (S - 1) / (a1 * a1 + b1)
        var_e = term1 * theta_for_var + term2 * theta_sq_for_var
        
        if var_e <= 0:
            if S == 1:
                zeng_e_values[w_idx] = (theta_l_val - theta_w_val) / np.sqrt(abs(term1))
            else:
                zeng_e_values[w_idx] = np.nan
        else:
            zeng_e_values[w_idx] = (theta_l_val - theta_w_val) / np.sqrt(var_e)
    
    # Create output dataset
    result = ds.copy()
    result["zeng_e"] = (["windows"], zeng_e_values)
    
    return result


def singletons(ds: xr.Dataset, call_genotype: str = "call_genotype", folded: bool = True) -> xr.Dataset:
    """
    Calculate singleton count for each variant.
    
    References:
    - Walsh and Lynch (2018) Equation 9.26b (folded) and 9.26c (unfolded)
    - Fu and Li (1993) Statistical tests of neutrality of mutations
    
    Textbook formulas (Walsh & Lynch 2018):
    
    Folded singletons (η₁, Equation 9.26b):
    η₁ = number of sites with minor allele count = 1
    
    Unfolded singletons (ζ₁, Equation 9.26c):
    ζ₁ = number of sites with derived allele count = 1
    
    Where:
    - Minor allele count = min(derived_count, n - derived_count)
    - Derived allele count = number of derived alleles in sample
    - n = sample size
    
    Singletons are sites where an allele appears only once in the sample.
    They are used in Fu and Li's tests and are sensitive to population size changes.
    
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

    # Get window information (windows are always defined)
    n_windows = len(ds.windows)
    window_starts = ds.window_start_idx.values
    window_stops = ds.window_stop_idx.values
    
    # Calculate singletons for each window
    singleton_values = np.zeros(n_windows)
    
    for w_idx in range(n_windows):
        # Extract variants in this window
        window_start = window_starts[w_idx]
        window_stop = window_stops[w_idx]
        window_variant_matrix = variant_matrix[window_start:window_stop, :]
        
        if folded:
            # Use folded SFS (η₁) for entire window
            sfs, n = get_folded_sfs(window_variant_matrix)
            singleton_count = sfs[0] if len(sfs) > 0 else 0
        else:
            # Use unfolded SFS (ζ₁) for entire window
            sfs, n = get_unfolded_sfs(window_variant_matrix)
            singleton_count = sfs[0] if len(sfs) > 0 else 0
            
        singleton_values[w_idx] = singleton_count
    
    # Create output dataset
    result = ds.copy()
    if folded:
        result["singletons_folded"] = (["windows"], singleton_values)
    else:
        result["singletons_unfolded"] = (["windows"], singleton_values)
    
    return result


def fay_wu_h(ds: xr.Dataset, call_genotype: str = "call_genotype") -> xr.Dataset:
    """
    Calculate Fay and Wu's H statistic for each variant.
    
    References:
    - Walsh and Lynch (2018) Equation 9.27b
    - Fay and Wu (2000) Statistical tests of neutrality of mutations
    - Zeng et al. (2006) Statistical Tests for Detecting Positive Selection
    
    Textbook formula (Walsh & Lynch 2018, Equation 9.27b):
    H = (θπ - θh) / sqrt(Var(θπ - θh))
    
    Where:
    - θπ = nucleotide diversity (pairwise differences)
    - θh = Fay and Wu's theta emphasizing high-frequency derived alleles
    - Var(θπ - θh) = u_Hθ + v_Hθ²
    
    Variance formula (Equation 9.27c):
    u_H = (n-2)/(6(n-1))
    v_H = [(18n²)(3n+2)b_{n+1} - (88n³ + 9n² - 13n + 6)] / 9n(n-1)²
    
    Where:
    - b_{n+1} = 2((n+1)²+(n+1)+3)/(9(n+1)(n+1-1)) = 2(n²+3n+5)/(9n(n+1))
    - n = sample size
    
    Fay and Wu's H tests for neutrality by comparing theta_pi and theta_h.
    
    Note: Following scikit-allel and sgkit conventions, harmonic numbers (a1, a2, etc.)
    are calculated using the maximum observed sample size across all variants.
    This maintains theoretical consistency with the assumption of constant n.
    
    Args:
        ds: sgkit Dataset containing genotype calls
        call_genotype: Name of the genotype variable
        
    Returns:
        Dataset with Fay and Wu's H values
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
    
    # Get window information (windows are always defined)
    n_windows = len(ds.windows)
    window_starts = ds.window_start_idx.values
    window_stops = ds.window_stop_idx.values
    
    # Calculate Fay and Wu's H for each window
    fay_wu_h_values = np.full(n_windows, np.nan)
    
    if n_windows > 0:
        print(f"DEBUG fay_wu_h: n_windows={n_windows}, first window variants: {window_stops[0] - window_starts[0]}", flush=True)
    
    for w_idx in range(n_windows):
        # Extract variants in this window
        window_start = window_starts[w_idx]
        window_stop = window_stops[w_idx]
        window_variant_matrix = variant_matrix[window_start:window_stop, :]
        
        # Calculate max sample size for this window
        window_max_n = 0
        for i in range(len(window_variant_matrix)):
            site_n = 0
            for j in range(n_samples):
                if window_variant_matrix[i, j] != -1:
                    site_n += 1
            if site_n > window_max_n:
                window_max_n = site_n
        
        if window_max_n <= 1:
            fay_wu_h_values[w_idx] = np.nan
            continue
        
        # Calculate harmonic numbers for this window
        a1 = calculate_a1(window_max_n)
        a2 = calculate_a2(window_max_n)
        u_h, v_h = calculate_v_h(window_max_n, a1, a2)
        
        # DEBUG: Print variance components for first window
        if w_idx == 0:
            print(f"DEBUG fay_wu_h variance components: n={window_max_n}, a1={a1:.4f}, a2={a2:.4f}, u_h={u_h:.4f}, v_h={v_h:.4f}", flush=True)
        
        # Get SFS for entire window
        sfs, n_max = get_unfolded_sfs(window_variant_matrix)
        S = calculate_S(sfs)
        
        if S == 0:
            fay_wu_h_values[w_idx] = np.nan
            continue
            
        pi = calculate_pi(window_variant_matrix)
        theta_h_val = calculate_theta_h(sfs, window_max_n)
        
        # Calculate variance using window harmonic numbers
        # From Walsh & Lynch 2018, Equation 9.27c (page 306):
        # Var(θ_π - θ_H) uses scaled theta from Equation 9.21b (page 301):
        # θ → S/a₁ and θ² → S(S-1)/(a₁² + b₁)
        # This is DIFFERENT from Tajima's D which uses S directly
        b1 = calculate_b1(window_max_n)
        theta_for_var = S / a1
        theta_sq_for_var = S * (S - 1) / (a1 * a1 + b1)
        var_h = u_h * theta_for_var + v_h * theta_sq_for_var
        
        # DEBUG: Check if theta_h is zero
        if w_idx < 3:  # Only print for first few windows
            print(f"DEBUG fay_wu_h window {w_idx}: S={S}, pi={pi}, theta_h={theta_h_val}, var_h={var_h}", flush=True)
        
        if var_h <= 0:
            if S == 1:
                fay_wu_h_values[w_idx] = (pi - theta_h_val) / np.sqrt(abs(u_h))
            else:
                fay_wu_h_values[w_idx] = np.nan
        else:
            fay_wu_h_values[w_idx] = (pi - theta_h_val) / np.sqrt(var_h)
    
    # Create output dataset
    result = ds.copy()
    result["fay_wu_h"] = (["windows"], fay_wu_h_values)
    
    return result

