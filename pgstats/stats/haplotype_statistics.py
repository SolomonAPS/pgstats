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
import warnings
from typing import Tuple, List


# =============================================================================
# HAPLOTYPE ESTIMATION AND HASHING
# =============================================================================

@numba.njit(nogil=True, fastmath=False)
def hash_haplotype(haplotype: np.ndarray, ignore_missing: bool = False) -> int:
    """
    Create a hash for a haplotype sequence using sgkit's DJBX33A hash function.
    
    Args:
        haplotype: Array of shape (n_variants,) with haplotype data
                  0 = reference allele, 1 = alternate allele, -1 = missing
        ignore_missing: If True, ignore missing data (-1) in hashing
                       If False, include missing data in hash (default sgkit behavior)
    
    Returns:
        Hash value for the haplotype
    """
    # DJBX33A hash function (matches sgkit exactly)
    hash_value = 5381
    for i in range(haplotype.shape[0]):
        if ignore_missing and haplotype[i] == -1:
            # Skip missing data when ignore_missing=True
            continue
        hash_value = hash_value * 33 + haplotype[i]
    
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
def hash_haplotypes(haplotypes: np.ndarray, ignore_missing: bool = False) -> np.ndarray:
    """
    Hash all haplotypes to create unique identifiers.
    
    This is the core of sgkit's approach - it creates unique identifiers
    for each haplotype by treating the sequence as a base-3 number.
    
    Args:
        haplotypes: Array of shape (n_variants, n_haplotypes) with haplotype data
        ignore_missing: If True, ignore missing data (-1) in hashing
                       If False, include missing data in hash (default sgkit behavior)
    
    Returns:
        Array of shape (n_haplotypes,) with hash values
    """
    n_variants, n_haplotypes = haplotypes.shape
    hash_values = np.zeros(n_haplotypes, dtype=np.int64)
    
    for haplotype_idx in range(n_haplotypes):
        hash_values[haplotype_idx] = hash_haplotype(haplotypes[:, haplotype_idx], ignore_missing)
    
    return hash_values


@numba.njit(nogil=True, fastmath=False)
def filter_haplotypes_by_missingness(haplotypes: np.ndarray, 
                                   max_missing: float,
                                   missing_is_percentage: bool = True) -> np.ndarray:
    """
    Filter haplotypes based on missing data threshold.
    
    This function removes entire haplotypes that exceed the missingness threshold,
    then statistics are calculated only on the remaining haplotypes.
    
    Args:
        haplotypes: Array of shape (n_variants, n_haplotypes) with haplotype data
        max_missing: Maximum allowed missing data (percentage 0-1 or absolute count)
        missing_is_percentage: If True, max_missing is percentage (0-1), 
                              if False, max_missing is absolute count
    
    Returns:
        Boolean array of shape (n_haplotypes,) indicating which haplotypes pass the filter
    """
    n_variants, n_haplotypes = haplotypes.shape
    
    # If no variants, return all False (no valid haplotypes)
    if n_variants == 0:
        return np.zeros(n_haplotypes, dtype=numba.boolean)
    
    valid_haplotypes = np.ones(n_haplotypes, dtype=numba.boolean)
    
    for i in range(n_haplotypes):
        missing_count = 0
        for j in range(n_variants):
            if haplotypes[j, i] == -1:
                missing_count += 1
        
        if missing_is_percentage:
            # max_missing is a percentage (0-1)
            missing_fraction = missing_count / n_variants
            if missing_fraction > max_missing:
                valid_haplotypes[i] = False
        else:
            # max_missing is an absolute count
            if missing_count > max_missing:
                valid_haplotypes[i] = False
    
    return valid_haplotypes


# =============================================================================
# PAIRWISE COMPARISON FOR MISSING DATA HANDLING
# =============================================================================

@numba.njit(nogil=True, fastmath=False)
def haplotypes_match_pairwise(haplotype1: np.ndarray, haplotype2: np.ndarray, 
                               min_sites_compared: int = 1) -> bool:
    """
    Compare two haplotypes, returning True only if they match at all positions
    where BOTH have non-missing data.
    
    This function implements pairwise comparison ignoring missing data:
    - Positions where either haplotype has missing data (-1) are excluded
    - Positions where both have non-missing data must match exactly
    - Two haplotypes are considered "the same" if they match at all compared positions
    - If there are no positions where both have non-missing data, returns False
      (cannot determine if they match, so treat as different)
    - If the number of compared sites is less than min_sites_compared, returns False
      (insufficient overlap to confidently call them the same)
    
    Args:
        haplotype1: Array of shape (n_variants,) with haplotype data
                   0 = reference allele, 1 = alternate allele, -1 = missing
        haplotype2: Array of shape (n_variants,) with haplotype data
                   0 = reference allele, 1 = alternate allele, -1 = missing
        min_sites_compared: Minimum number of sites that must be compared
                          (both non-missing) for a match to be considered valid.
                          Default is 1. Higher values prevent spurious matches
                          between haplotypes with little overlap.
    
    Returns:
        True if haplotypes match at all positions where both have non-missing data
        AND at least min_sites_compared positions were compared,
        False otherwise (including case where insufficient positions can be compared)
    """
    n_variants = haplotype1.shape[0]
    n_compared = 0
    
    for pos in range(n_variants):
        val1 = haplotype1[pos]
        val2 = haplotype2[pos]
        
        # Skip if either is missing
        if val1 == -1 or val2 == -1:
            continue
        
        # If both are non-missing, they must match
        n_compared += 1
        if val1 != val2:
            return False
    
    # Check if we compared enough sites
    if n_compared < min_sites_compared:
        return False
    
    # All compared positions matched and we had sufficient overlap
    return True


@numba.njit(nogil=True, fastmath=False)
def union_find_find(parent: np.ndarray, x: int) -> int:
    """
    Find the root of x (iterative version for Numba compatibility).
    
    Args:
        parent: Parent array for union-find structure
        x: Element to find root for
    
    Returns:
        Root of x
    """
    # Find root iteratively (simple version without path compression to avoid issues)
    root = x
    while parent[root] != root:
        root = parent[root]
    return root


@numba.njit(nogil=True, fastmath=False)
def union_find_union(parent: np.ndarray, size: np.ndarray, x: int, y: int):
    """
    Union two sets using union by size.
    
    Args:
        parent: Parent array for union-find structure
        size: Size array for union by size optimization
        x: First element
        y: Second element
    """
    root_x = union_find_find(parent, x)
    root_y = union_find_find(parent, y)
    
    if root_x == root_y:
        return  # Already in same set
    
    # Union by size: attach smaller tree to larger tree
    if size[root_x] < size[root_y]:
        parent[root_x] = root_y
        size[root_y] += size[root_x]
    else:
        parent[root_y] = root_x
        size[root_x] += size[root_y]


@numba.njit(nogil=True, fastmath=False)
def union_find_group_haplotypes(haplotypes: np.ndarray, 
                                min_sites_compared: int = 1) -> Tuple[np.ndarray, np.ndarray]:
    """
    Group haplotypes using union-find based on pairwise comparison.
    
    Two haplotypes are grouped together if they match at all positions where
    both have non-missing data. This handles transitive relationships correctly:
    if A matches B and A matches C, then A, B, and C are all in the same group,
    even if B and C don't directly match.
    
    This function is used when ignore_missing=True to correctly handle missing
    data without imputation. It is slower than hashing (O(n²) comparisons) but
    necessary for correctness.
    
    Args:
        haplotypes: Array of shape (n_variants, n_haplotypes) with haplotype data
                   0 = reference allele, 1 = alternate allele, -1 = missing
        min_sites_compared: Minimum number of sites that must be compared
                          (both non-missing) for two haplotypes to be grouped together.
                          Default is 1. Higher values prevent spurious grouping of
                          haplotypes with little overlap.
    
    Returns:
        Tuple of (group_ids, counts) where:
        - group_ids: Array of shape (n_haplotypes,) with group ID for each haplotype
        - counts: Array of shape (n_groups,) with count of haplotypes in each group
    """
    n_variants, n_haplotypes = haplotypes.shape
    
    if n_haplotypes == 0:
        return np.empty(0, dtype=np.int64), np.empty(0, dtype=np.int64)
    
    # Initialize union-find structure
    parent = np.arange(n_haplotypes, dtype=np.int64)
    size = np.ones(n_haplotypes, dtype=np.int64)
    
    # Compare all pairs of haplotypes
    for i in range(n_haplotypes):
        for j in range(i + 1, n_haplotypes):
            if haplotypes_match_pairwise(haplotypes[:, i], haplotypes[:, j], 
                                        min_sites_compared):
                union_find_union(parent, size, i, j)
    
    # Assign group IDs to each haplotype (compress all paths first)
    group_ids = np.zeros(n_haplotypes, dtype=np.int64)
    for i in range(n_haplotypes):
        group_ids[i] = union_find_find(parent, i)
    
    # Map group IDs to consecutive integers starting from 0
    unique_groups = np.unique(group_ids)
    n_groups = len(unique_groups)
    
    # Create mapping array (max group_id will be < n_haplotypes)
    max_group_id = np.max(group_ids) if n_haplotypes > 0 else 0
    group_map = np.zeros(max_group_id + 1, dtype=np.int64)
    for idx, group_id in enumerate(unique_groups):
        group_map[group_id] = idx
    
    # Remap group IDs
    for i in range(n_haplotypes):
        group_ids[i] = group_map[group_ids[i]]
    
    # Count frequencies for each group
    counts = np.zeros(n_groups, dtype=np.int64)
    for i in range(n_haplotypes):
        counts[group_ids[i]] += 1
    
    return group_ids, counts


def detect_bridge_events(haplotypes: np.ndarray, group_ids: np.ndarray, 
                         min_sites_compared: int = 1) -> Tuple[int, List[Tuple[int, int, int]]]:
    """
    Detect bridge events: cases where two haplotypes are in the same group but 
    don't directly match each other (they were bridged by a third haplotype).
    
    This function checks each group with more than one haplotype to see if all
    members directly match each other. If two members A and B are in the same
    group but don't directly match, it indicates they were bridged by some other
    haplotype C where A matches C and C matches B.
    
    Args:
        haplotypes: Array of shape (n_variants, n_haplotypes) with haplotype data
        group_ids: Array of shape (n_haplotypes,) with group ID for each haplotype
        min_sites_compared: Minimum number of sites for a valid match (same as used
                          in grouping)
    
    Returns:
        Tuple of (n_events, event_list) where:
        - n_events: Total number of bridge events detected (pairs that don't match)
        - event_list: List of tuples (group_id, hap_i, hap_j) for each non-matching pair
    """
    n_variants, n_haplotypes = haplotypes.shape
    
    if n_haplotypes == 0:
        return 0, []
    
    # Get unique groups and their members
    unique_groups = np.unique(group_ids)
    events = []
    
    for group_id in unique_groups:
        # Find all haplotypes in this group
        members = np.where(group_ids == group_id)[0]
        
        # Only check groups with more than one member
        if len(members) <= 1:
            continue
        
        # Check all pairs within this group
        for i in range(len(members)):
            for j in range(i + 1, len(members)):
                hap_i = members[i]
                hap_j = members[j]
                
                # Check if these two haplotypes directly match
                # We need to convert the Numba function call to work in regular Python
                n_compared = 0
                mismatch = False
                
                for pos in range(n_variants):
                    val1 = haplotypes[pos, hap_i]
                    val2 = haplotypes[pos, hap_j]
                    
                    if val1 == -1 or val2 == -1:
                        continue
                    
                    n_compared += 1
                    if val1 != val2:
                        mismatch = True
                        break
                
                # If they don't directly match (either mismatch or insufficient overlap)
                # but are in the same group, this is a bridge event
                if mismatch or n_compared < min_sites_compared:
                    events.append((int(group_id), int(hap_i), int(hap_j)))
    
    return len(events), events


# =============================================================================
# HAPLOTYPE STATISTICS CALCULATIONS
# =============================================================================

@numba.njit(nogil=True, fastmath=False)
def calculate_haplotype_diversity(haplotypes: np.ndarray, ignore_missing: bool = False, 
                                  min_sites_compared: int = 1) -> float:
    """
    Calculate haplotype diversity (H) - the probability that two randomly
    chosen haplotypes are different.
    
    H = (1 - sum(f_i^2)) * n / (n - 1)
    where f_i is the frequency of haplotype i and n is the number of haplotypes.
    
    When ignore_missing=True, uses pairwise comparison with union-find to correctly
    group haplotypes that match at all positions where both have non-missing data.
    This handles transitive relationships (A matches B and A matches C, but B and C
    don't directly match). When ignore_missing=False, uses faster hashing approach.
    
    Args:
        haplotypes: Array of shape (n_variants, n_haplotypes) with haplotype data
                  0 = reference allele, 1 = alternate allele, -1 = missing
        ignore_missing: If True, use pairwise comparison ignoring missing data
                        If False, include missing data in hash (default sgkit behavior)
        min_sites_compared: Minimum number of non-missing overlapping sites required
                          for two haplotypes to be considered matching (only used when
                          ignore_missing=True). Default is 1.
    
    Returns:
        Haplotype diversity value
    """
    n_variants, n_haplotypes = haplotypes.shape
    
    # Check for empty window or insufficient haplotypes
    if n_variants == 0 or n_haplotypes < 2:
        return np.nan
    
    # Use pairwise comparison when ignore_missing=True, hashing otherwise
    if ignore_missing:
        # Pairwise comparison with union-find (slower but correct for missing data)
        group_ids, counts = union_find_group_haplotypes(haplotypes, min_sites_compared)
    else:
        # Hash all haplotypes to identify unique haplotypes (faster)
        hash_values = hash_haplotypes(haplotypes, ignore_missing)
        # Count haplotype frequencies
        unique_hashes, counts = count_unique_values(hash_values)
    
    # Protect against division by zero
    if n_haplotypes == 0 or len(counts) == 0:
        return np.nan
    
    # Calculate frequencies
    frequencies = counts.astype(np.float64) / n_haplotypes
    
    # Calculate haplotype diversity using Nei's formula
    # H = (1 - sum(f_i^2)) * n / (n - 1)
    sum_squared_freq = 0.0
    for i in range(len(frequencies)):
        sum_squared_freq += frequencies[i] * frequencies[i]
    
    # Apply finite population correction: n / (n - 1)
    if n_haplotypes > 1:
        diversity = (1.0 - sum_squared_freq) * n_haplotypes / (n_haplotypes - 1.0)
    else:
        diversity = 0.0
    
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
def calculate_garud_h_statistics(haplotypes: np.ndarray,
                               ignore_missing: bool = False,
                               max_missing: float = 1.0,
                               missing_is_percentage: bool = True,
                               min_sites_compared: int = 1) -> np.ndarray:
    """
    Calculate Garud H1, H12, H123, and H2/H1 statistics.
    
    When ignore_missing=True, uses pairwise comparison with union-find to correctly
    group haplotypes that match at all positions where both have non-missing data.
    This handles transitive relationships (A matches B and A matches C, but B and C
    don't directly match). When ignore_missing=False, uses faster hashing approach.
    
    References:
    - Walsh and Lynch (2018) Equation 9.10-9.13
    - Garud et al. (2015) A selective sweep on the Drosophila X chromosome
    
    Args:
        haplotypes: Array of shape (n_variants, n_haplotypes) with haplotype data
        ignore_missing: If True, use pairwise comparison ignoring missing data
                       If False, include missing data in hash (default sgkit behavior)
        max_missing: Maximum allowed missing data (percentage 0-1 or absolute count)
                    Default: 1.0 (100% missing data allowed - no filtering)
        missing_is_percentage: If True, max_missing is percentage (0-1), 
                              if False, max_missing is absolute count
                              Default: True
        min_sites_compared: Minimum number of non-missing overlapping sites required
                          for two haplotypes to be considered matching (only used when
                          ignore_missing=True). Default is 1.
    
    Returns:
        Array of [H1, H12, H123, H2/H1] values
    """
    n_variants, n_haplotypes = haplotypes.shape
    
    # Check for empty window or insufficient haplotypes
    if n_variants == 0 or n_haplotypes < 2:
        return np.array([np.nan, np.nan, np.nan, np.nan])
    
    # Filter haplotypes by missingness threshold
    valid_mask = filter_haplotypes_by_missingness(haplotypes, max_missing, missing_is_percentage)
    valid_haplotypes = haplotypes[:, valid_mask]
    n_valid_haplotypes = np.sum(valid_mask)
    
    if n_valid_haplotypes < 2:
        return np.array([np.nan, np.nan, np.nan, np.nan])
    
    # Use pairwise comparison when ignore_missing=True, hashing otherwise
    if ignore_missing:
        # Pairwise comparison with union-find (slower but correct for missing data)
        group_ids, counts = union_find_group_haplotypes(valid_haplotypes, min_sites_compared)
    else:
        # Hash all haplotypes to create unique identifiers (faster)
        hash_values = hash_haplotypes(valid_haplotypes, ignore_missing)
        # Count haplotype frequencies
        unique_hashes, counts = count_unique_values(hash_values)
    
    # Protect against division by zero
    if n_valid_haplotypes == 0 or len(counts) == 0:
        return np.array([np.nan, np.nan, np.nan, np.nan])
    
    frequencies = counts / n_valid_haplotypes
    
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
                       call_genotype: str = "call_genotype",
                       ignore_missing: bool = False,
                       min_sites_compared: int = 1) -> xr.Dataset:
    """
    Calculate haplotype diversity for each window using Nei's gene diversity formula.
    
    Uses the formula: H = (1 - sum(f_i^2)) * n / (n - 1)
    where f_i is the frequency of haplotype i and n is the number of haplotypes.
    
    References:
    - Walsh and Lynch (2018) Equation 9.9
    - Nei (1987) Molecular Evolutionary Genetics
    
    Args:
        ds: sgkit Dataset containing genotype calls
        call_genotype: Name of the genotype variable
        ignore_missing: If True, ignore missing data (-1) in hashing
                        If False, include missing data in hash (default sgkit behavior)
        min_sites_compared: Minimum number of non-missing overlapping sites required
                          for two haplotypes to be considered matching (only used when
                          ignore_missing=True). Default is 1.
        
    Returns:
        Dataset with haplotype diversity values
    """
    genotypes = ds[call_genotype].values
    n_variants, n_samples, ploidy = genotypes.shape
    
    # Estimate haplotypes assuming phased data (sgkit approach)
    haplotypes = estimate_haplotypes_from_phased_data(genotypes)
    
    # Get window information (windows are always defined)
    n_windows = len(ds.windows)
    window_starts = ds.window_start_idx.values
    window_stops = ds.window_stop_idx.values
    
    # Get window positions for debug output
    positions = ds.variant_position.values
    window_start_positions = []
    window_stop_positions = []
    if 'window_start' in ds.coords:
        window_start_positions = ds.window_start.values
    if 'window_stop' in ds.coords:
        window_stop_positions = ds.window_stop.values
    
    # Calculate haplotype diversity for each window
    diversity_values = np.zeros(n_windows)
    
    # Write debug output to file (append mode for multiple runs, TSV format)
    import os
    from datetime import datetime
    debug_file = os.environ.get('PGSTATS_DEBUG_FILE', 'haplotype_diversity_debug.txt')
    run_id = os.environ.get('PGSTATS_RUN_ID', datetime.now().strftime('%Y%m%d_%H%M%S'))
    
    # Write header only if file doesn't exist (TSV format for easy pandas import)
    file_exists = os.path.exists(debug_file)
    with open(debug_file, 'a') as f:
        if not file_exists:
            f.write("Run_ID\tWindow\tWin_Size\tN_Var\tN_Hap\tMissing_Rate\tN_Var_Miss\tN_Hap_Miss\t"
                    "N_Unique\tN_Unique_Ig\tN_Sing\tN_Sing_Ig\tDiversity\tDiversity_Ig\t"
                    "Sum_Freq2\tSum_Freq2_Ig\tTop1_Freq\tTop2_Freq\tTop3_Freq\n")
    
    for w_idx in range(n_windows):
        window_start = window_starts[w_idx]
        window_stop = window_stops[w_idx]
        window_haplotypes = haplotypes[window_start:window_stop, :]
        diversity_values[w_idx] = calculate_haplotype_diversity(window_haplotypes, ignore_missing, 
                                                                 min_sites_compared)
        
        # Calculate debug info
        n_variants_win = window_haplotypes.shape[0]
        n_haplotypes_win = window_haplotypes.shape[1]
        
        # Calculate missing data metrics
        total_calls = n_variants_win * n_haplotypes_win
        missing_calls = np.sum(window_haplotypes == -1)
        missing_rate = missing_calls / total_calls if total_calls > 0 else 0.0
        
        # Count variants with any missing data
        variants_with_missing = np.sum(np.any(window_haplotypes == -1, axis=1))
        # Count haplotypes with any missing data
        haplotypes_with_missing = np.sum(np.any(window_haplotypes == -1, axis=0))
        
        # With ignore_missing=False
        hash_values_false = hash_haplotypes(window_haplotypes, False)
        unique_hashes_false, counts_false = count_unique_values(hash_values_false)
        n_unique_false = len(unique_hashes_false)
        n_singletons_false = np.sum(counts_false == 1)
        frequencies_false = counts_false.astype(np.float64) / n_haplotypes_win
        sum_freq_sq_false = np.sum(frequencies_false ** 2)
        diversity_false = calculate_haplotype_diversity(window_haplotypes, False)
        frequencies_false_sorted = np.sort(frequencies_false)[::-1]
        top3_freq_false = frequencies_false_sorted[:3] if len(frequencies_false_sorted) >= 3 else list(frequencies_false_sorted) + [0.0] * (3 - len(frequencies_false_sorted))
        
        # With ignore_missing=True
        # Use union-find for grouping (same as in calculate_haplotype_diversity)
        group_ids_true, counts_true = union_find_group_haplotypes(window_haplotypes, min_sites_compared)
        n_unique_true = len(counts_true)
        n_singletons_true = np.sum(counts_true == 1)
        
        # Detect bridge events (only when using ignore_missing)
        if ignore_missing and n_haplotypes_win > 1:
            n_bridge_events, bridge_event_list = detect_bridge_events(window_haplotypes, group_ids_true, 
                                                                       min_sites_compared)
            if n_bridge_events > 0:
                warnings.warn(
                    f"Window {w_idx}: Detected {n_bridge_events} bridge event(s) where haplotypes "
                    f"are grouped together but don't directly match. This indicates transitive grouping "
                    f"via intermediate 'bridge' haplotypes. Consider increasing --haplotype-min-sites "
                    f"(currently {min_sites_compared}) to require more overlap for grouping.",
                    UserWarning
                )
        frequencies_true = counts_true.astype(np.float64) / n_haplotypes_win
        sum_freq_sq_true = np.sum(frequencies_true ** 2)
        diversity_true = calculate_haplotype_diversity(window_haplotypes, True)
        frequencies_true_sorted = np.sort(frequencies_true)[::-1]
        top3_freq_true = frequencies_true_sorted[:3] if len(frequencies_true_sorted) >= 3 else list(frequencies_true_sorted) + [0.0] * (3 - len(frequencies_true_sorted))
        
        # Get window size in bp
        if len(window_start_positions) > w_idx and len(window_stop_positions) > w_idx:
            win_size_bp = window_stop_positions[w_idx] - window_start_positions[w_idx] + 1
        else:
            win_size_bp = 0
        
        # Write debug info to file (TSV format)
        with open(debug_file, 'a') as f:
            f.write(f"{run_id}\t{w_idx}\t{win_size_bp}\t{n_variants_win}\t{n_haplotypes_win}\t{missing_rate:.6f}\t"
                   f"{variants_with_missing}\t{haplotypes_with_missing}\t"
                   f"{n_unique_false}\t{n_unique_true}\t{n_singletons_false}\t{n_singletons_true}\t"
                   f"{diversity_false:.6f}\t{diversity_true:.6f}\t{sum_freq_sq_false:.6f}\t{sum_freq_sq_true:.6f}\t"
                   f"{top3_freq_false[0]:.6f}\t{top3_freq_false[1]:.6f}\t{top3_freq_false[2]:.6f}\n")
    
    print(f"[DEBUG] Haplotype diversity diagnostics written to {debug_file}", flush=True)
    
    # Create result dataset
    result = ds.copy()
    result["haplotype_diversity"] = (["windows"], diversity_values)
    
    return result


def garud_h_statistics(ds: xr.Dataset,
                      call_genotype: str = "call_genotype",
                      ignore_missing: bool = False,
                      max_missing: float = 1.0,
                      missing_is_percentage: bool = True,
                      min_sites_compared: int = 1) -> xr.Dataset:
    """
    Calculate Garud H1, H12, H123, and H2/H1 statistics using sgkit's hashing approach.
    
    References:
    - Walsh and Lynch (2018) Equations 9.10-9.13
    - Garud et al. (2015) A selective sweep on the Drosophila X chromosome
    
    Args:
        ds: sgkit Dataset containing genotype calls
        call_genotype: Name of the genotype variable
        ignore_missing: If True, ignore missing data (-1) in hashing
                       If False, include missing data in hash (default sgkit behavior)
        max_missing: Maximum allowed missing data (percentage 0-1 or absolute count)
                    Default: 1.0 (100% missing data allowed - no filtering)
        missing_is_percentage: If True, max_missing is percentage (0-1), 
                              if False, max_missing is absolute count
                              Default: True
        min_sites_compared: Minimum number of non-missing overlapping sites required
                          for two haplotypes to be considered matching (only used when
                          ignore_missing=True). Default is 1.
        
    Returns:
        Dataset with Garud H statistics
    """
    genotypes = ds[call_genotype].values
    n_variants, n_samples, ploidy = genotypes.shape
    
    # Estimate haplotypes assuming phased data (sgkit approach)
    haplotypes = estimate_haplotypes_from_phased_data(genotypes)
    
    # Get window information (windows are always defined)
    n_windows = len(ds.windows)
    window_starts = ds.window_start_idx.values
    window_stops = ds.window_stop_idx.values
    
    # Get window positions for debug output
    positions = ds.variant_position.values
    window_start_positions = []
    window_stop_positions = []
    if 'window_start' in ds.coords:
        window_start_positions = ds.window_start.values
    if 'window_stop' in ds.coords:
        window_stop_positions = ds.window_stop.values
    
    # Calculate Garud H statistics for each window
    h1_values = np.zeros(n_windows)
    h12_values = np.zeros(n_windows)
    h123_values = np.zeros(n_windows)
    h2_h1_values = np.zeros(n_windows)
    
    # Write debug output to file (append mode for multiple runs, TSV format)
    import os
    from datetime import datetime
    debug_file = os.environ.get('PGSTATS_DEBUG_FILE', 'garud_stats_debug.txt')
    run_id = os.environ.get('PGSTATS_RUN_ID', datetime.now().strftime('%Y%m%d_%H%M%S'))
    
    # Write header only if file doesn't exist (TSV format for easy pandas import)
    file_exists = os.path.exists(debug_file)
    with open(debug_file, 'a') as f:
        if not file_exists:
            f.write("Run_ID\tWindow\tWin_Size\tN_Var\tN_Hap\tMissing_Rate\tN_Var_Miss\tN_Hap_Miss\t"
                    "N_Unique\tN_Unique_Ig\tN_Sing\tN_Sing_Ig\t"
                    "H1\tH1_Ig\tH12\tH12_Ig\tH123\tH123_Ig\tH2_H1\tH2_H1_Ig\t"
                    "Top1_Freq\tTop2_Freq\tTop3_Freq\n")
    
    for w_idx in range(n_windows):
        window_start = window_starts[w_idx]
        window_stop = window_stops[w_idx]
        window_haplotypes = haplotypes[window_start:window_stop, :]
        stats = calculate_garud_h_statistics(window_haplotypes, 
                                           ignore_missing, max_missing, missing_is_percentage,
                                           min_sites_compared)
        h1_values[w_idx] = stats[0]
        h12_values[w_idx] = stats[1]
        h123_values[w_idx] = stats[2]
        h2_h1_values[w_idx] = stats[3]
        
        # Calculate debug info
        n_variants_win = window_haplotypes.shape[0]
        n_haplotypes_win = window_haplotypes.shape[1]
        
        # Calculate missing data metrics
        total_calls = n_variants_win * n_haplotypes_win
        missing_calls = np.sum(window_haplotypes == -1)
        missing_rate = missing_calls / total_calls if total_calls > 0 else 0.0
        
        # Count variants with any missing data
        variants_with_missing = np.sum(np.any(window_haplotypes == -1, axis=1))
        # Count haplotypes with any missing data
        haplotypes_with_missing = np.sum(np.any(window_haplotypes == -1, axis=0))
        
        # With ignore_missing=False
        valid_mask_false = filter_haplotypes_by_missingness(window_haplotypes, max_missing, missing_is_percentage)
        valid_haplotypes_false = window_haplotypes[:, valid_mask_false]
        n_valid_false = np.sum(valid_mask_false)
        if n_valid_false >= 2:
            hash_values_false = hash_haplotypes(valid_haplotypes_false, False)
            unique_hashes_false, counts_false = count_unique_values(hash_values_false)
            n_unique_false = len(unique_hashes_false)
            n_singletons_false = np.sum(counts_false == 1)
            frequencies_false = counts_false / n_valid_false
            frequencies_false_sorted = np.sort(frequencies_false)[::-1]
            stats_false = calculate_garud_h_statistics(window_haplotypes, False, max_missing, missing_is_percentage)
            top3_freq_false = frequencies_false_sorted[:3] if len(frequencies_false_sorted) >= 3 else list(frequencies_false_sorted) + [0.0] * (3 - len(frequencies_false_sorted))
        else:
            n_unique_false = 0
            n_singletons_false = 0
            stats_false = np.array([np.nan, np.nan, np.nan, np.nan])
            top3_freq_false = [0.0, 0.0, 0.0]
        
        # With ignore_missing=True
        valid_mask_true = filter_haplotypes_by_missingness(window_haplotypes, max_missing, missing_is_percentage)
        valid_haplotypes_true = window_haplotypes[:, valid_mask_true]
        n_valid_true = np.sum(valid_mask_true)
        if n_valid_true >= 2:
            # Use union-find for grouping (same as in calculate_garud_h_statistics)
            group_ids_true, counts_true = union_find_group_haplotypes(valid_haplotypes_true, min_sites_compared)
            n_unique_true = len(counts_true)
            n_singletons_true = np.sum(counts_true == 1)
            frequencies_true = counts_true.astype(np.float64) / n_valid_true
            frequencies_true_sorted = np.sort(frequencies_true)[::-1]
            stats_true = calculate_garud_h_statistics(window_haplotypes, True, max_missing, 
                                                      missing_is_percentage, min_sites_compared)
            top3_freq_true = frequencies_true_sorted[:3] if len(frequencies_true_sorted) >= 3 else list(frequencies_true_sorted) + [0.0] * (3 - len(frequencies_true_sorted))
            
            # Detect bridge events (only when using ignore_missing)
            if ignore_missing:
                n_bridge_events, bridge_event_list = detect_bridge_events(valid_haplotypes_true, 
                                                                          group_ids_true, 
                                                                          min_sites_compared)
                if n_bridge_events > 0:
                    warnings.warn(
                        f"Window {w_idx}: Detected {n_bridge_events} bridge event(s) where haplotypes "
                        f"are grouped together but don't directly match. This indicates transitive grouping "
                        f"via intermediate 'bridge' haplotypes. Consider increasing --haplotype-min-sites "
                        f"(currently {min_sites_compared}) to require more overlap for grouping.",
                        UserWarning
                    )
        else:
            n_unique_true = 0
            n_singletons_true = 0
            stats_true = np.array([np.nan, np.nan, np.nan, np.nan])
            top3_freq_true = [0.0, 0.0, 0.0]
        
        # Get window size in bp
        if len(window_start_positions) > w_idx and len(window_stop_positions) > w_idx:
            win_size_bp = window_stop_positions[w_idx] - window_start_positions[w_idx] + 1
        else:
            win_size_bp = 0
        
        # Write debug info to file (TSV format)
        with open(debug_file, 'a') as f:
            f.write(f"{run_id}\t{w_idx}\t{win_size_bp}\t{n_variants_win}\t{n_haplotypes_win}\t{missing_rate:.6f}\t"
                   f"{variants_with_missing}\t{haplotypes_with_missing}\t"
                   f"{n_unique_false}\t{n_unique_true}\t{n_singletons_false}\t{n_singletons_true}\t"
                   f"{stats_false[0]:.6f}\t{stats_true[0]:.6f}\t{stats_false[1]:.6f}\t{stats_true[1]:.6f}\t"
                   f"{stats_false[2]:.6f}\t{stats_true[2]:.6f}\t{stats_false[3]:.6f}\t{stats_true[3]:.6f}\t"
                   f"{top3_freq_false[0]:.6f}\t{top3_freq_false[1]:.6f}\t{top3_freq_false[2]:.6f}\n")
    
    print(f"[DEBUG] Garud stats diagnostics written to {debug_file}", flush=True)
    
    # Create result dataset
    result = ds.copy()
    result["garud_h1"] = (["windows"], h1_values)
    result["garud_h12"] = (["windows"], h12_values)
    result["garud_h123"] = (["windows"], h123_values)
    result["garud_h2_h1"] = (["windows"], h2_h1_values)
    
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
        diversity = calculate_haplotype_diversity(haplotypes, ignore_missing=False)
        garud_stats = calculate_garud_h_statistics(haplotypes, ignore_missing=False)
        
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
