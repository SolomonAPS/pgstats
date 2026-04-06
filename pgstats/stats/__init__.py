"""
Statistical functions for population genetics.

This module contains custom statistical functions implemented with Numba
for high-performance computation of population genetics statistics.

Organization:
- sfs_statistics: Site Frequency Spectrum based statistics
- (future) haplotype_statistics: Haplotype-based statistics
- (future) ld_statistics: Linkage disequilibrium statistics
"""

from .sfs_statistics import (
    # Theta estimators
    theta_pi,
    theta_w,
    theta_h,
    theta_l,
    # Neutrality tests
    tajima_d,
    fu_li_d,
    fu_li_f,
    fu_li_d_unfolded,
    fu_li_f_unfolded,
    zeng_e,
    fay_wu_h,
    # SFS components
    singletons,
    windowed_sfs,
)

from .ld_statistics import (
    # Linkage disequilibrium statistics
    calculate_ld_matrix,
    calculate_windowed_ld,
    ld_d,
    ld_d_prime,
    ld_r_squared,
    omega_statistic,
    ld_decay,
)

from .haplotype_statistics import (
    # Haplotype statistics
    haplotype_diversity,
    garud_h_statistics,
    calculate_windowed_haplotype_stats,
    garud_h1,
    garud_h12,
    garud_h123,
    garud_h2_h1,
)

__all__ = [
    # Theta estimators (diversity statistics)
    "theta_pi",
    "theta_w",
    "theta_h",
    "theta_l",
    
    # Neutrality tests
    "tajima_d",
    "fu_li_d",
    "fu_li_f",
    "fu_li_d_unfolded",
    "fu_li_f_unfolded",
    "zeng_e",
    "fay_wu_h",
    
    # SFS components
    "singletons",
    "windowed_sfs",
    
    # Linkage disequilibrium statistics
    "calculate_ld_matrix",
    "calculate_windowed_ld",
    "ld_d",
    "ld_d_prime",
    "ld_r_squared",
    "omega_statistic",
    "ld_decay",
    
    # Haplotype statistics
    "haplotype_diversity",
    "garud_h_statistics",
    "calculate_windowed_haplotype_stats",
    "garud_h1",
    "garud_h12",
    "garud_h123",
    "garud_h2_h1",
]
