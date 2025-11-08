#!/usr/bin/env python3
"""
Apply missingness profile to simulated VCF.
Designed to be fast for parallel processing of many VCFs.

Usage in your pipeline:
    python apply_missingness_profile.py \
        --sim-vcf sim_rep1.vcf.gz \
        --profile missingness_profile.npz \
        --output sim_rep1_matched.vcf.gz \
        --subsample 11

This will:
1. Subsample simulated individuals to match real data sample size
2. Apply exact missingness pattern from real data
3. Output matched VCF ready for pgstats
"""

import argparse
import numpy as np
import allel
import sys
import os


def load_missingness_profile(profile_path):
    """Load pre-computed missingness profile."""
    print(f"Loading missingness profile: {profile_path}", file=sys.stderr)
    data = np.load(profile_path, allow_pickle=True)
    
    profile = {
        'individual_ids': data['individual_ids'],
        'missingness_rate': data['missingness_rate'],
        'positions': data['positions'],
        'n_sites': int(data['n_sites']),
        'n_individuals': int(data['n_individuals']),
        'missing_sites_per_individual': data['missing_sites_per_individual']
    }
    
    print(f"  Real data individuals: {profile['n_individuals']}", file=sys.stderr)
    print(f"  Real data sites: {profile['n_sites']:,}", file=sys.stderr)
    print(f"  Mean missingness: {np.mean(profile['missingness_rate']):.4f}", file=sys.stderr)
    
    return profile


def apply_missingness_to_vcf(sim_vcf_path, profile, output_path, n_subsample=None, seed=None):
    """
    Apply missingness profile to simulated VCF.
    
    Args:
        sim_vcf_path: Path to simulated VCF
        profile: Missingness profile dict from load_missingness_profile()
        output_path: Output VCF path
        n_subsample: Number of individuals to subsample (must be <= profile['n_individuals'])
        seed: Random seed for reproducibility
    """
    if seed is not None:
        np.random.seed(seed)
    
    print(f"\nReading simulated VCF: {sim_vcf_path}", file=sys.stderr)
    
    # Read simulated VCF
    callset = allel.read_vcf(sim_vcf_path, fields=['samples', 'variants/POS', 'variants/CHROM', 
                                                     'variants/REF', 'variants/ALT', 'calldata/GT'])
    
    if callset is None:
        raise ValueError(f"Failed to read VCF: {sim_vcf_path}")
    
    sim_samples = callset['samples']
    sim_positions = callset['variants/POS']
    sim_gt = callset['calldata/GT']  # Shape: [n_variants, n_samples, ploidy]
    
    n_sim_individuals = len(sim_samples)
    n_sim_sites = len(sim_positions)
    
    print(f"  Simulated individuals: {n_sim_individuals}", file=sys.stderr)
    print(f"  Simulated sites: {n_sim_sites:,}", file=sys.stderr)
    
    # Determine subsampling
    n_real = profile['n_individuals']
    if n_subsample is None:
        n_subsample = n_real
    
    if n_subsample > n_real:
        raise ValueError(f"Cannot subsample {n_subsample} individuals from profile with only {n_real}")
    
    if n_subsample > n_sim_individuals:
        raise ValueError(f"Cannot subsample {n_subsample} individuals from sim VCF with only {n_sim_individuals}")
    
    print(f"\nSubsampling to {n_subsample} individuals", file=sys.stderr)
    
    # Randomly select individuals from sim and real
    sim_indices = np.random.choice(n_sim_individuals, size=n_subsample, replace=False)
    real_indices = np.random.choice(n_real, size=n_subsample, replace=False)
    
    print(f"  Sim individuals selected: {sim_indices}", file=sys.stderr)
    print(f"  Real individuals matched: {real_indices}", file=sys.stderr)
    
    # Subsample simulated genotypes
    sim_gt_sub = sim_gt[:, sim_indices, :]  # Shape: [n_sites, n_subsample, ploidy]
    sim_samples_sub = sim_samples[sim_indices]
    
    # Apply missingness pattern
    print(f"\nApplying missingness pattern...", file=sys.stderr)
    
    # Check if site counts match
    if n_sim_sites != profile['n_sites']:
        print(f"  WARNING: Site count mismatch! Sim: {n_sim_sites}, Real: {profile['n_sites']}", file=sys.stderr)
        print(f"  Will apply missingness by position matching", file=sys.stderr)
        
        # Match by position
        position_to_sim_idx = {pos: idx for idx, pos in enumerate(sim_positions)}
        
        for i, real_idx in enumerate(real_indices):
            missing_sites_real = profile['missing_sites_per_individual'][real_idx]
            missing_positions = profile['positions'][missing_sites_real]
            
            # Find corresponding sites in sim data
            for pos in missing_positions:
                if pos in position_to_sim_idx:
                    sim_site_idx = position_to_sim_idx[pos]
                    sim_gt_sub[sim_site_idx, i, :] = -1  # Mark as missing
    else:
        # Direct index matching (faster)
        for i, real_idx in enumerate(real_indices):
            missing_sites = profile['missing_sites_per_individual'][real_idx]
            sim_gt_sub[missing_sites, i, :] = -1  # Mark as missing
    
    # Calculate final missingness
    final_missing_mask = np.any(sim_gt_sub == -1, axis=2)
    final_missingness = np.mean(final_missing_mask, axis=0)
    
    print(f"\nFinal per-individual missingness:", file=sys.stderr)
    for i, (sample, rate) in enumerate(zip(sim_samples_sub, final_missingness)):
        real_sample = profile['individual_ids'][real_indices[i]]
        real_rate = profile['missingness_rate'][real_indices[i]]
        print(f"  {sample} (matched to {real_sample}): {rate:.4f} (target: {real_rate:.4f})", file=sys.stderr)
    
    # Write output VCF
    print(f"\nWriting output VCF: {output_path}", file=sys.stderr)
    
    # Create output dict for allel.write_vcf
    output_data = {
        'samples': sim_samples_sub,
        'variants/CHROM': callset['variants/CHROM'],
        'variants/POS': sim_positions,
        'variants/REF': callset['variants/REF'],
        'variants/ALT': callset['variants/ALT'],
        'calldata/GT': sim_gt_sub
    }
    
    allel.write_vcf(output_path, output_data, fill='.')
    
    print(f"Done! Output written to: {output_path}", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(
        description='Apply missingness profile to simulated VCF',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example:
    python apply_missingness_profile.py \\
        --sim-vcf sim_rep1.vcf.gz \\
        --profile missingness_profile.npz \\
        --output sim_rep1_matched.vcf.gz \\
        --subsample 11 \\
        --seed 42
        """
    )
    
    parser.add_argument('--sim-vcf', required=True,
                        help='Input simulated VCF file')
    parser.add_argument('--profile', required=True,
                        help='Missingness profile NPZ file')
    parser.add_argument('--output', required=True,
                        help='Output VCF file with matched missingness')
    parser.add_argument('--subsample', type=int, default=None,
                        help='Number of individuals to subsample (default: same as profile)')
    parser.add_argument('--seed', type=int, default=None,
                        help='Random seed for reproducibility')
    
    args = parser.parse_args()
    
    # Load profile
    profile = load_missingness_profile(args.profile)
    
    # Apply to simulated VCF
    apply_missingness_to_vcf(
        args.sim_vcf,
        profile,
        args.output,
        n_subsample=args.subsample,
        seed=args.seed
    )


if __name__ == '__main__':
    main()

