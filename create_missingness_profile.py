#!/usr/bin/env python3
"""
Create a missingness profile from real data VCF.
This profile can be rapidly applied to simulated data.

Usage:
    python create_missingness_profile.py \
        --vcf real_data.vcf.gz \
        --output missingness_profile.npz \
        --region chr2L:1000000-2000000

Output NPZ contains:
    - individual_ids: List of individual names
    - missingness_rate: Per-individual missingness rate [n_individuals]
    - missing_sites: List of site indices where each individual is missing [n_individuals, variable_length]
    - positions: Array of all variant positions
    - n_sites: Total number of sites
"""

import argparse
import numpy as np
import allel
import sys


def calculate_missingness_profile(vcf_path, region=None):
    """
    Calculate per-individual missingness profile from VCF.
    
    Returns:
        dict with keys:
            - individual_ids: list of sample names
            - missingness_rate: array of missingness rates per individual
            - missing_mask: boolean array [n_sites, n_individuals] where True = missing
            - positions: array of variant positions
            - contig: chromosome name
    """
    print(f"Reading VCF: {vcf_path}", file=sys.stderr)
    if region:
        print(f"  Region: {region}", file=sys.stderr)
    
    # Read VCF
    if region:
        callset = allel.read_vcf(vcf_path, region=region, fields=['samples', 'variants/POS', 'calldata/GT'])
    else:
        callset = allel.read_vcf(vcf_path, fields=['samples', 'variants/POS', 'calldata/GT'])
    
    if callset is None:
        raise ValueError(f"Failed to read VCF: {vcf_path}")
    
    samples = callset['samples']
    positions = callset['variants/POS']
    gt = callset['calldata/GT']  # Shape: [n_variants, n_samples, ploidy]
    
    n_sites = len(positions)
    n_individuals = len(samples)
    
    print(f"  Sites: {n_sites:,}", file=sys.stderr)
    print(f"  Individuals: {n_individuals}", file=sys.stderr)
    
    # Create missing mask: True where genotype is missing (-1)
    # For diploid: missing if either allele is -1
    missing_mask = np.any(gt == -1, axis=2)  # Shape: [n_sites, n_individuals]
    
    # Calculate per-individual missingness rate
    missingness_rate = np.mean(missing_mask, axis=0)  # Shape: [n_individuals]
    
    print("\nPer-individual missingness:", file=sys.stderr)
    for i, (sample, rate) in enumerate(zip(samples, missingness_rate)):
        print(f"  {sample}: {rate:.4f} ({int(rate * n_sites):,} / {n_sites:,} sites)", file=sys.stderr)
    
    print(f"\nOverall missingness: {np.mean(missingness_rate):.4f}", file=sys.stderr)
    
    # Extract chromosome name from region or VCF
    if region:
        contig = region.split(':')[0]
    else:
        # Try to get from VCF header or use first position's contig
        contig = "unknown"
    
    return {
        'individual_ids': samples,
        'missingness_rate': missingness_rate,
        'missing_mask': missing_mask,
        'positions': positions,
        'contig': contig,
        'n_sites': n_sites,
        'n_individuals': n_individuals
    }


def save_missingness_profile(profile, output_path):
    """Save missingness profile to compressed NPZ file."""
    print(f"\nSaving missingness profile to: {output_path}", file=sys.stderr)
    
    # Convert missing_mask to sparse format to save space
    # Store indices of missing sites for each individual
    missing_sites_per_individual = []
    for ind_idx in range(profile['n_individuals']):
        missing_sites = np.where(profile['missing_mask'][:, ind_idx])[0]
        missing_sites_per_individual.append(missing_sites)
    
    np.savez_compressed(
        output_path,
        individual_ids=profile['individual_ids'],
        missingness_rate=profile['missingness_rate'],
        positions=profile['positions'],
        contig=profile['contig'],
        n_sites=profile['n_sites'],
        n_individuals=profile['n_individuals'],
        # Save missing sites as object array (variable length arrays)
        missing_sites_per_individual=np.array(missing_sites_per_individual, dtype=object)
    )
    
    print(f"  File size: {np.os.path.getsize(output_path) / 1024 / 1024:.2f} MB", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(
        description='Create missingness profile from real data VCF',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Whole VCF
    python create_missingness_profile.py \\
        --vcf real_data.vcf.gz \\
        --output missingness_profile.npz
    
    # Specific region
    python create_missingness_profile.py \\
        --vcf real_data.vcf.gz \\
        --region chr2L:1000000-2000000 \\
        --output missingness_profile_chr2L_1-2Mb.npz
        """
    )
    
    parser.add_argument('--vcf', required=True,
                        help='Input VCF file (real data)')
    parser.add_argument('--output', required=True,
                        help='Output NPZ file with missingness profile')
    parser.add_argument('--region', default=None,
                        help='Genomic region (format: chr:start-end)')
    
    args = parser.parse_args()
    
    # Calculate profile
    profile = calculate_missingness_profile(args.vcf, args.region)
    
    # Save to file
    save_missingness_profile(profile, args.output)
    
    print("\nDone!", file=sys.stderr)


if __name__ == '__main__':
    main()

