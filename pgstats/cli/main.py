#!/usr/bin/env python3
"""
Main command-line interface for pgstats.

This module provides the main CLI entry point for calculating population genetics
statistics from VCF files.
"""

import sys
import argparse
import logging
from pathlib import Path
from typing import List, Tuple, Optional

from pgstats.core.dataset import GenomicDataset, WindowConfig, CallableSitesConfig
from pgstats import __version__


def setup_logging(verbose: int = 0):
    """
    Setup logging configuration.
    
    Args:
        verbose: Verbosity level (0=WARNING, 1=INFO, 2=DEBUG)
    """
    level = logging.WARNING
    if verbose == 1:
        level = logging.INFO
    elif verbose >= 2:
        level = logging.DEBUG
    
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )


def create_parser() -> argparse.ArgumentParser:
    """Create the main argument parser."""
    
    parser = argparse.ArgumentParser(
        prog='pgstats',
        description='Calculate population genetics statistics from VCF files',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Calculate genome-wide statistics
  pgstats stats input.vcf.gz --output results.csv
  
  # Windowed analysis with 100kb windows
  pgstats stats input.vcf.gz --output results.csv \\
      --window-size 100000 --step-size 50000
  
  # With callable sites masking
  pgstats stats input.vcf.gz --output results.csv \\
      --bed callable_sites.bed --max-missing 0.1
  
  # Region-specific analysis
  pgstats stats input.vcf.gz --output results.csv \\
      --region chr1:1000000-5000000 --window-size 50000
  
  # Multi-region analysis
  pgstats stats input.vcf.gz --output results.csv \\
      --regions-file candidate_genes.bed --window-size 10000
  
  # Multi-region with callable sites masking
  pgstats stats input.vcf.gz --output results.csv \\
      --bed callable_sites.bed --regions-file target_regions.bed \\
      --window-size 50000 --max-missing 0.1
  
  # With non-callable regions (default behavior)
  pgstats stats input.vcf.gz --output results.csv \\
      --bed repetitive_regions.bed --bed-format non_callable
  
  # With callable regions only
  pgstats stats input.vcf.gz --output results.csv \\
      --bed high_quality_regions.bed --bed-format callable
  
  # Keep Zarr files for faster re-runs (saved next to output file)
  pgstats stats input.vcf.gz --output results.csv \\
      --keep-zarr --window-size 100000
  
  # Specify custom Zarr directory (useful for clusters)
  pgstats stats input.vcf.gz --output results.csv \\
      --keep-zarr --zarr-dir /work/users/s/o/solsloat/zarr_cache

For more information, visit: https://github.com/SolomonAPS/pgstats
        """
    )
    
    parser.add_argument(
        '--version',
        action='version',
        version=f'%(prog)s {__version__}'
    )
    
    parser.add_argument(
        '-v', '--verbose',
        action='count',
        default=0,
        help='Increase verbosity (can be used multiple times: -v, -vv)'
    )
    
    subparsers = parser.add_subparsers(
        dest='command',
        title='commands',
        description='Available commands'
    )
    
    # Stats command
    stats_parser = subparsers.add_parser(
        'stats',
        help='Calculate population genetics statistics',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description='Calculate various population genetics statistics from VCF files'
    )
    add_stats_arguments(stats_parser)
    
    # Info command
    info_parser = subparsers.add_parser(
        'info',
        help='Display dataset information',
        description='Display information about a VCF file'
    )
    add_info_arguments(info_parser)
    
    return parser


def add_stats_arguments(parser: argparse.ArgumentParser):
    """Add arguments for the stats command."""
    
    # Required arguments
    required = parser.add_argument_group('required arguments')
    required.add_argument(
        'input',
        type=str,
        help='Input VCF file (.vcf, .vcf.gz, or .zarr)'
    )
    required.add_argument(
        '-o', '--output',
        type=str,
        required=True,
        help='Output file path'
    )
    
    # Statistics selection
    stats_group = parser.add_argument_group('statistics options')
    stats_group.add_argument(
        '-s', '--stats',
        type=str,
        nargs='+',
        default=['all'],
        choices=[
            # Neutrality tests
            'tajima_d', 'fu_li_d', 'fu_li_f', 'fu_li_d_unfolded', 'fu_li_f_unfolded', 'zeng_e', 'fay_wu_h',
            # Theta estimators (primary names)
            'theta_pi', 'theta_w', 'theta_h', 'theta_l',
            # Backward compatibility aliases
            'pi', 'nucleotide_diversity', 'watterson_theta', 'fay_wu_theta',
            # LD statistics
            'ld_d', 'ld_dprime', 'ld_r2', 'omega_statistic',
            # Haplotype statistics
            'haplotype_diversity', 'garud_h1', 'garud_h12', 'garud_h123', 'garud_h2_h1',
            # Singleton counts
            'singletons', 'singletons_unfolded',
            # Special option
            'all'
        ],
        help='Statistics to calculate (default: all)'
    )
    
    # Windowing options
    window_group = parser.add_argument_group('windowing options')
    window_group.add_argument(
        '-w', '--window-size',
        type=int,
        default=None,
        metavar='INT',
        help='Window size in base pairs (default: genome-wide)'
    )
    window_group.add_argument(
        '--step-size',
        type=int,
        default=None,
        metavar='INT',
        help='Step size for sliding windows (default: same as window-size)'
    )
    window_group.add_argument(
        '--min-variants',
        type=int,
        default=5,
        metavar='INT',
        help='Minimum number of variants per window (default: 5)'
    )
    # Region selection (mutually exclusive)
    region_group = parser.add_mutually_exclusive_group()
    region_group.add_argument(
        '-r', '--region',
        type=str,
        metavar='REGION',
        help='Single genomic region to analyze (format: chr:start-end or chr)'
    )
    region_group.add_argument(
        '--regions-file',
        type=str,
        metavar='FILE',
        help='BED file with multiple regions to analyze (one region per line)'
    )
    
    # Callable sites and filtering
    filter_group = parser.add_argument_group('filtering options')
    filter_group.add_argument(
        '-b', '--bed',
        type=str,
        metavar='FILE',
        help='BED file with callable sites (see --bed-format for interpretation)'
    )
    filter_group.add_argument(
        '--bed-format',
        type=str,
        choices=['callable', 'non_callable'],
        default='non_callable',
        help='BED file format: "callable" (regions to include) or "non_callable" (regions to exclude) (default: non_callable)'
    )
    filter_group.add_argument(
        '--max-missing',
        type=float,
        default=0.0,
        metavar='FLOAT',
        help='Maximum proportion of missing data per variant (default: 0.0 = no filtering; e.g., 0.2 = remove variants with >20%% missing)'
    )
    filter_group.add_argument(
        '--fixed-sample-size',
        action='store_true',
        help='Use fixed maximum sample size for theta_pi calculation (like SFS-based stats). Default uses per-site sample sizes.'
    )
    filter_group.add_argument(
        '--haplotype-ignore-missing',
        action='store_true',
        help='Ignore missing data when hashing haplotypes for Garud statistics. If True, haplotypes that match at all non-missing positions are grouped together. If False (default), missing data is included in the hash.'
    )
    filter_group.add_argument(
        '--keep-zarr',
        action='store_true',
        help='Keep intermediate Zarr files for faster re-runs. Files saved next to output file by default.'
    )
    filter_group.add_argument(
        '--zarr-dir',
        type=str,
        metavar='DIR',
        help='Directory to store Zarr files (overrides default: next to output file with --keep-zarr)'
    )
    
    # Output options
    output_group = parser.add_argument_group('output options')
    output_group.add_argument(
        '-f', '--format',
        type=str,
        default='csv',
        choices=['csv', 'tsv', 'parquet'],
        help='Output format (default: csv)'
    )
    output_group.add_argument(
        '--no-header',
        action='store_true',
        help='Do not write header row in output'
    )
    output_group.add_argument(
        '--profile',
        action='store_true',
        help='Enable detailed timing output for performance monitoring'
    )


def add_info_arguments(parser: argparse.ArgumentParser):
    """Add arguments for the info command."""
    
    parser.add_argument(
        'input',
        type=str,
        help='Input VCF file (.vcf, .vcf.gz, or .zarr)'
    )
    parser.add_argument(
        '-b', '--bed',
        type=str,
        metavar='FILE',
        help='BED file to check callable sites'
    )


def parse_regions_file(regions_file: str) -> List[Tuple[str, int, int]]:
    """
    Parse BED file with regions to analyze.
    
    Args:
        regions_file: Path to BED file with regions
        
    Returns:
        List of (contig, start, end) tuples (1-based, inclusive)
        
    Raises:
        ValueError: If file format is invalid
        FileNotFoundError: If file doesn't exist
    """
    import pandas as pd
    
    bed_path = Path(regions_file)
    if not bed_path.exists():
        raise FileNotFoundError(f"Regions file not found: {bed_path}")
    
    try:
        # Read BED file (0-based coordinates)
        bed_df = pd.read_csv(
            bed_path,
            sep='\t',
            names=['chrom', 'start', 'end'],
            usecols=[0, 1, 2],
            comment='#'
        )
        
        if bed_df.empty:
            raise ValueError("Regions file is empty")
        
        # Convert to 1-based coordinates for internal use
        regions = []
        for _, row in bed_df.iterrows():
            contig = str(row['chrom'])
            start = int(row['start']) + 1  # Convert 0-based to 1-based
            end = int(row['end'])          # BED end is exclusive, keep as-is
            
            if start > end:
                raise ValueError(f"Invalid region: {contig}:{start}-{end} (start > end)")
            
            regions.append((contig, start, end))
        
        return regions
        
    except Exception as e:
        raise ValueError(f"Error parsing regions file: {e}")


def parse_region(region_str: str) -> tuple:
    """
    Parse region string into contig, start, end.
    
    Args:
        region_str: Region string (e.g., "chr1:1000-5000" or "chr1")
        
    Returns:
        Tuple of (contig, start, end) or (contig, None, None)
    """
    if ':' not in region_str:
        return (region_str, None, None)
    
    contig, positions = region_str.split(':', 1)
    
    if '-' not in positions:
        raise ValueError(f"Invalid region format: {region_str}. Expected chr:start-end")
    
    start_str, end_str = positions.split('-', 1)
    start = int(start_str.replace(',', ''))
    end = int(end_str.replace(',', ''))
    
    return (contig, start, end)


def run_stats_command(args):
    """Execute the stats command."""
    
    logger = logging.getLogger(__name__)
    
    # Validate input file
    input_path = Path(args.input)
    if not input_path.exists():
        logger.error(f"Input file not found: {input_path}")
        sys.exit(1)
    
    # Validate BED file if provided
    if args.bed:
        bed_path = Path(args.bed)
        if not bed_path.exists():
            logger.error(f"BED file not found: {bed_path}")
            sys.exit(1)
    
    # Parse region(s) if provided
    region_contig = None
    region_start = None
    region_end = None
    regions_list = None
    
    if args.region:
        try:
            region_contig, region_start, region_end = parse_region(args.region)
            logger.info(f"Analyzing single region: {region_contig}:{region_start}-{region_end}")
        except ValueError as e:
            logger.error(f"Invalid region format: {e}")
            sys.exit(1)
    
    elif args.regions_file:
        try:
            regions_list = parse_regions_file(args.regions_file)
            logger.info(f"Analyzing {len(regions_list)} regions from file: {args.regions_file}")
            for i, (contig, start, end) in enumerate(regions_list[:3], 1):  # Show first 3
                logger.info(f"  Region {i}: {contig}:{start:,}-{end:,}")
            if len(regions_list) > 3:
                logger.info(f"  ... and {len(regions_list) - 3} more regions")
        except (ValueError, FileNotFoundError) as e:
            logger.error(f"Error with regions file: {e}")
            sys.exit(1)
    
    # Configure windowing
    window_config = WindowConfig(
        window_size=args.window_size,
        step_size=args.step_size,
        start=region_start,
        end=region_end,
        min_variants=args.min_variants
    )
    
    # Configure callable sites
    callable_config = CallableSitesConfig(
        bed_file=args.bed,
        bed_format=args.bed_format,
        max_missing=args.max_missing
    )
    
    # Extract output directory for Zarr file placement
    output_dir = None
    if args.output:
        output_dir = str(Path(args.output).parent)
    
    # Load dataset
    logger.info(f"Loading data from {input_path}")
    print(f"Loading VCF: {input_path}", flush=True)
    
    try:
        # Determine if we're doing region-based analysis (use lazy masking)
        # Use lazy masking for both single region (--region) and multi-region (--regions-file)
        use_lazy_masking = (region_contig is not None) or (regions_list is not None)
        
        # Temporarily disable BED loading in __init__ if we'll use lazy mode
        callable_config_for_init = callable_config
        if use_lazy_masking and callable_config.bed_file:
            # Create a config without BED file for __init__, we'll load it lazily after
            callable_config_for_init = CallableSitesConfig(
                bed_file=None,  # Don't load in __init__
                bed_format=callable_config.bed_format,
                max_missing=callable_config.max_missing
            )
        
        genomic_ds = GenomicDataset(
            data_source=str(input_path),
            callable_config=callable_config_for_init,
            window_config=window_config,
            keep_zarr=args.keep_zarr,
            zarr_dir=args.zarr_dir,
            output_dir=output_dir,
            enable_profiling=args.profile
        )
        
        # If using BED file with region-based analysis, use lazy loading
        if use_lazy_masking and callable_config.bed_file:
            # Restore the original callable config with BED file
            genomic_ds.callable_config = callable_config
            # Load BED in lazy mode (don't mask genome-wide, mask each region individually)
            genomic_ds._load_callable_sites(lazy_mode=True)
        
    except Exception as e:
        logger.error(f"Failed to load dataset: {e}")
        sys.exit(1)
    
    # Print dataset summary
    summary = genomic_ds.get_summary()
    print(f"\nDataset Summary:", flush=True)
    print(f"  Variants: {summary['n_variants']:,}", flush=True)
    print(f"  Samples: {summary['n_samples']:,}", flush=True)
    print(f"  Contigs: {summary['n_contigs']:,}", flush=True)
    
    # Filter missing data if needed
    if args.max_missing > 0:
        logger.info(f"Filtering variants with >{args.max_missing:.1%} missing data")
        print(f"\nFiltering variants with >{args.max_missing:.1%} missing data...", flush=True)
        genomic_ds.filter_missing_data()
        summary = genomic_ds.get_summary()
        print(f"  Retained variants: {summary['n_variants']:,}", flush=True)
    
    # Calculate statistics
    # Handle "all" option
    if args.stats == ['all']:
        stats_to_calculate = [
            # Neutrality tests
            'tajima_d', 'fu_li_d', 'fu_li_f', 'fu_li_d_unfolded', 'fu_li_f_unfolded', 'zeng_e', 'fay_wu_h',
            # Theta estimators
            'theta_pi', 'theta_w', 'theta_h', 'theta_l',
            # LD statistics
            'ld_d', 'ld_dprime', 'ld_r2', 'omega_statistic',
            # Haplotype statistics
            'haplotype_diversity', 'garud_h1', 'garud_h12', 'garud_h123', 'garud_h2_h1',
            # Singleton counts
            'singletons', 'singletons_unfolded'
        ]
    else:
        stats_to_calculate = args.stats
    
    logger.info(f"Calculating statistics: {', '.join(stats_to_calculate)}")
    print(f"\nCalculating statistics: {', '.join(stats_to_calculate)}", flush=True)
    
    try:
        if regions_list:
            # Multi-region analysis
            print(f"  Window size: {args.window_size:,} bp")
            print(f"  Step size: {args.step_size or args.window_size:,} bp")
            print(f"  Minimum variants: {args.min_variants}")
            print(f"\nStarting statistics calculation across {len(regions_list)} regions...")
            
            results_df = genomic_ds.calculate_stats_for_regions(
                regions=regions_list,
                window_size=args.window_size,
                step_size=args.step_size,
                stats=stats_to_calculate,
                min_variants=args.min_variants,
                use_callable_sites=True,
                use_fixed_n=args.fixed_sample_size,
                haplotype_ignore_missing=args.haplotype_ignore_missing
            )
            
            print(f"  Calculated statistics for {len(results_df)} total windows")
            
            # Save results
            logger.info(f"Saving results to {args.output}")
            print(f"\nSaving results to {args.output}...")
            results_df.to_csv(args.output, index=False)
            print(f"Results saved successfully")
            
        elif args.window_size or args.region:
            # Single region or genome-wide windowed analysis
            print(f"\nCreating windows...")
            genomic_ds.create_windows()
            summary = genomic_ds.get_summary()
            
            if args.window_size:
                print(f"  Window size: {args.window_size:,} bp")
                print(f"  Step size: {args.step_size or args.window_size:,} bp")
                print(f"  Minimum variants: {args.min_variants}")
                print(f"  Number of windows: {summary['n_windows']:,}")
            else:
                print(f"  Genome-wide analysis")
            
            print(f"\nStarting statistics calculation...")
            window_stats = genomic_ds.calculate_windowed_stats(
                stats=stats_to_calculate,
                use_callable_sites=True,
                use_fixed_n=args.fixed_sample_size,
                haplotype_ignore_missing=args.haplotype_ignore_missing
            )
            print(f"  Calculated statistics for {len(window_stats.windows)} windows")
            
            # Save results
            logger.info(f"Saving results to {args.output}")
            print(f"\nSaving results to {args.output}...")
            genomic_ds.save_results(args.output, format=args.format)
            print(f"Results saved successfully")
        else:
            # Genome-wide analysis
            print(f"\nStarting genome-wide statistics calculation...")
            genome_stats = genomic_ds.calculate_genome_wide_stats(stats=stats_to_calculate)
            print(f"  Genome-wide results:")
            for stat, value in genome_stats.items():
                print(f"    {stat}: {value:.6f}")
            
            # Save genome-wide results
            logger.info(f"Saving results to {args.output}")
            print(f"\nSaving results to {args.output}...")
            import pandas as pd
            df = pd.DataFrame([genome_stats])
            if args.format == 'csv':
                df.to_csv(args.output, index=False, header=not args.no_header)
            elif args.format == 'tsv':
                df.to_csv(args.output, sep='\t', index=False, header=not args.no_header)
            elif args.format == 'parquet':
                df.to_parquet(args.output, index=False)
            print(f"Results saved successfully")
    except Exception as e:
        logger.error(f"Failed to calculate/save statistics: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
    
    print(f"\nAnalysis complete!")


def run_info_command(args):
    """Execute the info command."""
    
    logger = logging.getLogger(__name__)
    
    # Validate input file
    input_path = Path(args.input)
    if not input_path.exists():
        logger.error(f"Input file not found: {input_path}")
        sys.exit(1)
    
    # Configure callable sites if BED provided
    callable_config = None
    if args.bed:
        bed_path = Path(args.bed)
        if not bed_path.exists():
            logger.error(f"BED file not found: {bed_path}")
            sys.exit(1)
        callable_config = CallableSitesConfig(bed_file=args.bed, bed_format=args.bed_format)
    
    # Extract output directory for Zarr file placement
    output_dir = None
    if args.output:
        output_dir = str(Path(args.output).parent)
    
    # Load dataset
    print(f"Loading dataset: {input_path}")
    
    try:
        genomic_ds = GenomicDataset(
            data_source=str(input_path),
            callable_config=callable_config,
            keep_zarr=args.keep_zarr,
            zarr_dir=args.zarr_dir,
            output_dir=output_dir
        )
    except Exception as e:
        logger.error(f"Failed to load dataset: {e}")
        sys.exit(1)
    
    # Get and display summary
    summary = genomic_ds.get_summary()
    
    print(f"\n{'='*60}")
    print(f"Dataset Information")
    print(f"{'='*60}")
    print(f"File: {input_path.name}")
    print(f"\nDimensions:")
    print(f"  Variants: {summary['n_variants']:,}")
    print(f"  Samples: {summary['n_samples']:,}")
    print(f"  Contigs: {summary['n_contigs']:,}")
    
    # Check for missing data
    from pgstats.utils.validation import check_missing_data
    missing_info = check_missing_data(genomic_ds.dataset)
    
    print(f"\nMissing Data:")
    print(f"  Total missing calls: {missing_info['total_missing']:,}")
    print(f"  Missing rate: {missing_info['missing_rate']:.3%}")
    print(f"  Variants with missing: {missing_info['variants_with_missing']:,}")
    print(f"  Samples with missing: {missing_info['samples_with_missing']:,}")
    
    print(f"\n{'='*60}")


def main(argv: Optional[List[str]] = None):
    """
    Main entry point for the CLI.
    
    Args:
        argv: Command-line arguments (default: sys.argv[1:])
    """
    parser = create_parser()
    
    if argv is None:
        argv = sys.argv[1:]
    
    # Show help if no arguments provided
    if len(argv) == 0:
        parser.print_help()
        sys.exit(0)
    
    args = parser.parse_args(argv)
    
    # Setup logging
    setup_logging(args.verbose)
    
    # Execute command
    if args.command == 'stats':
        run_stats_command(args)
    elif args.command == 'info':
        run_info_command(args)
    else:
        parser.print_help()
        sys.exit(1)


if __name__ == '__main__':
    main()

