#!/usr/bin/env python3
"""
Test FIXED PyRanges implementation of callable length calculation.
"""

import pandas as pd
import numpy as np
import pyranges as pr

# Get BED file path from previous run
bed_path = "/var/folders/v3/ybcp76dx0bx9gy8dznldx29c0000gp/T/tmp4p5mtw17.bed"

def calculate_callable_length_pyranges_fixed(window_start, window_stop, contig_name, bed_gr, bed_format="non_callable"):
    """
    FIXED PyRanges implementation of callable length calculation.
    """
    window_length = window_stop - window_start
    
    if bed_format == "non_callable":
        # BED defines non-callable regions
        # Create window PyRanges object
        window_gr = pr.PyRanges(pd.DataFrame({
            'Chromosome': [contig_name],
            'Start': [window_start],
            'End': [window_stop]
        }))
        
        # Find overlaps with BED regions
        overlaps = window_gr.overlap(bed_gr)
        
        if len(overlaps) > 0:
            # Calculate total masked length by finding intersections
            masked_length = 0
            
            # Get the overlapping BED regions
            for _, overlap_row in overlaps.df.iterrows():
                # The overlap_row contains the window coordinates, not the intersection
                # We need to find the actual BED regions that overlap and calculate intersections
                pass
            
            # Alternative approach: manually find overlapping BED regions
            bed_df = bed_gr.df
            overlapping_bed = bed_df[bed_df['Chromosome'] == contig_name]
            
            for _, bed_row in overlapping_bed.iterrows():
                bed_start = bed_row['Start']
                bed_end = bed_row['End']
                
                # Calculate intersection
                intersection_start = max(window_start, bed_start)
                intersection_end = min(window_stop, bed_end)
                
                if intersection_start < intersection_end:
                    masked_length += intersection_end - intersection_start
            
            callable_length = window_length - masked_length
        else:
            # No overlaps, entire window is callable
            callable_length = window_length
            
    elif bed_format == "callable":
        # BED defines callable regions
        # Create window PyRanges object
        window_gr = pr.PyRanges(pd.DataFrame({
            'Chromosome': [contig_name],
            'Start': [window_start],
            'End': [window_stop]
        }))
        
        # Find overlaps with BED regions
        overlaps = window_gr.overlap(bed_gr)
        
        if len(overlaps) > 0:
            # Calculate total callable length by finding intersections
            callable_length = 0
            
            # Get the overlapping BED regions
            bed_df = bed_gr.df
            overlapping_bed = bed_df[bed_df['Chromosome'] == contig_name]
            
            for _, bed_row in overlapping_bed.iterrows():
                bed_start = bed_row['Start']
                bed_end = bed_row['End']
                
                # Calculate intersection
                intersection_start = max(window_start, bed_start)
                intersection_end = min(window_stop, bed_end)
                
                if intersection_start < intersection_end:
                    callable_length += intersection_end - intersection_start
        else:
            # No overlaps, no callable sites
            callable_length = 0
    
    # Ensure at least 1 callable site to avoid division by zero
    return max(1, callable_length)

# Load BED file and convert to PyRanges
bed_df = pd.read_csv(bed_path, sep='\t', header=None, names=['chrom', 'start', 'end'])
bed_gr = pr.PyRanges(bed_df.rename(columns={'chrom': 'Chromosome', 'start': 'Start', 'end': 'End'}))

print("BED file loaded and converted to PyRanges:")
print(bed_gr)

# Load expected results from corrected implementation
expected_results = pd.read_csv('corrected_implementation_results_fixed.csv')
print(f"\nLoaded {len(expected_results)} expected results from corrected implementation")

# Test PyRanges implementation
print(f"\nTesting FIXED PyRanges implementation:")
print("Window\t\tContig\tExpected\tActual\tMatch")
print("-" * 60)

results = []
for _, row in expected_results.iterrows():
    window_start = int(row['window_start'])
    window_stop = int(row['window_stop'])
    contig_name = row['contig_name']
    expected = int(row['expected'])
    
    actual = calculate_callable_length_pyranges_fixed(window_start, window_stop, contig_name, bed_gr)
    match = "✓" if actual == expected else "✗"
    results.append({
        'window_start': window_start,
        'window_stop': window_stop,
        'contig_name': contig_name,
        'expected': expected,
        'actual': actual,
        'match': actual == expected
    })
    print(f"[{window_start:3d}, {window_stop:3d})\t{contig_name}\t{expected:3d}bp\t\t{actual:3d}bp\t{match}")

# Summary
total_tests = len(results)
passed_tests = sum(1 for r in results if r['match'])
print(f"\nSummary: {passed_tests}/{total_tests} tests passed")

if passed_tests == total_tests:
    print("✓ All tests passed! FIXED PyRanges implementation matches corrected implementation.")
    print("Ready to integrate into the main code!")
else:
    print("✗ Some tests failed. Need to debug PyRanges implementation further.")

# Save results for comparison
results_df = pd.DataFrame(results)
results_df.to_csv('pyranges_fixed_results.csv', index=False)
print("Results saved to: pyranges_fixed_results.csv")
