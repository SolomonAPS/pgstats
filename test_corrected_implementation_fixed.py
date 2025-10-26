#!/usr/bin/env python3
"""
Test corrected callable length implementation with FIXED expectations.
"""

import pandas as pd
import numpy as np

# Get BED file path from previous run
bed_path = "/var/folders/v3/ybcp76dx0bx9gy8dznldx29c0000gp/T/tmp4p5mtw17.bed"

def calculate_callable_length_corrected(window_start, window_stop, contig_name, bed_df, bed_format="non_callable"):
    """Corrected implementation that properly handles contig filtering."""
    window_length = window_stop - window_start
    
    if bed_format == "non_callable":
        # BED defines non-callable regions
        masked_length = 0
        
        for _, row in bed_df.iterrows():
            chrom = row['chrom']
            bed_start = row['start']
            bed_end = row['end']
            
            # Only process BED regions on the same contig
            if str(chrom) == str(contig_name):
                # Calculate overlap between BED region and window
                overlap_start = max(window_start, bed_start)
                overlap_end = min(window_stop, bed_end)
                
                if overlap_start < overlap_end:
                    masked_length += overlap_end - overlap_start
        
        # Calculate callable length
        callable_length = window_length - masked_length
        
    elif bed_format == "callable":
        # BED defines callable regions
        callable_length = 0
        
        for _, row in bed_df.iterrows():
            chrom = row['chrom']
            bed_start = row['start']
            bed_end = row['end']
            
            # Only process BED regions on the same contig
            if str(chrom) == str(contig_name):
                # Calculate overlap between BED region and window
                overlap_start = max(window_start, bed_start)
                overlap_end = min(window_stop, bed_end)
                
                if overlap_start < overlap_end:
                    callable_length += overlap_end - overlap_start
    
    # Ensure at least 1 callable site to avoid division by zero
    return max(1, callable_length)

# Load BED file
bed_df = pd.read_csv(bed_path, sep='\t', header=None, names=['chrom', 'start', 'end'])
print("BED file loaded:")
print(bed_df)

# CORRECTED test cases with proper expectations
test_cases = [
    # chr1 tests
    (50, 150, 'chr1', 50),    # Window [50,150): BED [100,200) → overlap [100,150) = 50bp masked → 100-50=50 callable
    (150, 250, 'chr1', 50),   # Window [150,250): BED [100,200) → overlap [150,200) = 50bp masked → 100-50=50 callable
    (250, 350, 'chr1', 50),   # Window [250,350): BED [300,400) → overlap [300,350) = 50bp masked → 100-50=50 callable
    (350, 450, 'chr1', 50),   # Window [350,450): BED [300,400) → overlap [350,400) = 50bp masked → 100-50=50 callable
    (450, 550, 'chr1', 50),   # Window [450,550): BED [500,600) → overlap [500,550) = 50bp masked → 100-50=50 callable
    (550, 650, 'chr1', 50),   # Window [550,650): BED [500,600) → overlap [550,600) = 50bp masked → 100-50=50 callable
    (650, 750, 'chr1', 100),  # Window [650,750): No BED overlap → 100bp callable
    
    # chr2 tests
    (0, 100, 'chr2', 50),     # Window [0,100): BED [50,150) → overlap [50,100) = 50bp masked → 100-50=50 callable
    (100, 200, 'chr2', 50),   # Window [100,200): BED [50,150) → overlap [100,150) = 50bp masked → 100-50=50 callable
    (200, 300, 'chr2', 50),   # Window [200,300): BED [250,350) → overlap [250,300) = 50bp masked → 100-50=50 callable
    (300, 400, 'chr2', 50),   # Window [300,400): BED [250,350) → overlap [300,350) = 50bp masked → 100-50=50 callable
    (400, 500, 'chr2', 100),  # Window [400,500): No BED overlap → 100bp callable
    
    # Edge cases
    (95, 105, 'chr1', 5),     # Window [95,105): BED [100,200) → overlap [100,105) = 5bp masked → 10-5=5 callable
    (195, 205, 'chr1', 5),    # Window [195,205): BED [100,200) → overlap [195,200) = 5bp masked → 10-5=5 callable
]

print(f"\nTesting corrected implementation with FIXED expectations:")
print("Window\t\tContig\tExpected\tActual\tMatch")
print("-" * 60)

results = []
for window_start, window_stop, contig_name, expected in test_cases:
    actual = calculate_callable_length_corrected(window_start, window_stop, contig_name, bed_df)
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
    print("✓ All tests passed! Corrected implementation works correctly.")
    print("Now we can implement the PyRanges version and verify it matches!")
else:
    print("✗ Some tests failed. Need to debug further.")

# Save results for comparison
results_df = pd.DataFrame(results)
results_df.to_csv('corrected_implementation_results_fixed.csv', index=False)
print("Results saved to: corrected_implementation_results_fixed.csv")
