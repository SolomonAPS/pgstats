#!/usr/bin/env python3
"""
Debug PyRanges overlap detection.
"""

import pandas as pd
import numpy as np
import pyranges as pr

# Get BED file path from previous run
bed_path = "/var/folders/v3/ybcp76dx0bx9gy8dznldx29c0000gp/T/tmp4p5mtw17.bed"

# Load BED file and convert to PyRanges
bed_df = pd.read_csv(bed_path, sep='\t', header=None, names=['chrom', 'start', 'end'])
bed_gr = pr.PyRanges(bed_df.rename(columns={'chrom': 'Chromosome', 'start': 'Start', 'end': 'End'}))

print("BED PyRanges object:")
print(bed_gr)

def debug_pyranges_overlap(window_start, window_stop, contig_name, bed_gr):
    """Debug PyRanges overlap detection step by step."""
    print(f"\n=== Debugging window [{window_start}, {window_stop}) on {contig_name} ===")
    window_length = window_stop - window_start
    print(f"Window length: {window_length}bp")
    
    # Create window PyRanges object
    window_gr = pr.PyRanges(pd.DataFrame({
        'Chromosome': [contig_name],
        'Start': [window_start],
        'End': [window_stop]
    }))
    
    print(f"Window PyRanges object:")
    print(window_gr)
    
    # Find overlaps with BED regions
    print(f"\nFinding overlaps...")
    overlaps = window_gr.overlap(bed_gr)
    print(f"Overlaps found: {len(overlaps)}")
    
    if len(overlaps) > 0:
        print(f"Overlap details:")
        print(overlaps.df)
        
        # Calculate total masked length
        masked_length = 0
        for _, row in overlaps.df.iterrows():
            print(f"  Processing overlap: BED [{row['Start']}, {row['End']})")
            # Calculate overlap between window and BED region
            overlap_start = max(window_start, row['Start'])
            overlap_end = min(window_stop, row['End'])
            print(f"    → Overlap: [{overlap_start}, {overlap_end})")
            if overlap_start < overlap_end:
                overlap_length = overlap_end - overlap_start
                masked_length += overlap_length
                print(f"    → Overlap length: {overlap_length}bp (added to masked)")
            else:
                print(f"    → No overlap")
        
        callable_length = window_length - masked_length
        print(f"Total masked: {masked_length}bp")
        print(f"Callable length: {window_length} - {masked_length} = {callable_length}bp")
    else:
        print("No overlaps found - entire window is callable")
        callable_length = window_length
    
    final_result = max(1, callable_length)
    print(f"Final result: max(1, {callable_length}) = {final_result}bp")
    
    return final_result

# Test a few cases
print("Testing PyRanges overlap detection:")

# Case 1: chr1 [150, 250) - should find overlap with BED [100, 200)
debug_pyranges_overlap(150, 250, 'chr1', bed_gr)

# Case 2: chr1 [650, 750) - should find no overlaps
debug_pyranges_overlap(650, 750, 'chr1', bed_gr)

# Case 3: chr2 [100, 200) - should find overlap with BED [50, 150)
debug_pyranges_overlap(100, 200, 'chr2', bed_gr)
