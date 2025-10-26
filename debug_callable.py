#!/usr/bin/env python3
"""
Debug the callable length calculation step by step.
"""

import pandas as pd
import numpy as np

# Load BED file
bed_path = "/var/folders/v3/ybcp76dx0bx9gy8dznldx29c0000gp/T/tmp4p5mtw17.bed"
bed_df = pd.read_csv(bed_path, sep='\t', header=None, names=['chrom', 'start', 'end'])
print("BED file:")
print(bed_df)

def debug_callable_length(window_start, window_stop, contig_name, bed_df):
    """Debug version that shows all steps."""
    print(f"\n=== Debugging window [{window_start}, {window_stop}) on {contig_name} ===")
    window_length = window_stop - window_start
    print(f"Window length: {window_length}bp")
    
    masked_length = 0
    
    for i, row in bed_df.iterrows():
        chrom = row['chrom']
        bed_start = row['start']
        bed_end = row['end']
        
        print(f"  BED region {i}: {chrom} [{bed_start}, {bed_end})")
        
        # Only process BED regions on the same contig
        if str(chrom) == str(contig_name):
            print(f"    → Same contig! Checking overlap...")
            
            # Calculate overlap between BED region and window
            overlap_start = max(window_start, bed_start)
            overlap_end = min(window_stop, bed_end)
            
            print(f"    → Overlap: [{overlap_start}, {overlap_end})")
            
            if overlap_start < overlap_end:
                overlap_length = overlap_end - overlap_start
                masked_length += overlap_length
                print(f"    → Overlap length: {overlap_length}bp (added to masked)")
            else:
                print(f"    → No overlap")
        else:
            print(f"    → Different contig, skipping")
    
    callable_length = window_length - masked_length
    print(f"Total masked: {masked_length}bp")
    print(f"Callable length: {window_length} - {masked_length} = {callable_length}bp")
    
    return max(1, callable_length)

# Test the failing cases
print("Testing failing cases:")

# Case 1: chr1 [150, 250) - should be 100bp
debug_callable_length(150, 250, 'chr1', bed_df)

# Case 2: chr1 [250, 350) - should be 100bp  
debug_callable_length(250, 350, 'chr1', bed_df)

# Case 3: chr1 [450, 550) - should be 100bp
debug_callable_length(450, 550, 'chr1', bed_df)

# Case 4: chr2 [100, 200) - should be 100bp
debug_callable_length(100, 200, 'chr2', bed_df)
