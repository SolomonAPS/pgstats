#!/usr/bin/env python3
"""
Test callable length calculation with known BED regions.
"""

import pandas as pd
import numpy as np
from pathlib import Path
import tempfile
import os

# Create test BED file
test_bed_content = """chr1	100	200
chr1	300	400
chr1	500	600
chr2	50	150
chr2	250	350
"""

# Create temporary BED file
with tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False) as f:
    f.write(test_bed_content)
    test_bed_path = f.name

print("Test BED file created:")
print(test_bed_content)
print(f"Saved to: {test_bed_path}")

# Test cases: (window_start, window_stop, expected_callable_length)
test_cases = [
    # chr1 tests
    (50, 150, 50),    # Window [50,150): BED [100,200) → overlap [100,150) = 50bp masked → 100-50=50 callable
    (150, 250, 100),  # Window [150,250): No BED overlap → 100bp callable
    (250, 350, 100),  # Window [250,350): No BED overlap → 100bp callable
    (350, 450, 50),   # Window [350,450): BED [300,400) → overlap [350,400) = 50bp masked → 100-50=50 callable
    (450, 550, 100),  # Window [450,550): No BED overlap → 100bp callable
    (550, 650, 50),   # Window [550,650): BED [500,600) → overlap [550,600) = 50bp masked → 100-50=50 callable
    (650, 750, 100),  # Window [650,750): No BED overlap → 100bp callable
    
    # chr2 tests  
    (0, 100, 50),     # Window [0,100): BED [50,150) → overlap [50,100) = 50bp masked → 100-50=50 callable
    (100, 200, 100),  # Window [100,200): No BED overlap → 100bp callable
    (200, 300, 50),   # Window [200,300): BED [250,350) → overlap [250,300) = 50bp masked → 100-50=50 callable
    (300, 400, 50),   # Window [300,400): BED [250,350) → overlap [300,350) = 50bp masked → 100-50=50 callable
    (400, 500, 100),  # Window [400,500): No BED overlap → 100bp callable
    
    # Edge cases
    (95, 105, 5),     # Window [95,105): BED [100,200) → overlap [100,105) = 5bp masked → 10-5=5 callable
    (195, 205, 5),    # Window [195,205): BED [100,200) → overlap [195,200) = 5bp masked → 10-5=5 callable
]

print(f"\nTest cases ({len(test_cases)} total):")
print("Window\t\tExpected Callable")
print("-" * 30)
for start, stop, expected in test_cases:
    print(f"[{start:3d}, {stop:3d})\t{expected:3d}bp")

# Save test cases for later use
test_cases_df = pd.DataFrame(test_cases, columns=['window_start', 'window_stop', 'expected_callable'])
test_cases_df.to_csv('test_callable_cases.csv', index=False)
print(f"\nTest cases saved to: test_callable_cases.csv")

print(f"\nBED file path: {test_bed_path}")
print("Use this path in your GenomicDataset test!")
