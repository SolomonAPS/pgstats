#!/bin/bash

# Concatenate results from array job
# Usage: ./concat_results.sh <results_dir> <output_file>

RESULTS_DIR="${1:-/proj/johrilab/projects/DFEpos/dpgp3/stats/array_results}"
OUTPUT_FILE="${2:-/proj/johrilab/projects/DFEpos/dpgp3/stats/dpgp3_100_regions_102925.tsv}"

# Check if results directory exists
if [ ! -d "$RESULTS_DIR" ]; then
    echo "Error: Results directory not found: $RESULTS_DIR"
    exit 1
fi

# Count CSV files
N_FILES=$(ls -1 ${RESULTS_DIR}/region_*.csv 2>/dev/null | wc -l)
if [ "$N_FILES" -eq 0 ]; then
    echo "Error: No region_*.csv files found in $RESULTS_DIR"
    exit 1
fi

echo "Found $N_FILES result files"

# Get header from first file (sorted numerically)
FIRST_FILE=$(ls -1v ${RESULTS_DIR}/region_*.csv | head -n 1)
head -n 1 "$FIRST_FILE" > "$OUTPUT_FILE"

# Append all data (skip headers)
for file in $(ls -1v ${RESULTS_DIR}/region_*.csv); do
    tail -n +2 "$file" >> "$OUTPUT_FILE"
done

# Count lines
N_LINES=$(wc -l < "$OUTPUT_FILE")
echo "Combined $N_FILES files into $OUTPUT_FILE"
echo "Total lines (including header): $N_LINES"
echo "Total windows: $((N_LINES - 1))"

