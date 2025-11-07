#!/bin/bash
# Script to create a non-callable sites BED file from an all-sites VCF
# This identifies sites that are missing or have low quality in the all-sites VCF

# Usage: ./create_noncallable_bed.sh <all_sites.vcf.gz> <output.bed>

set -e

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <all_sites_vcf.gz> <output_noncallable.bed>"
    echo ""
    echo "This script creates a BED file of non-callable sites by:"
    echo "  1. Extracting sites with missing genotypes (./.) or low quality"
    echo "  2. Converting to 0-based BED format"
    echo ""
    echo "Example:"
    echo "  $0 /path/to/all_sites.vcf.gz non_callable_sites.bed"
    exit 1
fi

ALL_SITES_VCF="$1"
OUTPUT_BED="$2"

echo "=========================================="
echo "Creating Non-Callable Sites BED File"
echo "=========================================="
echo ""
echo "Input VCF: $ALL_SITES_VCF"
echo "Output BED: $OUTPUT_BED"
echo ""

# Check if input file exists
if [ ! -f "$ALL_SITES_VCF" ]; then
    echo "ERROR: Input VCF file not found: $ALL_SITES_VCF"
    exit 1
fi

# Check if bcftools is available
if ! command -v bcftools &> /dev/null; then
    echo "ERROR: bcftools not found. Please load bcftools module or install it."
    exit 1
fi

echo "Step 1: Extracting non-callable sites from VCF..."
echo "  (Sites with missing genotypes, FILTER != PASS, or low quality)"
echo ""

# Extract sites where:
# - FILTER is not PASS (failed quality filters)
# - OR genotypes contain missing data (./. or .|.)
# Convert to BED format (0-based, half-open)

bcftools view -H "$ALL_SITES_VCF" | \
awk -v OFS="\t" '
BEGIN {
    print "# Non-callable sites BED file (0-based, half-open intervals)" > "/dev/stderr"
    print "# Created from: " ARGV[1] > "/dev/stderr"
    print "# Criteria: FILTER != PASS OR missing genotypes" > "/dev/stderr"
}
{
    chrom = $1
    pos = $2  # 1-based VCF position
    filter = $7
    genotypes = $10  # Assuming single sample, adjust if multi-sample
    
    # Check if site is non-callable
    is_noncallable = 0
    
    # Check FILTER field
    if (filter != "PASS" && filter != ".") {
        is_noncallable = 1
    }
    
    # Check for missing genotypes in any sample
    for (i = 10; i <= NF; i++) {
        if ($i ~ /^\.\/\./ || $i ~ /^\.\|\./ || $i ~ /^\.[:\/\|]/ || $i ~ /[:\/\|]\.$/) {
            is_noncallable = 1
            break
        }
    }
    
    # Output as BED (0-based, half-open)
    if (is_noncallable) {
        bed_start = pos - 1  # Convert 1-based to 0-based
        bed_end = pos        # Half-open interval
        print chrom, bed_start, bed_end
    }
}' | \
sort -k1,1 -k2,2n | \
bedtools merge -i - > "$OUTPUT_BED"

echo ""
echo "Step 2: Merging adjacent non-callable sites..."
echo ""

# Count sites
TOTAL_SITES=$(bcftools view -H "$ALL_SITES_VCF" | wc -l)
NONCALLABLE_REGIONS=$(wc -l < "$OUTPUT_BED")
NONCALLABLE_BP=$(awk '{sum += $3-$2} END {print sum+0}' "$OUTPUT_BED")

echo "=========================================="
echo "Summary:"
echo "=========================================="
echo "Total sites in VCF: $TOTAL_SITES"
echo "Non-callable regions: $NONCALLABLE_REGIONS"
echo "Non-callable bases: $NONCALLABLE_BP"
echo ""
echo "Output written to: $OUTPUT_BED"
echo ""
echo "First 10 regions:"
head -10 "$OUTPUT_BED"
echo ""
echo "=========================================="
echo "Done!"
echo "=========================================="

