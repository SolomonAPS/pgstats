#!/bin/bash
# Script to create a non-callable sites BED file from DPGP3 all-sites VCFs
# Processes 2L, 2R, 3L, 3R separately and concatenates them

set -e

# Configuration
VCF_DIR="/proj/johrilab/projects/DFEpos/dpgp3/vcf"
OUTPUT_BED="/proj/johrilab/projects/DFEpos/dpgp3/masking/non_callable_sites_corrected.bed"
TEMP_DIR=$(mktemp -d)

echo "=========================================="
echo "Creating Non-Callable Sites BED File"
echo "=========================================="
echo ""
echo "VCF directory: $VCF_DIR"
echo "Output BED: $OUTPUT_BED"
echo "Temp directory: $TEMP_DIR"
echo ""

# Check if bcftools is available
if ! command -v bcftools &> /dev/null; then
    echo "ERROR: bcftools not found. Please load bcftools module or install it."
    exit 1
fi

# Check if bedtools is available
if ! command -v bedtools &> /dev/null; then
    echo "ERROR: bedtools not found. Please load bedtools module or install it."
    exit 1
fi

# Process each chromosome arm
for CHR in 2L 2R 3L 3R; do
    VCF_FILE="${VCF_DIR}/dpgp3_${CHR}_mono_biallelic.vcf.gz"
    TEMP_BED="${TEMP_DIR}/${CHR}_noncallable.bed"
    
    echo "=========================================="
    echo "Processing chromosome ${CHR}"
    echo "=========================================="
    echo "Input: $VCF_FILE"
    
    if [ ! -f "$VCF_FILE" ]; then
        echo "WARNING: VCF file not found: $VCF_FILE"
        echo "Skipping..."
        continue
    fi
    
    echo "Extracting non-callable sites..."
    
    # Extract sites where any sample has missing genotype (.)
    # In this VCF format, missing genotypes are just "." not "./."
    bcftools view -H "$VCF_FILE" | \
    awk -v OFS="\t" '
    {
        chrom = $1
        pos = $2  # 1-based VCF position
        
        # Check if any sample (columns 10+) has missing genotype
        is_noncallable = 0
        for (i = 10; i <= NF; i++) {
            # Split on : to get GT field
            split($i, gt_fields, ":")
            gt = gt_fields[1]
            
            # Check if genotype is missing (just "." or contains ".")
            if (gt == "." || gt ~ /^\./) {
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
    }' > "$TEMP_BED"
    
    # Count sites
    NONCALLABLE_SITES=$(wc -l < "$TEMP_BED")
    echo "  Non-callable sites found: $NONCALLABLE_SITES"
    
    # Merge adjacent sites
    echo "  Merging adjacent sites..."
    bedtools merge -i "$TEMP_BED" > "${TEMP_BED}.merged"
    mv "${TEMP_BED}.merged" "$TEMP_BED"
    
    NONCALLABLE_REGIONS=$(wc -l < "$TEMP_BED")
    NONCALLABLE_BP=$(awk '{sum += $3-$2} END {print sum+0}' "$TEMP_BED")
    
    echo "  Non-callable regions (merged): $NONCALLABLE_REGIONS"
    echo "  Non-callable bases: $NONCALLABLE_BP"
    echo ""
done

echo "=========================================="
echo "Concatenating all chromosomes..."
echo "=========================================="

# Concatenate all chromosome BED files and sort
cat "${TEMP_DIR}"/*_noncallable.bed | \
    sort -k1,1 -k2,2n > "$OUTPUT_BED"

# Summary statistics
TOTAL_REGIONS=$(wc -l < "$OUTPUT_BED")
TOTAL_BP=$(awk '{sum += $3-$2} END {print sum+0}' "$OUTPUT_BED")

echo ""
echo "=========================================="
echo "Summary:"
echo "=========================================="
echo "Total non-callable regions: $TOTAL_REGIONS"
echo "Total non-callable bases: $TOTAL_BP"
echo ""
echo "Output written to: $OUTPUT_BED"
echo ""
echo "First 10 regions:"
head -10 "$OUTPUT_BED"
echo ""
echo "Last 10 regions:"
tail -10 "$OUTPUT_BED"
echo ""

# Cleanup
rm -rf "$TEMP_DIR"
echo "Cleaned up temporary files"
echo ""
echo "=========================================="
echo "Done!"
echo "=========================================="

