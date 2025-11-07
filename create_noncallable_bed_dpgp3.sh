#!/bin/bash
# Script to create a non-callable sites BED file from DPGP3 all-sites VCFs
# Non-callable sites are positions NOT present in the all-sites VCF (gaps between variants)
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
echo "Strategy: Non-callable = gaps between positions in all-sites VCF"
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

# Get chromosome lengths from FASTA files
declare -A CHR_LENGTHS

echo "Reading chromosome lengths from FASTA files..."
for CHR in 2L 2R 3L 3R; do
    FASTA_FILE="${VCF_DIR}/dpgp3_${CHR}.fasta"
    CHROM_NAME="chr${CHR}"
    
    if [ -f "$FASTA_FILE" ]; then
        # Get length from FASTA (count all non-header lines)
        LENGTH=$(grep -v "^>" "$FASTA_FILE" | tr -d '\n' | wc -c | tr -d ' ')
        CHR_LENGTHS[$CHROM_NAME]=$LENGTH
        echo "  ${CHROM_NAME}: $LENGTH bp"
    else
        echo "  WARNING: FASTA not found: $FASTA_FILE"
    fi
done
echo ""

# Process each chromosome arm
for CHR in 2L 2R 3L 3R; do
    VCF_FILE="${VCF_DIR}/dpgp3_${CHR}_mono_biallelic.vcf.gz"
    TEMP_BED="${TEMP_DIR}/${CHR}_noncallable.bed"
    CHROM_NAME="chr${CHR}"
    CHROM_LENGTH=${CHR_LENGTHS[$CHROM_NAME]}
    
    if [ -z "$CHROM_LENGTH" ]; then
        echo "ERROR: Chromosome length not found for ${CHROM_NAME}"
        echo "Skipping..."
        continue
    fi
    
    echo "=========================================="
    echo "Processing chromosome ${CHR}"
    echo "=========================================="
    echo "Input: $VCF_FILE"
    echo "Chromosome length: $CHROM_LENGTH bp"
    
    if [ ! -f "$VCF_FILE" ]; then
        echo "WARNING: VCF file not found: $VCF_FILE"
        echo "Skipping..."
        continue
    fi
    
    echo "Extracting callable positions..."
    
    # Extract all positions from VCF (these are callable)
    # Convert to 0-based BED format for each position
    bcftools view -H "$VCF_FILE" | \
    awk -v OFS="\t" '{
        pos = $2  # 1-based VCF position
        bed_start = pos - 1  # 0-based
        bed_end = pos        # half-open
        print $1, bed_start, bed_end
    }' | sort -k1,1 -k2,2n > "${TEMP_DIR}/${CHR}_callable.bed"
    
    CALLABLE_SITES=$(wc -l < "${TEMP_DIR}/${CHR}_callable.bed")
    echo "  Callable sites in VCF: $CALLABLE_SITES"
    
    # Create a BED file for the entire chromosome
    echo -e "${CHROM_NAME}\t0\t${CHROM_LENGTH}" > "${TEMP_DIR}/${CHR}_whole.bed"
    
    # Subtract callable sites from whole chromosome to get non-callable regions
    echo "  Finding gaps (non-callable regions)..."
    bedtools subtract -a "${TEMP_DIR}/${CHR}_whole.bed" -b "${TEMP_DIR}/${CHR}_callable.bed" > "$TEMP_BED"
    
    # Count non-callable regions
    NONCALLABLE_REGIONS=$(wc -l < "$TEMP_BED")
    NONCALLABLE_BP=$(awk '{sum += $3-$2} END {print sum+0}' "$TEMP_BED")
    CALLABLE_BP=$((CHROM_LENGTH - NONCALLABLE_BP))
    
    echo "  Callable bases: $CALLABLE_BP"
    echo "  Non-callable regions: $NONCALLABLE_REGIONS"
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
