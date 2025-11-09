#!/bin/bash
# Script to create a non-callable sites BED file from DPGP3 all-sites VCFs
# Non-callable sites are positions NOT present in the all-sites VCF (gaps between variants)
# Processes 2L, 2R, 3L, 3R separately and concatenates them

set -e

# Configuration
VCF_DIR="/proj/johrilab/projects/DFEpos/dpgp3/vcf"
REF_FAI="/proj/johrilab/projects/DFEpos/dpgp3/reference/dmel_r5.fasta.fai"
OUTPUT_BED="/proj/johrilab/projects/DFEpos/dpgp3/masking/non_callable_sites_corrected.bed"
TEMP_DIR=$(mktemp -d)

echo "=========================================="
echo "Creating Non-Callable Sites BED File"
echo "=========================================="
echo ""
echo "Strategy: Non-callable = gaps between positions in all-sites VCF"
echo "VCF directory: $VCF_DIR"
echo "Reference FAI: $REF_FAI"
echo "Output BED: $OUTPUT_BED"
echo "Temp directory: $TEMP_DIR"
echo ""

# Check if FAI file exists
if [ ! -f "$REF_FAI" ]; then
    echo "ERROR: Reference FAI file not found: $REF_FAI"
    exit 1
fi

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

# Get chromosome lengths from FAI file
declare -A CHR_LENGTHS

echo "Reading chromosome lengths from FAI file..."
while read -r chrom length rest; do
    # Map chromosome names (FAI uses "2L", VCF uses "chr2L")
    case "$chrom" in
        2L) CHR_LENGTHS["chr2L"]=$length; echo "  chr2L: $length bp" ;;
        2R) CHR_LENGTHS["chr2R"]=$length; echo "  chr2R: $length bp" ;;
        3L) CHR_LENGTHS["chr3L"]=$length; echo "  chr3L: $length bp" ;;
        3R) CHR_LENGTHS["chr3R"]=$length; echo "  chr3R: $length bp" ;;
    esac
done < "$REF_FAI"
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
    
    # Extract all positions from VCF (these are callable)
    # Convert to 0-based BED format for each position, then merge consecutive positions
    echo "  Extracting and merging callable positions..."
    bcftools view -H "$VCF_FILE" | \
    awk -v OFS="\t" '{
        pos = $2  # 1-based VCF position
        bed_start = pos - 1  # 0-based
        bed_end = pos        # half-open
        print $1, bed_start, bed_end
    }' | sort -k1,1 -k2,2n | \
    bedtools merge -i - > "${TEMP_DIR}/${CHR}_callable.bed"
    
    CALLABLE_REGIONS=$(wc -l < "${TEMP_DIR}/${CHR}_callable.bed")
    CALLABLE_BP=$(awk '{sum += $3-$2} END {print sum+0}' "${TEMP_DIR}/${CHR}_callable.bed")
    echo "  Callable regions (after merging): $CALLABLE_REGIONS"
    echo "  Callable bases: $CALLABLE_BP"
    
    # Create a BED file for the entire chromosome
    echo -e "${CHROM_NAME}\t0\t${CHROM_LENGTH}" > "${TEMP_DIR}/${CHR}_whole.bed"
    
    # Subtract callable sites from whole chromosome to get non-callable regions
    echo "  Finding gaps (non-callable regions)..."
    bedtools subtract -a "${TEMP_DIR}/${CHR}_whole.bed" -b "${TEMP_DIR}/${CHR}_callable.bed" > "$TEMP_BED"
    
    # Count non-callable regions
    NONCALLABLE_REGIONS=$(wc -l < "$TEMP_BED")
    NONCALLABLE_BP=$(awk '{sum += $3-$2} END {print sum+0}' "$TEMP_BED")
    
    echo "  Non-callable regions: $NONCALLABLE_REGIONS"
    echo "  Non-callable bases: $NONCALLABLE_BP ($(awk "BEGIN {printf \"%.2f\", ($NONCALLABLE_BP/$CHROM_LENGTH)*100}")%)"
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
