#!/bin/bash
#SBATCH --job-name=analyze_regions
#SBATCH --output=analyze_regions_%j.out
#SBATCH --error=analyze_regions_%j.err
#SBATCH --time=2:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --partition=general

# Analyze regions for callable sites and AA_prob quality
# For each region, output:
# - Region coordinates
# - Total bp
# - Callable bp (from VCF)
# - Non-callable bp (from BED mask)
# - AA_prob distribution (missing, <0.5, 0.5-0.8, 0.8-0.95, >=0.95)

set -e

# Load modules
module load bcftools
module load bedtools

# Configuration
REGIONS_BED="/proj/johrilab/projects/DFEpos/dpgp3/stats/regions.bed"
VCF_DIR="/proj/johrilab/projects/DFEpos/dpgp3/test_adaptive_priors/processed_11"
NONCALLABLE_BED="/proj/johrilab/projects/DFEpos/dpgp3/test_adaptive_priors/non_callable_sites_11.bed"
OUTPUT_TSV="/proj/johrilab/projects/DFEpos/dpgp3/stats/regions_quality.tsv"

echo "=========================================="
echo "Analyzing Region Quality"
echo "=========================================="
echo ""
echo "Regions file: $REGIONS_BED"
echo "VCF directory: $VCF_DIR"
echo "Non-callable BED: $NONCALLABLE_BED"
echo "Output: $OUTPUT_TSV"
echo ""

# Check files exist
if [ ! -f "$REGIONS_BED" ]; then
    echo "ERROR: Regions file not found: $REGIONS_BED"
    exit 1
fi

if [ ! -f "$NONCALLABLE_BED" ]; then
    echo "ERROR: Non-callable BED not found: $NONCALLABLE_BED"
    exit 1
fi

# Create output header
echo -e "region_id\tchrom\tstart\tend\ttotal_bp\tcallable_bp\tnoncallable_bp\tcallable_pct\tsites_in_vcf\tAA_missing\tAA_missing_pct\tAA_low\tAA_low_pct\tAA_medium\tAA_medium_pct\tAA_high\tAA_high_pct\tAA_very_high\tAA_very_high_pct" > "$OUTPUT_TSV"

# Summary variables
TOTAL_REGIONS=0
TOTAL_BP=0
TOTAL_CALLABLE_BP=0
TOTAL_NONCALLABLE_BP=0
TOTAL_SITES=0
TOTAL_AA_MISSING=0
TOTAL_AA_LOW=0
TOTAL_AA_MEDIUM=0
TOTAL_AA_HIGH=0
TOTAL_AA_VERY_HIGH=0

echo "Processing regions..."
echo ""

# Process each region
REGION_ID=0
while IFS=$'\t' read -r CHROM START END; do
    REGION_ID=$((REGION_ID + 1))
    TOTAL_REGIONS=$((TOTAL_REGIONS + 1))
    
    # Calculate total bp
    TOTAL_REGION_BP=$((END - START))
    TOTAL_BP=$((TOTAL_BP + TOTAL_REGION_BP))
    
    # Get chromosome arm for VCF file (chr3L -> 3L)
    CHR_ARM=${CHROM#chr}
    VCF_FILE="${VCF_DIR}/dpgp3_${CHR_ARM}_polarized_clean.vcf.gz"
    
    if [ ! -f "$VCF_FILE" ]; then
        echo "WARNING: VCF not found for $CHROM: $VCF_FILE"
        continue
    fi
    
    # Extract sites in this region from VCF
    # Count sites and get AA_prob distribution
    AA_STATS=$(bcftools view -H -r "${CHROM}:${START}-${END}" "$VCF_FILE" | \
    bcftools query -f '%AA_prob\n' | \
    awk '
    {
        total++
        aa_prob = $1
        if (aa_prob == ".") missing++
        else if (aa_prob < 0.5) low++
        else if (aa_prob < 0.8) medium++
        else if (aa_prob < 0.95) high++
        else very_high++
    }
    END {
        printf "%d\t%d\t%d\t%d\t%d\t%d", total, missing, low, medium, high, very_high
    }
    ')
    
    read SITES_IN_VCF AA_MISSING AA_LOW AA_MEDIUM AA_HIGH AA_VERY_HIGH <<< "$AA_STATS"
    
    # Calculate callable bp (sites present in VCF)
    # Each site in VCF is 1bp callable
    CALLABLE_BP=$SITES_IN_VCF
    
    # Calculate non-callable bp in this region
    NONCALLABLE_BP=$(bedtools intersect -a <(echo -e "${CHROM}\t${START}\t${END}") -b "$NONCALLABLE_BED" | \
    awk '{sum += $3-$2} END {print sum+0}')
    
    # Calculate percentages
    CALLABLE_PCT=$(awk "BEGIN {printf \"%.2f\", ($CALLABLE_BP/$TOTAL_REGION_BP)*100}")
    
    if [ $SITES_IN_VCF -gt 0 ]; then
        AA_MISSING_PCT=$(awk "BEGIN {printf \"%.2f\", ($AA_MISSING/$SITES_IN_VCF)*100}")
        AA_LOW_PCT=$(awk "BEGIN {printf \"%.2f\", ($AA_LOW/$SITES_IN_VCF)*100}")
        AA_MEDIUM_PCT=$(awk "BEGIN {printf \"%.2f\", ($AA_MEDIUM/$SITES_IN_VCF)*100}")
        AA_HIGH_PCT=$(awk "BEGIN {printf \"%.2f\", ($AA_HIGH/$SITES_IN_VCF)*100}")
        AA_VERY_HIGH_PCT=$(awk "BEGIN {printf \"%.2f\", ($AA_VERY_HIGH/$SITES_IN_VCF)*100}")
    else
        AA_MISSING_PCT=0
        AA_LOW_PCT=0
        AA_MEDIUM_PCT=0
        AA_HIGH_PCT=0
        AA_VERY_HIGH_PCT=0
    fi
    
    # Write to output file
    echo -e "${REGION_ID}\t${CHROM}\t${START}\t${END}\t${TOTAL_REGION_BP}\t${CALLABLE_BP}\t${NONCALLABLE_BP}\t${CALLABLE_PCT}\t${SITES_IN_VCF}\t${AA_MISSING}\t${AA_MISSING_PCT}\t${AA_LOW}\t${AA_LOW_PCT}\t${AA_MEDIUM}\t${AA_MEDIUM_PCT}\t${AA_HIGH}\t${AA_HIGH_PCT}\t${AA_VERY_HIGH}\t${AA_VERY_HIGH_PCT}" >> "$OUTPUT_TSV"
    
    # Update totals
    TOTAL_CALLABLE_BP=$((TOTAL_CALLABLE_BP + CALLABLE_BP))
    TOTAL_NONCALLABLE_BP=$((TOTAL_NONCALLABLE_BP + NONCALLABLE_BP))
    TOTAL_SITES=$((TOTAL_SITES + SITES_IN_VCF))
    TOTAL_AA_MISSING=$((TOTAL_AA_MISSING + AA_MISSING))
    TOTAL_AA_LOW=$((TOTAL_AA_LOW + AA_LOW))
    TOTAL_AA_MEDIUM=$((TOTAL_AA_MEDIUM + AA_MEDIUM))
    TOTAL_AA_HIGH=$((TOTAL_AA_HIGH + AA_HIGH))
    TOTAL_AA_VERY_HIGH=$((TOTAL_AA_VERY_HIGH + AA_VERY_HIGH))
    
    # Print progress every 10 regions
    if [ $((REGION_ID % 10)) -eq 0 ]; then
        echo "  Processed $REGION_ID regions..."
    fi
    
done < "$REGIONS_BED"

echo ""
echo "=========================================="
echo "Summary Statistics"
echo "=========================================="
echo ""
echo "Total regions analyzed: $TOTAL_REGIONS"
echo "Total bp in regions: $TOTAL_BP"
echo "Total callable bp: $TOTAL_CALLABLE_BP ($(awk "BEGIN {printf \"%.2f\", ($TOTAL_CALLABLE_BP/$TOTAL_BP)*100}")%)"
echo "Total non-callable bp: $TOTAL_NONCALLABLE_BP ($(awk "BEGIN {printf \"%.2f\", ($TOTAL_NONCALLABLE_BP/$TOTAL_BP)*100}")%)"
echo ""
echo "Total sites in VCFs: $TOTAL_SITES"
echo ""
echo "AA_prob Distribution:"
echo "  Missing: $TOTAL_AA_MISSING ($(awk "BEGIN {printf \"%.2f\", ($TOTAL_AA_MISSING/$TOTAL_SITES)*100}")%)"
echo "  < 0.5 (low): $TOTAL_AA_LOW ($(awk "BEGIN {printf \"%.2f\", ($TOTAL_AA_LOW/$TOTAL_SITES)*100}")%)"
echo "  0.5-0.8 (medium): $TOTAL_AA_MEDIUM ($(awk "BEGIN {printf \"%.2f\", ($TOTAL_AA_MEDIUM/$TOTAL_SITES)*100}")%)"
echo "  0.8-0.95 (high): $TOTAL_AA_HIGH ($(awk "BEGIN {printf \"%.2f\", ($TOTAL_AA_HIGH/$TOTAL_SITES)*100}")%)"
echo "  >= 0.95 (very high): $TOTAL_AA_VERY_HIGH ($(awk "BEGIN {printf \"%.2f\", ($TOTAL_AA_VERY_HIGH/$TOTAL_SITES)*100}")%)"
echo ""
echo "Output written to: $OUTPUT_TSV"
echo ""
echo "First 10 regions:"
head -11 "$OUTPUT_TSV" | column -t
echo ""
echo "=========================================="
echo "Done!"
echo "=========================================="
