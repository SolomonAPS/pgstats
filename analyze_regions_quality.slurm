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
echo -e "region_id\tchrom\tstart\tend\ttotal_bp\tcallable_bp\tnoncallable_bp\tcallable_pct\tall_sites\tbiallelic_sites\tAA_missing\tAA_missing_pct\tAA_low\tAA_low_pct\tAA_medium\tAA_medium_pct\tAA_high\tAA_high_pct\tAA_very_high\tAA_very_high_pct" > "$OUTPUT_TSV"

# Summary variables
TOTAL_REGIONS=0
TOTAL_BP=0
TOTAL_CALLABLE_BP=0
TOTAL_NONCALLABLE_BP=0
TOTAL_ALL_SITES=0
TOTAL_BIALLELIC_SITES=0
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
    
    # Get chromosome arm for VCF files (chr3L -> 3L)
    CHR_ARM=${CHROM#chr}
    ALL_SITES_VCF="${VCF_DIR}/dpgp3_${CHR_ARM}_polarized_clean.vcf.gz"
    BIALLELIC_VCF="${VCF_DIR}/dpgp3_${CHR_ARM}_polarized_biallelic.vcf.gz"
    
    if [ ! -f "$ALL_SITES_VCF" ]; then
        echo "WARNING: All-sites VCF not found for $CHROM: $ALL_SITES_VCF"
        continue
    fi
    
    if [ ! -f "$BIALLELIC_VCF" ]; then
        echo "WARNING: Biallelic VCF not found for $CHROM: $BIALLELIC_VCF"
        continue
    fi
    
    # Count total sites in region (from all-sites VCF for callable length)
    SITES_IN_VCF=$(bcftools view -H -r "${CHROM}:${START}-${END}" "$ALL_SITES_VCF" | wc -l)
    
    # Get AA_prob distribution from BIALLELIC sites only (variant sites used for stats)
    AA_STATS=$(bcftools query -r "${CHROM}:${START}-${END}" -f '%AA_prob\n' "$BIALLELIC_VCF" | \
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
    
    read BIALLELIC_SITES AA_MISSING AA_LOW AA_MEDIUM AA_HIGH AA_VERY_HIGH <<< "$AA_STATS"
    
    # Calculate callable bp (sites present in all-sites VCF)
    # Each site in VCF is 1bp callable
    CALLABLE_BP=$SITES_IN_VCF
    
    # Calculate non-callable bp in this region
    # Use temp file instead of process substitution for bedtools compatibility
    TEMP_REGION=$(mktemp)
    echo -e "${CHROM}\t${START}\t${END}" > "$TEMP_REGION"
    NONCALLABLE_BP=$(bedtools intersect -a "$TEMP_REGION" -b "$NONCALLABLE_BED" | \
    awk '{sum += $3-$2} END {print sum+0}')
    rm -f "$TEMP_REGION"
    
    # Calculate percentages
    CALLABLE_PCT=$(awk "BEGIN {printf \"%.2f\", ($CALLABLE_BP/$TOTAL_REGION_BP)*100}")
    
    if [ $BIALLELIC_SITES -gt 0 ]; then
        AA_MISSING_PCT=$(awk "BEGIN {printf \"%.2f\", ($AA_MISSING/$BIALLELIC_SITES)*100}")
        AA_LOW_PCT=$(awk "BEGIN {printf \"%.2f\", ($AA_LOW/$BIALLELIC_SITES)*100}")
        AA_MEDIUM_PCT=$(awk "BEGIN {printf \"%.2f\", ($AA_MEDIUM/$BIALLELIC_SITES)*100}")
        AA_HIGH_PCT=$(awk "BEGIN {printf \"%.2f\", ($AA_HIGH/$BIALLELIC_SITES)*100}")
        AA_VERY_HIGH_PCT=$(awk "BEGIN {printf \"%.2f\", ($AA_VERY_HIGH/$BIALLELIC_SITES)*100}")
    else
        AA_MISSING_PCT=0
        AA_LOW_PCT=0
        AA_MEDIUM_PCT=0
        AA_HIGH_PCT=0
        AA_VERY_HIGH_PCT=0
    fi
    
    # Write to output file
    echo -e "${REGION_ID}\t${CHROM}\t${START}\t${END}\t${TOTAL_REGION_BP}\t${CALLABLE_BP}\t${NONCALLABLE_BP}\t${CALLABLE_PCT}\t${SITES_IN_VCF}\t${BIALLELIC_SITES}\t${AA_MISSING}\t${AA_MISSING_PCT}\t${AA_LOW}\t${AA_LOW_PCT}\t${AA_MEDIUM}\t${AA_MEDIUM_PCT}\t${AA_HIGH}\t${AA_HIGH_PCT}\t${AA_VERY_HIGH}\t${AA_VERY_HIGH_PCT}" >> "$OUTPUT_TSV"
    
    # Update totals
    TOTAL_CALLABLE_BP=$((TOTAL_CALLABLE_BP + CALLABLE_BP))
    TOTAL_NONCALLABLE_BP=$((TOTAL_NONCALLABLE_BP + NONCALLABLE_BP))
    TOTAL_ALL_SITES=$((TOTAL_ALL_SITES + SITES_IN_VCF))
    TOTAL_BIALLELIC_SITES=$((TOTAL_BIALLELIC_SITES + BIALLELIC_SITES))
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
echo "Total sites (all): $TOTAL_ALL_SITES"
echo "Total biallelic sites: $TOTAL_BIALLELIC_SITES"
echo ""
echo "AA_prob Distribution (biallelic sites only):"
echo "  Missing: $TOTAL_AA_MISSING ($(awk "BEGIN {printf \"%.2f\", ($TOTAL_AA_MISSING/$TOTAL_BIALLELIC_SITES)*100}")%)"
echo "  < 0.5 (low): $TOTAL_AA_LOW ($(awk "BEGIN {printf \"%.2f\", ($TOTAL_AA_LOW/$TOTAL_BIALLELIC_SITES)*100}")%)"
echo "  0.5-0.8 (medium): $TOTAL_AA_MEDIUM ($(awk "BEGIN {printf \"%.2f\", ($TOTAL_AA_MEDIUM/$TOTAL_BIALLELIC_SITES)*100}")%)"
echo "  0.8-0.95 (high): $TOTAL_AA_HIGH ($(awk "BEGIN {printf \"%.2f\", ($TOTAL_AA_HIGH/$TOTAL_BIALLELIC_SITES)*100}")%)"
echo "  >= 0.95 (very high): $TOTAL_AA_VERY_HIGH ($(awk "BEGIN {printf \"%.2f\", ($TOTAL_AA_VERY_HIGH/$TOTAL_BIALLELIC_SITES)*100}")%)"
echo ""
echo "Output written to: $OUTPUT_TSV"
echo ""
echo "First 10 regions:"
head -11 "$OUTPUT_TSV" | column -t
echo ""
echo "=========================================="
echo "Done!"
echo "=========================================="

