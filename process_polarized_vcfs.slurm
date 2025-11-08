#!/bin/bash
#SBATCH --job-name=process_polarized_vcfs
#SBATCH --output=process_polarized_vcfs_%j.out
#SBATCH --error=process_polarized_vcfs_%j.err
#SBATCH --time=6:00:00
#SBATCH --mem=16G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --partition=general

# Process polarized VCFs:
# 1. Remove outgroup samples (droSim1, droYak2)
# 2. Recalculate AC/AN
# 3. Keep AA and AA_prob, remove AA_info
# 4. Create separate VCFs: all sites, mono+biallelic only
# 5. Generate BED masks from all sites

set -e

# Load modules
module load bcftools
module load bedtools

# Configuration
BASE_DIR="/proj/johrilab/projects/DFEpos/dpgp3/test_adaptive_priors"
REF_FAI="/proj/johrilab/projects/DFEpos/dpgp3/reference/dmel_r5.fasta.fai"
OUTPUT_DIR="${BASE_DIR}/processed_11"
BED_OUTPUT="${BASE_DIR}/non_callable_sites_11.bed"

mkdir -p "$OUTPUT_DIR"

echo "=========================================="
echo "Processing Polarized VCFs (n=11)"
echo "=========================================="
echo ""
echo "Base directory: $BASE_DIR"
echo "Output directory: $OUTPUT_DIR"
echo "BED output: $BED_OUTPUT"
echo ""

# Process each chromosome
for CHR in 2L 2R 3L 3R; do
    INPUT_VCF="${BASE_DIR}/results_${CHR}_11/polarized_${CHR}_n11.vcf.gz"
    
    # Output files
    ALL_SITES_VCF="${OUTPUT_DIR}/dpgp3_${CHR}_polarized_clean.vcf.gz"
    MONO_BI_VCF="${OUTPUT_DIR}/dpgp3_${CHR}_polarized_mono_bi.vcf.gz"
    
    echo "=========================================="
    echo "Processing chromosome ${CHR}"
    echo "=========================================="
    echo "Input: $INPUT_VCF"
    
    if [ ! -f "$INPUT_VCF" ]; then
        echo "ERROR: Input VCF not found: $INPUT_VCF"
        continue
    fi
    
    # Step 1: Remove outgroup samples, clean INFO, and trim unused ALT alleles
    echo "  Step 1: Removing outgroups, cleaning INFO, and trimming unused ALT alleles..."
    bcftools view -s ^droSim1,droYak2 "$INPUT_VCF" | \
    bcftools annotate -x INFO/AA_info | \
    bcftools view -a -Ou | \
    bcftools +fill-tags -Oz -o "$ALL_SITES_VCF" -- -t AC,AN,AF
    
    bcftools index -f "$ALL_SITES_VCF"
    
    # Step 2: Create mono+biallelic only version
    echo "  Step 2: Creating mono+biallelic VCF..."
    bcftools view -m2 -M2 -Oz -o "$MONO_BI_VCF" "$ALL_SITES_VCF"
    bcftools index -f "$MONO_BI_VCF"
    
    # Check site counts
    ALL_SITES=$(bcftools view -H "$ALL_SITES_VCF" | wc -l)
    MONO_BI_SITES=$(bcftools view -H "$MONO_BI_VCF" | wc -l)
    MULTIALLELIC=$((ALL_SITES - MONO_BI_SITES))
    
    echo "  Total sites: $ALL_SITES"
    echo "  Mono+biallelic sites: $MONO_BI_SITES"
    echo "  Multiallelic sites (excluded): $MULTIALLELIC"
    
    # Check AA_prob distribution
    echo "  Checking AA_prob distribution..."
    bcftools query -f '%AA_prob\n' "$ALL_SITES_VCF" | \
    awk '{
        if ($1 == ".") missing++
        else if ($1 < 0.5) low++
        else if ($1 < 0.8) medium++
        else if ($1 < 0.95) high++
        else very_high++
        total++
    }
    END {
        print "    Missing AA_prob: " missing " (" (missing/total)*100 "%)"
        print "    AA_prob < 0.5: " low " (" (low/total)*100 "%)"
        print "    AA_prob 0.5-0.8: " medium " (" (medium/total)*100 "%)"
        print "    AA_prob 0.8-0.95: " high " (" (high/total)*100 "%)"
        print "    AA_prob >= 0.95: " very_high " (" (very_high/total)*100 "%)"
    }'
    
    echo ""
done

echo "=========================================="
echo "Creating Non-Callable Sites BED File"
echo "=========================================="
echo ""
echo "Using ALL sites (including multiallelic) to define callable regions"
echo ""

# Get chromosome lengths
declare -A CHR_LENGTHS
while read -r chrom length rest; do
    case "$chrom" in
        2L) CHR_LENGTHS["chr2L"]=$length ;;
        2R) CHR_LENGTHS["chr2R"]=$length ;;
        3L) CHR_LENGTHS["chr3L"]=$length ;;
        3R) CHR_LENGTHS["chr3R"]=$length ;;
    esac
done < "$REF_FAI"

TEMP_DIR=$(mktemp -d)

# Process each chromosome for BED creation
for CHR in 2L 2R 3L 3R; do
    ALL_SITES_VCF="${OUTPUT_DIR}/dpgp3_${CHR}_polarized_clean.vcf.gz"
    TEMP_BED="${TEMP_DIR}/${CHR}_noncallable.bed"
    CHROM_NAME="chr${CHR}"
    CHROM_LENGTH=${CHR_LENGTHS[$CHROM_NAME]}
    
    echo "Processing $CHROM_NAME for BED..."
    
    # Extract all positions (callable sites)
    bcftools view -H "$ALL_SITES_VCF" | \
    awk -v OFS="\t" '{
        pos = $2  # 1-based VCF position
        bed_start = pos - 1  # 0-based
        bed_end = pos        # half-open
        print $1, bed_start, bed_end
    }' | sort -k1,1 -k2,2n > "${TEMP_DIR}/${CHR}_callable.bed"
    
    # Create whole chromosome BED
    echo -e "${CHROM_NAME}\t0\t${CHROM_LENGTH}" > "${TEMP_DIR}/${CHR}_whole.bed"
    
    # Subtract callable from whole to get non-callable
    bedtools subtract -a "${TEMP_DIR}/${CHR}_whole.bed" -b "${TEMP_DIR}/${CHR}_callable.bed" > "$TEMP_BED"
    
    CALLABLE_SITES=$(wc -l < "${TEMP_DIR}/${CHR}_callable.bed")
    NONCALLABLE_REGIONS=$(wc -l < "$TEMP_BED")
    NONCALLABLE_BP=$(awk '{sum += $3-$2} END {print sum+0}' "$TEMP_BED")
    CALLABLE_BP=$((CHROM_LENGTH - NONCALLABLE_BP))
    
    echo "  Callable sites: $CALLABLE_SITES"
    echo "  Callable bp: $CALLABLE_BP ($(awk "BEGIN {printf \"%.2f\", ($CALLABLE_BP/$CHROM_LENGTH)*100}")%)"
    echo "  Non-callable regions: $NONCALLABLE_REGIONS"
    echo "  Non-callable bp: $NONCALLABLE_BP"
    echo ""
done

# Concatenate all chromosomes
cat "${TEMP_DIR}"/*_noncallable.bed | sort -k1,1 -k2,2n > "$BED_OUTPUT"

# Summary
TOTAL_REGIONS=$(wc -l < "$BED_OUTPUT")
TOTAL_NONCALLABLE_BP=$(awk '{sum += $3-$2} END {print sum+0}' "$BED_OUTPUT")

echo "=========================================="
echo "Summary"
echo "=========================================="
echo ""
echo "Total non-callable regions: $TOTAL_REGIONS"
echo "Total non-callable bases: $TOTAL_NONCALLABLE_BP"
echo ""
echo "Output files created:"
echo "  Clean VCFs (all sites): ${OUTPUT_DIR}/dpgp3_*_polarized_clean.vcf.gz"
echo "  Mono+biallelic VCFs: ${OUTPUT_DIR}/dpgp3_*_polarized_mono_bi.vcf.gz"
echo "  Non-callable BED: $BED_OUTPUT"
echo ""
echo "First 10 non-callable regions:"
head -10 "$BED_OUTPUT"
echo ""
echo "Last 10 non-callable regions:"
tail -10 "$BED_OUTPUT"
echo ""

# Cleanup
rm -rf "$TEMP_DIR"

echo "=========================================="
echo "Done!"
echo "=========================================="

