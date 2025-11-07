#!/bin/bash
#SBATCH --job-name=compare_vcf_coverage
#SBATCH --output=compare_vcf_coverage_%j.out
#SBATCH --error=compare_vcf_coverage_%j.err
#SBATCH --time=2:00:00
#SBATCH --mem=4G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=general

# Compare site counts between original and filtered VCFs

set -e

# Load required modules
module load bcftools

VCF_DIR="/proj/johrilab/projects/DFEpos/dpgp3/vcf"
REF_FAI="/proj/johrilab/projects/DFEpos/dpgp3/reference/dmel_r5.fasta.fai"

echo "=========================================="
echo "Comparing VCF Coverage"
echo "=========================================="
echo ""

# Get chromosome lengths
declare -A CHR_LENGTHS
while read -r chrom length rest; do
    case "$chrom" in
        2L) CHR_LENGTHS["2L"]=$length ;;
        2R) CHR_LENGTHS["2R"]=$length ;;
        3L) CHR_LENGTHS["3L"]=$length ;;
        3R) CHR_LENGTHS["3R"]=$length ;;
    esac
done < "$REF_FAI"

echo "Chromosome | Length (bp) | Original Sites | Filtered Sites | Original % | Filtered %"
echo "-----------|-------------|----------------|----------------|------------|------------"

for CHR in 2L 2R 3L 3R; do
    CHROM_LENGTH=${CHR_LENGTHS[$CHR]}
    
    ORIGINAL_VCF="${VCF_DIR}/dpgp3_${CHR}.vcf.gz"
    FILTERED_VCF="${VCF_DIR}/dpgp3_${CHR}_mono_biallelic.vcf.gz"
    
    if [ -f "$ORIGINAL_VCF" ]; then
        ORIGINAL_SITES=$(bcftools view -H "$ORIGINAL_VCF" | wc -l)
        ORIGINAL_PCT=$(awk "BEGIN {printf \"%.2f\", ($ORIGINAL_SITES / $CHROM_LENGTH) * 100}")
    else
        ORIGINAL_SITES="N/A"
        ORIGINAL_PCT="N/A"
    fi
    
    if [ -f "$FILTERED_VCF" ]; then
        FILTERED_SITES=$(bcftools view -H "$FILTERED_VCF" | wc -l)
        FILTERED_PCT=$(awk "BEGIN {printf \"%.2f\", ($FILTERED_SITES / $CHROM_LENGTH) * 100}")
    else
        FILTERED_SITES="N/A"
        FILTERED_PCT="N/A"
    fi
    
    printf "%-10s | %11s | %14s | %14s | %9s%% | %9s%%\n" \
        "$CHR" "$CHROM_LENGTH" "$ORIGINAL_SITES" "$FILTERED_SITES" "$ORIGINAL_PCT" "$FILTERED_PCT"
done

echo ""
echo "=========================================="
echo "Site Retention After Filtering"
echo "=========================================="
echo ""

for CHR in 2L 2R 3L 3R; do
    ORIGINAL_VCF="${VCF_DIR}/dpgp3_${CHR}.vcf.gz"
    FILTERED_VCF="${VCF_DIR}/dpgp3_${CHR}_mono_biallelic.vcf.gz"
    
    if [ -f "$ORIGINAL_VCF" ] && [ -f "$FILTERED_VCF" ]; then
        ORIGINAL_SITES=$(bcftools view -H "$ORIGINAL_VCF" | wc -l)
        FILTERED_SITES=$(bcftools view -H "$FILTERED_VCF" | wc -l)
        
        if [ "$ORIGINAL_SITES" -gt 0 ]; then
            RETENTION=$(awk "BEGIN {printf \"%.2f\", ($FILTERED_SITES / $ORIGINAL_SITES) * 100}")
            REMOVED=$((ORIGINAL_SITES - FILTERED_SITES))
            echo "$CHR: Kept $FILTERED_SITES / $ORIGINAL_SITES sites (${RETENTION}% retained, $REMOVED removed)"
        fi
    fi
done

echo ""
echo "=========================================="
echo "First 20 sites in original vs filtered"
echo "=========================================="
echo ""

CHR="2L"
ORIGINAL_VCF="${VCF_DIR}/dpgp3_${CHR}.vcf.gz"
FILTERED_VCF="${VCF_DIR}/dpgp3_${CHR}_mono_biallelic.vcf.gz"

echo "Original VCF (first 20 positions):"
bcftools view -H "$ORIGINAL_VCF" | head -20 | awk '{print $1, $2, $4, $5, $7}'

echo ""
echo "Filtered VCF (first 20 positions):"
bcftools view -H "$FILTERED_VCF" | head -20 | awk '{print $1, $2, $4, $5, $7}'

echo ""
echo "=========================================="
echo "Done!"
echo "=========================================="

