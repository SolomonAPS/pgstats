#!/bin/bash
#SBATCH --job-name=check_ancestral
#SBATCH --output=check_ancestral_%j.out
#SBATCH --error=check_ancestral_%j.err
#SBATCH --time=1:00:00
#SBATCH --mem=4G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=general

# Check how ancestral allele (AA) info is distributed in VCFs

set -e

# Load required modules
module load bcftools

VCF_DIR="/proj/johrilab/projects/DFEpos/dpgp3/vcf"

echo "=========================================="
echo "Checking Ancestral Allele (AA) Info"
echo "=========================================="
echo ""

for CHR in 2L 2R 3L 3R; do
    ORIGINAL_VCF="${VCF_DIR}/dpgp3_${CHR}.vcf.gz"
    FILTERED_VCF="${VCF_DIR}/dpgp3_${CHR}_mono_biallelic.vcf.gz"
    
    echo "=========================================="
    echo "Chromosome $CHR"
    echo "=========================================="
    
    if [ -f "$ORIGINAL_VCF" ]; then
        echo ""
        echo "Original VCF (dpgp3_${CHR}.vcf.gz):"
        
        TOTAL=$(bcftools view -H "$ORIGINAL_VCF" | wc -l)
        WITH_AA=$(bcftools view -H "$ORIGINAL_VCF" | grep -c "AA=" || true)
        WITHOUT_AA=$((TOTAL - WITH_AA))
        
        echo "  Total sites: $TOTAL"
        echo "  Sites with AA info: $WITH_AA ($(awk "BEGIN {printf \"%.2f\", ($WITH_AA / $TOTAL) * 100}")%)"
        echo "  Sites without AA info: $WITHOUT_AA ($(awk "BEGIN {printf \"%.2f\", ($WITHOUT_AA / $TOTAL) * 100}")%)"
        
        # Check if sites without AA are monomorphic
        echo ""
        echo "  Checking first 10 sites WITHOUT AA info:"
        bcftools view -H "$ORIGINAL_VCF" | grep -v "AA=" | head -10 | awk '{print "    " $1, $2, $4, $5, $8}'
        
        echo ""
        echo "  Checking first 10 sites WITH AA info:"
        bcftools view -H "$ORIGINAL_VCF" | grep "AA=" | head -10 | awk '{print "    " $1, $2, $4, $5, $8}'
    fi
    
    if [ -f "$FILTERED_VCF" ]; then
        echo ""
        echo "Filtered VCF (dpgp3_${CHR}_mono_biallelic.vcf.gz):"
        
        TOTAL=$(bcftools view -H "$FILTERED_VCF" | wc -l)
        WITH_AA=$(bcftools view -H "$FILTERED_VCF" | grep -c "AA=" || true)
        WITHOUT_AA=$((TOTAL - WITH_AA))
        
        echo "  Total sites: $TOTAL"
        echo "  Sites with AA info: $WITH_AA ($(awk "BEGIN {printf \"%.2f\", ($WITH_AA / $TOTAL) * 100}")%)"
        echo "  Sites without AA info: $WITHOUT_AA ($(awk "BEGIN {printf \"%.2f\", ($WITHOUT_AA / $TOTAL) * 100}")%)"
    fi
    
    echo ""
done

echo "=========================================="
echo "Summary"
echo "=========================================="
echo ""
echo "Key questions to answer:"
echo "1. Do monomorphic sites (ALT='.') have AA info?"
echo "2. Do polymorphic sites always have AA info?"
echo "3. If polymorphic sites lack AA, can they still be polarized?"
echo ""
echo "Based on the output above:"
echo "- If monomorphic sites lack AA but polymorphic sites have it,"
echo "  then filtering for AA removed only invariant sites (expected)"
echo "- If some polymorphic sites also lack AA, you may have removed"
echo "  unpolarizable variants (which is correct for DFE analysis)"
echo ""
echo "=========================================="

