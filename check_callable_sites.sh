#!/bin/bash
# Script to manually verify callable sites calculation for specific windows

VCF="/proj/johrilab/projects/DFEpos/dpgp3/vcf/merged_biallelic_snps.vcf.gz"
BED="/proj/johrilab/projects/DFEpos/dpgp3/masking/non_callable_sites.bed"

echo "=========================================="
echo "Verifying Callable Sites Calculation"
echo "=========================================="
echo ""
echo "Testing 3 windows from region_1.csv:"
echo "  Window 1 (reasonable): chr3L:3594414-3594999 (54 variants, 155 callable sites reported)"
echo "  Window 2 (suspicious): chr3L:3595004-3595997 (176 variants, 176 callable sites reported)"
echo "  Window 3 (suspicious): chr3L:3596000-3596971 (155 variants, 155 callable sites reported)"
echo ""

# Window sizes
echo "=== Window Sizes (1-based inclusive) ==="
echo "Window 1: 3594414-3594999 = $((3594999 - 3594414 + 1)) bp"
echo "Window 2: 3595004-3595997 = $((3595997 - 3595004 + 1)) bp"
echo "Window 3: 3596000-3596971 = $((3596971 - 3596000 + 1)) bp"
echo ""

# Variant counts from VCF
echo "=== Variant Counts in VCF ==="
echo -n "Window 1 (chr3L:3594414-3594999): "
bcftools view -H -r chr3L:3594414-3594999 "$VCF" | wc -l

echo -n "Window 2 (chr3L:3595004-3595997): "
bcftools view -H -r chr3L:3595004-3595997 "$VCF" | wc -l

echo -n "Window 3 (chr3L:3596000-3596971): "
bcftools view -H -r chr3L:3596000-3596971 "$VCF" | wc -l
echo ""

# Callable sites calculation (BED is 0-based, so convert)
echo "=== Callable Sites Calculation ==="
echo "Window 1 (0-based: chr3L:3594413-3594999, length=586 bp):"
NONCALL=$(bedtools intersect -a <(echo -e "chr3L\t3594413\t3594999") -b "$BED" | awk '{sum += $3-$2} END {print sum+0}')
CALLABLE=$((586 - NONCALL))
echo "  Non-callable: $NONCALL bp"
echo "  Callable: $CALLABLE bp (reported: 155)"
echo ""

echo "Window 2 (0-based: chr3L:3595003-3595997, length=994 bp):"
NONCALL=$(bedtools intersect -a <(echo -e "chr3L\t3595003\t3595997") -b "$BED" | awk '{sum += $3-$2} END {print sum+0}')
CALLABLE=$((994 - NONCALL))
echo "  Non-callable: $NONCALL bp"
echo "  Callable: $CALLABLE bp (reported: 176)"
echo ""

echo "Window 3 (0-based: chr3L:3595999-3596971, length=972 bp):"
NONCALL=$(bedtools intersect -a <(echo -e "chr3L\t3595999\t3596971") -b "$BED" | awk '{sum += $3-$2} END {print sum+0}')
CALLABLE=$((972 - NONCALL))
echo "  Non-callable: $NONCALL bp"
echo "  Callable: $CALLABLE bp (reported: 155)"
echo ""

# Show actual BED overlaps
echo "=== BED Regions Overlapping Windows ==="
echo "Window 1 (chr3L:3594413-3594999):"
bedtools intersect -a <(echo -e "chr3L\t3594413\t3594999") -b "$BED" -wa -wb | head -5
COUNT=$(bedtools intersect -a <(echo -e "chr3L\t3594413\t3594999") -b "$BED" -wa -wb | wc -l)
if [ "$COUNT" -gt 5 ]; then
    echo "  ... and $((COUNT - 5)) more regions"
fi
echo ""

echo "Window 2 (chr3L:3595003-3595997):"
bedtools intersect -a <(echo -e "chr3L\t3595003\t3595997") -b "$BED" -wa -wb | head -5
COUNT=$(bedtools intersect -a <(echo -e "chr3L\t3595003\t3595997") -b "$BED" -wa -wb | wc -l)
if [ "$COUNT" -gt 5 ]; then
    echo "  ... and $((COUNT - 5)) more regions"
fi
echo ""

echo "Window 3 (chr3L:3595999-3596971):"
bedtools intersect -a <(echo -e "chr3L\t3595999\t3596971") -b "$BED" -wa -wb | head -5
COUNT=$(bedtools intersect -a <(echo -e "chr3L\t3595999\t3596971") -b "$BED" -wa -wb | wc -l)
if [ "$COUNT" -gt 5 ]; then
    echo "  ... and $((COUNT - 5)) more regions"
fi
echo ""

echo "=========================================="
echo "Done! Compare 'Callable bp' to 'reported' values."
echo "If they don't match, there's a bug in the Python code."
echo "=========================================="

