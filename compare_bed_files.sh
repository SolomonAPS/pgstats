#!/bin/bash
# Script to compare old vs new BED file callable sites calculations

VCF="/proj/johrilab/projects/DFEpos/dpgp3/vcf/merged_biallelic_snps.vcf.gz"
OLD_BED="/proj/johrilab/projects/DFEpos/dpgp3/masking/non_callable_sites.bed"
NEW_BED="/proj/johrilab/projects/DFEpos/dpgp3/masking/non_callable_sites_corrected.bed"

echo "=========================================="
echo "Comparing Old vs New BED Files"
echo "=========================================="
echo ""
echo "Old BED: $OLD_BED"
echo "New BED: $NEW_BED"
echo ""

# Test windows from region_1.csv
echo "Testing 3 windows from region_1.csv:"
echo "  Window 1: chr3L:3594414-3594999 (54 variants, 155 callable reported)"
echo "  Window 2: chr3L:3595004-3595997 (176 variants, 176 callable reported)"
echo "  Window 3: chr3L:3596000-3596971 (155 variants, 155 callable reported)"
echo ""

# Function to calculate callable sites for a window
calc_callable() {
    local chrom=$1
    local start_1based=$2
    local end_1based=$3
    local bed_file=$4
    local window_size=$((end_1based - start_1based + 1))
    
    # Convert to 0-based for BED
    local start_0based=$((start_1based - 1))
    local end_0based=$end_1based
    
    # Calculate non-callable bases
    local noncall=$(bedtools intersect -a <(echo -e "${chrom}\t${start_0based}\t${end_0based}") -b "$bed_file" | \
                    awk '{sum += $3-$2} END {print sum+0}')
    
    local callable=$((window_size - noncall))
    
    echo "$callable"
}

echo "=== Window 1: chr3L:3594414-3594999 (586 bp) ==="
echo "Reported callable: 155 bp"
OLD_CALL=$(calc_callable "chr3L" 3594414 3594999 "$OLD_BED")
NEW_CALL=$(calc_callable "chr3L" 3594414 3594999 "$NEW_BED")
echo "Old BED callable: $OLD_CALL bp"
echo "New BED callable: $NEW_CALL bp"
if [ "$NEW_CALL" -eq 155 ]; then
    echo "✓ NEW BED MATCHES!"
else
    echo "✗ Still doesn't match"
fi
echo ""

echo "=== Window 2: chr3L:3595004-3595997 (994 bp) ==="
echo "Reported callable: 176 bp"
OLD_CALL=$(calc_callable "chr3L" 3595004 3595997 "$OLD_BED")
NEW_CALL=$(calc_callable "chr3L" 3595004 3595997 "$NEW_BED")
echo "Old BED callable: $OLD_CALL bp"
echo "New BED callable: $NEW_CALL bp"
if [ "$NEW_CALL" -eq 176 ]; then
    echo "✓ NEW BED MATCHES!"
else
    echo "✗ Still doesn't match"
fi
echo ""

echo "=== Window 3: chr3L:3596000-3596971 (972 bp) ==="
echo "Reported callable: 155 bp"
OLD_CALL=$(calc_callable "chr3L" 3596000 3596971 "$OLD_BED")
NEW_CALL=$(calc_callable "chr3L" 3596000 3596971 "$NEW_BED")
echo "Old BED callable: $OLD_CALL bp"
echo "New BED callable: $NEW_CALL bp"
if [ "$NEW_CALL" -eq 155 ]; then
    echo "✓ NEW BED MATCHES!"
else
    echo "✗ Still doesn't match"
fi
echo ""

echo "=========================================="
echo "BED File Statistics"
echo "=========================================="
echo ""
echo "Old BED:"
echo "  Regions: $(wc -l < "$OLD_BED")"
echo "  Total non-callable bp: $(awk '{sum += $3-$2} END {print sum+0}' "$OLD_BED")"
echo ""
echo "New BED:"
echo "  Regions: $(wc -l < "$NEW_BED")"
echo "  Total non-callable bp: $(awk '{sum += $3-$2} END {print sum+0}' "$NEW_BED")"
echo ""
echo "=========================================="

