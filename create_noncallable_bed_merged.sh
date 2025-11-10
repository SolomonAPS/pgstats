#!/bin/bash

# Alternative approach: Merge consecutive callable positions first
# This will show if the gaps are real or an artifact of the method

VCF="/proj/johrilab/projects/DFEpos/dpgp3/vcf/dpgp3_2L_mono_biallelic.vcf.gz"
CHR_LENGTH=23011544
OUTPUT="/proj/johrilab/projects/DFEpos/dpgp3/masking/test_merged_approach.bed"

echo "Testing merged approach for chr2L..."

# Extract positions and convert to BED, then MERGE consecutive positions
bcftools view -H "$VCF" | \
awk -v OFS="\t" '{
    pos = $2
    bed_start = pos - 1
    bed_end = pos
    print $1, bed_start, bed_end
}' | \
bedtools merge -i - > /tmp/callable_merged.bed

echo "Callable regions after merging:"
wc -l /tmp/callable_merged.bed
head -20 /tmp/callable_merged.bed

echo ""
echo "Size distribution of callable regions:"
awk '{print $3-$2}' /tmp/callable_merged.bed | sort -n | awk '
{
    size = $1
    if (size == 1) one++
    else if (size <= 10) small++
    else if (size <= 100) medium++
    else if (size <= 1000) large++
    else huge++
    total++
    sum += size
}
END {
    printf "1 bp: %d\n", one
    printf "2-10 bp: %d\n", small
    printf "11-100 bp: %d\n", medium
    printf "101-1000 bp: %d\n", large
    printf ">1000 bp: %d\n", huge
    printf "Total callable regions: %d\n", total
    printf "Average size: %.2f bp\n", sum/total
}
'

# Now subtract to get non-callable
echo -e "chr2L\t0\t${CHR_LENGTH}" > /tmp/whole_chr.bed
bedtools subtract -a /tmp/whole_chr.bed -b /tmp/callable_merged.bed > /tmp/noncallable_merged.bed

echo ""
echo "Non-callable regions after merging callable sites:"
wc -l /tmp/noncallable_merged.bed
NONCALLABLE_BP=$(awk '{sum += $3-$2} END {print sum}' /tmp/noncallable_merged.bed)
echo "Non-callable bp: $NONCALLABLE_BP"

echo ""
echo "First 20 non-callable regions:"
head -20 /tmp/noncallable_merged.bed

echo ""
echo "Largest 20 non-callable regions:"
awk '{print $1, $2, $3, $3-$2}' /tmp/noncallable_merged.bed | sort -k4,4nr | head -20
