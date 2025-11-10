#!/bin/bash
#SBATCH --job-name=diagnose-noncallable
#SBATCH --output=diagnose_noncallable_%j.out
#SBATCH --error=diagnose_noncallable_%j.err
#SBATCH --time=00:30:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --partition=general

# Diagnose potential issues with non-callable site calculation

VCF="/proj/johrilab/projects/DFEpos/dpgp3/vcf/dpgp3_2L_mono_biallelic.vcf.gz"
BED="/proj/johrilab/projects/DFEpos/dpgp3/masking/non_callable_sites_corrected.bed"

echo "=========================================="
echo "Diagnostic checks for non-callable sites"
echo "=========================================="
echo ""

# Check 1: Are positions actually consecutive in the VCF?
echo "CHECK 1: Position spacing in VCF"
echo "Looking at first 100 positions to see if they're consecutive..."
bcftools query -f '%POS\n' "$VCF" | head -100 | awk '
NR==1 {prev=$1; first=$1; print "First position:", $1; next}
{
    gap = $1 - prev
    if (gap > 1) {
        gaps++
        if (gaps <= 5) {
            printf "  Gap found: pos %d to %d (gap size: %d bp)\n", prev, $1, gap-1
        }
    }
    prev = $1
}
END {
    printf "Total gaps in first 100 sites: %d\n", gaps
    printf "Position range: %d to %d\n", first, prev
}
'

echo ""
echo "CHECK 2: VCF header - is this really an all-sites VCF?"
echo "Looking for evidence in VCF header..."
bcftools view -h "$VCF" | grep -i "all\|gvcf\|invariant\|monomorphic" | head -5

echo ""
echo "CHECK 3: Sample AC=0 sites (monomorphic)"
echo "These should be present in an all-sites VCF..."
bcftools query -f '%AC\n' "$VCF" | head -10000 | awk '
{
    if ($1 == 0) mono++
    else if ($1 > 0) variant++
    total++
}
END {
    printf "Monomorphic (AC=0): %d (%.1f%%)\n", mono, (mono/total)*100
    printf "Variant (AC>0): %d (%.1f%%)\n", variant, (variant/total)*100
    printf "Total sites checked: %d\n", total
}
'

echo ""
echo "CHECK 4: BED file - size distribution of non-callable regions"
awk '{print $3-$2}' "$BED" | sort -n | awk '
{
    size = $1
    sizes[size]++
    total++
    sum += size
    
    if (size == 1) one_bp++
    else if (size <= 5) small++
    else if (size <= 50) medium++
    else if (size <= 500) large++
    else huge++
}
END {
    printf "Non-callable region size distribution:\n"
    printf "  1 bp: %d (%.1f%%)\n", one_bp, (one_bp/total)*100
    printf "  2-5 bp: %d (%.1f%%)\n", small, (small/total)*100
    printf "  6-50 bp: %d (%.1f%%)\n", medium, (medium/total)*100
    printf "  51-500 bp: %d (%.1f%%)\n", large, (large/total)*100
    printf "  >500 bp: %d (%.1f%%)\n", huge, (huge/total)*100
    printf "\nTotal regions: %d\n", total
    printf "Average size: %.2f bp\n", sum/total
}
'

echo ""
echo "CHECK 5: Large non-callable regions (likely real biology)"
echo "Top 20 largest non-callable regions:"
awk '{print $1, $2, $3, $3-$2}' "$BED" | sort -k4,4nr | head -20

echo ""
echo "CHECK 6: Chromosome ends (telomeres)"
echo "First and last callable positions per chromosome:"
for CHR in chr2L chr2R chr3L chr3R; do
    FIRST=$(bcftools query -f '%POS\n' "/proj/johrilab/projects/DFEpos/dpgp3/vcf/dpgp3_${CHR#chr}_mono_biallelic.vcf.gz" | head -1)
    LAST=$(bcftools query -f '%POS\n' "/proj/johrilab/projects/DFEpos/dpgp3/vcf/dpgp3_${CHR#chr}_mono_biallelic.vcf.gz" | tail -1)
    echo "  $CHR: first=$FIRST, last=$LAST"
done

echo ""
echo "CHECK 7: Are there duplicate positions in VCF?"
echo "Checking for duplicate positions (would cause artificial gaps)..."
bcftools query -f '%POS\n' "$VCF" | head -10000 | sort | uniq -d | head -10
if [ $? -eq 0 ]; then
    DUP_COUNT=$(bcftools query -f '%POS\n' "$VCF" | head -10000 | sort | uniq -d | wc -l)
    echo "Duplicates found in first 10k sites: $DUP_COUNT"
else
    echo "No duplicates found"
fi

echo ""
echo "=========================================="
echo "Done!"
echo "=========================================="
