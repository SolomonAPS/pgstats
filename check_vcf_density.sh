#!/bin/bash
#SBATCH --job-name=check-vcf-density
#SBATCH --output=check_vcf_density_%j.out
#SBATCH --error=check_vcf_density_%j.err
#SBATCH --time=00:30:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --partition=general

# Check the density of sites in the VCF to understand what we're dealing with

VCF="/proj/johrilab/projects/DFEpos/dpgp3/vcf/dpgp3_2L_mono_biallelic.vcf.gz"

echo "=========================================="
echo "Analyzing VCF site density"
echo "=========================================="
echo ""
echo "VCF: $VCF"
echo ""

# Sample first 1000 sites and check spacing
echo "First 1000 sites - checking position spacing:"
bcftools query -f '%POS\n' "$VCF" | head -1000 | awk '
NR==1 {prev=$1; next}
{
    gap = $1 - prev - 1
    gaps[gap]++
    total_gap += gap
    prev = $1
}
END {
    print "Gap size distribution (first 1000 sites):"
    for (g in gaps) {
        printf "  Gap %d bp: %d occurrences (%.1f%%)\n", g, gaps[g], (gaps[g]/999)*100
    }
    printf "\nAverage gap: %.2f bp\n", total_gap/999
}
'

echo ""
echo "Checking AC (allele count) distribution:"
bcftools query -f '%AC\n' "$VCF" | head -10000 | awk '
{
    ac = $1
    if (ac == 0) zero++
    else if (ac == 1) one++
    else if (ac >= 2 && ac <= 5) low++
    else if (ac >= 6 && ac <= 10) mid++
    else high++
    total++
}
END {
    printf "AC=0 (monomorphic): %d (%.1f%%)\n", zero, (zero/total)*100
    printf "AC=1 (singleton): %d (%.1f%%)\n", one, (one/total)*100
    printf "AC=2-5: %d (%.1f%%)\n", low, (low/total)*100
    printf "AC=6-10: %d (%.1f%%)\n", mid, (mid/total)*100
    printf "AC>10: %d (%.1f%%)\n", high, (high/total)*100
}
'

echo ""
echo "First 20 sites (POS, REF, ALT, AC, AN):"
bcftools query -f '%POS\t%REF\t%ALT\t%AC\t%AN\n' "$VCF" | head -20

echo ""
echo "=========================================="
echo "Done!"
echo "=========================================="
