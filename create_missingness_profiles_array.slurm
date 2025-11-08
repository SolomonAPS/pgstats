#!/bin/bash
#SBATCH --job-name=missingness_profiles
#SBATCH --output=logs/missingness_profile_%A_%a.out
#SBATCH --error=logs/missingness_profile_%A_%a.err
#SBATCH --time=2:00:00
#SBATCH --mem=8G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=general
#SBATCH --array=1-100

# Create missingness profiles for all 100 regions
# Each array job processes one region

set -e

# Configuration
REGIONS_FILE="regions.bed"
VCF_DIR="/proj/johrilab/projects/DFEpos/dpgp3/test_adaptive_priors/processed_11"
OUTPUT_DIR="/proj/johrilab/projects/DFEpos/dpgp3/test_adaptive_priors/missingness_profiles"
SCRIPT_DIR="$HOME/bin/many-stats"

# Create output directory
mkdir -p "$OUTPUT_DIR"
mkdir -p logs

# Get the region for this array task
REGION_LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$REGIONS_FILE")

if [ -z "$REGION_LINE" ]; then
    echo "ERROR: Could not read region line ${SLURM_ARRAY_TASK_ID} from ${REGIONS_FILE}"
    exit 1
fi

# Parse BED line (0-based, half-open)
CHROM=$(echo "$REGION_LINE" | awk '{print $1}')
START_0BASED=$(echo "$REGION_LINE" | awk '{print $2}')
END_0BASED=$(echo "$REGION_LINE" | awk '{print $3}')

# Convert to 1-based inclusive for VCF region query
START_1BASED=$((START_0BASED + 1))
END_1BASED=$((END_0BASED))

# Create region string for VCF query (1-based, inclusive)
REGION_STR="${CHROM}:${START_1BASED}-${END_1BASED}"

# Determine which VCF file to use based on chromosome
case "$CHROM" in
    chr2L)
        VCF_FILE="${VCF_DIR}/dpgp3_2L_polarized_mono_bi.vcf.gz"
        ;;
    chr2R)
        VCF_FILE="${VCF_DIR}/dpgp3_2R_polarized_mono_bi.vcf.gz"
        ;;
    chr3L)
        VCF_FILE="${VCF_DIR}/dpgp3_3L_polarized_mono_bi.vcf.gz"
        ;;
    chr3R)
        VCF_FILE="${VCF_DIR}/dpgp3_3R_polarized_mono_bi.vcf.gz"
        ;;
    *)
        echo "ERROR: Unknown chromosome: $CHROM"
        exit 1
        ;;
esac

# Check if VCF exists
if [ ! -f "$VCF_FILE" ]; then
    echo "ERROR: VCF file not found: $VCF_FILE"
    exit 1
fi

# Output file with region info
OUTPUT_FILE="${OUTPUT_DIR}/region_${SLURM_ARRAY_TASK_ID}_${CHROM}_${START_1BASED}_${END_1BASED}.npz"

echo "=========================================="
echo "Creating Missingness Profile"
echo "=========================================="
echo "Array task: ${SLURM_ARRAY_TASK_ID} / 100"
echo "Region: ${REGION_STR}"
echo "  BED (0-based): ${CHROM}:${START_0BASED}-${END_0BASED}"
echo "  VCF (1-based): ${REGION_STR}"
echo "VCF file: $VCF_FILE"
echo "Output: $OUTPUT_FILE"
echo ""

# Run the script
python3 "${SCRIPT_DIR}/create_missingness_profile.py" \
    --vcf "$VCF_FILE" \
    --region "$REGION_STR" \
    --output "$OUTPUT_FILE"

echo ""
echo "=========================================="
echo "Done!"
echo "=========================================="

