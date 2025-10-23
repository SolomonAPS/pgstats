#!/bin/bash
#SBATCH --job-name=many-stats-cluster
#SBATCH --output=many-stats_%j.out
#SBATCH --error=many-stats_%j.err
#SBATCH --time=04:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --partition=general

# Load modules
module load python
module load conda

# Initialize conda for this shell session
eval "$(conda shell.bash hook)"

# Activate conda environment
conda activate many-stats

# Change to working directory
cd /nas/longleaf/home/solsloat/stats

# Create Zarr cache directory in work space
ZARR_DIR="/work/users/s/o/solsloat/zarr_cache"
mkdir -p "$ZARR_DIR"

# Run many-stats analysis with persistent Zarr files
echo "Starting many-stats analysis..."
echo "Zarr files will be saved to: $ZARR_DIR"
echo "This enables fast re-runs without VCF conversion!"

many-stats stats \
    ../vcf/merged_biallelic_snps.vcf.gz \
    --bed uncallable_sites.bed \
    --bed-format callable \
    --regions-file regions.bed \
    --window-size 1000 \
    --step-size 1000 \
    --keep-zarr \
    --zarr-dir "$ZARR_DIR" \
    --output results.csv

echo "Analysis completed successfully!"
echo "Zarr files preserved in: $ZARR_DIR"
echo "Future runs will reuse these files for much faster processing!"
echo ""
echo "To reuse these Zarr files in future runs, use:"
echo "many-stats stats ../vcf/merged_biallelic_snps.vcf.gz --keep-zarr --zarr-dir $ZARR_DIR [other options]"
