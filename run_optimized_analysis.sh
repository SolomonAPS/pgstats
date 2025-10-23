#!/bin/bash
#SBATCH --job-name=many-stats-optimized
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

# Run many-stats analysis with optimizations
echo "Starting optimized many-stats analysis..."
echo "Using --keep-zarr for faster future runs"

many-stats stats \
    ../vcf/merged_biallelic_snps.vcf.gz \
    --bed uncallable_sites.bed \
    --bed-format callable \
    --regions-file regions.bed \
    --window-size 1000 \
    --step-size 1000 \
    --keep-zarr \
    --output results.csv

echo "Analysis completed successfully!"
echo "Zarr files preserved for faster future runs"
