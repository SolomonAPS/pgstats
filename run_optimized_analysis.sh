#!/bin/bash
#SBATCH --job-name=many-stats
#SBATCH --output=many-stats_%j.out
#SBATCH --error=many-stats_%j.err
#SBATCH --time=06:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --partition=general

# Initialize conda for this shell session
eval "$(conda shell.bash hook)"

# Activate conda environment
conda activate many-stats

# Change to working directory
cd /proj/johrilab/projects/DFEpos/dpgp3/stats/

# Run many-stats analysis
many-stats stats \
    ../vcf/merged_biallelic_snps.vcf.gz \
    --bed uncallable_sites.bed \
    --bed-format non_callable \
    --regions-file regions.bed \
    --window-size 1000 \
    --step-size 1000 \
    --output results.csv \
    --keep-zarr \
    --zarr-dir /proj/johrilab/projects/DFEpos/dpgp3/stats

echo "Analysis completed successfully!"
