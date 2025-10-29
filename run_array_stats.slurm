#!/bin/bash
#SBATCH --job-name=many-stats
#SBATCH --output=many-stats_%A_%a.out
#SBATCH --error=many-stats_%A_%a.err
#SBATCH --time=02:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1
#SBATCH --partition=general
#SBATCH --array=1-100

# Initialize conda
eval "$(conda shell.bash hook)"
conda activate many-stats

# Set up directories
STATS_DIR="/proj/johrilab/projects/DFEpos/dpgp3/stats"
TEMP_DIR="/work/users/s/o/solsloat/job_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"

# Create temp directory
mkdir -p ${TEMP_DIR}
cd ${TEMP_DIR}

# Extract this region
sed -n "${SLURM_ARRAY_TASK_ID}p" ${STATS_DIR}/regions.bed > region.bed

# Run analysis
many-stats stats \
    ${STATS_DIR}/../vcf/merged_biallelic_snps.vcf.gz \
    --bed ${STATS_DIR}/non_callable_sites.bed \
    --bed-format non_callable \
    --regions-file region.bed \
    --window-size 1000 \
    --step-size 1000 \
    --output results.csv \
    --keep-zarr \
    --stats all \
    --profile

# Save results
mkdir -p ${STATS_DIR}/array_results
cp results.csv ${STATS_DIR}/array_results/region_${SLURM_ARRAY_TASK_ID}.csv

# Cleanup
rm -rf ${TEMP_DIR}

