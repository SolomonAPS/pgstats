#!/bin/bash
#SBATCH --job-name=create-zarr
#SBATCH --output=create-zarr_%j.out
#SBATCH --error=create-zarr_%j.err
#SBATCH --time=02:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1
#SBATCH --partition=general

# Initialize conda
eval "$(conda shell.bash hook)"
conda activate pgstats

# Set up directories
STATS_DIR="/proj/johrilab/projects/DFEpos/dpgp3/stats"
POLAR_DIR="/proj/johrilab/projects/DFEpos/dpgp3/test_adaptive_priors"
ZARR_DIR="${STATS_DIR}/zarr_cache"

# Create zarr directory
mkdir -p ${ZARR_DIR}

# Create a dummy region file (just to trigger Zarr creation)
echo -e "chr2L\t1\t1000" > /tmp/dummy_region.bed

echo "Creating Zarr file (this will take ~30 minutes)..."
pgstats stats \
    ${POLAR_DIR}/processed_11/dpgp3_polarized_biallelic.vcf.gz \
    --bed ${POLAR_DIR}/non_callable_sites_11.bed \
    --bed-format non_callable \
    --regions-file /tmp/dummy_region.bed \
    --window-size 1000 \
    --step-size 1000 \
    --output /tmp/dummy_results.csv \
    --stats theta_pi \
    --keep-zarr \
    --zarr-dir ${ZARR_DIR}

echo "Zarr created at: ${ZARR_DIR}/dpgp3_polarized_biallelic.vcz"
ls -lh ${ZARR_DIR}/

# Clean up dummy files
rm -f /tmp/dummy_region.bed /tmp/dummy_results.csv

