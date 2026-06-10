#!/bin/bash -l
#SBATCH --job-name=bac_abundance
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=24:00:00
#SBATCH --mem=64gb
#SBATCH --output=logs/step10_bacteria_abundance_%j.out
#SBATCH --error=logs/step10_bacteria_abundance_%j.err
#SBATCH --partition=batch,guest
#SBATCH --mail-user=echandrasekara2@unl.edu
#SBATCH --mail-type=END,FAIL

set -euo pipefail

module load coverm

WORKDIR="/work/samodha/sachin/ShotgunM/Final_Analysis/data/combined_data/day1_25healthy_25infected/Output"

CATALOG="${WORKDIR}/11_NR_catalogs/bacteria/bacteria_NR_catalog.fa"
BAM_DIR="${WORKDIR}/12_Read_mapping/bacteria2"
OUTDIR="${WORKDIR}/13_Abundance_matrices"
THREADS=${SLURM_NTASKS_PER_NODE:-16}

mkdir -p "$OUTDIR" logs

echo "======================================"
echo "Bacteria-only abundance matrix rerun"
echo "Catalog : $CATALOG"
echo "BAM dir : $BAM_DIR"
echo "Output  : $OUTDIR"
echo "======================================"

BAM_FILES=($(ls ${BAM_DIR}/*_vs_bacteria.sorted.bam 2>/dev/null))

echo "BAM files found: ${#BAM_FILES[@]}"

if [[ ${#BAM_FILES[@]} -eq 0 ]]; then
    echo "ERROR: No bacteria BAM files found."
    exit 1
fi

if [[ ${#BAM_FILES[@]} -ne 45 ]]; then
    echo "WARNING: Expected 45 BAM files, found ${#BAM_FILES[@]}"
fi

echo "Building bacteria TPM matrix..."

coverm contig \
    --bam-files "${BAM_FILES[@]}" \
    --methods tpm \
    --min-read-percent-identity 95 \
    --min-read-aligned-percent 75 \
    --threads "$THREADS" \
    > "${OUTDIR}/bacteria_abundance_matrix.tsv"

echo "Building bacteria count matrix..."

coverm contig \
    --bam-files "${BAM_FILES[@]}" \
    --methods count \
    --min-read-percent-identity 95 \
    --min-read-aligned-percent 75 \
    --threads "$THREADS" \
    > "${OUTDIR}/bacteria_count_matrix.tsv"

echo "======================================"
echo "Done."
echo "TPM matrix:"
wc -l "${OUTDIR}/bacteria_abundance_matrix.tsv"
echo "Count matrix:"
wc -l "${OUTDIR}/bacteria_count_matrix.tsv"
echo "======================================"
