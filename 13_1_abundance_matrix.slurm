#!/bin/bash -l
#SBATCH --job-name=abundance_matrix
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=24:00:00
#SBATCH --mem=64gb
#SBATCH --output=logs/step10_abundance_%j.out
#SBATCH --error=logs/step10_abundance_%j.err
#SBATCH --partition=batch,guest
#SBATCH --mail-user=echandrasekara2@unl.edu
#SBATCH --mail-type=END,FAIL

set -euo pipefail

module load coverm     # or: conda activate coverm

# =========================
# PATHS
# =========================
WORKDIR="/work/samodha/sachin/ShotgunM/Final_Analysis/data/combined_data/day1_25healthy_25infected/Output"
CATALOG_DIR="${WORKDIR}/11_NR_catalogs"
BAM_DIR="${WORKDIR}/12_Read_mapping"
OUTDIR="${WORKDIR}/13_Abundance_matrices"
THREADS=${SLURM_NTASKS_PER_NODE:-16}

mkdir -p "${OUTDIR}" logs

echo "======================================"
echo "STEP 10: Abundance matrices with CoverM"
echo "======================================"

# Function: build abundance matrix for one domain
build_matrix() {
    local DOMAIN=$1
    local CATALOG="${CATALOG_DIR}/${DOMAIN}/${DOMAIN}_NR_catalog.fa"
    local BAM_GLOB="${BAM_DIR}/${DOMAIN}/*_vs_${DOMAIN}.sorted.bam"
    local OUT_MATRIX="${OUTDIR}/${DOMAIN}_abundance_matrix.tsv"

    if [[ ! -f "$CATALOG" ]]; then
        echo "WARNING: No catalog for ${DOMAIN}. Skipping."
        return
    fi

    # Collect BAM files
    BAM_FILES=($(ls ${BAM_GLOB} 2>/dev/null))
    if [[ ${#BAM_FILES[@]} -eq 0 ]]; then
        echo "WARNING: No BAM files found for ${DOMAIN}. Skipping."
        return
    fi

    echo ""
    echo "--- ${DOMAIN}: ${#BAM_FILES[@]} samples ---"

    # CoverM produces:
    #   rows    = genes (from catalog)
    #   columns = samples (one per BAM)
    #   values  = TPM (transcripts per million) — best for cross-sample comparison
    coverm contig \
        --bam-files "${BAM_FILES[@]}" \
        --methods tpm \
        --min-read-percent-identity 95 \
        --min-read-aligned-percent 75 \
        --threads "${THREADS}" \
        > "$OUT_MATRIX"

    echo "  Matrix: $(wc -l < $OUT_MATRIX) rows × $(head -1 $OUT_MATRIX | awk '{print NF}') columns"
    echo "  Saved : ${OUT_MATRIX}"

    # Also produce read counts (useful for DESeq2/edgeR later)
    coverm contig \
        --bam-files "${BAM_FILES[@]}" \
        --methods count \
        --min-read-percent-identity 95 \
        --min-read-aligned-percent 75 \
        --threads "${THREADS}" \
        > "${OUTDIR}/${DOMAIN}_count_matrix.tsv"

    echo "  Count matrix saved: ${OUTDIR}/${DOMAIN}_count_matrix.tsv"
}

build_matrix "bacteria"
build_matrix "virus"
build_matrix "eukaryote"

echo ""
echo "======================================"
echo "STEP 10 COMPLETE"
echo "Output matrices:"
for f in "${OUTDIR}"/*.tsv; do
    echo "  $(basename $f): $(wc -l < $f) genes × $(( $(head -1 $f | awk '{print NF}') - 1 )) samples"
done
echo "======================================"
