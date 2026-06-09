#!/bin/bash -l
#SBATCH --job-name=megahit_twopass
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=168:00:00
#SBATCH --mem=64gb
#SBATCH --output=logs/megahit_%A_%a.out
#SBATCH --error=logs/megahit_%A_%a.err
#SBATCH --partition=yinlab,batch,guest
#SBATCH --mail-user=echandrasekara2@unl.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --array=0-44%10          # adjust 49 to (total_samples - 1)

set -euo pipefail

module load megahit/1.2

# =========================
# PATHS
# =========================
BASE="/work/samodha/sachin/ShotgunM/Final_Analysis/data/combined_data/day1_25healthy_25infected/Output/04_Bovine_removal"
OUTDIR="/work/samodha/sachin/ShotgunM/Final_Analysis/data/combined_data/day1_25healthy_25infected/Output/05_Assembly_per_sample"
THREADS=${SLURM_NTASKS_PER_NODE:-16}

mkdir -p "${OUTDIR}" logs

# =========================
# PICK SAMPLE BY ARRAY INDEX
# =========================
mapfile -t R1_FILES < <(ls "${BASE}"/*_final_clean_paired.1.fastq 2>/dev/null | sort)

TOTAL=${#R1_FILES[@]}
if [[ "$TOTAL" -eq 0 ]]; then
    echo "ERROR: No R1 files found in ${BASE}"
    exit 1
fi

TASK_ID=${SLURM_ARRAY_TASK_ID}
if [[ "$TASK_ID" -ge "$TOTAL" ]]; then
    echo "Task ID ${TASK_ID} exceeds sample count ${TOTAL}. Nothing to do."
    exit 0
fi

# =========================
# DEFINE FILE PATHS
# =========================
R1="${R1_FILES[$TASK_ID]}"
SAMPLE=$(basename "$R1" _final_clean_paired.1.fastq)

R2="${BASE}/${SAMPLE}_final_clean_paired.2.fastq"
SINGLE="${BASE}/${SAMPLE}_final_clean_single.fastq"

# Per-sample working directory
SAMPLE_DIR="${OUTDIR}/${SAMPLE}"
mkdir -p "${SAMPLE_DIR}"

# Pass 1 paths
PASS1_DIR="${SAMPLE_DIR}/pass1_paired"
PAIRED_CONTIGS="${SAMPLE_DIR}/paired_contigs.fa"

# Pass 2 paths
PASS2_DIR="${SAMPLE_DIR}/pass2_final"
FINAL_CONTIGS="${SAMPLE_DIR}/final_contigs.fa"

# =========================
# LOGGING
# =========================
echo "======================================"
echo "SLURM Array Task : ${TASK_ID} / $((TOTAL - 1))"
echo "Sample           : ${SAMPLE}"
echo "R1               : ${R1}"
echo "R2               : ${R2}"
echo "Singleton        : $([ -f "$SINGLE" ] && echo "$SINGLE" || echo "NOT FOUND")"
echo "Sample dir       : ${SAMPLE_DIR}"
echo "Threads          : ${THREADS}"
echo "======================================"

# =========================
# VALIDATE INPUTS
# =========================
if [[ ! -f "$R2" ]]; then
    echo "ERROR: Missing R2 for ${SAMPLE}. Exiting."
    exit 1
fi

# =========================
# RESUME CHECK — skip entire sample if final output exists
# =========================
if [[ -f "$FINAL_CONTIGS" ]]; then
    echo "Final contigs already exist for ${SAMPLE}. Skipping."
    exit 0
fi

# ==================================================
# PASS 1 — Paired-end assembly only
# ==================================================
echo ""
echo "--- PASS 1: Paired-end assembly ---"

if [[ -f "$PAIRED_CONTIGS" ]]; then
    echo "Pass 1 already complete. Reusing ${PAIRED_CONTIGS}"
else
    # Clean up any incomplete pass 1 directory
    [[ -d "$PASS1_DIR" ]] && rm -rf "$PASS1_DIR"

    megahit \
        -1 "$R1" \
        -2 "$R2" \
        --k-min 21 \
        --k-max 141 \
        --k-step 12 \
        --min-contig-len 500 \
        -t "$THREADS" \
        -o "$PASS1_DIR"

    if [[ ! -f "${PASS1_DIR}/final.contigs.fa" ]]; then
        echo "ERROR: Pass 1 MEGAHIT failed for ${SAMPLE}."
        exit 1
    fi

    cp "${PASS1_DIR}/final.contigs.fa" "$PAIRED_CONTIGS"
    PASS1_COUNT=$(grep -c ">" "$PAIRED_CONTIGS")
    echo "Pass 1 complete. Contigs: ${PASS1_COUNT}"
fi

# ==================================================
# PASS 2 — paired_contigs.fa + singleton reads → final_contigs.fa
# ==================================================
echo ""
echo "--- PASS 2: paired_contigs + singletons ---"

if [[ ! -f "$SINGLE" ]]; then
    echo "WARNING: No singleton file found for ${SAMPLE}."
    echo "Copying Pass 1 contigs directly as final output."
    cp "$PAIRED_CONTIGS" "$FINAL_CONTIGS"
else
    # Clean up any incomplete pass 2 directory
    [[ -d "$PASS2_DIR" ]] && rm -rf "$PASS2_DIR"

    # Feed Pass 1 contigs as -r (unpaired/read input) alongside singletons
    # Both are treated as single-end reads in this pass
    megahit \
        -r "${PAIRED_CONTIGS},${SINGLE}" \
        --k-min 21 \
        --k-max 141 \
        --k-step 12 \
        --min-contig-len 500 \
        -t "$THREADS" \
        -o "$PASS2_DIR"

    if [[ ! -f "${PASS2_DIR}/final.contigs.fa" ]]; then
        echo "ERROR: Pass 2 MEGAHIT failed for ${SAMPLE}."
        exit 1
    fi

    cp "${PASS2_DIR}/final.contigs.fa" "$FINAL_CONTIGS"
fi

# =========================
# SUMMARY
# =========================
PASS1_COUNT=$(grep -c ">" "$PAIRED_CONTIGS")
FINAL_COUNT=$(grep -c ">" "$FINAL_CONTIGS")

echo ""
echo "======================================"
echo "DONE: ${SAMPLE}"
echo "Pass 1 contigs : ${PASS1_COUNT}  →  ${PAIRED_CONTIGS}"
echo "Final contigs  : ${FINAL_COUNT}  →  ${FINAL_CONTIGS}"
echo "======================================"
