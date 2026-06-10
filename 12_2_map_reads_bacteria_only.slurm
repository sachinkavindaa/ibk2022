#!/bin/bash -l
#SBATCH --job-name=bwa_bacteria
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=48:00:00
#SBATCH --mem=64gb
#SBATCH --output=logs/step9_bacteria_%A_%a.out
#SBATCH --error=logs/step9_bacteria_%A_%a.err
#SBATCH --partition=batch,guest
#SBATCH --mail-user=echandrasekara2@unl.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --array=0-49%10

set -euo pipefail

module load bwa
module load samtools

READS_DIR="/work/samodha/sachin/ShotgunM/Final_Analysis/data/combined_data/day1_25healthy_25infected/Output/04_Bovine_removal"
WORKDIR="/work/samodha/sachin/ShotgunM/Final_Analysis/data/combined_data/day1_25healthy_25infected/Output"

CATALOG="${WORKDIR}/11_NR_catalogs/bacteria/bacteria_NR_catalog.fa"
OUTDIR="${WORKDIR}/12_Read_mapping/bacteria2"
THREADS=${SLURM_NTASKS_PER_NODE:-16}

mkdir -p "$OUTDIR" logs

mapfile -t R1_FILES < <(ls "${READS_DIR}"/*_final_clean_paired.1.fastq 2>/dev/null | sort)

TOTAL=${#R1_FILES[@]}
TASK_ID=${SLURM_ARRAY_TASK_ID}

if [[ "$TASK_ID" -ge "$TOTAL" ]]; then
    echo "Task ID ${TASK_ID} exceeds sample count ${TOTAL}. Exiting."
    exit 0
fi

R1="${R1_FILES[$TASK_ID]}"
SAMPLE=$(basename "$R1" _final_clean_paired.1.fastq)
R2="${READS_DIR}/${SAMPLE}_final_clean_paired.2.fastq"

BAM_OUT="${OUTDIR}/${SAMPLE}_vs_bacteria.sorted.bam"

echo "Sample: $SAMPLE"
echo "Catalog: $CATALOG"

if [[ ! -f "$R2" ]]; then
    echo "ERROR: R2 not found: $R2"
    exit 1
fi

if [[ ! -f "${CATALOG}.bwt" ]]; then
    bwa index "$CATALOG"
fi

bwa mem \
    -t "$THREADS" \
    -M \
    -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ILLUMINA" \
    "$CATALOG" \
    "$R1" "$R2" \
| samtools view -@ "$THREADS" -bS -F 4 \
| samtools sort -@ "$THREADS" -o "$BAM_OUT"

samtools index "$BAM_OUT"

samtools flagstat "$BAM_OUT"

echo "Done: $BAM_OUT"
