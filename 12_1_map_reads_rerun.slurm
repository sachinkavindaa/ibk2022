#!/bin/bash -l
#SBATCH --job-name=bwa_map
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=48:00:00
#SBATCH --mem=64gb
#SBATCH --output=logs/step9_map_%A_%a.out
#SBATCH --error=logs/step9_map_%A_%a.err
#SBATCH --partition=batch,guest
#SBATCH --mail-user=echandrasekara2@unl.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --array=0-49%10     # one task per sample; adjust 49 to (total_samples - 1)

set -euo pipefail

module load bwa
module load samtools

# =========================
# PATHS
# =========================
READS_DIR="/work/samodha/sachin/ShotgunM/Final_Analysis/data/combined_data/day1_25healthy_25infected/Output/04_Bovine_removal"
WORKDIR="/work/samodha/sachin/ShotgunM/Final_Analysis/data/combined_data/day1_25healthy_25infected/Output"
CATALOG_DIR="${WORKDIR}/11_NR_catalogs"
OUTDIR="${WORKDIR}/12_Read_mapping"
THREADS=${SLURM_NTASKS_PER_NODE:-16}

mkdir -p "${OUTDIR}"/{bacteria,virus,eukaryote} logs

# =========================
# PICK SAMPLE BY ARRAY INDEX
# =========================
mapfile -t R1_FILES < <(ls "${READS_DIR}"/*_final_clean_paired.1.fastq 2>/dev/null | sort)

TOTAL=${#R1_FILES[@]}
TASK_ID=${SLURM_ARRAY_TASK_ID}

echo "Total samples found : ${TOTAL}"
echo "This task ID        : ${TASK_ID}"

if [[ "$TASK_ID" -ge "$TOTAL" ]]; then
    echo "Task ID ${TASK_ID} exceeds sample count ${TOTAL}. Exiting."
    exit 0
fi

R1="${R1_FILES[$TASK_ID]}"
SAMPLE=$(basename "$R1" _final_clean_paired.1.fastq)
R2="${READS_DIR}/${SAMPLE}_final_clean_paired.2.fastq"

echo "======================================"
echo "STEP 9: BWA-MEM read mapping"
echo "Sample  : ${SAMPLE}"
echo "R1      : ${R1}"
echo "R2      : ${R2}"
echo "Threads : ${THREADS}"
echo "======================================"

if [[ ! -f "$R2" ]]; then
    echo "ERROR: R2 not found for ${SAMPLE}: ${R2}"
    exit 1
fi

# =========================
# MAPPING FUNCTION
# =========================
map_to_catalog() {
    local DOMAIN=$1
    local CATALOG="${CATALOG_DIR}/${DOMAIN}/${DOMAIN}_NR_catalog.fa"
    local DOMAIN_OUT="${OUTDIR}/${DOMAIN}"
    local BAM_OUT="${DOMAIN_OUT}/${SAMPLE}_vs_${DOMAIN}.sorted.bam"

    if [[ ! -f "$CATALOG" ]]; then
        echo "WARNING: No catalog found for ${DOMAIN} at ${CATALOG}. Skipping."
        return
    fi

    # Index catalog (only needs to happen once; safe to call every time)
    if [[ ! -f "${CATALOG}.bwt" ]]; then
        echo "Indexing ${DOMAIN} catalog (first time only)..."
        bwa index "$CATALOG"
    fi

    # Remove old zero-read BAMs from previous failed runs
    if [[ -f "$BAM_OUT" ]]; then
        EXISTING=$(samtools view -c "$BAM_OUT" 2>/dev/null || echo 0)
        if [[ "$EXISTING" -eq 0 ]]; then
            echo "${DOMAIN}: Removing empty BAM from previous run: ${BAM_OUT}"
            rm -f "$BAM_OUT" "${BAM_OUT}.bai"
        else
            echo "${DOMAIN}: Valid BAM already exists (${EXISTING} reads). Skipping."
            return
        fi
    fi

    echo "Mapping ${SAMPLE} → ${DOMAIN} catalog..."

    # Map reads with BWA-MEM
    # -F 4  : discard unmapped reads
    # keeping all mapped reads (paired and singleton) for gene catalog abundance
    bwa mem \
        -t "${THREADS}" \
        -M \
        -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ILLUMINA" \
        "$CATALOG" \
        "$R1" "$R2" \
    | samtools view -@ "${THREADS}" -bS -F 4 \
    | samtools sort -@ "${THREADS}" -o "$BAM_OUT"

    samtools index "$BAM_OUT"

    # Report mapping stats
    TOTAL_READS=$(samtools flagstat "$BAM_OUT" | grep "in total" | awk '{print $1}')
    MAPPED=$(samtools flagstat "$BAM_OUT" | grep "mapped (" | head -1 | awk '{print $1}')
    PCT=$(samtools flagstat "$BAM_OUT" | grep "mapped (" | head -1 | grep -oP '\(\K[^)]+')

    echo "  ${DOMAIN}: ${MAPPED} / ${TOTAL_READS} reads mapped (${PCT})"
    echo "  BAM: ${BAM_OUT}"

    if [[ "$MAPPED" -eq 0 ]]; then
        echo "  WARNING: 0 reads mapped for ${DOMAIN}. Check catalog is nucleotide:"
        echo "  head -2 ${CATALOG}"
    fi
}

# Map to all three domain catalogs
map_to_catalog "bacteria"
map_to_catalog "virus"
map_to_catalog "eukaryote"

echo "======================================"
echo "STEP 9 COMPLETE: ${SAMPLE}"
echo "======================================"
