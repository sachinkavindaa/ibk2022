#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=6
#SBATCH --mem=64gb
#SBATCH --partition=guest,batch
#SBATCH --time=168:00:00
#SBATCH --job-name=HostRemoval_day1
#SBATCH --error=logs/HostRemoval_%j.err
#SBATCH --output=logs/HostRemoval_%j.out
#SBATCH --mail-user=echandrasekara2@unl.edu
#SBATCH --mail-type=END,FAIL

set -euo pipefail

module load bowtie/2.5

# ── Paths ──────────────────────────────────────────────────────────────────────
BASE_DIR="/work/samodha/sachin/ShotgunM/Final_Analysis/data/combined_data"
DATA_DIR="$BASE_DIR/day1_25healthy_25infected"
OUT_DIR="$DATA_DIR/Output"
TRIM_DIR="$OUT_DIR/02_Trimming"
HUMAN_DIR="$OUT_DIR/03_Human_removal"
BOVINE_DIR="$OUT_DIR/04_Bovine_removal"
STATS_DIR="$OUT_DIR/host_removal_stats"

HG38_INDEX="/work/samodha/sachin/hg38_index/genomeo"        # ← fixed typo (genomeo)
BOVINE_INDEX="/work/samodha/sachin/bovine_index/bosTau9"

mkdir -p "$HUMAN_DIR" "$BOVINE_DIR" "$STATS_DIR" logs

echo "Starting host removal..."

# ── Loop ───────────────────────────────────────────────────────────────────────
for FW in "$TRIM_DIR"/*_1_trimmed.fq; do

    [[ -e "$FW" ]] || { echo "No trimmed files found in $TRIM_DIR"; exit 1; }

    sample=$(basename "$FW" _1_trimmed.fq)
    RV="$TRIM_DIR/${sample}_2_trimmed.fq"
    SINGLE="$TRIM_DIR/${sample}_trimmed_single.fq"

    [[ -f "$RV" ]] || { echo "WARNING: No reverse read for $sample. Skipping."; continue; }

    echo "════════════════════════════════════"
    echo "Processing: $sample"

    # Step 1 — Human paired
    echo "  [1/4] Human removal (paired)..."
    bowtie2 -x "$HG38_INDEX" -1 "$FW" -2 "$RV" \
        --un-conc "$HUMAN_DIR/${sample}_nonhuman_paired.fastq" \
        -p 6 --very-sensitive --no-unal -S /dev/null \
        2> "$STATS_DIR/${sample}_human_paired.txt"

    # Step 2 — Bovine paired
    echo "  [2/4] Bovine removal (paired)..."
    bowtie2 -x "$BOVINE_INDEX" \
        -1 "$HUMAN_DIR/${sample}_nonhuman_paired.1.fastq" \
        -2 "$HUMAN_DIR/${sample}_nonhuman_paired.2.fastq" \
        --un-conc "$BOVINE_DIR/${sample}_final_clean_paired.fastq" \
        -p 6 --very-sensitive --no-unal -S /dev/null \
        2> "$STATS_DIR/${sample}_bovine_paired.txt"

    # Steps 3 & 4 — Singles (only if file exists)
    if [[ -f "$SINGLE" ]]; then

        echo "  [3/4] Human removal (single)..."
        bowtie2 -x "$HG38_INDEX" -U "$SINGLE" \
            --un "$HUMAN_DIR/${sample}_nonhuman_single.fastq" \
            -p 6 --very-sensitive --no-unal -S /dev/null \
            2> "$STATS_DIR/${sample}_human_single.txt"

        echo "  [4/4] Bovine removal (single)..."
        bowtie2 -x "$BOVINE_INDEX" \
            -U "$HUMAN_DIR/${sample}_nonhuman_single.fastq" \
            --un "$BOVINE_DIR/${sample}_final_clean_single.fastq" \
            -p 6 --very-sensitive --no-unal -S /dev/null \
            2> "$STATS_DIR/${sample}_bovine_single.txt"

    else
        echo "  No single reads for $sample — skipping steps 3 & 4."
    fi

    echo "  Completed: $sample"

done

echo "All host removal completed."
echo "Clean paired reads : $BOVINE_DIR/*_final_clean_paired.1/.2.fastq"
echo "Clean single reads : $BOVINE_DIR/*_final_clean_single.fastq"
