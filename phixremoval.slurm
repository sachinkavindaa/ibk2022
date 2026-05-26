#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=6
#SBATCH --mem=32gb
#SBATCH --partition=guest,batch
#SBATCH --time=168:00:00
#SBATCH --job-name=PhiX_day1
#SBATCH --error=logs/PhiX_%j.err
#SBATCH --output=logs/PhiX_%j.out
#SBATCH --mail-user=echandrasekara2@unl.edu
#SBATCH --mail-type=END,FAIL

set -euo pipefail

module load bbmap/39.06

# ── Paths ──────────────────────────────────────────────────────────────────────
BASE_DIR="/work/samodha/sachin/ShotgunM/Final_Analysis/data/combined_data"
DATA_DIR="$BASE_DIR/day1_25healthy_25infected"
OUT_DIR="$DATA_DIR/Output"
PHIX_DIR="$OUT_DIR/01_PhiX_contamination"
STATS_DIR="$PHIX_DIR/stats"

mkdir -p "$PHIX_DIR" "$STATS_DIR" logs

cd "$DATA_DIR"

echo "Starting PhiX removal..."
echo "Data directory: $DATA_DIR"

# ── Loop over forward reads only ───────────────────────────────────────────────
for FW in *_1.fq.gz *_R1.fq.gz; do

    [[ -e "$FW" ]] || continue

    # Strip suffix to get sample name
    sample="${FW%_1.fq.gz}"
    sample="${sample%_R1.fq.gz}"

    # Find matching reverse read
    if [[ -f "${sample}_2.fq.gz" ]]; then
        RV="${sample}_2.fq.gz"
    elif [[ -f "${sample}_R2.fq.gz" ]]; then
        RV="${sample}_R2.fq.gz"
    else
        echo "WARNING: Reverse read not found for $sample. Skipping."
        continue
    fi

    echo "Processing sample: $sample"
    echo "  Forward : $FW"
    echo "  Reverse : $RV"

    bbduk.sh \
        in1="$FW" \
        in2="$RV" \
        out1="$PHIX_DIR/${sample}_1_phiX_clean.fq.gz" \
        out2="$PHIX_DIR/${sample}_2_phiX_clean.fq.gz" \
        ref=artifacts,phix \
        k=31 \
        hdist=1 \
        threads=6 \
        stats="$STATS_DIR/${sample}_phiX_stats.txt"

    echo "Completed: $sample"
    echo "────────────────────────────────────"

done

echo "All PhiX removal jobs completed."
