#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=6
#SBATCH --mem=32gb
#SBATCH --partition=guest,batch
#SBATCH --time=168:00:00
#SBATCH --job-name=Trim_day1
#SBATCH --error=logs/Trim_%j.err
#SBATCH --output=logs/Trim_%j.out
#SBATCH --mail-user=echandrasekara2@unl.edu
#SBATCH --mail-type=END,FAIL

set -euo pipefail

module load sickle

# ── Paths ──────────────────────────────────────────────────────────────────────
BASE_DIR="/work/samodha/sachin/ShotgunM/Final_Analysis/data/combined_data"
DATA_DIR="$BASE_DIR/day1_25healthy_25infected"
OUT_DIR="$DATA_DIR/Output"
PHIX_DIR="$OUT_DIR/01_PhiX_contamination"
TRIM_DIR="$OUT_DIR/02_Trimming"
STATS_DIR="$TRIM_DIR/stats"

mkdir -p "$TRIM_DIR" "$STATS_DIR" logs

echo "Starting trimming..."
echo "Reading PhiX-cleaned files from: $PHIX_DIR"

# ── Loop over cleaned forward reads ───────────────────────────────────────────
for FW_GZ in "$PHIX_DIR"/*_1_phiX_clean.fq.gz; do

    [[ -e "$FW_GZ" ]] || continue

    filename=$(basename "$FW_GZ")
    sample="${filename%_1_phiX_clean.fq.gz}"

    RV_GZ="$PHIX_DIR/${sample}_2_phiX_clean.fq.gz"

    if [[ ! -f "$RV_GZ" ]]; then
        echo "WARNING: Reverse read not found for $sample. Skipping."
        continue
    fi

    echo "Processing: $sample"

    # Decompress temporarily — sickle cannot read .fq.gz
    FW_TMP="$TRIM_DIR/${sample}_1_tmp.fq"
    RV_TMP="$TRIM_DIR/${sample}_2_tmp.fq"

    zcat "$FW_GZ" > "$FW_TMP"
    zcat "$RV_GZ" > "$RV_TMP"

    # Run sickle
    sickle pe -t sanger \
        -f "$FW_TMP" \
        -r "$RV_TMP" \
        -o "$TRIM_DIR/${sample}_1_trimmed.fq" \
        -p "$TRIM_DIR/${sample}_2_trimmed.fq" \
        -s "$TRIM_DIR/${sample}_trimmed_single.fq" \
        -q 30 \
        -l 50 \
        2> "$STATS_DIR/${sample}_trim_stats.txt"

    # Clean up temp files immediately to save disk space
    rm "$FW_TMP" "$RV_TMP"

    echo "  Done: $sample"
    echo "────────────────────────────────────"

done

echo "All trimming completed."
echo "Output in: $TRIM_DIR"
