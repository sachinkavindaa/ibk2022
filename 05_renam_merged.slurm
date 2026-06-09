#!/bin/bash
#SBATCH --job-name=merge_contigs
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=2:00:00
#SBATCH --mem=16gb
#SBATCH --output=logs/merge.%J.out
#SBATCH --error=logs/merge.%J.err
#SBATCH --partition=yinlab,batch,guest

set -euo pipefail
shopt -s nullglob

ASSEM="/work/samodha/sachin/ShotgunM/Final_Analysis/data/combined_data/day1_25healthy_25infected/Output/05_Assembly_per_sample"

OUTDIR="/work/samodha/sachin/ShotgunM/Final_Analysis/data/combined_data/day1_25healthy_25infected/Output/06_Merged"

mkdir -p "${OUTDIR}" logs

ALL_CONTIGS="${OUTDIR}/all_contigs.fa"
> "${ALL_CONTIGS}"

COUNT=0

for FINAL in "${ASSEM}"/*/final_contigs.fa; do

    SAMPLE=$(basename "$(dirname "$FINAL")")

    echo "Processing: ${SAMPLE}"

    awk -v s="${SAMPLE}" '
        /^>/ { print ">" s "_" substr($0,2); next }
        { print }
    ' "${FINAL}" >> "${ALL_CONTIGS}"

    COUNT=$((COUNT+1))

done

TOTAL=$(grep -c "^>" "${ALL_CONTIGS}")

echo "Merged samples: ${COUNT}"
echo "Total contigs: ${TOTAL}"
echo "Output file: ${ALL_CONTIGS}"
