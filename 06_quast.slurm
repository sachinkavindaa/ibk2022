#!/bin/bash
#SBATCH --job-name=quast
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=4:00:00
#SBATCH --mem=32gb
#SBATCH --output=logs/quast.%J.out
#SBATCH --error=logs/quast.%J.err
#SBATCH --partition=yinlab,batch,guest

set -euo pipefail

ml quast

CONTIGS=/work/samodha/sachin/ShotgunM/Final_Analysis/data/combined_data/day1_25healthy_25infected/Output/06_Merged/all_contigs.fa
OUTDIR=/work/samodha/sachin/ShotgunM/Final_Analysis/data/combined_data/day1_25healthy_25infected/Output/07_QUAST
mkdir -p "${OUTDIR}" logs

quast.py \
    "${CONTIGS}" \
    --output-dir "${OUTDIR}" \
    --threads 16 \
    --min-contig 500 \
    --plots-format pdf

echo "QUAST complete. Report: ${OUTDIR}/report.html"
