#!/bin/bash
#SBATCH --job-name=dram_annotate
#SBATCH --partition=guest
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=48:00:00
#SBATCH --output=/work/samodha/sachin/ShotgunM/Final_Analysis/data/combined_data/day1_25healthy_25infected/Output/logs/dram_%j.out
#SBATCH --error=/work/samodha/sachin/ShotgunM/Final_Analysis/data/combined_data/day1_25healthy_25infected/Output/logs/dram_%j.err

set -euo pipefail

INPUT_FAA="/work/samodha/sachin/ShotgunM/Final_Analysis/data/combined_data/day1_25healthy_25infected/Output/11_NR_catalogs/bacteria/bacteria_NR_catalog.faa"
OUTDIR="/work/samodha/sachin/ShotgunM/Final_Analysis/data/combined_data/day1_25healthy_25infected/Output/12_functional_annotation/dram"

mkdir -p "${OUTDIR}"

echo "=== Checking input ==="
if [ ! -s "${INPUT_FAA}" ]; then
    echo "ERROR: ${INPUT_FAA} missing or empty" >&2
    exit 1
fi
echo "Sequences in input: $(grep -c '^>' ${INPUT_FAA})"

module purge
module load dram/1.2

echo "=== Running DRAM annotate_genes ==="
DRAM.py annotate_genes \
    -i "${INPUT_FAA}" \
    -o "${OUTDIR}/annotation" \
    --threads 16

echo "=== Running DRAM distill ==="
DRAM.py distill \
    -i "${OUTDIR}/annotation/annotations.tsv" \
    -o "${OUTDIR}/distillate"

echo "=== DRAM done. Check ${OUTDIR} ==="
