#!/bin/bash -l
#SBATCH --job-name=kraken2_classify
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=24:00:00
#SBATCH --mem=150gb
#SBATCH --output=logs/step4_kraken2_new_%j.out
#SBATCH --error=logs/step4_kraken2__new_%j.err
#SBATCH --partition=batch,guest
#SBATCH --mail-user=echandrasekara2@unl.edu
#SBATCH --mail-type=END,FAIL

set -euo pipefail

module load kraken2

# =========================
# PATHS
# =========================
WORKDIR="/work/samodha/sachin/ShotgunM/Final_Analysis/data/combined_data/day1_25healthy_25infected/Output"
ASSEMBLY="${WORKDIR}/06_Merged/all_contigs.fa"
OUTDIR="${WORKDIR}/08_Kraken2_classification_new"
KRAKEN2_DB="${KRAKEN2_PLUSPF_DB}"
THREADS=${SLURM_NTASKS_PER_NODE:-32}

mkdir -p "${OUTDIR}" logs

# Confirm we have the right contig count before running
CONTIG_COUNT=$(grep -c ">" "$ASSEMBLY")
echo "======================================"
echo "STEP 4: Kraken2 classification (rerun)"
echo "Input    : ${ASSEMBLY}"
echo "Contigs  : ${CONTIG_COUNT}"
echo "Database : ${KRAKEN2_DB}"
echo "Threads  : ${THREADS}"
echo "======================================"

if [[ "$CONTIG_COUNT" -lt 30000 ]]; then
    echo "ERROR: Only ${CONTIG_COUNT} contigs found — expected ~33719."
    echo "Check that all samples are merged into ${ASSEMBLY}"
    exit 1
fi

# Remove old incomplete output so we start fresh
rm -f "${OUTDIR}/kraken2_output.txt" "${OUTDIR}/kraken2_report.txt"

kraken2 \
    --db "${KRAKEN2_DB}" \
    --threads "${THREADS}" \
    --output "${OUTDIR}/kraken2_output.txt" \
    --report "${OUTDIR}/kraken2_report.txt" \
    --use-names \
    --confidence 0.2 \
    "${ASSEMBLY}"

# Summary
TOTAL=$(wc -l < "${OUTDIR}/kraken2_output.txt")
CLASSIFIED=$(awk '$1=="C"' "${OUTDIR}/kraken2_output.txt" | wc -l)
UNCLASSIFIED=$(awk '$1=="U"' "${OUTDIR}/kraken2_output.txt" | wc -l)

echo ""
echo "--- Classification Summary ---"
echo "Total contigs     : ${TOTAL}"
echo "Classified        : ${CLASSIFIED}"
echo "Unclassified      : ${UNCLASSIFIED}"
echo ""
echo "Top domain breakdown:"
grep -E "Bacteria|Viruses|Eukaryota|Fungi|Protozoa" "${OUTDIR}/kraken2_report.txt" | \
    awk '{printf "  %-30s %s%%\n", $6, $1}' | head -10

echo "======================================"
echo "STEP 4 RERUN COMPLETE"
echo "======================================"
