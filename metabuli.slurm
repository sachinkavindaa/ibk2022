#!/bin/bash -l
#SBATCH --job-name=metabuli_taxonomy
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=24:00:00
#SBATCH --mem=120gb
#SBATCH --output=logs/metabuli_%j.out
#SBATCH --error=logs/metabuli_%j.err
#SBATCH --partition=batch

set -euo pipefail

mkdir -p logs

BASE="/work/samodha/sachin/ShotgunM/Final_Analysis/data/combined_data/day1_25healthy_25infected/Output"

CONTIGS="${BASE}/06_Merged/all_contigs.fa"
DB="/work/samodha/sachin/databases/metabuli/refseq_standard"
OUT="${BASE}/07_Metabuli_taxonomy"
JOBID="all_contigs"

mkdir -p "$OUT"

module load metabuli/1.2

echo "Starting Metabuli at $(date)"
echo "Contigs: $CONTIGS"
echo "Database: $DB"
echo "Output: $OUT"

if [ ! -f "$CONTIGS" ]; then
    echo "ERROR: contig file not found: $CONTIGS"
    exit 1
fi

if [ ! -d "$DB" ]; then
    echo "ERROR: database folder not found: $DB"
    exit 1
fi

metabuli classify \
    "$CONTIGS" \
    "$DB" \
    "$OUT" \
    "$JOBID" \
    --threads 16 \
    --seq-mode 1 \
    --min-score 0.5 \
    --lineage 1

echo "Metabuli finished at $(date)"

CLASS="${OUT}/${JOBID}_classifications.tsv"
REPORT="${OUT}/${JOBID}_report.tsv"

if [ ! -f "$CLASS" ]; then
    echo "ERROR: classification file not created"
    ls -lh "$OUT"
    exit 1
fi

echo "Creating domain summaries..."

awk -F'\t' '$0 ~ /Bacteria|d__Bacteria/' "$CLASS" > "${OUT}/bacteria_classified.tsv"
awk -F'\t' '$0 ~ /Archaea|d__Archaea/' "$CLASS" > "${OUT}/archaea_classified.tsv"
awk -F'\t' '$0 ~ /Viruses|Virus|d__Viruses/' "$CLASS" > "${OUT}/viruses_classified.tsv"
awk -F'\t' '$0 ~ /Eukaryota|Eukarya|d__Eukaryota/' "$CLASS" > "${OUT}/eukaryota_classified.tsv"

echo ""
echo "====== Metabuli Classification Summary ======"
echo "Total lines: $(wc -l < "$CLASS")"
echo "Bacteria: $(wc -l < "${OUT}/bacteria_classified.tsv")"
echo "Archaea: $(wc -l < "${OUT}/archaea_classified.tsv")"
echo "Viruses: $(wc -l < "${OUT}/viruses_classified.tsv")"
echo "Eukaryota: $(wc -l < "${OUT}/eukaryota_classified.tsv")"
echo "============================================="
echo ""
echo "Classification file: $CLASS"
echo "Report file: $REPORT"
