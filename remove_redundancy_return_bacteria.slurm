#!/bin/bash -l
#SBATCH --job-name=bacteria_nr
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=24:00:00
#SBATCH --mem=100gb
#SBATCH --partition=batch,guest
#SBATCH --output=logs/bacteria_nr_%j.out
#SBATCH --error=logs/bacteria_nr_%j.err

set -euo pipefail

module load mmseqs2

WORKDIR="/work/samodha/sachin/ShotgunM/Final_Analysis/data/combined_data/day1_25healthy_25infected/Output"
ORFDIR="${WORKDIR}/10_ORF_prediction/bacteria"
OUTDIR="${WORKDIR}/11_NR_catalogs/bacteria"
TMPDIR="${WORKDIR}/tmp_mmseqs/bacteria"

THREADS=32
IDENTITY=0.95
COVERAGE=0.90

# Clean old bacteria catalog
rm -rf "$OUTDIR"
rm -rf "$TMPDIR"

mkdir -p "$OUTDIR"
mkdir -p "$TMPDIR"

INPUT="${ORFDIR}/bacteria_orfs.ffn"

echo "Input sequences:"
grep -c "^>" "$INPUT"

# Create database
mmseqs createdb \
    "$INPUT" \
    "${TMPDIR}/bacteria_db"

# Cluster
mmseqs cluster \
    "${TMPDIR}/bacteria_db" \
    "${TMPDIR}/bacteria_cluster_db" \
    "${TMPDIR}/tmp" \
    --min-seq-id ${IDENTITY} \
    -c ${COVERAGE} \
    --cov-mode 2 \
    --cluster-mode 0 \
    --threads ${THREADS}

# Extract representatives
mmseqs createsubdb \
    "${TMPDIR}/bacteria_cluster_db" \
    "${TMPDIR}/bacteria_db" \
    "${TMPDIR}/bacteria_rep_db"

mmseqs convert2fasta \
    "${TMPDIR}/bacteria_rep_db" \
    "${OUTDIR}/bacteria_NR_catalog.fa"

echo "NR genes:"
grep -c "^>" "${OUTDIR}/bacteria_NR_catalog.fa"
