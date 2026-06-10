#!/bin/bash -l
#SBATCH --job-name=euk_mmseqs
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=24:00:00
#SBATCH --mem=200gb
#SBATCH --output=logs/euk_mmseqs_%j.out
#SBATCH --error=logs/euk_mmseqs_%j.err
#SBATCH --partition=batch,guest

module load mmseqs2

WORKDIR="/work/samodha/sachin/ShotgunM/Final_Analysis/data/combined_data/day1_25healthy_25infected/Output"

INPUT="${WORKDIR}/10_ORF_prediction/eukaryote/eukaryote_orfs.codon.fas"

OUTDIR="${WORKDIR}/11_NR_catalogs/eukaryote"
TMP="${WORKDIR}/tmp_mmseqs/eukaryote"

mkdir -p "$OUTDIR" "$TMP"

mmseqs createdb \
    "$INPUT" \
    "${TMP}/euk_db"

mmseqs cluster \
    "${TMP}/euk_db" \
    "${TMP}/euk_cluster_db" \
    "${TMP}/tmp" \
    --min-seq-id 0.95 \
    -c 0.90 \
    --cov-mode 2 \
    --cluster-mode 0 \
    --threads 32

mmseqs createsubdb \
    "${TMP}/euk_cluster_db" \
    "${TMP}/euk_db" \
    "${TMP}/euk_rep_db"

mmseqs convert2fasta \
    "${TMP}/euk_rep_db" \
    "${OUTDIR}/eukaryote_NR_catalog.fa"

mmseqs createtsv \
    "${TMP}/euk_db" \
    "${TMP}/euk_db" \
    "${TMP}/euk_cluster_db" \
    "${OUTDIR}/eukaryote_cluster_membership.tsv" \
    --threads 32

echo "DONE"
