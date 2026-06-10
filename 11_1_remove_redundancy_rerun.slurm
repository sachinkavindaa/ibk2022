#!/bin/bash -l
#SBATCH --job-name=mmseqs2_cluster
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=48:00:00
#SBATCH --mem=200gb
#SBATCH --output=logs/step7_mmseqs_%j.out
#SBATCH --error=logs/step7_mmseqs_%j.err
#SBATCH --partition=batch,guest
#SBATCH --mail-user=echandrasekara2@unl.edu
#SBATCH --mail-type=END,FAIL

set -euo pipefail

module load mmseqs2

# =========================
# PATHS
# =========================
WORKDIR="/work/samodha/sachin/ShotgunM/Final_Analysis/data/combined_data/day1_25healthy_25infected/Output"
ORFDIR="${WORKDIR}/10_ORF_prediction"
OUTDIR="${WORKDIR}/11_NR_catalogs"
MMSEQS_TMP="${WORKDIR}/tmp_mmseqs"
THREADS=${SLURM_NTASKS_PER_NODE:-32}

# Clustering parameters (standard for gene catalogs)
IDENTITY=0.95        # 95% sequence identity
COVERAGE=0.90        # 90% coverage of shorter sequence
COV_MODE=2           # bidirectional coverage (correct for cluster-mode 0)

mkdir -p "${OUTDIR}" "${MMSEQS_TMP}" logs

echo "======================================"
echo "STEP 7: Redundancy removal — virus + eukaryote RERUN"
echo "NOTE: Using NUCLEOTIDE files for BWA compatibility"
echo "  Virus     : virus_orfs.ffn       (Prodigal nucleotide)"
echo "  Eukaryote : eukaryote_orfs.codon.fas (MetaEuk codon/nucleotide)"
echo "  Bacteria  : SKIPPED (already correct)"
echo "Identity  : ${IDENTITY}"
echo "Coverage  : ${COVERAGE}"
echo "Cov-mode  : ${COV_MODE} (bidirectional)"
echo "Threads   : ${THREADS}"
echo "======================================"

# Clean up old wrong protein-based catalogs before rerunning
echo "Removing old protein-based catalogs for virus and eukaryote..."
rm -rf "${OUTDIR}/virus"
rm -rf "${OUTDIR}/eukaryote"
rm -rf "${MMSEQS_TMP}/virus"
rm -rf "${MMSEQS_TMP}/eukaryote"
echo "Done. Rebuilding from nucleotide ORFs..."

# -------------------------------------------------------
# Function to run MMseqs2 clustering for one domain
# -------------------------------------------------------
run_mmseqs_cluster() {
    local DOMAIN=$1
    local INPUT_FA=$2
    local DOMAIN_OUT="${OUTDIR}/${DOMAIN}"
    local DOMAIN_TMP="${MMSEQS_TMP}/${DOMAIN}"

    mkdir -p "${DOMAIN_OUT}" "${DOMAIN_TMP}"

    if [[ ! -f "$INPUT_FA" ]]; then
        echo "ERROR: Input file not found: ${INPUT_FA}"
        echo "Please verify the file path and rerun."
        exit 1
    fi

    SEQ_COUNT=$(grep -c ">" "$INPUT_FA")
    if [[ "$SEQ_COUNT" -eq 0 ]]; then
        echo "WARNING: No sequences found in ${INPUT_FA}. Skipping ${DOMAIN}."
        return
    fi

    echo ""
    echo "--- ${DOMAIN}: ${SEQ_COUNT} input sequences ---"
    echo "    Input: ${INPUT_FA}"

    # Verify sequences are nucleotide (first non-header line should have ATCG)
    FIRST_SEQ=$(grep -v ">" "$INPUT_FA" | head -1)
    if echo "$FIRST_SEQ" | grep -qiP "^[ATCGN]+$"; then
        echo "    Sequence type: NUCLEOTIDE ✓"
    else
        echo "ERROR: ${INPUT_FA} does not appear to be nucleotide!"
        echo "First sequence: ${FIRST_SEQ:0:60}"
        echo "BWA requires nucleotide sequences. Please check your input file."
        exit 1
    fi

    # Step 1 — create MMseqs2 database
    # NOTE: createdb does NOT accept --threads
    echo "  [1/4] Creating MMseqs2 database..."
    mmseqs createdb \
        "$INPUT_FA" \
        "${DOMAIN_TMP}/${DOMAIN}_db"

    # Step 2 — cluster at defined identity + coverage
    echo "  [2/4] Clustering at ${IDENTITY} identity / ${COVERAGE} coverage..."
    mmseqs cluster \
        "${DOMAIN_TMP}/${DOMAIN}_db" \
        "${DOMAIN_TMP}/${DOMAIN}_cluster_db" \
        "${DOMAIN_TMP}/tmp" \
        --min-seq-id "${IDENTITY}" \
        -c "${COVERAGE}" \
        --cov-mode "${COV_MODE}" \
        --cluster-mode 0 \
        --threads "${THREADS}"

    # Step 3 — extract representative sequences (NR catalog)
    echo "  [3/4] Extracting representative sequences..."
    mmseqs createsubdb \
        "${DOMAIN_TMP}/${DOMAIN}_cluster_db" \
        "${DOMAIN_TMP}/${DOMAIN}_db" \
        "${DOMAIN_TMP}/${DOMAIN}_rep_db"

    mmseqs convert2fasta \
        "${DOMAIN_TMP}/${DOMAIN}_rep_db" \
        "${DOMAIN_OUT}/${DOMAIN}_NR_catalog.fa"

    # Step 4 — write cluster membership TSV
    echo "  [4/4] Writing cluster membership TSV..."
    mmseqs createtsv \
        "${DOMAIN_TMP}/${DOMAIN}_db" \
        "${DOMAIN_TMP}/${DOMAIN}_db" \
        "${DOMAIN_TMP}/${DOMAIN}_cluster_db" \
        "${DOMAIN_OUT}/${DOMAIN}_cluster_membership.tsv" \
        --threads "${THREADS}"

    NR_COUNT=$(grep -c ">" "${DOMAIN_OUT}/${DOMAIN}_NR_catalog.fa")
    echo ""
    echo "  ${DOMAIN} NR catalog: ${NR_COUNT} representative genes"
    echo "  Catalog saved : ${DOMAIN_OUT}/${DOMAIN}_NR_catalog.fa"
    echo "  Membership TSV: ${DOMAIN_OUT}/${DOMAIN}_cluster_membership.tsv"

    # Quick sanity check — confirm output is nucleotide
    FIRST_OUT=$(grep -v ">" "${DOMAIN_OUT}/${DOMAIN}_NR_catalog.fa" | head -1)
    if echo "$FIRST_OUT" | grep -qiP "^[ATCGN]+$"; then
        echo "  Output check  : NUCLEOTIDE ✓ — ready for BWA mapping"
    else
        echo "  WARNING: Output does not look like nucleotide. Check catalog."
    fi
}

# -------------------------------------------------------
# BACTERIA — already correct, skip
# -------------------------------------------------------
echo ""
echo "--- bacteria: SKIPPED (already built correctly from bacteria_orfs.faa) ---"
echo "    Existing catalog: ${OUTDIR}/bacteria/bacteria_NR_catalog.fa"
BAC_CAT="${OUTDIR}/bacteria/bacteria_NR_catalog.fa"
[[ -f "$BAC_CAT" ]] && echo "    Bacteria NR genes: $(grep -c '>' $BAC_CAT)" || echo "    WARNING: bacteria catalog not found!"

# -------------------------------------------------------
# VIRUS — use .ffn (nucleotide) not .faa (protein)
# -------------------------------------------------------
run_mmseqs_cluster "virus" "${ORFDIR}/virus/virus_orfs.ffn"

# -------------------------------------------------------
# EUKARYOTE — use .codon.fas (nucleotide CDS) not .fas (protein)
# -------------------------------------------------------
run_mmseqs_cluster "eukaryote" "${ORFDIR}/eukaryote/eukaryote_orfs.codon.fas"

# -------------------------------------------------------
# SUMMARY
# -------------------------------------------------------
echo ""
echo "======================================"
echo "STEP 7 RERUN COMPLETE — NR Catalog summary:"
for domain in bacteria virus eukaryote; do
    CAT="${OUTDIR}/${domain}/${domain}_NR_catalog.fa"
    if [[ -f "$CAT" ]]; then
        COUNT=$(grep -c ">" "$CAT")
        echo "  ${domain}: ${COUNT} NR genes → ${CAT}"
    else
        echo "  ${domain}: NOT FOUND — check logs"
    fi
done
echo ""
echo "Next step: delete old virus/eukaryote BAMs and rerun step9_map_reads.sh"
echo "======================================"
