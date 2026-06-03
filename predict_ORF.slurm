#!/bin/bash -l
#SBATCH --job-name=predict_orfs
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=168:00:00
#SBATCH --mem=120gb
#SBATCH --output=logs/step6_orfs_%j.out
#SBATCH --error=logs/step6_orfs_%j.err
#SBATCH --partition=batch,guest
#SBATCH --mail-user=echandrasekara2@unl.edu
#SBATCH --mail-type=END,FAIL

set -euo pipefail

module load prodigal
module load vibrant/1.2
module load metaeuk

WORKDIR="/work/samodha/sachin/ShotgunM/Final_Analysis/data/combined_data/day1_25healthy_25infected/Output"
SPLITDIR="${WORKDIR}/09_Domain_split"
OUTDIR="${WORKDIR}/10_ORF_prediction"
THREADS=${SLURM_NTASKS_PER_NODE:-16}

mkdir -p "${OUTDIR}"/{bacteria,virus,eukaryote} logs

echo "======================================"
echo "STEP 6: ORF Prediction per domain"
echo "Input dir : ${SPLITDIR}"
echo "Output dir: ${OUTDIR}"
echo "Threads   : ${THREADS}"
echo "======================================"

# -------------------------------------------------------
# BACTERIA — Prodigal meta mode
# -------------------------------------------------------
echo ""
echo "--- Bacteria: Prodigal ---"

BAC_IN="${SPLITDIR}/bacteria_contigs.fa"
BAC_OUT="${OUTDIR}/bacteria"

if [[ ! -s "$BAC_IN" ]] || [[ $(grep -c ">" "$BAC_IN") -eq 0 ]]; then
    echo "WARNING: No bacteria contigs found. Skipping bacteria."
else
    prodigal \
        -i "$BAC_IN" \
        -a "${BAC_OUT}/bacteria_orfs.faa" \
        -d "${BAC_OUT}/bacteria_orfs.ffn" \
        -f gff \
        -o "${BAC_OUT}/bacteria_orfs.gff" \
        -p meta \
        -q

    echo "Bacteria ORFs: $(grep -c '>' ${BAC_OUT}/bacteria_orfs.faa)"
fi


# -------------------------------------------------------
# VIRUS — VIBRANT
# -------------------------------------------------------
echo ""
echo "--- Virus: VIBRANT ---"

VIR_IN="${SPLITDIR}/virus_contigs.fa"
VIR_OUT="${OUTDIR}/virus"

if [[ ! -s "$VIR_IN" ]] || [[ $(grep -c ">" "$VIR_IN") -eq 0 ]]; then
    echo "WARNING: No virus contigs found. Skipping virus."
else
    VIBRANT_run.py \
        -i "$VIR_IN" \
        -folder "$VIR_OUT" \
        -t "$THREADS"

    echo "VIBRANT virus prediction finished."
    echo "Check output folder:"
    echo "${VIR_OUT}"
fi


# -------------------------------------------------------
# EUKARYOTE — MetaEuk
# -------------------------------------------------------
echo ""
echo "--- Eukaryote: MetaEuk ---"

EUK_IN="${SPLITDIR}/eukaryote_contigs.fa"
EUK_OUT="${OUTDIR}/eukaryote"
EUK_TMP="${OUTDIR}/eukaryote/tmp"

# IMPORTANT: update this database path
METAEUK_DB="/work/samodha/sachin/databases/metaeuk_uniref50/uniref50"

if [[ ! -s "$EUK_IN" ]] || [[ $(grep -c ">" "$EUK_IN") -eq 0 ]]; then
    echo "WARNING: No eukaryote contigs found. Skipping eukaryote."
else
    mkdir -p "$EUK_TMP"

    metaeuk easy-predict \
        "$EUK_IN" \
        "$METAEUK_DB" \
        "${EUK_OUT}/eukaryote_orfs" \
        "$EUK_TMP" \
        --threads "$THREADS"

    echo "MetaEuk eukaryotic ORF prediction finished."
    echo "Check output prefix:"
    echo "${EUK_OUT}/eukaryote_orfs"
fi


echo ""
echo "======================================"
echo "STEP 6 COMPLETE — ORF summary:"
echo "======================================"

if [[ -f "${OUTDIR}/bacteria/bacteria_orfs.faa" ]]; then
    echo "Bacteria ORFs: $(grep -c '>' ${OUTDIR}/bacteria/bacteria_orfs.faa)"
else
    echo "Bacteria: skipped"
fi

echo "Virus: check VIBRANT output folder:"
echo "${OUTDIR}/virus"

echo "Eukaryote: check MetaEuk output prefix:"
echo "${OUTDIR}/eukaryote/eukaryote_orfs"

echo "======================================"
