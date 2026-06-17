#!/bin/bash
#SBATCH --job-name=metaphlan4_taxonomy
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=24:00:00
#SBATCH --mem=80gb
#SBATCH --output=logs/metaphlan4_%j.out
#SBATCH --error=logs/metaphlan4_%j.err
#SBATCH --partition=batch

mkdir -p logs

BASE="/work/samodha/sachin/ShotgunM/Final_Analysis/data/combined_data/day1_25healthy_25infected/Output"
CONTIGS="${BASE}/06_Merged/all_contigs.fa"
OUT="${BASE}/07_Taxonomy"

mkdir -p "$OUT"

ml metaphlan/4.1

echo "Starting MetaPhlAn4 at $(date)"
echo "Input: $CONTIGS"
echo "Output: $OUT"

if [ ! -f "$CONTIGS" ]; then
    echo "ERROR: Input file not found"
    exit 1
fi

metaphlan "$CONTIGS" \
    --input_type fasta \
    --nproc 16 \
    --output_file "${OUT}/metaphlan4_profile.txt"

echo ""
echo "MetaPhlAn4 completed at $(date)"

if [ -f "${OUT}/metaphlan4_profile.txt" ]; then
    echo "Profile created successfully"
    
    grep "s__" "${OUT}/metaphlan4_profile.txt" > "${OUT}/species_abundance.txt"
    grep "k__Bacteria" "${OUT}/metaphlan4_profile.txt" > "${OUT}/bacteria_profile.txt"
    grep "k__Archaea" "${OUT}/metaphlan4_profile.txt" > "${OUT}/archaea_profile.txt"
    grep "k__Viruses" "${OUT}/metaphlan4_profile.txt" > "${OUT}/viruses_profile.txt"
    
    echo ""
    echo "Summary:"
    echo "Total lines: $(wc -l < ${OUT}/metaphlan4_profile.txt)"
    echo "Bacteria: $(wc -l < ${OUT}/bacteria_profile.txt)"
    echo "Archaea: $(wc -l < ${OUT}/archaea_profile.txt)"
    echo "Viruses: $(wc -l < ${OUT}/viruses_profile.txt)"
    echo "Species: $(wc -l < ${OUT}/species_abundance.txt)"
else
    echo "ERROR: Profile not created"
    exit 1
fi

echo ""
echo "Done at $(date)"
