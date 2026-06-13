#!/bin/bash -l
#SBATCH --job-name=mapped_reads_summary
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64gb
#SBATCH --time=24:00:00
#SBATCH --partition=guest,batch
#SBATCH --output=mapped_reads_summary_%j.out
#SBATCH --error=mapped_reads_summary_%j.err

set -euo pipefail

module load samtools

BASE="/work/samodha/sachin/ShotgunM/Final_Analysis/data/combined_data/day1_25healthy_25infected/Output/12_Read_mapping"

OUT="$BASE/mapped_reads_summary.tsv"

echo -e "Domain\tBAM_files\tUnique_Mapped_Reads" > "$OUT"

for domain in bacteria2 eukaryote virus
do
    DIR="$BASE/$domain"

    bam_count=$(find "$DIR" -name "*.bam" | wc -l)

    echo "Processing $domain ..."

    unique_reads=$(
        find "$DIR" -name "*.bam" \
        -exec samtools view -F 4 {} \; \
        | cut -f1 \
        | sort -u \
        | wc -l
    )

    echo -e "${domain}\t${bam_count}\t${unique_reads}" >> "$OUT"
done

echo
echo "================================"
echo "Mapping Summary"
echo "================================"
column -t "$OUT"
