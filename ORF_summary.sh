#!/bin/bash

BASE="/work/samodha/sachin/ShotgunM/Final_Analysis/data/combined_data/day1_25healthy_25infected/Output"

SPLIT_DIR="${BASE}/09_Domain_split"
ORF_DIR="${BASE}/10_ORF_prediction"

echo -e "Domain\tContigs\tPredicted_ORFs\tORFs_per_Contig" > ORF_summary.tsv

############################
# BACTERIA
############################
bac_contigs=$(grep -c "^>" "${SPLIT_DIR}/bacteria_contigs.fa" 2>/dev/null || echo 0)
bac_orfs=$(grep -c "^>" "${ORF_DIR}/bacteria/bacteria_orfs.faa" 2>/dev/null || echo 0)

bac_ratio=$(awk -v a="$bac_orfs" -v b="$bac_contigs" \
'BEGIN{if(b>0) printf "%.2f",a/b; else print 0}')

echo -e "Bacteria\t${bac_contigs}\t${bac_orfs}\t${bac_ratio}" >> ORF_summary.tsv

############################
# EUKARYOTE
############################
euk_contigs=$(grep -c "^>" "${SPLIT_DIR}/eukaryote_contigs.fa" 2>/dev/null || echo 0)

euk_file=$(find "${ORF_DIR}/eukaryote" -name "*.fas" | head -1)

if [[ -n "$euk_file" ]]; then
    euk_orfs=$(grep -c "^>" "$euk_file")
else
    euk_orfs=0
fi

euk_ratio=$(awk -v a="$euk_orfs" -v b="$euk_contigs" \
'BEGIN{if(b>0) printf "%.2f",a/b; else print 0}')

echo -e "Eukaryote\t${euk_contigs}\t${euk_orfs}\t${euk_ratio}" >> ORF_summary.tsv

############################
# VIRUS
############################
virus_contigs=$(grep -c "^>" "${SPLIT_DIR}/virus_contigs.fa" 2>/dev/null || echo 0)

virus_file=$(find "${ORF_DIR}/virus" -name "*.faa" | head -1)

if [[ -n "$virus_file" ]]; then
    virus_orfs=$(grep -c "^>" "$virus_file")
else
    virus_orfs=0
fi

virus_ratio=$(awk -v a="$virus_orfs" -v b="$virus_contigs" \
'BEGIN{if(b>0) printf "%.2f",a/b; else print 0}')

echo -e "Virus\t${virus_contigs}\t${virus_orfs}\t${virus_ratio}" >> ORF_summary.tsv

echo
echo "================================="
echo "ORF Summary"
echo "================================="
column -t ORF_summary.tsv
