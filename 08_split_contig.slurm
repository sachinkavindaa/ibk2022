#!/bin/bash -l
#SBATCH --job-name=split_contigs
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --time=4:00:00
#SBATCH --mem=32gb
#SBATCH --output=logs/step5_split_%j.out
#SBATCH --error=logs/step5_split_%j.err
#SBATCH --partition=batch,guest
#SBATCH --mail-user=echandrasekara2@unl.edu
#SBATCH --mail-type=END,FAIL

set -euo pipefail

module load python/3.9

# =========================
# PATHS
# =========================
WORKDIR="/work/samodha/sachin/ShotgunM/Final_Analysis/data/combined_data/day1_25healthy_25infected/Output"
ASSEMBLY="${WORKDIR}/06_Merged/all_contigs.fa"
KRAKEN_OUT="${WORKDIR}/08_Kraken2_classification/kraken2_output.txt"
KRAKEN_REPORT="${WORKDIR}/08_Kraken2_classification/kraken2_report.txt"
OUTDIR="${WORKDIR}/09_Domain_split"

mkdir -p "${OUTDIR}" logs

echo "======================================"
echo "STEP 5: Split contigs by domain"
echo "Assembly      : ${ASSEMBLY}"
echo "Kraken2 out   : ${KRAKEN_OUT}"
echo "Kraken2 report: ${KRAKEN_REPORT}"
echo "Output dir    : ${OUTDIR}"
echo "======================================"

PYSCRIPT="${OUTDIR}/split_by_domain.py"

cat > "$PYSCRIPT" << 'PYEOF'
import sys, re
from collections import Counter

kraken_out    = sys.argv[1]
kraken_report = sys.argv[2]
assembly      = sys.argv[3]
outdir        = sys.argv[4]

# -------------------------------------------------------
# Step A: Build taxid → domain from Kraken2 report
#
# Key insight from report inspection:
#   Rank D  + name "Eukaryota" → eukaryote
#   Rank D  + name "Bacteria"  → bacteria
#   Rank R1 + name "Viruses"   → virus  ← different rank!
#
# Strategy: track current_domain as we walk the report.
# Reset current_domain whenever we hit a top-level clade.
# Use indentation depth to detect when we leave a clade.
# Simpler approach: just check rank+name for known domains,
# and propagate until we hit another domain-level entry.
# -------------------------------------------------------

# These rank+name combinations signal a new domain
DOMAIN_TRIGGERS = {
    ("D",  "eukaryota") : "eukaryote",
    ("D",  "bacteria")  : "bacteria",
    ("D",  "archaea")   : "unclassified",
    ("R1", "viruses")   : "virus",       # Viruses sits at R1 in NCBI taxonomy
}

print("Building taxid -> domain map from report...")
taxid_to_domain = {}
current_domain  = "unclassified"

with open(kraken_report) as fh:
    for line in fh:
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 6:
            continue
        rank_code = parts[3].strip()
        taxid     = parts[4].strip()
        name      = parts[5].strip().lower()

        # Check if this line is a domain trigger
        key = (rank_code, name)
        if key in DOMAIN_TRIGGERS:
            current_domain = DOMAIN_TRIGGERS[key]

        taxid_to_domain[taxid] = current_domain

print("  Taxids mapped: " + str(len(taxid_to_domain)))

# Verify domains were found
domain_check = Counter(taxid_to_domain.values())
print("  Domain coverage in report:")
for d in ["bacteria", "virus", "eukaryote", "unclassified"]:
    print("    " + d.ljust(15) + ": " + str(domain_check.get(d, 0)) + " taxids")

# -------------------------------------------------------
# Step B: Parse Kraken2 output → contig_id : domain
# Extract taxid from column 3 with regex
# Use full contig ID as-is (ibk236_k141_122)
# -------------------------------------------------------
print("\nParsing Kraken2 output...")
contig_domain = {}

with open(kraken_out) as fh:
    for line in fh:
        parts     = line.strip().split("\t")
        status    = parts[0]
        contig_id = parts[1]

        if status == "U":
            contig_domain[contig_id] = "unclassified"
        else:
            m      = re.search(r'\(taxid (\d+)\)', parts[2])
            taxid  = m.group(1) if m else "0"
            domain = taxid_to_domain.get(taxid, "unclassified")
            contig_domain[contig_id] = domain

counts = Counter(contig_domain.values())
print("Domain counts:")
print("  Bacteria     : " + str(counts.get("bacteria", 0)))
print("  Virus        : " + str(counts.get("virus", 0)))
print("  Eukaryote    : " + str(counts.get("eukaryote", 0)))
print("  Unclassified : " + str(counts.get("unclassified", 0)))
print("  Total        : " + str(len(contig_domain)))

# -------------------------------------------------------
# Step C: Write domain-split FASTA using full contig IDs
# -------------------------------------------------------
out_handles = {
    "bacteria"     : open(outdir + "/bacteria_contigs.fa",     "w"),
    "virus"        : open(outdir + "/virus_contigs.fa",        "w"),
    "eukaryote"    : open(outdir + "/eukaryote_contigs.fa",    "w"),
    "unclassified" : open(outdir + "/unclassified_contigs.fa", "w"),
}

written        = Counter()
current_domain = "unclassified"
not_found      = 0

print("\nWriting FASTA files...")
with open(assembly) as fh:
    for line in fh:
        if line.startswith(">"):
            contig_id = line[1:].split()[0]
            if contig_id in contig_domain:
                current_domain = contig_domain[contig_id]
            else:
                current_domain = "unclassified"
                not_found += 1
            out_handles[current_domain].write(line)
            written[current_domain] += 1
        else:
            out_handles[current_domain].write(line)

for fh in out_handles.values():
    fh.close()

print("\nContigs written:")
for domain in ["bacteria", "virus", "eukaryote", "unclassified"]:
    print("  " + domain.ljust(15) + ": " + str(written.get(domain, 0)) + " contigs")

if not_found > 0:
    print("\nNote: " + str(not_found) + " contigs not in Kraken2 output (sent to unclassified)")

print("\nDone.")
PYEOF

python3 "$PYSCRIPT" "$KRAKEN_OUT" "$KRAKEN_REPORT" "$ASSEMBLY" "$OUTDIR"

echo ""
echo "--- Output file sizes ---"
for f in "${OUTDIR}"/*.fa; do
    COUNT=$(grep -c ">" "$f" 2>/dev/null || echo 0)
    SIZE=$(du -sh "$f" | cut -f1)
    echo "  $(basename $f): ${COUNT} contigs, ${SIZE}"
done

echo ""
echo "======================================"
echo "STEP 5 COMPLETE"
echo "======================================"
