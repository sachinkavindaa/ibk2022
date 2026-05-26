#!/bin/bash -l
#SBATCH --job-name=make_metaeuk_db
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=24:00:00
#SBATCH --mem=120gb
#SBATCH --output=logs/make_metaeuk_db_%j.out
#SBATCH --error=logs/make_metaeuk_db_%j.err
#SBATCH --partition=batch,guest

set -euo pipefail

module load metaeuk
module load biodata

mkdir -p logs
mkdir -p /work/samodha/sachin/databases/metaeuk_uniref50

metaeuk createdb \
  "$UNIPROT_UNIREF50" \
  /work/samodha/sachin/databases/metaeuk_uniref50/uniref50

echo "MetaEuk database created:"
ls -lh /work/samodha/sachin/databases/metaeuk_uniref50/
