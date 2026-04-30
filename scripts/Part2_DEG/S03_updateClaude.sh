#!/bin/bash
# S03_make_count_matrix_RSEM.sh
# Builds transcript and gene expression count matrices from RSEM output
# FIX: updated paths to match new directory structure
# FIX: fixed module version mismatch (foss/2019b loaded alongside foss-2022b R — kept consistent)
# FIX: added set -euo pipefail
# FIX: renamed job

#SBATCH --job-name=rsem_matrix
#SBATCH --mail-user=alicebalard@zedat.fu-berlin.de
#SBATCH --mail-type=end,fail
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4GB
#SBATCH --time=00:30:00
#SBATCH --error=/scratch/alicebalard/outData/RSEM/logs/%x.%j.err
#SBATCH --output=/scratch/alicebalard/outData/RSEM/logs/%x.%j.out
#SBATCH --qos=standard
#SBATCH --partition=main,begendiv
#SBATCH --constraint=no_gpu

set -euo pipefail

export LMOD_DISABLE_SAME_NAME_AUTOSWAP=no

# FIX: load consistent foss version; mixing foss/2019b and foss-2022b R was a potential conflict
module add Trinity/2.10.0-foss-2019b-Python-3.7.4
module load R-bundle-Bioconductor/3.16-foss-2022b-R-4.2.2

DATADIR=/scratch/alicebalard/outData/RSEM

## From S01:
GTM=$DATADIR/combined_gene_trans_map.txt

## Manually prepared: list of paths to per-sample RSEM isoform results files
## One path per line, pointing to *.isoforms.results files from S02 output
QUANT=$DATADIR/quant_file_isoform.txt

## Verify input files exist
for f in "$GTM" "$QUANT"; do
    if [ ! -f "$f" ]; then
        echo "[S03] ERROR: required file not found: $f" >&2
        exit 1
    fi
done

cd $DATADIR

echo "[S03] Building count matrices..."

$TRINITY_HOME/util/abundance_estimates_to_matrix.pl \
    --est_method RSEM \
    --gene_trans_map $GTM \
    --out_prefix RSEM \
    --name_sample_by_basedir \
    --quant_files $QUANT \
    --min_tpm 0

echo "[S03] Done. Count matrices written to $DATADIR/RSEM.gene.counts.matrix etc."
