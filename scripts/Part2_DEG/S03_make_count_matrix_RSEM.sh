#!/bin/bash

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

# Enable autoswapping of same-name modules
export LMOD_DISABLE_SAME_NAME_AUTOSWAP=no

# Load the Trinity, R and foss modules
module load foss/2019b
module add Trinity/2.10.0-foss-2019b-Python-3.7.4
module load R-bundle-Bioconductor/3.16-foss-2022b-R-4.2.2

## Created in S01:
DATADIR=/scratch/alicebalard/outData/RSEM

cd $DATADIR

## Automatically make listing S02 results:
find $DATADIR/ -name "RSEM.genes.results" | sort > $DATADIR/quant_file_genes.txt

# Build transcript and gene expression matrices
$TRINITY_HOME/util/abundance_estimates_to_matrix.pl \
    --est_method RSEM \
    --gene_trans_map 'none' \
    --out_prefix RSEM_new_hope \
    --name_sample_by_basedir \
    --quant_files $DATADIR/quant_file_genes.txt \
    --min_tpm 0

## no need for gene trans map at this stage as it's already genes agregated

## rename for clarity
for f in $DATADIR/RSEM_new_hope.isoform.*; do
    newname=$(echo $f | sed 's/isoform/gene/')
    mv $f $newname
    echo "$f → $newname"
done

## Move input for DESeq2 where useful
cp $DATADIR/RSEM_new_hope.gene.counts.matrix /home/alicebalard/Scripts/AliceScripts/cyanochytridMET/data/RSEM_new_hope.gene.counts.matrix
