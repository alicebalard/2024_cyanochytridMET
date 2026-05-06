#!/bin/bash
# S02_mappingBowtie.sh
# Maps reads to combined transcriptome and estimates abundance with RSEM
#SBATCH --job-name=align_rsem
#SBATCH --mail-user=alicebalard@zedat.fu-berlin.de
#SBATCH --mail-type=end,fail
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=2GB
#SBATCH --time=24:00:00
#SBATCH --output=/home/alicebalard/Scripts/AliceScripts/cyanochytridMET/scripts/logs_dir/%x.%j.out
#SBATCH --error=/home/alicebalard/Scripts/AliceScripts/cyanochytridMET/scripts/logs_dir/%x.%j.err
#SBATCH --qos=standard
#SBATCH --partition=main,begendiv
#SBATCH --constraint=no_gpu

set -euo pipefail

export LMOD_DISABLE_SAME_NAME_AUTOSWAP=no

module add Trinity/2.10.0-foss-2019b-Python-3.7.4
module add Bowtie/1.3.1-GCC-11.3.0
module add RSEM/1.3.3-foss-2022a

## Outputs from S01 (now all under same directory)
DATADIR=/scratch/alicebalard/outData/RSEM
ASSEMBLY_BOTH=$DATADIR/assembly_both.fna
GTM=$DATADIR/combined_gene_trans_map.txt

cd $DATADIR

## Manually prepared sample file (tab-separated: condition, replicate, fq1, [fq2])
cp /scratch/erikamr/cyano_chytrid_met/data/samples_file.txt .
SAMPLE_FILE=$DATADIR/samples_file.txt

OUTDIR=$DATADIR/out_trinity_align_rsem
mkdir -p $OUTDIR

echo "[S02] Starting alignment and abundance estimation..."
echo "[S02] Assembly: $ASSEMBLY_BOTH"
echo "[S02] GTM: $GTM"
echo "[S02] Samples: $SAMPLE_FILE"

$TRINITY_HOME/util/align_and_estimate_abundance.pl \
    --transcripts $ASSEMBLY_BOTH \
    --seqType fq \
    --samples_file $SAMPLE_FILE \
    --est_method RSEM \
    --aln_method bowtie \
    --bowtie_RSEM "--no-unal -k 20" \
    --thread_count 8 \
    --output_dir $OUTDIR \
    --gene_trans_map $GTM \
    --prep_reference

echo "[S02] Done. Results in $OUTDIR"
