#!/bin/bash

#SBATCH --job-name=7.4_trinotate
#SBATCH --error=/home/alicebalard/Scripts/AliceScripts/cyanochytridMET/scripts/logs_dir/%x.%j.err
#SBATCH --output=/home/alicebalard/Scripts/AliceScripts/cyanochytridMET/scripts/logs_dir/%x.%j.out
#SBATCH --mail-user=alicebalard@zedat.fu-berlin.de  
#SBATCH --mail-type=end
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=5GB 
#SBATCH --time=3-24:00:00
#SBATCH --qos=standard              
#SBATCH --partition=main,begendiv
#SBATCH --constraint=no_gpu

module purge
module load Anaconda3

source ~/.bashrc
conda init --all

## Install Trinotate and all dependencies in TRINOTATE_DATA_DIR
source /home/alicebalard/Scripts/AliceScripts/cyanochytridMET/scripts/Part1_transcriptome/7.1_finalTransAnnotation_prepareTrinotate.sh

## transcripts.fasta : your target transcriptome in fasta format
transcripts=/scratch/alicebalard/outData/assemblyMergedFungi/trinity_out_dir/Trinity.fasta

## coding_seqs.pep : coding regions translated in fasta format (specific header formatting required - see below. Most use TransDecoder to generate this)
coding_seqs=/scratch/alicebalard/outData/assemblyMergedFungi/annotation/transdecoder/Trinity.fasta.transdecoder.pep

## gene_to_trans_map.tsv : pairwise mappings between gene and transcript isoform identifiers
gene_to_trans_map=/scratch/alicebalard/outData/assemblyMergedFungi/trinity_out_dir/Trinity.fasta.gene_trans_map

OUTDIR=/scratch/alicebalard/outData/assemblyMergedFungi/annotation

cd $OUTDIR

## Initiate database in OUTDIR by dl important databases
DB=/scratch/alicebalard/outData/assemblyMergedFungi/annotation/assemblyMergedFungi.sqlite

$TRINOTATE_HOME/Trinotate-Trinotate-v4.0.2/Trinotate --create --db $DB --trinotate_data_dir $TRINOTATE_HOME/DATADIR --use_diamond

## Initialize Trinotate sqlite database:
$TRINOTATE_HOME/Trinotate-Trinotate-v4.0.2/Trinotate --db $DB --init \
            --gene_trans_map $gene_to_trans_map \
            --transcript_fasta $transcripts \
            --transdecoder_pep $coding_seqs

## Running sequence analyses

## Run all the other analyses:
$TRINOTATE_HOME/Trinotate-Trinotate-v4.0.2/Trinotate --db $DB --CPU 20 \
           --transcript_fasta $transcripts \
           --transdecoder_pep $coding_seqs \
	   --trinotate_data_dir $TRINOTATE_DATA_DIR \
           --run "swissprot_blastp swissprot_blastx pfam infernal" \
           --use_diamond

## Run signalp6 (in previous script 7.3)
## Add to SQLite
$TRINOTATE_HOME/Trinotate-Trinotate-v4.0.2/Trinotate --db $DB --LOAD_signalp sigP6outdir/output.gff3

## Run tmhmm outside because of python versions issues:
conda activate myannot
tmhmm-2.0c/bin/tmhmm --short $coding_seqs  > tmhmm.v2.out
conda deactivate

## Add to SQLite
$TRINOTATE_HOME/Trinotate-Trinotate-v4.0.2/Trinotate --db $DB --LOAD_tmhmmv2 tmhmm.v2.out

## Generate Trinotate report:
$TRINOTATE_HOME/Trinotate-Trinotate-v4.0.2/Trinotate --db $DB --report --incl_pep --incl_trans > assemblyMergedFungi.tsv

## Simplify output for later use with DESeq2:
cat assemblyMergedFungi.tsv | cut -f 1,2,3,13 > assemblyMergedFungi_simplified.tsv

awk 'BEGIN {OFS="\t"} 
     NR==1 {print $0, "gene_name"} 
     NR > 1 {split($3, a, "^"); print $0, a[1]}' assemblyMergedFungi_simplified.tsv > temp

mv temp assemblyMergedFungi_simplified.tsv

mv assemblyMergedFungi_simplified.tsv /home/alicebalard/Scripts/AliceScripts/cyanochytridMET/data/.
