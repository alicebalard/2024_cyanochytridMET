#!/bin/bash

#SBATCH --job-name=8.4_trinotate
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

source /home/alicebalard/Scripts/AliceScripts/cyanochytridMET/scripts/8.1_finalTransAnnotation_prepareTrinotate.sh

## transcripts.fasta : your target transcriptome in fasta forma
transcripts=/scratch/alicebalard/outData/assemblyFinal/trinity_out_dir/Trinity.filtered.Fungi.fasta

## coding_seqs.pep : coding regions translated in fasta format (specific header formatting required - see below. Most use TransDecoder to generate this)
coding_seqs=/scratch/alicebalard/outData/annotation/Trinotate/TransDecoder-TransDecoder-v5.7.1/Trinity.filtered.Fungi.fasta.transdecoder.pep

## gene_to_trans_map.tsv : pairwise mappings between gene and transcript isoform identifiers
gene_to_trans_map=/scratch/alicebalard/outData/assemblyFinal/trinity_out_dir/Trinity.fasta.gene_trans_map

## Initialize Trinotate sqlite database:
## $TRINOTATE_HOME/Trinotate-Trinotate-v4.0.2/Trinotate --db chytridTrinotate.sqlite --init \
##            --gene_trans_map $gene_to_trans_map \
##            --transcript_fasta $transcripts \
##            --transdecoder_pep $coding_seqs

## Running sequence analyses

## Run all the other analyses:
$TRINOTATE_HOME/Trinotate-Trinotate-v4.0.2/Trinotate --db chytridTrinotate.sqlite --CPU 20 \
           --transcript_fasta $transcripts \
           --transdecoder_pep $coding_seqs \
	   --trinotate_data_dir $TRINOTATE_DATA_DIR \
           --run "swissprot_blastp swissprot_blastx pfam infernal" \
           --use_diamond

## Run signalp6 (in previous script) and tmhmm outside because of python versions issues:
conda activate myannot
tmhmm-2.0c/bin/tmhmm --short $coding_seqs  > tmhmm.v2.out
conda deactivate

## Add to SQLite
$TRINOTATE_HOME/Trinotate-Trinotate-v4.0.2/Trinotate --db chytridTrinotate.sqlite --LOAD_tmhmmv2 tmhmm.v2.out

## Generate Trinotate report:
$TRINOTATE_HOME/Trinotate-Trinotate-v4.0.2/Trinotate --db chytridTrinotate.sqlite --report --incl_pep --incl_trans > $TRINOTATE_HOME/allFungiTrinotate.tsv

## Simplify output for later use with DESeq2:
cat $TRINOTATE_HOME/allFungiTrinotate.tsv | cut -f 1,2,3,13 > $TRINOTATE_HOME/allFungiTrinot_simplified.tsv

awk 'BEGIN {OFS="\t"} 
     NR==1 {print $0, "gene_name"} 
     NR > 1 {split($3, a, "^"); print $0, a[1]}' $TRINOTATE_HOME/allFungiTrinot_simplified.tsv > $TRINOTATE_HOME/temp

mv $TRINOTATE_HOME/temp $TRINOTATE_HOME/allFungiTrinot_simplified.tsv

cp TRINOTATE_HOME/allFungiTrinot_simplified.tsv /home/alicebalard/Scripts/AliceScripts/cyanochytridMET/figTab/.

