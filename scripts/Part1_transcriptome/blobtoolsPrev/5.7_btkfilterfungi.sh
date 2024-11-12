#!/bin/bash
#SBATCH --job-name=btk_filter_keep_fungiReads
#SBATCH --output=/home/alicebalard/Scripts/AliceScripts/cyanochytridMET/scripts/logs_dir/%x.%j.out
#SBATCH --error=/home/alicebalard/Scripts/AliceScripts/cyanochytridMET/scripts/logs_dir/%x.%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=5G
#SBATCH --time=2-24:00:00
#SBATCH --qos=standard

## Enable same name autoswapping
LMOD_DISABLE_SAME_NAME_AUTOSWAP=no

module purge
module load Anaconda3

source ~/.bashrc
conda init --all
conda activate btk

module load BEDTools
module load seqtk
module load GCCcore/12.3.0

##############################################################
## First assembly: Z1 to Z12, chytrids
ASSEMBLYDIR=/scratch/alicebalard/outData/assembly
READS1=$ASSEMBLYDIR/combined_left.fq
READS2=$ASSEMBLYDIR/combined_right.fq
ASSEMBLY=$ASSEMBLYDIR/trinity_out_dir/Trinity.fasta 

BTKOUT=/scratch/alicebalard/outData/blobtools/Z1Z12
cd $BTKOUT

## https://github.com/blobtoolkit/blobtoolkit/issues/49
##echo "blobtools filter make table"
#blobtools filter --param bestsumorder_kingdom--Keys=Fungi --invert --table tableFungi.tsv $BTKOUT
##echo "table done"

## need to run in a job
##echo "bedtools make list"
#bedtools bamtobed -i $BTKOUT/assembly.reads.bam | grep -wFf <(tail -n +2 tableFungi.tsv | cut -f 2) | cut -f 4 > nameFungi.lst
##echo "list made"

## Sort and keep unique reads names
sort -u nameFungi.lst > nameFungi_unique.lst

## rm /1 or /2 at the end of each line: we will keep both sides
sed 's/..$//' < nameFungi_unique.lst > nameFungiUnpaired.lst

## Sort and keep unique reads names
sort -u nameFungiUnpaired.lst > nameFungiUniqueUnpaired.lst

## Keep only the reads associated with fungi
echo "seqtk filter reads"
seqtk subseq $READS1 nameFungiUniqueUnpaired.lst > filtered/out_1.fq
seqtk subseq $READS2 nameFungiUniqueUnpaired.lst > filtered/out_2.fq
echo "Reads filtered!"

##############################################################
## Second assembly: In1 to In12, chytrids
ASSEMBLYDIR=/scratch/alicebalard/outData/assembly_In
READS1=$ASSEMBLYDIR/combined_left.fq
READS2=$ASSEMBLYDIR/combined_right.fq
ASSEMBLY=$ASSEMBLYDIR/trinity_out_dir/Trinity.fasta 

BTKOUT=/scratch/alicebalard/outData/blobtools/In1In12
cd $BTKOUT

#### https://github.com/blobtoolkit/blobtoolkit/issues/49
##echo "blobtools filter make table"
##blobtools filter --param bestsumorder_kingdom--Keys=Fungi --invert --table tableFungi.tsv $BTKOUT
##echo "table done"

#### need to run in a job
##echo "bedtools make list"
##bedtools bamtobed -i $BTKOUT/assembly.reads.bam | grep -wFf <(tail -n +2 tableFungi.tsv | cut -f 2) | cut -f 4 > nameFungi.lst
##echo "list made"

## Sort and keep unique reads names
sort -u nameFungi.lst > nameFungi_unique.lst

## rm /1 or /2 at the end of each line: we will keep both sides
sed 's/..$//' < nameFungi_unique.lst > nameFungiUnpaired.lst

## Sort and keep unique reads names
sort -u nameFungiUnpaired.lst > nameFungiUniqueUnpaired.lst

## Keep only the reads associated with fungi
echo "seqtk filter reads"
seqtk subseq $READS1 nameFungiUniqueUnpaired.lst > filtered/out_1.fq
seqtk subseq $READS2 nameFungiUniqueUnpaired.lst > filtered/out_2.fq
echo "Reads filtered!"

