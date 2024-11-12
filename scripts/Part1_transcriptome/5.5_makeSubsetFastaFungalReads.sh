#!/bin/bash
#SBATCH --job-name=5.5_makeSubsetReads
#SBATCH --output=/home/alicebalard/Scripts/AliceScripts/cyanochytridMET/scripts/logs_dir/%x.%j.out
#SBATCH --error=/home/alicebalard/Scripts/AliceScripts/cyanochytridMET/scripts/logs_dir/%x.%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=5G
#SBATCH --time=04:00:00
#SBATCH --qos=standard

## Enable same name autoswapping
LMOD_DISABLE_SAME_NAME_AUTOSWAP=no

module purge
source ~/.bashrc
module load BEDTools
module load seqtk
module load GCCcore/12.3.0

Z_FUNGIREADS_NAMES=/scratch/alicebalard/outData/alignments/assemblyZ.reads_filteredFungi
Z_READS1=/scratch/alicebalard/outData/assembly/combined_left.fq
Z_READS2=/scratch/alicebalard/outData/assembly/combined_right.fq

In_FUNGIREADS_NAMES=/scratch/alicebalard/outData/alignments/assemblyIn.reads_filteredFungi
In_READS1=/scratch/alicebalard/outData/assembly_In/combined_left.fq
In_READS2=/scratch/alicebalard/outData/assembly_In/combined_right.fq

cd /scratch/alicebalard/outData/assemblyMergedFungi

## Keep reads that align paired-end to our fungi transcripts
sed 's/\/[12]$//' $Z_FUNGIREADS_NAMES | uniq -c | awk '$1 == 2 {print $2}' > Z_fungireads_uniquepairedreads
sed 's/\/[12]$//' $In_FUNGIREADS_NAMES | uniq -c | awk '$1 == 2 {print $2}' > In_fungireads_uniquepairedreads

## Keep only the reads associated with fungi
echo "seqtk filter reads for Z..."
seqtk subseq $Z_READS1 Z_fungireads_uniquepairedreads > assemblyZ.reads_filteredFungi_left.fq
seqtk subseq $Z_READS2 Z_fungireads_uniquepairedreads > assemblyZ.reads_filteredFungi_right.fq
echo "Reads filtered!"

echo "seqtk filter reads for In..."
seqtk subseq $In_READS1 In_fungireads_uniquepairedreads > assemblyIn.reads_filteredFungi_left.fq
seqtk subseq $In_READS2 In_fungireads_uniquepairedreads > assemblyIn.reads_filteredFungi_right.fq
echo "Reads filtered!"
