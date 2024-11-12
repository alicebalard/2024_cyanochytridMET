#!/bin/bash
#SBATCH --job-name=5.4_filterReads
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
module load GCCcore/12.3.0

Z_FUNGTRANSCRIPTS=/scratch/alicebalard/outData/diamondBlastX/assZ_Fungi_transcripts
Z_READS1=/scratch/alicebalard/outData/assembly/combined_left.fq
Z_READS2=/scratch/alicebalard/outData/assembly/combined_right.fq
Z_BAM=/scratch/alicebalard/outData/alignments/assemblyZ.reads.bam
Z_OUT=/scratch/alicebalard/outData/alignments/assemblyZ.reads_filteredFungi

In_FUNGTRANSCRIPTS=/scratch/alicebalard/outData/diamondBlastX/assIn_Fungi_transcripts
In_READS1=/scratch/alicebalard/outData/assembly_In/combined_left.fq
In_READS2=/scratch/alicebalard/outData/assembly_In/combined_right.fq
In_BAM=/scratch/alicebalard/outData/alignments/assemblyIn.reads.bam
In_OUT=/scratch/alicebalard/outData/alignments/assemblyIn.reads_filteredFungi

bedtools bamtobed -i $Z_BAM | grep -wFf <(tail -n +2 $Z_FUNGTRANSCRIPTS) | cut -f 4 | sort -u > $Z_OUT

bedtools bamtobed -i $In_BAM | grep -wFf <(tail -n +2 $In_FUNGTRANSCRIPTS) | cut -f 4 | sort -u > $In_OUT
