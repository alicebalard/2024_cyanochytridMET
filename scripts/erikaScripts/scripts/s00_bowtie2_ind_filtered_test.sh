#!/bin/bash

#SBATCH --mail-user=erikamr@zedat.fu-berlin.de
#SBATCH --mail-type=end
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=16GB
#SBATCH --time=2:00:00
#SBATCH --error=/scratch/erikamr/cyano_chytrid_met/scripts/%x.%j.err
#SBATCH --output=/scratch/erikamr/cyano_chytrid_met/scripts/%x.%j.out
#SBATCH --qos=standard
#SBATCH --partition=main,begendiv
#SBATCH --constraint=no_gpu

# Enable autoswapping of same-name modules
export LMOD_DISABLE_SAME_NAME_AUTOSWAP=no

module add Bowtie2/2.4.5-GCC-11.3.0
module add SAMtools/1.16.1-GCC-11.3.0

# Define path for the input FASTA file
FASTA_FILE="/scratch/erikamr/cyano_chytrid_met/data/assembly_both_cds.fna"

# Build Bowtie2 index
bowtie2-build ${FASTA_FILE} ${FASTA_FILE}

#  Alignment
bowtie2 -p 2 -q --no-unal -k 20 -x $FASTA_FILE -1 /scratch/alicebalard/outData/sortmerna/Z1_non_rRNA_fwd.fq.gz -2 /scratch/alicebalard/outData/sortmerna/Z1_non_rRNA_rev.fq.gz  2>/scratch/erikamr/cyano_chytrid_met/data/aligned_bowtie2_ind/align_stats_test.txt | samtools view -@2 -Sb -o /scratch/erikamr/cyano_chytrid_met/data/aligned_bowtie2_ind/bowtie2_test.bam
