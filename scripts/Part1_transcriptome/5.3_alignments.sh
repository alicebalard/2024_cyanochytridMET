#!/bin/bash
#SBATCH --job-name=alignments
#SBATCH --error=/home/alicebalard/Scripts/AliceScripts/cyanochytridMET/scripts/logs_dir/%x.%j.err
#SBATCH --output=/home/alicebalard/Scripts/AliceScripts/cyanochytridMET/scripts/logs_dir/%x.%j.out
#SBATCH --mail-user=alicebalard@zedat.fu-berlin.de
#SBATCH --mail-type=end
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=40GB 
#SBATCH --time=1-24:00:00
#SBATCH --qos=standard              
#SBATCH --partition=main,begendiv
#SBATCH --constraint=no_gpu

## Enable same name autoswapping
LMOD_DISABLE_SAME_NAME_AUTOSWAP=no

module purge

source ~/.bashrc

cd /scratch/alicebalard/outData/alignments

num_threads=$SLURM_CPUS_PER_TASK

module load minimap2/2.29-GCCcore-13.3.0
module load SAMtools/1.21-GCC-13.3.0

###########################################
## First assembly: Z1 to Z12, only chytrids
ASSEMBLY=/scratch/alicebalard/outData/assemblies/assembly_Z/Trinity.fasta 
READS1=/scratch/alicebalard/outData/assemblies/assembly_Z/combined_left.fq
READS2=/scratch/alicebalard/outData/assemblies/assembly_Z/combined_right.fq

minimap2 -ax sr -t $num_threads $ASSEMBLY $READS1 $READS2 | samtools sort -@$num_threads -m 1G -O BAM -o assemblyZ.reads.bam -
## -m 1G aps samtools sort at about 1 GB per thread, which is much less likely to trigger OOM

##########################################################################
## Second assembly: In1 to In12, bacteria infected by chytrids + coculture
ASSEMBLY=/scratch/alicebalard/outData/assemblies/assembly_In_coculture/Trinity.fasta 
READS1=/scratch/alicebalard/outData/assemblies/assembly_In_coculture/combined_left.fq
READS2=/scratch/alicebalard/outData/assemblies/assembly_In_coculture/combined_right.fq

minimap2 -ax sr -t $num_threads $ASSEMBLY $READS1 $READS2 | samtools sort -@$num_threads -m 1G -O BAM -o assemblyIn.reads.bam -
