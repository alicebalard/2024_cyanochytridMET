#!/bin/sh
#SBATCH --job-name=checkQualityAssembly
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --mail-user=alicebalard@zedat.fu-berlin.de
#SBATCH --mail-type=end
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=60GB 
#SBATCH --time=03-24:00:00
#SBATCH --qos=standard              
#SBATCH --partition=main,begendiv
#SBATCH --constraint=no_gpu

module purge
module add Trinity/2.10.0-foss-2019b-Python-3.7.4

ASSEMBLY1="/scratch/alicebalard/outData/assembly/trinity_out_dir/Trinity.fasta"

cd /scratch/alicebalard/outData/qualityAssembly

echo "basic contig statistics"

# $TRINITY_HOME/util/TrinityStats.pl $ASSEMBLY1

echo "read content: Typical Trinity transcriptome assembly will have the vast majority of all reads mapping back to the assembly, and ~70-80%
of the mapped fragments found mapped as proper pairs" 

# First, build the bowtie2 index of the assembly
# bowtie2-build $ASSEMBLY1 Trinity_index

# Align reads to the assembly
bowtie2 --local --no-unal -x Trinity_index -q -1 ../../RawData/Z2_1_trimm_paired.fastq.gz -2 ../../RawData/Z2_2_trimm_paired.fastq.gz | samtools view -Sb -@ 10 | samtools sort -@ 10 -m 5G -n -o bowtie2.nameSorted.bam

# Produce alignment statistics
## $TRINITY_HOME/util/SAM_nameSorted_to_uniq_count_stats.pl bowtie2_out.nameSorted.bam

############################################
## Same for the decontaminated transcriptome

DECONT="/scratch/alicebalard/outData/MCSCdecontamination/Trinity_decont.fasta"
