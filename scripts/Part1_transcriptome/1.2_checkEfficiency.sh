#!/bin/bash
#SBATCH --job-name=1.2_checkEfficiency
#SBATCH --output=/home/alicebalard/Scripts/AliceScripts/cyanochytridMET/scripts/logs_dir/%x.%j.out
#SBATCH --error=/home/alicebalard/Scripts/AliceScripts/cyanochytridMET/scripts/logs_dir/%x.%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=5G
#SBATCH --time=2-24:00:00
#SBATCH --qos=standard

echo "number of reads Z trimmed:"
zcat /scratch/alicebalard/RawData/Z*trimm_paired.fastq.gz | echo $((`wc -l`/4))

echo "Z post rRNA filter:"
zcat /scratch/alicebalard/outData/sortmerna/Z*non_rRNA*fq.gz | echo $((`wc -l`/4))

echo "number of reads C trimmed:"
zcat /scratch/alicebalard/RawData/C*trimm_paired.fastq.gz | echo $((`wc -l`/4))

echo "C post rRNA filter:"
zcat /scratch/alicebalard/outData/sortmerna/C*non_rRNA*fq.gz | echo $((`wc -l`/4))

echo "number of reads In trimmed:"
zcat /scratch/alicebalard/RawData/In*trimm_paired.fastq.gz | echo $((`wc -l`/4))

echo "In post rRNA filter:"
zcat /scratch/alicebalard/outData/sortmerna/In*non_rRNA*fq.gz | echo $((`wc -l`/4))
