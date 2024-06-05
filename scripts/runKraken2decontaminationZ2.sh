#!/bin/sh --login
#SBATCH --job-name=runKraken
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=60GB 
#SBATCH --time=03-24:00:00
#SBATCH --qos=standard              
#SBATCH --partition=main,begendiv
#SBATCH --constraint=no_gpu
#SBATCH --mail-user=alicebalard@zedat.fu-berlin.de
#SBATCH --mail-type=end

## https://laramieakozbek.com/using-kraken2-bbduk-to-check-for-contamination-in-short-read-datasets/

module purge
module load Anaconda3

source ~/.bashrc
conda init --all
conda activate kraken

DBNAME=krakenDB
R1=/scratch/alicebalard/RawData/Z2_1_trimm_paired.fastq.gz
R2=/scratch/alicebalard/RawData/Z2_2_trimm_paired.fastq.gz

cd /scratch/alicebalard/outData/kraken2decontamination/

## create DB with RefSeq bacterial, archaeal, and viral domains plus the human genome and a handful of other important vectors
kraken2-build --standard --db $DBNAME --threads 20

## Run Kraken2
kraken2 --db $DBNAME --threads 20 --unclassified-out unclassified_raw_reads\#.fq --classified-out classified_raw_reads\#.fq --paired --report --use-names $R1 $R2

