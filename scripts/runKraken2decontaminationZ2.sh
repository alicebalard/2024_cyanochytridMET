#!/bin/sh --login
#SBATCH --job-name=runKraken
#SBATCH --error=/home/alicebalard/Scripts/AliceScripts/cyanochytridMET/scripts/logs_dir/%x.%j.err
#SBATCH --array=1-12
#SBATCH --cpus-per-task=3
#SBATCH --mem-per-cpu=20GB 
#SBATCH --partition=main,begendiv
#SBATCH --constraint=no_gpu
#SBATCH --time=03-24:00:00
#SBATCH --qos=standard              
#SBATCH --mail-user=alicebalard@zedat.fu-berlin.de
#SBATCH --mail-type=end

## https://laramieakozbek.com/using-kraken2-bbduk-to-check-for-contamination-in-short-read-datasets/

module purge
module load Anaconda3

source ~/.bashrc
conda init --all
conda activate kraken

DBNAME=krakenDB

## We want to clean samples Z1 to Z12 (only chytrid)
rawDir=/scratch/alicebalard/RawData/
extension="_trimm_paired.fastq.gz"

cd /scratch/alicebalard/outData/kraken2decontamination/
## Create the list to array over
for i in {1..12}; do echo "Z$i" >> listfilestemp ; done  

## Get the sample (e.g. Z1)
sample=$(sed "${SLURM_ARRAY_TASK_ID}q;d" listfilestemp)

## Get the full path of the sample
R1="${rawDir}${sample}_1${extension}"    
R2="${rawDir}${sample}_2${extension}"

num_threads=$SLURM_CPUS_PER_TASK
echo "Running with $num_threads threads"

cd /scratch/alicebalard/outData/kraken2decontamination/

# Standard DB library downloaded it 20240112
# DB with RefSeq bacterial, archaeal, and viral domains plus the human genome and a handful of other important vectors

## Run Kraken2
echo "run Kraken2"
kraken2 --db $DBNAME --threads $num_threads --unclassified-out out/${sample}_unclassified_raw_reads\#.fq --classified-out out/${sample}_classified_raw_reads\#.fq --paired --out "${sample}_kraken2out" --report "${sample}_kraken2report" --use-names --report-zero-counts $R1 $R2
echo "run Kraken2 FINISHED!"
