#!/bin/sh --login
#SBATCH --job-name=runSortmeRNA
#SBATCH --output=/home/alicebalard/Scripts/AliceScripts/cyanochytridMET/scripts/logs_dir/%x.%j.out
#SBATCH --error=/home/alicebalard/Scripts/AliceScripts/cyanochytridMET/scripts/logs_dir/%x.%j.err
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=20GB 
#SBATCH --time=05-24:00:00
#SBATCH --qos=standard              
#SBATCH --mail-user=alicebalard@zedat.fu-berlin.de
#SBATCH --mail-type=end
#SBATCH --partition=main,begendiv
#SBATCH --constraint=no_gpu
#SBATCH --array=1-24

module purge
module load Anaconda3

source ~/.bashrc
conda init --all
conda activate sortmerna_env

cd /scratch/alicebalard/outData/sortmerna

## sortmeRNA requires a database that I downloaded beforehand
#wget https://github.com/biocore/sortmerna/releases/download/v4.3.4/database.tar.gz
#mkdir rRNA_databases_v4
#tar -xvf database.tar.gz -C rRNA_databases_v4

REF=/scratch/alicebalard/outData/sortmerna/rRNA_databases_v4/smr_v4.3_default_db.fasta 

## references of the files to clean
## We want to clean samples Z1 to Z12 (only chytrid) 
rawDir=/scratch/alicebalard/RawData/
extension="_trimm_paired.fastq.gz"

## Create the list to array over OUTSIDE of the script
## done before for Zxx
# for i in {1..12}; do echo "Z$i" >> listfilestemp ; done

# for i in {1..12}; do echo "C$i" >> listfilestemp ; done
# for i in {1..12}; do echo "In$i" >> listfilestemp ; done

## Get the sample (e.g. Z1)
sample=$(sed "${SLURM_ARRAY_TASK_ID}q;d" listfilestemp)

## Get the full path of the sample
FILEF="${rawDir}${sample}_1${extension}"    
FILER="${rawDir}${sample}_2${extension}"

num_threads=$SLURM_CPUS_PER_TASK
echo "Running with $num_threads threads"

sortmerna -ref $REF -reads $FILEF -reads $FILER --aligned ${sample}_rRNA_hits --other ${sample}_non_rRNA --paired_in --fastx --out2 --threads $num_threads --workdir /scratch/alicebalard/outData/sortmerna/${sample}

## --paired_in will keep the read pairs together by putting both reads (aligned and unaligned) into the aligned.fasta output file
## --paired_out will keep the read pairs together by putting both reads (aligned and unaligned) into the other.fasta output file 
## --out2 makes 2 files

