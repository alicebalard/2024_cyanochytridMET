#!/bin/sh --login
#SBATCH --job-name=runSortmeRNAmissing
#SBATCH --output=/home/alicebalard/Scripts/AliceScripts/cyanochytridMET/scripts/logs_dir/%x.%j.out
#SBATCH --error=/home/alicebalard/Scripts/AliceScripts/cyanochytridMET/scripts/logs_dir/%x.%j.err
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=20GB 
#SBATCH --time=24:00:00
#SBATCH --qos=standard              
#SBATCH --mail-user=alicebalard@zedat.fu-berlin.de
#SBATCH --mail-type=end
#SBATCH --partition=main,begendiv
#SBATCH --constraint=no_gpu

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

## missing file
FILEF=/scratch/alicebalard/RawData/Z12_1_trimm_paired.fastq.gz
FILER=/scratch/alicebalard/RawData/Z12_2_trimm_paired.fastq.gz

num_threads=$SLURM_CPUS_PER_TASK
echo "Running with $num_threads threads"

sortmerna -ref $REF -reads $FILEF -reads $FILER --aligned Z12_rRNA_hits --other $Z12_non_rRNA --paired_in --fastx --out2 --threads $num_threads --workdir /scratch/alicebalard/outData/sortmerna/Z12

## --paired_in will keep the read pairs together by putting both reads (aligned and unaligned) into the aligned.fasta output file
## --paired_out will keep the read pairs together by putting both reads (aligned and unaligned) into the other.fasta output file 
## --out2 makes 2 files

