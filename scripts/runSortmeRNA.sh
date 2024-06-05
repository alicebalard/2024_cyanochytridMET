#!/bin/sh --login
#SBATCH --job-name=runSortmeRNA
#SBATCH --output=/scratch/alicebalard/code/logs_dir/%x.%j.out
#SBATCH --error=/scratch/alicebalard/code/logs_dir/%x.%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=10GB 
#SBATCH --time=03-24:00:00
#SBATCH --qos=standard              
#SBATCH --mail-user=alicebalard@zedat.fu-berlin.de
#SBATCH --mail-type=end

module purge
module load Anaconda3

source ~/.bashrc
conda init --all
conda activate sortmerna_env

REF=/scratch/alicebalard/outData/sortmerna/rRNA_databases_v4/smr_v4.3_default_db.fasta 
FILEF=/scratch/alicebalard/RawData/Z2_1_trimm_paired.fastq.gz 
FILER=/scratch/alicebalard/RawData/Z2_2_trimm_paired.fastq.gz 

cd /scratch/alicebalard/outData/sortmerna

sortmerna -ref $REF -reads $FILEF -reads $FILER --aligned --other --paired_in --fastx --out2 --threads 20 

## --paired_in will keep the read pairs together by putting both reads (aligned and unaligned) into the aligned.fasta output file .
## --paired_out will keep the read pairs together by putting both reads (aligned and unaligned) into the other.fasta output file 
## --out2 makes 2 files

