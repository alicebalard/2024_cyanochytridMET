#!/bin/bash
#SBATCH --job-name=fastqc
#SBATCH --error=/home/alicebalard/Scripts/AliceScripts/cyanochytridMET/scripts/logs_dir/%x.%j.err
#SBATCH --output=/home/alicebalard/Scripts/AliceScripts/cyanochytridMET/scripts/logs_dir/%x.%j.out
#SBATCH --array=1-24
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=3:00:00
#SBATCH --qos=standard 

# Get the current file based on the array task ID
FILE=$(ls /scratch/alicebalard/outData/sortmerna/*non_rRNA*fwd.fq.gz | sed -n "${SLURM_ARRAY_TASK_ID}p")

/scratch/alicebalard/FastQC/fastqc -o /scratch/alicebalard/outData/fastqc/afterSortmerna $FILE

