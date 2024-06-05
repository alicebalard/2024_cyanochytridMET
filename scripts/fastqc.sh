#!/bin/bash
#SBATCH --job-name=fastqc
#SBATCH --output=fastqc_%A_%a.out
#SBATCH --error=fastqc_%A_%a.err
#SBATCH --array=1-73
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=3:00:00
#SBATCH --qos=standard 

# Get the current file based on the array task ID
FILE=$(ls /scratch/alicebalard/RawData/*fastq.gz | sed -n "${SLURM_ARRAY_TASK_ID}p")

/scratch/alicebalard/FastQC/fastqc -o /scratch/alicebalard/outData/fastqc $FILE

