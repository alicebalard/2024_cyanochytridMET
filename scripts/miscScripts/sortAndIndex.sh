#!/bin/bash
#SBATCH --job-name=sort_index_bam
#SBATCH --output=/home/alicebalard/Scripts/AliceScripts/cyanochytridMET/scripts/logs_dir/%x.%j.out
#SBATCH --error=/home/alicebalard/Scripts/AliceScripts/cyanochytridMET/scripts/logs_dir/%x.%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=5G
#SBATCH --time=1-24:00:00
#SBATCH --qos=standard

## take > 5h for 24G bam

# Load necessary modules on Curta
module load SAMtools

# Assign input arguments to variables
INPUT_BAM="$1"

# Sort the BAM file
echo "Sorting BAM file..."
samtools sort -@ 10 -m 4G -o $INPUT_BAM.sorted.bam $INPUT_BAM
echo "BAM file sorted!"

# Index the sorted BAM file
echo "Indexing sorted BAM file..."
samtools index $INPUT_BAM.sorted.bam
echo "BAM file indexed!"

echo "Job completed :)"
