#!/bin/bash
#SBATCH --job-name=multiqc
#SBATCH --mail-user=alicebalard@zedat.fu-berlin.de  
#SBATCH --mail-type=end
#SBATCH --error=/home/alicebalard/Scripts/AliceScripts/cyanochytridMET/scripts/logs_dir/%x.%j.err
#SBATCH --output=/home/alicebalard/Scripts/AliceScripts/cyanochytridMET/scripts/logs_dir/%x.%j.out
#SBATCH --ntasks=1
#SBATCH --mem=5GB 
#SBATCH --time=01:00:00
#SBATCH --qos=standard              

module load MultiQC/1.9-foss-2020a-Python-3.8.2
## for reads after trimming
# multiqc /scratch/alicebalard/outData/fastqc/afterTrim -o /scratch/alicebalard/outData/multiQC/trimmedReads

## for reads after kraken2 decontamination
#multiqc /scratch/alicebalard/outData/fastqc/afterKraken2 -o /scratch/alicebalard/outData/multiQC/decontKrakenReads

## for reads after sortmerna decontamination
multiqc /scratch/alicebalard/outData/fastqc/afterSortmerna -o /scratch/alicebalard/outData/multiQC/sortmeRNA


