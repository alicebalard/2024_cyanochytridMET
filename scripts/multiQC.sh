#!/bin/bash

#SBATCH --job-name=multiqc
#SBATCH --mail-user=alicebalard@zedat.fu-berlin.de  
#SBATCH --mail-type=end
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8                             
#SBATCH --mem=60GB 
#SBATCH --time=01:00:00
#SBATCH --qos=standard              

module load MultiQC/1.9-foss-2020a-Python-3.8.2
multiqc /scratch/alicebalard/outData/fastqc
