#!/bin/bash

#SBATCH --job-name=mvrawfasta
#SBATCH --mail-user=alicebalard@zedat.fu-berlin.de  
#SBATCH --mail-type=end
#SBATCH --ntasks=1                             
#SBATCH --mem-per-cpu=1000
#SBATCH --time=20:00:00                       
#SBATCH --qos=standard              

scp /scratch/strassert/EN00002253_hdd1/erika/*_paired_*fastq.gz /scratch/alicebalard/RawData/.
