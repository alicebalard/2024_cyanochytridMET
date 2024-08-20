#!/bin/bash

#SBATCH --job-name=rmolddir
#SBATCH --mail-user=alicebalard@zedat.fu-berlin.de  
#SBATCH --mail-type=end
#SBATCH --ntasks=1                             
#SBATCH --mem-per-cpu=10                      
#SBATCH --time=08:00:00                       
#SBATCH --qos=standard              

rm -rf /scratch/alicebalard/outData/assembly_v1/

