#!/bin/bash

#SBATCH --job-name=rmolddir
#SBATCH --mail-user=alicebalard@zedat.fu-berlin.de  
#SBATCH --mail-type=end
#SBATCH --ntasks=1                             
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=10                      
#SBATCH --time=08:00:00                       
#SBATCH --qos=standard              

cd /scratch/alicebalard/outData/assembly_In/

find trinity_out_dir/ -maxdepth 3 -type d | xargs -P 20 -n 1 rm -rf
