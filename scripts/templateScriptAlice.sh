#!/bin/bash

#SBATCH --job-name=my_serial_job                # replace name
#SBATCH --mail-user=alicebalard@zedat.fu-berlin.de  
#SBATCH --mail-type=end
#SBATCH --ntasks=1                             
#SBATCH --mem-per-cpu=10                      
#SBATCH --time=08:00:00                       
#SBATCH --qos=standard              

## module add ExampleProg/1.2.3-foss-2018b         # load module if needed

cd /scratch/alicebalard/
exampleprog_serial                              # replace with your program

