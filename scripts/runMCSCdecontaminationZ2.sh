#!/bin/sh
#SBATCH --job-name=runMCSC 
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=3-24:00:00
#SBATCH --mail-user=alicebalard@zedat.fu-berlin.de
#SBATCH --mail-type=end
#SBATCH --ntasks=1                             
#SBATCH --qos=standard    

module purge
module load DIAMOND/2.0.13-GCC-11.2.0

## Enable same name autoswapping
LMOD_DISABLE_SAME_NAME_AUTOSWAP=no
module load R/4.3.2-gfbf-2023a

module load Ghostscript

sh /scratch/alicebalard/MCSC/MCSC_Decontamination/MCSC_decontamination.sh /scratch/alicebalard/code/param.MCSC.ini

