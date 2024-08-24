#!/bin/bash

#SBATCH --job-name=blobtk_makeblobdir
#SBATCH --error=/home/alicebalard/Scripts/AliceScripts/cyanochytridMET/scripts/logs_dir/%x.%j.err
#SBATCH --output=/home/alicebalard/Scripts/AliceScripts/cyanochytridMET/scripts/logs_dir/%x.%j.out
#SBATCH --mail-user=alicebalard@zedat.fu-berlin.de  
#SBATCH --mail-type=end
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=5GB 
#SBATCH --time=24:00:00
#SBATCH --qos=standard              
#SBATCH --partition=main,begendiv
#SBATCH --constraint=no_gpu

## Enable same name autoswapping
LMOD_DISABLE_SAME_NAME_AUTOSWAP=no

module purge
module load Anaconda3

source ~/.bashrc
conda init --all
conda activate btk

cd /scratch/alicebalard/outData/blobtools

ASSEMBLY=/scratch/alicebalard/outData/assemblyFinal/trinity_out_dir/Trinity.fasta 
OUT="FINALTran"

## Create blobdir
echo "Create a Blobdir with a fasta assembly..."
blobtools create --fasta $ASSEMBLY $OUT

# creates a new directory in the location specified by the last argument that contains a set of files containing values for GC-content (gc.json), length (length.json), number of Ns (ncount.json) and sequence names (identifiers.json) for each sequence in the assembly. A final file (meta.json) contains metadata for the dataset describing the datatypes of the available fields and the ranges of values for each of these fields
