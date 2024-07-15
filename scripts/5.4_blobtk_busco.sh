#!/bin/sh --login
#SBATCH --job-name=BUSCO
#SBATCH --output=/home/alicebalard/Scripts/AliceScripts/cyanochytridMET/scripts/logs_dir/%x.%j.out
#SBATCH --error=/home/alicebalard/Scripts/AliceScripts/cyanochytridMET/scripts/logs_dir/%x.%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=4GB 
#SBATCH --time=24:00:00
#SBATCH --partition=main,begendiv
#SBATCH --constraint=no_gpu
#SBATCH --qos=standard              
#SBATCH --mail-user=alicebalard@zedat.fu-berlin.de
#SBATCH --mail-type=end

## Enable same name autoswapping
LMOD_DISABLE_SAME_NAME_AUTOSWAP=no

module purge
module load Anaconda3

source ~/.bashrc
conda init --all
conda activate btk

module add BUSCO/5.1.2-foss-2020b

## First assembly Z1Z12 only chytrids:
BTKOUT=/scratch/alicebalard/outData/blobtools/Z1Z12assembly
TRANSC=/scratch/alicebalard/outData/assembly/trinity_out_dir/Trinity.fasta

## First assembly: Z1 to Z12, only chytrids
#ASSEMBLY=/scratch/alicebalard/outData/assembly/trinity_out_dir/Trinity.fasta 
#BTKOUT=/scratch/alicebalard/outData/blobtools/Z1Z12assembly

## Second assembly: In1 to In12, chytrids infected by bacteria
ASSEMBLY=/scratch/alicebalard/outData/assembly_In/trinity_out_dir/Trinity.fasta
BTKOUT=/scratch/alicebalard/outData/blobtools/In1In12assembly

cd $BTKOUT

## With more samples
#busco -i $ASSEMBLY -l fungi_odb10 -m transcriptome -c 10 -o "BUSCO_Z1Z12"    
busco -i $ASSEMBLY -l fungi_odb10 -m transcriptome -c 10 -o "BUSCO_In1In12"

## l = current fungi database (will be downloaded automatically)
## m = genome or transcriptome (in your case transcriptome)
## NOTE: BUSCO doesn't like slash signs in the fasta header, you may need to replace them before running BUSCO.

## Add to blob dir
# These files can be imported to add BUSCO annotations to the assembly contigs:

blobtools add \
#    --busco $BTKOUT/BUSCO_Z1Z12/run_fungi_odb10/full_table.tsv \
     --busco $BTKOUT/BUSCO_In1In12/run_fungi_odb10/full_table.tsv
    $BTKOUT

