#!/bin/bash

#SBATCH --job-name=7.3_blobtk_addcov
#SBATCH --error=/home/alicebalard/Scripts/AliceScripts/cyanochytridMET/scripts/logs_dir/%x.%j.err
#SBATCH --output=/home/alicebalard/Scripts/AliceScripts/cyanochytridMET/scripts/logs_dir/%x.%j.out
#SBATCH --mail-user=alicebalard@zedat.fu-berlin.de
#SBATCH --mail-type=end
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=21GB 
#SBATCH --time=1-24:00:00
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

## https://blobtoolkit.genomehubs.org/blobtools2/blobtools2-tutorials/getting-started-with-blobtools2/#open_dataset

num_threads=$SLURM_CPUS_PER_TASK

module load minimap2/2.24-GCCcore-11.3.0
module load SAMtools/1.17-GCC-12.2.0

## Final assembly
ZfungiL=/scratch/alicebalard/outData/blobtools/Z1Z12/filtered/out_1.fq
ZfungiR=/scratch/alicebalard/outData/blobtools/Z1Z12/filtered/out_2.fq
InfungiL=/scratch/alicebalard/outData/blobtools/In1In12/filtered/out_1.fq
InfungiR=/scratch/alicebalard/outData/blobtools/In1In12/filtered/out_2.fq

cat $ZfungiL $InfungiL > /scratch/alicebalard/outData/assemblyFinal/combined_left.fq
cat $ZfungiR $InfungiR > /scratch/alicebalard/outData/assemblyFinal/combined_right.fq

ASSEMBLY=/scratch/alicebalard/outData/assemblyFinal/trinity_out_dir/Trinity.fasta 
READS1=/scratch/alicebalard/outData/assemblyFinal/combined_left.fq
READS2=/scratch/alicebalard/outData/assemblyFinal/combined_right.fq
BTKOUT=/scratch/alicebalard/outData/blobtools/FINALTran

echo "Add mapping coverage..."
minimap2 -ax sr -t $num_threads $ASSEMBLY \
         $READS1 $READS2 \
| samtools sort -@$num_threads -O BAM -o $BTKOUT/assembly.reads.bam -

# The resulting BAM file(s) can then be imported into a BlobDir:
echo "Import results cov in blobdir..."
blobtools add --cov $BTKOUT/assembly.reads.bam $BTKOUT
echo "Done!"
