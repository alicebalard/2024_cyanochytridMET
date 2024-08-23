#!/bin/bash

#SBATCH --job-name=6_de_novo_assembly_final
#SBATCH --error=/home/alicebalard/Scripts/AliceScripts/cyanochytridMET/scripts/logs_dir/%x.%j.err
#SBATCH --output=/home/alicebalard/Scripts/AliceScripts/cyanochytridMET/scripts/logs_dir/%x.%j.out
#SBATCH --mail-user=alicebalard@zedat.fu-berlin.de  
#SBATCH --mail-type=end
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=300GB 
#SBATCH --time=5-24:00:00
#SBATCH --qos=standard              
#SBATCH --partition=main,begendiv
#SBATCH --constraint=no_gpu

module purge
module add Trinity/2.10.0-foss-2019b-Python-3.7.4

cd /scratch/alicebalard/outData/assemblyFinal/
## the output will be in 'trinity_out_dir/Trinity.fasta'

## Reads to be assembled:
ZfungiL=/scratch/alicebalard/outData/blobtools/Z1Z12/filtered/out_1.fq
ZfungiR=/scratch/alicebalard/outData/blobtools/Z1Z12/filtered/out_2.fq
InfungiL=/scratch/alicebalard/outData/blobtools/In1In12/filtered/out_1.fq
InfungiR=/scratch/alicebalard/outData/blobtools/In1In12/filtered/out_2.fq

cat $ZfungiL $InfungiL > combined_left.fq
cat $ZfungiR $InfungiR > combined_right.fq

echo "Number of left reads assembled in total:"
cat combined_left.fq | echo $((`wc -l`/4))
echo "Number of right reads assembled in total:"
cat combined_right.fq | echo $((`wc -l`/4))

num_threads=$SLURM_CPUS_PER_TASK
echo "Transcriptome assembly running with $num_threads threads"
MAX_MEM=300G

Trinity --seqType fq --normalize_by_read_set --max_memory $MAX_MEM --left combined_left.fq --right combined_right.fq --CPU $num_threads
