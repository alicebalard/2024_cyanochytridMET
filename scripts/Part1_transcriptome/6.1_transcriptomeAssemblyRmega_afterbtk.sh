#!/bin/bash

#SBATCH --job-name=6_de_novo_assembly_final
#SBATCH --error=/home/alicebalard/Scripts/AliceScripts/cyanochytridMET/scripts/logs_dir/%x.%j.err
#SBATCH --output=/home/alicebalard/Scripts/AliceScripts/cyanochytridMET/scripts/logs_dir/%x.%j.out
#SBATCH --mail-user=alicebalard@zedat.fu-berlin.de  
#SBATCH --mail-type=end
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=770GB 
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

echo "Number of left reads assembled from Z and In transcriptomes:"
cat $ZfungiL | echo $((`wc -l`/4))
cat $InfungiL | echo $((`wc -l`/4))

echo "Number of right reads assembled from Z and In transcriptomes:"
cat $ZfungiR | echo $((`wc -l`/4))
cat $InfungiR | echo $((`wc -l`/4))

num_threads=$SLURM_CPUS_PER_TASK
echo "Transcriptome assembly running with $num_threads threads"
MAX_MEM=768G

Trinity --seqType fq --normalize_by_read_set --max_memory $MAX_MEM --left $ZfungiL,$InfungiL --right $ZfungiR,$InfungiR --CPU $num_threads
