#!/bin/bash

#SBATCH --job-name=de_novo_assembly
#SBATCH --error=/home/alicebalard/Scripts/AliceScripts/cyanochytridMET/scripts/logs_dir/%x.%j.err
#SBATCH --output=/home/alicebalard/Scripts/AliceScripts/cyanochytridMET/scripts/logs_dir/%x.%j.out
#SBATCH --mail-user=alicebalard@zedat.fu-berlin.de  
#SBATCH --mail-type=end
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=250GB 
#SBATCH --time=7-24:00:00
#SBATCH --qos=standard              
#SBATCH --partition=main,begendiv
#SBATCH --constraint=no_gpu

module purge
module add Trinity/2.10.0-foss-2019b-Python-3.7.4

cd /scratch/alicebalard/outData/assembly/
## the output will be in 'trinity_out_dir/Trinity.fasta'

## Store the names of the forward and reverse variables: we can use Z1 to Z12 samples 
for i in $(seq 12); do var_name="F$i"; eval "$var_name='/scratch/alicebalard/outData/sortmerna/Z'$i'_non_rRNA_fwd.fq.gz'"; done
for i in $(seq 12); do var_name="R$i"; eval "$var_name='/scratch/alicebalard/outData/sortmerna/Z'$i'_non_rRNA_rev.fq.gz'"; done

zcat $F1 $F2 $F3 $F4 $F5 $F6 $F7 $F8 $F9 $F10 $F11 $F12 > combined_left.fq
zcat $R1 $R2 $R3 $R4 $R5 $R6 $R7 $R8 $R9 $R10 $R11 $R12 > combined_right.fq

num_threads=$SLURM_CPUS_PER_TASK
echo "Running with $num_threads threads"
MAX_MEM=250G

Trinity --seqType fq --normalize_by_read_set --max_memory $MAX_MEM --left combined_left.fq --right combined_right.fq --CPU $num_threads
