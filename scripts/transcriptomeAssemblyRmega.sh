#!/bin/bash

#SBATCH --job-name=de_novo_assembly                # replace name
#SBATCH --mail-user=alicebalard@zedat.fu-berlin.de  
#SBATCH --mail-type=end
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=60GB 
#SBATCH --time=03-24:00:00
#SBATCH --qos=standard              
#SBATCH --partition=main,begendiv
#SBATCH --constraint=no_gpu

module purge
module add Trinity/2.10.0-foss-2019b-Python-3.7.4

cd /scratch/alicebalard/outData/assembly/

## Store the names of the forward and reverse variables: we can use Z1 to Z6 samples but Z2 seems ideal 
## NB: merging sounds cool for completeness but could lead to chimeras
for i in $(seq 6); do var_name="F$i"; eval "$var_name='/scratch/alicebalard/RawData/Z'$i'_1_trimm_paired.fastq.gz'"; done
for i in $(seq 6); do var_name="R$i"; eval "$var_name='/scratch/alicebalard/RawData/Z'$i'_2_trimm_paired.fastq.gz'"; done

Trinity --seqType fq --normalize_by_read_set --max_memory 58G --left $F2 --right $R2 --CPU 20

## the output will be in 'trinity_out_dir/Trinity.fasta'
