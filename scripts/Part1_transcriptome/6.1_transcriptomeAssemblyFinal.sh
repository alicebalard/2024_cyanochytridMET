#!/bin/bash

#SBATCH --job-name=6_de_novo_assembly_final
#SBATCH --error=/home/alicebalard/Scripts/AliceScripts/cyanochytridMET/scripts/logs_dir/%x.%j.err
#SBATCH --output=/home/alicebalard/Scripts/AliceScripts/cyanochytridMET/scripts/logs_dir/%x.%j.out
#SBATCH --mail-user=alicebalard@zedat.fu-berlin.de
#SBATCH --mail-type=end
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=220GB
#SBATCH --time=10-24:00:00
#SBATCH --qos=standard
#SBATCH --partition=main,begendiv
#SBATCH --constraint=no_gpu

module purge
module add Trinity/2.15.2-foss-2023a

num_threads=$SLURM_CPUS_PER_TASK

DIR=/scratch/alicebalard/outData/assemblies/assemblyMergedFungi

mkdir -p $DIR; cd $DIR
## the output will be in 'trinity_out_dir/Trinity.fasta'

## Reads to be assembled:
ZfungiL=assemblyZ.reads_filteredFungi_left.fq
ZfungiR=assemblyZ.reads_filteredFungi_right.fq
InfungiL=assemblyIn.reads_filteredFungi_left.fq
InfungiR=assemblyIn.reads_filteredFungi_right.fq

echo "Number of left reads assembled from Z and In transcriptomes:"
#cat $ZfungiL | echo $((`wc -l`/4))
#cat $InfungiL | echo $((`wc -l`/4))

echo "Number of right reads assembled from Z and In transcriptomes:"
#cat $ZfungiR | echo $((`wc -l`/4))
#cat $InfungiR | echo $((`wc -l`/4))

num_threads=$SLURM_CPUS_PER_TASK
echo "Transcriptome assembly running with $num_threads threads"
MAX_MEM=200G   # Trinity ≥2.14 caps to ~200G internally, don’t request more here

Trinity --FORCE --seqType fq --normalize_by_read_set --max_memory $MAX_MEM --left $ZfungiL,$InfungiL --right $ZfungiR,$InfungiR --CPU $num_threads
## Retry with --FORCE at the end if enters a neverending loop!
