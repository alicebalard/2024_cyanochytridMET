#!/bin/bash

#SBATCH --job-name=4.1_de_novo_assembly
#SBATCH --error=/home/alicebalard/Scripts/AliceScripts/cyanochytridMET/scripts/logs_dir/%x.%j.err
#SBATCH --output=/home/alicebalard/Scripts/AliceScripts/cyanochytridMET/scripts/logs_dir/%x.%j.out
#SBATCH --mail-user=alicebalard@zedat.fu-berlin.de
#SBATCH --mail-type=end
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=220GB ## was 250
#SBATCH --time=7-24:00:00
#SBATCH --qos=standard
#SBATCH --partition=main,begendiv
#SBATCH --constraint=no_gpu

module purge
module add Trinity/2.15.2-foss-2023a

DIR=/scratch/alicebalard/outData/assembly_Z

mkdir -p $DIR; cd $DIR
## the output will be 'trinity_out_dir/Trinity.fasta' and 'Trinity.fasta.gene_trans_map'

if [ ! -f $DIR/combined_right.fq ]; then

    ## Store the names of the forward and reverse variables: we can use Z1 to Z12 samples 
    for i in $(seq 12); do var_name="F$i"; eval "$var_name='/scratch/alicebalard/outData/sortmerna/Z'$i'_non_rRNA_fwd.fq.gz'"; done
    for i in $(seq 12); do var_name="R$i"; eval "$var_name='/scratch/alicebalard/outData/sortmerna/Z'$i'_non_rRNA_rev.fq.gz'"; done

    zcat $F1 $F2 $F3 $F4 $F5 $F6 $F7 $F8 $F9 $F10 $F11 $F12 > $DIR/combined_left.fq
    zcat $R1 $R2 $R3 $R4 $R5 $R6 $R7 $R8 $R9 $R10 $R11 $R12 > $DIR/combined_right.fq

fi

num_threads=$SLURM_CPUS_PER_TASK
echo "Running with $num_threads threads"
MAX_MEM=200G # was 240

Trinity --seqType fq --normalize_by_read_set --max_memory $MAX_MEM --left $DIR/combined_left.fq --right $DIR/combined_right.fq --CPU $num_threads

## --normalize_by_read_set     → normalizes per library (good for mixed samples)
