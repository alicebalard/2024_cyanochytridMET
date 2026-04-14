#!/bin/bash

#SBATCH --job-name=4.2_de_novo_assembly_In_coculture
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
MAX_MEM=200G   # Trinity ≥2.14 caps to ~200G internally, don’t request more here

DIR=/scratch/alicebalard/outData/assemblies/assembly_In_coculture

mkdir -p $DIR; cd $DIR
## the output will be 'trinity_out_dir/Trinity.fasta' and 'Trinity.fasta.gene_trans_map'

if [ ! -f $DIR/combined_right.fq ]; then
    ## Store the names of the forward and reverse variables: we can use Z1 to Z12 samples
    for i in $(seq 12); do var_name="F$i"; eval "$var_name='/scratch/alicebalard/outData/sortmerna/In'$i'_non_rRNA_fwd.fq.gz'"; done
    for i in $(seq 12); do var_name="R$i"; eval "$var_name='/scratch/alicebalard/outData/sortmerna/In'$i'_non_rRNA_rev.fq.gz'"; done

    ## NEW 2026: add new data from Jurgen
    F13=/scratch/strassert/chytrids_e/02_rRNA_removal/CyanoP_non_rRNA_fwd.fq.gz
    R13=/scratch/strassert/chytrids_e/02_rRNA_removal/CyanoP_non_rRNA_rev.fq.gz

    zcat $F1 $F2 $F3 $F4 $F5 $F6 $F7 $F8 $F9 $F10 $F11 $F12 $F13 > combined_left.fq
    zcat $R1 $R2 $R3 $R4 $R5 $R6 $R7 $R8 $R9 $R10 $R11 $R12 $R13 > combined_right.fq
fi

echo "Running assembly with $num_threads threads"

#Trinity --seqType fq --normalize_by_read_set --max_memory $MAX_MEM --left combined_left.fq --right combined_right.fq --CPU $num_threads 

## enter a loop, I manually remove 3 broken parts at the end then relaunch:
## (base) [alicebalard@login assembly_In_coculture]$ rm -rf trinity_out_dir/read_partitions/Fb_3/CBin_3916/c397550.trinity.reads.fa.out
## (base) [alicebalard@login assembly_In_coculture]$ rm -rf trinity_out_dir/read_partitions/Fb_1/CBin_1129/c113589.trinity.reads.fa.out
## (base) [alicebalard@login assembly_In_coculture]$ rm -rf trinity_out_dir/read_partitions/Fb_0/CBin_751/c75285.trinity.reads.fa.out
## (base) [alicebalard@login assembly_In_coculture]$ rm -rf trinity_out_dir/read_partitions/Fb_3/CBin_3908/c396796.trinity.reads.fa.out

## Retry with --FORCE at the end if enters a neverending loop!
Trinity --FORCE --seqType fq --normalize_by_read_set --max_memory $MAX_MEM --left combined_left.fq --right combined_right.fq --CPU $num_threads 




echo "Assembly finished for the infected filaments!"
