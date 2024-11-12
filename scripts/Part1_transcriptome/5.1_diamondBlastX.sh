#!/bin/bash

#SBATCH --job-name=5.1_diamondblastx_e-5
#SBATCH --error=/home/alicebalard/Scripts/AliceScripts/cyanochytridMET/scripts/logs_dir/%x.%j.err
#SBATCH --output=/home/alicebalard/Scripts/AliceScripts/cyanochytridMET/scripts/logs_dir/%x.%j.out
#SBATCH --mail-user=alicebalard@zedat.fu-berlin.de  
#SBATCH --mail-type=end
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=500GB 
#SBATCH --time=2-24:00:00
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

cd /scratch/alicebalard/outData/diamondBlastX

num_threads=$SLURM_CPUS_PER_TASK

module load BLAST+/2.13.0-gompi-2022a
module load DIAMOND/2.0.13-GCC-11.2.0

## NCBI BLAST searches must be run with an appropriate output format in order to be imported dl database in nt repo: wget "ftp://ftp.ncbi.nlm.nih.gov/blast/db/nt.*.tar.gz" -P nt/ && for file in nt/*.tar.gz; do tar xf $file -C nt && rm $file; done

## echo "Make diamond db..."
## diamond makedb -p $num_threads --in /scratch/alicebalard/resources/uniprot/reference_proteomes.fasta.gz --taxonmap /scratch/alicebalard/resources/uniprot/reference_proteomes.taxid_map --taxonnodes /scratch/alicebalard/resources/taxdump/nodes.dmp --taxonnames /scratch/alicebalard/resources/taxdump/names.dmp --db /scratch/alicebalard/resources/uniprot/reference_proteomes.dmnd
## echo "diamond db done"

OUT=/scratch/alicebalard/outData/diamondBlastX
DB=/scratch/alicebalard/resources/uniprot/reference_proteomes.dmnd

###########################################
echo "First assembly: Z1 to Z12, only chytrids"
ASSEMBLY=/scratch/alicebalard/outData/assembly/trinity_out_dir/Trinity.fasta

echo "Create DIAMOND hits..."
## Run diamond
diamond blastx --query $ASSEMBLY \
	--db $DB \
        --outfmt 6 qseqid staxids bitscore qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sscinames sskingdoms skingdoms sphylums \
	--sensitive \
	--max-target-seqs 1 \
        --evalue 1e-5 \
        --threads $num_threads \
        > $OUT/assemblyZ_diamond_1e-5pval.out
echo "done."

echo "done."
	
##############################################################
echo "Second assembly: In1 to In12, chytrids infected by bacteria"
ASSEMBLY=/scratch/alicebalard/outData/assembly_In/trinity_out_dir/Trinity.fasta

echo "Create DIAMOND hits..."
## Run diamond
diamond blastx --query $ASSEMBLY \
	--db $DB \
        --outfmt 6 qseqid staxids bitscore qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sscinames sskingdoms skingdoms sphylums \
	--sensitive \
	--max-target-seqs 1 \
        --evalue 1e-5 \
        --threads $num_threads \
        > $OUT/assemblyIn_diamond_1e-5pval.out
echo "done."
