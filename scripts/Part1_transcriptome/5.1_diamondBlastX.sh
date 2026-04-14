#!/bin/bash

#SBATCH --job-name=5.1_diamondblastx_e-5
#SBATCH --error=/home/alicebalard/Scripts/AliceScripts/cyanochytridMET/scripts/logs_dir/%x.%j.err
#SBATCH --output=/home/alicebalard/Scripts/AliceScripts/cyanochytridMET/scripts/logs_dir/%x.%j.out
#SBATCH --mail-user=alicebalard@zedat.fu-berlin.de  
#SBATCH --mail-type=end
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
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
module load DIAMOND/2.1.8-GCC-12.3.0

echo "Download databases:"
wget -P /scratch/alicebalard/resources -c http://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.FULL.gz ## can be used zipped

wget -P /scratch/alicebalard/resources -c http://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.zip
gzip -d /scratch/alicebalard/resources/taxdmp.zip ## needs unzipped to access sub files

wget -P /scratch/alicebalard/resources -c https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz ## can be used zipped

echo "Make diamond db..."
diamond makedb -p $num_threads\
  --in /scratch/alicebalard/resources/nr.gz \
  --db /scratch/alicebalard/resources/nr \
  --taxonmap /scratch/alicebalard/resources/prot.accession2taxid.FULL.gz \
  --taxonnodes /scratch/alicebalard/resources/taxdmp/nodes.dmp \
  --taxonnames /scratch/alicebalard/resources/taxdmp/names.dmp
echo "diamond db done"

echo "Start diamond blastX search..."
OUT=/scratch/alicebalard/outData/diamondBlastX
DB=/scratch/alicebalard/resources/nr.dmnd

###########################################
echo "First assembly: Z1 to Z12, only chytrids"
ASSEMBLY=/scratch/alicebalard/outData/assemblies/assembly_Z/Trinity.fasta

echo "Create DIAMOND hits..."
## Run diamond
diamond blastx --query $ASSEMBLY \
	--db $DB \
        --outfmt 6 qseqid staxids bitscore qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sscinames sskingdoms skingdoms sphylums \
	--sensitive \
	--max-target-seqs 1 \
        --evalue 1e-5 \
        --threads $num_threads \
        > $OUT/assemblyZ_diamondNR_1e-5pval.out
echo "done."
	
##############################################################
echo "Second assembly: In1 to In12, chytrids infected by bacteria, & coculture"
ASSEMBLY=/scratch/alicebalard/outData/assemblies/assembly_In_coculture/Trinity.fasta

echo "Create DIAMOND hits..."
## Run diamond
diamond blastx --query $ASSEMBLY \
	--db $DB \
        --outfmt 6 qseqid staxids bitscore qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sscinames sskingdoms skingdoms sphylums \
	--sensitive \
	--max-target-seqs 1 \
        --evalue 1e-5 \
        --threads $num_threads \
        > $OUT/assemblyIn_cocult_diamondNR_1e-5pval.out
echo "done."
