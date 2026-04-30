#!/bin/bash

#SBATCH --job-name=7.4_diamondblastp
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

num_threads=$SLURM_CPUS_PER_TASK

module load BLAST+/2.13.0-gompi-2022a
module load DIAMOND/2.1.8-GCC-12.3.0

##echo "Download databases:"
##wget -P /scratch/alicebalard/resources -c http://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.FULL.gz ## can be used zipped
##
##wget -P /scratch/alicebalard/resources -c http://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.zip
##unzip -d taxdmp taxdmp.zip; rm taxdmp.zip ## needs unzipped to access sub files

## Manual fix: NCBI's newer taxonomy ranks need editing for DIAMOND 2.1.8 compatibility (will be fixed in later diamond version)
# Fix problematic ranks
# Simple string replacement (no regex escaping needed)
# sed -i 's/acellular root/superkingdom/g' nodes.dmp
# sed -i 's/cellular root/superkingdom/g' nodes.dmp
# sed -i 's/ domain / superkingdom /g' nodes.dmp
# sed -i 's/ realm / superkingdom /g' nodes.dmp

##wget -P /scratch/alicebalard/resources -c https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz ## can be used zipped

## echo "Make diamond db..."
## diamond makedb -p $num_threads --in /scratch/alicebalard/resources/nr.gz --db /scratch/alicebalard/resources/nr --taxonmap /scratch/alicebalard/resources/prot.accession2taxid.FULL.gz --taxonnodes /scratch/alicebalard/resources/taxdmp/nodes.dmp --taxonnames /scratch/alicebalard/resources/taxdmp/names.dmp
## echo "diamond db done"

cd /scratch/alicebalard/outData/assemblies/assemblyMergedFungi/annotation
mkdir -p blastp
cd blastp

echo "Start diamond blastP search..."
OUT=/scratch/alicebalard/outData/assemblies/assemblyMergedFungi/annotation/blastp
DB=/scratch/alicebalard/resources/nr.dmnd

## blastp: protein query (TransDecoder peptides) against nr
coding_seqs=/scratch/alicebalard/outData/assemblies/assemblyMergedFungi/annotation/transdecoder/Trinity.fasta.transdecoder.pep

echo "Create DIAMOND hits..."
## Run diamond
diamond blastp -q $coding_seqs \
	--db $DB \
        --outfmt 6 \
	--sensitive \
	--max-target-seqs 10 \
        --max-hsps 1 \
        --evalue 1e-5 \
        --threads $num_threads \
	-o $OUT/diamond_nr.outfmt6
echo "done."
