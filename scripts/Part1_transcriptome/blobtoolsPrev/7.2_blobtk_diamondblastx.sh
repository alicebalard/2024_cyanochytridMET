#!/bin/bash

#SBATCH --job-name=7.2_blobtk_diamondblastx_Final_e-5
#SBATCH --error=/home/alicebalard/Scripts/AliceScripts/cyanochytridMET/scripts/logs_dir/%x.%j.err
#SBATCH --output=/home/alicebalard/Scripts/AliceScripts/cyanochytridMET/scripts/logs_dir/%x.%j.out
#SBATCH --mail-user=alicebalard@zedat.fu-berlin.de  
#SBATCH --mail-type=end
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=200GB 
#SBATCH --time=1-24:00:00
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

cd /scratch/alicebalard/outData/blobtools

## https://blobtoolkit.genomehubs.org/blobtools2/blobtools2-tutorials/getting-started-with-blobtools2/#open_dataset

num_threads=$SLURM_CPUS_PER_TASK

module load BLAST+/2.13.0-gompi-2022a
module load DIAMOND/2.0.13-GCC-11.2.0

## The BlobTools approach uses Diamond hits to provide taxonomic annotation for each sequence in an assembly. NCBI BLAST searches must be run with an appropriate output format in order to be imported dl database in nt repo: wget "ftp://ftp.ncbi.nlm.nih.gov/blast/db/nt.*.tar.gz" -P nt/ && for file in nt/*.tar.gz; do tar xf $file -C nt && rm $file; done

## Make diamond db NB LOOOONG PROCESS!!
## diamond makedb -p $num_threads --in uniprot/reference_proteomes.fasta.gz --taxonmap uniprot/reference_proteomes.taxid_map --taxonnodes taxdump/nodes.dmp -d uniprot/reference_proteomes.dmnd

## Final assembly
ASSEMBLY=/scratch/alicebalard/outData/assemblyFinal/trinity_out_dir/Trinity.fasta 
BTKOUT=/scratch/alicebalard/outData/blobtools/FINALTran

echo "Create DIAMOND hits..."
## Run diamond: a few days runtime
diamond blastx --query $ASSEMBLY \
        --db /scratch/alicebalard/outData/blobtools/uniprot/reference_proteomes.dmnd \
        --outfmt 6 qseqid staxids bitscore qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore \
        --sensitive \
        --max-target-seqs 1 \
        --evalue 1e-5 \
        --threads $num_threads \
        > $BTKOUT/diamond_1e-5pval.out

# This file can be imported to assign taxonomic labels to the assembly scaffolds using taxrules:
echo "Add Diamond BLASTx hits to blobtools dir..."
blobtools add --hits $BTKOUT/diamond_1e-5pval.out \
	  --taxrule bestsumorder \
	  --taxdump /scratch/alicebalard/outData/blobtools/taxdump \
	  $BTKOUT
echo "Diamond BLASTx hits added to blobtools dir"
