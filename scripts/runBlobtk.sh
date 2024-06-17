#!/bin/bash

#SBATCH --job-name=blobtk
#SBATCH --error=/home/alicebalard/Scripts/AliceScripts/cyanochytridMET/scripts/logs_dir/%x.%j.err
#SBATCH --output=/home/alicebalard/Scripts/AliceScripts/cyanochytridMET/scripts/logs_dir/%x.%j.out
#SBATCH --mail-user=alicebalard@zedat.fu-berlin.de  
#SBATCH --mail-type=end
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-core=8GB 
#SBATCH --time=7-24:00:00
#SBATCH --qos=standard              
#SBATCH --partition=main,begendiv
#SBATCH --constraint=no_gpu

module purge
module load Anaconda3

source ~/.bashrc
conda init --all
conda activate btk

cd /scratch/alicebalard/outData/blobtools

## https://blobtoolkit.genomehubs.org/blobtools2/blobtools2-tutorials/getting-started-with-blobtools2/#open_dataset

ASSEMBLY=/scratch/alicebalard/outData/assembly/trinity_out_dir/Trinity.fasta 
READS1=/scratch/alicebalard/outData/assembly/combined_left.fq
READS2=/scratch/alicebalard/outData/assembly/combined_right.fq
num_threads=$SLURM_CPUS_PER_TASK

echo "Create a Blobdir with a fasta assembly..."

blobtools create --fasta $ASSEMBLY Z1Z12assembly

# creates a new directory in the location specified by the last argument that contains a set of files containing
# values for GC-content (gc.json), length (length.json), number of Ns (ncount.json) and sequence names (identifiers.json)
# for each sequence in the assembly. A final file (meta.json) contains metadata for the dataset describing the datatypes
# of the available fields and the ranges of values for each of these fields

echo "Add BLAST hits..."
# The BlobTools approach uses BLAST hits to provide taxonomic annotation for each sequence in an assembly. NCBI BLAST 
# searches must be run with an appropriate output format in order to be imported
## dl database in nt repo:
## wget "ftp://ftp.ncbi.nlm.nih.gov/blast/db/nt.*.tar.gz" -P nt/ && for file in nt/*.tar.gz; do tar xf $file -C nt && rm $file; done

module load BLAST+/2.13.0-gompi-2022a

blastn -db nt \
       -query $ASSEMBLY \
       -outfmt "6 qseqid staxids bitscore std" \
       -max_target_seqs 10 \
       -max_hsps 1 \
       -evalue 1e-25 \
       -num_threads $num_threads \
       -out blast.out

# This file can be imported to assign taxonomic labels to the assembly scaffolds using taxrules:

blobtools add \
    --hits blast.out \
    --taxrule bestsumorder \
    --taxdump ~/taxdump \
    Z1Z2assembly

echo "Add mapping coverage..."
# The final data required to generate a standard blob plot is coverage information from mapping sequencing reads back 
# to the assembly. BlobTools2 parses sorted SAM/BAM format files to calculate coverage information using the PySAM library.
# To map Illumina reads using minimap2:

module load minimap2/2.24-GCCcore-11.3.0

minimap2 -ax sr \
         -t $num_threads $ASSEMBLY \
         READS1 READS2 \
| samtools sort -@$num_threads -O BAM -o assembly.reads.bam -

# The resulting BAM file(s) can then be imported into a BlobDir:

./blobtools2/blobtools add \
    --cov assembly.reads.bam \
    Z1Z2assembly

echo "Add BUSCO scores..."
module add BUSCO/5.1.2-foss-2020b

## With more samples
busco -i $ASSEMBLY -l fungi_odb10 -o output_directory -m transcriptome -c $num_threads -o "BUSCO_transcZ1to12"

## l = current fungi database (will be downloaded automatically)
## m = genome or transcriptome (in your case transcriptome)

# Add the BUSCO full_table output to a BlobDir:

blobtools add \
    --busco BUSCO_transcZ1to12/run_fungi_odb10/full_table.tsv \
    $ASSEMBLY

echo "END BLOBTK"
