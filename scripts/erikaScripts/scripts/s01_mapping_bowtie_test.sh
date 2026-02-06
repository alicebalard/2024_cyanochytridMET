#!/bin/bash

#SBATCH --job-name=erika_blast
#SBATCH --mail-user=erikamr@zedat.fu-berlin.de
#SBATCH --mail-type=end
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2GB
#SBATCH --time=1:00:00
#SBATCH --error=/scratch/erikamr/cyano_chytrid_met/scripts/%x.%j.err
#SBATCH --output=/scratch/erikamr/cyano_chytrid_met/scripts/%x.%j.out
#SBATCH --qos=standard
#SBATCH --partition=main,begendiv
#SBATCH --constraint=no_gpu

module add Trinity/2.10.0-foss-2019b-Python-3.7.4
module add Bowtie/1.3.1-GCC-11.3.0
module add RSEM/1.3.3-foss-2022a

cd /scratch/erikamr/cyano_chytrid_met

# Merge both transcriptomes
# ASSEMBLY_CHYTRID=/data/Trinity.filtered.Fungi.fasta
# ASSEMBLY_CYANO=/data/GCF_904830935.1_P._agardhii_No.976_rna_from_genomic.fna
# cat $ASSEMBLY_CHYTRID $ASSEMBLY_CYANO > /data/assembly_both.fna
ASSEMBLY_BOTH=/scratch/erikamr/cyano_chytrid_met/data/assembly_bothv2.fna


$TRINITY_HOME/util/align_and_estimate_abundance.pl --transcripts $ASSEMBLY_BOTH --seqType fq --samples_file /scratch/erikamr/cyano_chytrid_met/data/samples_file.txt --est_method RSEM --aln_method bowtie2 --bowtie2_RSEM "--no-unal -k 20 -q" --thread_count 1 --output_dir /data/out_trinity_align_rsem --gene_trans_map /scratch/erikamr/cyano_chytrid_met/data/combined_gene_trans_map.txt --prep_reference
