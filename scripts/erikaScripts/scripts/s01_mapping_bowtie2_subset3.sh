#!/bin/bash

#SBATCH --job-name=erika_blast
#SBATCH --mail-user=erikamr@zedat.fu-berlin.de
#SBATCH --mail-type=end
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --mem-per-cpu=25GB
#SBATCH --time=12:00:00
#SBATCH --error=/scratch/erikamr/cyano_chytrid_met/scripts/%x.%j.err
#SBATCH --output=/scratch/erikamr/cyano_chytrid_met/scripts/%x.%j.out
#SBATCH --qos=standard
#SBATCH --partition=main,begendiv
#SBATCH --constraint=no_gpu

# Enable autoswapping of same-name modules
export LMOD_DISABLE_SAME_NAME_AUTOSWAP=no

module add Trinity/2.10.0-foss-2019b-Python-3.7.4
module add Bowtie2/2.4.5-GCC-11.3.0
module add RSEM/1.3.3-foss-2019b
module add SAMtools/1.16.1-GCC-11.3.0

#cd /scratch/erikamr/cyano_chytrid_met

ASSEMBLY_BOTH=/scratch/erikamr/cyano_chytrid_met/data/assembly_both_cds.fna
OUTPUT_DIR=/scratch/erikamr/cyano_chytrid_met/data/out_trinity_align_rsem_concatenated_bowtie2_subset3

$TRINITY_HOME/util/align_and_estimate_abundance.pl --transcripts $ASSEMBLY_BOTH --seqType fq --samples_file /scratch/erikamr/cyano_chytrid_met/data/samples_file_subset3.txt --est_method RSEM --aln_method bowtie2 --bowtie2_RSEM "--no-unal -k 4" --thread_count 5 --output_dir $OUTPUT_DIR --gene_trans_map /scratch/erikamr/cyano_chytrid_met/data/filtered_gene_trans_map_cds.txt --prep_reference 
