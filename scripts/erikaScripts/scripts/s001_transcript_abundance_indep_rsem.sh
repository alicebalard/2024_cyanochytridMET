#!/bin/bash

#SBATCH --job-name=erika_blast
#SBATCH --mail-user=erikamr@zedat.fu-berlin.de
#SBATCH --mail-type=end
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8 
#SBATCH --mem-per-cpu=4GB 
#SBATCH --time=15:00:00
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

#cd /scratch/erikamr/cyano_chytrid_met

ASSEMBLY_BOTH=/scratch/erikamr/cyano_chytrid_met/data/assembly_both_cds.fna
OUTPUT_DIR=/scratch/erikamr/cyano_chytrid_met/data/aligned_bowtie2_ind

$TRINITY_HOME/util/align_and_estimate_abundance.pl --transcripts $ASSEMBLY_BOTH --seqType fq --samples_file /scratch/erikamr/cyano_chytrid_met/data/samples_file_bam.txt --est_method RSEM --aln_method bam --thread_count 8 --output_dir $OUTPUT_DIR --gene_trans_map /scratch/erikamr/cyano_chytrid_met/data/filtered_gene_trans_map_cds.txt --prep_reference
