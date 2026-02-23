#!/bin/bash

#SBATCH --job-name=erika_blast
#SBATCH --mail-user=erikamr@zedat.fu-berlin.de
#SBATCH --mail-type=end
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=2GB
#SBATCH --time=24:00:00
#SBATCH --error=/scratch/erikamr/cyano_chytrid_met/scripts/%x.%j.err
#SBATCH --output=/scratch/erikamr/cyano_chytrid_met/scripts/%x.%j.out
#SBATCH --qos=standard
#SBATCH --partition=main,begendiv
#SBATCH --constraint=no_gpu

# Enable autoswapping of same-name modules
export LMOD_DISABLE_SAME_NAME_AUTOSWAP=no

module add Trinity/2.10.0-foss-2019b-Python-3.7.4
module add Bowtie/1.3.1-GCC-11.3.0
module add RSEM/1.3.3-foss-2022a

cd /scratch/erikamr/cyano_chytrid_met

## prepare in part2 S01:
ASSEMBLY_BOTH=/scratch/erikamr/cyano_chytrid_met/data/assembly_both_cds_final_hope.fna
GTM=/scratch/erikamr/cyano_chytrid_met/data/combined_gene_trans_map_cds_final_hope.txt

## manually prepared:
SAMPLE_FILE=/scratch/erikamr/cyano_chytrid_met/data/samples_file_remove.txt

OUTDIR=/scratch/erikamr/cyano_chytrid_met/data/out_trinity_align_rsem_final_hope

$TRINITY_HOME/util/align_and_estimate_abundance.pl --transcripts $ASSEMBLY_BOTH --seqType fq --samples_file $SAMPLE_FILE --est_method RSEM --aln_method bowtie --bowtie_RSEM "--no-unal -k 20" --thread_count 8 --output_dir $OUTDIR --gene_trans_map $GTM --prep_reference
