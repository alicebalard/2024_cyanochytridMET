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


# Build transcript and gene expression matrices
$TRINITY_HOME/util/abundance_estimates_to_matrix.pl --est_method RSEM --gene_trans_map /scratch/erikamr/cyano_chytrid_met/data/combined_gene_trans_map_cds.txt --out_prefix RSEM --quant_files /scratch/erikamr/cyano_chytrid_met/data/quant_file_isoform.txt --name_sample_by_basedir /scratch/erikamr/cyano_chytrid_met/data/out_trinity_align_rsem
