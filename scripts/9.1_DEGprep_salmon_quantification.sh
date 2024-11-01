#!/bin/bash

#SBATCH --job-name=9.1_DEGSalmon_cyano_In
#SBATCH --error=/home/alicebalard/Scripts/AliceScripts/cyanochytridMET/scripts/logs_dir/%x.%j.err
#SBATCH --output=/home/alicebalard/Scripts/AliceScripts/cyanochytridMET/scripts/logs_dir/%x.%j.out
#SBATCH --mail-user=alicebalard@zedat.fu-berlin.de
#SBATCH --mail-type=end
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=5GB
#SBATCH --time=1-24:00:00
#SBATCH --qos=standard
#SBATCH --partition=main,begendiv
#SBATCH --constraint=no_gpu

## https://combine-lab.github.io/salmon/getting_started/

module load Salmon/1.9.0-GCC-11.3.0

cd /scratch/alicebalard/outData/DEGSalmon_Alice

## load cyanobacteria transcriptome
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/904/830/935/GCF_904830935.1_P._agardhii_No.976/GCF_904830935.1_P._agardhii_No.976_rna_from_genomic.fna.gz

## index transcriptome
salmon index -t GCF_904830935.1_P._agardhii_No.976_cds_from_genomic.fna.gz -i GCF_904830935.1_P._agardhii_No.976_cds_from_genomic.fna.gz_index

INDEX="GCF_904830935.1_P._agardhii_No.976_cds_from_genomic.fna.gz_index"

## Quantifying the samples

# Define the input directory and output directory
input_dir="/scratch/alicebalard/outData/sortmerna"
output_dir="/scratch/alicebalard/outData/DEGSalmon_Alice/quants"

# Create output directory if it doesn't exist
mkdir -p ${output_dir}

# Loop through each sample pair
for fwd in ${input_dir}/In*_non_rRNA_fwd.fq.gz; do
    # Get the sample name by removing the "_fwd.fq.gz" part
    samp=$(basename ${fwd} _non_rRNA_fwd.fq.gz)
    
    # Define the reverse read file path
    rev="${input_dir}/${samp}_non_rRNA_rev.fq.gz"
    
    # Check if the reverse file exists
    if [[ -f ${rev} ]]; then
        echo "Processing sample ${samp}"
        salmon quant -i $INDEX -l A \
            -1 ${fwd} \
            -2 ${rev} \
            -p 8 --validateMappings --gcBias -o ${output_dir}/${samp}_quant
    else
        echo "Warning: Reverse read file for sample ${samp} not found: ${rev}"
    fi
done

## NB: We recommend using the --gcBias flag which estimates a correction factor for systematic biases commonly present in RNA-seq data (Love, Hogenesch, and Irizarry 2016; Patro et al. 2017), unless you are certain that your data do not contain such bias.

cd quants
for i in In*; do echo $i; cat $i/quant.sf | cut -f 5 | grep -v -e "0.000" | wc -l ; done 
