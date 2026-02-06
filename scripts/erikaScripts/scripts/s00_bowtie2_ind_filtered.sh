#!/bin/bash

#SBATCH --job-name=erika_blast
#SBATCH --mail-user=erikamr@zedat.fu-berlin.de
#SBATCH --mail-type=end
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=2GB
#SBATCH --time=24:00:00
#SBATCH --error=/scratch/erikamr/cyano_chytrid_met/scripts/%x.%j.err
#SBATCH --output=/scratch/erikamr/cyano_chytrid_met/scripts/%x.%j.out
#SBATCH --qos=standard
#SBATCH --partition=main,begendiv
#SBATCH --constraint=no_gpu

# Enable autoswapping of same-name modules
export LMOD_DISABLE_SAME_NAME_AUTOSWAP=no

module add Bowtie2/2.4.5-GCC-11.3.0
module add SAMtools/1.16.1-GCC-11.3.0

BOWTIE2_INDEX="/scratch/erikamr/cyano_chytrid_met/data/assembly_both_cds.fna"
INPUT_DIR="/scratch/alicebalard/outData/sortmerna"
OUTPUT_DIR="/scratch/erikamr/cyano_chytrid_met/data/aligned_bowtie2_ind"

# Find all forward read files (_non_rRNA_fwd.fq.gz) in the input directory
for FWD_READS in ${INPUT_DIR}/*_non_rRNA_fwd.fq.gz
do
    # Extract the base sample name by removing the forward read suffix
    SAMPLE_BASE=$(basename ${FWD_READS} _non_rRNA_fwd.fq.gz)

    # Define the corresponding reverse read file based on the sample base name
    REV_READS="${INPUT_DIR}/${SAMPLE_BASE}_non_rRNA_rev.fq.gz"

    # Define output filenames for alignment stats and BAM files
    ALIGN_STATS="${OUTPUT_DIR}/${SAMPLE_BASE}_align_stats.txt"
    BAM_FILE="${OUTPUT_DIR}/${SAMPLE_BASE}_bowtie2.bam"

    # Check if the reverse read file exists
    if [ -f "${REV_READS}" ]; then
        echo "Processing sample ${SAMPLE_BASE}..."

        # Run Bowtie2 and pipe the output to Samtools for BAM conversion
        bowtie2 -p 16 -q --no-unal -k 20 -x ${BOWTIE2_INDEX} -1 ${FWD_READS} -2 ${REV_READS} 2>${ALIGN_STATS} | samtools view -@16 -Sb -o ${BAM_FILE}

        echo "Finished sample ${SAMPLE_BASE}"
    else
        echo "Warning: Reverse read file for sample ${SAMPLE_BASE} not found! Skipping..."
    fi
done
