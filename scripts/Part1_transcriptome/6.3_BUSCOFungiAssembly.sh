#!/bin/sh --login
#SBATCH --job-name=6.3_BUSCO
#SBATCH --output=/home/alicebalard/Scripts/AliceScripts/cyanochytridMET/scripts/logs_dir/%x.%j.out
#SBATCH --error=/home/alicebalard/Scripts/AliceScripts/cyanochytridMET/scripts/logs_dir/%x.%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=4GB 
#SBATCH --time=24:00:00
#SBATCH --partition=main,begendiv
#SBATCH --constraint=no_gpu
#SBATCH --qos=standard
#SBATCH --mail-user=alicebalard@zedat.fu-berlin.de
#SBATCH --mail-type=end

## Enable same name autoswapping
LMOD_DISABLE_SAME_NAME_AUTOSWAP=no

module purge
source ~/.bashrc
module add BUSCO/5.8.2-foss-2023a

cd /scratch/alicebalard/outData/assemblies

## l = current fungi database (will be downloaded automatically)

echo "run BUSCO assembly_Z..."
#busco -i assembly_Z/Trinity.fasta -l fungi_odb10 -m transcriptome -c 10 -o "BUSCO_assembly_Z"

echo "run BUSCO assembly_In 2024..."
#busco -i assembly_In/Trinity.fasta -l fungi_odb10 -m transcriptome -c 10 -o "BUSCO_assembly_In_2024"

echo "run BUSCO assembly_In new with coculture..."
#busco -i assembly_In_coculture/Trinity.fasta -l fungi_odb10 -m transcriptome -c 10 -o "BUSCO_assembly_In_coculture"

echo "run BUSCO assembly_merged 2024 (Z + In)..."
busco -i assemblyMergedFungi/BAK/trinity_out_dir/Trinity.fasta -l fungi_odb10 -m transcriptome -c 10 -o "BUSCO_assembly_merged_fungi_2024"

echo "run BUSCO assembly_merged (Z + In_coculture)..."
busco -i assemblyMergedFungi/Trinity.fasta -l fungi_odb10 -m transcriptome -c 10 -o "BUSCO_assembly_merged_fungi"


