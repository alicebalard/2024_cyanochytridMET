#!/bin/sh
#SBATCH --job-name=6.2_checkQualityAssembly
#SBATCH --output=/home/alicebalard/Scripts/AliceScripts/cyanochytridMET/scripts/logs_dir/%x.%j.out
#SBATCH --error=/home/alicebalard/Scripts/AliceScripts/cyanochytridMET/scripts/logs_dir/%x.%j.err
#SBATCH --mail-user=alicebalard@zedat.fu-berlin.de
#SBATCH --mail-type=end
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=20GB 
#SBATCH --time=01:00:00
#SBATCH --qos=standard              
#SBATCH --partition=main

module purge
module add Trinity/2.15.2-foss-2023a

## Assembly Final

ASSEMBLY="/scratch/alicebalard/outData/assemblies/assemblyMergedFungi/Trinity.fasta"
echo "basic contig statistics"
$TRINITY_HOME/util/TrinityStats.pl $ASSEMBLY
