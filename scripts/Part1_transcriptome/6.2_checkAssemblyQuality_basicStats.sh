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
module add Trinity/2.10.0-foss-2019b-Python-3.7.4

## Assembly Final

ASSEMBLY="/scratch/alicebalard/outData/assemblyMergedFungi/trinity_out_dir/Trinity.fasta"
echo "basic contig statistics"
$TRINITY_HOME/util/TrinityStats.pl $ASSEMBLY
