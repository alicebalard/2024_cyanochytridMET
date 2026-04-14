#!/bin/sh
#SBATCH --job-name=checkQualityAssembly
#SBATCH --output=/home/alicebalard/Scripts/AliceScripts/cyanochytridMET/scripts/logs_dir/%x.%j.out
#SBATCH --error=/home/alicebalard/Scripts/AliceScripts/cyanochytridMET/scripts/logs_dir/%x.%j.err
#SBATCH --mail-user=alicebalard@zedat.fu-berlin.de
#SBATCH --mail-type=end
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=20GB 
#SBATCH --time=01:00:00
#SBATCH --qos=standard              
#SBATCH --partition=main

module purge
module add Trinity/2.10.0-foss-2019b-Python-3.7.4

echo "Z1Z12:"
ASSEMBLY="/scratch/alicebalard/outData/assemblies/assembly_Z/Trinity.fasta"
echo "basic contig statistics"
$TRINITY_HOME/util/TrinityStats.pl $ASSEMBLY


echo "In1In12 + coculture:"
ASSEMBLY="/scratch/alicebalard/outData/assemblies/assembly_In_coculture/Trinity.fasta"
echo "basic contig statistics"
$TRINITY_HOME/util/TrinityStats.pl $ASSEMBLY

echo "In1In12 only (previous):"
ASSEMBLY="/scratch/alicebalard/outData/assemblies/assembly_In/Trinity.fasta"
echo "basic contig statistics"
$TRINITY_HOME/util/TrinityStats.pl $ASSEMBLY

assembliesEarlyBAK/assembly_In/
