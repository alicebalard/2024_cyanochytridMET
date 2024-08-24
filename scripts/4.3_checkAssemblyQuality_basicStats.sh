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

## Z1Z12

ASSEMBLY="/scratch/alicebalard/outData/assembly/trinity_out_dir/Trinity.fasta"
echo "basic contig statistics"
$TRINITY_HOME/util/TrinityStats.pl $ASSEMBLY


## In1In12

ASSEMBLY="/scratch/alicebalard/outData/assembly_In/trinity_out_dir/Trinity.fasta"
echo "basic contig statistics"
$TRINITY_HOME/util/TrinityStats.pl $ASSEMBLY
