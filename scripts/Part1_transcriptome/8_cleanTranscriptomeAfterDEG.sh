#!/bin/sh
#SBATCH --job-name=9_afterDEGdecont
#SBATCH --output=/home/alicebalard/Scripts/AliceScripts/cyanochytridMET/scripts/logs_dir/%x.%j.out
#SBATCH --error=/home/alicebalard/Scripts/AliceScripts/cyanochytridMET/scripts/logs_dir/%x.%j.err
#SBATCH --mail-user=alicebalard@zedat.fu-berlin.de
#SBATCH --mail-type=end
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=20GB 
#SBATCH --time=06:00:00
#SBATCH --qos=standard              
#SBATCH --partition=main

## Transcriptome clean:

ASSEMBLY=/scratch/alicebalard/outData/assemblyMergedFungi/trinity_out_dir/Trinity_eukaryoteHits.fasta

## Transcripts to remove after transcripts count (see R scripts)
# /home/alicebalard/Scripts/AliceScripts/cyanochytridMET/data/listOfTranscriptContaminant_toRmFromChytridTranscriptome
## rm first line's "x" if needed
# sed -i '1d' /home/alicebalard/Scripts/AliceScripts/cyanochytridMET/data/listOfTranscriptContaminant_toRmFromChytridTranscriptome

## Filter them out
awk 'BEGIN {
    while ((getline < "/home/alicebalard/Scripts/AliceScripts/cyanochytridMET/data/listOfTranscriptContaminant_toRmFromChytridTranscriptome") > 0) 
        l[$1] = 1 
} 
/^>/ { 
    header = $1; 
    skip = 0;  # Initialize skip flag
    for (id in l) {
        if (header ~ "^>" id "_i[0-9]*$" || header ~ "^>" id "$") {
            skip = 1;  # Set skip flag if match found
            break;     # Exit loop if match is found
        }
    } 
    if (skip) { 
        getline;  # Skip the next line (the sequence)
        next;     # Skip this header
    }
} 
{ print }' "$ASSEMBLY" > "$ASSEMBLY.rmDEGconta.fasta"

########################
## Check assembly stats
module purge
module add Trinity/2.10.0-foss-2019b-Python-3.7.4

ASSEMBLY2="$ASSEMBLY.rmDEGconta.fasta"
echo "basic contig statistics for Trinity_eukaryoteHits.fasta:"
$TRINITY_HOME/util/TrinityStats.pl $ASSEMBLY2

##############
## Check BUSCO
module purge
source ~/.bashrc
module add BUSCO/5.1.2-foss-2020b

cd /scratch/alicebalard/outData/assemblyMergedFungi

echo "Eukaryote blastx hits"

busco -i $ASSEMBLY2 -l fungi_odb10 -m transcriptome -c 10 -o "BUSCO_MergedFungiTranscriptome_eukaryoteHits.rmDEGconta"
## l = current fungi database (will be downloaded automatically)  
