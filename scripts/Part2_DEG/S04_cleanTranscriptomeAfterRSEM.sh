#!/bin/bash
# S04_decontFinal
# Step 1: Find column indices of C and In samples from header
# Step 2: For TRINITY rows, check if any of those columns > 0
# Step 3: Output contaminating gene IDs
# Step 4: Filter them from the fasta

#SBATCH --job-name=S04
#SBATCH --mail-user=alicebalard@zedat.fu-berlin.de
#SBATCH --mail-type=end,fail
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=4GB 
#SBATCH --time=24:00:00
#SBATCH --output=/home/alicebalard/Scripts/AliceScripts/cyanochytridMET/scripts/logs_dir/%x.%j.out
#SBATCH --error=/home/alicebalard/Scripts/AliceScripts/cyanochytridMET/scripts/logs_dir/%x.%j.err
#SBATCH --qos=standard
#SBATCH --partition=main,begendiv
#SBATCH --constraint=no_gpu

set -euo pipefail                                                                                                                                                                       
export LMOD_DISABLE_SAME_NAME_AUTOSWAP=no

ASSEMBLYDIR=/scratch/alicebalard/outData/assemblies/assemblyMergedFungi

cd $ASSEMBLYDIR

MATRIX=/home/alicebalard/Scripts/AliceScripts/cyanochytridMET/data/RSEM_new_hope.gene.counts.matrix
ASSEMBLY=$ASSEMBLYDIR/Trinity_eukaryoteHits.fasta

# Get column indices of samples starting with C or In (control_cyano_C*, met_cyano_C*, control_both_In*, met_both_In*)
awk '
NR==1 {
    for (i=2; i<=NF; i++) {
        # Extract sample suffix after last underscore: C* or In*
        split($i, a, "_")
        suffix = a[length(a)]
        if (suffix ~ /^C/ || suffix ~ /^In/) {
            col[i] = 1
        }
    }
    next
}
# Only TRINITY rows (chytrid assembly)
/^TRINITY/ {
    for (i in col) {
        if ($i > 0) {
            print $1
            next
        }
    }
}
' "$MATRIX" > $ASSEMBLYDIR/contaminant_genes.txt

echo "Number of contaminant genes found:"
wc -l $ASSEMBLYDIR/contaminant_genes.txt

# use the awk approach but simplified to remove contaminants:
awk 'NR==FNR {l[$1]=1; next}
     /^>/ {
         id = substr($1, 2)          # remove ">"
         # match gene ID (strip _i* suffix if present)
         gsub(/_i[0-9]+$/, "", id)
         skip = (id in l)
     }
     !skip' $ASSEMBLYDIR/contaminant_genes.txt "$ASSEMBLY" > "${ASSEMBLY%.fasta}.rmConta.fasta"

# Check counts before and after
echo "N transcripts before clean:"
grep -c "^>" "$ASSEMBLY"

echo "N transcripts after clean:"
grep -c "^>" "${ASSEMBLY%.fasta}.rmConta.fasta"

echo "Sanity check: none of the contaminants remain"
head -3 contaminant_genes.txt | while read id; do
    grep -c "$id" "${ASSEMBLY%.fasta}.rmConta.fasta" && echo "WARNING: $id still present!" || echo "OK: $id removed"
done

echo "Check assembly stats:"
module purge
module add Trinity

$TRINITY_HOME/util/TrinityStats.pl "${ASSEMBLY%.fasta}.rmConta.fasta"

echo "Check BUSCO:"
module add BUSCO

busco -i "${ASSEMBLY%.fasta}.rmConta.fasta" -l fungi_odb10 -m transcriptome -c 10 -o "BUSCO_MergedFungiTranscriptome_eukaryoteHits.rmRSEMconta"
## l = current fungi database (will be downloaded automatically)

