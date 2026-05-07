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
THRESHOLD=1  # adjust if needed

# ============================================================
# STEP 1: Get contaminant gene IDs
# TRINITY genes with any count in cyano-alone samples (cols 14-19, 32-37)
# NOTE: header has leading tab causing NF=36 vs data NF=37, so hardcode columns
# control_cyano_C1-C6 = data cols 14-19
# met_cyano_C10-C9   = data cols 32-37
# ============================================================
awk '
NR==1 {
    # Store column index by name - no position ambiguity
    for(i=1;i<=NF;i++) {
        if($i ~ /cyano/) cyanocol[i]=1
    }
    # Print which columns were found for verification
    for(i in cyanocol) print "CYANO col " i ": " $i > "/dev/stderr"
    next
}
/^TRINITY/ {
    cyano_sum=0
    for(i in cyanocol) cyano_sum += $(i+1)  # +1 offset for leading tab
    if(cyano_sum > 0) print $1
}' "$MATRIX" > "$ASSEMBLYDIR/contaminant_genes.txt"

echo "Contaminant genes: $(wc -l < $ASSEMBLYDIR/contaminant_genes.txt)"

# ============================================================
# STEP 2: Expand gene IDs to isoform IDs using FASTA headers
# Gene IDs:   TRINITY_DN888_c0_g1
# Isoform IDs: TRINITY_DN888_c0_g1_i18, _i2, etc.
# ============================================================
grep "^>" "$ASSEMBLY" | sed 's/>//' | awk '
NR==FNR { l[$1]=1; next }
{
    gene=$1
    gsub(/_i[0-9]+$/, "", gene)
    if(gene in l) print $1
}' "$ASSEMBLYDIR/contaminant_genes.txt" - > "$ASSEMBLYDIR/contaminant_isoforms.txt"

echo "Contaminant isoforms: $(wc -l < $ASSEMBLYDIR/contaminant_isoforms.txt)"
echo "Total isoforms:       $(grep -c '^>' $ASSEMBLY)"

# ============================================================
# STEP 3: Filter contaminating isoforms from FASTA
# ============================================================
awk 'NR==FNR {drop[$1]=1; next}
     /^>/ { id=substr($1,2); skip=(id in drop) }
     !skip' "$ASSEMBLYDIR/contaminant_isoforms.txt" "$ASSEMBLY" \
     > "${ASSEMBLY%.fasta}.rmConta.fasta"

echo "Isoforms after cleaning: $(grep -c '^>' ${ASSEMBLY%.fasta}.rmConta.fasta)"
echo "Expected:                $(( $(grep -c '^>' $ASSEMBLY) - $(wc -l < $ASSEMBLYDIR/contaminant_isoforms.txt) ))"

# ============================================================
# STEP 3b: Filter contaminant genes from count matrix
# ============================================================
awk 'NR==FNR {drop[$1]=1; next}
     NR==1 {print; next}          # keep header
     !($1 in drop)                # keep non-contaminant rows
' "$ASSEMBLYDIR/contaminant_genes.txt" "$MATRIX" \
> "${MATRIX%.matrix}.rmConta.matrix"

echo "Genes in original matrix: $(wc -l < $MATRIX)"
echo "Genes in clean matrix:    $(wc -l < ${MATRIX%.matrix}.rmConta.matrix)"
echo "Difference (contaminants removed): $(( $(wc -l < $MATRIX) - $(wc -l < ${MATRIX%.matrix}.rmConta.matrix) ))"

# ============================================================
# STEP 4: Sanity check - verify no contaminants remain
# ============================================================
echo "Sanity check (should all be 0):"
shuf -n 5 "$ASSEMBLYDIR/contaminant_genes.txt" | while read gene; do
    count=$(grep -c "^>${gene}_i" "${ASSEMBLY%.fasta}.rmConta.fasta" 2>/dev/null || echo 0)
    echo "  $gene: $count sequences remaining"
done

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

