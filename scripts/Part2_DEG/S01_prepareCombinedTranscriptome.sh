#!/bin/bash
# S01_prepareCombinedTranscriptome.sh	

# Prepares combined chytrid + cyanobacteria transcriptome for RSEM mapping
# FIX: removed _g1_i1 suffix (was misleading; gene_trans_map assigns genes by NR anyway)
# FIX: fixed original_headers_cds.txt to keep full original headers for traceability
# FIX: added set -euo pipefail for safer execution
# FIX: updated output paths to match new working directory
#SBATCH --job-name=S01
#SBATCH --mail-user=alicebalard@zedat.fu-berlin.de
#SBATCH --mail-type=end,fail
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=2GB
#SBATCH --time=24:00:00
#SBATCH --output=/home/alicebalard/Scripts/AliceScripts/cyanochytridMET/scripts/logs_dir/%x.%j.out
#SBATCH --error=/home/alicebalard/Scripts/AliceScripts/cyanochytridMET/scripts/logs_dir/%x.%j.err
#SBATCH --qos=standard
#SBATCH --partition=main,begendiv
#SBATCH --constraint=no_gpu

set -euo pipefail                                                                                                                                                                       
export LMOD_DISABLE_SAME_NAME_AUTOSWAP=no

T_CHY=/scratch/alicebalard/outData/assemblies/assemblyMergedFungi/Trinity_eukaryoteHits.fasta
T_CHY_GTM=/scratch/alicebalard/outData/assemblies/assemblyMergedFungi/Trinity_eukaryoteHits.fasta.gene_trans_map
T_CYA=/scratch/alicebalard/outData/mergedTransc/GCF_904830935.1_P._agardhii_No.976_cds_from_genomic.fna

OUTDIR=/scratch/alicebalard/outData/RSEM
mkdir -p $OUTDIR
cd $OUTDIR

echo "[S01] Starting transcriptome preparation..."

## -------------------------------------------------------
## 1. Prepare the gene_trans_map for Planktothrix transcriptome
## -------------------------------------------------------

## Rename cyanobacterium CDS headers:
sed '/^>/ s/^>lcl|/>cyano_/; /^>/ s/ .*//; /^>/ s/$/_g1_i1/' $T_CYA > cyano_transcriptome_renamed_cds.fna

echo "[S01] Renamed cyano headers."

## Create a file to store the original Fasta headers and the new name to be able to trace back to genes and proteins later
grep '^>' $T_CYA | sed 's/^lcl|/>cyano_/' > original_headers.txt
grep '^>' cyano_transcriptome_renamed_cds.fna > renamed_headers.txt
paste renamed_headers.txt original_headers.txt > header_lookup_table.txt

echo "[S01] Created header lookup table."

## Extract transcript IDs (no '>' symbol) for gene_trans_map generation
grep "^>" cyano_transcriptome_renamed_cds.fna | sed 's/>//' > transcript_ids_cds.txt

## Generate the gene_trans_map file by appending a unique gene identifier to each transcript
awk '{print "cyano_gene" NR "\t" $1}' transcript_ids_cds.txt > gene_trans_map_cds.txt

echo "[S01] Generated cyano gene_trans_map ($(wc -l < gene_trans_map_cds.txt) entries)."

## -------------------------------------------------------
## 2. Combine gene_trans_maps
## -------------------------------------------------------
cat $T_CHY_GTM gene_trans_map_cds.txt > combined_gene_trans_map.txt

## Sanity check: no duplicate transcript IDs across the two transcriptomes
N_TOTAL=$(wc -l < combined_gene_trans_map.txt)
N_UNIQUE=$(awk '{print $2}' combined_gene_trans_map.txt | sort -u | wc -l)
if [ "$N_TOTAL" -ne "$N_UNIQUE" ]; then
    echo "[S01] WARNING: duplicate transcript IDs detected in combined gene_trans_map!" >&2
    echo "[S01] Total: $N_TOTAL, Unique: $N_UNIQUE" >&2
else
    echo "[S01] gene_trans_map OK: $N_TOTAL unique transcript entries."
fi

## -------------------------------------------------------
## 3. Combine FASTA transcriptomes
## -------------------------------------------------------
cat $T_CHY cyano_transcriptome_renamed_cds.fna > assembly_both.fna

echo "[S01] Combined FASTA created ($(grep -c '^>' assembly_both.fna) sequences)."

## -------------------------------------------------------
## Outputs
## -------------------------------------------------------
cp $OUTDIR/combined_gene_trans_map.txt /home/alicebalard/Scripts/AliceScripts/cyanochytridMET/data/.
cp $OUTDIR/header_lookup_table.txt /home/alicebalard/Scripts/AliceScripts/cyanochytridMET/data/.

echo "[S01] Done. Outputs:"
echo "  $OUTDIR/assembly_both.fna"
echo "  /home/alicebalard/Scripts/AliceScripts/cyanochytridMET/data/combined_gene_trans_map.txt"
echo "  /home/alicebalard/Scripts/AliceScripts/cyanochytridMET/data/header_lookup_table.txt"
