#!/bin/bash
# S01_prepareCombinedTranscriptome.sh	

# Prepares combined chytrid + cyanobacteria transcriptome for RSEM mapping
# FIX: removed _g1_i1 suffix (was misleading; gene_trans_map assigns genes by NR anyway)
# FIX: fixed original_headers_cds.txt to keep full original headers for traceability
# FIX: added set -euo pipefail for safer execution
# FIX: updated output paths to match new working directory

set -euo pipefail

T_CHY=/scratch/alicebalard/outData/assemblies/assemblyMergedFungi/Trinity_eukaryoteHits.fasta
T_CHY_GTM=/scratch/alicebalard/outData/assemblies/assemblyMergedFungi/Trinity_eukaryoteHits.fasta.gene_trans_map
T_CYA=/scratch/alicebalard/outData/mergedTransc/GCF_904830935.1_P._agardhii_No.976_rna_from_genomic.fna

OUTDIR=/scratch/alicebalard/outData/RSEM
mkdir -p $OUTDIR
cd $OUTDIR

echo "[S01] Starting transcriptome preparation..."

## -------------------------------------------------------
## 1. Prepare the gene_trans_map for Planktothrix transcriptome
## -------------------------------------------------------

## Rename cyanobacterium CDS headers:
##   - replace ">lcl|" with ">cyano_"
##   - strip everything after the first space (description)
## FIX: removed the "_g1_i1" suffix — it implied Trinity isoform structure
##      but gene assignment is done by NR in awk below, making the suffix misleading
sed '/^>/ s/^>lcl|/>cyano_/; /^>/ s/ .*//' $T_CYA > cyano_transcriptome_renamed_cds.fna

echo "[S01] Renamed cyano headers."

## Create lookup table to trace back renamed IDs to original FASTA headers
## FIX: original_headers_cds.txt now keeps the full original header (no partial sed transform)
##      so the two columns in the lookup table are directly comparable
grep '^>' $T_CYA | sed 's/^>//' > original_headers_cds.txt
grep '^>' cyano_transcriptome_renamed_cds.fna | sed 's/>//' > renamed_headers_cds.txt
paste renamed_headers_cds.txt original_headers_cds.txt > header_lookup_table.txt

echo "[S01] Created header lookup table."

## Extract transcript IDs (no '>' symbol) for gene_trans_map generation
grep "^>" cyano_transcriptome_renamed_cds.fna | sed 's/>//' > transcriptome_ids_cds.txt

## Generate gene_trans_map: one unique gene per CDS entry (appropriate for a reference genome)
## Format: gene_id <TAB> transcript_id  (required by Trinity/RSEM utilities)
## NB: each CDS gets its own cyano_geneN since this is a reference genome, not a de novo assembly
awk '{print "cyano_gene" NR "\t" $1}' transcriptome_ids_cds.txt > gene_trans_map_cds.txt

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
echo "[S01] Done. Outputs:"
echo "  $OUTDIR/assembly_both.fna"
echo "  $OUTDIR/combined_gene_trans_map.txt"
echo "  $OUTDIR/header_lookup_table.txt"
