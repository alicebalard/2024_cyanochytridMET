T_CHY=/scratch/alicebalard/outData/assemblyMergedFungi/trinity_out_dir/Trinity_eukaryoteHits.fasta
T_CHY_GTM=/scratch/alicebalard/outData/assemblyMergedFungi/trinity_out_dir/Trinity_eukaryoteHits.fasta.gene_trans_map

T_CYA=/scratch/alicebalard/outData/mergedTransc/GCF_904830935.1_P._agardhii_No.976_rna_from_genomic.fna

## To be processed by Erika, we create data in her folder
cd /scratch/erikamr/cyano_chytrid_met/data/

# 1.Prepare the gene_trans_map for Planktothrix transcriptome

## adjust the names of the cds cyanobacterium
sed '/^>/ s/^>lcl|/>cyano_/; /^>/ s/ .*//; /^>/ s/$/_g1_i1/' $T_CYA > cyano_transcriptome_renamed_cds.fna

## Create a file to store the original Fasta headers and the new name to be able to trace back to genes and proteins later
grep '^>' $T_CYA | sed 's/^lcl|/>cyano_/' > original_headers_cds.txt
grep '^>' cyano_transcriptome_renamed_cds.fna > renamed_headers_cds.txt
paste renamed_headers.txt original_headers.txt > header_lookup_table.txt

## extract headers. I will have a list of all transcript IDs without the > symbol
grep "^>" cyano_transcriptome_renamed_cds.fna | sed 's/>//' > transcriptome_ids_cds.txt

## Generate the gene_trans_map file by appending a unique gene identifier to each transcript
awk '{print "cyano_gene" NR "\t" $1}' transcript_ids_cds.txt > gene_trans_map_cds.txt

# 2. Combined gene trans map
cat $T_CHY_GTM gene_trans_map_cds.txt > combined_gene_trans_map_cds_final_hope.txt

# 3. Create a new combined file of transcriptomes
cat $T_CHY cyano_transcriptome_renamed_cds.fna > assembly_both_cds_final_hope.fna

## Outputs:
## /scratch/erikamr/cyano_chytrid_met/data/combined_gene_trans_map_cds_final_hope.txt
## /scratch/erikamr/cyano_chytrid_met/data/assembly_both_cds_final_hope.fna
