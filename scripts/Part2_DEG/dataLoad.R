## 11th of December 2024
## A. Balard

## load all files needed for the project

## design table
samples_data <- read.table(
  "../../data/sample_data_remove_r.txt", header = TRUE, sep = "\t")

## Chytrid

# assembly: /scratch/alicebalard/outData/assemblyMergedFungi/trinity_out_dir/Trinity_eukaryoteHits.fasta
# gene_trans_map: /scratch/alicebalard/outData/assemblyMergedFungi/trinity_out_dir/Trinity_eukaryoteHits.fasta.gene_trans_map

# annotation: /scratch/alicebalard/outData/assemblyMergedFungi/annotation/assemblyMergedFungi_filterEuk_simplified.tsv
annotationChytridFULL <- read.csv("../../gitignore/assemblyMergedFungi_filterEuk_simplified.tsv", sep = "\t")
## extract GO terms
annotationChytridFULL$GO.accession <- str_extract_all(annotationChytridFULL$gene_ontology_BLASTX, "GO:\\d+")

########################################################################
## make an annotation df with: customGeneName, gene_name, GOaccession ##
annotationChytrid <- annotationChytridFULL[c("X.gene_id", "gene_name", "GO.accession")] %>% 
  unnest(GO.accession, keep_empty = T) %>%
  dplyr::rename("custom_gene_name" = "X.gene_id") %>% data.frame()

annotationChytrid <- unique(annotationChytrid)

############################################################
## make a GOdf table with GO accession, name and ontology ##
split_by_backtick <- lapply(annotationChytridFULL$gene_ontology_BLASTX,
                            function(x) strsplit(x, "`")[[1]])
split_by_caret <- lapply(split_by_backtick, function(go_list) {
  lapply(go_list, function(x) strsplit(x, "\\^")[[1]])
})
result_list <- unlist(split_by_caret, recursive = FALSE)

# Create a data frame while ensuring each entry has 3 elements
result_df <- do.call(rbind, lapply(result_list, function(x) {
  length(x) <- 3  # Ensure each entry has 3 elements (fill with NA if necessary)
  return(x)
}))
# Convert to data frame and set column names
GO_chytrid <- as.data.frame(result_df, stringsAsFactors = FALSE)
GO_chytrid=na.omit(GO_chytrid)
GO_chytrid=unique(GO_chytrid)

names(GO_chytrid) <- c("GO.accession", "GO.ontology", "GO.name")

## Cyanobacteria

########################################################################
## make an annotation df with: customGeneName, gene_name, GOaccession ##
# 1. header_lookup_table.txt This file contains the shorter names I tried for the cyanobacterium genes and the full name separated by a tab. For the short name I add a number that represents the gene (1...4428) and the g1 and i1 even though each gene has only one isoform. For example
# cyano_NZ_LR882952.1_cds_WP_027249510.1_1_g1_i1 >lcl|NZ_LR882952.1_cds_WP_027249510.1_1 [locus_tag=NMG88_RS00005] [db_xref=GeneID:77286138] [protein=TIGR01548 family HAD-type hydrolase] [protein_id=WP_027249510.1] [location=135..920] [gbkey=CDS]
header_lookup_table = read.csv("../../data/run_DESEQ2_Erika/header_lookup_table.txt", sep="\t", header = F)
header_lookup_table$V1 = gsub(">", "", header_lookup_table$V1)  
names(header_lookup_table) <- c("gene_name_long", "info")

# 2. gene_trans_map_cds.txt This file contains the short name and the shortest name separated by a tab. For example
# cyano_gene1 cyano_NZ_LR882952.1_cds_WP_027249510.1_1_g1_i1
gene_trans_map_cyano = read.csv("../../data/run_DESEQ2_Erika/gene_trans_map_cds.txt", sep="\t", header = F)
names(gene_trans_map_cyano) <- c("custom_gene_name", "gene_name_long")

annotationCyanoFULL <- merge(gene_trans_map_cyano, header_lookup_table)

# extract protein name
annotationCyanoFULL$protein = sub(".*\\[protein=(.*?)\\].*", "\\1", annotationCyanoFULL$info)
annotationCyanoFULL$protein_id = sub(".*\\[protein.id=(.*?)\\].*", "\\1", annotationCyanoFULL$info)
## not all have gene name
annotationCyanoFULL$gene_name = sub(".*\\[db_xref=(.*?)\\].*", "\\1", annotationCyanoFULL$info)
annotationCyanoFULL$gene_name[grep("gene=", annotationCyanoFULL$info)] = 
  sub(".*\\[gene=(.*?)\\].*", "\\1", annotationCyanoFULL$info[grep("gene=", annotationCyanoFULL$info)])

## add GO
annotationCyano2 <- import("../../data/annotations/GCF_904830935.1_P._agardhii_No.976_genomic.gtf") 
annotationCyano2 <- annotationCyano2[!is.na(annotationCyano2$Ontology_term),
                                   c("protein_id", "product", "gene", "Ontology_term", "go_function", "go_process", "go_component")] %>%
  data.frame()
annotationCyano2 <- annotationCyano2 %>% melt(id.vars = names(annotationCyano2)[1:9])
annotationCyanoFULL <- merge(annotationCyanoFULL, annotationCyano2, all = T)
                          
annotationCyano = annotationCyanoFULL[c("custom_gene_name", "gene_name", "Ontology_term")] %>%
                      dplyr::rename("GO.accession" = "Ontology_term")
                     
############################################################
## make a GOdf table with GO accession, name and ontology ##
GO_cyano = data.frame(GO.ID = annotationCyanoFULL$Ontology_term,
                    GO.cat = ifelse(annotationCyanoFULL$variable == "go_function", "molecular_biological_function", ifelse(
                      annotationCyanoFULL$variable == "go_process", "biological_process", ifelse(
                        annotationCyanoFULL$variable == "go_component", "cellular_component", NA))),
                    GO.term = sapply(str_split(annotationCyanoFULL$value, "\\|"), function(x) x[1]))

GO_cyano=na.omit(GO_cyano)
GO_cyano=unique(GO_cyano)
names(GO_cyano) <- c("GO.accession", "GO.ontology", "GO.name")

source("functions.R")
