## 11th of December 2024
## A. Balard

## load all files needed for the project

## design table
samples_data <- read.table(here("data/samples_file.txt"), header = F, sep = "\t")

############################################################
## Chytrid
# assembly: /scratch/alicebalard/outData/assemblies/assemblyMergedFungi/Trinity_eukaryoteHits.rmConta.fasta
# gene_trans_map: /scratch/alicebalard/outData/assemblies/assemblyMergedFungi/Trinity_eukaryoteHits.fasta.gene_trans_map

# annotation: /scratch/alicebalard/outData/assemblies/assemblyMergedFungi/annotation/assemblyMergedFungi_filterEuk_simplified_GOKegg.tsv
annotationChytridFULL <- read.csv(here("gitignore/assemblyMergedFungi_filterEuk_simplified_GOKegg.tsv"), sep = "\t")
## extract GO terms
annotationChytridFULL$GO.accession <- str_extract_all(annotationChytridFULL$gene_ontology_BLASTX, "GO:\\d+")

## make an annotation df with: customGeneName, gene_name, GOaccession, Kegg ##
annotationChytrid <- annotationChytridFULL[c("X.gene_id", "gene_name", "GO.accession", "Kegg")] %>% 
  unnest(GO.accession, keep_empty = T) %>%
  dplyr::rename("custom_gene_name" = "X.gene_id") %>% data.frame()

annotationChytrid <- unique(annotationChytrid)

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

message("=== Chytrid annotation summary ===")
message("Total genes:               ", n_distinct(annotationChytridFULL$X.gene_id))
message("With GO terms:             ", annotationChytrid %>% filter(!is.na(GO.accession)) %>% distinct(custom_gene_name) %>% nrow())
message("Without GO terms:          ", annotationChytrid %>% filter(is.na(GO.accession))  %>% distinct(custom_gene_name) %>% nrow())
message("Unique gene-GO pairs:      ", nrow(annotationChytrid))
message("Unique GO terms:           ", nrow(GO_chytrid))

## to keep:
annotationChytrid_final <- annotationChytridFULL
names(annotationChytrid_final)[names(annotationChytrid_final) %in% "X.gene_id"] <- "gene_id"

rm(split_by_caret, split_by_backtick, result_list, result_df, annotationChytridFULL, annotationChytrid)

############################################################
## Cyanobacteria annotation

# 1. Load header lookup table (full info & long name for each gene)
header_lookup_table <- read.csv(here("data/header_lookup_table.txt"), sep="\t", header=FALSE)
header_lookup_table$V1 <- gsub(">", "", header_lookup_table$V1)
names(header_lookup_table) <- c("gene_name_long", "info")

# 2. Load gene/transcript map (long and short name key table)
gene_trans_map_cyano <- read.csv(here("data/combined_gene_trans_map.txt"), sep=" ", header=FALSE)
gene_trans_map_cyano <- gene_trans_map_cyano[!grepl("TRINITY", gene_trans_map_cyano$V1),]
gene_trans_map_cyano <- tidyr::separate(gene_trans_map_cyano, V1, into=c("V1","V2"), sep="\\t")
names(gene_trans_map_cyano) <- c("custom_gene_name", "gene_name_long")

# 3. Load GTF once and build both bridge and GO table
gtf_full <- import(here("data/GCF_904830935.1_P._agardhii_No.976_genomic.gtf")) %>%
  data.frame()

# Bridge: locus_tag → protein_id + product name
locus_to_protein <- gtf_full %>%
  filter(!is.na(protein_id), !is.na(locus_tag)) %>%
  distinct(locus_tag, protein_id, product)  # no gene here, comes from annotationCyanoFULL

# GO terms in long format
annotationCyano_GO <- gtf_full %>%
  filter(!is.na(Ontology_term)) %>%
  dplyr::select(protein_id, product, gene, Ontology_term,
                go_function, go_process, go_component) %>%
  distinct() %>%
  pivot_longer(
    cols = c(go_function, go_process, go_component),
    names_to = "variable",
    values_to = "value"
  ) %>%
  filter(!is.na(value))

# 4. Merge and extract fields from info column
annotationCyanoFULL <- merge(gene_trans_map_cyano, header_lookup_table) %>%
  mutate(
    full_id   = str_extract(info, "(?<=>lcl\\|)\\S+"),
    gene      = str_extract(info, "(?<=\\[gene=)[^\\]]+"),
    locus_tag = str_extract(info, "(?<=\\[locus_tag=)[^\\]]+"),
    gene_id   = str_extract(info, "(?<=\\[db_xref=GeneID:)[^\\]]+"),
    protein   = str_extract(info, "(?<=\\[protein=)[^\\]]+")
  ) %>%
  left_join(
    locus_to_protein %>%
      dplyr::select(locus_tag, protein_id, product) %>%
      distinct(),
    by = "locus_tag"
  ) %>%
  left_join(
    annotationCyano_GO %>% distinct(protein_id, Ontology_term, variable, value),
    by = "protein_id",
    relationship = "many-to-many"
  )

# 5. Slim annotation for DESeq2 (one row per gene + GO accession)
annotationCyano <- annotationCyanoFULL %>%
  distinct(custom_gene_name, gene, Ontology_term) %>%
  dplyr::rename(GO.accession = Ontology_term)

# 6. GO reference table (accession + ontology category + term name)
GO_cyano <- annotationCyanoFULL %>%
  filter(!is.na(Ontology_term), !is.na(variable), !is.na(value)) %>%
  transmute(
    GO.accession = Ontology_term,
    GO.ontology = case_when(
      variable == "go_function"  ~ "molecular_function",
      variable == "go_process"   ~ "biological_process",
      variable == "go_component" ~ "cellular_component",
      TRUE ~ NA_character_
    ),
    GO.name = sapply(str_split(value, "\\|"), function(x) x[1])
  ) %>%
  distinct() %>%
  na.omit()

message("=== Cyano annotation summary ===")
message("Total genes:               ", n_distinct(annotationCyanoFULL$custom_gene_name))
message("With protein_id:           ", annotationCyanoFULL %>% filter(!is.na(protein_id)) %>% distinct(custom_gene_name) %>% nrow())
message("With GO terms:             ", annotationCyanoFULL %>% filter(!is.na(Ontology_term)) %>% distinct(custom_gene_name) %>% nrow())
message("Unmatched (no protein_id): ", annotationCyanoFULL %>% filter(is.na(protein_id)) %>% distinct(custom_gene_name) %>% nrow())
message("Unique GO terms:           ", nrow(GO_cyano))
message("Unique gene-GO pairs: ", nrow(annotationCyano))
message("Unique GO terms:      ", nrow(GO_cyano))

## to keep:
annotationCyano_final <- annotationCyanoFULL

rm(gtf_full, header_lookup_table, locus_to_protein, annotationCyano_GO,
   annotationCyano,annotationCyanoFULL, gene_trans_map_cyano)
