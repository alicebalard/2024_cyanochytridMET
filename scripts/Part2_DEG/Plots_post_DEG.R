#################
## Libraries load
## BiocManager::install('EnhancedVolcano') # package for pretty volcano plots
library(EnhancedVolcano)
library(ggplot2)
library(stringr)

############
## Data load

## 1) new_gene_trans_map_fungi.txt - this file contains the gene number and its name without any additional information (e.g. gene 1 TRINITY_DN43456_c0_g1_i1).
new_gene_trans_map_fungi <- read.table("../../data/new_gene_trans_map_fungi.txt")

## 2) transcripts_from_filtered_fungi.txt - this file contains the full names of the genes with additional information (TRINITY_DN43456_c0_g1_i1 len=541 path=[0:0-540]).
## not needed now

## 3) Chytrid annotation
annotChyt <- read.csv("../../data/allFungiTrinot_simplified.tsv", sep ="\t")

## extract GO terms in their own column for GO analysis downstream
annotChyt$GO <- 
  sapply(annotChyt$gene_ontology_BLASTX, function(x) {
    go_terms <- str_extract_all(x, "GO:\\d+")
    paste(unlist(go_terms), collapse = ",")
  })

## Read in DESeq2 results:
## 1. Chytrid gene expression

## MET control, chytrid infecting cyanobacteria vs not infecting (zoospores)
resr_inf_effect_control <- readRDS("../../data/resr_inf_effect_control.rds")

## Match correct names
resr_inf_effect_control$transcript_id <- new_gene_trans_map_fungi$V2[
  match(rownames(resr_inf_effect_control), new_gene_trans_map_fungi$V1)]

## Match gene name based on blastx top hit
resr_inf_effect_control$gene_name_blastx <- annotChyt$gene_name[
  match(resr_inf_effect_control$transcript_id, annotChyt$transcript_id)]


head(resr_inf_effect_control)

resr_inf_effect_control[100,]
baseMean log2FoldChange     lfcSE      stat    pvalue      padj
<numeric>      <numeric> <numeric> <numeric> <numeric> <numeric>
  gene10087   3.85973        3.54772   3.70592  0.957312   0.33841  0.953593
transcript_id gene_name_blastx
<character>      <character>
  gene10087 TRINITY_DN836_c1_g2_i3       METE_SCHPO

annotChyt[annotChyt$transcript_id %in% "TRINITY_DN836_c1_g2_i3",]

makeVolcano <- function(res, title, subtitle){
  ### Extracting significant differentially expressed genes
  
  ## Convert Ensembl names to gene names (if they are known)
  results_df = as.data.frame(res)
  # results_df$Ensembl_id = row.names(results_df)
  results_df = results_df[order(results_df$padj),]
  
  # results_genes = gconvert(row.names(res), organism = "mmusculus",
  # target = "ENTREZGENE_ACC", filter_na = FALSE)
  # add the gene names
  # results_df = merge(results_df,
  #                    results_genes[,c("input", "target", "name", "description")],
  #                    by.x = "Ensembl_id", by.y = "input")
  
  # results_df$Name <- ifelse(is.na(results_df$name), results_df$Ensembl_id, results_df$name)
  
  # Subset the results to keep only significant genes
  ressig = results_df[results_df$padj < 0.05 & !is.na(results_df$padj),]
  
  ## Volcano plot
  plot = EnhancedVolcano(results_df,
                         lab = rownames(results_df),
                         x = 'log2FoldChange', title = title, subtitle = subtitle,
                         y = 'padj', pCutoff = 0.05,
                         drawConnectors = TRUE)
  
  return(list(signifGenes = ressig, plot = plot))
}

V_inf_effect_control <- makeVolcano(
  res = resr_inf_effect_control,
  title = "infected vs non infected P. agardhii",
  subtitle = "MET unexposed")

## open bigger window
dev.new(width = 15, height = 12)
V_inf_effect_control$plot
