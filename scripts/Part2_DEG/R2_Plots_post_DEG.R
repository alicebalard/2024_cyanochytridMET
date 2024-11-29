#################
## Libraries load
## BiocManager::install('EnhancedVolcano') # package for pretty volcano plots
library(EnhancedVolcano)
library(ggplot2)
library(stringr)

## Files used for the DESeq2 analysis:
# assembly: /scratch/alicebalard/outData/assemblyMergedFungi/trinity_out_dir/Trinity_eukaryoteHits.fasta
# gene_trans_map: /scratch/alicebalard/outData/assemblyMergedFungi/trinity_out_dir/Trinity_eukaryoteHits.fasta.gene_trans_map
# annotation: /scratch/alicebalard/outData/assemblyMergedFungi/annotation/assemblyMergedFungi_filterEuk_simplified.tsv

## Read in DESeq2 results:
chytrid_inf_effect_control <- readRDS("../../data/run_DESEQ2_Erika/chytrid_inf_effect_control.rds")
chytrid_inf_effect_met <- readRDS("../../data/run_DESEQ2_Erika/chytrid_inf_effect_met.rds")
chytrid_met_effect_1org <- readRDS("../../data/run_DESEQ2_Erika/chytrid_met_effect_1org.rds")
chytrid_met_effect_2orgs <- readRDS("../../data/run_DESEQ2_Erika/chytrid_met_effect_2orgs.rds")
cyano_inf_effect_control <- readRDS("../../data/run_DESEQ2_Erika/cyano_inf_effect_control.rds")
cyano_inf_effect_met <- readRDS("../../data/run_DESEQ2_Erika/cyano_resr_inf_effect_met.rds")
cyano_met_effect_1org <- readRDS("../../data/run_DESEQ2_Erika/cyano_met_effect_1org.rds")
cyano_met_effect_2org <- readRDS("../../data/run_DESEQ2_Erika/cyano_met_effect_2orgs.rds")


## TBC

### gene number and name without any additional information (e.g. gene 1 TRINITY_DN43456_c0_g1_i1).
new_gene_trans_map_fungi <- read.table("../../data/new_gene_trans_map_fungi.txt")

## 1. Chytrid gene expression

# addInfoChyt <- function(res_chyt){

res_chyt = chytrid_inf_effect_control

  ## Match correct names
  res_chyt$transcript_id <- new_gene_trans_map_fungi$V2[
    match(rownames(res_chyt), new_gene_trans_map_fungi$V1)]
  
  ## Match gene name based on blastx top hit
  res_chyt$gene_name_blastx <- annotChyt$gene_name[
    match(res_chyt$transcript_id, annotChyt$transcript_id)]
  return(res_chyt)
}

resr_inf_effect_control <- addInfoChyt(res_chyt = chytrid_inf_effect_control)

makeVolcano <- function(res, title, subtitle, chyt=T){#by default for chytrid
  if (chyt == T){
    res$geneName = ifelse(res$gene_name_blastx != ".", res$gene_name_blastx, "")
  }
  
  ### Extracting significant differentially expressed genes
  
  ## Convert Ensembl names to gene names (if they are known)
  results_df = as.data.frame(res)
  results_df = results_df[order(results_df$padj),]
  
  # Subset the results to keep only significant genes
  ressig = results_df[results_df$padj < 0.05 & !is.na(results_df$padj),]
  
  ## Volcano plot
  plot = EnhancedVolcano(results_df,
                         lab = results_df$geneName,
                         x = 'log2FoldChange', title = title, subtitle = subtitle,
                         y = 'padj', pCutoff = 0.05,
                         drawConnectors = TRUE, labSize = 3)
  
  return(list(signifGenes = ressig, plot = plot))
}

V_inf_effect_control <- makeVolcano(
  res = resr_inf_effect_control,
  title = "infected vs non infected P. agardhii gene expression",
  subtitle = "MET unexposed")

## open bigger window
dev.new(width = 15, height = 12)
V_inf_effect_control$plot

## Extreme log2 fold changes could also arise from issues related to:
### Low counts: If the gene has very low counts in one condition (e.g., zero or near-zero),
# even a small count in the other condition can lead to a high fold change.
### Normalization methods: Ensure appropriate normalization techniques are applied during 
# the analysis. Poor normalization can distort fold change calculations.

## Idea: explore raw data for outliers

## If we remove OL
hist(resr_inf_effect_control$log2FoldChange, breaks = 100)

V_inf_effect_control_noOL <- makeVolcano(
  res = resr_inf_effect_control[!is.na(resr_inf_effect_control$log2FoldChange) &
                                  abs(resr_inf_effect_control$log2FoldChange) < 10,],
  title = "infected vs non infected P. agardhii gene expression",
  subtitle = "MET unexposed")

dev.new(width = 15, height = 12)
V_inf_effect_control_noOL$plot

