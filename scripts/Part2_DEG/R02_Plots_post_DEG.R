source("R01_Erika_DESeq2.R")

##################################################
## Rename DESEq2 files with meaningfull gene names
## TO DO annotate longuest isoform (each isoform has an annotation otherwise, so we'll take the first for now)
giveGoodGeneNameChytrid <- function(chytridDESEq2){
  row.names(chytridDESEq2) = annotationChytrid$gene_name[
    match(row.names(chytridDESEq2), annotationChytrid$X.gene_id)]
  return(chytridDESEq2)
}
chytrid_inf_effect_control <- giveGoodGeneNameChytrid(contrast_chytridgenome$resr_inf_effect_control)
chytrid_inf_effect_met <- giveGoodGeneNameChytrid(contrast_chytridgenome$resr_inf_effect_met)
chytrid_met_effect_1org <- giveGoodGeneNameChytrid(contrast_chytridgenome$resr_met_effect_1org)
chytrid_met_effect_2orgs <- giveGoodGeneNameChytrid(contrast_chytridgenome$resr_met_effect_2orgs)

giveGoodGeneNameCyano <- function(cyanoDESEq2){
  row.names(cyanoDESEq2) = gene_trans_map_cyano$protein[
    match(row.names(cyanoDESEq2), gene_trans_map_cyano$V1)]
  return(cyanoDESEq2)
}

cyano_inf_effect_control <- giveGoodGeneNameCyano(contrast_cyanogenome$resr_inf_effect_control)
cyano_inf_effect_met <- giveGoodGeneNameCyano(contrast_cyanogenome$resr_inf_effect_met)
cyano_met_effect_1org <- giveGoodGeneNameCyano(contrast_cyanogenome$resr_met_effect_1org)
cyano_met_effect_2orgs <- giveGoodGeneNameCyano(contrast_cyanogenome$resr_met_effect_2orgs)

mylistResDESEQ2 <- list(chytrid_inf_effect_control=chytrid_inf_effect_control,
                        chytrid_inf_effect_met=chytrid_inf_effect_met,
                        chytrid_met_effect_1org=chytrid_met_effect_1org,
                        chytrid_met_effect_2orgs=chytrid_met_effect_2orgs,
                        cyano_inf_effect_control=cyano_inf_effect_control,
                        cyano_inf_effect_met=cyano_inf_effect_met,
                        cyano_met_effect_1org=cyano_met_effect_1org,
                        cyano_met_effect_2orgs=cyano_met_effect_2orgs)

###############
## Run plots ##
###############

makeVolcano <- function(res, title){
  results_df = as.data.frame(res)
  results_df = results_df[order(results_df$padj),]
  
  # Subset the results to keep only significant genes
  ressig = results_df[results_df$padj < 0.05 & !is.na(results_df$padj),]
  
  ## Volcano plot
  plot = EnhancedVolcano(results_df,
                         lab = row.names(results_df),
                         x = 'log2FoldChange', title = title,
                         y = 'padj', pCutoff = 0.05,
                         drawConnectors = TRUE, labSize = 3)
  
  return(list(signifGenes = ressig, plot = plot))
}

## open bigger window
V_chytrid_inf_effect_control <- makeVolcano(
  res = chytrid_inf_effect_control,
  title = "infecting vs non infecting chytrid gene expression, no MET")

V_chytrid_inf_effect_met <- makeVolcano(
  res = chytrid_inf_effect_met,
  title = "infecting vs non infecting chytrid gene expression, MET")

V_chytrid_met_effect_1org <- makeVolcano(
  res = chytrid_met_effect_1org,
  title = "MET effect on chytrid gene expression")

V_chytrid_met_effect_2orgs <- makeVolcano(
  res = chytrid_met_effect_2orgs,
  title = "MET effect on chytrid infecting bacteria gene expression")

V_cyano_inf_effect_control <- makeVolcano(
  res = cyano_inf_effect_control,
  title = "infected vs non infected cyano gene expression, no MET")

V_cyano_inf_effect_met <- makeVolcano(
  res = cyano_inf_effect_met,
  title = "infected vs non infected cyano gene expression, MET")

V_cyano_met_effect_1org <- makeVolcano(
  res = cyano_met_effect_1org,
  title = "MET effect on cyano gene expression")

V_cyano_met_effect_2orgs <- makeVolcano(
  res = cyano_met_effect_2orgs,
  title = "MET effect on cyano infecting bacteria gene expression")

## Extreme log2 fold changes could also arise from issues related to:
### Low counts: If the gene has very low counts in one condition (e.g., zero or near-zero),
# even a small count in the other condition can lead to a high fold change.
### Normalization methods: Ensure appropriate normalization techniques are applied during 
# the analysis. Poor normalization can distort fold change calculations.

## open bigger window
dev.new(width = 15, height = 12)
pdf("../../figures/Fig1_chytrid_volc.pdf", width = 15, height = 15)
cowplot::plot_grid(V_chytrid_inf_effect_control$plot,
                   V_chytrid_inf_effect_met$plot,
                   V_chytrid_met_effect_1org$plot,
                   V_chytrid_met_effect_2orgs$plot)
dev.off()

dev.new(width = 15, height = 12)
pdf("../../figures/Fig2_cyano_volc.pdf", width = 15, height = 15)
cowplot::plot_grid(V_cyano_inf_effect_control$plot,
                   V_cyano_inf_effect_met$plot,
                   V_cyano_met_effect_1org$plot,
                   V_cyano_met_effect_2orgs$plot)
dev.off()

fullDEGTable = rbind(
  V_chytrid_inf_effect_control$signifGenes %>% mutate(group = "infection effect on chytrid gene expression in absence of MET"),
  V_chytrid_inf_effect_met$signifGenes %>% mutate(group = "infection effect on chytrid gene expression in presence of MET"),
  V_chytrid_met_effect_1org$signifGenes %>% mutate(group = "MET effect on chytrid gene expression in absence of cyanobacteria"),
  V_chytrid_met_effect_2orgs$signifGenes %>% mutate(group = "MET effect on chytrid gene expression in presence of cyanobacteria"),
  V_cyano_inf_effect_control$signifGenes %>% mutate(group = "infection effect on cyanobacteria gene expression in absence of MET"),
  V_cyano_inf_effect_met$signifGenes %>% mutate(group = "infection effect on cyanobacteria gene expression in presence of MET"),
  V_cyano_met_effect_1org$signifGenes %>% mutate(group = "MET effect on cyanobacteria gene expression in absence of chytrid"),
  V_cyano_met_effect_2orgs$signifGenes %>% mutate(group = "MET effect on cyanobacteria gene expression in presence of chytrid"))

fullDEGTable$geneName = rownames(fullDEGTable)
fullDEGTable=fullDEGTable[c("geneName", names(fullDEGTable)[!names(fullDEGTable) %in% "geneName"])]

write.table(x = fullDEGTable, "../../figures/TableS1_fullDEGTable.tsv", quote = F, sep = "\t", row.names = F)
