library(DESeq2)
library(reshape2)
library("viridis")
library("pheatmap")
library("RColorBrewer")
library(ggplot2)

##############
# Data loading

## Chytrid annotation
annotChyt <- read.csv("../../data/allFungiTrinot_simplified.tsv", sep ="\t")
## rm multiple lines due to multiple sources of annotations
annotChyt <- unique(annotChyt)

## extract GO terms in their own column for GO analysis downstream
annotChyt$GO <- 
  sapply(annotChyt$gene_ontology_BLASTX, function(x) {
    go_terms <- str_extract_all(x, "GO:\\d+")
    paste(unlist(go_terms), collapse = ",")
  })

annotChytFungi <- annotChyt[grep("Fungi", annotChyt$sprot_Top_BLASTX_hit),]

### Load all the samples
samples_data <- read.table(
  "../../data/sample_data_remove_r.txt", header = TRUE, sep = "\t")

### Load the count matrix
count_matrixr <- read.csv(
  "../../data/RSEM_remove.gene.counts.matrix", sep = "\t")

#############
## Split by: everything with chytrid + everything with cyano
samples_chytrid_genome <- samples_data[grep("chy|both", samples_data$condition),]

### Select chytrid genes
count_matrixr_chytridgenome <- count_matrixr[grep("cyano", count_matrixr$X, invert = T),]

## rename chytrid genes
### gene number and name without any additional information (e.g. gene 1 TRINITY_DN43456_c0_g1_i1).
new_gene_trans_map_fungi <- read.table(
  "../../data/new_gene_trans_map_fungi.txt")
new_gene_trans_map_fungi$V3 <- annotChyt$gene_name[
  match(new_gene_trans_map_fungi$V2, annotChyt$transcript_id)]

count_matrixr_chytridgenome$isoform <- new_gene_trans_map_fungi$V2[
  match(count_matrixr_chytridgenome$X, new_gene_trans_map_fungi$V1)]

count_matrixr_chytridgenome$geneName <- new_gene_trans_map_fungi$V3[
  match(count_matrixr_chytridgenome$X, new_gene_trans_map_fungi$V1)]

head(count_matrixr_chytridgenome)

### SUBSELECT ONLY FUNGI (to be corrected earlier)
count_matrixr_chytridgenome <- count_matrixr_chytridgenome[
  count_matrixr_chytridgenome$X %in% annotChytFungi$gene_name,]


table(count_matrixr_chytridgenome$X %in% ".")

count_matrixr_chytridgenome[
  count_matrixr_chytridgenome$X %in% ".",]


## Select samples chytrid and both
count_matrixr_chytridgenome <- count_matrixr_chytridgenome[grep("chy|both", names(count_matrixr_chytridgenome))]

###################
## 2. cyanobacteria (tbc)
samples_cyano_genome <- samples_data[grep("cyano|both", samples_data$condition),]
count_matrixr_cyanogenome <- count_matrixr[grep("cyano", row.names(count_matrixr)),]
count_matrixr_cyanogenome <- count_matrixr_cyanogenome[grep("cyano|both", names(count_matrixr_cyanogenome))]

#################
# DESeq2 analysis- for this analysis I need un-normalized counts. 
# I need to have integers. My values are numeric, but According to a blog in which someone asked about the same problem, Michael Love said: Just round the counts to integer. There isn't a loss in precision to this operation (sampling variance from the experiment is much larger). https://support.bioconductor.org/p/105964/

calculateContrasts <- function(my_samples, my_countsmatrix, my_org){
  ## Define our 4 conditions, we are comparing 2 by 2
  my_groups=c("control_both","met_both", paste0("control_",my_org), paste0("met_",my_org))
  
  ## make a ddsr object comparing between conditions
  ddsr <- DESeqDataSetFromMatrix(countData = round(my_countsmatrix),
                                 colData = my_samples,
                                 design = ~ condition)
  ddsr <- DESeq(ddsr)
  
  # Variance stabilising transformation
  vstr <- as.data.frame(assay(vst(ddsr)))
  vstr$Gene <- rownames(vstr)
  
  ## Calculate the contrasts (pairwise comparisons of interest)
  print(paste0("comparison of groups: ", my_groups[1], " vs ", my_groups[2]))
  resr_met_effect_2orgs <- results(ddsr, contrast=c("condition", my_groups[1], my_groups[2]), alpha=0.05) 
  print(summary(resr_met_effect_2orgs))
  print(sum(resr_met_effect_2orgs$padj < 0.05, na.rm=TRUE))
  
  print(paste0("comparison of groups: ", my_groups[3]," vs ", my_groups[4]))
  resr_met_effect_1org <- results(ddsr, contrast=c("condition", my_groups[3], my_groups[4]), alpha=0.05) 
  print(summary(resr_met_effect_1org))
  print(sum(resr_met_effect_1org$padj < 0.05, na.rm=TRUE))
  
  print(paste0("comparison of groups: ", my_groups[3], " vs ", my_groups[1]))
  resr_inf_effect_control <- results(ddsr, contrast=c("condition", my_groups[3], my_groups[1]), alpha=0.05) 
  print(summary(resr_inf_effect_control))
  print(sum(resr_inf_effect_control$padj < 0.05, na.rm=TRUE))
  
  print(paste0("comparison of groups: ", my_groups[4], " vs ", my_groups[2]))
  resr_inf_effect_met <- results(ddsr, contrast=c("condition", my_groups[4], my_groups[2]), alpha=0.05) 
  print(summary(resr_inf_effect_met))
  print(sum(resr_inf_effect_met$padj < 0.05, na.rm=TRUE))
  
  # plotMA
  plotMA(resr_met_effect_2orgs)
  plotMA(resr_met_effect_1org)
  plotMA(resr_inf_effect_control)
  plotMA(resr_inf_effect_met)
  
  return(list(resr_met_effect_2orgs=resr_met_effect_2orgs, 
              resr_met_effect_1org=resr_met_effect_1org,
              resr_inf_effect_control=resr_inf_effect_control,
              resr_inf_effect_met=resr_inf_effect_met,
              vstr=vstr))
}

contrast_cyanogenome <- calculateContrasts(
  my_samples=samples_cyano_genome,
  my_countsmatrix=count_matrixr_cyanogenome,
  my_org="cyano")

contrast_chytridgenome <- calculateContrasts(
  my_samples=samples_chytrid_genome,
  my_countsmatrix=count_matrixr_chytridgenome,
  my_org="chy")

plotHeatmap <- function(vstr){
  # Overwrite our original data frame with the long format
  deseq2VSTr <- melt(vstr, id.vars=c("Gene"))
  
  # Plot heatmap
  heatmap2 <- ggplot(deseq2VSTr, aes(x=variable, y=Gene, fill=value)) + 
    geom_raster() + 
    scale_fill_viridis(trans="sqrt") + 
    theme(axis.text.x=element_text(angle=65, hjust=1),
          axis.text.y=element_blank(), axis.ticks.y=element_blank())
  heatmap2
  return(heatmap2)
}

plotHeatmap(contrast_cyanogenome_Case1$vstr)
plotHeatmap(contrast_chytridgenome_Case1$vstr)


# heatmap of the sample-to-sample distances

makeCorrelationPlot <- function(my_vstr){
  my_vstr=my_vstr[!names(my_vstr)%in%"Gene"]
  
  sampleDists <- dist(t(my_vstr))
  
  sampleDistMatrix <- as.matrix(sampleDists)
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  pheatmap(sampleDistMatrix,
           clustering_distance_rows=sampleDists,
           clustering_distance_cols=sampleDists,
           col=colors)
}

makeCorrelationPlot(contrast_cyanogenome_Case1$vstr)
#makeCorrelationPlot(contrast_cyanogenome_Case2$vstr)
makeCorrelationPlot(contrast_chytridgenome_Case1$vstr)
#makeCorrelationPlot(contrast_chytridgenome_Case2$vstr)