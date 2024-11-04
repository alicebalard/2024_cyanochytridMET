library(DESeq2)
library(reshape2)
library("viridis")
library("pheatmap")
library("RColorBrewer")
library(ggplot2)

# Set path
setwd("/Users/akire/Documents/IGB-tasks/Transcriptome_assay/diff_exp/data")

##############
# Data loading
## case 1- one alignment 
### Load the count matrix
count_matrixr <- read.table("RSEM.gene.counts.matrix", header=T)

### Load all the samples
samples_data <- read.table("/Users/akire/Documents/IGB-tasks/Transcriptome_assay/diff_exp/required_tables/sample_data_r.txt", header = TRUE, sep = "\t")

## Split by: everything with chytrid + everything with cyano
samples_chytrid_genome <- samples_data[grep("chy|both", samples_data$condition),]
count_matrixr_chytridgenome <- count_matrixr[grep("cyano", row.names(count_matrixr),invert = TRUE ),]
count_matrixr_chytridgenome <- count_matrixr_chytridgenome[grep("chy|both", names(count_matrixr_chytridgenome))]

samples_cyano_genome <- samples_data[grep("cyano|both", samples_data$condition),]
count_matrixr_cyanogenome <- count_matrixr[grep("cyano", row.names(count_matrixr)),]
count_matrixr_cyanogenome <- count_matrixr_cyanogenome[grep("cyano|both", names(count_matrixr_cyanogenome))]

## case 2- four alignments
# RSEM_cyano _chytri _both_cyano _both_chytrid
datar_chytrid <- read.table("RSEM_chytrid.gene.counts.matrix", header=T)
datar_both_chytrid <- read.table("RSEM_both_chytrid.gene.counts.matrix", header=T)
count_matrixr_chytridgenome_CASE2 <-cbind(datar_chytrid, datar_both_chytrid)

datar_cyano <- read.table("RSEM_cyano.gene.counts.matrix", header=T)
datar_both_cyano <- read.table("RSEM_both_cyano.gene.counts.matrix", header=T)
count_matrixr_cyanogenome_CASE2 <- cbind(datar_cyano, datar_both_cyano)

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

contrast_cyanogenome_Case1 <- calculateContrasts(
  my_samples=samples_cyano_genome,
  my_countsmatrix=count_matrixr_cyanogenome,
  my_org="cyano")

contrast_chytridgenome_Case1 <- calculateContrasts(
  my_samples=samples_chytrid_genome,
  my_countsmatrix=count_matrixr_chytridgenome,
  my_org="chy")

contrast_cyanogenome_Case2 <- calculateContrasts(
  my_samples=samples_cyano_genome,
  my_countsmatrix=count_matrixr_cyanogenome_CASE2,
  my_org="cyano")

contrast_chytridgenome_Case2 <- calculateContrasts(
  my_samples=samples_chytrid_genome,
  my_countsmatrix=count_matrixr_chytridgenome_CASE2,
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
plotHeatmap(contrast_cyanogenome_Case2$vstr)
plotHeatmap(contrast_chytridgenome_Case1$vstr)
plotHeatmap(contrast_chytridgenome_Case2$vstr)

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