##############
# Data loading
source("R00_analysisRSEMcountmatrix.R")

##### Preparation for DE analysis #####

# (c) chytrid genes expressed in chytrid alone and infecting 

## Select only chytrid and both genes
RSEM_final_hope.gene_chytrid <- RSEM_final_hope.gene[RSEM_final_hope.gene$X %in% RSEM_final_hope.gene[
  RSEM_final_hope.gene$whichOrg %in% "in_chytrid_alone in_both_organisms" &
    RSEM_final_hope.gene$whichTranscriptome %in% "chytrid","X"],]

## Select only chytrid and both samples 
RSEM_final_hope.gene_chytrid = 
  RSEM_final_hope.gene_chytrid[grep("cyano", names(RSEM_final_hope.gene_chytrid),invert = T)]

## Clean
names(RSEM_final_hope.gene_chytrid)[1]=""
RSEM_final_hope.gene_chytrid=RSEM_final_hope.gene_chytrid[
  !names(RSEM_final_hope.gene_chytrid) %in% c("whichOrg", "whichTranscriptome")]

## Investigate outliers (my function)
makeClusterWGCNA(t(RSEM_final_hope.gene_chytrid[-1]))

## met_both_In11 and met_chy_Z10 outliers
RSEM_final_hope.gene_chytrid <- 
  RSEM_final_hope.gene_chytrid[!names(RSEM_final_hope.gene_chytrid) %in% c("met_chy_Z10", "met_both_In11")]

makeClusterWGCNA(t(RSEM_final_hope.gene_chytrid[-1])) # good

# Z6 is a technical replicate of Z5
# Z12 is a technical replicate of Z7
# In6 is a technical replicate of In5
# In12 is a technical replicate of In11 --> no cluster by replicates, so these are sequencing differences

## Remove lowly expressed genes for robustness across experiments 
# (the both experiments may have lower depth and not even detect these genes), 
## At least gene exp non null in half of each group inf/not inf
RSEM_final_hope.gene_chytrid = RSEM_final_hope.gene_chytrid[
  (rowSums(RSEM_final_hope.gene_chytrid[grep("chy", names(RSEM_final_hope.gene_chytrid))]>0)>
     ncol(RSEM_final_hope.gene_chytrid[grep("chy", names(RSEM_final_hope.gene_chytrid))])/2) &
    (rowSums(RSEM_final_hope.gene_chytrid[grep("both", names(RSEM_final_hope.gene_chytrid))]>0)>
       ncol(RSEM_final_hope.gene_chytrid[grep("both", names(RSEM_final_hope.gene_chytrid))])/2),]

makeClusterWGCNA(t(RSEM_final_hope.gene_chytrid[-1])) # more homogeneous

write.table(RSEM_final_hope.gene_chytrid, 
            "../../data/run_DESEQ2_Erika/RSEM_final_hope.gene_chytrid_03122024.matrix", 
            sep="\t", quote = F, row.names = F)

rownames(RSEM_final_hope.gene_chytrid) <- RSEM_final_hope.gene_chytrid[[1]]
RSEM_final_hope.gene_chytrid=RSEM_final_hope.gene_chytrid[-1]

names(RSEM_final_hope.gene_chytrid)
nrow(RSEM_final_hope.gene_chytrid) # 1192

# (f) cyanobacteria genes expressed in cyano alone and infected --> give to Erika for DESEq2 3420 genes
## Select only cyano and both genes
RSEM_final_hope.gene_cyano <- RSEM_final_hope.gene[RSEM_final_hope.gene$X %in% RSEM_final_hope.gene[
  RSEM_final_hope.gene$whichOrg %in% "in_both_organisms in_cyano_alone" &
    RSEM_final_hope.gene$whichTranscriptome %in% "cyano","X"],]

## Select only cyano and both samples 
RSEM_final_hope.gene_cyano = 
  RSEM_final_hope.gene_cyano[grep("chy", names(RSEM_final_hope.gene_cyano),invert = T)]

## Clean
names(RSEM_final_hope.gene_cyano)[1]=""
RSEM_final_hope.gene_cyano=RSEM_final_hope.gene_cyano[
  !names(RSEM_final_hope.gene_cyano) %in% c("whichOrg", "whichTranscriptome")]

## Investigate outliers (my function)
makeClusterWGCNA(t(RSEM_final_hope.gene_cyano[-1]))

## no outliers!

## Remove lowly expressed genes for robustness across experiments 
# (the both experiments may have lower depth and not even detect these genes), 
## At least gene exp non null in half of each group inf/not inf
RSEM_final_hope.gene_cyano = RSEM_final_hope.gene_cyano[
  (rowSums(RSEM_final_hope.gene_cyano[grep("cyano", names(RSEM_final_hope.gene_cyano))]>0)>
     ncol(RSEM_final_hope.gene_cyano[grep("cyano", names(RSEM_final_hope.gene_cyano))])/2) &
    (rowSums(RSEM_final_hope.gene_cyano[grep("both", names(RSEM_final_hope.gene_cyano))]>0)>
       ncol(RSEM_final_hope.gene_cyano[grep("both", names(RSEM_final_hope.gene_cyano))])/2),]

makeClusterWGCNA(t(RSEM_final_hope.gene_cyano[-1])) # more homogeneous

write.table(RSEM_final_hope.gene_cyano, 
            "../../data/run_DESEQ2_Erika/RSEM_final_hope.gene_cyano_03122024.matrix", 
            sep="\t", quote = F, row.names = F)

rownames(RSEM_final_hope.gene_cyano) <- RSEM_final_hope.gene_cyano[[1]]
RSEM_final_hope.gene_cyano=RSEM_final_hope.gene_cyano[-1]

names(RSEM_final_hope.gene_cyano)

nrow(RSEM_final_hope.gene_cyano) # 2166

#################
# DESeq2 analysis- for this analysis I need un-normalized counts. 
# I need to have integers. My values are numeric, but According to a blog in which 
# someone asked about the same problem, Michael Love said: Just round the counts to integer. 
# There isn't a loss in precision to this operation (sampling variance from the experiment
# is much larger). https://support.bioconductor.org/p/105964/

calculateContrasts <- function(my_countsmatrix, my_org, my_samples=samples_data){
  
  ## Define our 4 conditions, we are comparing 2 by 2
  my_groups=c("control_both","met_both", paste0("control_",my_org), paste0("met_",my_org))
  
  my_samples=my_samples[my_samples$sample %in% names(my_countsmatrix),]
  
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
  my_countsmatrix=RSEM_final_hope.gene_cyano[-1],
  my_org="cyano")

contrast_chytridgenome <- calculateContrasts(
  my_countsmatrix=RSEM_final_hope.gene_chytrid[-1],
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

plotHeatmap(contrast_cyanogenome$vstr)
plotHeatmap(contrast_chytridgenome$vstr)

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

makeCorrelationPlot(contrast_cyanogenome$vstr)
makeCorrelationPlot(contrast_chytridgenome$vstr)

##################################################
## Rename DESEq2 files with meaningfull gene names
renameDESeq <- function(DESEq2, annotation){
  DESEq2$custom_gene_name <- row.names(DESEq2)
  row.names(DESEq2) = annotation$gene_name[
    match(row.names(DESEq2), annotation$custom_gene_name)]
  return(DESEq2)
}

DESEq2 = contrast_cyanogenome$resr_inf_effect_control
annotation = annotationCyano

DESEq2$custom_gene_name <- row.names(DESEq2)
annotation$gene_name[
  match(row.names(DESEq2), annotation$custom_gene_name)]

chytrid_inf_effect_control <- renameDESeq(contrast_chytridgenome$resr_inf_effect_control, 
                                          annotation = annotationChytrid)
chytrid_inf_effect_met <- renameDESeq(contrast_chytridgenome$resr_inf_effect_met, 
                                      annotation = annotationChytrid)
chytrid_met_effect_1org <- renameDESeq(contrast_chytridgenome$resr_met_effect_1org, 
                                       annotation = annotationChytrid)
chytrid_met_effect_2orgs <- renameDESeq(contrast_chytridgenome$resr_met_effect_2orgs, 
                                        annotation = annotationChytrid)

cyano_inf_effect_control <- renameDESeq(contrast_cyanogenome$resr_inf_effect_control, 
                                        annotation = annotationCyano)
cyano_inf_effect_met <- renameDESeq(contrast_cyanogenome$resr_inf_effect_met, 
                                    annotation = annotationCyano)
cyano_met_effect_1org <- renameDESeq(contrast_cyanogenome$resr_met_effect_1org, 
                                     annotation = annotationCyano)
cyano_met_effect_2orgs <- renameDESeq(contrast_cyanogenome$resr_met_effect_2orgs, 
                                      annotation = annotationCyano)

mylistResDESEQ2 <- list(chytrid_inf_effect_control=chytrid_inf_effect_control,
                        chytrid_inf_effect_met=chytrid_inf_effect_met,
                        chytrid_met_effect_1org=chytrid_met_effect_1org,
                        chytrid_met_effect_2orgs=chytrid_met_effect_2orgs,
                        cyano_inf_effect_control=cyano_inf_effect_control,
                        cyano_inf_effect_met=cyano_inf_effect_met,
                        cyano_met_effect_1org=cyano_met_effect_1org,
                        cyano_met_effect_2orgs=cyano_met_effect_2orgs)


################
## group d: 124 cyano genes only expressed outside of infection,
## at least in half of each group met/not met

## make a ddsr object comparing between conditions
ddsr <- DESeqDataSetFromMatrix(countData = round(RSEM_final_hope.gene_d),
                               colData = samples_data[samples_data$sample %in% names(RSEM_final_hope.gene_d),],
                               design = ~ condition)
ddsr <- DESeq(ddsr)

# Variance stabilising transformation
vstr <- as.data.frame(assay(varianceStabilizingTransformation(ddsr)))
vstr$Gene <- rownames(vstr)

## Calculate the contrasts (pairwise comparisons of interest)
resd <- results(ddsr, alpha = .05)
table(resd$padj < 0.05) # no DEG
