## 29th of November 2024
source("libLoad.R")
source("dataLoad.R")

## Let's observe the count matrix calculated by Trinity
RSEM_final_hope.gene <- 
  read.csv("../../data/run_DESEQ2_Erika/RSEM_final_hope.gene.counts.matrix", sep="\t")
rownames(RSEM_final_hope.gene)=RSEM_final_hope.gene$X
## 9766 genes

#######################################################################
## Split by group depending on which gene is expressed in which case ##
#######################################################################
a = ifelse(rowSums(RSEM_final_hope.gene[grep("chy", names(RSEM_final_hope.gene))]) !=0,
           "in_chytrid_alone", "")
b = ifelse(rowSums(RSEM_final_hope.gene[grep("both", names(RSEM_final_hope.gene))]) !=0,
           "in_both_organisms", "")
c = ifelse(rowSums(RSEM_final_hope.gene[grep("cyano", names(RSEM_final_hope.gene))]) !=0,
           "in_cyano_alone", "")

RSEM_final_hope.gene$whichOrg <- trimws(paste(a,b,c, sep = " "))

RSEM_final_hope.gene$whichTranscriptome <- ifelse(
  grepl("TRINITY", RSEM_final_hope.gene$X), "chytrid", "cyano")

table(RSEM_final_hope.gene$whichTranscriptome,
      RSEM_final_hope.gene$whichOrg)
#            0  in_both_organisms in_both_organisms in_cyano_alone in_chytrid_alone
# chytrid  679               1764                              751              360
# cyano    629                  1                             3420                0
# 
#         in_chytrid_alone  in_cyano_alone in_chytrid_alone in_both_organisms
# chytrid                                5                               1623
# cyano                                  2                                  0
# 
#         in_chytrid_alone in_both_organisms in_cyano_alone in_cyano_alone
# chytrid                                               140             16
# cyano                                                 161            215

## Figure X ##

##### 1. rm contamination in my chytrid assembly #####
# chytrid genes also found in when cyanobacteria is alone:
# chytrid + in_both_organisms in_cyano_alone
# chytrid + in_chytrid_alone  in_cyano_alone
# chytrid + in_chytrid_alone in_both_organisms in_cyano_alone
# chytrid + in_cyano_alone

listOfTranscriptContaminant_toRmFromChytridTranscriptome <-
  RSEM_final_hope.gene[
    RSEM_final_hope.gene$whichTranscriptome %in% "chytrid" & RSEM_final_hope.gene$whichOrg %in% 
      c("in_both_organisms in_cyano_alone",
        "in_chytrid_alone  in_cyano_alone",
        "in_chytrid_alone in_both_organisms in_cyano_alone",
        "in_cyano_alone"),"X"]

write.csv(listOfTranscriptContaminant_toRmFromChytridTranscriptome,
  "../../data/listOfTranscriptContaminant_toRmFromChytridTranscriptome", quote = F, row.names = F)

#### 2. Which genes have been sequenced for chytrid and cyano (and are not contamination)?
#############
## Chytrid ##
sequencedChytridGenes <- RSEM_final_hope.gene[
  RSEM_final_hope.gene$whichOrg %in% c("in_chytrid_alone in_both_organisms",
                                       "in_chytrid_alone", "in_both_organisms") &
    RSEM_final_hope.gene$whichTranscriptome %in% "chytrid","X"]
length(sequencedChytridGenes) # 3747 sequenced genes

## Select only "chytrid and both" genes
# 3747 genes investigated
RSEM_final_hope.gene_chytrid <- RSEM_final_hope.gene[
  RSEM_final_hope.gene$X %in% sequencedChytridGenes,]

## Select only "chytrid and both" samples 
RSEM_final_hope.gene_chytrid = RSEM_final_hope.gene_chytrid[
  grep("cyano", names(RSEM_final_hope.gene_chytrid),invert = T)] 

## Clean
RSEM_final_hope.gene_chytrid=RSEM_final_hope.gene_chytrid[
  !names(RSEM_final_hope.gene_chytrid) %in% c("X", "whichOrg", "whichTranscriptome")]

## Rename based on annotations
rownames(RSEM_final_hope.gene_chytrid) = make.unique(annotationChytrid$gene_name[
  match(row.names(RSEM_final_hope.gene_chytrid), annotationChytrid$custom_gene_name)])

nrow(RSEM_final_hope.gene_chytrid) # 3747
## Merge identical proteins in only one row, suming the counts
RSEM_final_hope.gene_chytrid = RSEM_final_hope.gene_chytrid %>%
  mutate(base_name = sub("\\.\\d+$", "", rownames(RSEM_final_hope.gene_chytrid))) %>%
  group_by(base_name) %>%
  summarise(across(everything(), sum)) %>%
  tibble::column_to_rownames("base_name") %>% data.frame()

nrow(RSEM_final_hope.gene_chytrid) # 3156 genes

# 3747 - 3156 = 591 genes that had multiple transcript names

###################
## Cyanobacteria ##
sequencedCyanoGenes <- RSEM_final_hope.gene[
  RSEM_final_hope.gene$whichOrg %in% c("in_both_organisms in_cyano_alone",
                                       "in_cyano_alone", "in_both_organisms") &
    RSEM_final_hope.gene$whichTranscriptome %in% "cyano","X"]
length(sequencedCyanoGenes) # 3636 sequenced genes

## Select only "cyano and both" genes
RSEM_final_hope.gene_cyano <- RSEM_final_hope.gene[
  RSEM_final_hope.gene$X %in% sequencedCyanoGenes,]

## Select only "cyano and both" samples 
RSEM_final_hope.gene_cyano = RSEM_final_hope.gene_cyano[
  grep("chy", names(RSEM_final_hope.gene_cyano),invert = T)] 

## Clean
RSEM_final_hope.gene_cyano=RSEM_final_hope.gene_cyano[
  !names(RSEM_final_hope.gene_cyano) %in% c("X", "whichOrg", "whichTranscriptome")]

## Rename based on annotations
rownames(RSEM_final_hope.gene_cyano) = make.unique(annotationCyano$gene_name[
  match(row.names(RSEM_final_hope.gene_cyano), annotationCyano$custom_gene_name)])

nrow(RSEM_final_hope.gene_cyano) # 3636 genes

## Merge identical proteins in only one row, suming the counts
RSEM_final_hope.gene_cyano = RSEM_final_hope.gene_cyano %>%
  mutate(base_name = sub("\\.\\d+$", "", rownames(RSEM_final_hope.gene_cyano))) %>%
  group_by(base_name) %>%
  summarise(across(everything(), sum)) %>%
  tibble::column_to_rownames("base_name") %>% data.frame()

nrow(RSEM_final_hope.gene_cyano) # 3589 genes

# 3636-3589=47 genes with duplicated counts

###################################
## Low quality samples filtering ##
###################################

## Investigate outliers (my function)
makeClusterWGCNA(t(RSEM_final_hope.gene_chytrid))
## In11 clear outlier
makeClusterWGCNA(t(RSEM_final_hope.gene_cyano))
# no clear outlier

#####!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
## Rm In11 from the chytrid dataset
RSEM_final_hope.gene_chytrid = RSEM_final_hope.gene_chytrid[
  !names(RSEM_final_hope.gene_chytrid) %in% "met_both_In11"]
## rm genes with zero counts (196 genes of chytrid were only in In11!!)
RSEM_final_hope.gene_chytrid = RSEM_final_hope.gene_chytrid[
  rowSums(RSEM_final_hope.gene_chytrid) != 0,]

# To determine if a sample has been sufficiently sequenced, a saturation curve
# can be generated with software like vegan in R

# if (!requireNamespace("BiocManager", quietly = TRUE))
  # install.packages("BiocManager")
# BiocManager::install("NOISeq")

library(NOISeq)

# Sequencing depth & Expression Quantification
# The plots in this section can be generated by only providing the expression data, since no other biological
# information is required. Their purpose is to assess if the sequencing depth of the samples is enough to detect the
# features of interest and to get a good quantification of their expression.

# Saturation plot
# The “Saturation" plot shows the number of features in the genome detected with more than k counts with
# the sequencing depth of the sample, and with higher and lower simulated sequencing depth
# Create an ExpressionSet object
counts_chy <- as.matrix(RSEM_final_hope.gene_chytrid)
eset_chy <- ExpressionSet(assayData = counts_chy)
mysaturation_chy = dat(eset_chy, k = 0, ndepth = 7, type = "saturation")
## Chytrid alone:
explo.plot(mysaturation_chy, toplot = 1, samples = 1:11)
## Saturation of chytrid low (between 29.7% and 44.1%)
## Both organisms:
explo.plot(mysaturation_chy, toplot = 1, samples = 12:18)
## 34.4 to 50.1%

counts_cyano <- as.matrix(RSEM_final_hope.gene_cyano)
eset_cyano <- ExpressionSet(assayData = counts_cyano)
mysaturation_cyano = dat(eset_cyano, k = 0, ndepth = 7, type = "saturation")

## Both organisms:
explo.plot(mysaturation_cyano, toplot = 1, samples = 1:8) 
# good (60+) except for In2 and In11 (<35%)
## Cyano alone
explo.plot(mysaturation_cyano, toplot = 1, samples = 9:20) 
# very good saturation level >97.3%

#####!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
## Rm In2 and In11 from the cyano dataset
RSEM_final_hope.gene_cyano = RSEM_final_hope.gene_cyano[
  !names(RSEM_final_hope.gene_cyano) %in% c("control_both_In2", "met_both_In11")]
## rm genes with zero counts
RSEM_final_hope.gene_cyano = RSEM_final_hope.gene_cyano[
  rowSums(RSEM_final_hope.gene_cyano) != 0,]

###################################
## Low quality genes filtering ##
###################################

RSEM_final_hope.gene_chytrid <- RSEM_final_hope.gene_chytrid[
  rowSums(RSEM_final_hope.gene_chytrid >= 10) >= 3,]
RSEM_final_hope.gene_cyano <- RSEM_final_hope.gene_cyano[
  rowSums(RSEM_final_hope.gene_cyano >= 10) >= 3,]

nrow(RSEM_final_hope.gene_chytrid) # 1336
nrow(RSEM_final_hope.gene_cyano) # 3497

#####################
## DESeq2 jan 2025 ##
#####################

## I. Chytrid transcripts
contrast_chytridgenome<- calculateContrasts(
  my_countsmatrix=RSEM_final_hope.gene_chytrid,
  my_org="chy")

## II. Cyano transcripts
contrast_cyanogenome <- calculateContrasts(
  my_countsmatrix=RSEM_final_hope.gene_cyano,
  my_org="cyano")

## Volcano plots
V_chytrid_inf_effect_control <- makeVolcano(
  res = contrast_chytridgenome$resr_inf_effect_control,
  title = "infecting vs non infecting chytrid gene expression, no MET")

V_chytrid_inf_effect_met <- makeVolcano(
  res = contrast_chytridgenome$resr_inf_effect_met,
  title = "infecting vs non infecting chytrid gene expression, MET")

V_chytrid_met_effect_1org <- makeVolcano(
  res = contrast_chytridgenome$resr_met_effect_1org,
  title = "MET effect on chytrid gene expression")

V_chytrid_met_effect_2orgs <- makeVolcano(
  res = contrast_chytridgenome$resr_met_effect_2orgs,
  title = "MET effect on chytrid infecting bacteria gene expression")

V_cyano_inf_effect_control <- makeVolcano(
  res = contrast_cyanogenome$resr_inf_effect_control,
  title = "infected vs non infected cyano gene expression, no MET")

V_cyano_inf_effect_met <- makeVolcano(
  res = contrast_cyanogenome$resr_inf_effect_met,
  title = "infected vs non infected cyano gene expression, MET")

V_cyano_met_effect_1org <- makeVolcano(
  res = contrast_cyanogenome$resr_met_effect_1org,
  title = "MET effect on cyano gene expression")

V_cyano_met_effect_2orgs <- makeVolcano(
  res = contrast_cyanogenome$resr_met_effect_2orgs,
  title = "MET effect on cyano infecting bacteria gene expression")

## open bigger window
dev.new(width = 15, height = 12)
pdf("../../figures/Fig1-part1_chytrid_volc.pdf", width = 15, height = 15)
cowplot::plot_grid(V_chytrid_inf_effect_control$plot,
                   V_chytrid_inf_effect_met$plot,
                   V_chytrid_met_effect_1org$plot,
                   V_chytrid_met_effect_2orgs$plot)
dev.off()

dev.new(width = 15, height = 12)
pdf("../../figures/Fig2-part1_cyano_volc.pdf", width = 15, height = 15)
cowplot::plot_grid(V_cyano_inf_effect_control$plot,
                   V_cyano_inf_effect_met$plot,
                   V_cyano_met_effect_1org$plot,
                   V_cyano_met_effect_2orgs$plot)
dev.off()

## Venn diagrams
# devtools::install_github("yanlinlin82/ggvenn")
library(ggvenn)

## Select differentially expressed genes in each comparison
getGenes <- function(x) {rownames(x[x$padj < 0.05 & !is.na(x$padj),])}

# Combine the lists into a named list for the Venn diagram
venn_data <- list("infection effect chytrid - control"= 
                    getGenes(contrast_chytridgenome$resr_inf_effect_control),
                  "met effect chytrid - free living"= 
                    getGenes(contrast_chytridgenome$resr_met_effect_1org),
                  "met effect chytrid - infecting"= 
                    getGenes(contrast_chytridgenome$resr_met_effect_2orgs),
                  "infection effect chytrid - met"= 
                    getGenes(contrast_chytridgenome$resr_inf_effect_met))
  
pdf("../../figures/Fig1-part2_Venn.pdf", width = 5, height = 5)
ggvenn(
  venn_data, show_stats = "c", fill_color = rep("white", 4),
  stroke_size = 0.5, set_name_size = 4
)
dev.off()

# Combine the lists into a named list for the Venn diagram
venn_data <- list("infection effect cyano - control"= 
                    getGenes(contrast_cyanogenome$resr_inf_effect_control),
                  "met effect cyano - alone"= 
                    getGenes(contrast_cyanogenome$resr_met_effect_1org),
                  "met effect cyano - infected"= 
                    getGenes(contrast_cyanogenome$resr_met_effect_2orgs),
                  "infection effect cyano - met"= 
                    getGenes(contrast_cyanogenome$resr_inf_effect_met))

pdf("../../figures/Fig2-part2_Venn.pdf", width = 5, height = 5)
ggvenn(
  venn_data, show_stats = "c", fill_color = rep("white", 4),
  stroke_size = 0.5, set_name_size = 4
)
dev.off()

## Save results in tables
selectDEGenes <- function(x) {x[x$padj < 0.05 & !is.na(x$padj),]}

## 1. chytrid
contrast_chytridgenome_DEG <- lapply(contrast_chytridgenome[1:4], selectDEGenes)

x <- list()
for (i in 1:4) {
  if (nrow(contrast_chytridgenome_DEG[[i]]) > 0) {  # Check if the data frame is not empty
    x[[i]] <- data.frame(
      geneName = rownames(contrast_chytridgenome_DEG[[i]]),
      padj = contrast_chytridgenome_DEG[[i]]$padj,
      log2FoldChange = contrast_chytridgenome_DEG[[i]]$log2FoldChange
    )
    x[[i]]$comparison <- names(contrast_chytridgenome_DEG)[i]
  } else {
    x[[i]] <- data.frame()  # Create an empty data frame if no rows are present
  }
}

contrast_chytridgenome_DEG <- do.call(rbind, x)

contrast_chytridgenome_DEG$comparison[
  contrast_chytridgenome_DEG$comparison %in% "resr_inf_effect_control"] <-
  "Infection effect on chytrid gene expression, in the absence of MET"
contrast_chytridgenome_DEG$comparison[
  contrast_chytridgenome_DEG$comparison %in% "resr_inf_effect_met"] <-
  "Infection effect on chytrid gene expression, in the presence of MET "
contrast_chytridgenome_DEG$comparison[
  contrast_chytridgenome_DEG$comparison %in% "resr_met_effect_1org"] <-
  "MET effect on chytrid gene expression, in free-living chytrid zoospores"
contrast_chytridgenome_DEG$comparison[
  contrast_chytridgenome_DEG$comparison %in% "resr_met_effect_2orgs"] <-
  "MET effect on chytrid gene expression, in infecting chytrid"

table(contrast_chytridgenome_DEG$comparison)

## 2. cyanobacteria
contrast_cyanogenome_DEG <- lapply(contrast_cyanogenome[1:4], selectDEGenes)

x <- list()
for (i in 1:4) {
  if (nrow(contrast_cyanogenome_DEG[[i]]) > 0) {  # Check if the data frame is not empty
    x[[i]] <- data.frame(
      geneName = rownames(contrast_cyanogenome_DEG[[i]]),
      padj = contrast_cyanogenome_DEG[[i]]$padj,
      log2FoldChange = contrast_cyanogenome_DEG[[i]]$log2FoldChange
    )
    x[[i]]$comparison <- names(contrast_cyanogenome_DEG)[i]
  } else {
    x[[i]] <- data.frame()  # Create an empty data frame if no rows are present
  }
}

contrast_cyanogenome_DEG <- do.call(rbind, x)

contrast_cyanogenome_DEG$comparison[
  contrast_cyanogenome_DEG$comparison %in% "resr_inf_effect_control"] <-
  "Infection effect on cyanobacteria gene expression, in the absence of MET"
contrast_cyanogenome_DEG$comparison[
  contrast_cyanogenome_DEG$comparison %in% "resr_inf_effect_met"] <-
  "Infection effect on cyanobacteria gene expression, in the presence of MET "
contrast_cyanogenome_DEG$comparison[
  contrast_cyanogenome_DEG$comparison %in% "resr_met_effect_1org"] <-
  "MET effect on cyanobacteria gene expression, in uninfected cyanobacteria"
contrast_cyanogenome_DEG$comparison[
  contrast_cyanogenome_DEG$comparison %in% "resr_met_effect_2orgs"] <-
  "MET effect on cyanobacteria gene expression, in infected bacteria"

table(contrast_cyanogenome_DEG$comparison)

fullDEGTable <- rbind(contrast_chytridgenome_DEG, contrast_cyanogenome_DEG)
  
fullDEGTable$padj <- signif(fullDEGTable$padj, 2)
fullDEGTable$log2FoldChange <- signif(fullDEGTable$log2FoldChange, 2)

write.csv(fullDEGTable, "../../figures/TableS1_fullDEGTable.tsv", row.names = F)        

## TBC from here


#############
## Chytrid ##
a=getGenes(contrast_chytridgenome$resr_met_effect_1org)
b=getGenes(contrast_chytridgenome$resr_met_effect_2orgs)
c=getGenes(contrast_chytridgenome$resr_inf_effect_control)
d=getGenes(contrast_chytridgenome$resr_inf_effect_met)

## Metolachlor effect on free-living AND infecting chytrid: 3 genes ##
intersect(a, b)[!intersect(a, b) %in% union(c, d)]
# 3 catabolism linked protein
# "PCCB_PIG"    "PDX1_DICDI"  "PRP16_ARATH"

##################################################################
## Infection effect on control AND metolachlor chytrid: 6 genes ##
intersect(c, d)[!intersect(c, d) %in% union(a, b)]
# "APC1_DICDI"  "ERT1_USTMA"  "LORF2_MOUSE" "NEP1_YEAST"  "SCP36_ARATH" "XDH_DICDI"  

#################
## GO analysis ##

## take into account filtration
universe_chytrid = rownames(RSEM_final_hope.gene_chytrid[rowSums(RSEM_final_hope.gene_chytrid >= 10) >= 3,])
universe_cyano = rownames(RSEM_final_hope.gene_cyano[rowSums(RSEM_final_hope.gene_cyano >= 10) >= 3,])
## done before!


getGOBubbleZ(universe = universe_chytrid, annotation = annotationChytrid, isDE = F, 
             genelist = a, GO_df = GO_chytrid, isbubble = F)
getGOBubbleZ(universe = universe_chytrid, annotation = annotationChytrid, isDE = F, 
             genelist = b, GO_df = GO_chytrid, isbubble = F)
getGOBubbleZ(universe = universe_chytrid, annotation = annotationChytrid, isDE = F, 
             genelist = c, GO_df = GO_chytrid, isbubble = F)
getGOBubbleZ(universe = universe_chytrid, annotation = annotationChytrid, isDE = F, 
             genelist = d, GO_df = GO_chytrid, isbubble = F)

# no significant GO terms

###################
## Cyanobacteria ##
a=getGenes(contrast_cyanogenome$resr_met_effect_1org) # empty
b=getGenes(contrast_cyanogenome$resr_met_effect_2orgs)
c=getGenes(contrast_cyanogenome$resr_inf_effect_control)
d=getGenes(contrast_cyanogenome$resr_inf_effect_met)

## metolachlor effect on cyanobacteria, infected of not -> nul
intersect(a, b)[!intersect(a, b) %in% union(c, d)]

## Infection effect on control AND metolachlor cyano: 217 genes ##
intersect(c, d)[!intersect(c, d) %in% union(a, b)]

#################
## GO analysis ##
getGOBubbleZ(universe = universe_cyano, annotation = annotationCyano, isDE = F, 
             genelist = b, GO_df = GO_cyano, isbubble = F)
getGOBubbleZ(universe = universe_cyano, annotation = annotationCyano, isDE = F, 
             genelist = c, GO_df = GO_cyano, isbubble = F)
# no significant GO terms

d_cyano_GO = getGOBubbleZ(universe = universe_cyano, annotation = annotationCyano, isDE = F, 
             genelist = d, GO_df = GO_cyano, isbubble = F)

d_cyano_GO$enrichment[d_cyano_GO$enrichment$p.adjust < 0.05,]
# DNA-binding transcription factor activity GO:0003700 gene ratio 0.02807018 12/1450 p.adjust 0.01254802

##########################################
## Top genes expressed in each category ##
contrast_cyanogenome$resr_met_effect_2orgs$padj < 0.05

x=contrast_cyanogenome$resr_met_effect_2orgs
x=x[!is.na(x$padj),]
head(x[order(x$log2FoldChange),], n=10)
tail(x[order(x$log2FoldChange),], n=10)

############ PREV
RSEM_final_hope.gene_chytrid[grepl("YI31B_YEAST", rownames(RSEM_final_hope.gene_chytrid)),]


head(contrast_chytridgenome$resr_met_effect_2orgs) #let's look at the results table


table(contrast_chytridgenome$resr_met_effect_2orgs$padj < 0.05 &
        contrast_chytridgenome$resr_met_effect_2orgs$log2FoldChange)


result <- sapply(contrast_chytridgenome[2:5], function(x) 
  renameDESeq(x, annotation = annotationChytrid))

chytrid_met_effect_2orgs <- renameDESeq(contrast_chytridgenome$resr_met_effect_2orgs, 
                                        annotation = annotationChytrid)
chytrid_met_effect_1org <- renameDESeq(contrast_chytridgenome$resr_met_effect_1org, 
                                       annotation = annotationChytrid)
chytrid_inf_effect_control <- renameDESeq(contrast_chytridgenome$resr_inf_effect_control, 
                                          annotation = annotationChytrid)
chytrid_inf_effect_met <- renameDESeq(contrast_chytridgenome$resr_inf_effect_met, 
                                      annotation = annotationChytrid)

contrast_chytridgenome <- list(chytrid_met_effect_2orgs=chytrid_met_effect_2orgs,
                               chytrid_met_effect_1org=chytrid_met_effect_1org,
                               chytrid_inf_effect_control=chytrid_inf_effect_control,
                        chytrid_inf_effect_met=chytrid_inf_effect_met)
                        


table(contrast_chytridgenome$chytrid_met_effect_2orgs$padj < 0.05)

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

makeVolcano(
  res = chytrid_met_effect_2orgs,
  title =  "MET effect on chytrid infecting bacteria gene expression")

## rm 

table(chytrid_met_effect_2orgs$padj < 0.05)

# Can you plot the shrunken LFC using lfcShrink.


# My guess is that the groups are too heterogenous to combine them and you should just run pairs, like so
# https://support.bioconductor.org/p/9150170/

### Genes of interest only present in one condition:
rownames(RSEM_final_hope.gene) <- RSEM_final_hope.gene$X

dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ condition)



> dds <- DESeq(dds)
> res<-results(dds,independentFiltering = F)



################################################################################
# (a) 360 chytrid genes in chytrid alone only
RSEM_final_hope.gene_a <- RSEM_final_hope.gene[
  RSEM_final_hope.gene$whichOrg %in% "in_chytrid_alone" &
    RSEM_final_hope.gene$whichTranscriptome %in% "chytrid",]
## Select only chytrid samples
RSEM_final_hope.gene_a<- RSEM_final_hope.gene_a[
  grep("chy", names(RSEM_final_hope.gene_a))]




# 
# 
# ## Keep genes with expression levels in at least 3 samples/trt (met/nomet)
# RSEM_final_hope.gene_d = RSEM_final_hope.gene_d[
#   (rowSums(RSEM_final_hope.gene_d[grep("control", names(RSEM_final_hope.gene_d))]>0)>
#      ncol(RSEM_final_hope.gene_d[grep("control", names(RSEM_final_hope.gene_d))])/2) &
#     (rowSums(RSEM_final_hope.gene_d[grep("met", names(RSEM_final_hope.gene_d))]>0)>
#        ncol(RSEM_final_hope.gene_d[grep("met", names(RSEM_final_hope.gene_d))])/2),]

################################################################################


### Genes for DESeq2 -> see next script R01 R02






# (d)
# 215 --> cyano genes only expressed outside of infection?
## TO DO
# in how many samples are they found (>half)
# top expressed
# --> GO on these to find the big functions
# DEG MET/not MET --> does met affects the virulence?
RSEM_final_hope.gene_d <- RSEM_final_hope.gene[
  RSEM_final_hope.gene$whichOrg %in% "in_cyano_alone" &
    RSEM_final_hope.gene$whichTranscriptome %in% "cyano",]
## Select only cyano samples
RSEM_final_hope.gene_d <- RSEM_final_hope.gene_d[
  grep("cyano", names(RSEM_final_hope.gene_d))]

## At least gene exp non null in half of each group met/not met
RSEM_final_hope.gene_d = RSEM_final_hope.gene_d[
  (rowSums(RSEM_final_hope.gene_d[grep("control", names(RSEM_final_hope.gene_d))]>0)>
     ncol(RSEM_final_hope.gene_d[grep("control", names(RSEM_final_hope.gene_d))])/2) &
    (rowSums(RSEM_final_hope.gene_d[grep("met", names(RSEM_final_hope.gene_d))]>0)>
       ncol(RSEM_final_hope.gene_d[grep("met", names(RSEM_final_hope.gene_d))])/2),]

nrow(RSEM_final_hope.gene_d) # 124 remaining!!
# GO on those TBC --> see R03 script: no significant GO term

# DE on met? --> R01 no DEG found

df = melt(RSEM_final_hope.gene_d)
df$trt = ifelse(grepl("met", df$variable), "met", "control")
t.test(df$value[df$trt %in% "met"], df$value[df$trt %in% "control"])
## no difference between metolachlor values

## describe top 10

# (e)
## one gene cyano only expressed during infection
RSEM_final_hope.gene_e <- RSEM_final_hope.gene[
  RSEM_final_hope.gene$whichOrg %in% "in_both_organisms" &
    RSEM_final_hope.gene$whichTranscriptome %in% "cyano",]
## Select only cyano samples
RSEM_final_hope.gene_e <- RSEM_final_hope.gene_e[
  grep("both", names(RSEM_final_hope.gene_e))]

## At least gene exp non null in half of each group met/not met
RSEM_final_hope.gene_e[
  (rowSums(RSEM_final_hope.gene_e[grep("control", names(RSEM_final_hope.gene_e))]>0)>
     ncol(RSEM_final_hope.gene_e[grep("control", names(RSEM_final_hope.gene_e))])/2) &
    (rowSums(RSEM_final_hope.gene_e[grep("met", names(RSEM_final_hope.gene_e))]>0)>
       ncol(RSEM_final_hope.gene_e[grep("met", names(RSEM_final_hope.gene_e))])/2),]
## no
