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
## Low quality genes filtering ##
###################################
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
## Both organisms:
explo.plot(mysaturation_chy, toplot = 1, samples = 12:19)

counts_cyano <- as.matrix(RSEM_final_hope.gene_cyano)
eset_cyano <- ExpressionSet(assayData = counts_cyano)
mysaturation_cyano = dat(eset_cyano, k = 0, ndepth = 7, type = "saturation")

## Both organisms:
explo.plot(mysaturation_cyano, toplot = 1, samples = 1:8) 
## Cyano alone
explo.plot(mysaturation_cyano, toplot = 1, samples = 9:20) 

## Due to the discrepancies between the groups, we remove genes with
## zeros in at least 3/group (met+inf)
RSEM_final_hope.gene_chytrid <- filterRSEMno3Nullpergp(RSEM_final_hope.gene_chytrid)
RSEM_final_hope.gene_cyano <- filterRSEMno3Nullpergp(RSEM_final_hope.gene_cyano)

## check saturation
counts_chy <- as.matrix(RSEM_final_hope.gene_chytrid)
eset_chy <- ExpressionSet(assayData = counts_chy)
mysaturation_chy = dat(eset_chy, k = 0, ndepth = 7, type = "saturation")

ncol(RSEM_final_hope.gene_chytrid)
## all
explo.plot(mysaturation_chy, toplot = 1, samples = 1:19)
## 86.9 to 100% genes detected

counts_cyano <- as.matrix(RSEM_final_hope.gene_cyano)
eset_cyano <- ExpressionSet(assayData = counts_cyano)
mysaturation_cyano = dat(eset_cyano, k = 0, ndepth = 7, type = "saturation")

names(RSEM_final_hope.gene_cyano)
explo.plot(mysaturation_cyano, toplot = 1, samples = 1:20) 
## Both organisms:
explo.plot(mysaturation_cyano, toplot = 1, samples = 1:9) 
## Cyano alone
explo.plot(mysaturation_cyano, toplot = 1, samples = 10:20) 
## meth_both_In11 58.9% everything else 92.5 to 100%

RSEM_final_hope.gene_cyano = RSEM_final_hope.gene_cyano[
  !names(RSEM_final_hope.gene_cyano) %in% c("met_both_In11")]

nrow(RSEM_final_hope.gene_chytrid) # 835
names(RSEM_final_hope.gene_chytrid)

nrow(RSEM_final_hope.gene_cyano) # 555
names(RSEM_final_hope.gene_cyano)

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
    title = "no MET, chytrids zoospores vs chytrids infecting",
    mylogo="logos/logo1.png", positionLogo = "left")

V_chytrid_inf_effect_met <- makeVolcano(
  res = contrast_chytridgenome$resr_inf_effect_met,
  title = "MET, chytrids zoospores vs chytrids infecting ",
  mylogo="logos/logo2.png", positionLogo = "left")

V_chytrid_met_effect_1org <- makeVolcano(
  res = contrast_chytridgenome$resr_met_effect_1org,
  title = "free-living zoospores, no MET vs MET",
  mylogo="logos/logo3.png", positionLogo = "left")

V_chytrid_met_effect_2orgs <- makeVolcano(
  res = contrast_chytridgenome$resr_met_effect_2orgs,
  title = "chytrids infecting, no MET vs MET",
  mylogo="logos/logo4.png", positionLogo = "left")

V_cyano_inf_effect_control <- makeVolcano(
  res = contrast_cyanogenome$resr_inf_effect_control,
  title = "no MET, uninfected vs infected cyanobacteria",
  mylogo="logos/logo5.png", positionLogo = "left")

V_cyano_inf_effect_met <- makeVolcano(
  res = contrast_cyanogenome$resr_inf_effect_met,
  title = "MET, uninfected vs infected cyanobacteria",
  mylogo="logos/logo6.png", positionLogo = "left")

V_cyano_met_effect_1org <- makeVolcano(
  res = contrast_cyanogenome$resr_met_effect_1org,
  title = "uninfected cyanobacteria, no MET vs MET",
  mylogo="logos/logo7.png", positionLogo = "right")

V_cyano_met_effect_2orgs <- makeVolcano(
  res = contrast_cyanogenome$resr_met_effect_2orgs,
  title = "infected cyanobacteria, no MET vs MET",
  mylogo="logos/logo8.png", positionLogo = "left")

## open bigger window
W=10; H=8

dev.new(width = W, height = H)
pdf("../../figures/Fig2_chytrid_volc.pdf", width = W, height = H)
cowplot::plot_grid(V_chytrid_inf_effect_control$plot,
                   V_chytrid_inf_effect_met$plot,
                   V_chytrid_met_effect_1org$plot,
                   V_chytrid_met_effect_2orgs$plot,
                   labels = c("a", "b", "c", "d"), label_size = 20)
dev.off()

dev.new(width = W, height = H)
pdf("../../figures/Fig3_cyano_volc.pdf", width = W, height = H)
cowplot::plot_grid(V_cyano_inf_effect_control$plot,
                   V_cyano_inf_effect_met$plot,
                   V_cyano_met_effect_1org$plot,
                   V_cyano_met_effect_2orgs$plot, 
                   labels = c("a", "b", "c", "d"), label_size = 20)
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

ggvenn(
  venn_data, show_percentage = F, fill_color = rep("white", 4),
  stroke_size = 0.5, set_name_size = 4
)

# Combine the lists into a named list for the Venn diagram
venn_data <- list("infection effect cyano - control"= 
                    getGenes(contrast_cyanogenome$resr_inf_effect_control),
                  "met effect cyano - alone"= 
                    getGenes(contrast_cyanogenome$resr_met_effect_1org),
                  "met effect cyano - infected"= 
                    getGenes(contrast_cyanogenome$resr_met_effect_2orgs),
                  "infection effect cyano - met"= 
                    getGenes(contrast_cyanogenome$resr_inf_effect_met))

ggvenn(
  venn_data, show_percentage = F, fill_color = rep("white", 4),
  stroke_size = 0.5, set_name_size = 4
)

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

#############################################
## GO of the four groups per transcriptome ##
#############################################

universe_chytrid = rownames(RSEM_final_hope.gene_chytrid)
universe_cyano = rownames(RSEM_final_hope.gene_cyano)

getGOBubbleZ(universe = universe_chytrid, annotation = annotationChytrid, 
             genelist = getGenes(contrast_chytridgenome$resr_inf_effect_control), 
             GO_df = GO_chytrid, isbubble = F)
## axoneme (cellular component) padj=0.003898744

getGOBubbleZ(universe = universe_chytrid, annotation = annotationChytrid, 
             genelist = getGenes(contrast_chytridgenome$resr_met_effect_1org), 
             GO_df = GO_chytrid, isbubble = F)
# no significant GO terms
# getGOBubbleZ(universe = universe_chytrid, annotation = annotationChytrid, 
#              genelist = getGenes(contrast_chytridgenome$resr_met_effect_1org[
#                contrast_chytridgenome$resr_met_effect_1org$log2FoldChange < 0,]), 
#              GO_df = GO_chytrid, isbubble = F)

getGOBubbleZ(universe = universe_cyano, annotation = annotationCyano, 
             genelist = getGenes(contrast_cyanogenome$resr_met_effect_1org[
               contrast_cyanogenome$resr_met_effect_1org$log2FoldChange > 0,]), 
             GO_df = GO_cyano, isbubble = F)

getGOBubbleZ(universe = universe_chytrid, annotation = annotationChytrid, 
             genelist = getGenes(contrast_chytridgenome$resr_met_effect_2orgs), 
             GO_df = GO_chytrid, isbubble = F)
# axoneme p.adj 0.007747258

# getGOBubbleZ(universe = universe_chytrid, annotation = annotationChytrid, 
#              genelist = getGenes(contrast_chytridgenome$resr_inf_effect_met), 
#              GO_df = GO_chytrid, isbubble = F)

## 2. Cyano
# getGOBubbleZ(universe = universe_cyano, annotation = annotationCyano,
#              genelist = getGenes(contrast_cyanogenome$resr_inf_effect_control),
#              GO_df = GO_cyano, isbubble = F)
# no DEG

getGOBubbleZ(universe = universe_cyano, annotation = annotationCyano, 
              genelist = getGenes(contrast_cyanogenome$resr_met_effect_1org), 
              GO_df = GO_cyano, isbubble = F)
# "no significant GO terms"

# getGOBubbleZ(universe = universe_cyano, annotation = annotationCyano,
#              genelist = getGenes(contrast_cyanogenome$resr_met_effect_2orgs),
#              GO_df = GO_cyano, isbubble = F)
# no DEG

## infection effect on cyanobacteria in presence of metolachlor
# getGOBubbleZ(universe = universe_cyano, annotation = annotationCyano, 
#              genelist = getGenes(contrast_cyanogenome$resr_inf_effect_met), 
#              GO_df = GO_cyano, isbubble = F)
# no signif
