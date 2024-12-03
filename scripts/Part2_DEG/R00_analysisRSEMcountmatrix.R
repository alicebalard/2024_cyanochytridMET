## 29th of November 2024
source("dataLoad.R")

## Let's observe the count matrix calculated by Trinity
RSEM_final_hope.gene <- 
  read.csv("../../data/run_DESEQ2_Erika/RSEM_final_hope.gene.counts.matrix", sep="\t")

##########################################################################
## Split by type of genes depending on in which samples they are expressed
# (1) chytrid genes expressed in chytrid alone and infecting --> give to Erika for DESEq2 1623 genes
# (2) cyanobacteri genes expressed in cyano alone and infected --> give to Erika for DESEq2 3420 genes
# (3) a mix that we need to understand --> GO analysis
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
#              in_both_organisms in_both_organisms in_cyano_alone in_chytrid_alone
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

# (1) chytrid genes expressed in chytrid alone and infecting --> give to Erika for DESEq2 1623 genes

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

# (2) cyanobacteria genes expressed in cyano alone and infected --> give to Erika for DESEq2 3420 genes
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

##########################################
## What can we eliminate as contamination? TBC later
# # all the "chytrid" genes found when cyanobacteria is alone are contaminants! N=912
# contaChytridTrans <- RSEM_final_hope.gene[RSEM_final_hope.gene$whichTranscriptome %in% "chytrid" & 
#                               grepl("in_cyano_alone", RSEM_final_hope.gene$whichOrg),"X"]
# 
# # all the "cyano" genes found when chytrid is alone are contaminants! N=163
# contaCyanoTrans <- RSEM_final_hope.gene[RSEM_final_hope.gene$whichTranscriptome %in% "cyano" & 
#                                                    grepl("in_chytrid_alone", RSEM_final_hope.gene$whichOrg),"X"]
# 
# ## Explore contaminants
# gene_trans_map_cyano[gene_trans_map_cyano$V1 %in% contaCyanoTrans,]
# 
# table(annotationChytrid$blastxPhylum[
#   annotationChytrid$X.gene_id %in% contaChytridTrans])
# 
# ## Remove contaminants
# RSEM_final_hope.gene <- RSEM_final_hope.gene[!RSEM_final_hope.gene$X %in% c(contaChytridTrans, contaCyanoTrans),]
# 
# table(RSEM_final_hope.gene$whichTranscriptome,
#       RSEM_final_hope.gene$whichOrg) # much better
# 
# #         in_both_organisms in_both_organisms in_cyano_alone in_chytrid_alone
# # chytrid              1764                                0              360
# # cyano                   1                             3420                0
# # 
# #         in_chytrid_alone in_both_organisms in_cyano_alone
# # chytrid                               1623              0
# # cyano                                    0            215
# 
# 
# 
# ## Same question, but for genes found in ALL SAMPLES
# df = RSEM_final_hope.gene_chytrid[grep("cyano", names(RSEM_final_hope.gene_chytrid), 
#                                        invert = T)]
# chytridGenesCompleteCase <-
#   RSEM_final_hope.gene_chytrid[row.names(df[apply(df[-1], 1, function(x) all(x > 0)),]),]
# 
# df = RSEM_final_hope.gene_cyano[grep("chytrid", names(RSEM_final_hope.gene_cyano), 
#                                      invert = T)]
# cyanoGenesCompleteCase <-
#   RSEM_final_hope.gene_cyano[row.names(df[apply(df[-1], 1, function(x) all(x > 0)),]),]
# 
# table(apply(df[-1], 1, function(x) all(x > 0)))
# 
# 
# hist(RSEM_final_hope.gene_chytrid$control_chy_Z1, breaks = 100)
# hist(RSEM_final_hope.gene_chytrid$control_chy_Z2, breaks = 100)
# hist(RSEM_final_hope.gene_chytrid$control_chy_Z3, breaks = 100)
# hist(RSEM_final_hope.gene_chytrid$control_chy_Z4, breaks = 100)
# hist(RSEM_final_hope.gene_chytrid$control_chy_Z5, breaks = 100)
# hist(RSEM_final_hope.gene_chytrid$control_chy_Z6, breaks = 100)
# hist(RSEM_final_hope.gene_chytrid$met_chy_Z7, breaks = 100)
# hist(RSEM_final_hope.gene_chytrid$met_chy_Z8, breaks = 100)
# hist(RSEM_final_hope.gene_chytrid$met_chy_Z10, breaks = 100)
# hist(RSEM_final_hope.gene_chytrid$met_chy_Z11, breaks = 100)
# hist(RSEM_final_hope.gene_chytrid$met_chy_Z12, breaks = 100)
# 
# hist(RSEM_final_hope.gene_cyano$control_both_In5, breaks = 100)
# 
# hist(RSEM_final_hope.gene_cyano$met_cyano_C8, breaks = 100)
# 
# range(RSEM_final_hope.gene_cyano$met_cyano_C8, na.rm = T)
# 
# # (3) a mix that we need to understand --> GO analysis
# 
# #### TBC
# 
# # 1764 chytrid genes are expressed only when both organisms are together
# DF <- RSEM_final_hope.gene[
#   RSEM_final_hope.gene$whichTranscriptome %in% "chytrid" &
#   RSEM_final_hope.gene$whichOrg %in% "in_both_organisms",]
# 
# 
# # annotationChytrid$gene_ontology_BLASTX[
#   # annotationChytrid$X.gene_id %in% DF$X]
# 
# 
# # in_both_organisms, chytrid transcriptome, in a lot of samples
# table(rowSums(DF[grep("both", names(DF))] !=0))
# 
# DF[rowSums(DF[grep("both", names(DF))] !=0) == 5,]
# 
# "TRINITY_DN2472_c0_g1" # virulence factor in a plant fungi
# "TRINITY_DN34253_c0_g1" # in arabidopsis
# "TRINITY_DN611_c0_g1" # RNA metabolic process
# 
# assemblyMergedFungi_filterEuk_simplified[
#   assemblyMergedFungi_filterEuk_simplified$X.gene_id %in% "TRINITY_DN611_c0_g1",]
# 
# DF[rowSums(DF[grep("both", names(DF))] !=0) == 4,"X"]
# 
# assemblyMergedFungi_filterEuk_simplified[
#   assemblyMergedFungi_filterEuk_simplified$X.gene_id %in% "TRINITY_DN4458_c0_g1",]
# 
# DF[rowSums(DF[grep("both", names(DF))] !=0) == 4,][
#   grep("both",DF[rowSums(DF[grep("both", names(DF))] !=0) == 4,])]
# 
# ## Highest mean of expresssion, only in "both", reads in at least 4 samples
# DF[rowSums(DF[grep("both", names(DF))] !=0) == 4,][
#   order(rowMeans(DF[rowSums(DF[grep("both", names(DF))] !=0) == 4,][
#   grep("both", names(DF[rowSums(DF[grep("both", names(DF))] !=0) == 4,]))]), decreasing = T),]
# 
# assemblyMergedFungi_filterEuk_simplified[
#   assemblyMergedFungi_filterEuk_simplified$X.gene_id %in% "TRINITY_DN30306_c0_g1",]
# 
# 
# 
# library(dplyr)
# 
# assemblyMergedFungi_filterEuk_simplified[
#   assemblyMergedFungi_filterEuk_simplified$X.gene_id %in% RSEM_final_hope.gene[
#     RSEM_final_hope.gene$whichTranscriptome %in% "chytrid" &
#       RSEM_final_hope.gene$whichOrg %in% "in_both_organisms","X"] & 
#     assemblyMergedFungi_filterEuk_simplified$blastxPhylum %in% "Fungi","sprot_Top_BLASTX_hit"] %>% unique() 
# 
# 
# RSEM_final_hope.gene[
#   RSEM_final_hope.gene$X %in% "TRINITY_DN891_c0_g1",]
# 
# RSEM_final_hope.gene[1:10,]



