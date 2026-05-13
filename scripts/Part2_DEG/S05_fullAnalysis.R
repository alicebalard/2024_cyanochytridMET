## Updated 2026
library(here)
source("libLoad.R")
source("dataLoad.R")
# === Chytrid annotation summary ===
#   Total genes:               4860
# With GO terms:             4842
# Without GO terms:          46
# Unique gene-GO pairs:      142801
# Unique GO terms:           9419
# === Cyano annotation summary ===
#   Total genes:               4451
# With protein_id:           4356
# With GO terms:             1842
# Unmatched (no protein_id): 95
# Unique GO terms:           1314
# Unique gene-GO pairs: 5876
# Unique GO terms:      1314
source("functions.R")

## Load the count matrix calculated by Trinity
RSEM_newhope.gene <- read.csv(here("data/RSEM_new_hope.gene.counts.rmConta.matrix"), sep="\t")
RSEM_newhope.gene <- column_to_rownames(RSEM_newhope.gene, "X")
nrow(RSEM_newhope.gene)
## 8345 genes (after removing conta)

#######################################################################
## Split by group depending on which gene is expressed in which case ##
#######################################################################
a = ifelse(rowSums(RSEM_newhope.gene[grep("chy", names(RSEM_newhope.gene))]) !=0,
           "in_chytrid_alone", "")
b = ifelse(rowSums(RSEM_newhope.gene[grep("both", names(RSEM_newhope.gene))]) !=0,
           "in_both_organisms", "")
c = ifelse(rowSums(RSEM_newhope.gene[grep("cyano", names(RSEM_newhope.gene))]) !=0,
           "in_cyano_alone", "")

RSEM_newhope.gene$whichOrg <- trimws(paste(a,b,c, sep = " "))

RSEM_newhope.gene$whichTranscriptome <- ifelse(
  grepl("TRINITY", rownames(RSEM_newhope.gene)), "chytrid", "cyano")

table(RSEM_newhope.gene$whichOrg,
      RSEM_newhope.gene$whichTranscriptome)
#                                                    chytrid cyano
#                                                       286   603
# in_both_organisms                                    1702     2
# in_both_organisms in_cyano_alone                        0  3401
# in_chytrid_alone                                      330     0
# in_chytrid_alone  in_cyano_alone                        0     2
# in_chytrid_alone in_both_organisms                   1576     0
# in_chytrid_alone in_both_organisms in_cyano_alone       0   186
# in_cyano_alone                                          0   257

## All good, contamination removed

#### CHYTRID ####

RSEM_chytrid <- RSEM_newhope.gene[RSEM_newhope.gene$whichTranscriptome %in% "chytrid",] %>%
  dplyr::select(!c(whichOrg, whichTranscriptome))

## Rename based on annotations
table(is.na(annotationChytrid_final$gene_name))
table(is.na(annotationChytrid_final$gene_id))
## good, all annotated

row.names(RSEM_chytrid) <- make.unique(annotationChytrid_final$gene_name[
  match(row.names(RSEM_chytrid), annotationChytrid_final$gene_id)])

## Merge identical proteins in only one row, suming the counts
RSEM_chytrid = RSEM_chytrid %>%
  mutate(base_name = sub("\\.\\d+$", "", rownames(RSEM_chytrid))) %>%
  group_by(base_name) %>%
  summarise(across(everything(), sum)) %>%
  tibble::column_to_rownames("base_name") %>% data.frame()

## remove cyano only samples
RSEM_chytrid <- RSEM_chytrid[!grepl("cyano", names(RSEM_chytrid))]

nrow(RSEM_chytrid) # 3240 genes

#### CYANO ####
RSEM_cyano <- RSEM_newhope.gene[RSEM_newhope.gene$whichTranscriptome %in% "cyano",]%>%
  dplyr::select(!c(whichOrg, whichTranscriptome))

## Rename based on annotations
table(is.na(annotationCyano_final$gene_id))
table(is.na(annotationCyano_final$custom_gene_name))

rownames(RSEM_cyano) = make.unique(annotationCyano_final$gene_id[
  match(row.names(RSEM_cyano), annotationCyano_final$custom_gene_name)])

## remove chytrid only samples
RSEM_cyano <- RSEM_cyano[!grepl("chy", names(RSEM_cyano))]

nrow(RSEM_cyano) # 4451 genes

## Merge identical proteins in only one row, suming the counts
RSEM_cyano = RSEM_cyano %>%
  mutate(base_name = sub("\\.\\d+$", "", rownames(RSEM_cyano))) %>%
  group_by(base_name) %>%
  summarise(across(everything(), sum)) %>%
  tibble::column_to_rownames("base_name") %>% data.frame()

nrow(RSEM_cyano) # 4451 genes

#############################################
## Low quality genes & samples filtering   ##
#############################################

## Chytrid
counts_chy <- as.matrix(RSEM_chytrid)
counts_chy <- counts_chy[rowSums(counts_chy) > 0, ]  # remove all-zero rows

## Cyano
counts_cyano <- as.matrix(RSEM_cyano)
counts_cyano <- counts_cyano[rowSums(counts_cyano) > 0, ]  # remove all-zero rows

# Saturation plot: number of genes detected at observed and simulated
# sequencing depths, and marginal gain in gene detection per million additional reads
eset_chy <- ExpressionSet(assayData = counts_chy)
mysaturation_chy <- dat(eset_chy, k = 0, ndepth = 7, type = "saturation")

eset_cyano <- ExpressionSet(assayData = counts_cyano)
mysaturation_cyano <- dat(eset_cyano, k = 0, ndepth = 7, type = "saturation")

## Plots:
pdf(here("figures/saturation_plots.pdf"), width=13, height=9)

## Chytrid
par(mfrow=c(2,2))  # 2x2 layout
explo.plot(mysaturation_chy,   toplot=1, samples=grep("chy",  colnames(counts_chy)))
title(sub="Chytrid alone")
explo.plot(mysaturation_chy,   toplot=1, samples=grep("both", colnames(counts_chy)))
title(sub="Chytrid in co-culture")

## Cyano
explo.plot(mysaturation_cyano, toplot=1, samples=grep("cyano", colnames(counts_cyano)))
title(sub="Cyano alone")
explo.plot(mysaturation_cyano, toplot=1, samples=grep("both",  colnames(counts_cyano)))
title(sub="Cyano in co-culture")

dev.off()
par(mfrow=c(1,1))  # reset layout

# Ratio chytrid/cyano counts per sample
both_chy   <- colSums(counts_chy)[grep("both", colnames(counts_chy))]
both_cyano <- colSums(counts_cyano)[grep("both", colnames(counts_cyano))]
# Make sure same samples in both
common_both <- intersect(names(both_chy), names(both_cyano))
sort(both_chy[common_both] / both_cyano[common_both])

# control_both_In3    met_both_In12 control_both_In5     met_both_In9 
# 0.9724276        1.0419057        1.2084985        3.7342325 
# met_both_In8     met_both_In7 control_both_In2 control_both_In1 
# 4.2212903       20.4887370       81.2082762      146.3645498 
# met_both_In10 control_both_In6 control_both_In4    met_both_In11 
# 180.9254192      184.6268227      366.2886701      446.6517864 

## For samples alone, Z9 is an outlier for chytrid, we can keep all samples for cyano

## For co-infection, certain samples have an unbalanced cyano/chytrid read depth
# Ratio ~1:     In12, In3, In5  → balanced cyano/chytrid ✅
# Ratio <100:   In2, In7, In9, In8   → moderate chytrid dominance
# Ratio >100:    In1, In4, In6, In10, In11   → high chytrid dominance ⚠️

## Based on (1) saturation plots and (2) ratio chytrid/cyano counts per sample,
## we remove the following samples:

# Chytrid matrix: met_chy_Z9 and control_both_In1 are clearly outliers (14.3%, 69.3%)
remove_chy <- c("met_chy_Z9", # 14.3% saturation vs 30%+ for all others
                "control_both_In1",
                "met_both_In11") # big outlier

# Cyano matrix: 2 different for cyano alone (keep all) vs co-culture (remove outliers)
remove_cyano <- c("control_both_In1",   # 0.9% saturation + high chytrid dominance
                  "control_both_In4",   # 1.3% saturation + high chytrid dominance
                  "control_both_In6",   # 1.9% saturation + high chytrid dominance
                  "met_both_In10", # 0.9% saturation + high chytrid dominance
                  "met_both_In11") # 36.8% saturation + high chytrid dominance

## --- 1. Infection effect in chytrid, no MET ---
## --- 2. Infection effect in chytrid, with MET ---
## --- 3. MET effect in chytrid alone ---
## --- 4. MET effect in chytrid infecting cyanobacteria
## --- 5. Infection effect in cyano, no MET ---
## --- 6. Infection effect in cyano, with MET ---
## --- 7. MET effect in cyano alone ---
## --- 8. MET effect in cyano infected by chytrid ---

# Group-aware filtering: require >= 10 counts in >= 3 samples
# within all groups. This prevents genes expressed
# in only 1-2 samples from passing the filter due to cross-group pooling,
# which would otherwise cause extreme LFC inflation.
filter_by_group <- function(counts, min_count=10, min_samples=3) {
  groups <- sub("_[^_]+$", "", colnames(counts))
  keep <- sapply(rownames(counts), function(gene) {
    all(sapply(unique(groups), function(grp) {
      grp_counts <- counts[gene, groups == grp]
      sum(grp_counts >= min_count) >= min_samples
    }))
  })
  counts[keep, ]
}

## Subset count matrices
counts1_infeffect_noMET_chy <- counts_chy[ ,grep("control", colnames(counts_chy))]
counts1_infeffect_noMET_chy <- counts1_infeffect_noMET_chy[
  , !colnames(counts1_infeffect_noMET_chy) %in% remove_chy]
colnames(counts1_infeffect_noMET_chy) # 5 vs 6

counts2_infeffect_MET_chy <- counts_chy[ ,grep("met", colnames(counts_chy))]
counts2_infeffect_MET_chy <- counts2_infeffect_MET_chy[
  , !colnames(counts2_infeffect_MET_chy) %in% remove_chy]
colnames(counts2_infeffect_MET_chy) # 5 vs 6

counts3_METeffect_alone_chy <- counts_chy[ ,grep("chy", colnames(counts_chy))]
counts3_METeffect_alone_chy <- counts3_METeffect_alone_chy[
  , !colnames(counts3_METeffect_alone_chy) %in% remove_chy]
colnames(counts3_METeffect_alone_chy) # 6 vs 5

counts4_METeffect_infecting_chy <- counts_chy[ ,grep("both", colnames(counts_chy))]
counts4_METeffect_infecting_chy <- counts4_METeffect_infecting_chy[
  , !colnames(counts4_METeffect_infecting_chy) %in% remove_chy]
colnames(counts4_METeffect_infecting_chy) # 5 vs 6

counts5_infeffect_noMET_cyano <- counts_cyano[ ,grep("control", colnames(counts_cyano))]
counts5_infeffect_noMET_cyano <- counts5_infeffect_noMET_cyano[
  , !colnames(counts5_infeffect_noMET_cyano) %in% remove_cyano]
colnames(counts5_infeffect_noMET_cyano) # 3 vs 6 LESS RELIABLE!!

counts6_infeffect_MET_cyano <- counts_cyano[ ,grep("met", colnames(counts_cyano))]
counts6_infeffect_MET_cyano <- counts6_infeffect_MET_cyano[
  , !colnames(counts6_infeffect_MET_cyano) %in% remove_cyano]
colnames(counts6_infeffect_MET_cyano) # 4 vs 6 LESS RELIABLE!!

counts7_METeffect_alone_cyano <- counts_cyano[ ,grep("cyano", colnames(counts_cyano))]
counts7_METeffect_alone_cyano <- counts7_METeffect_alone_cyano[
  , !colnames(counts7_METeffect_alone_cyano) %in% remove_cyano]
colnames(counts7_METeffect_alone_cyano) # 6 vs 6

counts8_METeffect_infecting_cyano <- counts_cyano[ ,grep("both", colnames(counts_cyano))]
counts8_METeffect_infecting_cyano <- counts8_METeffect_infecting_cyano[
  , !colnames(counts8_METeffect_infecting_cyano) %in% remove_cyano]
colnames(counts8_METeffect_infecting_cyano) # 3 vs 4

## Reapply filtering
# Create named list automatically from object names
listCounts <- mget(c("counts1_infeffect_noMET_chy",
                     "counts2_infeffect_MET_chy",
                     "counts3_METeffect_alone_chy",
                     "counts4_METeffect_infecting_chy",
                     "counts5_infeffect_noMET_cyano",
                     "counts6_infeffect_MET_cyano",
                     "counts7_METeffect_alone_cyano",
                     "counts8_METeffect_infecting_cyano"))

listCounts <- lapply(listCounts, filter_by_group)

# Automated summary for all count matrices in the list
invisible(lapply(names(listCounts), function(name) {
  x <- listCounts[[name]]
  message(name, ":  ", ncol(x), " samples, ", nrow(x), " genes")
}))
## invisible() suppresses the NULL list output that lapply would otherwise print to console.

# counts1_infeffect_noMET_chy:  11 samples, 837 genes
# counts2_infeffect_MET_chy:  11 samples, 1030 genes
# counts3_METeffect_alone_chy:  11 samples, 941 genes
# counts4_METeffect_infecting_chy:  11 samples, 856 genes
# counts5_infeffect_noMET_cyano:  9 samples, 95 genes ***** low trust
# counts6_infeffect_MET_cyano:  10 samples, 988 genes ***** low trust
# counts7_METeffect_alone_cyano:  12 samples, 3553 genes
# counts8_METeffect_infecting_cyano:  7 samples, 66 genes

#####################
## DESeq2          ##
#####################

## rebuild to match filtered samples
samples_data <- data.frame(
  sample    = union(colnames(counts_chy), colnames(counts_cyano)),
  condition = sub("_[^_]+$", "", union(colnames(counts_chy), colnames(counts_cyano)))
)
rownames(samples_data) <- samples_data$sample 
samples_data$condition <- as.factor(samples_data$condition)

###############
## Functions ##
###############

testEffectSize <- function(count){
  ddsr <- DESeqDataSetFromMatrix(
    countData = round(count),
    colData   = samples_data[colnames(count), ],
    design    = ~ condition)
  ddsr_check <- estimateSizeFactors(ddsr)  # faster than full DESeq()
  message("Size factors (should be ~1):")
  print(sort(sizeFactors(ddsr_check)))
}

myDESeq2_withlfcShrink <- function(count, a, b){
  ddsr <- DESeqDataSetFromMatrix(
    countData = round(count),
    colData   = samples_data[colnames(count), ],
    design    = ~ condition)
  
  ddsr <- DESeq(ddsr, fitType="local")
  res <- results(ddsr,
                 contrast = c("condition", b, a),
                 alpha    = 0.05)
  res <- lfcShrink(ddsr,
                   contrast = c("condition", b, a),
                   res      = res,
                   type     = "ashr")
  print(summary(res))
  return(res)
}

testIfCrazyFoldChange <- function(res){
  res %>% as.data.frame() %>%
    filter(!is.na(padj) & padj < 0.05) %>%
    arrange(desc(abs(log2FoldChange))) %>%
    head(10) %>%
    dplyr::select(baseMean, log2FoldChange, padj)
}
############################

lapply(listCounts, testEffectSize) ## In11 flagged as high chytrid

# Run all contrasts using the filtered listCounts
## using 'ashr' for LFC shrinkage. If used in published research, please cite:
# Stephens, M. (2016) False discovery rates: a new deal. Biostatistics, 18:2.
# https://doi.org/10.1093/biostatistics/kxw041

## --- 1. Infection effect in chytrid, no MET ---

res_counts1_infeffect_noMET_chy <- myDESeq2_withlfcShrink(
  listCounts$counts1_infeffect_noMET_chy, a="control_chy", b="control_both")
# out of 837 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 3, 0.36%
# LFC < 0 (down)     : 3, 0.36%
# outliers [1]       : 24, 2.9%
# low counts [2]     : 0, 0%
# (mean count < 10)

testIfCrazyFoldChange(res_counts1_infeffect_noMET_chy) # OK

## --- 2. Infection effect in chytrid, with MET ---
res_counts2_infeffect_MET_chy <- myDESeq2_withlfcShrink(
  listCounts$counts2_infeffect_MET_chy, a="met_chy", b="met_both")
# out of 1030 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 7, 0.68%
# LFC < 0 (down)     : 3, 0.29%
# outliers [1]       : 27, 2.6%
# low counts [2]     : 60, 5.8%
# (mean count < 27)

testIfCrazyFoldChange(res_counts2_infeffect_MET_chy) # OK

## --- 3. MET effect in chytrid alone ---
res_counts3_METeffect_alone_chy <- myDESeq2_withlfcShrink(
  listCounts$counts3_METeffect_alone_chy, a="control_chy", b="met_chy")

# out of 941 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 2, 0.21%
# LFC < 0 (down)     : 4, 0.43%
# outliers [1]       : 19, 2%
# low counts [2]     : 0, 0%
# (mean count < 10)

testIfCrazyFoldChange(res_counts3_METeffect_alone_chy) # ok

## --- 4. MET effect in chytrid infecting cyanobacteria
res_counts4_METeffect_infecting_chy <- myDESeq2_withlfcShrink(
  listCounts$counts4_METeffect_infecting_chy, a="control_both", b="met_both")
# out of 856 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 6, 0.7%
# LFC < 0 (down)     : 7, 0.82%
# outliers [1]       : 31, 3.6%
# low counts [2]     : 132, 15%
# (mean count < 62)

testIfCrazyFoldChange(res_counts4_METeffect_infecting_chy)

## --- 5. Infection effect in cyano, no MET --- ***** UNRELIABLE
res_counts5_infeffect_noMET_cyano <- myDESeq2_withlfcShrink(
  listCounts$counts5_infeffect_noMET_cyano, a="control_cyano", b="control_both")
# out of 95 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 15, 16%
# LFC < 0 (down)     : 25, 26%
# outliers [1]       : 0, 0%
# low counts [2]     : 0, 0%
# (mean count < 136)

testIfCrazyFoldChange(res_counts5_infeffect_noMET_cyano)

## --- 6. Infection effect in cyano, with MET --- ***** UNRELIABLE
res_counts6_infeffect_MET_cyano <- myDESeq2_withlfcShrink(
  listCounts$counts6_infeffect_MET_cyano, a="met_cyano", b="met_both")
# out of 988 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 39, 3.9%
# LFC < 0 (down)     : 117, 12%
# outliers [1]       : 1, 0.1%
# low counts [2]     : 0, 0%
# (mean count < 47)
testIfCrazyFoldChange(res_counts6_infeffect_MET_cyano)

## --- 7. MET effect in cyano alone ---
res_counts7_METeffect_alone_cyano <- myDESeq2_withlfcShrink(
  listCounts$counts7_METeffect_alone_cyano, a="control_cyano", b="met_cyano")
# out of 3553 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 1, 0.028%
# LFC < 0 (down)     : 0, 0%
# outliers [1]       : 0, 0%
# low counts [2]     : 0, 0%
# (mean count < 9)

testIfCrazyFoldChange(res_counts7_METeffect_alone_cyano)

# An environmentally relevant concentration of metolachlor induced minimal transcriptional 
# response in P. agardhii, with only one gene (GeneID:77286457) significantly differentially
# expressed (padj=0.04, LFC=0.29). This weak response is consistent with the sub-lethal, 
# environmentally relevant exposure concentration used, and suggests that at this concentration, 
# metolachlor does not strongly perturb cyanobacterial gene expression under single-species conditions."

# What is that gene?
annotationCyano_final[annotationCyano_final$gene_id == "77286457", 
                      c("custom_gene_name", "locus_tag", "gene", "protein", "product")]
# 1 DEG: iron uptake porin (GeneID:77286457, padj=0.04, LFC=0.28)
# Minimal transcriptional response consistent with sub-lethal environmentally relevant dose

## p-value distribution check (uniform = well-calibrated)
hist(res_counts7_METeffect_alone_cyano$pvalue, breaks=50,
     main="p-value distribution: control_cyano vs met_cyano")

## --- 8. MET effect in cyano infected by chytrid ---
res_counts8_METeffect_infecting_cyano <- myDESeq2_withlfcShrink(
  listCounts$counts8_METeffect_infecting_cyano, a="control_both", b="met_both")
# out of 66 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 0, 0%
# LFC < 0 (down)     : 0, 0%
# outliers [1]       : 5, 7.6%
# low counts [2]     : 0, 0%
# (mean count < 20)

## --- Combine into named lists ---
contrast_chytridgenome <- mget(c(
  "res_counts1_infeffect_noMET_chy",  
  "res_counts2_infeffect_MET_chy",
  "res_counts3_METeffect_alone_chy",
  "res_counts4_METeffect_infecting_chy"
))

contrast_cyanogenome <- mget(c(
  "res_counts5_infeffect_noMET_cyano",  # ⚠️ interpret cautiously
  "res_counts6_infeffect_MET_cyano",  # ⚠️ interpret cautiously
  "res_counts7_METeffect_alone_cyano",  # MET effect, cyano alone
  "res_counts8_METeffect_infecting_cyano"   # MET effect, co-culture
))

# Automated summary for all count matrices in the list
invisible(lapply(names(contrast_chytridgenome), function(name) {
  x <- listCounts[[name]]
}))
invisible(lapply(names(contrast_cyanogenome), function(name) {
  x <- listCounts[[name]]
}))

########## PLOTS ########## 

## Volcano plots
V_chytrid_inf_effect_control <- makeVolcano(
  res          = contrast_chytridgenome$res_counts1_infeffect_noMET_chy,
  title        = "no MET, chytrids zoospores vs chytrids infecting",
  positionLogoStart = 0, positionLogoStop = 0,
  
  mylogo       = "logos/logo1.png")

V_chytrid_inf_effect_met <- makeVolcano(
  res          = contrast_chytridgenome$res_counts2_infeffect_MET_chy,
  title        = "MET, chytrids zoospores vs chytrids infecting",
  positionLogoStart = 0, positionLogoStop = 0,
  
  mylogo       = "logos/logo2.png")

V_chytrid_met_effect_1org <- makeVolcano(
  res          = contrast_chytridgenome$res_counts3_METeffect_alone_chy,
  title        = "free-living zoospores, no MET vs MET",
  positionLogoStart = 0, positionLogoStop = 0,
  
  mylogo       = "logos/logo3.png")

V_chytrid_met_effect_2orgs <- makeVolcano(
  res          = contrast_chytridgenome$res_counts4_METeffect_infecting_chy,
  title        = "chytrids infecting, no MET vs MET",
  positionLogoStart = 0, positionLogoStop = 0,
  mylogo       = "logos/logo4.png")

# V_cyano_inf_effect_control <- makeVolcano(
  # res          = contrast_cyanogenome$res_counts5_infeffect_noMET_cyano,
#   title        = "no MET, uninfected vs infected cyanobacteria [interpret cautiously]", #***
#   mylogo       = "logos/logo5.png", positionLogoStart = -1, positionLogoStop = 3)
# 
# V_cyano_inf_effect_met <- makeVolcano(
  # res          = contrast_cyanogenome$res_counts6_infeffect_MET_cyano,
#   title        = "MET, uninfected vs infected cyanobacteria [interpret cautiously]", #***
#   mylogo       = "logos/logo6.png", positionLogoStart = -1, positionLogoStop = 3)

V_cyano_met_effect_1org <- makeVolcano(
  res          = contrast_cyanogenome$res_counts7_METeffect_alone_cyano,
  title        = "uninfected cyanobacteria, no MET vs MET",
  mylogo       = "logos/logo7.png")

V_cyano_met_effect_2orgs <- makeVolcano(
  res          = contrast_cyanogenome$res_counts8_METeffect_infecting_cyano,
  title        = "infected cyanobacteria, no MET vs MET",
  mylogo       = "logos/logo8.png")

## Save figures
W <- 10; H <- 8

pdf(here("figures/Fig2_chytrid_volc.pdf"), width=W, height=H)
cowplot::plot_grid(V_chytrid_inf_effect_control$plot,
                   V_chytrid_inf_effect_met$plot,
                   V_chytrid_met_effect_1org$plot,
                   V_chytrid_met_effect_2orgs$plot,
                   labels=c("a","b","c","d"), label_size=20)
dev.off()

pdf(here("figures/Fig3_cyano_volc.pdf"), width=W, height=H/2)
cowplot::plot_grid(V_cyano_met_effect_1org$plot,
                   V_cyano_met_effect_2orgs$plot,
                   labels=c("a","b"), label_size=20)
dev.off()

## Venn diagrams
getGenes <- function(x) rownames(x[!is.na(x$padj) & x$padj < 0.05, ])

p1 <- ggVennDiagram(x = list("Infection effect\nabsence of metolachlor" = 
                               getGenes(contrast_chytridgenome$res_counts1_infeffect_noMET_chy),
                             "Metolachlor effect\nfree living zoospores" = 
                               getGenes(contrast_chytridgenome$res_counts2_infeffect_MET_chy),
                             "Metolachlor effect\nduring infection" = 
                               getGenes(contrast_chytridgenome$res_counts3_METeffect_alone_chy),
                             "Infection effect\npresence of metolachlor"     = 
                               getGenes(contrast_chytridgenome$res_counts4_METeffect_infecting_chy)),
                    label = "both", label_alpha = 0, ) + 
  scale_fill_gradient(low="grey90",high = "red")+
  theme(legend.position = "none")+
  coord_sf(xlim = c(-.1, 1.1))

## Sanity check
lapply(names(contrast_chytridgenome), function(name) {
  x <- contrast_chytridgenome[[name]]  # access by name
  count_name <- names(listCounts)[grepl(gsub("res_", "", name), names(listCounts))]
  co <- listCounts[[count_name]]
  
  if(is.null(x)) return(NULL)  # skip NULL contrasts
  
  genes <- getGenes(x)
  if(length(genes) == 0) {
    message(name, ": no significant genes, skipping")
    return(NULL)
  }
  
  # Prepare data
  # Define desired order globally
  condition_levels <- c("control_chy", "control_both", "met_chy", "met_both",
                        "control_cyano", "met_cyano")  # include cyano too
  
  deg_counts <- co[rownames(co) %in% genes, ] %>%
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    pivot_longer(-gene, names_to="sample", values_to="count") %>%
    mutate(
      condition = sub("_[^_]+$", "", sample),
      replicate = sub(".*_", "", sample),
      condition = factor(condition, 
                         levels = condition_levels[condition_levels %in% condition])
      # only keeps levels present in the data, but in the right order #***
    )

  # Plot
  ggplot(deg_counts, aes(x=condition, y=count+1, color=condition)) +
    ggtitle(name) +                   
    geom_violin(alpha=0.3) +
    geom_boxplot(width=0.2, alpha=0.5) +
    geom_jitter(width=0.2, size=2, alpha=0.7) +
    geom_text(aes(label = sub(".*_", "", sample)))+
    scale_y_log10() +                  # log scale recommended for counts
    facet_wrap(~gene) +
    theme_bw() +
    theme(axis.text.x = element_text(angle=45, hjust=1),
          legend.position = "none")
})

p2 <- ggVennDiagram(x = list("Metolachlor effect\nuninfected cyanobacteria" = getGenes(contrast_cyanogenome$resr_met_effect_1org),
                             "Metolachlor effect\nduring infection" = getGenes(contrast_cyanogenome$resr_met_effect_2orgs)),
                    label = "both", label_alpha = 0, ) + 
  scale_fill_gradient(low="grey90",high = "red")+
  theme(legend.position = "none")+
  coord_sf(xlim = c(-8, 8))

pdf(here("figures/Fig4_venns.pdf"), width=6, height=9)
p1 + ggtitle("Eﬀect on chytrid gene expression")
dev.off()

## Save DEG tables
selectDEGenes <- function(x) x[!is.na(x$padj) & x$padj < 0.05, ]

format_DEG_table <- function(contrast_list, comparison_labels) {
  x <- lapply(seq_along(contrast_list[1:length(comparison_labels)]), function(i) {
    df <- contrast_list[[i]]
    deg <- selectDEGenes(df)
    if (nrow(deg) > 0) {
      data.frame(
        geneName       = rownames(deg),
        padj           = deg$padj,
        log2FoldChange = deg$log2FoldChange
      ) %>% mutate(comparison = comparison_labels[i])
    }
  })
  do.call(rbind, x)
}

## 1. Chytrid
names(contrast_chytridgenome)
contrast_chytridgenome_DEG <- format_DEG_table(
  contrast_chytridgenome,
  comparison_labels = c(
    "Infection effect on chytrid gene expression, in the absence of MET",
    "Infection effect on chytrid gene expression, in the presence of MET",
    "MET effect on chytrid gene expression, in free-living chytrid zoospores",
    "MET effect on chytrid gene expression, in infecting chytrid"))

## 2. Cyano
names(contrast_cyanogenome)
contrast_cyanogenome_DEG <- format_DEG_table(
  contrast_cyanogenome,
  comparison_labels = c(
    "Infection effect on cyanobacteria gene expression, in the absence of MET",
    "Infection effect on cyanobacteria gene expression, in the presence of MET",
    "MET effect on cyanobacteria gene expression, in uninfected cyanobacteria",
    "MET effect on cyanobacteria gene expression, in infected cyanobacteria"))

contrast_cyanogenome_DEG <- contrast_cyanogenome_DEG[
  grep("Infection effect", contrast_cyanogenome_DEG$comparison, invert = T),]

fullDEGTable <- rbind(contrast_chytridgenome_DEG, contrast_cyanogenome_DEG) %>%
  mutate(padj           = signif(padj, 2),
         log2FoldChange = signif(log2FoldChange, 2)) %>%
  arrange(geneName)

table(fullDEGTable$comparison)
# Infection effect on chytrid gene expression, in the absence of MET 
# 6 
# Infection effect on chytrid gene expression, in the presence of MET 
# 3 
# MET effect on chytrid gene expression, in free-living chytrid zoospores 
# 6 
# MET effect on chytrid gene expression, in infecting chytrid 
# 9 
# MET effect on cyanobacteria gene expression, in uninfected cyanobacteria 
# 1 

write.csv(fullDEGTable, "../../figures/TableS1_fullDEGTable.tsv", row.names=FALSE)

#############################################
## GO enrichment                           ##
#############################################

## Chytrid
getGOBubbleZ(universe  = rownames(counts1_infeffect_noMET_chy), 
             annotation = annotationChytrid_final,
             genelist  = getGenes(contrast_chytridgenome$res_counts1_infeffect_noMET_chy),
             GO_df     = GO_chytrid, isbubble = FALSE)

getGOBubbleZ(universe  = rownames(counts2_infeffect_MET_chy), 
             annotation = annotationChytrid_final,
             genelist  = getGenes(contrast_chytridgenome$res_counts2_infeffect_MET_chy),
             GO_df     = GO_chytrid, isbubble = FALSE)

getGOBubbleZ(universe  = rownames(counts3_METeffect_alone_chy), 
             annotation = annotationChytrid_final,
             genelist  = getGenes(contrast_chytridgenome$res_counts3_METeffect_alone_chy),
             GO_df     = GO_chytrid, isbubble = FALSE)

getGOBubbleZ(universe  = rownames(counts4_METeffect_infecting_chy), 
             annotation = annotationChytrid_final,
             genelist  = getGenes(contrast_chytridgenome$res_counts4_METeffect_infecting_chy),
             GO_df     = GO_chytrid, isbubble = FALSE)
# no significant GO terms (not enough DEG)