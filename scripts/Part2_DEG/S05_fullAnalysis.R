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

# Remove genes detected in only 1 to 3 sample - likely artifacts
n_samples_detected_chy <- rowSums(counts_chy > 0)
message("Chytrid genes by n samples detected:")            
print(table(n_samples_detected_chy))                      

counts_chy <- counts_chy[n_samples_detected_chy >= 4, ]
message("Chytrid genes after removing genes covered in 3 samples or less: ",  
        nrow(counts_chy))               
## 1349

## Same for cyano
counts_cyano <- as.matrix(RSEM_cyano)
counts_cyano <- counts_cyano[rowSums(counts_cyano) > 0, ]

n_samples_detected_cyano <- rowSums(counts_cyano > 0)    
message("Cyano genes by n samples detected:")             
print(table(n_samples_detected_cyano))                      

counts_cyano <- counts_cyano[n_samples_detected_cyano >= 4, ]
message("Cyano genes after removing genes covered in 3 samples or less: ",  
        nrow(counts_cyano))   
## 3753

# Saturation plot: number of genes detected at observed and simulated
# sequencing depths, and marginal gain in gene detection per million additional reads
eset_chy <- ExpressionSet(assayData = counts_chy)
mysaturation_chy <- dat(eset_chy, k = 0, ndepth = 7, type = "saturation")

eset_cyano <- ExpressionSet(assayData = counts_cyano)
mysaturation_cyano <- dat(eset_cyano, k = 0, ndepth = 7, type = "saturation")

## Plots:
pdf(here("figures/Fig_S2a.saturation_plots.pdf"), width=13, height=9)

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

# met_both_In12 control_both_In3 control_both_In5     met_both_In9     met_both_In8 
# 0.8774226        0.9162581        1.1230317        3.6461793        4.1771954 
# met_both_In7 control_both_In2 control_both_In1    met_both_In10 control_both_In6 
# 20.3524841       80.1941478      122.2305252      165.6350000      182.5740257 
# control_both_In4    met_both_In11 
# 361.4186156      446.0448654 

## ============================================================
## Visualise sequencing composition across all samples
## ============================================================

# For each sample: how many chytrid genes have >0 counts / total chytrid genes
# and same for cyano
table(rowSums(counts_chy) == 0)
table(rowSums(counts_cyano) == 0)

total_chy_genes   <- nrow(counts_chy)
total_cyano_genes <- nrow(counts_cyano)

# % of chytrid genes detected per sample
pct_chy_genes <- apply(counts_chy, 2, function(x) 
  100 * sum(x > 0) / total_chy_genes)

# % of cyano genes detected per sample  
pct_cyano_genes <- apply(counts_cyano, 2, function(x)
  100 * sum(x > 0) / total_cyano_genes)

all_samples <- union(names(pct_chy_genes), names(pct_cyano_genes))

# Add total genes detected per sample
n_chy_genes   <- apply(counts_chy,   2, function(x) sum(x > 0))
n_cyano_genes <- apply(counts_cyano, 2, function(x) sum(x > 0))

composition_df <- data.frame(
  sample      = all_samples,
  pct_chy     = pct_chy_genes[all_samples],
  pct_cyano   = pct_cyano_genes[all_samples],
  n_chy       = n_chy_genes[all_samples],    
  n_cyano     = n_cyano_genes[all_samples]   
) %>%
  replace(is.na(.), 0) %>%
  mutate(
    n_total     = n_chy + n_cyano,           
    sample_type = case_when(
      grepl("both",  sample) ~ "co-culture",
      grepl("chy",   sample) ~ "chytrid alone",
      grepl("cyano", sample) ~ "cyano alone"
    ),
    treatment = ifelse(grepl("^control", sample), "control", "MET"),
    label     = sub(".*_", "", sample)
  )

p_composition <- ggplot(composition_df,
                        aes(x=pct_cyano, y=pct_chy,
                            color=sample_type, shape=treatment)) +
  geom_point(size = 5) +
  ggrepel::geom_label_repel(aes(label=label), size=3, show.legend=FALSE) +
  scale_alpha_identity() +
  scale_color_manual(values=c("co-culture"    = "purple",
                              "chytrid alone" = "darkgreen",
                              "cyano alone"   = "steelblue")) +
  scale_shape_manual(values=c("control"=16, "MET"=17)) +
  annotate("text", x=90, y=90, label="balanced co-culture",
           color="grey50", size=3, fontface="italic") +
  annotate("text", x=10,  y=102, label="chytrid dominated",
           color="grey50", size=3, fontface="italic") +
  annotate("text", x=85, y=5,  label="cyano dominated",
           color="grey50", size=3, fontface="italic") +
  labs(x        = "% cyanobacterial genes detected",
       y        = "% chytrid genes detected",
       title    = "Sequencing composition per sample",
       color    = "Sample type",
       shape    = "Treatment") +
  theme_bw() +
  theme(legend.position="top")

composition_long <- composition_df %>%
  dplyr::select(sample, sample_type, treatment, label,
                n_chy, n_cyano) %>%
  pivot_longer(cols      = c(n_chy, n_cyano),
               names_to  = "organism",
               values_to = "n_genes") %>%
  mutate(
    organism = case_when(
      organism == "n_chy"   ~ "Chytrid",
      organism == "n_cyano" ~ "Cyanobacteria"
    ),
    # Remove rows where organism wasn't present in that sample type
    relevant = case_when(
      organism == "Chytrid"       & sample_type == "cyano alone"   ~ FALSE,
      organism == "Cyanobacteria" & sample_type == "chytrid alone" ~ FALSE,
      TRUE ~ TRUE
    )
  ) %>%
  filter(relevant)

p_ngenes <- ggplot(composition_long,
       aes(x    = reorder_within(label, n_genes, organism),  #***
           y    = n_genes,
           fill = sample_type)) +
  geom_col() +
  scale_x_reordered() +                                       #*** strips facet suffix
  scale_fill_manual(values = c("co-culture"    = "purple",
                               "chytrid alone" = "darkgreen",
                               "cyano alone"   = "steelblue")) +
  facet_wrap(~ organism, scales="free", ncol=1) +
  labs(x     = "Sample",
       y     = "Number of genes detected (count > 0)",
       title = "Genes detected per sample per organism",
       fill  = "Sample type") +
  theme_bw() +
  theme(axis.text.x  = element_text(angle=45, hjust=1),
        legend.position = "top")

pdf(here("figures/Fig_S3.genes_detected.pdf"), width=11, height=6)
cowplot::plot_grid(p_composition, p_ngenes)
dev.off()

# Build presence/absence matrix for each organism
pa_chy   <- (counts_chy   > 0) * 1  # 1 if expressed, 0 if not
pa_cyano <- (counts_cyano > 0) * 1

# Convert to list of genes present per sample
gene_lists_chy <- apply(pa_chy, 2, function(x)
  rownames(pa_chy)[x == 1])

upset(fromList(gene_lists_chy),
      nsets       = ncol(pa_chy),
      order.by    = "freq",
      main.bar.color = "darkgreen",
      sets.bar.color = "darkgreen",
      mainbar.y.label = "Genes in intersection",
      sets.x.label    = "Genes per sample")

gene_lists_cyano <- apply(pa_cyano, 2, function(x)
  rownames(pa_cyano)[x == 1])

upset(fromList(gene_lists_cyano),
      nsets       = ncol(pa_cyano),
      order.by    = "freq",
      main.bar.color = "steelblue",
      sets.bar.color = "steelblue",
      mainbar.y.label = "Genes in intersection",
      sets.x.label    = "Genes per sample")

# Sample          Chytrid matrix    Cyano matrix    Reason
# ────────────────────────────────────────────────────────────────
# met_chy_Z9      ❌ REMOVE         -               OClearr outlier
# control_both_In1 ❌ REMOVE        ❌ REMOVE       <100 cyano genes, low chytrid
# control_both_In4 keep             ❌ REMOVE       <100 cyano genes  
# control_both_In6 keep             ❌ REMOVE       <100 cyano genes
# met_both_In10   keep              ❌ REMOVE       <100 cyano genes
# met_both_In11   ❌ REMOVE         ❌ REMOVE       Outlier in both (high chytrid ratio)
# control_both_In2 keep             ❌ REMOVE   
# All others      ✅ KEEP           ✅ KEEP

## Based on (1) saturation plots and (2) ratio chytrid/cyano counts per sample,
## we remove the following samples:

# Chytrid matrix: met_chy_Z9 and control_both_In1 are clearly outliers (14.3%, 69.3%)
remove_chy <- c("met_chy_Z9", "control_both_In1", "met_both_In11")

# Cyano matrix: 2 different for cyano alone (keep all) vs co-culture (remove outliers)
remove_cyano <- c("control_both_In1", "control_both_In4",
                  "control_both_In6", "met_both_In10",
                  "met_both_In11", "control_both_In2") 

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
colnames(counts2_infeffect_MET_chy) # 5 vs 5

counts3_METeffect_alone_chy <- counts_chy[ ,grep("chy", colnames(counts_chy))]
counts3_METeffect_alone_chy <- counts3_METeffect_alone_chy[
  , !colnames(counts3_METeffect_alone_chy) %in% remove_chy]
colnames(counts3_METeffect_alone_chy) # 6 vs 5

counts4_METeffect_infecting_chy <- counts_chy[ ,grep("both", colnames(counts_chy))]
counts4_METeffect_infecting_chy <- counts4_METeffect_infecting_chy[
  , !colnames(counts4_METeffect_infecting_chy) %in% remove_chy]
colnames(counts4_METeffect_infecting_chy) # 5 vs 5

counts5_infeffect_noMET_cyano <- counts_cyano[ ,grep("control", colnames(counts_cyano))]
counts5_infeffect_noMET_cyano <- counts5_infeffect_noMET_cyano[
  , !colnames(counts5_infeffect_noMET_cyano) %in% remove_cyano]
colnames(counts5_infeffect_noMET_cyano) # 2 vs 6 NOT RELIABLE!!

counts6_infeffect_MET_cyano <- counts_cyano[ ,grep("met", colnames(counts_cyano))]
counts6_infeffect_MET_cyano <- counts6_infeffect_MET_cyano[
  , !colnames(counts6_infeffect_MET_cyano) %in% remove_cyano]
colnames(counts6_infeffect_MET_cyano) # 4 vs 6

counts7_METeffect_alone_cyano <- counts_cyano[ ,grep("cyano", colnames(counts_cyano))]
counts7_METeffect_alone_cyano <- counts7_METeffect_alone_cyano[
  , !colnames(counts7_METeffect_alone_cyano) %in% remove_cyano]
colnames(counts7_METeffect_alone_cyano) # 6 vs 6

counts8_METeffect_infecting_cyano <- counts_cyano[ ,grep("both", colnames(counts_cyano))]
counts8_METeffect_infecting_cyano <- counts8_METeffect_infecting_cyano[
  , !colnames(counts8_METeffect_infecting_cyano) %in% remove_cyano]
colnames(counts8_METeffect_infecting_cyano) # 2 vs 4 NOT RELIABLE!!

## Reapply filtering
# Create named list automatically from object names
listCounts <- mget(c("counts1_infeffect_noMET_chy",
                     "counts2_infeffect_MET_chy",
                     "counts3_METeffect_alone_chy",
                     "counts4_METeffect_infecting_chy",
                     "counts7_METeffect_alone_cyano"))

listCounts <- lapply(listCounts, filter_by_group)

# Automated summary for all count matrices in the list
invisible(lapply(names(listCounts), function(name) {
  x <- listCounts[[name]]
  message(name, ":  ", ncol(x), " samples, ", nrow(x), " genes")
}))
## invisible() suppresses the NULL list output that lapply would otherwise print to console.

# counts1_infeffect_noMET_chy:  11 samples, 837 genes
# counts2_infeffect_MET_chy:  10 samples, 945 genes
# counts3_METeffect_alone_chy:  11 samples, 941 genes
# counts4_METeffect_infecting_chy:  10 samples, 829 genes
# counts7_METeffect_alone_cyano:  12 samples, 3553 genes

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

lapply(listCounts, testEffectSize)

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
# out of 945 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 2, 0.21%
# LFC < 0 (down)     : 1, 0.11%
# outliers [1]       : 34, 3.6%
# low counts [2]     : 0, 0%
# (mean count < 6)

testIfCrazyFoldChange(res_counts2_infeffect_MET_chy) # OK

## --- 3. MET effect in chytrid alone ---
res_counts3_METeffect_alone_chy <- myDESeq2_withlfcShrink(
  listCounts$counts3_METeffect_alone_chy, a="control_chy", b="met_chy")

# out of 941 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 4, 0.43%
# LFC < 0 (down)     : 2, 0.21%
# outliers [1]       : 19, 2%
# low counts [2]     : 0, 0%
# (mean count < 10)

testIfCrazyFoldChange(res_counts3_METeffect_alone_chy) # ok

## --- 4. MET effect in chytrid infecting cyanobacteria
res_counts4_METeffect_infecting_chy <- myDESeq2_withlfcShrink(
  listCounts$counts4_METeffect_infecting_chy, a="control_both", b="met_both")
# out of 829 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 5, 0.6%
# LFC < 0 (down)     : 4, 0.48%
# outliers [1]       : 34, 4.1%
# low counts [2]     : 0, 0%
# (mean count < 12)

testIfCrazyFoldChange(res_counts4_METeffect_infecting_chy)

## --- 5. Infection effect in cyano, no MET --- ***** UNRELIABLE

## --- 6. Infection effect in cyano, with MET --- ***** UNRELIABLE

## --- 7. MET effect in cyano alone ---
res_counts7_METeffect_alone_cyano <- myDESeq2_withlfcShrink(
  listCounts$counts7_METeffect_alone_cyano, a="control_cyano", b="met_cyano")
# out of 3553 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 0, 0%
# LFC < 0 (down)     : 1, 0.028%
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

## --- 8. MET effect in cyano infected by chytrid ---***** UNRELIABLE

## --- Combine into named lists ---
contrast_chytridgenome <- mget(c(
  "res_counts1_infeffect_noMET_chy",  
  "res_counts2_infeffect_MET_chy",
  "res_counts3_METeffect_alone_chy",
  "res_counts4_METeffect_infecting_chy"
))

contrast_cyanogenome <- mget(c(
  "res_counts7_METeffect_alone_cyano"  # MET effect, cyano alone
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
  positionLogoStart = -2, positionLogoStop = 2,
  mylogo       = "logos/logo1.png")

V_chytrid_inf_effect_met <- makeVolcano(
  res          = contrast_chytridgenome$res_counts2_infeffect_MET_chy,
  title        = "MET, chytrids zoospores vs chytrids infecting",
  positionLogoStart = -2, positionLogoStop = 2,
  mylogo       = "logos/logo2.png")

V_chytrid_met_effect_1org <- makeVolcano(
  res          = contrast_chytridgenome$res_counts3_METeffect_alone_chy,
  title        = "free-living zoospores, no MET vs MET",
  positionLogoStart = -2, positionLogoStop = 2,
  mylogo       = "logos/logo3.png")

V_chytrid_met_effect_2orgs <- makeVolcano(
  res          = contrast_chytridgenome$res_counts4_METeffect_infecting_chy,
  title        = "chytrids infecting, no MET vs MET",
  positionLogoStart = -2, positionLogoStop = 2,
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

# V_cyano_met_effect_2orgs <- makeVolcano(
#   res          = contrast_cyanogenome$res_counts8_METeffect_infecting_cyano,
#   title        = "infected cyanobacteria, no MET vs MET",
#   mylogo       = "logos/logo8.png")

## Venn diagrams
getGenes <- function(x) rownames(x[!is.na(x$padj) & x$padj < 0.05, ])

p1 <- ggVennDiagram(x = list("Infection effect\nabsence of metolachlor" = 
                               getGenes(contrast_chytridgenome$res_counts1_infeffect_noMET_chy),
                             "Infection effect\npresence of metolachlor" = 
                               getGenes(contrast_chytridgenome$res_counts2_infeffect_MET_chy),
                             "Metolachlor effect\nfree living zoospores" = 
                               getGenes(contrast_chytridgenome$res_counts3_METeffect_alone_chy),
                             "Metolachlor effect\nduring infection"     = 
                               getGenes(contrast_chytridgenome$res_counts4_METeffect_infecting_chy)),
                    label = "both", label_alpha = 0, ) + 
  scale_fill_gradient(low="grey90",high = "red")+
  theme(legend.position = "none")+
  coord_sf(xlim = c(-.1, 1.1))

pdf(here("figures/Fig2.pdf"), width=10, height=10)
cowplot::plot_grid(V_chytrid_inf_effect_control$plot,
                   V_chytrid_inf_effect_met$plot,
                   V_chytrid_met_effect_1org$plot,
                   V_chytrid_met_effect_2orgs$plot,
                   p1 + ggtitle("Eﬀect on chytrid gene expression"),
                   V_cyano_met_effect_1org$plot,
                   labels=c("a","b","c","d", "e", "f"), nrow = 3,
                   label_size=20)
dev.off()

####################################################################
## Raw count plots for DEGs - one plot per contrast, saved to PDF ##
####################################################################
condition_levels <- c("control_chy", "control_both", "met_chy", "met_both",
                      "control_cyano", "met_cyano")

# Step 1: Collect ALL DEGs across all contrasts with their significance info
all_deg_data <- lapply(names(contrast_chytridgenome), function(name) {
  x          <- contrast_chytridgenome[[name]]
  count_name <- names(listCounts)[grepl(gsub("res_", "", name), names(listCounts))]
  co         <- listCounts[[count_name]]
  if(is.null(x)) return(NULL)
  genes <- getGenes(x)
  if(length(genes) == 0) return(NULL)

  # Which two conditions are compared in this contrast?
  contrast_conditions <- list(
    res_counts1_infeffect_noMET_chy     = c("control_chy",  "control_both"),
    res_counts2_infeffect_MET_chy       = c("met_chy",      "met_both"),
    res_counts3_METeffect_alone_chy     = c("control_chy",  "met_chy"),
    res_counts4_METeffect_infecting_chy = c("control_both", "met_both")
  )
  conds <- contrast_conditions[[name]]
  
  # Get padj for each significant gene
  sig_info <- as.data.frame(x) %>%
    filter(!is.na(padj) & padj < 0.05) %>%
    rownames_to_column("gene") %>%
    dplyr::select(gene, padj) %>%
    mutate(cond1 = conds[1], cond2 = conds[2])
  
  # Get counts for ALL 4 conditions for these genes
  all_counts <- counts_chy[rownames(counts_chy) %in% genes, ] %>%
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    pivot_longer(-gene, names_to="sample", values_to="count") %>%
    mutate(condition = sub("_[^_]+$", "", sample))
  list(counts=all_counts, sig=sig_info)
}) %>% Filter(Negate(is.null), .)

# Combine counts and significance across all contrasts
all_counts_combined <- bind_rows(lapply(all_deg_data, `[[`, "counts")) %>%
  distinct() %>%
  mutate(condition = factor(condition, levels=condition_levels))

all_sig_combined <- bind_rows(lapply(all_deg_data, `[[`, "sig")) %>%
  mutate(stars = case_when(
    padj < 0.001 ~ "***",
    padj < 0.01  ~ "**",
    padj < 0.05  ~ "*"
  ))

# Order genes alphabetically
all_genes <- sort(unique(all_counts_combined$gene))

message("Total unique DEGs across all contrasts: ", length(all_genes))

# Step 2: One plot per gene, all conditions, with significance bars
plots_pergene <- lapply(all_genes, function(g) {
  
  df <- all_counts_combined %>% filter(gene == g)
  sig_df <- all_sig_combined %>% filter(gene == g)
  
  # Build significance annotations for ggsignif
  # Need y positions staggered
  y_max <- max(log10(df$count + 1), na.rm=TRUE)
  sig_df <- sig_df %>%
    mutate(y_pos = y_max + 0.3 * row_number())
  
  p <- ggplot(df, aes(x=condition, y=count+1, color=condition)) +
    geom_violin(alpha=0.3, aes(fill=condition)) +
    geom_boxplot(width=0.15, alpha=0.7, outlier.shape=NA) +
    geom_jitter(width=0.15, size=1.5, alpha=0.8) +
    scale_y_log10() +
    scale_color_manual(values=c("control_chy"  = "darkgreen",
                                "met_chy"      = "lightgreen",
                                "control_both" = "purple",
                                "met_both"     = "plum")) +
    scale_fill_manual(values=c("control_chy"  = "darkgreen",
                               "met_chy"      = "lightgreen",
                               "control_both" = "purple",
                               "met_both"     = "plum")) +
    ggtitle(sub("_.*", "", g))+ #  simplify gene ID name for cyano) +
    labs(x=NULL, y="Count + 1 (log10)") +
    coord_cartesian(clip="off") +          #*** allow drawing outside plot area
    theme_minimal() +
    theme(axis.text.x  = element_text(angle=45, hjust=1),
          legend.position = "none",
          plot.title     = element_text(size=9, face="bold"))
  
  # Add significance bars
  if(nrow(sig_df) > 0) {
    for(i in seq_len(nrow(sig_df))) {
      p <- p + ggsignif::geom_signif(
        comparisons = list(c(sig_df$cond1[i], sig_df$cond2[i])),
        annotations = sig_df$stars[i],
        y_position  = sig_df$y_pos[i],
        tip_length  = 0.02,
        color       = "black",
        size        = 0.4,
        textsize    = 4)
    }
  }
  p
})
names(plots_pergene) <- all_genes

# Step 3: Save - auto-calculate grid dimensions
n_plots <- length(plots_pergene)
ncols   <- 6
nrows   <- ceiling(n_plots / ncols)

pdf(here("figures/Fig_rawcounts_DEG.pdf"),
    width  = ncols * 2.5,
    height = nrows * 3.5)
cowplot::plot_grid(plotlist = plots_pergene,
                   ncol     = ncols)
dev.off()

message("Saved ", n_plots, " gene plots (", nrows, " rows x ", ncols, " cols)")

####################
## Save DEG table ##
####################
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
    "MET effect on cyanobacteria gene expression, in uninfected cyanobacteria"))

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

## Add annotation
fullDEGTable$geneFull <- annotationChytrid_final$sprot_Top_BLASTX_hit[
  match(fullDEGTable$geneName, annotationChytrid_final$gene_name)]

fullDEGTable$geneFull <- str_extract(fullDEGTable$geneFull, "(?<=Full=)[^;]+")

x <- annotationCyano_final$protein[
  match(fullDEGTable$geneName, annotationCyano_final$gene_id)]

fullDEGTable$geneFull[which(is.na(fullDEGTable$geneFull))] <- 
  x[which(is.na(fullDEGTable$geneFull))]
rm(x)

fullDEGTable$potential_explanation <- case_when(
  fullDEGTable$geneName == "77286457" ~
    "Iron uptake porin upregulated under MET - iron acquisition may be disrupted by MET-induced oxidative stress (MET inhibits fatty acid synthesis, indirectly affecting membrane integrity and iron transport)",
  
  fullDEGTable$geneName == "AB6C_ARATH" ~
    "ABC-C transporter (MRP family) - fungal ABC transporters transport xenobiotics including herbicides and fungicides across membranes; upregulation likely reflects active efflux of MET as a detoxification response",
  
  fullDEGTable$geneName == "ABCB6_HUMAN" ~
    "Mitochondrial porphyrin/heme transporter - ABCB6 regulates heme biosynthesis and reduces reactive oxygen species; downregulation during infection may reflect redirection of mitochondrial resources toward parasite growth",
  
  fullDEGTable$geneName %in% c("ACA1_SCHPO", "ACA1_SCHPO.1") ~
    "N-acetyltransferase - involved in histone acetylation and gene regulation; differential expression under MET suggests epigenetic reprogramming in response to herbicide stress",
  
  fullDEGTable$geneName == "AL1A3_HUMAN" ~
    "Retinaldehyde dehydrogenase 3 (ALDH) - oxidises aldehydes generated by oxidative stress; upregulation under MET consistent with oxidative stress response; SIK3-deficient mice also show low ALDH1a expression, linking this gene to the same energy-stress pathway as SIK3",
  
  fullDEGTable$geneName == "AMYA_DROMA" ~
    "Alpha-amylase - starch/glycogen hydrolysis enzyme; upregulation during infection may reflect increased carbohydrate mobilisation to fuel parasite growth and zoospore production",
  
  fullDEGTable$geneName == "ARF_CRYNB" ~
    "ADP-ribosylation factor - key regulator of vesicle trafficking and cytoskeletal dynamics; downregulation during MET-treated infection suggests impaired vesicle-mediated secretion, consistent with MET disrupting membrane lipid composition required for vesicle formation",
  
  fullDEGTable$geneName == "ATC2_CRYNH" ~
    "Calcium-transporting ATPase - maintains intracellular Ca²⁺ homeostasis; downregulation during infection may reflect suppression of calcium signalling required for zoospore motility and host recognition",
  
  fullDEGTable$geneName == "BIP_ASPAW" ~
    "ER chaperone BiP (Hsp70 family) - canonical unfolded protein response (UPR) marker; upregulation under MET indicates protein folding stress, consistent with MET disrupting lipid synthesis and membrane protein biogenesis in the ER",
  
  fullDEGTable$geneName == "CAMKI_MACNP" ~
    "Calcium/calmodulin-dependent protein kinase I - mediates calcium-dependent signalling; upregulation during infection suggests activation of Ca²⁺ signalling cascades during host-parasite interaction, possibly linked to zoospore encystment and rhizoid development",
  
  fullDEGTable$geneName %in% c("CBPC1_HUMAN", "CBPC1_HUMAN.1") ~
    "Cytosolic carboxypeptidase 1 - deglutamylase involved in microtubule modification; downregulation is consistent with reduced tubulin post-translational modification, potentially impairing flagellar/ciliary function critical for zoospore motility",
  
  fullDEGTable$geneName == "DBP4_CRYNJ" ~
    "ATP-dependent RNA helicase DBP4 - involved in rRNA processing and ribosome biogenesis; upregulation under MET may reflect compensatory increase in translational machinery under stress conditions",
  
  fullDEGTable$geneName == "FEN1_NEMVE" ~
    "Flap endonuclease 1 - DNA repair enzyme involved in Okazaki fragment processing and base excision repair; downregulation during MET-treated infection suggests reduced DNA replication/repair capacity, possibly reflecting cell cycle arrest",
  
  fullDEGTable$geneName == "MOGT1_HUMAN" ~
    "2-acylglycerol O-acyltransferase 1 - catalyses monoacylglycerol re-acylation in lipid metabolism; downregulation during MET infection aligns with MET's primary mode of action as a VLCFAS inhibitor, further suppressing lipid metabolic capacity",
  
  fullDEGTable$geneName == "NPAL2_HUMAN" ~
    "NIPA-like protein 2 - magnesium transporter; upregulation under MET may reflect compensatory ion transport in response to membrane disruption caused by VLCFAS inhibition",
  
  fullDEGTable$geneName == "NU4M_USTMA" ~
    "NADH-ubiquinone oxidoreductase chain 4 (Complex I) - core mitochondrial respiratory chain subunit; upregulation during infection suggests increased oxidative phosphorylation to meet energy demands of parasite growth and host penetration",
  
  fullDEGTable$geneName == "PAQR1_HUMAN" ~
    "Adiponectin receptor 1 - membrane receptor linked to fatty acid oxidation and AMPK activation; upregulation under MET may represent a compensatory response to MET-induced disruption of fatty acid synthesis, activating alternative lipid catabolism pathways",
  
  fullDEGTable$geneName == "PHOP1_MOUSE" ~
    "Phosphoethanolamine/phosphocholine phosphatase - involved in phospholipid head group metabolism; upregulation during infection suggests active remodelling of membrane phospholipid composition during host-parasite interaction",
  
  fullDEGTable$geneName == "PPCE_MOUSE" ~
    "Prolyl endopeptidase - serine protease cleaving proline-containing peptides; upregulation under MET in infecting chytrid may reflect increased proteolytic activity during host cell degradation",
  
  fullDEGTable$geneName == "SIK3_HUMAN" ~
    "Salt-inducible kinase 3 (AMPK-related kinase) - central regulator of lipid and energy metabolism via LKB1-SIK3-HDAC4 axis; strongly upregulated under MET during infection (LFC=3.9, padj=6.7e-5), consistent with MET's VLCFAS inhibition triggering energy stress signalling; SIK3-deficient mice show lipodystrophy and reduced fatty acid synthesis, mirroring MET's mechanism of action",
  
  fullDEGTable$geneName %in% c("SYKC_SCHPO", "SYKC_SCHPO.1") ~
    "Lysine-tRNA ligase, cytoplasmic - aminoacyl-tRNA synthetase essential for translation; differential expression may reflect altered protein synthesis capacity; also a SIK family homolog in S. pombe, suggesting conserved kinase-related stress response",
  
  fullDEGTable$geneName == "TRNL_SCHPO" ~
    "tRNA ligase 1 - involved in tRNA splicing and the unfolded protein response (UPR); upregulation during infection suggests activation of the tRNA ligase-dependent UPR branch, consistent with BIP_ASPAW upregulation indicating ER stress",
  
  TRUE ~ NA_character_
)

write.csv(fullDEGTable, here("figures/TableS1_fullDEGTable_annotated.csv"),
          row.names=FALSE)

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
