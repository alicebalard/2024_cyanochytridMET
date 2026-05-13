## ============================================================
## Part 3: Co-expression & Cross-species analysis
## ============================================================
library(here)
source(here("scripts/Part2_DEG/S05_fullAnalysis.R"))
setwd(here("scripts/Part2_WGCNA/"))

library(WGCNA)
library(igraph); library(tidygraph); library(ggraph)
library(purrr); library(patchwork); library(openxlsx)

options(stringsAsFactors=FALSE)
allowWGCNAThreads()

## ============================================================
## 1. Compute VST matrices
## ============================================================

getVST <- function(count) {
  ddsr <- DESeqDataSetFromMatrix(
    countData = round(count),
    colData   = samples_data[colnames(count), ],
    design    = ~ condition)
  assay(varianceStabilizingTransformation(ddsr, blind=TRUE))
}

# VST for single-species analyses (MET effect contrasts only)
all_chy_genes <- Reduce(union, lapply(
  listCounts[intersect(grep("chy",      names(listCounts)),
                       grep("METeffect", names(listCounts)))], rownames))

all_cyano_genes <- Reduce(union, lapply(
  listCounts[intersect(grep("cyano",    names(listCounts)),
                       grep("METeffect", names(listCounts)))], rownames))

# Chytrid alone samples
counts_chy_forVST <- counts_chy[
  rownames(counts_chy) %in% all_chy_genes,
  grep("^control_chy|^met_chy", colnames(counts_chy))]
counts_chy_forVST <- counts_chy_forVST[rowSums(counts_chy_forVST) > 0, ]
vst_chy <- getVST(counts_chy_forVST)
message("VST chy:   ", nrow(vst_chy), " genes x ", ncol(vst_chy), " samples")
# VST chy:   958 genes x 12 samples

# Cyano alone samples
counts_cyano_forVST <- counts_cyano[
  rownames(counts_cyano) %in% all_cyano_genes,
  grep("cyano", colnames(counts_cyano))]
counts_cyano_forVST <- counts_cyano_forVST[rowSums(counts_cyano_forVST) > 0, ]
vst_cyano <- getVST(counts_cyano_forVST)
message("VST cyano: ", nrow(vst_cyano), " genes x ", ncol(vst_cyano), " samples")
# VST cyano: 3553 genes x 12 samples

# VST for cross-species analysis (co-culture "both" samples)
counts_chy_both_forVST <- counts_chy[
  rownames(counts_chy) %in% all_chy_genes,
  grep("both", colnames(counts_chy))]
counts_chy_both_forVST <- counts_chy_both_forVST[
  , !colnames(counts_chy_both_forVST) %in% remove_chy]
counts_chy_both_forVST <- counts_chy_both_forVST[
  rowSums(counts_chy_both_forVST) > 0, ]

counts_cyano_both_forVST <- counts_cyano[
  rownames(counts_cyano) %in% all_cyano_genes,
  grep("both", colnames(counts_cyano))]
counts_cyano_both_forVST <- counts_cyano_both_forVST[
  , !colnames(counts_cyano_both_forVST) %in% remove_cyano]
counts_cyano_both_forVST <- counts_cyano_both_forVST[
  rowSums(counts_cyano_both_forVST) > 0, ]

vst_chy_both   <- getVST(counts_chy_both_forVST)
vst_cyano_both <- getVST(counts_cyano_both_forVST)

# Keep only shared co-culture samples
shared_samples <- intersect(colnames(vst_chy_both), colnames(vst_cyano_both))
message("Shared co-culture samples: ", length(shared_samples))
# Shared co-culture samples: 7

vst_chy_both   <- vst_chy_both[,   shared_samples]
vst_cyano_both <- vst_cyano_both[, shared_samples]
stopifnot(all(colnames(vst_chy_both) == colnames(vst_cyano_both)))

## ============================================================
## 2. Single-species WGCNA
## ============================================================

pickSoftPow <- function(vst) {
  par(mfrow=c(1, 2))
  sth <- pickSoftThreshold(t(vst), powerVector=1:20, verbose=5)
  plot(sth$fitIndices$Power, sth$fitIndices$SFT.R.sq, type="b",
       xlab="Power", ylab="Scale Free R²", main="Scale-Free Fit")
  abline(h=0.8, col="red", lty=2); abline(h=0.9, col="red", lty=1)
  plot(sth$fitIndices$Power, sth$fitIndices$mean.k., type="b",
       xlab="Power", ylab="Mean connectivity", main="Mean Connectivity")
  abline(h=5, col="red", lty=2)
  par(mfrow=c(1, 1))
  invisible(sth)
}

pickSoftPow(vst_chy)
softPower_chy <- 9
# Power 9:  SFT.R.sq=0.828, mean.k=20.5 → R²>0.8 ✅, good connectivity ✅
#
pickSoftPow(vst_cyano)
softPower_cyano <- 5  # R²=0.837, mean.k=59.5, slope=-2.76 ✅

cor <- WGCNA::cor  # override base cor

net_chy <- blockwiseModules(
  datExpr        = t(vst_chy),
  power          = softPower_chy,
  TOMType        = "signed", minModuleSize=30,
  reassignThreshold=0, mergeCutHeight=0.25,
  saveTOMs=TRUE,  saveTOMFileBase="chyTOM", verbose=3)

net_cyano <- blockwiseModules(
  datExpr        = t(vst_cyano),
  power          = softPower_cyano,
  TOMType        = "signed", minModuleSize=30,
  reassignThreshold=0, mergeCutHeight=0.25,
  saveTOMs=TRUE,  saveTOMFileBase="cyanoTOM", verbose=3)

net_chy$colors[net_chy$colors     == "grey"] <- "white"
net_cyano$colors[net_cyano$colors == "grey"] <- "white"

message("Chytrid modules: ");  print(table(net_chy$colors))
message("Cyano modules: ");    print(table(net_cyano$colors))

# Dendrograms
pdf(here("figures/Fig5_dendrograms.pdf"), width=15, height=5)
plotDendroAndColors(net_chy$dendrograms[[1]],
                    net_chy$colors[net_chy$blockGenes[[1]]], "Module colors",
                    dendroLabels=FALSE, hang=0.03, addGuide=TRUE, guideHang=0.05,
                    main="Chytrid co-expression modules")
plotDendroAndColors(net_cyano$dendrograms[[1]],
                    net_cyano$colors[net_cyano$blockGenes[[1]]], "Module colors",
                    dendroLabels=FALSE, hang=0.03, addGuide=TRUE, guideHang=0.05,
                    main="Cyano co-expression modules")
dev.off()

## ============================================================
## 3. Single-species: module-trait correlation with MET
## ============================================================

moduleTraitTest <- function(net, vst, trait_label="MET") {
  treatment <- ifelse(grepl("^control", colnames(vst)), 0, 1)
  MEs  <- net$MEs
  cors <- cor(MEs, treatment, use="p")
  pval <- corPvalueStudent(cors, length(treatment))
  
  sig <- data.frame(
    module  = rownames(cors),
    cor     = round(cors[,1], 3),
    pvalue  = round(pval[,1], 4)
  ) %>% filter(pvalue < 0.05) %>% arrange(pvalue)
  
  message(trait_label, " - significant modules:")
  print(sig)
  
  labeledHeatmap(Matrix=cors, xLabels=trait_label,
                 yLabels=names(MEs), colors=blueWhiteRed(50),
                 main=paste("Module-trait:", trait_label))
  invisible(list(cors=cors, pval=pval, sig=sig))
}

pdf(here("figures/Fig5_moduleTraitHeatmaps.pdf"), width=4, height=6)
res_chy_MET   <- moduleTraitTest(net_chy,   vst_chy,   "MET (chytrid)")
res_cyano_MET <- moduleTraitTest(net_cyano, vst_cyano, "MET (cyano)")
dev.off()

# Eigengene plots for significant modules
plotSigEigengenes <- function(net, vst, sig_df, title) {
  if(nrow(sig_df) == 0) { message("No significant modules for: ", title); return(NULL) }
  meta <- data.frame(
    Sample    = colnames(vst),
    Treatment = sub("_.*", "", colnames(vst)))
  plot_df <- net$MEs %>%
    as.data.frame() %>%
    mutate(Sample=rownames(.)) %>%
    left_join(meta, by="Sample") %>%
    pivot_longer(starts_with("ME"), names_to="Module", values_to="Eigengene") %>%
    filter(Module %in% sig_df$module)
  ggplot(plot_df, aes(x=Treatment, y=Eigengene, fill=Treatment)) +
    geom_boxplot(alpha=0.7) +
    geom_jitter(width=0.2, size=2) +
    facet_wrap(~Module, scales="free_y") +
    theme_minimal() + theme(legend.position="none") +
    labs(title=title, y="Module Eigengene")
}

pdf(here("figures/Fig5_eigengenes.pdf"), width=8, height=4)
print(plotSigEigengenes(net_chy,   vst_chy,
                        res_chy_MET$sig,   "Chytrid modules correlated with MET"))
print(plotSigEigengenes(net_cyano, vst_cyano,
                        res_cyano_MET$sig, "Cyano modules correlated with MET"))
dev.off()

# Which DEGs fall in which module?
message("Chytrid DEGs in modules:")
print(data.frame(module=na.omit(
  net_chy$colors[unique(fullDEGTable$geneName)])) %>% arrange(module))

message("Cyano DEGs in modules:")
print(data.frame(module=na.omit(
  net_cyano$colors[unique(fullDEGTable$geneName)])) %>% arrange(module))
# module
# AB6C_ARATH       blue
# AMYA_DROMA       blue
# BIP_ASPAW        blue
# NU4M_USTMA       blue
# ABCB6_HUMAN turquoise
# ACA1_SCHPO  turquoise
# AL1A3_HUMAN turquoise
# ARF_CRYNB   turquoise
# CAMKI_MACNP turquoise
# CBPC1_HUMAN turquoise
# DBP4_CRYNJ  turquoise
# FEN1_NEMVE  turquoise
# MOGT1_HUMAN turquoise
# NPAL2_HUMAN turquoise
# PAQR1_HUMAN turquoise
# PHOP1_MOUSE turquoise
# PPCE_MOUSE  turquoise
# SIK3_HUMAN  turquoise
# SYKC_SCHPO  turquoise
# TRNL_SCHPO  turquoise
# ATC2_CRYNH      white

## ============================================================
## 4. Cross-species correlation (WGCNA not feasible: n=7 samples,
##    R² < 0.57 at all powers)
## ============================================================

# Spearman correlation matrix: chytrid genes x cyano genes
cor_matrix <- cor(t(vst_chy_both), t(vst_cyano_both), method="spearman")
message("Correlation matrix: ", nrow(cor_matrix), " x ", ncol(cor_matrix))
# Correlation matrix: 958 x 3500

# Extract strongly correlated pairs
high_cor <- which(abs(cor_matrix) > 0.9, arr.ind=TRUE)
cross_cor_df <- data.frame(
  chytrid_gene = rownames(cor_matrix)[high_cor[,1]],
  cyano_gene   = colnames(cor_matrix)[high_cor[,2]],
  spearman_r   = cor_matrix[high_cor]
) %>% arrange(desc(abs(spearman_r)))

message("Cross-species pairs |r|>0.9: ", nrow(cross_cor_df))

# Flag zero-inflated spurious correlations
cross_cor_df <- cross_cor_df %>%
  rowwise() %>%
  mutate(
    n_nonzero_chy   = sum(vst_chy_both[chytrid_gene, ]   > 0),
    n_nonzero_cyano = sum(vst_cyano_both[cyano_gene, ]   > 0),
    min_nonzero     = min(n_nonzero_chy, n_nonzero_cyano),
    suspicious      = (abs(spearman_r) >= 0.99 & min_nonzero <= 3)
  ) %>%
  ungroup()

message("Suspicious (zero-inflated) pairs: ", sum(cross_cor_df$suspicious))
cross_cor_df <- cross_cor_df %>% filter(!suspicious)
message("Clean pairs remaining: ", nrow(cross_cor_df))

# Add cyano gene annotation
cross_cor_df <- cross_cor_df %>%
  left_join(annotationCyano_final %>%
              distinct(gene_id, locus_tag, protein, product),
            by=c("cyano_gene"="gene_id"))

# Filter to pairs involving DEGs
cross_cor_DEG <- cross_cor_df %>%
  filter(chytrid_gene %in% fullDEGTable$geneName |
           cyano_gene   %in% fullDEGTable$geneName)

message("DEG-involving cross-species pairs: ", nrow(cross_cor_DEG))
print(cross_cor_DEG)

## ============================================================
## 5. Cross-species: correlation with MET treatment
## ============================================================

# For each chytrid DEG, test if its cyano partners differ by MET
treatment_both <- ifelse(grepl("^control", shared_samples), 0, 1)

# Correlate each chytrid DEG's expression with MET treatment
chy_DEGs_in_both <- intersect(fullDEGTable$geneName, rownames(vst_chy_both))

chy_MET_cor <- data.frame(
  chytrid_gene = chy_DEGs_in_both,
  cor_with_MET = sapply(chy_DEGs_in_both, function(g)
    cor(vst_chy_both[g, ], treatment_both, method="spearman")),
  pval_MET = sapply(chy_DEGs_in_both, function(g)
    cor.test(vst_chy_both[g, ], treatment_both,
             method="spearman")$p.value)
) %>% arrange(pval_MET)

message("Chytrid DEGs correlated with MET in co-culture (p<0.05):")
print(chy_MET_cor %>% filter(pval_MET < 0.05))

# Same for cyano DEGs
cyano_DEGs_in_both <- intersect(fullDEGTable$geneName, rownames(vst_cyano_both))

cyano_MET_cor <- data.frame(
  cyano_gene   = cyano_DEGs_in_both,
  cor_with_MET = sapply(cyano_DEGs_in_both, function(g)
    cor(vst_cyano_both[g, ], treatment_both, method="spearman")),
  pval_MET = sapply(cyano_DEGs_in_both, function(g)
    cor.test(vst_cyano_both[g, ], treatment_both,
             method="spearman")$p.value)
) %>% arrange(pval_MET)

message("Cyano DEGs correlated with MET in co-culture (p<0.05):")
print(cyano_MET_cor %>% filter(pval_MET < 0.05))

## ============================================================
## 6. Visualise key cross-species pairs
## ============================================================

plot_cross_cor <- function(gene_chy, gene_cyano) {
  cyano_info <- annotationCyano_final %>%
    filter(gene_id==gene_cyano) %>%
    slice(1) %>%
    pull(product)
  cyano_label <- ifelse(is.na(cyano_info), gene_cyano,
                        paste0(gene_cyano, "\n(", cyano_info, ")"))
  df <- data.frame(
    chytrid   = vst_chy_both[gene_chy, ],
    cyano     = vst_cyano_both[gene_cyano, ],
    sample    = shared_samples,
    treatment = ifelse(grepl("^control", shared_samples), "control", "MET")
  )
  ggplot(df, aes(x=chytrid, y=cyano, color=treatment,
                 label=sub(".*_", "", sample))) +
    geom_point(size=3) +
    ggrepel::geom_text_repel(size=3) +
    geom_smooth(method="lm", se=TRUE, color="grey40", alpha=0.2) +
    scale_color_manual(values=c("control"="blue", "MET"="red")) +
    labs(title = paste(gene_chy, "↔", gene_cyano),
         subtitle = paste("Spearman r =",
                          round(cor(df$chytrid, df$cyano,
                                    method="spearman"), 3)),
         x=paste("VST:", gene_chy),
         y=cyano_label) +
    theme_bw()
}

# Plot top DEG-involving pairs
top_DEG_pairs <- head(cross_cor_DEG, 12)

pdf(here("figures/Fig6_cross_species_correlations.pdf"), width=12, height=8)
plots <- lapply(1:min(9, nrow(top_DEG_pairs)), function(i)
  plot_cross_cor(top_DEG_pairs$chytrid_gene[i],
                 top_DEG_pairs$cyano_gene[i]))
print(cowplot::plot_grid(plotlist=plots, ncol=3))
dev.off()

## ============================================================
## 7. Save results
## ============================================================

write.csv(cross_cor_df,    here("figures/cross_species_correlations.csv"),
          row.names=FALSE)
write.csv(cross_cor_DEG,   here("figures/cross_species_DEG_correlations.csv"),
          row.names=FALSE)
write.csv(chy_MET_cor,     here("figures/chytrid_DEG_MET_correlation.csv"),
          row.names=FALSE)
write.csv(cyano_MET_cor,   here("figures/cyano_DEG_MET_correlation.csv"),
          row.names=FALSE)