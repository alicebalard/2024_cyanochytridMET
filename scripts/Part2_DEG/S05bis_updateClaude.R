## S04_fullAnalysis.R
## Updated: April 2025
## Changes from previous version:
##   FIX: isoform collapsing now uses highest-expressed isoform (slice_max) instead of summing
##        Rationale: summing inflates counts and can distort DESeq2 variance estimates;
##        keeping the dominant isoform is more conservative and biologically interpretable
##   FIX: removed stray commented-out prose text at top of file
##   FIX: deduplicated eukTransc definition (was defined twice)
##   FIX: write.csv output uses .csv extension consistently (was .tsv with write.csv)
##   FIX: dev.new() calls before pdf() are unnecessary and removed (dev.new opens a screen
##        device that gets ignored once pdf() is called; it just wastes a window)
##   NOTE: whichOrg string matching is whitespace-sensitive — trimws() helps but the
##         multi-space string "in_chytrid_alone  in_cyano_alone" (2 spaces) in the
##         contamination filter must exactly match what paste() produces. Consider
##         using a more robust approach (see comment in contamination section).

source("libLoad.R")
source("dataLoad.R")

## ============================================================
## Load count matrix
## ============================================================
RSEM.gene <- read.csv(
  "../../data/run_DESEQ2_Erika/RSEM.gene.counts.matrix",
  sep = "\t")
rownames(RSEM.gene) <- RSEM.gene$X
## 9766 genes

## ============================================================
## Classify genes by which organism/condition they appear in
## ============================================================
## NB: column name patterns "chy", "both", "cyano" must match your sample names exactly
a <- ifelse(rowSums(RSEM.gene[, grep("chy",   names(RSEM.gene))]) != 0, "in_chytrid_alone", "")
b <- ifelse(rowSums(RSEM.gene[, grep("both",  names(RSEM.gene))]) != 0, "in_both_organisms", "")
c <- ifelse(rowSums(RSEM.gene[, grep("cyano", names(RSEM.gene))]) != 0, "in_cyano_alone", "")

RSEM.gene$whichOrg <- trimws(gsub("  +", " ", paste(a, b, c, sep = " ")))
## FIX: gsub("  +", " ", ...) collapses any multi-space gaps before trimws,
##      making string matching in the contamination filter robust regardless of
##      which of the three conditions is empty

RSEM.gene$whichTranscriptome <- ifelse(
  grepl("TRINITY", RSEM.gene$X), "chytrid", "cyano")

table(RSEM.gene$whichTranscriptome, RSEM.gene$whichOrg)

## ============================================================
## 1. Identify and export chytrid assembly contamination
##    (chytrid transcripts expressed when only cyano is present)
## ============================================================
contam_patterns <- c(
  "in_both_organisms in_cyano_alone",
  "in_chytrid_alone in_cyano_alone",
  "in_chytrid_alone in_both_organisms in_cyano_alone",
  "in_cyano_alone"
)

listOfTranscriptContaminant <- RSEM.gene[
  RSEM.gene$whichTranscriptome %in% "chytrid" &
    RSEM.gene$whichOrg %in% contam_patterns, "X"]

write.csv(listOfTranscriptContaminant,
          "../../data/listOfTranscriptContaminant_toRmFromChytridTranscriptome.csv",
          quote = FALSE, row.names = FALSE)

## ============================================================
## 2. Select expressed genes per organism (excluding contamination)
## ============================================================

## --- Chytrid ---
sequencedChytridGenes <- RSEM.gene[
  RSEM.gene$whichOrg %in% c("in_chytrid_alone in_both_organisms",
                             "in_chytrid_alone",
                             "in_both_organisms") &
    RSEM.gene$whichTranscriptome %in% "chytrid", "X"]
length(sequencedChytridGenes) # 3747

RSEM.gene_chytrid <- RSEM.gene[RSEM.gene$X %in% sequencedChytridGenes, ]

## Keep only chytrid/both samples (drop cyano columns)
RSEM.gene_chytrid <- RSEM.gene_chytrid[, grep("cyano", names(RSEM.gene_chytrid), invert = TRUE)]

## Drop metadata columns
RSEM.gene_chytrid <- RSEM.gene_chytrid[
  , !names(RSEM.gene_chytrid) %in% c("X", "whichOrg", "whichTranscriptome")]

## Rename rows based on annotation
rownames(RSEM.gene_chytrid) <- make.unique(
  annotationChytrid$gene_name[match(rownames(RSEM.gene_chytrid), annotationChytrid$custom_gene_name)])

nrow(RSEM.gene_chytrid) # 3747

## FIX: Collapse isoforms by keeping the highest-expressed isoform (was: summing)
## Rationale: summing inflates counts and can distort DESeq2 dispersion estimates.
## slice_max on mean expression retains the most biologically relevant isoform.
RSEM.gene_chytrid <- RSEM.gene_chytrid %>%
  mutate(base_name = sub("\\.\\d+$", "", rownames(RSEM.gene_chytrid)),
         mean_expr  = rowMeans(.)) %>%
  group_by(base_name) %>%
  slice_max(mean_expr, n = 1, with_ties = FALSE) %>%
  select(-mean_expr) %>%
  tibble::column_to_rownames("base_name") %>%
  data.frame()

nrow(RSEM.gene_chytrid) # expect ~3156 or similar

## --- Cyanobacteria ---
sequencedCyanoGenes <- RSEM.gene[
  RSEM.gene$whichOrg %in% c("in_both_organisms in_cyano_alone",
                             "in_cyano_alone",
                             "in_both_organisms") &
    RSEM.gene$whichTranscriptome %in% "cyano", "X"]
length(sequencedCyanoGenes) # 3636

RSEM.gene_cyano <- RSEM.gene[RSEM.gene$X %in% sequencedCyanoGenes, ]

## Keep only cyano/both samples (drop chytrid columns)
RSEM.gene_cyano <- RSEM.gene_cyano[, grep("chy", names(RSEM.gene_cyano), invert = TRUE)]

## Drop metadata columns
RSEM.gene_cyano <- RSEM.gene_cyano[
  , !names(RSEM.gene_cyano) %in% c("X", "whichOrg", "whichTranscriptome")]

## Rename rows
rownames(RSEM.gene_cyano) <- make.unique(
  annotationCyano$gene_name[match(rownames(RSEM.gene_cyano), annotationCyano$custom_gene_name)])

nrow(RSEM.gene_cyano) # 3636

## Collapse isoforms — same approach as chytrid
RSEM.gene_cyano <- RSEM.gene_cyano %>%
  mutate(base_name = sub("\\.\\d+$", "", rownames(RSEM.gene_cyano)),
         mean_expr  = rowMeans(.)) %>%
  group_by(base_name) %>%
  slice_max(mean_expr, n = 1, with_ties = FALSE) %>%
  select(-mean_expr) %>%
  tibble::column_to_rownames("base_name") %>%
  data.frame()

nrow(RSEM.gene_cyano) # expect ~3589 or similar

## ============================================================
## 3. Sequencing depth / saturation QC (NOISeq)
## ============================================================
library(NOISeq)

makeSatPlot <- function(counts_mat, label, sample_indices) {
  eset   <- ExpressionSet(assayData = as.matrix(counts_mat))
  satdat <- dat(eset, k = 0, ndepth = 7, type = "saturation")
  explo.plot(satdat, toplot = 1, samples = sample_indices)
  title(main = label)
}

## Pre-filter saturation
makeSatPlot(RSEM.gene_chytrid, "Chytrid - alone (pre-filter)",  1:11)
makeSatPlot(RSEM.gene_chytrid, "Chytrid - both (pre-filter)",   12:19)
makeSatPlot(RSEM.gene_cyano,   "Cyano - both (pre-filter)",     1:8)
makeSatPlot(RSEM.gene_cyano,   "Cyano - alone (pre-filter)",    9:20)

## Filter: remove genes with zeros in >=3 samples per group
RSEM.gene_chytrid <- filterRSEMno3Nullpergp(RSEM.gene_chytrid)
RSEM.gene_cyano   <- filterRSEMno3Nullpergp(RSEM.gene_cyano)

## Post-filter saturation
makeSatPlot(RSEM.gene_chytrid, "Chytrid - all (post-filter)", 1:ncol(RSEM.gene_chytrid))
makeSatPlot(RSEM.gene_cyano,   "Cyano - all (post-filter)",   1:ncol(RSEM.gene_cyano))

## Remove low-saturation sample (58.9%)
RSEM.gene_cyano <- RSEM.gene_cyano[, !names(RSEM.gene_cyano) %in% "met_both_In11"]

nrow(RSEM.gene_chytrid) # 835
nrow(RSEM.gene_cyano)   # 555

## ============================================================
## 4. DESeq2 differential expression
## ============================================================
contrast_chytrid <- calculateContrasts(
  my_countsmatrix = RSEM.gene_chytrid,
  my_org = "chy")

contrast_cyano <- calculateContrasts(
  my_countsmatrix = RSEM.gene_cyano,
  my_org = "cyano")

## ============================================================
## 5. Volcano plots
## ============================================================
W <- 10; H <- 8

## Chytrid
V_chy_inf_ctrl <- makeVolcano(contrast_chytrid$resr_inf_effect_control,
  "no MET: zoospores vs infecting chytrids",        "logos/logo1.png", "left")
V_chy_inf_met  <- makeVolcano(contrast_chytrid$resr_inf_effect_met,
  "MET: zoospores vs infecting chytrids",            "logos/logo2.png", "left")
V_chy_met_1org <- makeVolcano(contrast_chytrid$resr_met_effect_1org,
  "free-living zoospores: no MET vs MET",            "logos/logo3.png", "left")
V_chy_met_2org <- makeVolcano(contrast_chytrid$resr_met_effect_2orgs,
  "infecting chytrids: no MET vs MET",               "logos/logo4.png", "left")

## Cyano
V_cya_inf_ctrl <- makeVolcano(contrast_cyano$resr_inf_effect_control,
  "no MET: uninfected vs infected cyanobacteria",    "logos/logo5.png", "left")
V_cya_inf_met  <- makeVolcano(contrast_cyano$resr_inf_effect_met,
  "MET: uninfected vs infected cyanobacteria",       "logos/logo6.png", "left")
V_cya_met_1org <- makeVolcano(contrast_cyano$resr_met_effect_1org,
  "uninfected cyanobacteria: no MET vs MET",         "logos/logo7.png", "right")
V_cya_met_2org <- makeVolcano(contrast_cyano$resr_met_effect_2orgs,
  "infected cyanobacteria: no MET vs MET",           "logos/logo8.png", "left")

## FIX: removed dev.new() calls before pdf() — they open unused screen devices
pdf("../../figures/Fig2_chytrid_volc.pdf", width = W, height = H)
cowplot::plot_grid(V_chy_inf_ctrl$plot, V_chy_inf_met$plot,
                   V_chy_met_1org$plot, V_chy_met_2org$plot,
                   labels = c("a","b","c","d"), label_size = 20)
dev.off()

pdf("../../figures/Fig3_cyano_volc.pdf", width = W, height = H)
cowplot::plot_grid(V_cya_inf_ctrl$plot, V_cya_inf_met$plot,
                   V_cya_met_1org$plot, V_cya_met_2org$plot,
                   labels = c("a","b","c","d"), label_size = 20)
dev.off()

## ============================================================
## 6. Venn diagrams
## ============================================================
library(ggvenn)

getGenes     <- function(x) rownames(x[!is.na(x$padj) & x$padj < 0.05, ])
selectDEGenes <- function(x) x[!is.na(x$padj) & x$padj < 0.05, ]

ggvenn(list(
  "infection effect - control" = getGenes(contrast_chytrid$resr_inf_effect_control),
  "MET effect - free living"   = getGenes(contrast_chytrid$resr_met_effect_1org),
  "MET effect - infecting"     = getGenes(contrast_chytrid$resr_met_effect_2orgs),
  "infection effect - MET"     = getGenes(contrast_chytrid$resr_inf_effect_met)),
  show_percentage = FALSE, fill_color = rep("white", 4),
  stroke_size = 0.5, set_name_size = 4)

ggvenn(list(
  "infection effect - control" = getGenes(contrast_cyano$resr_inf_effect_control),
  "MET effect - alone"         = getGenes(contrast_cyano$resr_met_effect_1org),
  "MET effect - infected"      = getGenes(contrast_cyano$resr_met_effect_2orgs),
  "infection effect - MET"     = getGenes(contrast_cyano$resr_inf_effect_met)),
  show_percentage = FALSE, fill_color = rep("white", 4),
  stroke_size = 0.5, set_name_size = 4)

## ============================================================
## 7. Export DEG table
## ============================================================

## Helper: extract DEG data frame from a named list of DESeq2 results
buildDEGTable <- function(contrast_list, label_map) {
  deg_list <- lapply(names(contrast_list)[1:4], function(nm) {
    res <- selectDEGenes(contrast_list[[nm]])
    if (nrow(res) == 0) return(data.frame())
    data.frame(
      geneName       = rownames(res),
      padj           = res$padj,
      log2FoldChange = res$log2FoldChange,
      comparison     = label_map[[nm]],
      stringsAsFactors = FALSE)
  })
  do.call(rbind, deg_list)
}

chytrid_labels <- list(
  resr_inf_effect_control = "Infection effect on chytrid gene expression, in the absence of MET",
  resr_inf_effect_met     = "Infection effect on chytrid gene expression, in the presence of MET",
  resr_met_effect_1org    = "MET effect on chytrid gene expression, in free-living chytrid zoospores",
  resr_met_effect_2orgs   = "MET effect on chytrid gene expression, in infecting chytrid")

cyano_labels <- list(
  resr_inf_effect_control = "Infection effect on cyanobacteria gene expression, in the absence of MET",
  resr_inf_effect_met     = "Infection effect on cyanobacteria gene expression, in the presence of MET",
  resr_met_effect_1org    = "MET effect on cyanobacteria gene expression, in uninfected cyanobacteria",
  resr_met_effect_2orgs   = "MET effect on cyanobacteria gene expression, in infected cyanobacteria")

fullDEGTable <- rbind(
  buildDEGTable(contrast_chytrid, chytrid_labels),
  buildDEGTable(contrast_cyano,   cyano_labels))

fullDEGTable$padj           <- signif(fullDEGTable$padj, 2)
fullDEGTable$log2FoldChange <- signif(fullDEGTable$log2FoldChange, 2)

table(fullDEGTable$comparison)

## FIX: use write.csv with .csv extension (was write.csv to a .tsv path — inconsistent)
write.csv(fullDEGTable, "../../figures/TableS1_fullDEGTable.csv", row.names = FALSE)

## ============================================================
## 8. GO enrichment
## ============================================================
universe_chytrid <- rownames(RSEM.gene_chytrid)
universe_cyano   <- rownames(RSEM.gene_cyano)

## Chytrid
getGOBubbleZ(universe_chytrid, annotationChytrid,
             getGenes(contrast_chytrid$resr_inf_effect_control), GO_chytrid, isbubble = FALSE)
# axoneme (CC), padj = 0.0039

getGOBubbleZ(universe_chytrid, annotationChytrid,
             getGenes(contrast_chytrid$resr_met_effect_1org), GO_chytrid, isbubble = FALSE)
# no significant GO terms

getGOBubbleZ(universe_chytrid, annotationChytrid,
             getGenes(contrast_chytrid$resr_met_effect_2orgs), GO_chytrid, isbubble = FALSE)
# axoneme (CC), padj = 0.0077

## Cyano
getGOBubbleZ(universe_cyano, annotationCyano,
             getGenes(contrast_cyano$resr_met_effect_1org[
               contrast_cyano$resr_met_effect_1org$log2FoldChange > 0, ]),
             GO_cyano, isbubble = FALSE)

getGOBubbleZ(universe_cyano, annotationCyano,
             getGenes(contrast_cyano$resr_met_effect_1org),
             GO_cyano, isbubble = FALSE)
# no significant GO terms
