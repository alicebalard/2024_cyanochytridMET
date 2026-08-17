## Updated June 2026
library(here)
source(here("scripts/Part2_DEG/libLoad.R"))
source(here("scripts/Part2_DEG/dataLoad.R"))
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
source(here("scripts/Part2_DEG/functions.R"))

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

################################################################
## removing repeated elements and known spurious-hit families ##
################################################################
annot <- read.delim(here("gitignore/assemblyMergedFungi_filterEuk_simplified_GOKegg.tsv"),
                    sep = "\t", stringsAsFactors = FALSE)

# Safe long tokens — substring match on BOTH fields (low false-positive risk)
te_terms <- c(
  "transpos",                 # transposon, transposable, transposase, retrotransposon
  "retroelement", "retrovir", # retrovirus, retrovirus-related, retroviral
  "integrase", "reverse transcriptase",
  "mobile element", "helitron", "mariner", "harbinger",
  "piggybac", "\\bcopia\\b", "\\bgypsy\\b", "\\bTy[0-9]",
  "satellite", "low complexity",
  "intron[- ]encoded", "homing endonuclease", "maturase")   # <- mitochondrial intron ORFs

# Short TE phrasings — only with word boundaries
te_phrases <- c("\\bgag-?pol\\b", "\\bpol polyprotein\\b",
                "\\bgag polyprotein\\b", "\\bLINE-?1\\b", "\\bSINE\\b", "\\bLTR\\b")

desc_pattern <- paste(c(te_terms, te_phrases), collapse = "|")   # remove rna_terms

# 4) UniProt TE entry names — anchored, gene_name only (safe vs polymerases)
entry_pattern <- "^(POL|GAG|ENV|GAGPOL|LORF[0-9]|LINE1?|L1RE[0-9]|RTBV|RTJK|TF2|TY[0-9]|AI[0-9])_"

keep_override <- grepl("telomerase", annot$sprot_Top_BLASTX_hit, ignore.case = TRUE)

is_bad <-
  grepl(desc_pattern,  annot$gene_name,            ignore.case = TRUE, perl = TRUE) |
  grepl(desc_pattern,  annot$sprot_Top_BLASTX_hit, ignore.case = TRUE, perl = TRUE) |
  grepl(entry_pattern, annot$gene_name,            ignore.case = TRUE, perl = TRUE)

is_bad <- is_bad & !keep_override

# Stricter than your version: drop a gene if ANY of its annotation rows looks like a TE
bad_genes  <- unique(annot$gene_name[is_bad])
keep_genes <- setdiff(unique(annot$gene_name), bad_genes)
counts_chy <- counts_chy[rownames(counts_chy) %in% keep_genes, ]

# remove all-zero rows
counts_chy <- counts_chy[rowSums(counts_chy) > 0, ] 

# Remove genes detected in only 1 to 3 sample - likely artifacts
n_samples_detected_chy <- rowSums(counts_chy > 0)
message("Chytrid genes by n samples detected:")            
ggplot(data.frame(table(n_samples_detected_chy)), 
       aes(x = n_samples_detected_chy, y = Freq))+
  geom_bar(stat = "identity")+
  theme_bw()

counts_chy <- counts_chy[n_samples_detected_chy >= 4, ]
message("Chytrid genes after removing genes covered in 3 samples or less: ",  
        nrow(counts_chy))               
## 1317

## Check
removed <- sort(setdiff(unique(annot$gene_name), keep_genes))
head(removed, 50)                      # eyeball for false positives
"SYKC_SCHPO"  %in% keep_genes          # TRUE
"LORF2_HUMAN" %in% keep_genes          # FALSE
"AI7_USTMA"  %in% keep_genes   # FALSE

grep("telomerase", annot$sprot_Top_BLASTX_hit[!is_bad], ignore.case = TRUE, value = TRUE)[1] 
# should now return TERT, while 
"TC1A_CAEEL" %in% keep_genes # stays FALSE.

grep("rRNA|fibrillarin|methyltransferase",
     annot$sprot_Top_BLASTX_hit[!is_bad], ignore.case = TRUE, value = TRUE)[1:10]
# rRNA-modification enzymes should now SURVIVE
grep("polymerase|repeat", 
     annot$sprot_Top_BLASTX_hit[!is_bad], ignore.case = TRUE, value = TRUE)[1:10]
# ^ confirm polymerases / repeat-domain proteins SURVIVED

####################
## Same for cyano ##
####################
counts_cyano <- as.matrix(RSEM_cyano)

# Mobile / repetitive elements (substring match is safe for these)
te_terms <- c("transpos", "retroelement", "retrotranspos", "retron",
              "integrase", "reverse transcriptase", "group II intron",
              "mobile element", "insertion sequence", "\\bIS element\\b",
              "low complexity")

# Phage / viral structural — but NOT 'phage shock protein'
phage_terms <- c("prophage", "bacteriophage", "\\bphage\\b",
                 "capsid", "terminase", "portal protein",
                 "baseplate", "tail fiber", "tail protein")

# Structural / non-coding RNA *features* — NOT tRNA enzymes or ribosomal proteins
aa <- paste("Ala|Arg|Asn|Asp|Cys|Gln|Glu|Gly|His|Ile|Leu|Lys|Met",
            "Phe|Pro|Ser|Thr|Trp|Tyr|Val|Sec|Pyl|fMet|OTHER|Xxx|Undet", sep = "|")
# structural RNA: match the feature TYPE, not the product name
is_rna_feature <- grepl("\\[gbkey=(rRNA|tRNA|tmRNA|ncRNA|misc_RNA|antisense_RNA)\\]",
                        annotationCyano_final$info, ignore.case = TRUE)

bad_pattern_cyano <- paste(c(te_terms, phage_terms), collapse = "|")   # drop rna_terms
is_bad <- grepl(bad_pattern_cyano, annotationCyano_final$info,
                ignore.case = TRUE, perl = TRUE) | is_rna_feature

keep_override <- grepl("phage shock", annotationCyano_final$info, ignore.case = TRUE)
is_bad <- is_bad & !keep_override

# Drop-if-any-bad-row (consistent with the fungal filter)
bad_genes  <- unique(annotationCyano_final$gene_id[is_bad])
keep_genes <- setdiff(unique(annotationCyano_final$gene_id), bad_genes)
counts_cyano <- counts_cyano[rownames(counts_cyano) %in% keep_genes, ]

## manual check
removed2 <- unique(annotationCyano_final$gene_id[is_bad])
grep("methyltransferase|pseudouridine|ribosomal protein|cell envelope|pentapeptide",
     annotationCyano_final$protein[!is_bad], ignore.case = TRUE, value = TRUE) |> head()  # should SURVIVE
grep("\\[gbkey=rRNA\\]", annotationCyano_final$info[is_bad])[1:3]                          # rRNA features still removed

removed_cyano <- unique(annotationCyano_final[is_bad, c("gene_id", "gene", "protein")])
removed_cyano[order(removed_cyano$protein), ]        # eyeball THIS for false positives

grep("methyltransferase|pseudouridine|rRNA|ribosomal RNA",
     removed_cyano$protein, ignore.case = TRUE, value = TRUE)
# any hits here = real enzymes wrongly removed

grep("pentapeptide repeat|cell envelope|--tRNA ligase|ribosomal protein|phage shock",
     annotationCyano_final$info[!is_bad], ignore.case = TRUE, value = TRUE)[1:10]

# remove all-zero rows
counts_cyano <- counts_cyano[rowSums(counts_cyano) > 0, ]

n_samples_detected_cyano <- rowSums(counts_cyano > 0)    
message("Cyano genes by n samples detected:")             
ggplot(data.frame(table(n_samples_detected_cyano)), 
       aes(x = n_samples_detected_cyano, y = Freq))+
  geom_bar(stat = "identity")+
  theme_bw()

counts_cyano <- counts_cyano[n_samples_detected_cyano >= 4, ]
message("Cyano genes after removing genes covered in 3 samples or less: ",  
        nrow(counts_cyano))   
## 3706

# Saturation plot: number of genes detected at observed and simulated
# sequencing depths, and marginal gain in gene detection per million additional reads
eset_chy <- ExpressionSet(assayData = counts_chy)
mysaturation_chy <- dat(eset_chy, k = 0, ndepth = 7, type = "saturation")

eset_cyano <- ExpressionSet(assayData = counts_cyano)
mysaturation_cyano <- dat(eset_cyano, k = 0, ndepth = 7, type = "saturation")

## Plots:
pdf(here("figures/Fig_S2.saturation_plots.pdf"), width=13, height=9)

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

# met_both_In12 control_both_In3 control_both_In5     met_both_In9     met_both_In8     met_both_In7 
# 0.8445958        0.8831965        1.0774280        3.4430569        3.9483074       19.8092919 
# control_both_In2 control_both_In1    met_both_In10 control_both_In6 control_both_In4    met_both_In11 
# 77.0348627      118.8968005      157.5408108      177.1222394      344.7722327      429.9527636 

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

pdf(here("figures/Fig2.compositionPlot.pdf"), width=11, height=6)
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

# Based on saturation plots and ratio chytrid/cyano counts per sample, 
# we remove the following samples:
#   
# Chytrid matrix: 
remove_chy <- c("met_chy_Z9")
counts_chy <- counts_chy[,!colnames(counts_chy) %in% remove_chy]

# Cyano matrix: 
remove_cyano <- c("control_both_In1", "control_both_In4", "control_both_In6", "met_both_In10") 
counts_cyano <- counts_cyano[,!colnames(counts_cyano) %in% remove_cyano]

# NB: strong difference infected/non infected → are the cyano missing genes in infected random, or infection related?

#####################
## Part 1. Chytrid ##
#####################

## --- 1. Infection effect in chytrid, no MET ---
## --- 2. Infection effect in chytrid, with MET ---
## --- 3. MET effect in chytrid alone ---
## --- 4. MET effect in chytrid infecting cyanobacteria

## Subset count matrices — explicit patterns to avoid cross-contamination

# --- 1. Infection effect in chytrid, no MET ---
# control_chy (alone) vs control_both (infecting), no MET
counts1_infeffect_noMET_chy <- counts_chy[
  , grep("^control_chy|^control_both", colnames(counts_chy))]
colnames(counts1_infeffect_noMET_chy)
table(sub("_[^_]+$", "", colnames(counts1_infeffect_noMET_chy)))
# control_both control_chy
#            6           6

# --- 2. Infection effect in chytrid, with MET ---
# met_chy (alone) vs met_both (infecting), with MET
counts2_infeffect_MET_chy <- counts_chy[
  , grep("^met_chy|^met_both", colnames(counts_chy))]
colnames(counts2_infeffect_MET_chy)
table(sub("_[^_]+$", "", colnames(counts2_infeffect_MET_chy)))
# met_both met_chy
#        6       5

# --- 3. MET effect in chytrid alone ---
# control_chy vs met_chy, no infection
counts3_METeffect_alone_chy <- counts_chy[
  , grep("^control_chy|^met_chy", colnames(counts_chy))]
colnames(counts3_METeffect_alone_chy)
table(sub("_[^_]+$", "", colnames(counts3_METeffect_alone_chy)))
# control_chy met_chy
#           6       5

# --- 4. MET effect in chytrid infecting cyanobacteria ---
# control_both vs met_both, during infection
counts4_METeffect_infecting_chy <- counts_chy[
  , grep("^control_both|^met_both", colnames(counts_chy))]
colnames(counts4_METeffect_infecting_chy)
table(sub("_[^_]+$", "", colnames(counts4_METeffect_infecting_chy)))
# control_both met_both
#            6        6

myFilterByExpWithinComp <- function(subcount){
  ## EdgeR::filterByExpr to determine which genes have sufficiently large counts to be retained in a statistical analysis.
  group <- factor(sub("_[^_]+$", "", colnames(subcount)))
  y <- DGEList(counts = subcount, group = group)
  keep <- filterByExpr(y, group = group)
  return(subcount[keep, ])
}

counts1_infeffect_noMET_chy <- myFilterByExpWithinComp(counts1_infeffect_noMET_chy)
counts2_infeffect_MET_chy <- myFilterByExpWithinComp(counts2_infeffect_MET_chy)
counts3_METeffect_alone_chy <- myFilterByExpWithinComp(counts3_METeffect_alone_chy)
counts4_METeffect_infecting_chy <- myFilterByExpWithinComp(counts4_METeffect_infecting_chy)

# Create named list automatically from object names
listCounts <- mget(c("counts1_infeffect_noMET_chy",
                     "counts2_infeffect_MET_chy",
                     "counts3_METeffect_alone_chy",
                     "counts4_METeffect_infecting_chy"))

# Automated summary for all count matrices in the list
invisible(lapply(names(listCounts), function(name) {
  x <- listCounts[[name]]
  message(name, ":  ", ncol(x), " samples, ", nrow(x), " genes")
}))
## invisible() suppresses the NULL list output that lapply would otherwise print to console.

# counts1_infeffect_noMET_chy:  12 samples, 882 genes
# counts2_infeffect_MET_chy:  11 samples, 1050 genes
# counts3_METeffect_alone_chy:  11 samples, 1012 genes
# counts4_METeffect_infecting_chy:  12 samples, 930 genes

############
## DESeq2 ##
############

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
  
  ddsr <- DESeq(ddsr, fitType="mean")
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
# out of 882 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 5, 0.57%
# LFC < 0 (down)     : 2, 0.23%
# outliers [1]       : 80, 9.1%
# low counts [2]     : 0, 0%
# (mean count < 9)

testIfCrazyFoldChange(res_counts1_infeffect_noMET_chy) # OK
# baseMean log2FoldChange        padj
# GDIR_YEAST  149.19884       3.337832 0.007524027
# ABCB6_HUMAN 615.66917      -2.571616 0.007524027
# AMYA_DROMA   88.31743       2.515683 0.012311049
# RFC4_ARATH  138.09361       2.368097 0.029433416
# CBPC1_HUMAN  55.21466      -2.095739 0.029433416
# ACA1_SCHPO  197.41375       2.092877 0.030769149
# TRNL_SCHPO  358.89619       1.471060 0.029433416

## --- 2. Infection effect in chytrid, with MET ---
res_counts2_infeffect_MET_chy <- myDESeq2_withlfcShrink(
  listCounts$counts2_infeffect_MET_chy, a="met_chy", b="met_both")
# out of 1050 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 0, 0%
# LFC < 0 (down)     : 0, 0%
# outliers [1]       : 28, 2.7%
# low counts [2]     : 0, 0%
# (mean count < 6)

## --- 3. MET effect in chytrid alone ---
res_counts3_METeffect_alone_chy <- myDESeq2_withlfcShrink(
  listCounts$counts3_METeffect_alone_chy, a="control_chy", b="met_chy")

# out of 1012 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 2, 0.2%
# LFC < 0 (down)     : 1, 0.099%
# outliers [1]       : 27, 2.7%
# low counts [2]     : 0, 0%
# (mean count < 5)

testIfCrazyFoldChange(res_counts3_METeffect_alone_chy) # ok
# baseMean log2FoldChange         padj
# AATM_MOUSE   27.06906       3.297744 0.0000263773
# ABCG2_MACMU  44.40865       2.936434 0.0022230968
# CBPC1_HUMAN 103.36536      -2.056562 0.0336600434

## --- 4. MET effect in chytrid infecting cyanobacteria
res_counts4_METeffect_infecting_chy <- myDESeq2_withlfcShrink(
  listCounts$counts4_METeffect_infecting_chy, a="control_both", b="met_both")
# out of 930 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 3, 0.32%
# LFC < 0 (down)     : 1, 0.11%
# outliers [1]       : 96, 10%
# low counts [2]     : 0, 0%
# (mean count < 8)

testIfCrazyFoldChange(res_counts4_METeffect_infecting_chy)
# baseMean log2FoldChange        padj
# SIK3_HUMAN 281.17333       3.444385 0.004672581
# ATC2_CRYNH 291.21001      -2.833865 0.030398485
# YBA9_SCHPO  66.94740       2.669970 0.037918844
# FHIT_DICDI  23.54381       2.514341 0.024805725

contrast_chytridgenome <- mget(c(
  "res_counts1_infeffect_noMET_chy",  
  "res_counts2_infeffect_MET_chy",
  "res_counts3_METeffect_alone_chy",
  "res_counts4_METeffect_infecting_chy"
))

## Interpretation
lapply(contrast_chytridgenome, function(x){
  annotationChytrid_final[match(rownames(x[
    x$padj < 0.05 &
      !is.na(x$padj),]),
    annotationChytrid_final$gene_name),c("gene_name", "sprot_Top_BLASTX_hit")]
  
})

###################
## Volcano plots ##
###################

V_chytrid_inf_effect_control <- makeVolcano(
  res          = contrast_chytridgenome$res_counts1_infeffect_noMET_chy,
  title        = "no MET, chytrids zoospores vs chytrids infecting",
  positionLogoStart = -2, positionLogoStop = 2,
  mylogo       = here("scripts/Part2_DEG/logos/logo1.png"))

V_chytrid_inf_effect_met <- makeVolcano(
  res          = contrast_chytridgenome$res_counts2_infeffect_MET_chy,
  title        = "MET, chytrids zoospores vs chytrids infecting",
  positionLogoStart = -2, positionLogoStop = 2,
  mylogo       = here("scripts/Part2_DEG/logos/logo2.png"))

V_chytrid_met_effect_1org <- makeVolcano(
  res          = contrast_chytridgenome$res_counts3_METeffect_alone_chy,
  title        = "free-living zoospores, no MET vs MET",
  positionLogoStart = -2, positionLogoStop = 2,
  mylogo       = here("scripts/Part2_DEG/logos/logo3.png"))

V_chytrid_met_effect_2orgs <- makeVolcano(
  res          = contrast_chytridgenome$res_counts4_METeffect_infecting_chy,
  title        = "chytrids infecting, no MET vs MET",
  positionLogoStart = -2, positionLogoStop = 2,
  mylogo       = here("scripts/Part2_DEG/logos/logo4.png"))

####################################################################
## Raw count plots for DEGs - one plot per contrast, saved to PDF ##
####################################################################
rawplot <- plot_degs_raw(
  contrast_list       = contrast_chytridgenome,
  count_matrix        = counts_chy,
  contrast_conditions = list(
    res_counts1_infeffect_noMET_chy     = c("control_chy",  "control_both"),
    res_counts2_infeffect_MET_chy       = c("met_chy",      "met_both"),
    res_counts3_METeffect_alone_chy     = c("control_chy",  "met_chy"),
    res_counts4_METeffect_infecting_chy = c("control_both", "met_both")),
  condition_levels = c("control_chy", "control_both", "met_chy", "met_both"),
  colour_map = c(control_chy = "darkgreen", met_chy = "lightgreen",
                 control_both = "purple",  met_both = "plum"),
  ncol = 7)

pdf(here("figures/Fig3_DEGchytrid.pdf"), width=14, height=10)
cowplot::plot_grid(
  cowplot::plot_grid(V_chytrid_inf_effect_control$plot,
                     V_chytrid_inf_effect_met$plot,
                     V_chytrid_met_effect_1org$plot,
                     V_chytrid_met_effect_2orgs$plot,
                     labels=c("a","b","c","d"), ncol = 4,
                     label_size=20),
  rawplot, labels = c("","e"), label_size=20, rel_heights = c(.5,1), ncol = 1)
dev.off()

## Cyanobacteria

## 1/ MET effect in cyano alone 
# control_cyano vs met_cyano — both cyano-dominated → classic DESeq2 is valid here
counts1_METeffect_alone_cyano <- counts_cyano[
  , grep("^control_cyano|^met_cyano", colnames(counts_cyano))]
table(sub("_[^_]+$", "", colnames(counts1_METeffect_alone_cyano)))
# control_cyano     met_cyano 
# 6             6

# mirror the chytrid pipeline: re-filter for this subset 
counts1_METeffect_alone_cyano <- myFilterByExpWithinComp(counts1_METeffect_alone_cyano)

testEffectSize(counts1_METeffect_alone_cyano)   # size factors should be ~1

res_counts1_METeffect_alone_cyano <- myDESeq2_withlfcShrink(
  counts1_METeffect_alone_cyano, a = "control_cyano", b = "met_cyano")

## NO genes affected!!
# out of 3512 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 0, 0%
# LFC < 0 (down)     : 0, 0%
# outliers [1]       : 0, 0%
# low counts [2]     : 0, 0%
# (mean count < 9)

## 2/ MET effect in cyano infected 
# control_both vs met_both → classic DESeq2 is valid here
counts2_METeffect_both_cyano <- counts_cyano[
  , grep("^control_both|^met_both", colnames(counts_cyano))]
table(sub("_[^_]+$", "", colnames(counts2_METeffect_both_cyano)))
# control_both     met_both 
# 3            5 

# mirror the chytrid pipeline: re-filter for this subset 
# grp  <- factor(sub("_[^_]+$", "", colnames(counts2_METeffect_both_cyano)))
# keep <- filterByExpr(DGEList(counts = counts2_METeffect_both_cyano, group = grp), group = grp)
# counts2_METeffect_both_cyano <- counts2_METeffect_both_cyano[keep, ]
counts2_METeffect_both_cyano  <- myFilterByExpWithinComp(counts2_METeffect_both_cyano)

testEffectSize(counts2_METeffect_both_cyano)   # size factors x18 between extremes

res_counts2_METeffect_both_cyano <- myDESeq2_withlfcShrink(
  counts2_METeffect_both_cyano, a = "control_both", b = "met_both")

# out of 1943 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 10, 0.51%
# LFC < 0 (down)     : 3, 0.15%
# outliers [1]       : 102, 5.2%
# low counts [2]     : 113, 5.8%
# (mean count < 5)

testIfCrazyFoldChange(res_counts2_METeffect_both_cyano)
# baseMean log2FoldChange        padj
# 77287375  92.47731       5.454977 0.004956047
# 77286711  20.98320       4.856485 0.007723129
# 77287779 160.53516      -4.357565 0.009767132
# 77289563  18.98694      -4.065735 0.009767132
# 77287962  19.61255       3.713944 0.044944474
# 77287091  15.97261       3.652874 0.049479498
# 77287231  38.49329       3.625604 0.044944474
# 77288018  10.44370       3.521393 0.049479498
# 77286801  17.23683       3.432724 0.049479498
# 77288928  11.53882       3.428198 0.049479498

# Create named list automatically from object names
listCounts <- mget(c("counts1_METeffect_alone_cyano",
                     "counts2_METeffect_both_cyano"))

# Automated summary for all count matrices in the list
invisible(lapply(names(listCounts), function(name) {
  x <- listCounts[[name]]
  message(name, ":  ", ncol(x), " samples, ", nrow(x), " genes")
}))

contrast_cyanogenome <- mget(c(
  "res_counts1_METeffect_alone_cyano",  
  "res_counts2_METeffect_both_cyano"
))

## Interpretation
lapply(contrast_cyanogenome, function(x){
  annotationCyano_final[match(rownames(x[
    x$padj < 0.05 &
      !is.na(x$padj),]),
    annotationCyano_final$gene_id),c("gene_id", "gene", "protein")]
})

# $res_counts2_METeffect_both_cyano
# gene_id gene                                              protein
# 6016 77286364 rimO 30S ribosomal protein S12 methylthiotransferase RimO
# 2404 77286711 secY                  preprotein translocase subunit SecY
# 4455 77286717 <NA>            DNA-directed RNA polymerase subunit alpha
# 2113 77286801 <NA>                           recombinase family protein
# 4609 77287091 <NA>                        RuBisCO accumulation factor 1
# 1611 77287231 <NA>                      F0F1 ATP synthase subunit gamma
# 5834 77287375 gvpC                             gas vesicle protein GvpC
# 2348 77287779 <NA>                 glycosyltransferase family 2 protein
# 3622 77287962 purB                               adenylosuccinate lyase
# 6412 77288018 <NA>                      GNAT family N-acetyltransferase
# 2303 77288928 groL                                     chaperonin GroEL
# 2309 77289032 <NA>                               pseudouridine synthase
# 1420 77289563 <NA>                NADH-quinone oxidoreductase subunit J

###############################################
################ Table results ################ 
###############################################
## shared, organism-neutral labels so the columns line up
disp <- c("MET\nBoth organisms", "no MET\nBoth organisms",
          "MET\nAlone",          "no MET\nAlone")

M_chy <- build_comparison_matrix(
  res_list = contrast_chytridgenome,
  contrast_conditions = list(
    res_counts1_infeffect_noMET_chy     = c("control_chy",  "control_both"),
    res_counts2_infeffect_MET_chy       = c("met_chy",      "met_both"),
    res_counts3_METeffect_alone_chy     = c("control_chy",  "met_chy"),
    res_counts4_METeffect_infecting_chy = c("control_both", "met_both")),
  conditions = c("met_both", "control_both", "met_chy", "control_chy"),
  display    = disp, organism = "chytrid")

M_cyano <- build_comparison_matrix(
  res_list = contrast_cyanogenome,
  contrast_conditions = list(
    res_counts1_METeffect_alone_cyano = c("control_cyano", "met_cyano"),
    res_counts2_METeffect_both_cyano  = c("control_both",  "met_both")),
  conditions = c("met_both", "control_both", "met_cyano", "control_cyano"),
  display    = disp, organism = "cyano")

save_combined_matrix_csv(
  list("Chytrid transcripts"       = M_chy,
       "Cyanobacteria transcripts" = M_cyano),
  file = here("figures/Table2_comparison.csv"))       

###################
## Volcano plots ##
###################
# label DEGs by gene symbol, falling back to protein, then gene_id
cyano_label <- function(g) {
  i   <- match(g, annotationCyano_final$gene_id)
  lab <- annotationCyano_final$gene[i]
  if (is.na(lab) || lab == "") lab <- annotationCyano_final$protein[i]
  if (is.na(lab) || lab == "") lab <- g
  substr(lab, 1, 28)
}

# cyano — gene symbols instead of numeric gene_ids
V_cyano_met_effect_1org <- makeVolcano(
  res = contrast_cyanogenome$res_counts1_METeffect_alone_cyano,
  title = "free-living cyanobacteria, no MET vs MET",
  mylogo = here("scripts/Part2_DEG/logos/logo7.png"),
  label_fun = cyano_label)

V_cyano_met_effect_2org <- makeVolcano(
  res = contrast_cyanogenome$res_counts2_METeffect_both_cyano,
  title = "infected cyanobacteria, no MET vs MET",
  mylogo = here("scripts/Part2_DEG/logos/logo8.png"),
  label_fun = cyano_label, maxx = 5.6) 

####################################
## Cyano raw-count panel + figure ##
####################################

rawplot_cyano <- plot_degs_raw(
  contrast_list       = contrast_cyanogenome,
  count_matrix        = counts_cyano,
  contrast_conditions = list(
    res_counts1_METeffect_alone_cyano = c("control_cyano", "met_cyano"),
    res_counts2_METeffect_both_cyano  = c("control_both",  "met_both")),
  condition_levels    = c("control_cyano", "control_both", "met_cyano", "met_both"),
  colour_map          = c(control_cyano = "steelblue", met_cyano = "lightblue",
                          control_both  = "purple",    met_both  = "plum"),
  label_fun           = cyano_label,
  # label_samples       = c("control_both", "met_both"),   # <- only these get point labels
  ncol = 7)

pdf(here("figures/Fig4_DEGcyano.pdf"), width = 14, height = 10) 
cowplot::plot_grid(
  cowplot::plot_grid(V_cyano_met_effect_1org$plot,
                     V_cyano_met_effect_2org$plot,
                     labels = c("a", "b"), ncol = 2, label_size = 20),
  rawplot_cyano, labels = c("", "c"), label_size = 20, rel_heights = c(.5, 1),
  ncol = 1)
dev.off()

## ============================================================
## Infection-linked LOSS of cyanobacterial genes
## Do alone-expressed genes drop out of co-culture MORE than the
## reduced cyano biomass / sequencing depth alone predicts?
## ============================================================
counts        <- counts_cyano
alone_cols    <- grep("cyano", colnames(counts), value = TRUE)   # cyano alone
infected_cols <- grep("both",  colnames(counts), value = TRUE)   # co-culture
depth         <- colSums(counts)
k             <- 5                                               # detection threshold (reads)

## relative abundance of each gene in the cyano transcriptome when alone
a_g      <- rowMeans(sweep(counts[, alone_cols], 2, depth[alone_cols], "/"))
universe <- rowMeans(counts[, alone_cols] >= k) >= 0.8 & a_g > 0  # reliably ON when alone
message(sum(universe), " genes reliably expressed when alone")
# 3536 genes reliably expressed when alone

## biomass-only null: expected reads = relative abundance x library cyano depth
lambda   <- outer(a_g, depth[infected_cols])
p_detect <- 1 - pnbinom(k - 1, mu = lambda, size = 1)            # NB, overdispersed
obs_det  <- rowSums(counts[, infected_cols] >= k)

## P(detections <= observed | null), Poisson-binomial lower tail
pb_lower <- function(p, x) { d <- 1; for (pp in p) d <- c(d,0)*(1-pp) + c(0,d)*pp; sum(d[seq_len(x+1)]) }
p_loss <- rep(NA_real_, nrow(counts))
for (g in which(universe)) p_loss[g] <- pb_lower(p_detect[g, ], obs_det[g])
padj_loss <- rep(NA_real_, nrow(counts))
padj_loss[universe] <- p.adjust(p_loss[universe], method = "BH")

## compositional effect size (relative-abundance log2FC, infected vs alone)
b_g       <- rowMeans(sweep(counts[, infected_cols], 2, depth[infected_cols], "/"))
eps       <- min(a_g[a_g > 0]) / 10
l2fc_prop <- log2((b_g + eps) / (a_g + eps))

## candidates: lost more than depth predicts, relative abundance collapsed,
## detection genuinely expected, and absent in ~all co-culture libraries
res <- data.frame(gene = rownames(counts),
                  cpm_alone = round(a_g*1e6, 1), cpm_infect = round(b_g*1e6, 1),
                  exp_infect_det = round(rowSums(p_detect), 2), obs_infect_det = obs_det,
                  l2fc_prop = round(l2fc_prop, 2), padj_loss = padj_loss)[universe, ]

candidates <- subset(res, padj_loss < 0.05 & l2fc_prop < -1 &
                       exp_infect_det >= 0.7 * length(infected_cols) &
                       obs_infect_det <= 1)
candidates <- candidates[order(candidates$padj_loss), ]
candidates$label <- sapply(candidates$gene, cyano_label)

## sanity check: are candidates really ~0 even in the best-sequenced infected libraries?
deep_inf <- names(sort(depth[infected_cols], decreasing = TRUE))
counts[candidates$gene, deep_inf, drop = FALSE]

candidates

## Plot
n_inf <- length(infected_cols)   # number of co-culture libraries

plt <- res
plt$status <- with(plt, factor(
  ifelse(padj_loss < 0.05 & l2fc_prop < -1 &
           exp_infect_det >= 0.7 * n_inf & obs_infect_det <= 1, "candidate",
         ifelse(padj_loss < 0.05, "lost > biomass predicts (p<0.05)",
                "consistent with biomass")),
  levels = c("consistent with biomass",
             "lost > biomass predicts (p<0.05)", "candidate")))

ggplot(plt, aes(exp_infect_det, obs_infect_det)) +
  ## candidate zone: expected in >=70% of libraries, observed in <=1
  annotate("rect", xmin = 0.7 * n_inf, xmax = n_inf, ymin = -0.5, ymax = 1.5,
           fill = "firebrick", alpha = 0.07) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey40") +   # biomass-only expectation
  geom_jitter(aes(color = status, size = status, alpha = status),
              width = 0, height = 0.18) +
  ggrepel::geom_text_repel(
    data = subset(plt, status == "candidate"),
    aes(label = sapply(gene, cyano_label)),
    color = "firebrick", size = 3, min.segment.length = 0, box.padding = 0.6) +
  scale_color_manual(values = c("consistent with biomass"          = "grey70",
                                "lost > biomass predicts (p<0.05)"  = "#E69F00",
                                "candidate"                         = "firebrick")) +
  scale_size_manual(values  = c(0.6, 1.6, 3),  guide = "none") +
  scale_alpha_manual(values = c(0.25, 0.85, 1), guide = "none") +
  scale_x_continuous(breaks = 0:n_inf) +
  scale_y_continuous(breaks = 0:n_inf) +
  labs(x = "Expected co-culture libraries detecting the gene (biomass-only null)",
       y = "Observed co-culture libraries detecting the gene",
       color = NULL,
       title = "Gene loss in co-culture matches the biomass-only prediction",
       subtitle = sprintf("%d genes reliably expressed when alone — %d candidate after all criteria",
                          nrow(plt), sum(plt$status == "candidate"))) +
  theme_bw(base_size = 12) +
  theme(legend.position = "top")

ggsave(here("figures/FigS3_infection_loss_null.pdf"), width = 7.5, height = 6.5)

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
    "MET effect on cyanobacteria gene expression, in uninfected cyanobacteria",
    "MET effect on cyanobacteria gene expression, in infected cyanobacteria"))

contrast_cyanogenome_DEG <- contrast_cyanogenome_DEG[
  grep("Infection effect", contrast_cyanogenome_DEG$comparison, invert = T),]

fullDEGTable <- rbind(contrast_chytridgenome_DEG, contrast_cyanogenome_DEG) %>%
  mutate(padj           = signif(padj, 2),
         log2FoldChange = signif(log2FoldChange, 2)) %>%
  arrange(geneName)

table(fullDEGTable$comparison)

## Add annotation
fullDEGTable$geneFull <- annotationChytrid_final$sprot_Top_BLASTX_hit[
  match(fullDEGTable$geneName, annotationChytrid_final$gene_name)]

fullDEGTable$geneFull <- str_extract(fullDEGTable$geneFull, "(?<=Full=)[^;]+")

x <- annotationCyano_final$protein[
  match(fullDEGTable$geneName, annotationCyano_final$gene_id)]

fullDEGTable$geneFull[which(is.na(fullDEGTable$geneFull))] <- 
  x[which(is.na(fullDEGTable$geneFull))]
rm(x)

write.csv(fullDEGTable, here("figures/TableS1_fullDEGTable_annotated.csv"),
          row.names=FALSE)
# 
# #############################################
# ## GO enrichment                           ##
# #############################################
# 
# ## Chytrid
# getGOBubbleZ(universe  = rownames(counts1_infeffect_noMET_chy), 
#              annotation = annotationChytrid_final,
#              genelist  = getGenes(contrast_chytridgenome$res_counts1_infeffect_noMET_chy),
#              GO_df     = GO_chytrid, isbubble = FALSE)
# 
# getGOBubbleZ(universe  = rownames(counts2_infeffect_MET_chy), 
#              annotation = annotationChytrid_final,
#              genelist  = getGenes(contrast_chytridgenome$res_counts2_infeffect_MET_chy),
#              GO_df     = GO_chytrid, isbubble = FALSE)
# 
# getGOBubbleZ(universe  = rownames(counts3_METeffect_alone_chy), 
#              annotation = annotationChytrid_final,
#              genelist  = getGenes(contrast_chytridgenome$res_counts3_METeffect_alone_chy),
#              GO_df     = GO_chytrid, isbubble = FALSE)
# 
# getGOBubbleZ(universe  = rownames(counts4_METeffect_infecting_chy), 
#              annotation = annotationChytrid_final,
#              genelist  = getGenes(contrast_chytridgenome$res_counts4_METeffect_infecting_chy),
#              GO_df     = GO_chytrid, isbubble = FALSE)
# # no significant GO terms (not enough DEG)
