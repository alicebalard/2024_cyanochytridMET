## ================================
## Part 3: Co-expression analysis
## ================================
library(here)
source(here("scripts/Part2_DEG/S05_fullAnalysis.R"))
library(WGCNA); library(DESeq2)
options(stringsAsFactors = FALSE)
allowWGCNAThreads()

## =====================================================================
## Per-organism WGCNA  (traits = data.frame of numeric sample annotations)
## =====================================================================
runWGCNA <- function(counts, traits,
                     organism      = "organism",
                     networkType   = "signed",
                     minModuleSize = 30,
                     mergeCutHeight = 0.25,
                     deepSplit     = 2,
                     RsquaredCut   = 0.85,
                     label_fun     = identity) {
  
  samples <- rownames(traits)
  stopifnot(all(samples %in% colnames(counts)))
  
  ## 1. subset + VST (blind -> design irrelevant) ----------------------
  counts <- round(as.matrix(counts[, samples, drop = FALSE]))
  dds <- DESeqDataSetFromMatrix(
    counts, data.frame(row.names = samples), design = ~ 1)
  dds <- estimateSizeFactors(dds)
  datExpr <- t(assay(vst(dds, blind = TRUE)))
  
  gsg     <- goodSamplesGenes(datExpr, verbose = 0)
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
  traits  <- traits[rownames(datExpr), , drop = FALSE]
  message(organism, ": ", nrow(datExpr), " samples x ", ncol(datExpr), " genes")
  
  ## 2. soft power = lowest reaching scale-free R^2 >= RsquaredCut ------
  powers <- 1:20
  sft <- pickSoftThreshold(datExpr, powerVector = powers,
                           networkType = networkType, corFnc = "bicor",
                           corOptions = list(maxPOutliers = 0.1),
                           RsquaredCut = RsquaredCut, verbose = 0)
  fit <- -sign(sft$fitIndices$slope) * sft$fitIndices$SFT.R.sq
  pwr <- sft$powerEstimate
  if (is.na(pwr)) pwr <- powers[which(fit >= RsquaredCut)[1]]
  if (is.na(pwr)) pwr <- if (networkType == "signed") 12 else 6
  message(organism, ": power = ", pwr,
          " (R2 = ", round(fit[pwr], 3),
          ", median k = ", round(sft$fitIndices$median.k.[pwr], 1),
          ", max R2 reached = ", round(max(fit), 3), ")")
  
  ## 3. modules (bicor, signed, dynamic tree cut) ----------------------
  net <- blockwiseModules(datExpr, power = pwr,
                          networkType = networkType, TOMType = networkType,
                          corType = "bicor", maxPOutliers = 0.1,
                          deepSplit = deepSplit, minModuleSize = minModuleSize,
                          mergeCutHeight = mergeCutHeight,
                          numericLabels = TRUE, pamRespectsDendro = FALSE,
                          maxBlockSize = ncol(datExpr) + 100, saveTOMs = FALSE, verbose = 0)
  moduleColors <- labels2colors(net$colors)
  names(moduleColors) <- colnames(datExpr)
  
  ## 4. module eigengenes vs EVERY trait -------------------------------
  MEs   <- orderMEs(moduleEigengenes(datExpr, moduleColors)$eigengenes)
  meCor <- cor(MEs, traits, use = "p")
  meP   <- corPvalueStudent(meCor, nrow(datExpr))
  meAdj <- apply(meP, 2, p.adjust, method = "BH")        # BH within each trait
  me_trait <- do.call(rbind, lapply(colnames(traits), function(tr)
    data.frame(trait = tr,
               module  = sub("^ME", "", rownames(meCor)),
               cor     = round(meCor[, tr], 3),
               p       = signif(meP[, tr], 3),
               padj    = signif(meAdj[, tr], 3),
               n_genes = as.integer(table(moduleColors)[sub("^ME","",rownames(meCor))]),
               row.names = NULL)))
  me_trait <- me_trait[order(me_trait$trait, me_trait$p), ]
  
  ## 5. hub genes: top 10% kWithin per module --------------------------
  kIN <- intramodularConnectivity.fromExpr(datExpr, moduleColors,
                                           networkType = networkType, power = pwr, corFnc = "bicor")
  rownames(kIN) <- colnames(datExpr)
  hubs <- do.call(rbind, lapply(setdiff(unique(moduleColors), "grey"), function(m) {
    g  <- names(moduleColors)[moduleColors == m]
    kw <- kIN[g, "kWithin"]
    h  <- g[kw >= quantile(kw, 0.9, na.rm = TRUE)]
    data.frame(module = m, gene = h, label = vapply(h, label_fun, character(1)),
               kWithin = round(kIN[h, "kWithin"], 2), row.names = NULL)
  }))
  hubs <- hubs[order(hubs$module, -hubs$kWithin), ]
  
  list(power = pwr, sft = sft, net = net, moduleColors = moduleColors,
       MEs = MEs, meCor = meCor, meP = meP, me_trait = me_trait,
       hubs = hubs, traits = traits, datExpr = datExpr)
}

## per-sample organism library sizes (for the composition trait) -------
depth_chy   <- colSums(counts_chy)
depth_cyano <- colSums(counts_cyano)

## =====================================================================
## CYANO — free-living samples only.  Only MET varies (all are "alone").
## =====================================================================
cyano_alone <- setdiff(grep("cyano", colnames(counts_cyano), value = TRUE), remove_cyano)
traits_cyano <- data.frame(row.names = cyano_alone,
                           MET = as.numeric(grepl("^met", cyano_alone)))

wgcna_cyano <- runWGCNA(counts_cyano, traits_cyano,
                        organism = "cyano", label_fun = cyano_label)
# cyano: 12 samples x 3706 genes
# cyano: power = 15 (R2 = 0.852, median k = 11.4, max R2 reached = 0.901)

subset(wgcna_cyano$me_trait, trait == "MET")

## =====================================================================
## Are the DEGs concentrated in particular modules?
## =====================================================================
degInModules <- function(wg, degs) {
  mc   <- wg$moduleColors
  degs <- intersect(unique(as.character(degs)), names(mc))
  data.frame(gene = degs, module = mc[degs], row.names = NULL)
}

## cyano DEGs (numeric gene_ids; only the infected-MET set, flagged low-confidence)
cyano_degs <- fullDEGTable$geneName[grepl("cyanobacteria", fullDEGTable$comparison)]
degmodcyano <- degInModules(wgcna_cyano, cyano_degs)
df <- data.frame(geneName = sapply(degmodcyano$gene, cyano_label),
                 gene = degmodcyano$gene,
                 module = degmodcyano$module)%>% 
  arrange

## ---- panel a: module dendrogram --------------------------------------
## plotDendroAndColors() draws to the device and returns NULL, so it can't
## be assigned straight into a cowplot grid — wrap it as a ggplot/grob.
pA <- as.ggplot(function()
  plotDendroAndColors(wgcna_cyano$net$dendrograms[[1]],
                      wgcna_cyano$moduleColors[wgcna_cyano$net$blockGenes[[1]]],
                      "Modules", dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05, main = "Cyanobacteria alone modules"))

## ---- panel b: which module each cyano DEG falls in -------------------
df <- df %>%
  arrange(module, geneName) %>%
  mutate(geneName = factor(geneName, levels = rev(unique(geneName))))

pB <- ggplot(df, aes(x = 1, y = geneName)) +
  geom_point(aes(fill = module), shape = 21, size = 5, colour = "grey30") +
  geom_text(aes(x = 1.12, label = module), hjust = 0, size = 3) +
  scale_fill_identity() +                       # module names are valid R colours
  scale_x_continuous(limits = c(0.9, 2.4)) +
  labs(x = NULL, y = NULL,
       title = "Cyanobacterial DEGs grouped by modules") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        panel.grid  = element_blank(),
        plot.title  = element_text(size = 10, face = "bold"))

met <- subset(wgcna_cyano$me_trait, trait == "MET")
## No correlation with module

## =====================================================================
## CHYTRID — free-living samples only.  Only MET varies (all are "alone").
## =====================================================================
chy_alone <- setdiff(grep("chy", colnames(counts_chy), value = TRUE), remove_chy)
traits_chy <- data.frame(row.names = chy_alone,
                         MET = as.numeric(grepl("^met", chy_alone)))

wgcna_chy <- runWGCNA(counts_chy, traits_chy,
                      organism = "chy", label_fun = function(g) sub("_.*", "", g))
# chy: 11 samples x 1252 genes
# chy: power = 12 (R2 = 0.256, median k = 11.2, max R2 reached = 0.693)

subset(wgcna_chy$me_trait, trait == "MET")

## =====================================================================
## Are the DEGs concentrated in particular modules?
## =====================================================================
## chy DEGs (numeric gene_ids; only the infected-MET set, flagged low-confidence)
chy_degs <- fullDEGTable$geneName[grepl("chytrid", fullDEGTable$comparison)]
degmodchy <- degInModules(wgcna_chy, chy_degs)
df <- data.frame(geneName = degmodchy$gene,
                 module = degmodchy$module) %>% 
  arrange(module)

## ---- panel a: module dendrogram --------------------------------------
## plotDendroAndColors() draws to the device and returns NULL, so it can't
## be assigned straight into a cowplot grid — wrap it as a ggplot/grob.
pC <- as.ggplot(function()
  plotDendroAndColors(wgcna_chy$net$dendrograms[[1]],
                      wgcna_chy$moduleColors[wgcna_chy$net$blockGenes[[1]]],
                      "Modules", dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05, main = "Chytrid alone modules"))

## ---- panel b: which module each chy DEG falls in -------------------
df <- df %>%
  arrange(module, geneName) %>%
  mutate(geneName = factor(geneName, levels = rev(unique(geneName))))
pD <- ggplot(df, aes(x = 1, y = geneName)) +
  geom_point(aes(fill = module), shape = 21, size = 5, colour = "grey30") +
  geom_text(aes(x = 1.12, label = module), hjust = 0, size = 3) +
  scale_fill_identity() +                       # module names are valid R colours
  scale_x_continuous(limits = c(0.9, 2.4)) +
  labs(x = NULL, y = NULL,
       title = "Chytrid DEGs grouped by modules") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        panel.grid  = element_blank(),
        plot.title  = element_text(size = 10, face = "bold"))

## ----  module–MET correlations
met <- subset(wgcna_chy$me_trait, trait == "MET")
met # no correlation

## ---- assemble & save ------------------------------------------------
fig5 <- cowplot::plot_grid(
  cowplot::plot_grid(pA, pB, pC, pD, labels = c("a", "b", "c", "d"),
                     rel_widths = c(1, .4, 1, .4), ncol = 2))

ggsave(here("figures/Fig5_networks.pdf"), fig5, width = 15, height = 8)
