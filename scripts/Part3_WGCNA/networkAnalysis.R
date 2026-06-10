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
# trait          module    cor      p  padj n_genes
# 10   MET  darkolivegreen -0.683 0.0144 0.325      49
# 40   MET  palevioletred3  0.671 0.0170 0.325      35
# 11   MET     floralwhite -0.665 0.0183 0.325      43
# 39   MET     lightyellow  0.658 0.0200 0.325      59
# 8    MET       honeydew1 -0.578 0.0492 0.625      32
# 55   MET          grey60  0.561 0.0577 0.625      61
# 9    MET      orangered4 -0.534 0.0739 0.650      47
# 48   MET       steelblue  0.525 0.0799 0.650      52
# 14   MET          brown4 -0.465 0.1270 0.913      42
# 24   MET             red -0.433 0.1600 0.913      80
# 57   MET             tan  0.402 0.1960 0.913      70
# 28   MET          orange  0.401 0.1970 0.913      56
# 56   MET        skyblue3  0.388 0.2130 0.913      48
# 46   MET    mediumorchid  0.371 0.2360 0.913      30
# 15   MET       darkgreen -0.365 0.2430 0.913      57
# 26   MET           black -0.363 0.2460 0.913      79
# 41   MET     saddlebrown  0.358 0.2530 0.913      52
# 20   MET    navajowhite2 -0.358 0.2540 0.913      34
# 54   MET   darkseagreen4  0.325 0.3030 0.913      32
# 45   MET        thistle2  0.311 0.3250 0.913      39
# 59   MET         darkred -0.307 0.3320 0.913      57
# 35   MET            pink  0.305 0.3340 0.913      77
# 64   MET           white  0.293 0.3550 0.913      55
# 5    MET        thistle1 -0.291 0.3590 0.913      36
# 53   MET      lightcyan1 -0.289 0.3620 0.913      43
# 61   MET       turquoise  0.259 0.4150 0.913     242
# 38   MET   darkslateblue  0.257 0.4190 0.913      41
# 1    MET    midnightblue  0.257 0.4210 0.913      64
# 65   MET            grey  0.253 0.4280 0.913      69
# 44   MET     yellowgreen  0.242 0.4480 0.913      49
# 12   MET           green -0.240 0.4510 0.913      83
# 50   MET  lavenderblush3 -0.237 0.4570 0.913      32
# 27   MET          salmon -0.218 0.4960 0.913      69
# 29   MET   antiquewhite4 -0.218 0.4960 0.913      31
# 37   MET        skyblue2  0.213 0.5070 0.913      30
# 31   MET        darkgrey -0.206 0.5210 0.913      56
# 21   MET           brown -0.196 0.5410 0.913     114
# 47   MET          coral2  0.196 0.5410 0.913      30
# 6    MET            blue -0.187 0.5600 0.913     177
# 13   MET         skyblue -0.162 0.6140 0.913      54
# 36   MET lightsteelblue1 -0.153 0.6360 0.913      44
# 49   MET          coral1  0.149 0.6440 0.913      31
# 18   MET     greenyellow -0.148 0.6460 0.913      73
# 30   MET     darkmagenta  0.144 0.6540 0.913      49
# 34   MET      lightpink4 -0.137 0.6710 0.913      33
# 33   MET          maroon  0.132 0.6830 0.913      33
# 25   MET           plum1 -0.127 0.6950 0.913      47
# 19   MET            cyan -0.121 0.7080 0.913      64
# 16   MET   darkturquoise  0.120 0.7100 0.913      57
# 22   MET          yellow -0.119 0.7120 0.913      84
# 62   MET          purple  0.113 0.7260 0.913      75
# 7    MET         salmon4 -0.111 0.7300 0.913      35
# 17   MET         magenta -0.094 0.7710 0.916      76
# 60   MET         bisque4 -0.092 0.7770 0.916      41
# 4    MET   mediumpurple3 -0.091 0.7800 0.916      46
# 3    MET     darkorange2  0.080 0.8050 0.916      43
# 42   MET      lightgreen  0.078 0.8110 0.916      61
# 52   MET           ivory -0.075 0.8170 0.916      43
# 58   MET          violet -0.062 0.8470 0.934      50
# 43   MET   paleturquoise -0.045 0.8900 0.945      51
# 23   MET           plum2 -0.045 0.8910 0.945      40
# 32   MET      darkorange -0.040 0.9010 0.945      56
# 51   MET       lightcyan -0.024 0.9400 0.970      61
# 2    MET       royalblue -0.004 0.9910 0.997      58
# 63   MET         sienna3  0.001 0.9970 0.997      49

## =====================================================================
## CHYTRID — alone + co-culture. Test MET *and* Infection, with the
## composition fraction included so artifacts are visible side-by-side.
## =====================================================================
chy_samples <- setdiff(colnames(counts_chy), remove_chy)
traits_chy  <- data.frame(
  row.names = chy_samples,
  MET       = as.numeric(grepl("met",  chy_samples)),
  Infection = as.numeric(grepl("both", chy_samples)),
  ChyFrac   = depth_chy[chy_samples] /
    (depth_chy[chy_samples] + depth_cyano[chy_samples]))

wgcna_chytrid <- runWGCNA(counts_chy, traits_chy,
                          organism = "chytrid",
                          label_fun = function(g) sub("_.*", "", g))
# chytrid: 23 samples x 1317 genes
# chytrid: power = 12 (R2 = 0.734, median k = 3.9, max R2 reached = 0.764)

## ---- infection vs composition: is the signal real or an artifact? ----
wide <- reshape(wgcna_chytrid$me_trait[, c("module", "trait", "cor")],
                idvar = "module", timevar = "trait", direction = "wide")
names(wide) <- sub("^cor\\.", "", names(wide))
wide$flag <- ifelse(abs(wide$ChyFrac) > 0.5, "composition-driven", "clean")
wide[order(-abs(wide$Infection)), ]      # modules tracking infection, flagged
# module ChyFrac Infection   MET               flag
# 9      grey   0.830    -0.289 0.173 composition-driven
# 8 turquoise   0.090     0.163 0.356              clean
# 7      blue   0.908     0.008 0.297 composition-driven

## strict within-co-culture artifact check (infection is constant here) -
both   <- grep("both", rownames(wgcna_chytrid$datExpr), value = TRUE)
r_both <- depth_chy[both] / depth_cyano[both]
round(cor(wgcna_chytrid$MEs[both, ], log2(r_both), use = "p"), 2)
# [,1]
# MEblue      0.79
# MEturquoise 0.11
# MEgrey      0.64

## ---- figures -------------------------------------------------------
plotDendroAndColors(wgcna_cyano$net$dendrograms[[1]],
                    wgcna_cyano$moduleColors[wgcna_cyano$net$blockGenes[[1]]],
                    "Modules", dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05, main = "Cyano modules")

plotDendroAndColors(wgcna_chytrid$net$dendrograms[[1]],
                    wgcna_chytrid$moduleColors[wgcna_chytrid$net$blockGenes[[1]]],
                    "Modules", dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05, main = "Chytrid modules")

## Panel A — module × trait correlation heatmap
mt <- wgcna_chytrid$me_trait
mt$trait  <- factor(mt$trait, levels = c("MET", "Infection", "ChyFrac"))
ord       <- mt %>% filter(trait == "ChyFrac") %>% arrange(cor) %>% pull(module)
mt$module <- factor(mt$module, levels = ord)

pA <- ggplot(mt, aes(trait, module, fill = cor)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f\n(p=%.2g)", cor, p)), size = 3) +
  scale_fill_gradient2(limits = c(-1, 1), low = "#2166AC",
                       mid = "white", high = "#B2182B") +
  labs(x = NULL, y = "Module eigengene", fill = "r",
       title = "Chytrid modules track library composition, not infection") +
  theme_minimal(base_size = 12)

## Panel B — within co-culture: eigengenes vs biomass ratio
both  <- grep("both", rownames(wgcna_chytrid$datExpr), value = TRUE)
ratio <- log2(depth_chy[both] / depth_cyano[both])
pB_df <- data.frame(ratio,
                    blue      = wgcna_chytrid$MEs[both, "MEblue"],
                    turquoise = wgcna_chytrid$MEs[both, "MEturquoise"]) |>
  pivot_longer(-ratio, names_to = "module", values_to = "ME")

pB <- ggplot(pB_df, aes(ratio, ME, color = module)) +
  geom_point(size = 2.5) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_color_manual(values = c(blue = "#3B5BA5", turquoise = "#1B9E77")) +
  labs(x = "log2(chytrid : cyano read ratio)", y = "Module eigengene",
       title = "Correlation of module eigengenes \n (overall expression profile) with composition") +
  theme_bw(base_size = 12) +
  theme(legend.position = "inside", legend.position.inside =  c(.8, .8))

pB

met <- subset(wgcna_cyano$me_trait, trait == "MET")
n   <- nrow(wgcna_cyano$datExpr)        # 12 free-living samples
m   <- nrow(met)                        # number of modules (BH m)

## |r| giving nominal p = 0.05, and the FDR "floor" = .05/m for the single
## best module (if the top can't clear it, nothing can)
r_at_p  <- function(pv) { t <- qt(1 - pv/2, n - 2); t / sqrt((n - 2) + t^2) }
r_nom   <- r_at_p(0.05)
r_fdr   <- r_at_p(0.05 / m)

met <- met %>% mutate(module = reorder(module, cor),
                      nominal = p < 0.05)

pC <- ggplot(met, aes(cor, module)) +
  annotate("rect", xmin = -r_nom, xmax = r_nom, ymin = -Inf, ymax = Inf,
           fill = "grey93") +
  geom_segment(aes(x = 0, xend = cor, yend = module), color = "grey75") +
  geom_vline(xintercept = 0, color = "grey30", linewidth = .3) +
  geom_vline(xintercept = c(-r_nom, r_nom), linetype = 2, color = "grey45") +
  geom_vline(xintercept = c(-r_fdr, r_fdr), color = "firebrick") +
  geom_point(aes(fill = as.character(module), size = nominal),
             shape = 21, color = "grey20", stroke = .3) +
  scale_fill_identity() +
  scale_size_manual(values = c(`FALSE` = 1.8, `TRUE` = 4.5), guide = "none") +
  scale_x_continuous(limits = c(-1, 1)) +
  labs(x = "Correlation of module eigengene with MET (r)", y = NULL,
       title = "Cyanobacteria: no co-expression module correlates with metolachlor") +
  theme_minimal(base_size = 11) +
  theme(axis.text.y = element_text(size = 5.5), panel.grid.major.y = element_blank())

fig5 <- cowplot::plot_grid(
  cowplot::plot_grid(pA, pB, labels = c("a", "b"), 
                   rel_widths = c(1, .8), ncol = 2),
  pC, ncol = 1, labels = c("", "c"), rel_heights = c(.8,1))

ggsave(here("figures/Fig5_networks.pdf"), width = 11, height = 10)

## =====================================================================
## Are the DEGs concentrated in particular modules?
## =====================================================================
degInModules <- function(wg, degs) {
  mc   <- wg$moduleColors
  degs <- intersect(unique(as.character(degs)), names(mc))
  data.frame(gene = degs, module = mc[degs], row.names = NULL)
}

degEnrichment <- function(wg, degs) {                 # one-sided Fisher per module
  mc   <- wg$moduleColors
  degs <- intersect(unique(as.character(degs)), names(mc))
  res  <- do.call(rbind, lapply(setdiff(unique(mc), "grey"), function(m) {
    tab <- table(in_module = mc == m, is_deg = names(mc) %in% degs)
    data.frame(module = m, module_size = sum(mc == m),
               degs_in_module = sum(mc == m & names(mc) %in% degs),
               n_degs = length(degs),
               p = signif(fisher.test(tab, alternative = "greater")$p.value, 3))
  }))
  res$padj <- signif(p.adjust(res$p, "BH"), 3)
  res[order(res$p), ]
}

## chytrid DEGs (UniProt-style names)
chy_degs <- fullDEGTable$geneName[grepl("chytrid", fullDEGTable$comparison)]
degInModules(wgcna_chytrid, chy_degs)
# gene    module
# 1   AATM_MOUSE      grey
# 2  ABCB6_HUMAN      grey
# 3  ABCG2_MACMU      grey
# 4   ACA1_SCHPO      grey
# 5   AMYA_DROMA      grey
# 6   ATC2_CRYNH      grey
# 7  CBPC1_HUMAN      grey
# 8   FHIT_DICDI      grey
# 9   GDIR_YEAST      grey
# 10  RFC4_ARATH      grey
# 11  SIK3_HUMAN turquoise
# 12  TRNL_SCHPO      grey
# 13  YBA9_SCHPO      grey

degEnrichment(wgcna_chytrid, chy_degs)
# module module_size degs_in_module n_degs    p padj
# 1 turquoise         119              1     13 0.71    1
# 2      blue          81              0     13 1.00    1

## cyano DEGs (numeric gene_ids; only the infected-MET set, flagged low-confidence)
cyano_degs <- fullDEGTable$geneName[grepl("cyanobacteria", fullDEGTable$comparison)]
degInModules(wgcna_cyano, cyano_degs)
# gene       module
# 1  77286364     skyblue3
# 2  77286711    turquoise
# 3  77286717    turquoise
# 4  77286801   darkorange
# 5  77287091      magenta
# 6  77287231     thistle2
# 7  77287375          red
# 8  77287779          tan
# 9  77287962 midnightblue
# 10 77288018      salmon4
# 11 77288928       purple
# 12 77289032    turquoise
# 13 77289563         pink

sapply(degInModules(wgcna_cyano, cyano_degs)$gene, cyano_label)
# 77286364                       77286711                       77286717 
# "rimO"                         "secY" "DNA-directed RNA polymerase " 
# 77286801                       77287091                       77287231 
# "recombinase family protein" "RuBisCO accumulation factor " "F0F1 ATP synthase subunit ga" 
# 77287375                       77287779                       77287962 
# "gvpC" "glycosyltransferase family 2"                         "purB" 
# 77288018                       77288928                       77289032 
# "GNAT family N-acetyltransfer"                         "groL"       "pseudouridine synthase" 
# 77289563 
# "NADH-quinone oxidoreductase " 

degEnrichment(wgcna_cyano, cyano_degs)
# module module_size degs_in_module n_degs      p padj
# 22       turquoise         242              3     13 0.0483    1
# 62         salmon4          35              1     13 0.1160    1
# 24        thistle2          39              1     13 0.1290    1
# 43        skyblue3          48              1     13 0.1560    1
# 59      darkorange          56              1     13 0.1800    1
# 12    midnightblue          64              1     13 0.2030    1
# 23             tan          70              1     13 0.2200    1
# 55          purple          75              1     13 0.2340    1
# 21         magenta          76              1     13 0.2360    1
# 41            pink          77              1     13 0.2390    1
# 9              red          80              1     13 0.2470    1
# 1            white          55              0     13 1.0000    1
# 2   darkolivegreen          49              0     13 1.0000    1
# 3       lightcyan1          43              0     13 1.0000    1
# 4     navajowhite2          34              0     13 1.0000    1
# 5       lightpink4          33              0     13 1.0000    1
# 6     mediumorchid          30              0     13 1.0000    1
# 7    darkslateblue          41              0     13 1.0000    1
# 8         darkgrey          56              0     13 1.0000    1
# 10   darkseagreen4          32              0     13 1.0000    1
# 11           brown         114              0     13 1.0000    1
# 13      lightgreen          61              0     13 1.0000    1
# 14          orange          56              0     13 1.0000    1
# 15   darkturquoise          57              0     13 1.0000    1
# 16            cyan          64              0     13 1.0000    1
# 17     greenyellow          73              0     13 1.0000    1
# 18         darkred          57              0     13 1.0000    1
# 19  palevioletred3          35              0     13 1.0000    1
# 20           black          79              0     13 1.0000    1
# 25           plum2          40              0     13 1.0000    1
# 26      orangered4          47              0     13 1.0000    1
# 27         skyblue          54              0     13 1.0000    1
# 28   mediumpurple3          46              0     13 1.0000    1
# 29        thistle1          36              0     13 1.0000    1
# 30       lightcyan          61              0     13 1.0000    1
# 31     saddlebrown          52              0     13 1.0000    1
# 32          grey60          61              0     13 1.0000    1
# 33     lightyellow          59              0     13 1.0000    1
# 34           plum1          47              0     13 1.0000    1
# 35  lavenderblush3          32              0     13 1.0000    1
# 36            blue         177              0     13 1.0000    1
# 37 lightsteelblue1          44              0     13 1.0000    1
# 38     yellowgreen          49              0     13 1.0000    1
# 39         bisque4          41              0     13 1.0000    1
# 40           green          83              0     13 1.0000    1
# 42          brown4          42              0     13 1.0000    1
# 44     darkmagenta          49              0     13 1.0000    1
# 45          yellow          84              0     13 1.0000    1
# 46   paleturquoise          51              0     13 1.0000    1
# 47   antiquewhite4          31              0     13 1.0000    1
# 48        skyblue2          30              0     13 1.0000    1
# 49          violet          50              0     13 1.0000    1
# 50       darkgreen          57              0     13 1.0000    1
# 51          salmon          69              0     13 1.0000    1
# 52       royalblue          58              0     13 1.0000    1
# 53     darkorange2          43              0     13 1.0000    1
# 54     floralwhite          43              0     13 1.0000    1
# 56       steelblue          52              0     13 1.0000    1
# 57         sienna3          49              0     13 1.0000    1
# 58           ivory          43              0     13 1.0000    1
# 60          maroon          33              0     13 1.0000    1
# 61          coral1          31              0     13 1.0000    1
# 63          coral2          30              0     13 1.0000    1
# 64       honeydew1          32              0     13 1.0000    1