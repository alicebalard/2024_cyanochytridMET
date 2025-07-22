# Load necessary libraries

## Results of previous analysis
setwd("../Part2_DEG/")
source("fullAnalysis.R")

setwd("../Part3_WGCNA/")

## ============================================================
## 📦 1. Load Variance-stabilized matrices from previous script
## ============================================================
vst_cyano <- contrast_cyanogenome$vstr
vst_chy <- contrast_chytridgenome$vstr
vst_chy %>% nrow # 835 genes
vst_cyano %>% nrow # 555 genes

## For big data
options(stringsAsFactors = FALSE)
allowWGCNAThreads()

##############################################
## I. Single species co-expression analysis ##
##############################################

## ===============================
## ✅ 2.Pick soft power threshold
## ===============================
# The soft threshold (β) is a key parameter in WGCNA.
# It raises the correlation matrix to this power to emphasize strong correlations 
# and suppress weak ones, so your network approximates a scale-free topology — 
# the biological reality that a few genes are highly connected (hubs) while most aren’t.
# Scale-free topology fit index (SFT.R.sq): Good ≥ 0.8 (≥0.9 better); bad <0.6
# Mean connectivity between 5 and 100 is usually reasonable to not lose meaningful connections.

pickSoftPow <- function(vst){
  par(mfrow = c(2, 1))
  sth <- pickSoftThreshold(t(vst), powerVector = c(1:20), verbose = 5)
  plot(sth$fitIndices$Power, sth$fitIndices$SFT.R.sq, type="b", 
       xlab="Soft Threshold (power)", ylab="Scale Free Topology Model Fit, signed R^2",
       main="Scale-Free Topology Fit")
  abline(h=0.8, col="red", lty=2); abline(h=0.9, col="red", lty=1)
  plot(sth$fitIndices$Power, sth$fitIndices$mean.k., type="b", 
       xlab="Soft Threshold (power)", ylab="Mean connectivity",
       main="Scale-Free Topology Fit")
  abline(h=5, col="red", lty=2)
  par(mfrow = c(1, ))
}

pickSoftPow(vst_chy)
softPower_chy <- 8 # mean connectivity 1.5900 ok

pickSoftPow(vst_cyano)
softPower_chy <- 14 # mean connectivity 2.090 ok

## ===============================
## 🔗 3. Build modules
## ===============================
cor <- WGCNA::cor

net_chy <- blockwiseModules(
  datExpr = t(vst_chy),
  power = 8, 
  TOMType = "signed", minModuleSize = 30, reassignThreshold = 0,  
  mergeCutHeight = 0.25, saveTOMs = TRUE, saveTOMFileBase = "chyTOM",
  verbose = 3)

# Number of modules identified:
table(net_chy$colors)
# blue     brown     green      grey turquoise    yellow 
# 60        46        41       576        71        41 

# Change grey to NA to rm unassigned geneds
net_chy$colors[net_chy$colors == "grey"] <- "white"

net_cyano <- blockwiseModules(
  datExpr = t(vst_cyano),
  power = 12, 
  TOMType = "signed", minModuleSize = 30, reassignThreshold = 0,  
  mergeCutHeight = 0.25, saveTOMs = TRUE, saveTOMFileBase = "cyanoTOM",
  verbose = 3)

# Number of modules identified:
table(net_cyano$colors)
# blue     brown      grey turquoise    yellow 
# 107        92       126       180        50 

# Change grey to NA
net_cyano$colors[net_cyano$colors == "grey"] <- "white"

## ===============================
## 🎨 4. Dendrogram & colors
## ===============================
plotDendroAndColors(
  net_chy$dendrograms[[1]], 
  net_chy$colors[net_chy$blockGenes[[1]]], "Module colors", 
  dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)

plotDendroAndColors(
  net_cyano$dendrograms[[1]], 
  net_cyano$colors[net_cyano$blockGenes[[1]]], "Module colors", 
  dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)

## ==============================
## 5. Run enrichment on modules
## ==============================
res <- lapply(names(table(net_chy$colors)), function(i){
  getGOBubbleZ(universe = colnames(t(vst_chy)),
               annotation = annotationChytrid, 
               genelist = names(net_chy$colors)[net_chy$colors == i], 
               GO_df = GO_chytrid, isbubble = F)
})

res <- lapply(names(table(net_cyano$colors)), function(i){
  getGOBubbleZ(universe = colnames(t(vst_cyano)),
               annotation = annotationCyano, 
               genelist = names(net_cyano$colors)[net_cyano$colors == i], 
               GO_df = GO_cyano, isbubble = F)
})

## No significant GO term in any module

## =================
## 6. Find our DEG
## =================
table(na.omit(net_chy$colors[unique(fullDEGTable$geneName)]))
data.frame(modules=na.omit(net_chy$colors[unique(fullDEGTable$geneName)])) %>% 
  arrange(modules) ## white = no module

table(na.omit(net_cyano$colors[unique(fullDEGTable$geneName)]))
data.frame(modules=na.omit(net_cyano$colors[unique(fullDEGTable$geneName)])) %>% 
  arrange(modules) ## white = no module

## ======================================
## 7. Test association with treatment
## ======================================

treatment <- ifelse(grepl("^control", colnames(vst_chy)), 0, 1)
moduleTraitCor = cor(net_chy$MEs, treatment, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, length(treatment))
moduleTraitCor; moduleTraitPvalue
labeledHeatmap(Matrix = moduleTraitCor, xLabels = "Treatment",
               yLabels = names(net_chy$MEs), colors = blueWhiteRed(50))

treatment <- ifelse(grepl("chy", colnames(vst_chy)), 0, 1)
moduleTraitCor <- cor(net_chy$MEs, treatment, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, length(treatment))
moduleTraitCor; moduleTraitPvalue
labeledHeatmap(Matrix = moduleTraitCor, xLabels = "Treatment",
               yLabels = names(net_chy$MEs), colors = blueWhiteRed(50))

treatment <- ifelse(grepl("^control", colnames(vst_cyano)), 0, 1)
moduleTraitCor <- cor(net_cyano$MEs, treatment, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, length(treatment))
moduleTraitCor; moduleTraitPvalue
labeledHeatmap(Matrix = moduleTraitCor, xLabels = "Treatment",
               yLabels = names(net_cyano$MEs), colors = blueWhiteRed(50))

treatment <- ifelse(grepl("cyano", colnames(vst_cyano)), 0, 1)
moduleTraitCor <- cor(net_cyano$MEs, treatment, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, length(treatment))
moduleTraitCor; moduleTraitPvalue
labeledHeatmap(Matrix = moduleTraitCor, xLabels = "Treatment",
               yLabels = names(net_cyano$MEs), colors = blueWhiteRed(50))

#####################################
## II. Dual co-expression analysis ##
#####################################
# We keep only "both" (dual transcriptome)
# vst_cyano: genes x samples
# vst_chy: genes x samples
shared_samples <- intersect(colnames(vst_cyano), colnames(vst_chy))
# Keep same sample order
vst_cyano_both <- vst_cyano[, shared_samples]
vst_chy_both <- vst_chy[, shared_samples]

# Make sure samples are in same order:
vst_cyano_both <- vst_cyano_both[, colnames(vst_chy_both)]
stopifnot(all(colnames(vst_cyano_both) == colnames(vst_chy_both)))

## combined
vst_combined <- rbind(vst_chy_both, vst_cyano_both)

## ===============================
## ✅ 2.Pick soft power threshold
## ===============================
pickSoftPow(vst_combined)
softPower_combined <- 1 # to avoid too high connectivity

## ===============================
## 🔗 3. Build modules
## ===============================
cor <- WGCNA::cor

## Biweight midcorrelation (bicor) is a robust correlation that downweights outliers — it’s one of the best choices for cross-species WGCNA.

net_combined <- blockwiseModules(
  datExpr = t(vst_combined),
  power = softPower_combined, 
  corFnc = "bicor",
  networkType = "signed",
  TOMType = "signed",
  minModuleSize = 30,
  mergeCutHeight = 0.25,
  reassignThreshold = 0,
  verbose = 5
)

# Number of modules identified:
table(net_combined$colors)

## ===============================
## 🎨 4. Dendrogram & colors
## ===============================
pdf("../../figures/Fig5a.pdf", width = 15, height = 5)
plotDendroAndColors(
  net_combined$dendrograms[[1]], 
  net_combined$colors[net_combined$blockGenes[[1]]], "Module colors", 
  dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
dev.off()

## ==============================
## 5. Run enrichment on modules
## ==============================
res <- lapply(names(table(net_combined$colors)), function(i){
  getGOBubbleZ(universe = colnames(t(vst_combined)),
               annotation = rbind(annotationChytrid[
                 names(annotationChytrid) %in%
                   c("custom_gene_name", "gene_name", "GO.accession")],
                 annotationCyano), 
               genelist = names(net_combined$colors)[net_combined$colors == i], 
               GO_df = rbind(GO_chytrid,GO_cyano), isbubble = F)
})

## No significant GO term in any module

## =================
## 6. Find our DEG
## =================
table(na.omit(net_combined$colors[unique(fullDEGTable$geneName)]))
data.frame(modules=na.omit(net_combined$colors[unique(fullDEGTable$geneName)])) %>% 
  arrange(modules) 

## blue, yellow, and greenyellow have a lot of DEG

## ======================================
## 7. Test association with MET treatment
## ======================================

treatment <- ifelse(grepl("^control", colnames(vst_combined)), 0, 1)
moduleTraitCor = cor(net_combined$MEs, treatment,use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, length(treatment))
moduleTraitCor; moduleTraitPvalue
# MEgreenyellow: cor = 0.85060744; p-value= 0.01526171

# If you plot eigengene boxplots for these 3, you’ll see how they change by treatment.
meta <- data.frame(
  Sample = rownames(net_combined$MEs),
  Treatment = sub("_.*", "",  rownames(net_combined$MEs)),
  stringsAsFactors = TRUE)
plot_df <- net_combined$MEs %>%
  as.data.frame() %>%
  mutate(Sample = rownames(.)) %>%
  left_join(meta, by = "Sample")

pdf("../../figures/Fig5b.pdf", width = 3, height = 3)
ggplot(plot_df, aes(x = Treatment, y = MEgreenyellow, fill = Treatment)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, size = 2) +
  labs(title = "Module Eigengene for MEgreenyellow",
       y = "Module Eigengene Value") +
  theme_minimal() + theme(legend.position = "none")
dev.off()
# ## as well as the hub genes (core genes for a module). We could focus only on modules
# (groups of super co-expressed genes) that have a lot of host AND parasite genes together, 
# or plotting all. I'll prepare a plot with all and we can discuss to focus on the most interesting. 
# I'll try to find time to do that ASAP.
#
## ✅ (1) Find the ratio of both species per module
# Get gene names
chy_genes <- row.names(vst_chy)
cyano_genes <- row.names(vst_cyano)

# Combine into a data.frame: gene + species + module
gene_species_df <- data.frame(
  gene = names(net_combined$colors),
  module = net_combined$colors,
  species = ifelse(names(net_combined$colors) %in% chy_genes, "Chytrid",
                   ifelse(names(net_combined$colors) %in% cyano_genes, "Cyano", NA))
)

# Count table
table_per_module <- table(gene_species_df$module, gene_species_df$species)
print(table_per_module)

# Proportions per module
prop_table <- prop.table(table_per_module, margin = 1)
print(round(prop_table, 2))

## ✅ (2) Find top 10% hub genes in each module

# Calculate adjacency
adjacency <- adjacency(t(vst_combined), power = softPower_combined)

# Calculate intramodular connectivity (kWithin)
IMconn <- intramodularConnectivity(adjacency, net_combined$colors)

# Add module membership too:
ME <- moduleEigengenes(t(vst_combined), net_combined$colors)$eigengenes
kME <- signedKME(t(vst_combined), ME)

# Combine results
hub_df <- data.frame(
  gene = rownames(vst_combined),
  module = net_combined$colors,
  kWithin = IMconn$kWithin,
  kME = apply(kME, 1, max)
)

# Get top 10% per module by kWithin:
hub_genes <- hub_df %>%
  group_by(module) %>%
  slice_max(order_by = kWithin, prop = 0.1)

# Add a column for organism type (host vs. parasite)
hub_genes <- hub_genes %>%
  mutate(org = ifelse(gene %in% annotationChytrid$gene_name, "chytrid",
                      ifelse(gene %in% annotationCyano$gene_name, "cyano", NA)))

## Which hub genes are also DEG?
intersect(hub_genes$gene, fullDEGTable$geneName)
## "GeneID:77286325" "HDA1A_XENLA"     "RS15A_BOVIN" 

hub_genes[hub_genes$module %in% hub_genes[grep("GeneID:77286325", hub_genes$gene),"module"],]
hub_genes[hub_genes$module %in% hub_genes[grep("HDA1A_XENLA", hub_genes$gene),"module"],]
hub_genes[hub_genes$module %in% hub_genes[grep("RS15A_BOVIN", hub_genes$gene),"module"],]

## Save xls file for hub genes:
library(openxlsx)

summary_table <- hub_genes %>%
  group_by(module) %>%
  summarise(
    CyanoGenes = sum(org == "cyano"),
    ChytridGenes = sum(org == "chytrid"),
    CyanoHubGenes = paste(gene[org == "cyano"], collapse = ", "),
    ChytridHubGenes = paste(gene[org == "chytrid"], collapse = ", ")
  ) %>%
  arrange(desc(CyanoGenes + ChytridGenes))%>% # Order by total genes
  filter(CyanoGenes >0 & ChytridGenes > 0)  # rm if only one organism

modules2plot <- summary_table$module

# Create a workbook
wb <- createWorkbook()
addWorksheet(wb, "WGCNA Summary")

# Add WGCNA color codes (base R handles them as named colors)
summary_table <- summary_table %>%
  mutate(module_color = labels2colors(module))

# Write data (exclude the color column for display)
writeData(wb, sheet = 1, x = summary_table %>% select(-module_color), startCol = 1, startRow = 1)

# Highlight each cell in the module column with its module color
for (i in seq_len(nrow(summary_table))) {
  mod_col <- summary_table$module_color[i]
  
  # Set cell style with fill color
  style <- createStyle(fgFill = mod_col, fontColour = "#FFFFFF")  # White text
  addStyle(wb, sheet = 1, style = style, rows = i + 1, cols = 1, gridExpand = FALSE)
}

# Optional: make headers bold
headerStyle <- createStyle(textDecoration = "bold")
addStyle(wb, sheet = 1, style = headerStyle, rows = 1, cols = 1:5, gridExpand = TRUE)

# Save the Excel file
saveWorkbook(wb, "../../figures/WGCNA_summary_table.xlsx", overwrite = TRUE)

## ✅ (3)  Plot network highlighting hub genes
library(igraph)
library(tidygraph)
library(ggraph)
library(purrr)
library(dplyr)
library(rlang)
library(ggplot2)
library(tidyr)
library(patchwork)

# Generate enhanced network plots per module
plots <- hub_genes %>%  
  filter(module %in% modules2plot) %>%  # remove or adjust this filter as needed
  group_by(module) %>%
  group_split() %>%
  map(function(df) {
    genes <- df$gene
    mod <- unique(df$module)
    
    # Subset adjacency matrix to only these genes
    valid_genes <- intersect(genes, rownames(adjacency))
    if (length(valid_genes) < 2) return(NULL)
    
    adj_sub <- adjacency[valid_genes, valid_genes]
    
    # Build edge list (upper triangle)
    edge_df <- as.data.frame(as.table(adj_sub))
    colnames(edge_df) <- c("from", "to", "weight")
    edge_df <- edge_df %>%
      filter(from != to) %>%
      filter(as.numeric(factor(from)) < as.numeric(factor(to))) %>%
      filter(weight > 0.05)  # adjustable threshold
    
    if (nrow(edge_df) == 0) return(NULL)
    
    # Node metadata: org and kWithin
    nodes <- df %>%
      filter(gene %in% c(edge_df$from, edge_df$to)) %>%
      select(gene, org, kWithin) %>%
      distinct() %>%
      rename(name = gene)
    
    g <- graph_from_data_frame(edge_df, vertices = nodes, directed = FALSE)
    tg <- as_tbl_graph(g)
    
    # Plot using ggraph:
    ggraph(tg, layout = "fr") +
      geom_edge_link(aes(alpha = weight), color = "black", width = 0.2, show.legend = FALSE) +       # thin black edges
      geom_node_label(
        aes(label = gsub("GeneID:", "", sub("_.*", "", name)), 
            color = org,  # text color mapped to org
            fill = org), size = 4,
        label.r = unit(0.25, "lines"),       # rounded corners = oval feel
        label.padding = unit(0.25, "lines"),  # controls horizontal/vertical size
        fontface = "bold", label.size = 0
      ) +
      scale_fill_manual(values = c("cyano" = "white", "chytrid" = "black")) +
      scale_color_manual(values = c("cyano" = "black", "chytrid" = "white")) +
      guides(size = "none") +
      theme_void() +
      ggtitle(mod) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", colour = mod),
        legend.position = "none",
        plot.margin = unit(c(1,1,1,1), "cm")  # add margin
      )+   coord_cartesian(clip = "off")
  })

# Remove NULL plots if any
plots <- compact(plots) 

pdf("../../figures/Fig5c.network.pdf", width = 12, height = 1)
wrap_plots(plots, ncol = 4)
dev.off()


