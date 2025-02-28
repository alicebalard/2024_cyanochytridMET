# Load necessary libraries
library(WGCNA)
library(data.table)
library(DESeq2)

## Results of previous analysis
setwd("../Part2_DEG/")
source("fullAnalysis.R")

setwd("../Part3_DEG/")

# 'exprData' is the expression matrix (genes x samples)

findSoftPower <- function(exprData){
  # Perform WGCNA
  options(stringsAsFactors = FALSE)
  enableWGCNAThreads()
  set.seed(123) # For reproducibility
  
  # Choose a set of soft-thresholding powers
  powers = c(c(1:10), seq(from = 12, to=20, by=2))
  # Call the network topology analysis function
  sft = pickSoftThreshold(exprData, powerVector = powers, verbose = 5)
  
  # Plot the results:
  sizeGrWindow(9, 5)
  par(mfrow = c(1,2));
  cex1 = 0.9;
  
  # Scale-free topology fit index as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
       main = paste("Scale independence"));
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red");
  # this line corresponds to using an R^2 cut-off of h
  abline(h=0.90,col="red")
  
  # Mean connectivity as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], sft$fitIndices[,5],
       xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
       main = paste("Mean connectivity"))
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
}

findSoftPower(exprData = t(contrast_chytridgenome$vstr))
## Chytrid: N=12

findSoftPower(exprData = t(contrast_cyanogenome$vstr))
## Cyano: N=16

makeNetwork <- function(exprData, mypower, fileNamePlot, fileNameNet){
  ## Network
  net = blockwiseConsensusModules(
    fixDataStructure(exprData), power = mypower, minModuleSize = 30, deepSplit = 2,
    pamRespectsDendro = FALSE,
    mergeCutHeight = 0.25, numericLabels = TRUE,
    minKMEtoStay = 0,
    saveTOMs = TRUE, verbose = 5)
  
  consMEs = net$multiMEs;
  moduleLabels = net$colors;
  # Convert the numeric labels to color labels
  moduleColors = labels2colors(moduleLabels)
  consTree = net$dendrograms[[1]];
  
  ## Save for further
  save(consMEs, moduleLabels, moduleColors, consTree, file = fileNameNet)
  
  sizeGrWindow(8,6);
  pdf(file = "fileNamePlot", wi = 8, he = 6)
  plotDendroAndColors(consTree, moduleColors,
                      "Module colors",
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05,
                      main = "Consensus gene dendrogram and module colors")
  dev.off()
}

makeNetwork(exprData = t(contrast_chytridgenome$vstr), mypower = 12, 
            fileNamePlot = "../../data/ConsensusDendrogram-chytrid.pdf", 
            fileNameNet = "../../data/Network_chytrid.RData")

makeNetwork(exprData = t(contrast_cyanogenome$vstr), mypower = 16, 
            fileNamePlot = "../../data/ConsensusDendrogram-cyano.pdf", 
            fileNameNet = "../../data/Network_cyano.RData")

# Figures: Gene dendrogram obtained by clustering the dissimilarity based on consensus Topological Overlap with
# the corresponding module colors indicated by the color row.
# B. Zhang and S. Horvath. A general framework for weighted gene co-expression network analysis. Statistical
# Applications in Genetics and Molecular Biology, 4(1):Article 17, 2005.

## 1. Chytrids

## Quantifying module–trait associations
# Define numbers of genes and samples

datExpr = t(contrast_chytridgenome$vstr)
## load species specific consMEs, moduleLabels, moduleColors, consTree
load("../../data/Network_chytrid.RData")

# Get data traits
datTraits0 <- data.frame(sample_id = rownames(datExpr))

# Split the sample_id column into three new columns
datTraits <- data.frame(do.call(rbind, strsplit(as.character(datTraits0$sample_id), "_")))

# Rename the columns for clarity
colnames(datTraits) <- c("MetTrt", "InfTrt", "ID")
datTraits = datTraits[-3]

# Combine the original data with the split data
rownames(datTraits)=datTraits0$sample_id

datTraits$MetTrt <- as.numeric(as.factor(datTraits$MetTrt))
datTraits$InfTrt <- as.numeric(as.factor(datTraits$InfTrt))

nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
# Figure: Module-trait associations. Each row corresponds to a module eigengene, column to a trait. Each cell
# contains the corresponding correlation and p-value. The table is color-coded by correlation according to the color
# legend.

# Define variable weight containing the InfTrt column of datTrait
InfTrt = as.data.frame(datTraits$InfTrt);
names(InfTrt) = "InfTrt"
# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(datExpr, InfTrt, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(InfTrt), sep="");
names(GSPvalue) = paste("p.GS.", names(InfTrt), sep="");

## Chytrid: MEturquoise correlate with infection status (cor = 0.46, p = 0.05)

#############################################################################################
# Intramodular analysis: identifying genes with high Gene Significance and Module Membership

module = "turquoise"
column = match(module, modNames);
moduleGenes = moduleColors==module;
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for InfTrt",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
# Figure: A scatterplot of Gene Significance (GS) for weight vs. Module Membership (MM) in the module.
# There is a significant correlation between GS and MM in this module.

## GO on these: nothing
getGOBubbleZ(universe = universe_chytrid, annotation = annotationChytrid, 
             genelist = colnames(datExpr)[moduleColors=="turquoise"], 
             GO_df = GO_chytrid, isbubble = F)

## In which modules are DEG?
data.frame(moduleLabels, moduleColors) # 1 is turquoise, 2 blue

x = moduleLabels[
  names(moduleLabels)%in% fullDEGTable[grep("on chytrid gene expression", fullDEGTable$comparison), "geneName"]] 

x[x > 0]

## coexpressed DEG: c("RS15A_BOVIN", "ILVB_CRYNH")
# DEGs within associated modules are likely important for the condition-specific responses

contrast_chytridgenome$vstr[rownames(contrast_chytridgenome$vstr) %in% c("RS15A_BOVIN", "ILVB_CRYNH"),] %>%
  data.frame() %>%
  rownames_to_column("genes") %>%
  pivot_longer(
    cols = -genes,
    names_to = "sample",
    values_to = "value"
  ) %>%
  pivot_wider(
    id_cols = sample,
    names_from = genes,
    values_from = value
  ) %>% ggplot(aes(x=ILVB_CRYNH, y=RS15A_BOVIN)) + geom_point()

## 2. Cyano

## Quantifying module–trait associations
# Define numbers of genes and samples

datExpr = t(contrast_cyanogenome$vstr)
## load species specific consMEs, moduleLabels, moduleColors, consTree
load("../../data/Network_cyano.RData")

# Get data traits
datTraits0 <- data.frame(sample_id = rownames(datExpr))

# Split the sample_id column into three new columns
datTraits <- data.frame(do.call(rbind, strsplit(as.character(datTraits0$sample_id), "_")))

# Rename the columns for clarity
colnames(datTraits) <- c("MetTrt", "InfTrt", "ID")
datTraits = datTraits[-3]

# Combine the original data with the split data
rownames(datTraits)=datTraits0$sample_id

datTraits$MetTrt <- as.numeric(as.factor(datTraits$MetTrt))
datTraits$InfTrt <- as.numeric(as.factor(datTraits$InfTrt))

nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
# Figure: Module-trait associations. Each row corresponds to a module eigengene, column to a trait. Each cell
# contains the corresponding correlation and p-value. The table is color-coded by correlation according to the color
# legend.

# Define variable weight containing the InfTrt column of datTrait
InfTrt = as.data.frame(datTraits$InfTrt);
names(InfTrt) = "InfTrt"
# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(datExpr, InfTrt, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(InfTrt), sep="");
names(GSPvalue) = paste("p.GS.", names(InfTrt), sep="");

# Cyano: no significant correlation

#############################################################################################
## In which modules are DEG?
table(moduleColors)
data.frame(moduleLabels, moduleColors) # 1 is turquoise, 2 blue, 3 brown

x = moduleLabels[
  names(moduleLabels)%in% fullDEGTable[grep("on cyanobacteria gene expression", fullDEGTable$comparison), "geneName"]] 

x[x > 0]

## coexpressed DEG: 
# DEGs within associated modules are likely important for the condition-specific responses

## 1. 77286325, 77288421, 77288446, 77288523, 77288731, 77289571, 77289591, 77289792, 77290147
## 2. 77287852, folK
## 3. 77287257, 77287779

## 1.
# allophycocyanin, tetratricopeptide repeat protein, N-acetylmuramoyl-L-alanine amidase,
# Uma2 family endonuclease, Rne/Rng family ribonuclease, photosynthesis system II assembly factor Ycf48,
# hypothetical protein, hypothetical protein, aldo/keto reductase 
