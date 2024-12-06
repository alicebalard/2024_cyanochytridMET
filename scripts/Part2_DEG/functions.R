## List of functions used
# makeClusterWGCNA


makeClusterWGCNA <- function(datExpr){
  gsg = WGCNA::goodSamplesGenes(datExpr, verbose = 3)
  message(gsg$allOK)
  if (!gsg$allOK){
    # Optionally, print the gene and sample names that were removed:
    if (sum(!gsg$goodGenes)>0) 
      printFlush(paste("Removing genes:", paste(names(datExpr)[!gsg$goodGenes], collapse = ", ")));
    if (sum(!gsg$goodSamples)>0) 
      printFlush(paste("Removing samples:", paste(rownames(datExpr)[!gsg$goodSamples], collapse = ", ")));
    # Remove the offending genes and samples from the data:
    datExpr = datExpr[gsg$goodSamples, gsg$goodGenes]
  }
  
  # Cluster the samples to detect outliers
  sampleTree_chytrid = hclust(dist(datExpr), method = "average")
  
  # WGCNA::sizeGrWindow(12,9)
  par(cex = 0.6);par(mar = c(0,4,2,0))
  plot(sampleTree_chytrid, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
       cex.axis = 1.5, cex.main = 2)
}
