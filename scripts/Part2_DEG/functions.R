## List of functions used
# filterRSEMno2Nullpergp
# makeClusterWGCNA
# makeVolcano for volcano plot
# calculateContrasts calculate DESeq2 and contrast
# getGOBubbleZ for GO analysis

filterRSEMno3Nullpergp <- function(RSEM){
  # Reshape the data for easier grouping and filtering
  df_long <- RSEM %>%
    rownames_to_column("row_name") %>% # Keep row names as a column
    pivot_longer(-row_name, names_to = "sample", values_to = "value") %>%
    separate(sample, into = c("met", "inf", "sampleshort"), sep = "_")
  
  df_long$group <- paste(df_long$met, df_long$inf)
  df_long$sample <- paste(df_long$met, df_long$inf, df_long$sampleshort, sep ="_")
  
  # Filter rows where there are at least two samples with values >0 per group
  filtered_df <- df_long %>%
    group_by(row_name, group) %>%
    summarise(nonzero_count = sum(value > 0), .groups = "drop") %>%
    filter(nonzero_count >= 3) %>%
    inner_join(df_long, by = c("row_name", "group")) %>%
    select(-nonzero_count)
  
  # Reshape back to wide format, preserving zeros
  RSEM_2 <- filtered_df %>%
    select(row_name, value, sample) %>%
    pivot_wider(names_from = sample, values_from = value
    ) %>%
    column_to_rownames("row_name") %>%
    na.omit() %>% as.data.frame()
  
  RSEM_2 %>% select(all_of(names(RSEM)))
  
  return(RSEM_2)
}

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
  return(data.frame(t(datExpr)))
}

makeVolcano <- function(res, title, subtitle){
  results_df = as.data.frame(res)
  results_df = results_df[order(results_df$padj),]
  
  # Subset the results to keep only significant genes
  ressig = results_df[results_df$padj < 0.05 & !is.na(results_df$padj),]
  
  ## simplify gene ID name for cyano
  labs = gsub("GeneID:", "", row.names(results_df))
  
  ## Volcano plot
  plot = EnhancedVolcano(results_df,
                         lab = labs,
                         x = 'log2FoldChange', title = title,
                         y = 'padj', pCutoff = 0.05,
                         widthConnectors = .1,
                         arrowheads = TRUE,
                         max.overlaps = Inf,
                         boxedLabels = TRUE,
                         drawConnectors = TRUE, labSize = 4,
                         ylim = c(0, 4),
                         # ylim = c(0, max(-log10(results_df$padj)) * 1.2),
                         ylab = bquote(~-Log[10]~adjusted~italic(P))) +
    ggtitle(label = title, subtitle = subtitle) 
  
  return(list(signifGenes = ressig, plot = plot))
}

#################
# DESeq2 analysis- for this analysis I need un-normalized counts. 
# I need to have integers. My values are numeric, but According to a blog in which 
# someone asked about the same problem, Michael Love said: Just round the counts to integer. 
# There isn't a loss in precision to this operation (sampling variance from the experiment
# is much larger). https://support.bioconductor.org/p/105964/
calculateContrasts <- function(my_countsmatrix, my_org, my_samples=samples_data,
                               smallestGroupSize = 3){
  
  ## Define our 4 conditions, we are comparing 2 by 2
  my_groups=c("control_both","met_both", paste0("control_",my_org), paste0("met_",my_org))
  
  my_samples=my_samples[my_samples$sample %in% names(my_countsmatrix),]
  
  ## make a ddsr object comparing between conditions
  ddsr <- DESeqDataSetFromMatrix(countData = round(my_countsmatrix),
                                 colData = my_samples,
                                 design = ~ condition)
  ## Filter for low count (less than 10 in less than the smallest group size)
  ## see "prefiltering" in vignette
  keep = rowSums(counts(ddsr) >= 10) >= smallestGroupSize # already filtered before
  ddsr = ddsr[keep,]
  
  # Variance stabilising transformation for downstream analysis (done by default in DESeq)
  vstr <-  assay(varianceStabilizingTransformation(ddsr))
  
  ## Run the DE seq
  ddsr <- DESeq(ddsr)
  
  ## Calculate the contrasts (pairwise comparisons of interest)
  print(paste0("comparison of groups: ", my_groups[1], " vs ", my_groups[2]))
  resr_met_effect_2orgs <- results(ddsr, contrast=c("condition", my_groups[1], my_groups[2]), alpha=0.05) 
  print(summary(resr_met_effect_2orgs))
  print(sum(resr_met_effect_2orgs$padj < 0.05, na.rm=TRUE))
  
  print(paste0("comparison of groups: ", my_groups[3]," vs ", my_groups[4]))
  resr_met_effect_1org <- results(ddsr, contrast=c("condition", my_groups[3], my_groups[4]), alpha=0.05) 
  print(summary(resr_met_effect_1org))
  print(sum(resr_met_effect_1org$padj < 0.05, na.rm=TRUE))
  
  print(paste0("comparison of groups: ", my_groups[3], " vs ", my_groups[1]))
  resr_inf_effect_control <- results(ddsr, contrast=c("condition", my_groups[3], my_groups[1]), alpha=0.05) 
  print(summary(resr_inf_effect_control))
  print(sum(resr_inf_effect_control$padj < 0.05, na.rm=TRUE))
  
  print(paste0("comparison of groups: ", my_groups[4], " vs ", my_groups[2]))
  resr_inf_effect_met <- results(ddsr, contrast=c("condition", my_groups[4], my_groups[2]), alpha=0.05) 
  print(summary(resr_inf_effect_met))
  print(sum(resr_inf_effect_met$padj < 0.05, na.rm=TRUE))
  
  # plotMA
  plotMA(resr_met_effect_2orgs)
  plotMA(resr_met_effect_1org)
  plotMA(resr_inf_effect_control)
  plotMA(resr_inf_effect_met)
  
  return(list(resr_met_effect_2orgs=resr_met_effect_2orgs, 
              resr_met_effect_1org=resr_met_effect_1org,
              resr_inf_effect_control=resr_inf_effect_control,
              resr_inf_effect_met=resr_inf_effect_met, 
              vstr=vstr))
}

##### GO analysis: 4 files needed
## 1. Genes sequenced in our study that are not contamination
# universe_chytrid
# universe_cyano

## 2. an annotation file with "custom_gene_name", "gene_name", and "GO.accession" 
### annotationChytrid and annotationCyano

## 3. a GO file with "GO.accession", "GO.ontology", "GO.name"
### GO_chytrid and GO_cyano

## 4. a vector of genes of interest

# https://github.com/dadrasarmin/enrichment_analysis_for_non_model_organism
getGOBubbleZ <- function(universe, annotation, GO_df, group, isbubble=T, genelist=NA){
  # A table that matches Term to Gene ID
  term2gene = annotation[c("GO.accession", "gene_name")] %>%
    dplyr::rename("term" = "GO.accession", "gene" = "gene_name") %>% 
    unique %>% data.frame()
  
  ## Select for genes which have a GO term associated
  term2gene = term2gene[!is.na(term2gene$term),]
  genelist = genelist[genelist %in% term2gene$gene]
  universe = universe[universe %in% term2gene$gene]
  
  # A table that matches Term to names
  term2name = GO_df[c("GO.accession", "GO.name")] %>% 
    dplyr::rename("term" = "GO.accession", "name" = "GO.name") %>% unique
  
  enrichment <- clusterProfiler::enricher(gene = genelist,
                                          TERM2GENE = term2gene,
                                          TERM2NAME = term2name, 
                                          pvalueCutoff = 0.05,
                                          universe = universe,
                                          qvalueCutoff = 0.05, 
                                          pAdjustMethod = "fdr")
  
  ## Add GO info
  enrichmentRes = enrichment@result %>%
    rename("GO.name"="Description") 
  enrichmentRes = merge(enrichmentRes, GO_df)
  
  if (sum(enrichmentRes$p.adjust < 0.05)==0){
    message("no significant GO terms")
    return(list(enrichment=enrichmentRes, result="no significant GO terms"))
  } else {
    if (isbubble){
      ## Prepare data for plotting with GOplot
      # circledat takes two data frames as an input. The first one 
      # contains the results of the functional analysis and should have at least four 
      # columns (category, term, genes, adjusted p-value). Additionally, a data frame 
      # of the selected genes and their logFC is needed.
      dfGO=data.frame(ID=enrichment@result$ID,
                      Term = enrichment@result$Description,
                      Genes = gsub("/", ", ", enrichment@result$geneID),
                      adj_pval = enrichment@result$p.adjust) 
      dfGO$Category = ifelse(
        GO_df$GO.ontology[match(dfGO$Term, GO_df$GO.name)] == "molecular_function", "MF",
        ifelse(GO_df$GO.ontology[match(dfGO$Term, GO_df$GO.name)] == "biological_process", "BP",
               ifelse(GO_df$GO.ontology[match(dfGO$Term, GO_df$GO.name)] == "cellular_component", "CC", NA)))
      
      dfDE=data.frame(ID=rownames(DESeqDF),
                      logFC=DESeqDF$log2FoldChange,
                      adj.P.Val=DESeqDF$padj,
                      B=DESeqDF$baseMean)
      
      # Generate the plotting object
      circ <- circle_dat(dfGO, dfDE)
      
      if (nrow(circ[circ$adj_pval < 0.05,]) == 0){
        plot = GOBubble(circ, table.legend = F, title = group)
      } else {
        plot = GOBubble(circ, ID = T, labels = -log10(0.05), title = group)
      }
      return(list(enrichment=enrichmentRes, circ=circ, plot=plot))
    } else { # plot classic GO plot
      enrichmentRes$GeneRatio <- sapply(enrichmentRes$GeneRatio, function(x) eval(parse(text=x)))
      GOplot = enrichmentRes %>%
        dplyr::filter(p.adjust < 0.05) %>% 
        ggplot(aes(x=GO.ontology, y = factor(GO.name))) +
        geom_point(aes(color = p.adjust, size = GeneRatio)) +
        scale_color_gradient(
          name="adjusted\np-value", low = "red", high = "blue", 
          limits = c(0, 0.05), breaks = c(0, 0.02, 0.04), labels =c("0", "0.02", "0.04")) +
        scale_size_continuous(name = "% of genes")+
        theme_bw() + ylab("") + xlab("") +
        theme(legend.box.background = element_rect(fill = "#ebebeb", color = "#ebebeb"),
              legend.background = element_rect(fill = "#ebebeb", color = "#ebebeb"),
              legend.key = element_rect(fill = "#ebebeb", color = "#ebebeb"), # grey box for legend
              legend.position="top",
              axis.text.y = element_text(size = 8),  # Decrease y-axis text size
              axis.text.x = element_text(size = 8, hjust = 1)  # Increase x-axis text size 
        )+
        facet_wrap(.~fct_inorder(GO.ontology), scales = "free")
      return(list(enrichment=enrichmentRes, GOplot=GOplot))
    }
  }
}
