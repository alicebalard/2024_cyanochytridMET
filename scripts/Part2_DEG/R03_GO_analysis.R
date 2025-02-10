##############################
### Gene Ontology analysis ###
##############################
source("R02_VolcanoPlots.R")

###############
## GO in DEG ##
###############

##### GO analysis: 4 files needed
## 1. Genes sequenced in our study that are not contamination
universe_chytrid = annotationChytrid$gene_name[match(sequencedChytridGenes, annotationChytrid$custom_gene_name)]
universe_cyano = annotationCyano$gene_name[match(sequencedCyanoGenes, annotationCyano$custom_gene_name)]

## 2. an annotation file with "custom_gene_name", "gene_name", and "GO.accession" 
### annotationChytrid and annotationCyano

## 3. a GO file with "GO.accession", "GO.ontology", "GO.name"
### GO_chytrid and GO_cyano

## 4. a list of genes of interest in gene_name format
# mylistResDESEQ2 (results DESEq2, R01)
## rownames must be the genes names

# https://github.com/dadrasarmin/enrichment_analysis_for_non_model_organism
getGOBubbleZ <- function(universe, annotation, GO_df, group, isbubble=T, isDE=T, genelist=NA){
  if (isDE){ # by default after differential gene expression
    DESeqDF = mylistResDESEQ2[[group]][
    mylistResDESEQ2[[group]]$padj < 0.05 & 
      !is.na(mylistResDESEQ2[[group]]$padj),]
  
  # A list of genes of interest from DESeq
  genes = rownames(DESeqDF)
  } else {
    genes = genelist
  }
  
  # A table that matches Term to Gene ID
  term2gene = annotation[c("GO.accession", "gene_name")] %>%
    dplyr::rename("term" = "GO.accession", "gene" = "gene_name") %>% data.frame()
  
  # A table that matches Term to names
  term2name = GO_df[c("GO.accession", "GO.name")] %>% 
    dplyr::rename("term" = "GO.accession", "name" = "GO.name")
  
  enrichment <- clusterProfiler::enricher(gene = genes,
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




## For the comparison:
## Plot GO
resclassicGO <- list(getGOBubbleZ(universe = universe_chytrid, annotation = annotationChytrid, 
                                  GO_df = GO_chytrid, group = names(mylistResDESEQ2)[1], isbubble = F),
                     getGOBubbleZ(universe = universe_chytrid, annotation = annotationChytrid, 
                                  GO_df = GO_chytrid, group = names(mylistResDESEQ2)[2], isbubble = F),
                     getGOBubbleZ(universe = universe_chytrid, annotation = annotationChytrid, 
                                  GO_df = GO_chytrid, group = names(mylistResDESEQ2)[3], isbubble = F),
                     getGOBubbleZ(universe = universe_chytrid, annotation = annotationChytrid, 
                                  GO_df = GO_chytrid, group = names(mylistResDESEQ2)[4], isbubble = F),
                     getGOBubbleZ(universe = universe_cyano, annotation = annotationCyano, 
                                  GO_df = GO_cyano, group = names(mylistResDESEQ2)[5], isbubble = F),
                     getGOBubbleZ(universe = universe_cyano, annotation = annotationCyano, 
                                  GO_df = GO_cyano, group = names(mylistResDESEQ2)[6], isbubble = F),
                     getGOBubbleZ(universe = universe_cyano, annotation = annotationCyano, 
                                  GO_df = GO_cyano, group = names(mylistResDESEQ2)[8], isbubble = F))
names(resclassicGO) <- names(mylistResDESEQ2)[c(1:6,8)]

resBubbleZ <- list(getGOBubbleZ(universe = universe_chytrid, annotation = annotationChytrid, 
                                GO_df = GO_chytrid, group = names(mylistResDESEQ2)[3]))
names(resBubbleZ) <- names(mylistResDESEQ2)[c(3)]

resclassicGO$chytrid_met_effect_1org$enrichment
resBubbleZ$chytrid_met_effect_1org$enrichment

################
## group a: 124 cyano genes only expressed outside of infection,
## at least in half of each group met/not met
getGOBubbleZ(universe = universe_chytrid, annotation = annotationChytrid, 
             GO_df = GO_chytrid, group = names(mylistResDESEQ2)[1], isbubble = F, 
             isDE = F, genelist = annotationChytrid$gene_name[match(rownames(RSEM_final_hope.gene_a), 
                                                                  annotationChytrid$custom_gene_name)])


################
## group d: 124 cyano genes only expressed outside of infection,
## at least in half of each group met/not met
getGOBubbleZ(universe = universe_cyano, annotation = annotationCyano, 
             GO_df = GO_cyano, group = names(mylistResDESEQ2)[1], isbubble = F, 
             isDE = F, genelist = annotationCyano$gene_name[match(rownames(RSEM_final_hope.gene_d), 
                                                                  annotationCyano$custom_gene_name)])

