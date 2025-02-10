source("R03_VolcanoPlots.R")

# (a)
# 360 only zoospore transcripts
# in how many samples are they found (>2)
RSEM_final_hope.gene_a <- RSEM_final_hope.gene[
  RSEM_final_hope.gene$whichOrg %in% "in_chytrid_alone" &
    RSEM_final_hope.gene$whichTranscriptome %in% "chytrid",]
## Select only chytrid samples
RSEM_final_hope.gene_a<- RSEM_final_hope.gene_a[
  grep("chy", names(RSEM_final_hope.gene_a))]

## At least gene exp non null in 2
RSEM_final_hope.gene_a <- RSEM_final_hope.gene_a[(rowSums(RSEM_final_hope.gene_a>0)>=2),]
nrow(RSEM_final_hope.gene_a) # 17 genes at least expressed in 2 samples
annotationChytrid$gene_name[match(rownames(RSEM_final_hope.gene_a), annotationChytrid$custom_gene_name)]
# "TBC31_ORYLA" "ERG7_GANLU"  "VPS21_YEAST" "FAS1_SCHPO"  "SNU13_USTMA" "ETFD_HUMAN"  "HGS_HUMAN"   "GPA1_NEUCR" 
# "TAF6B_ARATH" "MRP2_MOUSE"  "HRQ1_SCHPO"  "ALG8_DICDI"  "SWI6_CRYNH"  "YAMB_SCHPO"  "ARP3_ACACA"  "SPAZ_ARTOA" 
# "MNS1_CANAX" 

## Group a
getGOBubbleZ(universe_chytrid, annotationChytrid, GO_df=GO_chytrid,
             group=NA, isbubble=T, isDE=F, genelist=
               annotationChytrid$gene_name[match(rownames(RSEM_final_hope.gene_a), 
                                                 annotationChytrid$custom_gene_name)])

# no significant GO terms