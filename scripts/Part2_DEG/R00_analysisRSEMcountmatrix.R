## 29th of November 2024
source("libLoad.R")
source("dataLoad.R")

## Let's observe the count matrix calculated by Trinity
RSEM_final_hope.gene <- 
  read.csv("../../data/run_DESEQ2_Erika/RSEM_final_hope.gene.counts.matrix", sep="\t")

##########################################################################
## Split by type of genes depending on in which samples they are expressed
# (1) chytrid genes expressed in chytrid alone and infecting --> give to Erika for DESEq2 1623 genes
# (2) cyanobacteria genes expressed in cyano alone and infected --> give to Erika for DESEq2 3420 genes
# (3) a mix that we need to understand --> GO analysis
a = ifelse(rowSums(RSEM_final_hope.gene[grep("chy", names(RSEM_final_hope.gene))]) !=0,
           "in_chytrid_alone", "")
b = ifelse(rowSums(RSEM_final_hope.gene[grep("both", names(RSEM_final_hope.gene))]) !=0,
           "in_both_organisms", "")
c = ifelse(rowSums(RSEM_final_hope.gene[grep("cyano", names(RSEM_final_hope.gene))]) !=0,
           "in_cyano_alone", "")

RSEM_final_hope.gene$whichOrg <- trimws(paste(a,b,c, sep = " "))

RSEM_final_hope.gene$whichTranscriptome <- ifelse(
  grepl("TRINITY", RSEM_final_hope.gene$X), "chytrid", "cyano")

table(RSEM_final_hope.gene$whichTranscriptome,
      RSEM_final_hope.gene$whichOrg)
#            0  in_both_organisms in_both_organisms in_cyano_alone in_chytrid_alone
# chytrid  679               1764                              751              360
# cyano    629                  1                             3420                0
# 
#         in_chytrid_alone  in_cyano_alone in_chytrid_alone in_both_organisms
# chytrid                                5                               1623
# cyano                                  2                                  0
# 
#         in_chytrid_alone in_both_organisms in_cyano_alone in_cyano_alone
# chytrid                                               140             16
# cyano                                                 161            215

###################
## Define groups ##
###################
## Figure X ##

##### 1. rm contamination #####
# chytrid genes also found in when cyanobacteria is alone:
# chytrid + in_both_organisms in_cyano_alone
# chytrid + in_chytrid_alone  in_cyano_alone
# chytrid + in_chytrid_alone in_both_organisms in_cyano_alone
# chytrid + in_cyano_alone

listOfTranscriptContaminant_toRmFromChytridTranscriptome <-
  RSEM_final_hope.gene[
    RSEM_final_hope.gene$whichTranscriptome %in% "chytrid" & RSEM_final_hope.gene$whichOrg %in% 
      c("in_both_organisms in_cyano_alone",
        "in_chytrid_alone  in_cyano_alone",
        "in_chytrid_alone in_both_organisms in_cyano_alone",
        "in_cyano_alone"),"X"]

write.csv(listOfTranscriptContaminant_toRmFromChytridTranscriptome,
  "../../data/listOfTranscriptContaminant_toRmFromChytridTranscriptome", quote = F, row.names = F)

#### 2. Which genes have been sequenced for chytrid and cyano (and are not contamination)?
sequencedChytridGenes <- RSEM_final_hope.gene[
  RSEM_final_hope.gene$whichOrg %in% c("in_chytrid_alone in_both_organisms",
                                       "in_chytrid_alone", "in_both_organisms") &
    RSEM_final_hope.gene$whichTranscriptome %in% "chytrid","X"]
length(sequencedChytridGenes) # 3747 sequenced genes

sequencedCyanoGenes <- RSEM_final_hope.gene[
  RSEM_final_hope.gene$whichOrg %in% c("in_both_organisms in_cyano_alone",
                                       "in_cyano_alone", "in_both_organisms") &
    RSEM_final_hope.gene$whichTranscriptome %in% "cyano","X"]
length(sequencedCyanoGenes) # 3636 sequenced genes

### Genes for DESeq2 -> see next script R01 R02

### Genes of interest only present in one condition:
rownames(RSEM_final_hope.gene) <- RSEM_final_hope.gene$X
# (a)
# 360 only zoospore transcripts
# in how many samples are they found (>half)
RSEM_final_hope.gene_a <- RSEM_final_hope.gene[
      RSEM_final_hope.gene$whichOrg %in% "in_chytrid_alone" &
      RSEM_final_hope.gene$whichTranscriptome %in% "chytrid",]
## Select only chytrid samples
RSEM_final_hope.gene_a<- RSEM_final_hope.gene_a[
  grep("chy", names(RSEM_final_hope.gene_a))]

## At least gene exp non null in half of each group met/not met
RSEM_final_hope.gene_a = RSEM_final_hope.gene_a[
  (rowSums(RSEM_final_hope.gene_a[grep("control", names(RSEM_final_hope.gene_a))]>0)>
     ncol(RSEM_final_hope.gene_a[grep("control", names(RSEM_final_hope.gene_a))])/2) &
    (rowSums(RSEM_final_hope.gene_a[grep("met", names(RSEM_final_hope.gene_a))]>0)>
       ncol(RSEM_final_hope.gene_a[grep("met", names(RSEM_final_hope.gene_a))])/2),]
## none remaining!

# (b)
# 1764 chytrid in_both_organisms --> virulence genes of chytrid!!! 10 filaments with TONS of infection
RSEM_final_hope.gene_b <- RSEM_final_hope.gene[
  RSEM_final_hope.gene$whichOrg %in% "in_both_organisms" &
    RSEM_final_hope.gene$whichTranscriptome %in% "chytrid",]
## Select only infection samples
RSEM_final_hope.gene_b <- RSEM_final_hope.gene_b[grep("both", names(RSEM_final_hope.gene_b))]

## At least gene exp non null in half of each group met/not met
RSEM_final_hope.gene_b = RSEM_final_hope.gene_b[
  (rowSums(RSEM_final_hope.gene_b[grep("control", names(RSEM_final_hope.gene_b))]>0)>
     ncol(RSEM_final_hope.gene_b[grep("control", names(RSEM_final_hope.gene_b))])/2) &
    (rowSums(RSEM_final_hope.gene_b[grep("met", names(RSEM_final_hope.gene_b))]>0)>
       ncol(RSEM_final_hope.gene_b[grep("met", names(RSEM_final_hope.gene_b))])/2),]
## 3 remaining!

df = melt(RSEM_final_hope.gene_b)
df$trt = ifelse(grepl("met", df$variable), "met", "control")
t.test(df$value[df$trt %in% "met"], df$value[df$trt %in% "control"])
## All three non significant effect of Met

annotationChytrid[
  match(row.names(RSEM_final_hope.gene_b), annotationChytrid$custom_gene_name),]
# RLR78_PLAVT" "4CLL8_ARATH" "MAK5_CRYNB"
# RLR78_PLAVT = code for the secreted RxLR effector protein 78, virulence factor identified in the Downy mildew of grapevine
# GO for extracellular region

# 4CLL8_ARATH = code for 4-coumarate--CoA ligase-like 8 Arabidopsis thaliana (Mouse-ear cress)
## Identification of a peroxisomal acyl-activating enzyme involved in the biosynthesis of jasmonic acid in Arabidopsis.
## positble roles
# https://www.frontiersin.org/journals/fungal-biology/articles/10.3389/ffunb.2021.708813/full
# Host Interaction and Penetration
# The process of host penetration and colonization by parasitic chytrids involves complex interactions with the host cell wall and cytoplasm. A 4CL homolog might be involved in producing enzymes or metabolites that help in degrading the host cell wall or in forming structures necessary for penetration, such as appressoria or infection pegs1.
# Stress Response and Virulence
# Reactive Oxygen Species (ROS) play a significant role in fungal pathogenicity and stress response. A 4CL-like enzyme could be involved in pathways that help the fungus cope with ROS generated by the host as a defense mechanism, thereby enhancing its virulence and survival within the host1.


# (d)
# 215 --> cyano genes only expressed outside of infection?
## TO DO
# in how many samples are they found (>half)
# top expressed
# --> GO on these to find the big functions
# DEG MET/not MET --> does met affects the virulence?
RSEM_final_hope.gene_d <- RSEM_final_hope.gene[
  RSEM_final_hope.gene$whichOrg %in% "in_cyano_alone" &
    RSEM_final_hope.gene$whichTranscriptome %in% "cyano",]
## Select only cyano samples
RSEM_final_hope.gene_d <- RSEM_final_hope.gene_d[
  grep("cyano", names(RSEM_final_hope.gene_d))]

## At least gene exp non null in half of each group met/not met
RSEM_final_hope.gene_d = RSEM_final_hope.gene_d[
  (rowSums(RSEM_final_hope.gene_d[grep("control", names(RSEM_final_hope.gene_d))]>0)>
     ncol(RSEM_final_hope.gene_d[grep("control", names(RSEM_final_hope.gene_d))])/2) &
    (rowSums(RSEM_final_hope.gene_d[grep("met", names(RSEM_final_hope.gene_d))]>0)>
       ncol(RSEM_final_hope.gene_d[grep("met", names(RSEM_final_hope.gene_d))])/2),]

nrow(RSEM_final_hope.gene_d) # 124 remaining!!
# GO on those TBC --> see R03 script: no significant GO term

# DE on met? --> R01 no DEG found

df = melt(RSEM_final_hope.gene_d)
df$trt = ifelse(grepl("met", df$variable), "met", "control")
t.test(df$value[df$trt %in% "met"], df$value[df$trt %in% "control"])
## no difference between metolachlor values

## describe top 10

# (e)
## one gene cyano only expressed during infection
RSEM_final_hope.gene_e <- RSEM_final_hope.gene[
  RSEM_final_hope.gene$whichOrg %in% "in_both_organisms" &
    RSEM_final_hope.gene$whichTranscriptome %in% "cyano",]
## Select only cyano samples
RSEM_final_hope.gene_e <- RSEM_final_hope.gene_e[
  grep("both", names(RSEM_final_hope.gene_e))]

## At least gene exp non null in half of each group met/not met
RSEM_final_hope.gene_e[
  (rowSums(RSEM_final_hope.gene_e[grep("control", names(RSEM_final_hope.gene_e))]>0)>
     ncol(RSEM_final_hope.gene_e[grep("control", names(RSEM_final_hope.gene_e))])/2) &
    (rowSums(RSEM_final_hope.gene_e[grep("met", names(RSEM_final_hope.gene_e))]>0)>
       ncol(RSEM_final_hope.gene_e[grep("met", names(RSEM_final_hope.gene_e))])/2),]
## no