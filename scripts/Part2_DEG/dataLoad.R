## 29th of Novemeber 2024
## A. Balard

## load all libraries needed for the project

#################
## Libraries load
list.of.packages <- c(
  "ggplot2", "reshape2","viridis", "pheatmap","RColorBrewer", "cowplot", "ggvenn",
  "WGCNA",
  # "goEnrichment",
  "stringr", # to modify characters
  "tidyverse")  # tidyverse will pull in ggplot2, readr, other useful libraries

message("Install CRAN packages if missing, and load CRAN packages...")

## install from CRAN and require all libraries from CRAN and github
install_if_missing <- function(packages, dependencies = TRUE) {
  new_packages <- packages[!(packages %in% installed.packages()[,"Package"])]
  if(length(new_packages)) {
    message("Installing missing packages: ", paste(new_packages, collapse = ", "))
    install.packages(new_packages, dependencies = dependencies)
  } else {
    message("All CRAN packages are already installed.\n")
  }
}

install_if_missing(list.of.packages)

message("if db error, downgrade dplyr")
# devtools::install_version("dbplyr", version = "2.3.4")

message("Loading CRAN packages...")
load_packages <- function(package_list) {
  not_loaded <- character()
  for (pkg in package_list) {
    if (suppressPackageStartupMessages(require(pkg, character.only = TRUE, quietly = TRUE))) {
      # Package loaded successfully, do nothing
    } else {
      not_loaded <- c(not_loaded, pkg)
    }
  }  
  if (length(not_loaded) > 0) {
    message("The following packages could not be loaded: ", paste(not_loaded, collapse = ", "))
  } else {
    message("All packages loaded successfully.\n")
  }
}

load_packages(list.of.packages)

##########################################
## install packages from github if not yet
packages_to_install <- c(
  "asishallab/goEnrichment")

install_and_load_github_packages <- function(packages) {
  # Ensure devtools is installed and loaded
  if (!requireNamespace("devtools", quietly = TRUE)) {
    install.packages("devtools")
  }
  library(devtools)
  
  not_loaded <- character()
  
  for (pkg in packages) {
    pkg_name <- strsplit(pkg, "/")[[1]][2]
    
    if (!requireNamespace(pkg_name, quietly = TRUE)) {
      message(paste("Installing", pkg, "from GitHub..."))
      tryCatch(
        devtools::install_github(pkg, quiet = TRUE),
        error = function(e) {
          message(paste("Error installing", pkg, ":", e$message))
        }
      )
    }
    
    if (suppressPackageStartupMessages(require(pkg_name, character.only = TRUE, quietly = TRUE))) {
      # Package loaded successfully, do nothing
    } else {
      not_loaded <- c(not_loaded, pkg_name)
    }
  }
  
  if (length(not_loaded) > 0) {
    message("The following packages could not be loaded: ", paste(not_loaded, collapse = ", "))
  } else {
    message("All packages loaded successfully.\n")
  }
}

message("Install github packages if missing, and load github packages...")
install_and_load_github_packages(packages_to_install)

#####################################################
## install from biocmanager and require all libraries
## Biocmanager packages 
bioc_packages <- c("Category", # for hypergeometric GO test
                   "EnhancedVolcano",
                   "GOstats", # for GO analysis
                   "GSEABase",  # for GO term GeneSetCollection
                   "GO.db", # for GO slim retrieving
                   "qvalue", # for FDR after PQLseq
                   "impute", "preprocessCore",#for WGCNA
                   "DESeq2"
) 
# devtools::install_version("dbplyr", version = "2.3.4") # if buggy UniProt.ws

install_and_load_bioc_packages <- function(package_list) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  
  not_loaded <- character()
  
  for (pkg in package_list) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      message(paste("Installing", pkg, "..."))
      BiocManager::install(pkg, update = FALSE, ask = FALSE)
    }
    
    if (suppressPackageStartupMessages(require(pkg, character.only = TRUE, quietly = TRUE))) {
      # Package loaded successfully, do nothing
    } else {
      not_loaded <- c(not_loaded, pkg)
    }
  }
  
  if (length(not_loaded) > 0) {
    message("The following packages could not be loaded: ", paste(not_loaded, collapse = ", "))
  } else {
    message("All packages loaded successfully.\n")
  }
}

message("Loading bioconductor packages...")
install_and_load_bioc_packages(bioc_packages)

######################################
## Files used for the DESeq2 analysis:

# --> On Alice local Thinkpad machine, files in the path: 
# /home/alice/Documents/GIT/2024_cyanochytridMET/ignoreThinkpad/assembly15nov24

## Chytrid
# assembly: /scratch/alicebalard/outData/assemblyMergedFungi/trinity_out_dir/Trinity_eukaryoteHits.fasta
# gene_trans_map: /scratch/alicebalard/outData/assemblyMergedFungi/trinity_out_dir/Trinity_eukaryoteHits.fasta.gene_trans_map
gene_trans_map_chytrid <- read.csv("../../gitignore/Trinity_eukaryoteHits.fasta.gene_trans_map", 
                                   header = F, sep=" ")

# annotation: /scratch/alicebalard/outData/assemblyMergedFungi/annotation/assemblyMergedFungi_filterEuk_simplified.tsv
annotationChytrid <- read.csv("../../gitignore/assemblyMergedFungi_filterEuk_simplified.tsv", sep = "\t")

## Cyanobacteria
# 1. header_lookup_table.txt This file contains the shorter names I tried for the cyanobacterium genes and the full name separated by a tab. For the short name I add a number that represents the gene (1...4428) and the g1 and i1 even though each gene has only one isoform. For example
# cyano_NZ_LR882952.1_cds_WP_027249510.1_1_g1_i1 >lcl|NZ_LR882952.1_cds_WP_027249510.1_1 [locus_tag=NMG88_RS00005] [db_xref=GeneID:77286138] [protein=TIGR01548 family HAD-type hydrolase] [protein_id=WP_027249510.1] [location=135..920] [gbkey=CDS]
header_lookup_table = read.csv("../../data/run_DESEQ2_Erika/header_lookup_table.txt", sep="\t", header = F)

# 2. gene_trans_map_cds.txt This file contains the short name and the shortest name separated by a tab. For example
# cyano_gene1 cyano_NZ_LR882952.1_cds_WP_027249510.1_1_g1_i1
gene_trans_map_cyano = read.csv("../../data/run_DESEQ2_Erika/gene_trans_map_cds.txt", sep="\t", header = F)

# extract protein name
header_lookup_table$protein = sub(".*\\[protein=(.*?)\\].*", "\\1", header_lookup_table$V2)
# Remove '>'
header_lookup_table$V1 = gsub(">", "", header_lookup_table$V1)  
# match with gene trans map
gene_trans_map_cyano$protein = header_lookup_table$protein[
  match(gene_trans_map_cyano$V2, header_lookup_table$V1)]

source("functions.R")
