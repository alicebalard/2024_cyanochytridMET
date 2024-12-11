## 29th of Novemeber 2024
## A. Balard

## load all libraries needed for the project

#################
## Libraries load
list.of.packages <- c(
  "ggplot2", "reshape2","viridis", "pheatmap","RColorBrewer", "cowplot", "ggvenn",
  "MASS",
  "WGCNA",
  "GOplot", # plot GO results
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
  "asishallab/goEnrichment", "GuangchuangYu/clusterProfiler")

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
                   "clusterProfiler", # more GO functions
                   "EnhancedVolcano",
                   "GOstats", # for GO analysis
                   "GSEABase",  # for GO term GeneSetCollection
                   "GOSim",# link GO terms and ID
                   "GO.db", # for GO slim retrieving
                   "qvalue", # for FDR after PQLseq
                   "impute", "preprocessCore",#for WGCNA
                   "DESeq2",
                   "rtracklayer"
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
