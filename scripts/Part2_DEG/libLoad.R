for (pkg in c("BiocManager", "remotes")) {  # remotes instead of devtools
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}

cran_packages <- c(
  "ggplot2", "ggrepel", "reshape2", "viridis", "pheatmap", "RColorBrewer",
  "cowplot", "ggvenn", "MASS", "png", "grid",
  "WGCNA", "data.table",
  "GOplot",
  "stringr",
  "tidyverse"
)

bioc_packages <- c(
  "Category",
  "clusterProfiler",
  "GOstats",
  "GSEABase",
  # "GOSim",  # deprecated in Bioconductor 3.23, removed
  "GO.db",
  "qvalue",
  "impute",
  "preprocessCore",
  "DESeq2",
  "rtracklayer",
  "openxlsx"
)

github_packages <- c(
  "asishallab/goEnrichment",
  "GuangchuangYu/clusterProfiler"
)

install_and_load <- function(packages, source = "cran") {
  not_loaded <- character()
  for (pkg in packages) {
    pkg_name <- basename(pkg)
    
    if (!requireNamespace(pkg_name, quietly = TRUE)) {
      message("Installing: ", pkg)
      tryCatch({
        switch(source,
               cran   = install.packages(pkg),
               bioc   = BiocManager::install(pkg, update = FALSE, ask = FALSE),
               github = remotes::install_github(pkg, quiet = TRUE)  # remotes not devtools
        )
      }, error = function(e) message("Failed to install ", pkg, ": ", e$message))
    }
    
    if (!suppressPackageStartupMessages(
      require(pkg_name, character.only = TRUE, quietly = TRUE))) {
      not_loaded <- c(not_loaded, pkg_name)
    }
  }
  
  if (length(not_loaded) > 0) {
    message("Could not load: ", paste(not_loaded, collapse = ", "))
  } else {
    message(source, " packages: OK")
  }
}

install_and_load(cran_packages,   source = "cran")
install_and_load(bioc_packages,   source = "bioc")
install_and_load(github_packages, source = "github")