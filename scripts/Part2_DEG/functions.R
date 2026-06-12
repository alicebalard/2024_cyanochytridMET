## List of functions used
# makeClusterWGCNA
# makeVolcano for volcano plot
# calculateContrasts calculate DESeq2 and contrast
# getGOBubbleZ for GO analysis
# build_comparison_matrix & render_comparison_matrix Build a "comparison matrix" from DESeq2 result objects

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
makeVolcano <- function(res, title, mylogo, positionLogoStart = -2,
                        positionLogoStop = 2, top_n = 100,
                        label_fun = function(g) sub("_.*", "", gsub("GeneID:", "", g))) {
  results_df <- as.data.frame(res)
  results_df <- results_df[order(results_df$padj), ]
  
  # Subset the results to keep only significant genes
  ressig <- results_df[results_df$padj < 0.05 & !is.na(results_df$padj), ]
  
  ## label each gene via label_fun (applied per-row, so scalar helpers work too)
  results_df$labs <- unname(vapply(row.names(results_df), label_fun, character(1)))
  
  ## Volcano plot
  results_df <- results_df %>%
    mutate(
      negLog10Padj = -log10(padj),
      sig = ifelse(padj < 0.05, "Significant", "Not significant")
    )
  
  # Top N significant genes to label (ordered by padj)
  top_genes <- results_df %>%
    filter(padj < 0.05) %>%
    arrange(padj) %>%
    head(top_n)
  
  # Basic volcano plot
  plot <- ggplot(results_df, aes(x = log2FoldChange, y = negLog10Padj)) +
    geom_hline(yintercept = -log10(0.05), colour = "grey", linetype = "dashed") +
    geom_vline(xintercept = -1, colour = "grey", linetype = "dashed") +
    geom_vline(xintercept = 1, colour = "grey", linetype = "dashed") +
    geom_point(aes(color = sig), size = 2) +
    scale_color_manual(values = c("Significant" = "red", "Not significant" = "grey")) +
    theme_minimal() +
    labs(title = title,
         x = expression(Log[2]~Fold~Change),
         y = expression(-Log[10]~adjusted~italic(P))) +
    theme(legend.position = "none")
  
  # Label only top N genes
  plot <- plot +
    geom_label_repel(
      data = top_genes,
      aes(label = labs),
      max.overlaps = Inf,
      size = 4,
      box.padding = 0.3,
      point.padding = 0.2,
      segment.size = 0.5, segment.colour = "grey",
      arrow = arrow(length = unit(0.02, "npc"))
    ) + ylim(0, 5) + xlim(-5, 5)
  
  # Load your image
  img <- png::readPNG(mylogo)
  logo_grob <- rasterGrob(img, interpolate = TRUE)
  plot <- plot +
    annotation_custom(logo_grob, xmin = positionLogoStart, xmax = positionLogoStop,
                      ymin = 3, ymax = 4)
  return(list(signifGenes = ressig, plot = plot))
}

getGOBubbleZ <- function(universe, annotation, GO_df, isbubble=FALSE, genelist=NA) { #*** removed unused 'group' parameter
  
  # Build term2gene - handle list column for GO.accession #***
  term2gene <- annotation %>%
    dplyr::select(gene_name, GO.accession) %>%
    unnest(GO.accession) %>%
    filter(!is.na(GO.accession)) %>%
    distinct() %>%
    dplyr::rename(term=GO.accession, gene=gene_name) %>%
    dplyr::select(term, gene) %>%  # term FIRST, gene SECOND #***
    data.frame()
  
  # Filter genelist and universe to genes with GO annotations
  genelist <- genelist[genelist %in% term2gene$gene]
  universe <- universe[universe %in% term2gene$gene]
  
  if(length(genelist) == 0) {
    message("No DEGs with GO annotations")
    return(NULL)
  }
  
  # term2name
  term2name <- GO_df %>%
    dplyr::select(GO.accession, GO.name) %>%
    dplyr::rename(term = GO.accession, name = GO.name) %>%
    distinct()
  
  # Run enrichment
  enrichment <- clusterProfiler::enricher(
    gene          = genelist,
    TERM2GENE     = term2gene,
    TERM2NAME     = term2name,
    pvalueCutoff  = 0.05,
    universe      = universe,
    qvalueCutoff  = 0.05,
    minGSSize     = 1,        #***
    pAdjustMethod = "fdr")
  
  if(is.null(enrichment)) {
    message("No enrichment results")
    return(NULL)
  }
  
  # Add GO info
  enrichmentRes <- enrichment@result %>%
    left_join(GO_df, by = c("ID" = "GO.accession"))
  
  if(sum(enrichmentRes$p.adjust < 0.05, na.rm=TRUE) == 0) {
    message("No significant GO terms")
    return(list(enrichment=enrichmentRes, result="no significant GO terms"))
  }
  
  if(isbubble) {
    dfGO <- data.frame(
      ID      = enrichment@result$ID,
      Term    = enrichment@result$Description,
      Genes   = gsub("/", ", ", enrichment@result$geneID),
      adj_pval = enrichment@result$p.adjust)
    
    dfGO$Category <- case_when(                                      #*** replaces nested ifelse
      GO_df$GO.ontology[match(dfGO$Term, GO_df$GO.name)] == "molecular_function"  ~ "MF",
      GO_df$GO.ontology[match(dfGO$Term, GO_df$GO.name)] == "biological_process"  ~ "BP",
      GO_df$GO.ontology[match(dfGO$Term, GO_df$GO.name)] == "cellular_component"  ~ "CC",
      TRUE ~ NA_character_)
    
    dfDE <- data.frame(
      ID        = rownames(DESeqDF),
      logFC     = DESeqDF$log2FoldChange,
      adj.P.Val = DESeqDF$padj,
      B         = DESeqDF$baseMean)
    
    circ <- circle_dat(dfGO, dfDE)
    plot <- if(nrow(circ[circ$adj_pval < 0.05, ]) == 0) {
      GOBubble(circ, table.legend=FALSE)
    } else {
      GOBubble(circ, ID=TRUE, labels=-log10(0.05))
    }
    return(list(enrichment=enrichmentRes, circ=circ, plot=plot))
    
  } else {
    enrichmentRes$GeneRatio <- sapply(
      enrichmentRes$GeneRatio, function(x) eval(parse(text=x)))
    
    GOplot <- enrichmentRes %>%
      filter(p.adjust < 0.05) %>%
      ggplot(aes(x=GO.ontology, y=factor(GO.name))) +
      geom_point(aes(color=p.adjust, size=GeneRatio)) +
      scale_color_gradient(
        name   = "adjusted\np-value",
        low    = "red", high = "blue",
        limits = c(0, 0.05),
        breaks = c(0, 0.02, 0.04),
        labels = c("0", "0.02", "0.04")) +
      scale_size_continuous(name="% of genes") +
      theme_bw() + ylab("") + xlab("") +
      theme(
        legend.box.background = element_rect(fill="#ebebeb", color="#ebebeb"),
        legend.background     = element_rect(fill="#ebebeb", color="#ebebeb"),
        legend.key            = element_rect(fill="#ebebeb", color="#ebebeb"),
        legend.position       = "top",
        axis.text.y           = element_text(size=8),
        axis.text.x           = element_text(size=8, hjust=1)) +
      facet_wrap(.~fct_inorder(GO.ontology), scales="free")
    
    return(list(enrichment=enrichmentRes, GOplot=GOplot))
  }
}

####################################################################
## Reusable raw-count DEG plotter (used for both chytrid & cyano) ##
####################################################################
plot_degs_raw <- function(contrast_list, count_matrix, contrast_conditions,
                          condition_levels, colour_map,
                          label_fun = function(g) sub("_.*", "", g), ncol = 5,
                          label_samples   = NULL,                       # <- NEW: conditions to label
                          sample_label_fun = function(s) sub(".*_", "", s)) {  # <- NEW: e.g. "control_both_In2" -> "In2"
  
  getGenes <- function(x) rownames(x)[!is.na(x$padj) & x$padj < 0.05]
  
  deg_data <- lapply(names(contrast_list), function(name) {
    x <- contrast_list[[name]]; if (is.null(x)) return(NULL)
    genes <- getGenes(x);       if (length(genes) == 0) return(NULL)
    conds <- contrast_conditions[[name]]
    sig <- as.data.frame(x) %>% filter(!is.na(padj) & padj < 0.05) %>%
      rownames_to_column("gene") %>% dplyr::select(gene, padj) %>%
      mutate(cond1 = conds[1], cond2 = conds[2])
    cnts <- count_matrix[rownames(count_matrix) %in% genes, , drop = FALSE] %>%
      as.data.frame() %>% rownames_to_column("gene") %>%
      pivot_longer(-gene, names_to = "sample", values_to = "count") %>%
      mutate(condition = sub("_[^_]+$", "", sample))
    list(counts = cnts, sig = sig)
  }) %>% Filter(Negate(is.null), .)
  
  counts_all <- bind_rows(lapply(deg_data, `[[`, "counts")) %>% distinct() %>%
    mutate(condition = factor(condition, levels = condition_levels))
  sig_all <- bind_rows(lapply(deg_data, `[[`, "sig")) %>%
    mutate(stars = case_when(padj < 0.001 ~ "***", padj < 0.01 ~ "**", padj < 0.05 ~ "*"))
  
  genes <- sort(unique(counts_all$gene))
  message("Total unique DEGs: ", length(genes))
  
  log_max  <- max(log10(counts_all$count + 1), na.rm = TRUE)
  max_bars <- if (nrow(sig_all)) max(table(sig_all$gene)) else 0
  y_top    <- 10^(log_max + 0.3 * (max_bars + 1))
  
  plots <- lapply(genes, function(g) {
    df <- counts_all %>% filter(gene == g)
    sig_df <- sig_all %>% filter(gene == g) %>%
      mutate(y_pos = log_max + 0.3 * row_number())
    
    pos <- position_jitter(width = 0.15, height = 0, seed = 42)  # <- shared jitter
    
    p <- ggplot(df, aes(condition, count + 1, colour = condition)) +
      geom_violin(aes(fill = condition), alpha = 0.3) +
      geom_boxplot(width = 0.3, alpha = 0.7, outlier.shape = NA) +
      geom_point(position = pos, size = 1.5, alpha = 0.8) +        # <- was geom_jitter
      scale_y_log10(limits = c(1, y_top)) +
      scale_colour_manual(values = colour_map) +
      scale_fill_manual(values = colour_map) +
      ggtitle(label_fun(g)) + labs(x = NULL, y = "Count + 1 (log10)") +
      coord_cartesian(clip = "off") + theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "none", plot.title = element_text(size = 9, face = "bold"))
    
    # ---- NEW: label individual sample points for the requested conditions ----
    if (!is.null(label_samples))
      p <- p + ggrepel::geom_label_repel(
        aes(label = ifelse(condition %in% label_samples,
                           sample_label_fun(sample), NA_character_)),
        position = pos, na.rm = TRUE,           # same seeded jitter -> labels match points
        colour = "grey20", fill = alpha("white", 0.7),
        size = 2.3, label.padding = 0.12, label.size = 0.2,
        min.segment.length = 0, max.overlaps = Inf, show.legend = FALSE)
    
    for (i in seq_len(nrow(sig_df)))
      p <- p + ggsignif::geom_signif(
        comparisons = list(c(sig_df$cond1[i], sig_df$cond2[i])),
        annotations = sig_df$stars[i], y_position = sig_df$y_pos[i],
        tip_length = 0.02, colour = "black", size = 0.4, textsize = 4)
    p
  })
  cowplot::plot_grid(plotlist = plots, ncol = ncol)
}

## ============================================================
## Build a "comparison matrix" from DESeq2 result objects.
##
## Conditions are a 2-factor design encoded as "<treatment>_<lifestyle>"
##   treatment: control | met        (no MET vs MET)
##   lifestyle: both | chy | cyano   (in co-culture vs alone)
##
## For every PAIR of conditions that differs in exactly ONE factor:
##   - lower triangle  -> captured summary(res)  (only if that contrast was run)
##   - upper triangle  -> auto-generated description of what the contrast means
## Pairs differing in BOTH factors, and the diagonal, stay "-".
## ============================================================

build_comparison_matrix <- function(res_list,
                                    contrast_conditions,
                                    conditions,
                                    display     = conditions,
                                    organism    = c("chytrid", "cyano"),
                                    sep         = "\n",
                                    not_tested  = "(not tested)") {
  organism <- match.arg(organism)
  
  ## --- factor labels (the only organism-specific wording) -------------
  treat_lab <- c(control = "no MET",  met = "MET")
  if (organism == "chytrid") {
    org_word  <- "chytrid"
    life_word <- c(both = "infecting chytrid", chy = "chytrid alone")
  } else {
    org_word  <- "cyanobacteria"
    life_word <- c(both = "cyanobacteria in co-culture", cyano = "cyanobacteria alone")
  }
  
  parse_cond <- function(x) {
    p <- strsplit(x, "_")[[1]]
    list(treatment = p[1], lifestyle = p[2])
  }
  pair_key <- function(a, b) paste(sort(c(a, b)), collapse = "|")
  
  ## --- lookup: unordered condition pair -> result object --------------
  res_list <- res_list[names(contrast_conditions)]          # align order
  keys     <- vapply(contrast_conditions,
                     function(p) pair_key(p[1], p[2]), character(1))
  res_lookup <- setNames(res_list, keys)
  
  n <- length(conditions)
  M <- matrix("-", n, n, dimnames = list(display, display))
  
  for (i in seq_len(n)) for (j in seq_len(n)) {
    if (i == j) next
    a <- conditions[i]; b <- conditions[j]
    ca <- parse_cond(a); cb <- parse_cond(b)
    same_treat <- ca$treatment == cb$treatment
    same_life  <- ca$lifestyle == cb$lifestyle
    if (same_treat == same_life) next           # differ in both -> leave "-"
    
    key <- pair_key(a, b)
    has_res <- key %in% names(res_lookup)
    
    if (i > j) {
      ## lower triangle: the actual DESeq2 summary
      M[i, j] <- if (has_res) {
        sm <- utils::capture.output(summary(res_lookup[[key]]))
        sm <- sm[!grepl("^\\[[12]\\]", sm)]   # drop the [1]/[2] footnote lines
        sm <- sm[nzchar(trimws(sm))]          # drop blank lines summary() adds
        paste(sm, collapse = sep)
      } else not_tested
    } else {
      ## upper triangle: auto-generated description
      M[i, j] <- if (!same_treat) {            # treatment differs -> MET effect
        sprintf("Effect of MET on %s", life_word[[ca$lifestyle]])
      } else {                                 # lifestyle differs -> infection effect
        qual <- if (ca$treatment == "met") "presence" else "absence"
        sprintf("Effect of infection on %s, in %s of MET", org_word, qual)
      }
      if (!has_res) M[i, j] <- paste0(M[i, j], " ", not_tested)
    }
  }
  M
}

## ------------------------------------------------------------
## Render the matrix as an HTML table (for an Rmd report).
## Multi-line cells are turned into <br>; summaries shown monospace.
## ------------------------------------------------------------
render_comparison_matrix <- function(M, caption = NULL) {
  br <- function(x) gsub("\n", "<br>", x, fixed = TRUE)
  body <- apply(M, c(1, 2), br)
  dimnames(body) <- lapply(dimnames(M), br)
  knitr::kable(body, format = "html", escape = FALSE, caption = caption,
               align = "l", table.attr = 'style="font-family:monospace;
                                          font-size:11px; white-space:pre;"')
}

## ------------------------------------------------------------
## Render SEVERAL comparison matrices stacked in ONE HTML table,
## one labelled section per organism (needs kableExtra).
##   `mats`    : named list of matrices; names become section headers
##   All matrices must share the SAME column (display) labels.
## ------------------------------------------------------------
render_comparison_matrices <- function(mats, caption = NULL) {
  br <- function(x) gsub("\n", "<br>", x, fixed = TRUE)
  
  combined <- do.call(rbind, mats)            # requires identical colnames
  combined <- apply(combined, c(1, 2), br)
  colnames(combined) <- br(colnames(combined))
  
  k <- knitr::kable(
    combined, format = "html", escape = FALSE, align = "l", caption = caption,
    table.attr = 'style="font-family:monospace; font-size:11px; white-space:pre;"')
  
  start <- 1L                                  # add one section header per matrix
  for (nm in names(mats)) {
    n <- nrow(mats[[nm]])
    k <- kableExtra::pack_rows(k, nm, start, start + n - 1L)
    start <- start + n
  }
  k
}


## ------------------------------------------------------------
## Stack several per-organism matrices into ONE rendered table.
## `mats` = named list of character matrices that share the SAME
## column names (build each with the same generic `display`).
## Group headers come from the list names. Renders inline and,
## if `file` is given, saves the *rendered* table (HTML / PNG / PDF).
## Requires the kableExtra package.
## ------------------------------------------------------------
render_combined_matrix <- function(mats, file = NULL, caption = NULL,
                                   font_size  = 11,
                                   col_widths = NULL,          # CSS width per column, incl. "Comparison"
                                   blank      = c("-", "(not tested)"),
                                   ...) {                      # ... -> save_kable -> webshot (vwidth, zoom, ...)
  br <- function(x) gsub("\n", "<br>", x, fixed = TRUE)
  
  ## one data frame per organism, row label moved into a column
  dfs <- lapply(mats, function(M) {
    d <- as.data.frame(M, stringsAsFactors = FALSE, check.names = FALSE)
    cbind(Comparison = rownames(M), d, row.names = NULL,
          stringsAsFactors = FALSE)
  })
  combined <- do.call(rbind, dfs)
  
  ## hide empty contrasts: blank out placeholder / not-run cells
  combined[] <- lapply(combined, function(col)
    ifelse(col %in% blank, "", col))
  
  combined[]      <- lapply(combined, br)          # newlines -> <br> in cells
  names(combined) <- br(names(combined))           # ... and in headers
  
  ## row ranges for each organism group
  sizes  <- vapply(mats, nrow, integer(1))
  ends   <- cumsum(sizes)
  starts <- ends - sizes + 1
  
  kb <- kableExtra::kbl(combined, escape = FALSE, align = "l",
                        caption = caption) |>
    kableExtra::kable_styling(full_width = FALSE,
                              bootstrap_options = c("bordered", "condensed"),
                              html_font = "monospace", font_size = font_size)
  
  ## width control: set an explicit CSS width per column if requested
  if (!is.null(col_widths))
    for (cc in seq_along(col_widths))
      kb <- kableExtra::column_spec(kb, cc, width = col_widths[cc])
  
  for (g in seq_along(mats))
    kb <- kableExtra::pack_rows(kb, names(mats)[g], starts[g], ends[g],
                                background = "#eef2f5")
  
  if (!is.null(file)) kableExtra::save_kable(kb, file, ...)  # .html / .png / .pdf
  kb
}


## ------------------------------------------------------------
## Save the COMBINED table as a flat CSV.
## CSV can't keep the visual rendering, so the organism grouping
## becomes an `Organism` column and each cell's line breaks are
## flattened with `flatten`. `mats` = same named list of matrices
## (built with identical column `display`) used for rendering.
## ------------------------------------------------------------
save_combined_matrix_csv <- function(mats, file, flatten = " | ") {
  dfs <- lapply(seq_along(mats), function(g) {
    M    <- mats[[g]]
    # flat <- apply(M, c(1, 2), function(x) gsub("\n", flatten, x, fixed = TRUE))
    cbind(Organism   = names(mats)[g],
          Comparison = rownames(M),
          as.data.frame(M, check.names = FALSE, stringsAsFactors = FALSE),
          row.names = NULL)
  })
  combined        <- do.call(rbind, dfs)
  names(combined) <- gsub("\n", flatten, names(combined), fixed = TRUE)
  write.csv(combined, file, row.names = FALSE)
  invisible(combined)
}