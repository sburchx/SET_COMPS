
### HELPER FUNCTIONS ----

# updated SB 3/26/25
make_subdir <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE)
    message(paste("📁 Folder created:", path))
  } else {
    message(paste("📁 Folder already exists:", path))
  }
}

save_csv <- function(data, path, row.names = FALSE, make_dir = TRUE) {
  dir <- dirname(path)
  
  if (make_dir && !dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
    message(paste("📁 Created directory:", dir))
  }
  
  write.csv(data, file = path, row.names = row.names)
  cat("✅ CSV saved to:", normalizePath(path), "\n")
}

## ---
## *DESeq2 FUNCTIONS* ----
## ---

# updated 3/26/25 SAB 
salmon_import <- function(star_salmon_path) {
  # Ensure required packages are installed and loaded
  pkgs <- c("readr", "stringr", "tximport")
  lapply(pkgs, function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg)
    }
    library(pkg, character.only = TRUE)
  })
  
  ##  Read tx2gene file
  tx2gene <- readr::read_tsv(file.path(star_salmon_path, "tx2gene.tsv"), col_names = FALSE) %>%
    as.data.frame()
  
  ##  Get list of quant.sf files
  files <- list.files(
    path = star_salmon_path,
    pattern = "quant.sf",
    recursive = TRUE,
    full.names = TRUE
  )

  # Clean names for the files (remove full path and "quant.sf")
  names(files) <- stringr::str_remove_all(files, paste(c(star_salmon_path, "/", "quant.sf"), collapse = "|"))
  
  ##  Import quantification data
  txi <- tximport::tximport(files, type = "salmon", tx2gene = tx2gene)
  
  return(c(txi, list(tx2gene = tx2gene)))
}

# updated 3/26/25 SAB
extract_samples <- function(files){
  ## capture sample in file name  or path
  x <- lapply( strsplit(files, "/"), rev)
  samples <- sapply(x, "[", 1 )
  # remove file extension
  samples <- gsub("\\..*", "", samples)
  if( any(duplicated(samples )) ){
    # use parent directory, 13555X2/Log.final.out
    samples <- sapply(x, "[", 2 )
  }
  if( any(duplicated(samples )) ){
    samples <- sapply(x, "[", 3 )
  }
  if( any(duplicated(samples )) ){
    stop("Sample names are not unique:  \n  ", paste(files, collapse="\n  "), call.=FALSE)
  }
  samples
}


## ---
## *PCA FUNCTIONS* ----
## ---

# helper to get ID for PCA
# 2026-03-18 12:30 SAB ---
get_ids <- function(meta, timepoints = NULL, treatments = NULL,
                    genotypes = NULL, pretreat_days = NULL, baseline_treatment = "CTRL") {
  
  m <- meta
  
  if (!is.null(timepoints))  m <- dplyr::filter(m, timepoint %in% timepoints)
  if (!is.null(treatments))  m <- dplyr::filter(m, treatment %in% treatments)
  if (!is.null(genotypes))   m <- dplyr::filter(m, genotype %in% genotypes)
  
  ids_main <- m$sampleID
  
  if (!is.null(pretreat_days)) {
    m_pre <- meta
    m_pre <- dplyr::filter(m_pre, timepoint %in% pretreat_days)
    m_pre <- dplyr::filter(m_pre, treatment %in% baseline_treatment)
    if (!is.null(genotypes)) m_pre <- dplyr::filter(m_pre, genotype %in% genotypes)
    
    ids_main <- unique(c(ids_main, m_pre$sampleID))
  }
  
  unique(ids_main)
}


# 2026-06-07 17:54 SAB ---
generatePCA <- function(vsd, metadata, sampleIDs, colby, shape, label, save_prefix,
                        removeVar = 0.9, flip_pc2 = FALSE, annotation_colors = NULL,
                        labSize = 3, pointSize = 3, legendPosition = "right",
                        save_inputs = FALSE, quiet = TRUE, show_caption = TRUE,
                        scale_to_axes = TRUE, pca_base = 10,
                        base_width = 10, min_height = 6,
                        fig_width = NULL, fig_height = NULL, fig_units = "mm") {
  
  # enforce order + alignment 
  sampleIDs <- intersect(sampleIDs, colnames(vsd))
  if (length(sampleIDs) == 0) stop("❌ No sampleIDs overlap with vsd colnames.")
  
  # metadata rownames = sampleID
  missing_meta <- setdiff(sampleIDs, rownames(metadata))
  if (length(missing_meta) > 0) {
    stop("❌ These sampleIDs are missing from metadata rownames:\n", paste(missing_meta, collapse = ", "))
  }
  
  meta_subset <- metadata[sampleIDs, , drop = FALSE]
  vsd_subset  <- vsd[, sampleIDs, drop = FALSE]
  
  stopifnot(identical(colnames(vsd_subset), rownames(meta_subset)))
  
  # label handling
  lab <- if (!is.na(label) && !is.null(label)) meta_subset[[label]] else NULL
  
  # build caption on bottom of plot 
  plot_caption <- stringr::str_wrap(
    paste(tail(strsplit(save_prefix, "/")[[1]], 4), collapse = "/")
  )
  
  # Generate PCA
  pca_result <- if (quiet) {
    suppressMessages(PCAtools::pca(vsd_subset, metadata = meta_subset, removeVar = removeVar))
  } else {
    PCAtools::pca(vsd_subset, metadata = meta_subset, removeVar = removeVar)
  }
  
  if (flip_pc2) {
    pca_result$rotated[, 2]  <- pca_result$rotated[, 2]  * -1
    pca_result$loadings[, 2] <- pca_result$loadings[, 2] * -1
  }
  
  # --- build PCA plot using biplot --- 
  pca_plot <- PCAtools::biplot(
    pca_result,
    colby = colby,
    shape = shape,
    lab = lab,
    labSize = labSize,
    pointSize = pointSize,
    legendPosition = legendPosition,
    returnPlot = TRUE
  ) 
  
  if (show_caption) {
    pca_plot <- pca_plot +
      ggplot2::labs(caption = plot_caption)
  }
  
  # standardize axis and font
  pca_plot <- pca_plot +
    ggplot2::theme(
      text = element_text(
        family = "helvetica_neue", colour = "black", size = pca_base),
      
      panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.4),
      axis.line = element_blank(), 
      
      plot.caption = if (show_caption) {
        element_text(size = pca_base * 0.5)
      } else {
        element_blank()
      },
      
      plot.margin = margin(8, 10, 10, 18),
      panel.grid.major = element_line(colour = "#e6e6e6",linewidth = 0.2),
      panel.grid.minor = element_line(colour = "#e6e6e6",linewidth = 0.1),
      
      axis.title.x = element_text(vjust = -0.2, size = pca_base * 1.2),
      axis.title.y = element_text(angle = 90, vjust = 2, size = pca_base* 1.2),
      axis.text.x = element_text(size = pca_base, colour = "black", margin = margin(t = 2)),
      axis.text.y = element_text(size = pca_base, colour = "black", margin = margin(r = 2)),
      axis.ticks = element_line(colour = "black",linewidth = 0.4),
      axis.ticks.length = unit(0.1, "cm"),
      
      legend.text = element_text(size = pca_base),
      legend.key.size = unit(0.25, "cm"),
      legend.box.spacing = unit(0.03, "cm"),
      
      strip.text = element_text(size = pca_base * 0.9),
      
      plot.title = element_text(size = pca_base * 1.2, face = "bold", hjust = 0.5)
    ) +
    guides(
      colour = guide_legend(override.aes = list(size = pointSize * 0.8), title = NULL),
      shape = guide_legend(override.aes = list(size = pointSize * 0.8), title = NULL)
    )
      
  # palette override
  if (!is.null(annotation_colors) && !is.null(annotation_colors[[colby]])) {
    pca_plot <- pca_plot + ggplot2::scale_colour_manual(
      values = annotation_colors[[colby]],
      breaks = names(annotation_colors[[colby]]),
      drop = FALSE
    )
  }
  
  # figure-sized export overrides auto sizing
  if (!is.null(fig_width) && !is.null(fig_height)) {
    
    plot_width  <- fig_width
    plot_height <- fig_height
    plot_units  <- fig_units
    
  } else if (scale_to_axes) {
    # auto-scale plot dimensions based on PCA axis ranges
    pc1_range <- diff(range(pca_result$rotated[, 1], na.rm = TRUE))
    pc2_range <- diff(range(pca_result$rotated[, 2], na.rm = TRUE))
    aspect_ratio <- pc2_range / pc1_range
    
    plot_width  <- base_width
    plot_height <- max(min_height, base_width * aspect_ratio)
    plot_units <- "in"
    
  } else {
    
    plot_width  <- 12
    plot_height <- 10
    plot_units  <- "in"
  }
  
  ggplot2::ggsave(
    paste0(save_prefix, "_PCA.pdf"),
    pca_plot,
    width = plot_width,
    height = plot_height,
    units = plot_units,
    dpi = 600,
    bg = "white"
  )
  
  
  # optional: save reproducibility inputs for figure-grade tracking
  if (save_inputs) {
    pca_inputs <- list(
      exp_name = if ("exp_name" %in% names(parent.frame())) get("exp_name", parent.frame()) else NA,
      save_prefix = save_prefix,
      sampleIDs = sampleIDs,
      removeVar = removeVar,
      flip_pc2 = flip_pc2,
      colby = colby,
      shape = shape,
      label = label,
      vsd_dim = dim(vsd_subset),
      timestamp = as.character(Sys.time())
    )
    saveRDS(pca_inputs, paste0(save_prefix, "_PCA_inputs.rds"))
  }
  
  invisible(NULL)
}

# 2026-03-03 17:40 SAB ---
pcaLoadings <- function(vsd_data, metadata, sample_ids, exp_name, name, results_dir,
                        removeVar = 0.9, flip_pc2 = FALSE, quiet = TRUE, scale_to_axes = TRUE,
                        base_width = 10, min_height = 6, max_height = 12) {
  
  sample_ids <- intersect(sample_ids, colnames(vsd_data))
  if (length(sample_ids) == 0) stop("❌ No sample_ids overlap with vsd_data colnames.")
  
  missing_meta <- setdiff(sample_ids, rownames(metadata))
  if (length(missing_meta) > 0) {
    stop("❌ These sample_ids are missing from metadata rownames:\n", paste(missing_meta, collapse = ", "))
  }
  
  metadata_subset <- metadata[sample_ids, , drop = FALSE]
  vsd_subset <- vsd_data[, sample_ids, drop = FALSE]
  stopifnot(identical(colnames(vsd_subset), rownames(metadata_subset)))
  
  pca_result <- if (quiet) {
    suppressMessages(PCAtools::pca(vsd_subset, metadata = metadata_subset, removeVar = removeVar))
  } else {
    PCAtools::pca(vsd_subset, metadata = metadata_subset, removeVar = removeVar)
  }
  
  if (flip_pc2) {
    pca_result$rotated[, 2]  <- pca_result$rotated[, 2]  * -1
    pca_result$loadings[, 2] <- pca_result$loadings[, 2] * -1
  }
  
  screeplot_path <- file.path(results_dir, paste0("screeplot_", name, "_", exp_name, ".pdf"))
  loadings_path  <- file.path(results_dir, paste0("loadings_", name, "_", exp_name, ".pdf"))
  
  scree <- PCAtools::screeplot(
    pca_result,
    components = PCAtools::getComponents(pca_result, 1:10),
    axisLabSize = 18,
    titleLabSize = 22
  )
  ggplot2::ggsave(screeplot_path, plot = scree, width = 12, height = 6, dpi = 600)
  
  p_load <- PCAtools::biplot(
    pca_result,
    showLoadings = TRUE,
    showLoadingsNames = TRUE,
    ntopLoadings = 10,
    fillBoxedLoadings = scales::alpha("white", 1/4),
    pointSize = 2,
    lab = NULL,
    drawConnectors = FALSE,
    legendPosition = "none",
    title = paste("Top 10 PCA Loadings:", name)
  ) +
    ggplot2::coord_fixed()
  
  # calculate proportional dimensions based on PC1/PC2 ranges
  if (scale_to_axes) {
    pc1_range <- diff(range(pca_result$rotated[, 1], na.rm = TRUE))
    pc2_range <- diff(range(pca_result$rotated[, 2], na.rm = TRUE))
    
    aspect_ratio <- pc2_range / pc1_range
    
    plot_width  <- base_width
    plot_height <- base_width * aspect_ratio
    
    plot_height <- max(min_height, plot_height)
  } else {
    plot_width  <- 12
    plot_height <- 10
  }
  
  ggplot2::ggsave(
    loadings_path,
    plot = p_load,
    width = plot_width,
    height = plot_height,
    dpi = 600
  )
  
  invisible(NULL)
}

## ---
## *CHROMOSOME MAPPING* ----
## ---

# -- Summarize chromosome counts for every comp/day as a ratio DEG/counts
summarize_chr_counts <- function(chromosome_results, gene_chrom, all_tested) {
  standard_chrs <- c(paste0("chr", 1:22), "chrX", "chrY")
  regs          <- c("Upregulated", "Downregulated")
  
  # tested counts per chromosome (denominator) + complete list
  tested_df <- tibble(gene_id = as.character(all_tested)) %>%
    left_join(gene_chrom, by = "gene_id") %>%
    dplyr::count(chrom, name = "Tested_Count") %>%
    tidyr::complete(chrom = standard_chrs, fill = list(Tested_Count = 0L)) %>%
    dplyr::rename(Chromosome = chrom) %>%
    dplyr::mutate(Chromosome = factor(Chromosome, levels = standard_chrs))
  
  #DEGS per chrom (numerator)
  chr_counts_df <- bind_rows(
    lapply(names(chromosome_results), function(main_comp) { #loop over comps 
      day_list <- chromosome_results[[main_comp]]
      bind_rows(
        lapply(names(day_list), function(day_name) { # loop over each day 
          df <- day_list[[day_name]]
          regulation <- ifelse(grepl("\\.up$", day_name), "Upregulated", "Downregulated")
          
          df %>%
            mutate(
              Main_Comparison = main_comp,
              Comparison      = day_name,
              Regulation      = regulation
            ) %>%
            dplyr::count(Main_Comparison, Comparison, Regulation, chrom, name = "Gene_Count")
        })
      )
    })
  )
  
  # zero-fill missing chr×reg per panel, add day, join denominator, compute ratio 
  chr_counts_df <- chr_counts_df %>%
    dplyr::mutate(
      Day        = as.integer(sub(".*_D(\\d+)_.*", "\\1", Comparison)),
      Chromosome = chrom,
      Regulation = factor(Regulation, levels = regs)
    ) %>%
    dplyr::group_by(Main_Comparison, Comparison, Day) %>%
    tidyr::complete(
      Chromosome = standard_chrs,
      Regulation = regs,
      fill       = list(Gene_Count = 0L)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(Chromosome = factor(Chromosome, levels = standard_chrs)) %>%
    dplyr::left_join(tested_df, by = "Chromosome") %>%
    dplyr::mutate(
      Tested_Count = tidyr::replace_na(Tested_Count, 0L),
      Ratio        = ifelse(Tested_Count > 0, Gene_Count / Tested_Count, NA_real_)
    )
  
  chr_counts_df
}

# -- for plot labels
prettify_chr_counts <- function(chr_counts_df) {
  chr_counts_df %>%
    dplyr::mutate(
      Comparison_label = gsub("\\.up$|\\.dn$", "", Comparison),
      Day_for_order = as.integer(sub(".*_D(\\d+)_.*", "\\1", Comparison_label))
    ) %>%
    dplyr::group_by(Main_Comparison) %>%
    dplyr::mutate(
      Comparison_label = factor(
        Comparison_label,
        levels = unique(Comparison_label[order(Day_for_order)])
      )
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(-Day_for_order)
}


plot_chr_updown <- function(chr_counts_df, save_dir,
                            mode = c("facet", "individual"),
                            use_ratio = TRUE,
                            ratio_limits = c(0, 0.5),
                            counts_limits = NULL) {
  
  mode <- match.arg(mode)
  standard_chrs <- c(paste0("chr", 1:22), "chrX", "chrY")
  yvar <- if (use_ratio) "Ratio" else "Gene_Count"
  ylab <- if (use_ratio) "Ratio of DEGs" else "Number of DEGs"
  
  yscale <- if (use_ratio) {
    if (is.null(ratio_limits)) {
      scale_y_continuous(expand = expansion(mult = c(0, 0.01)))
    } else {
      scale_y_continuous(limits = ratio_limits, expand = expansion(mult = c(0, 0)))
    }
  } else {
    if (is.null(counts_limits)) {
      scale_y_continuous(expand = expansion(mult = c(0, 0.02)))
    } else {
      scale_y_continuous(limits = counts_limits, expand = expansion(mult = c(0, 0)))
    }
  }
  
  plots <- list()
  
  for (cmp in unique(chr_counts_df$Main_Comparison)) {
    
    subdf <- chr_counts_df %>%
      dplyr::filter(Main_Comparison == cmp, Chromosome %in% standard_chrs) %>%
      dplyr::mutate(Chromosome = factor(Chromosome, levels = standard_chrs))
    
    if (mode == "facet") {
      
      p <- ggplot(subdf, aes(x = Chromosome, y = .data[[yvar]], fill = Regulation)) +
        geom_col(position = position_dodge(width = 0.8), width = 0.7) +
        scale_fill_manual(values = c("Upregulated" = "red3",
                                     "Downregulated" = "#1a95eeff")) +
        facet_wrap(~ Comparison_label, nrow = 1) +
        yscale +
        labs(
          title = paste0("DEG Chromosome Distribution — ", cmp),
          x = NULL, y = ylab, fill = ""
        ) +
        theme_pub() +
        theme(
          axis.text.x = element_text(angle = 90, hjust = 1, size = 8, vjust = 0.5),
          strip.text = element_text(size = 8, face = "bold"),
          axis.text.y = element_text(size = 8),
          axis.title.y = element_text(size = 10),
          panel.spacing = unit(1, "lines")
        )
      
      key <- paste0(cmp, "_chrom_by_day_", if (use_ratio) "ratio" else "counts")
      plots[[key]] <- p
      
      ggsave(file.path(save_dir, paste0(key, ".pdf")), p, width = 14, height = 3)
    }
    
    if (mode == "individual") {
      
      panels <- unique(as.character(subdf$Comparison_label))
      
      for (panel in panels) {
        
        df_d <- subdf %>%
          dplyr::filter(as.character(Comparison_label) == panel)
        
        p <- ggplot(df_d, aes(x = Chromosome, y = .data[[yvar]], fill = Regulation)) +
          geom_col(position = position_dodge(width = 0.8), width = 0.7) +
          scale_fill_manual(values = c("Upregulated" = "red3",
                                       "Downregulated" = "#1a95eeff")) +
          yscale +
          labs(
            title = paste0("DEG Chromosome Distribution — ", panel),
            x = NULL, y = ylab, fill = ""
          ) +
          theme_pub() +
          theme(
            legend.position = "top",
            axis.text.x = element_text(angle = 90, hjust = 1, size = 8, vjust = 0.5),
            axis.text.y = element_text(size = 8),
            axis.title.y = element_text(size = 10),
            plot.title = element_text(size = 11, face = "plain")
          )
        
        key <- paste0(panel, "_updown_", if (use_ratio) "ratio" else "counts")
        plots[[key]] <- p
        
        ggsave(file.path(save_dir, paste0(key, ".pdf")), p, width = 10, height = 3)
      }
    }
  }
  
  return(plots)
}

summarize_total_ratio <- function(chr_counts_df) {
  
  chr_counts_df %>%
    dplyr::group_by(Main_Comparison, Chromosome, Day) %>%
    dplyr::summarize(
      Total_DEGs   = sum(Gene_Count, na.rm = TRUE),
      Tested_Count = dplyr::first(Tested_Count),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      Total_Ratio = dplyr::if_else(
        Tested_Count > 0,
        Total_DEGs / Tested_Count,
        NA_real_
      )
    )
}

plot_total_ratio_by_chrom_and_day <- function(total_df, save_dir,
                                              ratio_limits = c(0, 0.5),
                                              annotation_colors = NULL) {
  standard_chrs <- c(paste0("chr", 1:22), "chrX", "chrY")
  comps <- unique(total_df$Main_Comparison)
  
  # y scale
  yscale <- if (is.null(ratio_limits)) {
    scale_y_continuous(expand = expansion(mult = c(0, 0.01)))
  } else {
    scale_y_continuous(limits = ratio_limits, expand = expansion(mult = c(0, 0)))
  }
  
  plots <- list()
  
  for (cmp in comps) {
    subdf <- total_df %>%
      dplyr::filter(
        Main_Comparison == cmp,
        Chromosome %in% standard_chrs
      ) %>%
      dplyr::mutate(
        Chromosome = factor(Chromosome, levels = standard_chrs),
        Day = factor(Day, levels = sort(unique(Day)))
      )
    
    #  palette for days
    pal_days <- NULL
    if (!is.null(annotation_colors) && "timepoint" %in% names(annotation_colors)) {
      pal_days <- annotation_colors$timepoint
      names(pal_days) <- sub("^D", "", names(pal_days))  # "D3" → "3"
    }
    
    p <- ggplot(subdf, aes(x = Chromosome, y = Total_Ratio, fill = Day)) +
      geom_col(position = position_dodge(width = 0.8), width = 0.7) +
      yscale +
      labs(
        title = paste0("DEG Ratio per Chromosome by Day — ", cmp),
        x = NULL, y = "Ratio of DEGs", fill = "Day"
      ) +
      theme_pub() +
      theme(
        axis.text.x  = element_text(angle = 90, hjust = 1, size = 8, vjust = 0.5),
        axis.text.y  = element_text(size = 8),
        axis.title.y = element_text(size = 9),
        strip.background = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background  = element_rect(fill = "transparent", color = NA)
      )
    
    # apply custom colors
    if (!is.null(pal_days)) {
      p <- p + scale_fill_manual(values = pal_days)
    }
    
    key <- paste0(cmp, "_total_ratio_by_day")
    plots[[key]] <- p 
    
    ggsave(file.path(save_dir, paste0(key, ".pdf")),
           p, width = 10, height = 4)
  }
  return(plots)
}

make_sex_qc <- function(vsd_mat, metadata, line_name = exp_name, y_threshold = 5) {
  
  y_markers <- c("RPS4Y1", "DDX3Y", "KDM5D", "EIF1AY", "UTY", "ZFY")
  y_markers <- intersect(y_markers, rownames(vsd_mat))
  
  sex_scores <- tibble::tibble(
    sampleID = colnames(vsd_mat),
    Y_score = colMeans(vsd_mat[y_markers, , drop = FALSE])
  ) %>%
    dplyr::left_join(metadata, by = "sampleID")
  
  line_score <- mean(sex_scores$Y_score, na.rm = TRUE)
  predicted_sex <- ifelse(line_score >= y_threshold, "XY-like / male", "XX-like / female")
  
  p <- ggplot(sex_scores, aes(x = line_name, y = Y_score)) +
    geom_jitter(aes(color = genotype, shape = treatment),
                width = 0.08, size = 3) +
    geom_hline(yintercept = y_threshold, linetype = "dashed") +
    annotate("text", x = 1, y = line_score, 
             label = predicted_sex, vjust = -1, fontface = "bold") +
    theme_pub() +
    labs(
      title = paste0("Sex QC: ", line_name),
      subtitle = paste0("Mean Y-marker score = ", round(line_score, 2)),
      x = NULL,
      y = "Mean chrY marker expression"
    )
  
  p <- p +
    theme(
      legend.position = "top",
      legend.direction = "horizontal",
      legend.box = "horizontal",
      legend.title = element_blank(),
      legend.text = element_text(size = 9),
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 10),
      axis.title = element_text(size = 11),
      axis.text = element_text(size = 9),
      plot.margin = margin(5, 5, 5, 5)
    )
  
  list(
    plot = p,
    scores = sex_scores,
    line_score = line_score,
    predicted_sex = predicted_sex
  )
}


## ---
## *VOLCANO PLOT* ----
## ---

# 2026-01-21 11:39 SAB ---
createVolcano <- function(geneList, log2FC_Cutoff, pval_cutoff, labelGenes, 
                          plotTitle, file_path, ylim_max = NULL, cap = 60, 
                          labelGenes_override = NULL) {
  

  df <- geneList
  
  # -- compatibility layer --
  if (!"gene_id" %in% colnames(df) && "gene" %in% colnames(df)) {
    df$gene_id <- df$gene
  }
  
  if (!"log2FoldChange_raw" %in% colnames(df) && "log2FoldChange" %in% colnames(df)) {
    df$log2FoldChange_raw <- df$log2FoldChange
  }
  
  if (!"log2FoldChange_shrunk" %in% colnames(df) && "log2FoldChange" %in% colnames(df)) {
    df$log2FoldChange_shrunk <- df$log2FoldChange
  }
  
  # basic checks
  req_cols <- c("gene_id", "padj", "log2FoldChange", "baseMean", "pvalue")
  missing_cols <- setdiff(req_cols, colnames(df))
  if (length(missing_cols)) {
    stop("createVolcano(): missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  # LFC & padj transforms
  df$log2FC <- df$log2FoldChange_shrunk
  df$logpadj <- -log10(df$padj)
  
  label_set <- if (!is.null(labelGenes_override)) {
    labelGenes_override
  } else {
    labelGenes
  }
  
  df$panel <- ifelse(df$gene_id %in% label_set, "label", "no label")
  
  # cap the padj axis and mark capped points
  df$logpadj_capped <- pmin(df$logpadj, cap)
  df$capped         <- df$logpadj > cap   # TRUE if we hit the cap
  
  # classify DEGs on raw LFC and padj
  df$differentialExpression <- case_when(
    df$log2FoldChange_raw >  log2FC_Cutoff & df$padj < pval_cutoff ~ "UP",
    df$log2FoldChange_raw < -log2FC_Cutoff & df$padj < pval_cutoff ~ "DOWN",
    TRUE ~ "NO"
  )
  
  #Axis limits
  x_max    <- 10
  xlim_range <- c(- x_max,  x_max)
  
  y_max    <- cap
  ylim_range <- c(0, y_max * 1.1)
  
  # DEG counts (based on classification from raw log2FC)
  up_count <- sum(df$differentialExpression == "UP", na.rm = TRUE)
  down_count <- sum(df$differentialExpression == "DOWN", na.rm = TRUE)
  
  # Parse title
  title_parts <- strsplit(plotTitle, "_vs_")[[1]]
  group1 <- title_parts[1]
  group2 <- title_parts[2]
  
  folder_short <- paste(tail(strsplit(dirname(file_path), "/")[[1]], 5), collapse = "/")
  
  
  volcano <- ggplot(df, aes(x = log2FC, y = logpadj_capped, color = differentialExpression)) +
    ggrastr::geom_point_rast(aes(shape = capped), size = 1.5, alpha = 1) +
    scale_shape_manual(values = c(`FALSE` = 16,  # circle
                                  `TRUE`  = 16), # triangle = 17
                       guide  = "none" ) +
    scale_color_manual(values = c("DOWN" = "#1a95eeff", "NO" = "grey", "UP" = "red3")) +
    geom_hline(yintercept = -log10(pval_cutoff), linetype = "dashed", color = "black", linewidth = 0.75) +
    geom_vline(xintercept = c(-log2FC_Cutoff, log2FC_Cutoff), linetype = "dashed", color = "black", linewidth = 0.75) +
    coord_cartesian(xlim = xlim_range, ylim = ylim_range, clip = "off") +
    scale_y_continuous( limits = ylim_range, breaks = scales::pretty_breaks(n = 5)) + #removed trans = "log1p",
    labs(
      x = "log2FoldChange",
      y = expression(-log[10]~"(FDR-adjusted p-value)"),
      title = bquote(bold(.(group1)) ~ "vs." ~ bold(.(group2)))
    ) +
    geom_label_repel(
      aes(label = ifelse(panel == "label" & differentialExpression != "NO", gene_id, NA)),
      box.padding = 0.25, point.padding = 0.5,
      segment.color = "gray30", fill = "white", alpha = 0.9,
      force_pull = 100, max.overlaps = Inf,
      size = 3, label.size = 0.25, na.rm = TRUE
    )  +
    annotate("text", x = -x_max * 0.8, y = y_max * 1.1, label = paste("Down:", down_count), size = 6, color = "#1a95eeff", fontface = "bold") +
    annotate("text", x = x_max * 0.8, y = y_max * 1.1, label = paste("Up:", up_count), size = 6, color = "red3", fontface = "bold") +
    annotate(
      "label",
      x = x_max + 0.5, y = y_max * 1.1,
      label = paste0("log2FoldChange > |", log2FC_Cutoff, "|\npadj < ", pval_cutoff),
      size = 3.5, color = "black", fill = "white", label.size = 0.5, hjust = 1,  vjust = -0.4 
    ) +
    theme_pub(base_size = 14) +
    theme(
      aspect.ratio = 1, #1 is square
      legend.position = "none",
      plot.margin = margin(t = 40,  b = 20),
      plot.title = element_text(margin = margin(b = 30))
    )
  
  volcano <- volcano +
    labs(
      caption = paste0(
        "DEGs defined using *raw* log2FoldChange (|log2FC| > ", log2FC_Cutoff,
        ", FDR < ", pval_cutoff, "); shrunken log2FoldChange used for plotting.\n",
        "log2FoldChange capped at ±10 and –log₁₀(padj) capped at ", cap, "\n", #(triangles = true –log₁₀>cap).
        "Source: ", paste(folder_short, collapse = "/")
      )
    ) 
  
  return(volcano)
}

## ---
## *GENE ONTOLOGY ENRICHMENT* WITH enrichGO() ----
## ---

# 2026-05-30 18:56 SAB ---
run_go_terms <- function(gene_list, all_genes, file_suffix, output_directory, organism, 
                         go_pval_cutoff, go_qval_cutoff) {
  
  org_pkg <- deparse(substitute(organism))
  if (!requireNamespace(org_pkg, quietly = TRUE)) {
    stop("Required annotation package ", org_pkg, " is not installed.")
  }
  
  gene_ids <- bitr(gene_list, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = organism) %>%
    distinct(SYMBOL, .keep_all = TRUE)
  
  background_ids <- bitr(all_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = organism) %>%
    distinct(SYMBOL, .keep_all = TRUE)
  
  unmapped <- setdiff(gene_list, gene_ids$SYMBOL)
  
  if (length(unmapped) > 0) {
    fail_pct <- round(length(unmapped) / length(gene_list) * 100, 2)
    message(sprintf(" %.2f%% of genes failed to map. First 10: %s",
                    fail_pct, paste(head(unmapped, 10), collapse = ", ")))
    
    unmapped_dir <- file.path(output_directory, "Unmapped_Genes")
    if (!dir.exists(unmapped_dir)) dir.create(unmapped_dir, recursive = TRUE)
    
    unmapped_path <- file.path(unmapped_dir, paste0("unmappedgenes_", file_suffix, ".csv"))
    write_csv(data.frame(Failed_Genes = unmapped), file = unmapped_path)
    message("Unmapped gene list saved to: ", normalizePath(unmapped_path))
  } else {
    message("All genes successfully mapped.")
  }
  
  go_enrich <- enrichGO(
    gene = gene_ids$ENTREZID,
    OrgDb = organism,
    keyType = "ENTREZID",
    ont = "BP",
    readable = TRUE,
    pvalueCutoff = go_pval_cutoff,
    qvalueCutoff = go_qval_cutoff,
    pAdjustMethod = "BH",
    universe = background_ids$ENTREZID
  )
  
  # handle NULL or empty enrichment 
  if (is.null(go_enrich) || nrow(as.data.frame(go_enrich)) == 0) {
    message("No GO terms enriched for ", file_suffix, " at the specified cutoffs.")
    # Save an empty CSV for bookkeeping
    go_enrich_df <- data.frame()
    output_path <- file.path(output_directory, paste0("GO_", file_suffix, ".csv"))
    write_csv(go_enrich_df, file = output_path)
    return(list(
      go       = go_enrich_df,
      unmapped = unmapped
    ))
  }
  
  
  go_enrich <- simplify(go_enrich, cutoff = 0.7, by = "p.adjust", select_fun = min)
  go_enrich_df <- as.data.frame(go_enrich)
  
  output_path <- file.path(output_directory, paste0("GO_", file_suffix, ".csv"))
  write_csv(go_enrich_df, file = output_path)
  message("GO results saved: ", normalizePath(output_path))
  
  return(list(
    go = go_enrich_df,
    unmapped = unmapped
  ))
}


# updated 8/27/25 SAB
bubbleplotGO <- function(upGO, downGO, n_terms = 10, source_path) {
  
  make_bubbleplot <- function(go_df, color, direction_label, flip_y = FALSE) {
    if (is.null(go_df) || nrow(go_df) == 0) return(NULL)
    
    # Clean and score
    go_df <- go_df[!is.na(go_df$p.adjust), ]
    go_df <- go_df[1:min(n_terms, nrow(go_df)), ]
    go_df$qscore <- -log10(go_df$p.adjust)
    go_df$Description <- as.factor(go_df$Description)
    go_df$title <- paste("GO BP:", direction_label)
    
    source_label <- paste("Source:", paste(tail(strsplit(source_path, "/")[[1]], 5), collapse = "/"))
    
    # Base plot: Description on x, qscore on y
    p <- ggplot(go_df, aes(x = forcats::fct_reorder(Description, qscore), y = qscore, size = Count)) +
      geom_point(color = color, alpha = 0.8) +
      coord_flip() +
      scale_size_continuous(name = "Gene count", range = c(4, 10)) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
      scale_x_discrete(
        position = ifelse(direction_label == "Up", "top", "bottom"),
        labels = scales::wrap_format(35),
        expand = expansion(mult = c(0.075, 0.075))
      ) +
      theme_pub() +
      theme(
        legend.position = if (direction_label == "Up") "right" else "left",
        legend.direction = "vertical",
        legend.justification = "center",
        panel.grid.major = element_line(colour = "#f0f0f0"),
        axis.text.y = element_text(size = 16, margin = margin(l = 8)),
        axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_blank(),
        strip.text.x = element_text(size = 18, face = "bold", color = "black"),
        plot.margin = margin(t = 20, r = 40, b = 20, l = 40),
        plot.caption = element_text(size = 8, color = "gray30", hjust = 0, vjust = 0, margin = margin(t = 12))
      ) +
      labs(
        caption = source_label,
        x = NULL,
        y = expression(-log[10]~"(FDR-adjusted p-value)")
      ) +
      facet_grid(. ~ title)
    
    # Flip y-axis for upregulated so bubbles go right
    if (flip_y) {
      p <- p + scale_y_reverse(expand = expansion(mult = c(0, 0.1)))
    }
    
    return(p)
  }
  
  list(
    upbarplot = make_bubbleplot(upGO, color = "red3", direction_label = "Up", flip_y = FALSE),
    downbarplot = make_bubbleplot(downGO, color = "#1a95eeff", direction_label = "Down", flip_y = TRUE)
  )
}


# 2026-03-18 14:44 SAB ---
plot_single_go_bubble <- function(go_df, title_text, source_path, n_terms = 10) {
  
  if (is.null(go_df) || nrow(go_df) == 0) return(NULL)
  
  go_df <- go_df[!is.na(go_df$p.adjust), ]
  go_df <- go_df[1:min(n_terms, nrow(go_df)), ]
  
  # parse GeneRatio like "22/85" -> numeric
  go_df$GeneRatio_num <- vapply(go_df$GeneRatio, function(x) {
    parts <- strsplit(x, "/")[[1]]
    as.numeric(parts[1]) / as.numeric(parts[2])
  }, numeric(1))
  
  # wrap labels first, then order by GeneRatio
  go_df$Description_wrapped <- scales::wrap_format(40)(go_df$Description)
  go_df$Description_wrapped <- forcats::fct_reorder(go_df$Description_wrapped, go_df$GeneRatio_num)
  
  source_label <- paste("Source:", paste(tail(strsplit(source_path, "/")[[1]], 5), collapse = "/"))
  
  ggplot(go_df, aes(x = GeneRatio_num, y = Description_wrapped, size = Count, fill = p.adjust)) +
    geom_point(color = "black", alpha = 0.9, stroke = 0.6, shape = 21) +
    scale_size_continuous(name = "Count", range = c(4, 10)) +
    scale_fill_gradient(
      low = "#d6604d", high = "#4393c3",
      trans = "reverse", name = "FDR"
    ) +
    scale_x_continuous(expand = expansion(mult = c(0.02, 0.08))) +
    labs(
      title = title_text,
      x = "Gene Ratio",
      y = NULL,
      caption = source_label
    ) +
    theme_pub() +
    theme(
      legend.position = "right",
      legend.direction = "vertical",
      legend.justification = "center",
      panel.grid.major = element_line(colour = "#f7f7f7"),
      axis.text.y = element_text(size = 12, margin = margin(r = 8)),
      axis.text.x = element_text(size = 11),
      axis.title.x = element_text(size = 13),
      axis.title.y = element_blank(),
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      axis.line = element_line(size = 0.5),
      plot.margin = margin(t = 20, r = 20, b = 20, l = 20),
      plot.caption = element_text(size = 8, color = "gray30", hjust = 0, margin = margin(t = 10)),
      legend.box = "vertical",
      legend.spacing.y = unit(0.3, "cm")
    )+ guides(
      fill = guide_colorbar(order = 1),
      size = guide_legend(order = 2)
    )
}

# helper to compare groups within a day
compare_go_terms_within_day <- function(go_df, day) {
  
  day_df <- go_df %>%
    dplyr::filter(timepoint == day) %>%
    dplyr::select(timepoint, response_group, Description, p.adjust, Count, GeneRatio)
  
  shared_terms <- day_df %>%
    dplyr::filter(response_group == "shared") %>%
    dplyr::pull(Description) %>%
    unique()
  
  d21_terms <- day_df %>%
    dplyr::filter(response_group == "D21_only") %>%
    dplyr::pull(Description) %>%
    unique()
  
  t21_terms <- day_df %>%
    dplyr::filter(response_group == "T21_only") %>%
    dplyr::pull(Description) %>%
    unique()
  
  list(
    shared_only = setdiff(shared_terms, union(d21_terms, t21_terms)),
    D21_only_unique = setdiff(d21_terms, union(shared_terms, t21_terms)),
    T21_only_unique = setdiff(t21_terms, union(shared_terms, d21_terms)),
    shared_vs_D21_overlap = intersect(shared_terms, d21_terms),
    shared_vs_T21_overlap = intersect(shared_terms, t21_terms),
    D21_vs_T21_overlap = intersect(d21_terms, t21_terms)
  )
}

## ---
## *STITCH GO + VOLCANO PLOT* ----
## ---
# 8/27/25 SAB
stitch_combined_plot <- function(down_plot, volcano_plot, up_plot, width_ratios = c(1, 1.3, 1), caption = NULL) {
  
  if (is.null(volcano_plot)) return(NULL)
  
  # Create blank plot to fill if needed
  blank_msg <- ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = "No GO terms found", size = 6, hjust = 0.5, vjust = 0.5) +
    theme_void()
  
  # Replace NULLs with blank plots
  if (is.null(down_plot)) down_plot <- blank_msg
  if (is.null(up_plot)) up_plot <- blank_msg
  
  
  # Apply consistent margins
  down_plot <- down_plot + theme(plot.margin = margin(60, -5, 30, 20))
  volcano_plot <- volcano_plot + theme(plot.margin = margin(20, -5, 15, -5))
  up_plot <- up_plot + theme(plot.margin = margin(60, 20, 30, -5))
  
  # Stitch plots side by side
  combo <- cowplot::plot_grid(
    down_plot, volcano_plot, up_plot,
    nrow = 1, rel_widths = width_ratios, align = "h"
  )
  
  
  
  return(combo)
}

## ---
##  *MODULE FUNCTIONS* ----
## ---

score_modules <- function(vsd_mat, module_gene_sets) {
  
  expr <- as.matrix(vsd_mat)
  
  # gene-wise z-score across selected samples
  expr_z <- t(scale(t(expr)))
  expr_z[is.na(expr_z)] <- 0
  
  module_scores <- lapply(names(module_gene_sets), function(m) {
    genes_m <- intersect(module_gene_sets[[m]], rownames(expr_z))
    
    if (length(genes_m) == 0) {
      return(rep(NA_real_, ncol(expr_z)))
    }
    
    colMeans(expr_z[genes_m, , drop = FALSE], na.rm = TRUE)
  })
  
  module_scores <- do.call(rbind, module_scores)
  rownames(module_scores) <- names(module_gene_sets)
  colnames(module_scores) <- colnames(expr_z)
  
  return(module_scores)
}

sym2ent <- function(genes){
  genes <- unique(as.character(genes))
  genes <- genes[!is.na(genes) & genes != ""]
  
  eg <- clusterProfiler::bitr(
    genes,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
  )
  
  # keep only valid mappings & unique ENTREZ
  eg <- eg %>%
    dplyr::filter(!is.na(ENTREZID), ENTREZID != "") %>%
    dplyr::distinct(ENTREZID) %>%
    dplyr::pull(ENTREZID)
  
  unique(eg)
}

collapse_terms <- function(df, top_n = 5){
  if (is.null(df) || nrow(df) == 0) return(NA_character_)
  out <- df %>%
    dplyr::arrange(p.adjust) %>%
    dplyr::slice_head(n = min(top_n, nrow(df))) %>%
    dplyr::transmute(txt = paste0(Description, " (FDR=", signif(p.adjust, 2), ")")) %>%
    dplyr::pull(txt)
  paste(out, collapse = " | ")
}

run_GO_BP <- function(genes){
  eg <- sym2ent(genes)
  if (length(eg) < 10) return(NA_character_)
  
  out <- tryCatch({
    ego <- clusterProfiler::enrichGO(
      gene = eg,
      OrgDb = org.Hs.eg.db,
      ont = "BP",
      pAdjustMethod = "BH",
      qvalueCutoff = 0.2,
      readable = TRUE
    )
    collapse_terms(as.data.frame(ego))
  }, error = function(e) NA_character_)
  
  out
}

run_Reactome <- function(genes){
  eg <- sym2ent(genes)
  if (length(eg) < 10) return(NA_character_)
  
  out <- tryCatch({
    er <- ReactomePA::enrichPathway(
      gene = eg,
      organism = "human",
      pAdjustMethod = "BH",
      qvalueCutoff = 0.2,
      readable = TRUE
    )
    collapse_terms(as.data.frame(er))
  }, error = function(e) NA_character_)
  
  out
}

check_map_fail <- function(genes){
  genes <- unique(genes[!is.na(genes) & genes != ""])
  m <- suppressMessages(
    clusterProfiler::bitr(
      genes,
      fromType = "SYMBOL",
      toType   = "ENTREZID",
      OrgDb    = org.Hs.eg.db
    )
  )
  mapped <- unique(m$SYMBOL)
  setdiff(genes, mapped)
}


## ---
##  *SINGSCORE* ----
## ---
# 2026-06-07 20:40 SAB ---

plot_module_score <- function(df, module_name, title = NULL, y_lab = "SingScore",
                              output_dir = NULL, file_prefix = NULL, base_size = 10,
                              width = 8, height = 5, dpi = 600, units = "in", facet_by = NULL) {
  
  if (is.null(title)) {
    title <- module_name
  }
  
  plot_df <- df %>%
    dplyr::filter(module == module_name)
  
  p <- plot_df %>%
    ggplot(aes(x = day_num, y = score, color = genotype, group = genotype)) +
    
    # raw sample points
    geom_point(aes(group = interaction(genotype, timepoint)),
      size = 1, alpha = 0.7) +
    
    # error bares (mean SE)
    stat_summary( fun.data = mean_se,
      geom = "errorbar", width = 0.2, linewidth = 0.5) +
    
    # mean trajectory
    stat_summary(fun = mean, geom = "line", linewidth = 1) +
    stat_summary(fun = mean, geom = "point", size = 1.5) +
    
    scale_x_continuous(
      breaks = c(3, 4, 6, 8, 10),
      labels = c("D3", "D4", "D6", "D8", "D10")
    ) +
    scale_y_continuous(
      breaks = scales::pretty_breaks(n = 4),
      expand = expansion(mult = c(0.08, 0.08))
    ) +
    scale_color_manual(values = annotation_colors$genotype) +
    labs(title = title,
      caption = paste0("Source module: ", module_name),
      x = NULL, y = y_lab, color = NULL
    ) +
    theme_pub(base_size = base_size) +
    guides(
      color = guide_legend(override.aes = list(linewidth = 1.5, size = 3, alpha = 1)))
  
  if (!is.null(facet_by)) {
    p <- p + facet_wrap(as.formula(paste("~", facet_by)))
  }
  
  if (!is.null(output_dir)) {
    
    if (is.null(file_prefix)) {
      file_prefix <- module_name
    }
    
    ggsave(
      filename = file.path(output_dir, paste0(file_prefix, "_trajectory_singscore.pdf")),
      plot = p,
      width = width,
      height = height,
      dpi = dpi,
      units = units,
      device = cairo_pdf,
      bg = "white",
      
    )
  }
  
  return(p)
}


# top driver genes
# take one module, on set of samples and looks at ranked matrix 
get_module_driver_table <- function(module_name, ranked_mat,
                                    gene_sets, sample_ids,
                                    top_n = 20) {
  
  if (!module_name %in% names(gene_sets)) {
    stop("Module not found: ", module_name)
  }
  
  # Pull  genes & sampleIDs in  module and keep genes present in ranked matrix
  sig_genes <- gene_sets[[module_name]] %>%
    intersect(rownames(ranked_mat))
  
  sample_ids <- intersect(sample_ids, colnames(ranked_mat))
  
  if (length(sig_genes) == 0) {
    warning("No genes found in ranked_mat for module: ", module_name)
    return(tibble())
  }
  
  if (length(sample_ids) == 0) {
    warning("No samples found in ranked_mat for module: ", module_name)
    return(tibble())
  }
  
  
  # Subset the ranked matrix
  rank_sub <- ranked_mat[sig_genes, sample_ids, drop = FALSE]
  
  # For each gene, calculate: mean rank across  samples, SD then sort highest to lowest
  driver_df <- tibble(
    gene = rownames(rank_sub),
    mean_rank = rowMeans(rank_sub, na.rm = TRUE),
    sd_rank = apply(rank_sub, 1, sd, na.rm = TRUE),
    n_samples = length(sample_ids)
  ) %>%
    arrange(desc(mean_rank)) %>%
    mutate(
      module = module_name,
      rank_order = row_number()
    )
  
  driver_df <- driver_df %>%
    slice_head(n = min(top_n, nrow(driver_df)))
  
  return(driver_df)
}

## ---
##  *GSEA* ----
## ---

# 2026-04-13 18:29 SAB ---
# Rank genes for GSEA
prep_gsea_ranks <- function(res_df) {
  
  rank_df <- res_df %>%
    dplyr::filter(
      !is.na(gene_id),
      !is.na(log2FoldChange_raw),
      !is.na(pvalue)
    ) %>%
    dplyr::mutate(
      pvalue = pmax(pvalue, 1e-300),
      rank_metric = sign(log2FoldChange_raw) * -log10(pvalue)
    ) %>%
    dplyr::select(gene_id, rank_metric) %>%
    dplyr::group_by(gene_id) %>%
    dplyr::summarise(
      rank_metric = rank_metric[which.max(abs(rank_metric))],
      .groups = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(rank_metric))
  
  ranks <- rank_df$rank_metric
  names(ranks) <- rank_df$gene_id
  sort(ranks, decreasing = TRUE)
}


## ---
##  *MASIGPRO* ----
## ---
# 2026-05-01 13:49 SAB -

plot_cluster_summary <- function(km_obj, mat_z, metadata_plot, k_label,
                                 annotation_colors,
                                 save_dir = NULL,
                                 prefix = "CTRL",
                                 show_missing_genotypes = FALSE) {
  
  metadata_plot <- metadata_plot %>%
    dplyr::mutate(
      genotype = factor(genotype, levels = c("D21", "T21"))
    )
  
  # gene-to-cluster assignment
  cluster_assignments <- tibble::tibble(
    gene = rownames(mat_z),
    cluster = factor(km_obj$cluster)
  )
  
  # calculate one cluster score per sample
  cluster_scores <- as.data.frame(mat_z) %>%
    tibble::rownames_to_column("gene") %>%
    dplyr::left_join(cluster_assignments, by = "gene") %>%
    tidyr::pivot_longer(
      cols = -c(gene, cluster),
      names_to = "sampleID",
      values_to = "expression_z"
    ) %>%
    dplyr::group_by(cluster, sampleID) %>%
    dplyr::summarise(
      cluster_score = mean(expression_z, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::left_join(metadata_plot, by = "sampleID")
  
  # summarize across biological replicates
  cluster_summary <- cluster_scores %>%
    dplyr::group_by(cluster, Time, timepoint, genotype) %>%
    dplyr::summarise(
      mean_z = mean(cluster_score, na.rm = TRUE),
      sem_z = sd(cluster_score, na.rm = TRUE) / sqrt(dplyr::n()),
      n = dplyr::n(),
      .groups = "drop"
    )
  
  p <- ggplot(
    cluster_summary,
    aes(x = Time, y = mean_z, color = genotype, fill = genotype, group = genotype)
  ) +
    geom_errorbar(
      aes(ymin = mean_z - sem_z, ymax = mean_z + sem_z),
      width = 0.25,
      linewidth = 0.7
    ) +
    geom_line(linewidth = 1.3) +
    geom_point(size = 3) +
    scale_color_manual(
      values = annotation_colors$genotype,
      drop = !show_missing_genotypes
    ) +
    scale_fill_manual(
      values = annotation_colors$genotype,
      drop = !show_missing_genotypes
    ) +
    facet_wrap(~ cluster, scales = "free_y") +
    theme_classic(base_size = 14) +
    labs(
      title = paste0(prefix, " MaSigPro temporal clusters"),
      subtitle = paste0("k = ", k_label),
      x = "Day",
      y = "Scaled expression",
      color = "Genotype",
      fill = "Genotype"
    ) +
    theme(
      strip.background = element_rect(fill = "white", color = "black"),
      strip.text = element_text(face = "bold"),
      plot.title = element_text(face = "bold"),
      legend.title = element_text(face = "bold")
    )
  
  if (!is.null(save_dir)) {
    ggsave(
      file.path(save_dir, paste0(prefix, "_cluster_summary_k", k_label, ".pdf")),
      p,
      width = 8, height = 5, bg = "white")
  }
  
  return(p)
}

# module projection function (plotting with pre-defined modules)
plot_module_projection <- function(cluster_assignments, mat_z, metadata_plot,
                                   annotation_colors,
                                   group_var = "treatment",
                                   save_dir = NULL,
                                   prefix = "Projected_modules") {
  
  # ensure grouping variable has consistent order 
  if (group_var %in% names(annotation_colors)) {
    metadata_plot[[group_var]] <- factor(
      metadata_plot[[group_var]],
      levels = names(annotation_colors[[group_var]])
    )
  }
  
  # Make cluster labels categorical (for plotting/faceting)
  cluster_assignments <- cluster_assignments %>%
    dplyr::mutate(cluster = factor(cluster))
  
  # Calculate module score per sample
  cluster_scores <- as.data.frame(mat_z) %>%
    tibble::rownames_to_column("gene") %>%
    dplyr::inner_join(cluster_assignments, by = "gene") %>%
    tidyr::pivot_longer(
      cols = -c(gene, cluster),
      names_to = "sampleID",
      values_to = "expression_z"
    ) %>%
    dplyr::group_by(cluster, sampleID) %>%
    dplyr::summarise(
      module_score = mean(expression_z, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::left_join(metadata_plot, by = "sampleID") %>%
    dplyr::filter(!is.na(.data[[group_var]]))  
  
  # Summarize across replicates 
  module_summary <- cluster_scores %>%
    dplyr::group_by(cluster, Time, timepoint, .data[[group_var]]) %>%
    dplyr::summarise(
      mean_z = mean(module_score, na.rm = TRUE),
      sem_z = sd(module_score, na.rm = TRUE) / sqrt(dplyr::n()),
      n = dplyr::n(),
      .groups = "drop"
    )
  
  # Build plot 
  p <- ggplot(module_summary,
    aes(x = Time,
      y = mean_z,
      color = .data[[group_var]],
      fill = .data[[group_var]],
      group = .data[[group_var]]
    )
  ) +
    
    # error bars = SEM
    geom_errorbar(aes(ymin = mean_z - sem_z, ymax = mean_z + sem_z),
      width = 0.25, linewidth = 0.7) +
    geom_line(linewidth = 1.3) +
    geom_point(size = 3) +
    
    # apply predefined color palette
    scale_color_manual(values = annotation_colors[[group_var]]) +
    scale_fill_manual(values = annotation_colors[[group_var]]) +
    
    facet_wrap(~ cluster, scales = "free_y") +
    
    # clean styling
    theme_classic(base_size = 14) +
    labs(title = prefix,
      x = "Day",
      y = "Module score",
      color = group_var,
      fill = group_var
    ) +
    theme(
      strip.background = element_rect(fill = "white", color = "black"),
      strip.text = element_text(face = "bold"),
      plot.title = element_text(face = "bold"),
      legend.title = element_text(face = "bold")
    )
  
  # save plot 
  if (!is.null(save_dir)) {
    ggsave(file.path(save_dir, paste0(prefix, "_module_projection.pdf")),
      p, width = 9, height = 5, bg = "white")
  }
  return(list(
    plot = p,
    cluster_scores = cluster_scores,
    module_summary = module_summary
  ))
}


## ---
##  *SCATTERPLOT* ----
## ---

# GENERAL SCATTERPLOT 

get_days_for_contrast <- function(all_list, contrast) {
  if (is.null(all_list[[contrast]])) return(character(0))
  nms <- names(all_list[[contrast]])
  if (is.null(nms)) return(character(0))
  
  # Day token occurs right after T21_ or D21_
  days <- unique(unlist(stringr::str_extract_all(nms, "(?<=T21_)D\\d+|(?<=D21_)D\\d+")))
  days <- days[!is.na(days)]
  
  days[order(as.integer(sub("^D", "", days)))]
}


build_inner_key <- function(contrast, day) {
  sides <- strsplit(contrast, "_vs_", fixed = TRUE)[[1]]
  
  make_side <- function(side) {
    parts <- strsplit(side, "_", fixed = TRUE)[[1]]
    paste0(parts[1], "_", day, "_", paste(parts[-1], collapse = "_"))
  }
  paste0(make_side(sides[1]), "_vs_", make_side(sides[2]))
}



run_scatter_within_experiment <- function(all_list, contrast_x, contrast_y, days = NULL, exp_label,
                                          log_cut = 0.5, p_cut = 0.05, color_by = "DEG", autoprin = TRUE, 
                                          out_dir, show_cor = TRUE, show_fit = TRUE) {
  
  # If days not provided, plot all days that exist in BOTH contrasts
  if (is.null(days)) {
    dx <- get_days_for_contrast(all_list, contrast_x)
    dy <- get_days_for_contrast(all_list, contrast_y)
    days <- intersect(dx, dy)
  }
  
  if (length(days) == 0) {
    message("No shared days between ", contrast_x, " and ", contrast_y)
    return(invisible(NULL))
  }
  
  
  for (d in days) {
    
    key_x <- build_inner_key(contrast_x, d)
    key_y <- build_inner_key(contrast_y, d)
    
    if (is.null(all_list[[contrast_x]][[key_x]]) || is.null(all_list[[contrast_y]][[key_y]])) {
      message("Skipping ", d, " (missing key for one contrast).")
      next
    }
    
    df_x <- prep_for_scatter(all_list[[contrast_x]][[key_x]])
    df_y <- prep_for_scatter(all_list[[contrast_y]][[key_y]])
    
    # keep your function's "df1/df2 must be object names" requirement
    var_x <- paste0(exp_label, "_", d, "_X")
    var_y <- paste0(exp_label, "_", d, "_Y")
    assign(var_x, df_x, envir = .GlobalEnv)
    assign(var_y, df_y, envir = .GlobalEnv)
    
    md <- list()
    md[[var_x]] <- paste0(d, " ", gsub("_", " ", contrast_x))
    md[[var_y]] <- paste0(d, " ", gsub("_", " ", contrast_y))
    
    
    p <- plotScatterAnalysis(
      df1 = var_x,
      df2 = var_y,
      ColorBy = color_by,
      AutoPrint = autoprin,
      label_genes = character(0),
      logCutoff = log_cut,
      pvalCutoff = p_cut,
      Metadata = md,
      out_dir = out_dir,
      ReturnDataOnly = FALSE,
      show_cor = TRUE, show_fit = TRUE
    )
    
    safe_x <- gsub("_vs_", "-vs-", contrast_x, fixed = TRUE)
    safe_y <- gsub("_vs_", "-vs-", contrast_y, fixed = TRUE)
    
    
    fn <- file.path(out_dir, paste0(exp_label, "_", d, "_X-", safe_x, "_Y-", safe_y, "_scatter.pdf"))
    ggplot2::ggsave(fn, p, width = 8, height = 5, dpi = 300, bg = "transparent")
    message("Saved: ", fn)
  }
  
  invisible(TRUE)
}

prep_for_scatter <- function(x) {
  df <- as.data.frame(x)
  if (!"gene_id" %in% names(df)) df <- tibble::rownames_to_column(df, "gene_id")
  
  padj_col <- dplyr::case_when(
    "padj"   %in% names(df) ~ "padj",
    "qvalue" %in% names(df) ~ "qvalue",
    TRUE ~ NA_character_
  )
  if (is.na(padj_col)) stop("prep_for_scatter(): Missing adjusted p-value column (padj/qvalue).")
  
  # We want BOTH raw + shrunk available.
  # Raw: prefer explicit *_raw if present; otherwise fall back to log2FoldChange.
  raw_col <- dplyr::case_when(
    "log2FoldChange_raw" %in% names(df) ~ "log2FoldChange_raw",
    "log2FoldChange"     %in% names(df) ~ "log2FoldChange",
    "log2FC"             %in% names(df) ~ "log2FC",
    TRUE ~ NA_character_
  )
  if (is.na(raw_col)) stop("prep_for_scatter(): Missing raw LFC (log2FoldChange_raw or log2FoldChange/log2FC).")
  
  # Shrunk: prefer *_shrunk if present; otherwise fall back to raw (so plotting still works)
  shrunk_col <- dplyr::case_when(
    "log2FoldChange_shrunk" %in% names(df) ~ "log2FoldChange_shrunk",
    TRUE ~ raw_col
  )
  
  dplyr::transmute(df,
                   gene_id    = as.character(.data$gene_id),
                   padj       = as.numeric(.data[[padj_col]]),
                   lfc_raw    = as.numeric(.data[[raw_col]]),
                   lfc_shrunk = as.numeric(.data[[shrunk_col]])
  )
}


plotScatterAnalysis <- function(df1, df2, Metadata, out_dir,
                                ColorBy = "DEG",         # "DEG" or "Significance"
                                AutoPrint = TRUE,        # TRUE, FALSE, or "Manual"
                                label_genes = character(0),
                                logCutoff = 0.5, pvalCutoff = 0.05,
                                ReturnDataOnly = FALSE,
                                show_cor = TRUE, show_fit = TRUE) {
  
  d1 <- get(df1)
  d2 <- get(df2)
  
  # Join on common genes
  common_genes <- intersect(d1$gene_id, d2$gene_id)
  
  final_df <- dplyr::inner_join(
    d1 |>
      dplyr::filter(gene_id %in% common_genes) |>
      dplyr::select(gene_id,
                    padj_df1 = padj,
                    lfc_raw_df1 = lfc_raw,
                    lfc_shrunk_df1 = lfc_shrunk),
    d2 |>
      dplyr::filter(gene_id %in% common_genes) |>
      dplyr::select(gene_id,
                    padj_df2 = padj,
                    lfc_raw_df2 = lfc_raw,
                    lfc_shrunk_df2 = lfc_shrunk),
    by = "gene_id"
  )
  
  # Use shrunken for DISPLAY
  final_df <- final_df |>
    dplyr::mutate(
      lfc_plot_df1 = lfc_shrunk_df1,
      lfc_plot_df2 = lfc_shrunk_df2
    )
  
  # Significance + DEG definitions (RAW for cutoffs)
  final_df <- final_df |>
    dplyr::mutate(
      Sig_df1 = padj_df1 < pvalCutoff,
      Sig_df2 = padj_df2 < pvalCutoff,
      
      Significance = dplyr::case_when(
        Sig_df1 & Sig_df2 ~ "Significant in both",
        Sig_df1 & !Sig_df2 ~ paste0("Significant in ", Metadata[[df1]], " only"),
        !Sig_df1 & Sig_df2 ~ paste0("Significant in ", Metadata[[df2]], " only"),
        TRUE ~ "Not significant"
      ),
      
      # DEG calling uses RAW (this is what you want)
      is_all = (Sig_df1 & Sig_df2 &
                  abs(lfc_raw_df1) >= logCutoff & abs(lfc_raw_df2) >= logCutoff),
      is_df1 = (Sig_df1 & !Sig_df2 &
                  abs(lfc_raw_df1) >= logCutoff),
      is_df2 = (Sig_df2 & !Sig_df1 &
                  abs(lfc_raw_df2) >= logCutoff),
      
      DEGcat = dplyr::case_when(
        is_all ~ "DEG in both",
        is_df1 ~ paste0("DEG in ", Metadata[[df1]], " only"),
        is_df2 ~ paste0("DEG in ", Metadata[[df2]], " only"),
        TRUE ~ "Non DEG"
      ),
      
      # threshold quadrant ONLY for significant-in-both (RAW)
      Quadrant = dplyr::case_when(
        is_all & lfc_raw_df1 < -logCutoff & lfc_raw_df2 >=  logCutoff ~ "Quad.TopL",
        is_all & lfc_raw_df1 >=  logCutoff & lfc_raw_df2 >=  logCutoff ~ "Quad.TopR",
        is_all & lfc_raw_df1 >=  logCutoff & lfc_raw_df2 <  -logCutoff ~ "Quad.BottomR",
        is_all & lfc_raw_df1 <  -logCutoff & lfc_raw_df2 <  -logCutoff ~ "Quad.BottomL",
        is_df1 ~ "df1only",
        is_df2 ~ "df2only",
        TRUE ~ "Not selected"
      ),
      
      # For CORNER PLACEMENT use DISPLAY signs (so numbers live where your dots visually are)
      SignQuad = dplyr::case_when(
        lfc_plot_df1 >= 0 & lfc_plot_df2 >= 0 ~ "TopR",
        lfc_plot_df1 <  0 & lfc_plot_df2 >= 0 ~ "TopL",
        lfc_plot_df1 <  0 & lfc_plot_df2 <  0 ~ "BottomL",
        lfc_plot_df1 >= 0 & lfc_plot_df2 <  0 ~ "BottomR",
        TRUE ~ NA_character_
      )
    )
  
  # Save DEG tables 
  safe1 <- gsub("[^A-Za-z0-9]+", "_", Metadata[[df1]])
  safe2 <- gsub("[^A-Za-z0-9]+", "_", Metadata[[df2]])
  
  excel_name <- paste0("scatter_tables__", safe1, "__vs__", safe2, ".xlsx")
  
  deg_union <- final_df |> 
    dplyr::filter(is_all | is_df1 | is_df2)
  
  writexl::write_xlsx(
    setNames(
      list(
        deg_union,
        dplyr::filter(final_df, is_all),
        dplyr::filter(final_df, is_df1),
        dplyr::filter(final_df, is_df2)
      ),
      c(
        "All_DEGs",
        paste0("Shared"),
        paste0("Only__", safe1),
        paste0("Only__", safe2)
      )
    ),
    path = file.path(out_dir, excel_name)
  )
  
  # Color palette / column
  if (ColorBy == "Significance") {
    colors <- setNames(
      c("#D81B60", "#1E88E5"),
      c(paste0("Significant in ", Metadata[[df1]], " only"),
        paste0("Significant in ", Metadata[[df2]], " only"))
    )
    colors <- c("Significant in both" = "#FFC107", colors, "Not significant" = "#7C7C7C")
    colColumn <- "Significance"
  } else {
    colors <- setNames(
      c("#D81B60", "#1E88E5"),
      c(paste0("DEG in ", Metadata[[df1]], " only"),
        paste0("DEG in ", Metadata[[df2]], " only"))
    )
    colors <- c("DEG in both" = "#FFC107", colors, "Non DEG" = "#7C7C7C")
    colColumn <- "DEGcat"
  }
  
  # Plot limits based on DISPLAY lfc (what you see)
  q <- 0.998
  lim_x <- quantile(final_df$lfc_plot_df1, probs = c(1 - q, q), na.rm = TRUE)
  lim_y <- quantile(final_df$lfc_plot_df2, probs = c(1 - q, q), na.rm = TRUE)
  M <- max(abs(lim_x), abs(lim_y), na.rm = TRUE)
  buffer <- 0.20 * M
  xLIM <- c(-(M + buffer), (M + buffer))
  yLIM <- xLIM
  
  # Top labels for DEGs in both: pick by RAW magnitude (honest)
  top_labels <- final_df |>
    dplyr::filter(is_all, Quadrant %in% c("Quad.TopL","Quad.TopR","Quad.BottomR","Quad.BottomL")) |>
    dplyr::group_by(Quadrant) |>
    dplyr::slice_max(order_by = abs(lfc_raw_df1) + abs(lfc_raw_df2), n = 5) |>
    dplyr::ungroup()
  
  # Corner positions
  quad_pos <- tibble::tibble(
    SignQuad = c("TopL", "TopR", "BottomR", "BottomL"),
    x = c(xLIM[1], xLIM[2], xLIM[2], xLIM[1]),
    y = c(yLIM[2], yLIM[2], yLIM[1], yLIM[1]),
    hjust = c(0, 1, 1, 0),
    vjust = c(1, 1, 0, 0)
  )
  
  # Counts per sign quadrant (RAW DEG definitions, placed by DISPLAY sign)
  count_all <- final_df |> dplyr::filter(is_all) |> dplyr::count(SignQuad, name = "n")
  count_df1 <- final_df |> dplyr::filter(is_df1) |> dplyr::count(SignQuad, name = "n")
  count_df2 <- final_df |> dplyr::filter(is_df2) |> dplyr::count(SignQuad, name = "n")
  
  dy1 <- 0.03 * diff(yLIM)
  dy2 <- 0.09 * diff(yLIM)
  dy3 <- 0.15 * diff(yLIM)
  
  quad_all <- quad_pos |>
    dplyr::left_join(count_all, by = "SignQuad") |>
    dplyr::mutate(
      n = tidyr::replace_na(n, 0),
      x = x + ifelse(hjust == 0, 0.03 * diff(xLIM), -0.03 * diff(xLIM)),
      y = y + ifelse(vjust == 1, -dy1, dy1)
    )
  
  quad_df1 <- quad_pos |>
    dplyr::left_join(count_df1, by = "SignQuad") |>
    dplyr::mutate(
      n = tidyr::replace_na(n, 0),
      x = x + ifelse(hjust == 0, 0.03 * diff(xLIM), -0.03 * diff(xLIM)),
      y = y + ifelse(vjust == 1, -dy2, dy2)
    )
  
  quad_df2 <- quad_pos |>
    dplyr::left_join(count_df2, by = "SignQuad") |>
    dplyr::mutate(
      n = tidyr::replace_na(n, 0),
      x = x + ifelse(hjust == 0, 0.03 * diff(xLIM), -0.03 * diff(xLIM)),
      y = y + ifelse(vjust == 1, -dy3, dy3)
    )
  
  # Plot order (so "both" points are drawn last = on top)
  plot_df <- final_df |>
    dplyr::arrange(factor(Significance, levels = c(
      paste0("Significant in ", Metadata[[df1]], " only"),
      paste0("Significant in ", Metadata[[df2]], " only"),
      "Not significant",
      "Significant in both"
    )))
  
  # Correlation (SHARED DEGs ONLY) ---
  cor_label <- NULL
  if (isTRUE(show_cor)) {
    df_shared <- plot_df |> dplyr::filter(is_all)
    
    if (nrow(df_shared) >= 3) {
      cor_stats <- stats::cor.test(df_shared$lfc_plot_df1, df_shared$lfc_plot_df2, method = "pearson")
      cor_label <- paste0(
        "Pearson r (shared DEGs) = ", round(unname(cor_stats$estimate), 3),
        "\np = ", format.pval(cor_stats$p.value, digits = 3, eps = .001)
      )
    } else {
      cor_label <- "Pearson r (shared DEGs): NA\n(n<3)"
    }
  }
  
  
  plot <- ggplot2::ggplot(plot_df,
                          ggplot2::aes(x = lfc_plot_df1, y = lfc_plot_df2, col = .data[[colColumn]])) +
    ggplot2::geom_point(shape = 16, size = 1.5, alpha = 0.8) +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::ggtitle(paste0("Comparison between ", Metadata[[df1]], " and ", Metadata[[df2]])) +
    ggplot2::xlab(paste0("log2FC (shrunken) - ", Metadata[[df1]])) +
    ggplot2::ylab(paste0("log2FC (shrunken) - ", Metadata[[df2]])) +
    ggplot2::coord_cartesian(xlim = xLIM, ylim = yLIM) +
    # cutoff lines drawn at raw cutoff value (still same numeric threshold)
    ggplot2::geom_vline(xintercept = c(-logCutoff, logCutoff), linetype = "dashed",
                        color = "#7C7C7C", linewidth = 0.5) +
    ggplot2::geom_hline(yintercept = c(-logCutoff, logCutoff), linetype = "dashed",
                        color = "#7C7C7C", linewidth = 0.5) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.border = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(colour = "#000000"),
      legend.position = "right",
      legend.title = ggplot2::element_blank()
    ) +
    # Corner numbers: yellow (both) + df1-only (pink) + df2-only (blue)
    ggplot2::geom_label(
      data = quad_all,
      ggplot2::aes(x = x, y = y, label = n, hjust = hjust, vjust = vjust),
      inherit.aes = FALSE, size = 4, color = "#FFC107", fill = "white",
      alpha = 0.85, label.size = 0,
      label.r = grid::unit(0.12, "lines"),
      label.padding = grid::unit(0.10, "lines")
    ) +
    ggplot2::geom_label(
      data = quad_df1,
      ggplot2::aes(x = x, y = y, label = n, hjust = hjust, vjust = vjust),
      inherit.aes = FALSE, size = 4, color = "#D81B60", fill = "white",
      alpha = 0.85, label.size = 0,
      label.r = grid::unit(0.12, "lines"),
      label.padding = grid::unit(0.10, "lines")
    ) +
    ggplot2::geom_label(
      data = quad_df2,
      ggplot2::aes(x = x, y = y, label = n, hjust = hjust, vjust = vjust),
      inherit.aes = FALSE, size = 4, color = "#1E88E5",   fill = "white",
      alpha = 0.85, label.size = 0,
      label.r = grid::unit(0.12, "lines"),
      label.padding = grid::unit(0.10, "lines")
    ) 
  
  # optional regression line
  if (isTRUE(show_fit)) {
    plot <- plot + ggplot2::geom_smooth(method = "lm", se = FALSE, color = "#444444",
                                        linetype = "solid", linewidth = 0.6)
  }
  
  # optional correlation annotation
  if (isTRUE(show_cor) && !is.null(cor_label)) {
    plot <- plot + ggplot2::annotate(
      "text",
      x = xLIM[1] + 0 * diff(xLIM),
      y = yLIM[2] + 0.05 * diff(yLIM),
      label = cor_label,
      hjust = 0, vjust = 1, size = 3, fontface = "italic", color = "#444444"
    )
  }
  
  
  if (isTRUE(AutoPrint)) {
    plot <- plot +
      ggrepel::geom_text_repel(
        data = top_labels,
        ggplot2::aes(label = gene_id, color = DEGcat),
        fontface = "bold",
        size = 3,
        show.legend = FALSE,
        box.padding = 0.6,
        bg.color = "white",
        bg.r = 0.15
      )
  } else if (identical(AutoPrint, "Manual")) {
    plot <- plot +
      ggrepel::geom_text_repel(
        data = plot_df |> dplyr::filter(gene_id %in% label_genes),
        ggplot2::aes(label = gene_id, color = DEGcat),
        fontface = "bold",
        size = 3,
        show.legend = FALSE,
        max.overlaps = Inf
      )
  }
  
  if (!ReturnDataOnly) return(plot) else return(plot_df)
}




## ---
##  PLOT THEME ----
## ---
# 2026-06-07 18:15 SAB ---
theme_pub <- function(base_size = 10, base_family = "helvetica_neue") {
  library(grid)
  library(ggthemes)
  
  theme_foundation(base_size = base_size, base_family = base_family) +
    theme(
      plot.title = element_text(
        face = "bold", size = rel(1.2), hjust = 0.5,
        margin = margin(6, 6, 8,  6)),
      text = element_text(color = "black"),
      panel.background = element_rect(colour = NA, fill = NA), #makes background transparent
      plot.background = element_rect(colour = NA, fill = NA),
      panel.border = element_rect(colour = NA),
      
      # Axis titles and text sizes 
      axis.title.x = element_text(vjust = -0.2, size = rel(1.1)),
      axis.title.y = element_text(angle = 90, vjust = 2, size = rel(1.1)),
      axis.text.x = element_text(size = rel(1.0), color = "black"),
      axis.text.y = element_text(size = rel(1.0), color = "black"),
      axis.ticks.length = unit(0.1, "cm"),
      
      axis.line.x = element_line(colour = "black", linewidth = 0.4),
      axis.line.y = element_line(colour = "black", linewidth = 0.4),
      axis.ticks = element_line(colour = "black", linewidth = 0.4),
      
      #panel.grid.major = element_line(colour = "#f0f0f0"), #if you want grid lines
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      
      legend.key = element_rect(colour = NA, fill = NA),
      legend.position = "top",
      legend.direction = "horizontal",
      legend.box = "vertical",
      legend.key.size = unit(0.35, "cm"),
      legend.title = element_text(size = rel(0.9), face = "plain"),
      legend.text = element_text(size = rel(0.8)),
      
      #caption settings
      plot.caption = element_text(
        size = rel(0.6), color = "gray30",
        hjust = 0,             # left-aligned
        vjust = 0,           # shift it down
        margin = margin(t = 8, b = 0)
      ),
      plot.caption.position = "plot",
      
      # for facet labelts 
      strip.background = element_rect(colour = "#f0f0f0", fill = "#f0f0f0"),
      strip.text = element_text(size = rel(0.9), face = "plain"),
      
      plot.margin = margin(4, 4, 4, 4)
      
    )
}

