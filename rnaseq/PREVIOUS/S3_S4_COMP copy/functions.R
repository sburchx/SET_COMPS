make_subdir <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE)
    message(paste("📁 Folder created:", path))
  } else {
    message(paste("📁 Folder already exists:", path))
  }
}

# ---
# Find the inner element for a given day 
get_day_df <- function(DEG_all, contrast, day) {
  stopifnot(contrast %in% names(DEG_all))
  
  inner <- DEG_all[[contrast]]
  inner_names <- names(inner)
  if (is.null(inner_names)) stop("Inner list has no names; can't find day keys.")
  
  # robust match: any key that contains "D4" as a token (handles "D4", "day_D4", "D4_something")
  hits <- inner_names[str_detect(inner_names, regex(paste0("\\b", day, "\\b"), ignore_case = TRUE))]
  
  # fallback: if day isn't a clean token, match just "D4" anywhere
  if (length(hits) == 0) {
    hits <- inner_names[str_detect(inner_names, fixed(day, ignore_case = TRUE))]
  }
  
  if (length(hits) == 0) {
    stop(paste0("No inner key matched day=", day, " for contrast=", contrast,
                ". Available keys: ", paste(inner_names, collapse = ", ")))
  }
  
  # If multiple matches (rare), choose first but warn
  if (length(hits) > 1) {
    warning(paste0("Multiple keys matched ", day, " for ", contrast, ": ",
                   paste(hits, collapse = ", "), ". Using: ", hits[1]))
  }
  
  as.data.frame(inner[[hits[1]]])
}

# ---
# filter DEGS
filter_deg_genes <- function(df,
                             padj_cut = 0.05,
                             lfc_cut = 0.5,
                             lfc_col = "log2FoldChange_raw",
                             gene_col_candidates = c("gene_id","gene","symbol")) {
  
  if (!("padj" %in% names(df))) stop("No 'padj' column found.")
  if (!(lfc_col %in% names(df))) stop(paste0("No '", lfc_col, "' column found."))
  
  gene_col <- gene_col_candidates[gene_col_candidates %in% names(df)][1]
  if (is.na(gene_col)) stop("No gene column found (expected gene_id/gene/symbol).")
  
  df_f <- df %>%
    filter(!is.na(padj), padj < padj_cut) %>%
    filter(!is.na(.data[[lfc_col]]), abs(.data[[lfc_col]]) > lfc_cut)
  
  genes_all <- unique(df_f[[gene_col]])
  genes_up  <- unique(df_f[[gene_col]][df_f[[lfc_col]] > 0])
  genes_dn  <- unique(df_f[[gene_col]][df_f[[lfc_col]] < 0])
  
  list(df = df_f, all = genes_all, up = genes_up, down = genes_dn, gene_col = gene_col)
}


# ---
# TPM Function

plot_TPM_for_genes <- function(genes_of_interest, tpm, metadata, out_dir,
                               file_prefix = exp_name, width = 5, height = 6,
                               log10_y = FALSE, pal_group = NULL) {
  
  for (gene_of_interest in genes_of_interest) {
    
    tpm_subset <- tpm[rownames(tpm) == gene_of_interest, , drop = FALSE]
    if (nrow(tpm_subset) == 0) next
    
    # Melt the TPM dataframe for easier plotting
    tpm_melt <- melt(tpm_subset)
    colnames(tpm_melt) <- c("sampleID", "TPM") 
    
    # Merge with metadata
    tpm_melt <- merge(tpm_melt, metadata, by.x = "sampleID", by.y = "sampleID") 
    tpm_melt$genotype_day <- as.factor(tpm_melt$genotype_timepoint) 
    
    y_scale <- if (log10_y) scale_y_log10() else
      scale_y_continuous(expand = expansion(mult = c(0, 0.1)), limits = c(0, NA))
    
    p <- ggplot(tpm_melt, aes(x = timepoint, y = TPM, fill = genotype_treatment, group = genotype_timepoint_treatment)) +
      geom_boxplot(position = position_dodge(width = 0), alpha = 0.9) +  
      geom_jitter(position = position_jitter(0.15), size = 1) +  
      stat_summary(fun = mean, geom = "line", aes(group = genotype_treatment, color = genotype_treatment), 
                   size = .8, linetype = "solid") + 
      labs(title = paste(gene_of_interest), 
           x = "Timepoint",
           y = "TPM (Transcripts Per Million)",
           fill = NULL,  
           color = NULL) +  
      y_scale +
      theme_pub() +  
      theme(
        plot.title.position = "panel",
        plot.title = element_text(face = "bold"),
        panel.grid.major = element_line(colour = "grey92"),
        panel.grid.minor = element_line(colour = "grey97"),
        legend.position = "right",
        legend.direction = "vertical")
    
    if (!is.null(pal_group)) {
      keep <- sort(unique(tpm_melt$genotype_treatment))
      pal_use <- pal_group[names(pal_group) %in% keep]
      
      p <- p +
        scale_fill_manual(values = pal_use, guide = guide_legend(override.aes = list(alpha = 1))) +
        scale_color_manual(values = pal_use)
    }
    
    ggsave(file.path(out_dir, sprintf("%s_%s_TPM_boxplot.pdf", file_prefix, gene_of_interest)),
           plot = p, width = width, height = height, device = cairo_pdf)
    
  }
}


# theme pub
theme_pub <- function(base_size = 12, base_family = "helvetica_neue") {
  library(grid)
  library(ggthemes)
  
  theme_foundation(base_size = base_size, base_family = base_family) +
    theme(
      plot.title = element_text(
        face = "bold", size = rel(1.3), hjust = 0.5,
        margin = margin(10, 10, 20,  10)),
      text = element_text(),
      panel.background = element_rect(colour = NA, fill = NA), #makes background transparent
      plot.background = element_rect(colour = NA, fill = NA),
      panel.border = element_rect(colour = NA),
      
      # Axis titles and text sizes 
      axis.title.x = element_text(vjust = -0.2, size = rel(1.2)),
      axis.title.y = element_text(angle = 90, vjust = 2, size = rel(1.2)),
      axis.text.x = element_text(size = rel(1), color = "black"),
      axis.text.y = element_text(size = rel(1), color = "black"),
      
      axis.line.x = element_line(colour = "black"),
      axis.line.y = element_line(colour = "black"),
      axis.ticks = element_line(),
      
      #panel.grid.major = element_line(colour = "#f0f0f0"), #if you want grid lines
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      
      legend.key = element_rect(colour = NA),
      legend.position = "top",
      legend.direction = "horizontal",
      legend.box = "vertical",
      legend.key.size = unit(0.5, "cm"),
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 12),
      
      #caption settings
      plot.caption = element_text(
        size = rel(0.6), color = "gray30",
        hjust = 0,             # left-aligned
        vjust = 0,           # shift it down
        margin = margin(t = 12, b = 0)
      ),
      plot.caption.position = "plot",
      
      # for facet labelts 
      strip.background = element_rect(colour = "#f0f0f0", fill = "#f0f0f0"),
      strip.text = element_text(size = 14, face = "bold")
    )
}
