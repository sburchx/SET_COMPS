# Scatter plot functions

make_subdir <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE)
    message(paste("📁 Folder created:", path))
  } else {
    message(paste("📁 Folder already exists:", path))
  }
}


plotScatterAnalysis <- function(df1, # Dataset 1 (x-axis). Must be a string of the exact same name as the variable containing the data
                                df2, # Dataset 2 (y-axis). Must be a string of the exact same name as the variable containing the data
                                ColorBy = "DEG", # "DEG" or "Significance" // A gene can be significant but non DEG because of a log2FC cutoff
                                AutoPrint = T, # T, F or "Manual" // T to auto print up to 5 gene labels DEG in both per quadrant, F to avoid plotting gene labels, Manual to use a specific gene list "label_genes"
                                label_genes = c("Gli1", "Hhip", "Ptch1"), # print specific names on the plot; need AutoPrint = "Manual"
                                logCutoff_vitro = 2, # Log2 fold-change cutoff
                                logCutoff_vivo = 0.5, # Log2 fold-change cutoff
                                pvalCutoff = 0.05, # Adjusted p-value cutoff
                                Metadata = Metadata, # list of dataset variable and annotation to be printed as specific label
                                ReturnDataOnly = F # T or F, return the pre-processed dataset prior to plotting (to get quadrant info, gene lists, ...)
) {
  
  common_genes <- intersect(get(df1)$gene_id, get(df2)$gene_id)
  
  final_df <- dplyr::inner_join(
    get(df1) |> dplyr::filter(gene_id %in% common_genes) |>
      dplyr::select(gene_id,
                    padj_df1 = padj,
                    log2FC_df1 = log2FoldChange),
    get(df2) |> dplyr::filter(gene_id %in% common_genes) |>
      dplyr::select(gene_id,
                    padj_df2 = padj,
                    log2FC_df2 = log2FoldChange),
    by = "gene_id"
  )
  
  final_df <- final_df |> 
    dplyr::mutate(
      Significance = dplyr::case_when(
        padj_df1 < pvalCutoff & padj_df2 < pvalCutoff ~ "Significant in both",
        padj_df1 < pvalCutoff & padj_df2 >= pvalCutoff ~ paste0("Significant in ",Metadata[[df1]]," only"),
        padj_df1 >= pvalCutoff & padj_df2 < pvalCutoff ~ paste0("Significant in ",Metadata[[df2]]," only"),
        TRUE ~ "Not significant"
      ),
      Quadrant = dplyr::case_when(
        Significance == "Significant in both" & log2FC_df1 < -logCutoff_vivo & log2FC_df2 >= logCutoff_vitro ~ "Quad.TopL",
        Significance == "Significant in both" & log2FC_df1 >= logCutoff_vivo & log2FC_df2 >= logCutoff_vitro ~ "Quad.TopR",
        Significance == "Significant in both" & log2FC_df1 >= logCutoff_vivo & log2FC_df2 < -logCutoff_vitro ~ "Quad.BottomR",
        Significance == "Significant in both" & log2FC_df1 < -logCutoff_vivo & log2FC_df2 < -logCutoff_vitro ~ "Quad.BottomL",
        Significance == paste0("Significant in ",Metadata[[df1]]," only") & (log2FC_df1 < -logCutoff_vivo | log2FC_df1 >= logCutoff_vivo) ~ "df1only",
        Significance == paste0("Significant in ",Metadata[[df2]]," only") & (log2FC_df2 < -logCutoff_vitro | log2FC_df2 >= logCutoff_vitro) ~ "df2only",
        TRUE ~ "Not selected"
      ),
      DEGcat = dplyr::case_when(
        Quadrant %in% c("Quad.TopL","Quad.TopR","Quad.BottomR","Quad.BottomL") ~ "DEG in both",
        Quadrant == "df1only" ~ paste0("DEG in ", Metadata[[df1]], " only"),
        Quadrant == "df2only" ~ paste0("DEG in ", Metadata[[df2]], " only"),
        TRUE ~ "Non DEG"
      )
    )
  
  if(ColorBy == "Significance"){
    colors <- setNames(
      c("#D81B60", "#1E88E5"),
      c(paste0("Significant in ", Metadata[[df1]], " only"),
        paste0("Significant in ", Metadata[[df2]], " only"))
    )
    colors <- c("Significant in both" = "#FFC107", colors, "Not significant" = "#7C7C7C")
    colColumn <- "Significance"
  } else{
    colors <- setNames(
      c("#D81B60", "#1E88E5"),
      c(paste0("DEG in ", Metadata[[df1]], " only"),
        paste0("DEG in ", Metadata[[df2]], " only"))
    )
    colors <- c("DEG in both" = "#FFC107", colors, "Non DEG" = "#7C7C7C")
    colColumn <- "DEGcat"
  }
  
  # plot limits
  q <- 0.998   # keep central 99% (adjust to 0.99, 0.98, etc. if you want tighter zoom)
  lim_x <- quantile(final_df$log2FC_df1, probs = c(1 - q, q), na.rm = TRUE)
  lim_y <- quantile(final_df$log2FC_df2, probs = c(1 - q, q), na.rm = TRUE)
  
  # choose the larger absolute bound so axes are symmetric
  M <- max(abs(lim_x), abs(lim_y), na.rm = TRUE)
  xLIM <- yLIM <- c(-M, M)
  
  # add a small buffer 
  buffer <- 0.20 * M
  xLIM <- yLIM <- c(-(M + buffer), M + buffer)
  
  top_labels <- final_df |> 
    filter(DEGcat == "DEG in both", Quadrant %in% c("Quad.TopL", "Quad.TopR", "Quad.BottomR", "Quad.BottomL")) |>
    group_by(Quadrant) |>
    slice_max(order_by = abs(log2FC_df1) + abs(log2FC_df2), n = 5) |>  # pick top 5 by log2FC
    ungroup()
  
  quad_labels <- tibble::tibble(
    Quadrant = c("Quad.TopL", "Quad.TopR", "Quad.BottomR", "Quad.BottomL"),
    x = c(xLIM[1], xLIM[2], xLIM[2], xLIM[1]),
    y = c(yLIM[2], yLIM[2], yLIM[1], yLIM[1])
  ) |> 
    left_join(
      final_df |>
        filter(DEGcat == "DEG in both") |>
        group_by(Quadrant) |>
        summarise(n = n(), .groups = "drop"),
      by = "Quadrant"
    ) |> 
    mutate(n = tidyr::replace_na(n, 0))
  
  plot_df <- final_df |> 
    dplyr::arrange(factor(Significance, levels = c(
      paste0("Significant in ", Metadata[[df1]], " only"),
      paste0("Significant in ", Metadata[[df2]], " only"),
      "Not significant",
      "Significant in both"  # render last = printed on top
    )))
  
  cor_stats <- cor.test(plot_df$log2FC_df1, plot_df$log2FC_df2, method = "pearson")
  cor_label <- paste0("Pearson r = ", round(cor_stats$estimate, 3),
                      "\np = ", format.pval(cor_stats$p.value, digits = 3, eps = .001))
  
  plot <- ggplot(plot_df, aes(x=log2FC_df1, y=log2FC_df2, col=get(colColumn))) +
    geom_point(shape=16) +
    scale_color_manual(values = colors) +
    ggtitle(paste0("Comparison between ",Metadata[[df1]]," and ",Metadata[[df2]])) +
    xlab(paste0("log2FC - ",Metadata[[df1]]))+
    ylab(paste0("log2FC - ",Metadata[[df2]]))+
    xlim(xLIM) +
    ylim(yLIM) +
    geom_vline(xintercept=-logCutoff_vivo, linetype="dashed", color = "#7C7C7C", linewidth=0.5) +
    geom_vline(xintercept=logCutoff_vivo, linetype="dashed", color = "#7C7C7C", linewidth=0.5) +
    geom_hline(yintercept=-logCutoff_vitro, linetype="dashed", color = "#7C7C7C", linewidth=0.5) +
    geom_hline(yintercept=logCutoff_vitro, linetype="dashed", color = "#7C7C7C", linewidth=0.5) +
    theme_bw() +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "#000000"),
          legend.position="right",legend.title = element_blank())+
    geom_text(data = quad_labels,
              aes(x = x, y = y, label = n),
              size = 5, color = "#FFC107")+
    geom_smooth(method = "lm", se = FALSE, color = "#444444", linetype = "solid", linewidth = 0.6) +
    #geom_smooth(method = "lm", se = FALSE, linewidth = 0.6) + # seperate line for all sets of dots
    #geom_smooth(data = dplyr::filter(plot_df, DEGcat != "Non DEG"), aes(x = log2FC_df1, y = log2FC_df2, group = 1), method = "lm", se = FALSE, color = "#444444", linewidth = 0.6) # line for only DEGS
    annotate("text", x = xLIM[1] + 0.05 * diff(xLIM), y = yLIM[2] - 0.05 * diff(yLIM),
             label = cor_label, hjust = 0, vjust = 1, size = 4, fontface = "italic", color = "#444444")
  
  
  if(AutoPrint == T){
    plot <- plot + geom_text_repel(data = top_labels, aes(label = gene_id, color = DEGcat), fontface = "bold", size = 3, show.legend = FALSE,box.padding = 0.6, bg.color = "white", bg.r = 0.15)
  } else if(AutoPrint == "Manual"){
    plot <- plot + ggrepel::geom_text_repel(data = plot_df |> dplyr::filter(gene_id %in% label_genes), aes(label = gene_id, color = DEGcat), fontface = "bold",size = 3, show.legend = FALSE, max.overlaps = Inf)
  }
  
  if(ReturnDataOnly == F){
    return(plot)
  } else{return(plot_df)}
}

## ---
##  *SCATTERPLOT* ----
## ---


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
