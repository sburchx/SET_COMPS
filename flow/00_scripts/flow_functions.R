## ---
## HELPER FUNCTIONS ----
## ---

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
## GG PLOT ----
## ---

make_flow_plot <- function(measurement_name, summary_data, raw_data,
                           treatment_val = NULL, cell_val = NULL, facet_by = NULL,
                           folder, exp_name, color_map) {
  
  # Get pretty y-axis label
  label_name <- get_label(measurement_name)
  
  # Subset data by measurement
  summary_data <- summary_data %>% filter(measurement == measurement_name)
  raw_data <- raw_data %>% filter(measurement == measurement_name)
  
  # Filter by treatment or cell if requested
  if (!is.null(treatment_val)) {
    summary_data <- summary_data %>% filter(treatment == treatment_val)
    raw_data     <- raw_data     %>% filter(treatment == treatment_val)
  }
  if (!is.null(cell_val)) {
    summary_data <- summary_data %>% filter(str_detect(cell, paste0(cell_val, "$")))
    raw_data     <- raw_data     %>% filter(str_detect(cell, paste0(cell_val, "$")))
  }
  
  # Ensure treatment_label is factor with same order
  if (!"treatment_label" %in% colnames(summary_data)) {
    stop("summary_data is missing `treatment_label` column.")
  }
  if (!"treatment_label" %in% colnames(raw_data)) {
    stop("raw_data is missing `treatment_label` column.")
  }
  
  treatment_levels <- unique(summary_data$treatment_label)
  summary_data$treatment_label <- factor(summary_data$treatment_label, levels = treatment_levels)
  raw_data$treatment_label     <- factor(raw_data$treatment_label,     levels = treatment_levels)
  
  folder_short <- tail(strsplit(folder, "/")[[1]], 2)
  treatment_label_map <- raw_data %>%
    mutate(
      treatment = as.character(treatment),
      treatment_label = as.character(treatment_label)
    ) %>%
    select(treatment, treatment_label) %>%
    distinct()
  
  # Now join with color_map safely
  legend_map <- color_map %>%
    mutate(treatment_label = treatment_label_map$treatment_label[match(Group, treatment_label_map$treatment)]) %>%
    filter(!is.na(treatment_label))
  
  # Final named vector: treatment_label named by Group
  label_vector <- setNames(legend_map$treatment_label, legend_map$Group)
  
  # Plot
  p <- ggplot(summary_data, aes(x = factor(day), y = mean, fill = treatment, color = treatment)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9),
             width = 0.75, size = 1, aes(alpha = genotype, group = interaction(cell, treatment))) +
    geom_errorbar(
      aes(ymin = mean - sem, ymax = mean + sem, group = interaction(cell, treatment)),
      width = 0.2, size = 0.5, position = position_dodge(width = 0.9), color = "black"
    ) +
    geom_point(data = raw_data, aes(x = factor(day), y = value, group = interaction(cell, treatment), shape = as.factor(exp_no)),
               position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.9),
               size = 1.5, alpha = 0.8, color = "black", show.legend = FALSE) +
    scale_y_continuous(limits = c(0, 105), breaks = seq(0, 100, 20),
                       labels = scales::percent_format(scale = 1)) +
    scale_fill_manual(name = "D21", values = setNames(color_map$Color, color_map$Group),
                      labels = label_vector) +
    scale_color_manual(name = "T21", values = setNames(color_map$Color, color_map$Group),
                       labels = label_vector) +
    scale_alpha_manual(values = c("D21" = 1, "T21" = 0)) +
    labs(
      title = paste0(exp_name, " - ", label_name, "\n",
                     if (!is.null(treatment_val)) paste(raw_data$treatment_label), 
                     if (!is.null(cell_val)) paste0(cell_val)),
      x = "Day",
      y = label_name,
      caption = paste("Source:\n", paste(folder_short, collapse = "/"))
    ) +
    guides(
      fill = guide_legend(order = 1, override.aes = list(alpha = 1)),
      color = guide_legend(order = 2, override.aes = list(fill = NA)),
      alpha = "none"
    ) +
    theme_pub() 
  
  if (!is.null(facet_by)) {
    p <- p + facet_wrap(as.formula(paste0("~", facet_by)))
  }
  
  return(p)
}

## -- PULL PRETTY LABELS FOR PLOT -- ##

get_label <- function(m) {
  measurements$measurement_label[measurements$measurement == m]
}

## --- 
## STATS FUNCTIONS ----
## ---

### assumption and stats 
run_pairwise_stats <- function(data,
                               measurements,
                               group_var,
                               group1, group2,
                               stratify_var,
                               id_var = "exp_no") {
  results <- list()
  
  for (m in measurements$measurement) {
    for (d in unique(data$day)) {
      for (level in unique(data[[stratify_var]])) {
        
        df <- data %>%
          filter(day == d,
                 measurement == m,
                 !!sym(group_var) %in% c(group1, group2),
                 !!sym(stratify_var) == level) %>%
          select(all_of(c(id_var, group_var, "value"))) %>%
          drop_na()
        
        if (nrow(df) < 3 || n_distinct(df[[group_var]]) < 2) next
        
        # Assumption checks
        group1_vals <- df %>% filter(!!sym(group_var) == group1) %>% pull(value)
        group2_vals <- df %>% filter(!!sym(group_var) == group2) %>% pull(value)
        
        mean1 <- mean(group1_vals, na.rm = TRUE)
        mean2 <- mean(group2_vals, na.rm = TRUE)
        
        # Skip stat test if both means are below 1%
        if (mean1 < 1 & mean2 < 2) next
        
        norm_g1 <- tryCatch(shapiro.test(group1_vals)$p.value, error = function(e) NA)
        norm_g2 <- tryCatch(shapiro.test(group2_vals)$p.value, error = function(e) NA)
        levene_p <- tryCatch(car::leveneTest(value ~ get(group_var), data = df)[["Pr(>F)"]][1], error = function(e) NA)
        
        equal_var <- if (!is.na(levene_p)) levene_p > 0.05 else NA
        normality <- if (!is.na(norm_g1) & !is.na(norm_g2)) {
          if (norm_g1 > 0.05 & norm_g2 > 0.05) "Normal" else "Not Normal"
        } else {
          NA
        }
        
        df_wide <- df %>%
          pivot_wider(names_from = !!sym(group_var), values_from = value) %>%
          drop_na()
        
        paired_t <- tryCatch(t.test(df_wide[[group1]], df_wide[[group2]], paired = TRUE)$p.value, error = function(e) NA)
        wilcox <- tryCatch(wilcox.test(df_wide[[group1]], df_wide[[group2]], paired = TRUE, exact = FALSE)$p.value, error = function(e) NA)
        unpaired_t <- tryCatch(t.test(value ~ get(group_var), data = df, var.equal = equal_var)$p.value, error = function(e) NA)
        
        # Get direction and effect
        eff_info <- get_effect_and_direction(df, group_var, value_var = "value", group1, group2)
        
        stars <- case_when(
          paired_t <= 0.001 ~ "***",
          paired_t <= 0.01  ~ "**",
          paired_t <= 0.05  ~ "*",
          paired_t <= 0.07  ~ ".",
          TRUE              ~ ""
        )
        
        results[[length(results) + 1]] <- tibble(
          measurement = m,
          day = d,
          group_col = group_var,
          group1 = group1,
          group2 = group2,
          stratify_by = stratify_var,
          stratify_level = level,
          p.value = paired_t,
          wilcox_p = wilcox,
          unpaired_t = unpaired_t,
          test_used = "Paired t-test",
          normality_group1 = norm_g1,
          normality_group2 = norm_g2,
          normality = normality,
          equal_variance = equal_var,
          levene_p = levene_p,
          stars = stars,
          fold_label = eff_info$fold_label,
          global_change = eff_info$global_change,
          effect_size = eff_info$effect_size 
        )
      }
    }
  }
  
  bind_rows(results)
}

multi_flow_test <- function(data, measurements, group_var, subject_id = "exp_no",
                            stratify_var, day_var = "day", value_var = "value",
                            save_dir, comparison) {
  
  required_packages <- c("afex", "emmeans", "broom", "dplyr", "tidyr", "ggplot2", "glue")
  invisible(lapply(required_packages, require, character.only = TRUE))
  
  results <- list()
  
  for (m in measurements$measurement) {
    for (d in unique(data[[day_var]])) {
      for (level in unique(data[[stratify_var]])) {
        
        df <- data %>%
          filter(
            measurement == m,
            !!sym(day_var) == d,
            !!sym(stratify_var) == level
          ) %>%
          select(all_of(c(subject_id, group_var, value_var))) %>%
          drop_na()
        
        means_by_group <- df %>%
          group_by(!!sym(group_var)) %>%
          summarise(mean_val = mean(!!sym(value_var), na.rm = TRUE), .groups = "drop")
        
        if (all(means_by_group$mean_val < 2)) next
        if (nrow(df) < 3 || n_distinct(df[[group_var]]) < 2) next
        
        model <- tryCatch({
          afex::aov_ez(
            id = subject_id,
            dv = value_var,
            data = df,
            within = group_var,
            type = 3,
            fun_aggregate = mean,
            factorize = FALSE,
            return = "afex_aov"
          )
        }, error = function(e) NULL)
        
        if (is.null(model)) next
        
        assumptions_ok <- shapiro.test(residuals(model$lm))$p.value > 0.05
        
        if (assumptions_ok) {
          em_res <- emmeans::emmeans(model, as.formula(paste("pairwise ~", group_var)), adjust = "tukey")
          contrast_df <- broom::tidy(em_res$contrasts) %>%
            separate(contrast, into = c("group1", "group2"), sep = " - ")
          
          effect_info <- purrr::map2_dfr(contrast_df$group1, contrast_df$group2, ~ {
            res <- get_effect_and_direction(df, group_var, value_var, .x, .y)
            tibble(group1 = .x, group2 = .y, !!!res)
          })
          
          df_res <- contrast_df %>%
            mutate(p.value = adj.p.value) %>%
            left_join(effect_info, by = c("group1", "group2"))
        } else {
          df_wide <- df %>%
            pivot_wider(names_from = group_var, values_from = value_var) %>%
            column_to_rownames(var = subject_id)
          
          friedman_result <- friedman.test(as.matrix(df_wide))
          
          df_res <- tibble(
            group1 = "overall", group2 = "overall",
            p.value = friedman_result$p.value
          )
          
          if (friedman_result$p.value <= 0.05) {
            pwc <- pairwise.wilcox.test(df[[value_var]], df[[group_var]],
                                        paired = TRUE, p.adjust.method = "bonferroni")
            pwc_df <- as.data.frame(as.table(pwc$p.value)) %>%
              filter(!is.na(Freq)) %>%
              rename(group1 = Var1, group2 = Var2, p.value = Freq)
            
            effect_info <- purrr::map2_dfr(pwc_df$group1, pwc_df$group2, ~ {
              res <- get_effect_and_direction(df, group_var, value_var, .x, .y)
              tibble(group1 = .x, group2 = .y, !!!res)
            })
            
            pwc_df <- pwc_df %>%
              left_join(effect_info, by = c("group1", "group2"))
            
            df_res <- bind_rows(df_res, pwc_df)
          }
        }
        
        df_res <- df_res %>%
          mutate(
            measurement = m,
            day = d,
            test_used = ifelse(assumptions_ok, "Paired ANOVA + Tukey", "Friedman + Wilcoxon"),
            stratify_by = stratify_var,
            stratify_level = level,
            stars = case_when(
              p.value <= 0.001 ~ "***",
              p.value <= 0.01  ~ "**",
              p.value <= 0.05  ~ "*",
              p.value <= 0.07  ~ ".",
              TRUE ~ ""
            )
          )
        
        results[[length(results) + 1]] <- df_res
      }
    }
  }
  
  final_results <- bind_rows(results)
  
  filename_prefix <- glue::glue("{comparison}_stat_results")
  save_csv(data = final_results, path = file.path(save_dir, paste0(filename_prefix, ".csv")), row.names = FALSE)
  
  return(final_results)
}

## effect size and direction of change 
get_effect_and_direction <- function(df, group_var, value_var = "value",
                                     group1, group2, min_abs_diff = 5) {
  # Pull values for each group
  g1_vals <- df %>% filter(as.character(!!sym(group_var)) == group1) %>% pull(!!sym(value_var))
  g2_vals <- df %>% filter(as.character(!!sym(group_var)) == group2) %>% pull(!!sym(value_var))
  
  g1_mean <- mean(g1_vals, na.rm = TRUE)
  g2_mean <- mean(g2_vals, na.rm = TRUE)
  abs_diff <- abs(g2_mean - g1_mean)
  
  # Direction based on paired differences (optional global trend)
  global_change <- if (length(g1_vals) == length(g2_vals)) {
    diffs <- g2_vals - g1_vals
    if (all(diffs > 0)) "↑"
    else if (all(diffs < 0)) "↓"
    else "↔"
  } else {
    "↔"
  }
  
  # Handle divide-by-zero or missing
  if (g1_mean == 0 || is.na(g1_mean)) {
    return(list(
      effect_size = NA,
      direction = "↔",
      fold_label = "NA",
      global_change = global_change
    ))
  }
  
  raw_fc <- g2_mean / g1_mean
  
  # Default outputs
  direction <- "↔"
  fc_display <- 1
  fold_label <- "1× ↔"
  
  # Suppress small differences even if fold change is high
  if (abs_diff < min_abs_diff) {
    fold_label <- "NS"
  } else if (abs(raw_fc - 1) >= 0.05) {
    if (raw_fc > 1) {
      direction <- "↑"
      fc_display <- round(raw_fc, 2)
      fold_label <- paste0(fc_display, "× ", direction)
    } else {
      direction <- "↓"
      fc_display <- round(1 / raw_fc, 2)
      fold_label <- paste0(fc_display, "× ", direction)
    }
  }
  
  return(list(
    global_change = global_change,
    fold_label = fold_label,
    effect_size = fc_display
  ))
}

#works wth pairwise 
add_pairwise_annotation <- function( plot, stats_data, measurement, group_col, 
                                 stratify_by, stratify_level, 
                                 show_pval = TRUE, show_effect = TRUE) {
  # Filter relevant stats for this measurement + group
  stat_subset <- stats_data %>%
    filter(
      measurement == !!measurement,
      group_col == !!group_col,
      stratify_by == !!stratify_by,
      stratify_level == !!stratify_level,
      !is.na(p.value)
    ) %>%
    mutate(
      p.value = as.numeric(p.value),
      stars = case_when(
        p.value <= 0.001 ~ "***",
        p.value <= 0.01  ~ "**",
        p.value <= 0.05  ~ "*",
        TRUE             ~ ""
      ),
      show_star = p.value <= 0.05,
      show_pval = p.value <= 0.072,
      show_effect = !is.na(effect_size) & effect_size >= 1.1,
      label_pval = paste0("p = ", formatC(p.value, format = "f", digits = 3)),
      label_effect = fold_label
    )
  
  # Early exit if no stars, pvals, or effects to show
  if (nrow(stat_subset %>% filter(show_star | show_pval | show_effect)) == 0) return(plot)
  
  cat(glue::glue("📊 {measurement} | {stratify_by}: {stratify_level} → {nrow(stat_subset)} annotated\n"))
  
  # X-axis positions
  x_levels <- levels(factor(plot$data$day))
  
  stat_subset <- stat_subset %>%
    mutate(
      x_day = as.numeric(factor(day, levels = x_levels)),
      x_offset = 0.2,
      x_start = x_day - x_offset,
      x_end   = x_day + x_offset,
      stars_y = 102,
      bar_y   = 100,
      label_y = 97,
      effect_y = 90
    )
  
  # Add test used to caption
  tests_used_all <- stats_data %>%
    filter(
      measurement == !!measurement,
      group_col == !!group_col,
      stratify_by == !!stratify_by,
      stratify_level == !!stratify_level
    ) %>%
    pull(test_used) %>%
    unique() %>%
    discard(is.na)
  
  # Count stars in filtered subset
  has_stars <- any(stat_subset$show_star, na.rm = TRUE)
  
  # Add appropriate caption
  if (length(tests_used_all) == 1) {
    caption_text <- if (has_stars) {
      paste0("Test used: ", tests_used_all)
    } else {
      paste0("Test used: ", tests_used_all, " — no significant differences detected")
    }
    plot <- plot + labs(caption = paste0(plot$labels$caption, "\n", caption_text))
  }
  
  p_out <- plot
  
  # add star + bar if p ≤ 0.05
  if (any(stat_subset$show_star)) {
    p_out <- p_out +
      geom_segment(
        data = stat_subset %>% filter(show_star),
        aes(x = x_start, xend = x_end, y = bar_y, yend = bar_y),
        inherit.aes = FALSE,
        linewidth = 0.7,
        color = "black"
      ) +
      geom_text(
        data = stat_subset %>% filter(show_star),
        aes(x = x_day, y = stars_y, label = stars),
        inherit.aes = FALSE,
        size = 10,
        fontface = "bold"
      )
  }
  
  # Add p-values if toggled on
  if (show_pval && any(stat_subset$show_pval)) {
    p_out <- p_out +
      geom_text(
        data = stat_subset %>% filter(show_pval),
        aes(x = x_day, y = label_y, label = label_pval),
        inherit.aes = FALSE,
        size = 3,
        fontface = "italic"
      )
  }
  
  # Add effect size and direction if toggled on
  if (show_effect && any(stat_subset$show_effect)) {
    p_out <- p_out +
      geom_text(
        data = stat_subset %>% filter(show_effect),
        aes(x = x_day, y = effect_y, label = label_effect),
        inherit.aes = FALSE,
        size = 3,
        fontface = "plain",
        family = "Arial Unicode MS"
      )
  }
  
  return(p_out)
}


## staggers yes 
add_multi_annotation <- function( plot, stats_data, measurement, group_col, 
                                 stratify_by, stratify_level, spacing_factor = 0.08,
                                 show_pval = TRUE, show_effect = TRUE) {
  # Filter relevant stats for this measurement + group
  stat_subset <- stats_data %>%
    filter(
      measurement == !!measurement,
      group_col == !!group_col,
      stratify_by == !!stratify_by,
      stratify_level == !!stratify_level,
      !is.na(p.value),
      !(group1 == "overall" & group2 == "overall")
    ) %>%
    mutate(
      p.value = as.numeric(p.value),
      stars = case_when(
        p.value <= 0.001 ~ "***",
        p.value <= 0.01  ~ "**",
        p.value <= 0.05  ~ "*",
        TRUE             ~ ""
      ),
      show_star = p.value <= 0.05,
      show_pval = p.value <= 0.072,
      show_effect = !is.na(effect_size) & effect_size >= 1.1,
      label_pval = paste0("p = ", formatC(p.value, format = "f", digits = 3)),
      label_effect = fold_label
    )
  
  # Early exit if no stars, pvals, or effects to show
  if (nrow(stat_subset %>% filter(show_star | show_pval | show_effect)) == 0) return(plot)
  
  n_stars   <- sum(stat_subset$show_star, na.rm = TRUE)
  n_pval    <- if (show_pval)   sum(stat_subset$show_pval, na.rm = TRUE) else 0
  n_effects <- if (show_effect) sum(stat_subset$show_effect, na.rm = TRUE) else 0
  
  cat(glue::glue(
    "🔎 {measurement} | {stratify_by}: {stratify_level} → ",
    "{n_stars} stars, {n_pval} p-values, {n_effects} effects\n"
  ))
  
  
  # X-axis positions
  dodge_width <- 0.9
  n_treatments <- length(unique(plot$data$treatment))
  
  # Get positions of each treatment within each day (for dodged bars)
  x_pos_lookup <- plot$data %>%
    distinct(day, treatment) %>%
    arrange(day, treatment) %>%
    mutate(
      day_num = as.numeric(factor(day)),
      treatment_num = as.numeric(factor(treatment)),
      xpos = day_num + ((treatment_num - 1) / (n_treatments - 1)) * dodge_width - (dodge_width / 2)
    )
  
  # Get max y-axis tick to scale spacing
  y_max_tick <- max(ggplot_build(plot)$layout$panel_params[[1]]$y$get_breaks(), na.rm = TRUE)
  spacing <- y_max_tick * spacing_factor
  
  # Join in to your stat_subset
  stat_subset <- stat_subset %>%
    left_join(x_pos_lookup %>% rename(group1 = treatment, x_start = xpos), by = c("day", "group1")) %>%
    left_join(x_pos_lookup %>% rename(group2 = treatment, x_end   = xpos), by = c("day", "group2")) %>%
    mutate(x_day = (x_start + x_end) / 2)
  
  # Apply staggering only to stars
  star_y_lookup <- stat_subset %>%
    filter(show_star) %>%
    arrange(day, x_start, x_end) %>%
    group_by(day) %>%
    mutate(
      row_num = row_number(),
      bar_y = y_max_tick - spacing * (row_num - 1),
      stars_y = bar_y + spacing * 0.2
    ) %>%
    ungroup() %>%
    select(day, group1, group2, bar_y, stars_y)
  
  # Merge back staggered y-positions
  stat_subset <- stat_subset %>%
    left_join(star_y_lookup, by = c("day", "group1", "group2")) %>%
    mutate(
      bar_y   = coalesce(bar_y, y_max_tick),
      stars_y = coalesce(stars_y, bar_y + spacing * 0.2),
      label_y = bar_y - spacing * 0.25,
      effect_y = bar_y - spacing * 0.9,
      
      # shorten bar width relative to center
      bar_shrink = 0.2,  # try 0.1–0.2 for subtle shortening
      x_start_adj = x_day - ((x_day - x_start) * (1 - bar_shrink)),
      x_end_adj   = x_day + ((x_end - x_day) * (1 - bar_shrink))
    )
  
  # Add test used to caption (if only one kind used)
  tests_used <- unique(stat_subset$test_used)
  if (length(tests_used) == 1 && !is.na(tests_used)) {
    plot <- plot +
      labs(caption = paste0(plot$labels$caption, "\nTest used: ", tests_used))
  }
  
  p_out <- plot
  
  # add star + bar if p ≤ 0.05
  if (any(stat_subset$show_star)) {
    p_out <- p_out +
      geom_segment(
        data = stat_subset %>% filter(show_star),
        aes(x = x_start_adj, xend = x_end_adj, y = bar_y, yend = bar_y),
        inherit.aes = FALSE,
        linewidth = 0.7,
        color = "black"
      ) +
      geom_text(
        data = stat_subset %>% filter(show_star),
        aes(x = x_day, y = stars_y, label = stars),
        inherit.aes = FALSE,
        size = 10,
        fontface = "bold"
      )
  }
  
  # Add p-values if toggled on
  if (show_pval && any(stat_subset$show_pval)) {
    p_out <- p_out +
      geom_text(
        data = stat_subset %>% filter(show_pval),
        aes(x = x_day, y = label_y, label = label_pval),
        inherit.aes = FALSE,
        size = 3,
        fontface = "italic"
      )
  }
  
  # Add effect size and direction if toggled on
  if (show_effect && any(stat_subset$show_effect)) {
    p_out <- p_out +
      geom_text(
        data = stat_subset %>% filter(show_effect),
        aes(x = x_day, y = effect_y, label = label_effect),
        inherit.aes = FALSE,
        size = 4,
        fontface = "plain"
      )
  }
  
  return(p_out)
}

add_friedman_annotation <- function(plot, stats_data, measurement, stratify_by, stratify_level,
                                    symbol = "†", spacing_factor = 0.05) {
  # Filter for Friedman-only significant results
  stat_subset <- stats_data %>%
    filter(
      measurement == !!measurement,
      stratify_by == !!stratify_by,
      stratify_level == !!stratify_level,
      group1 == "overall", group2 == "overall",
      !is.na(p.value), p.value <= 0.05
    )
  
  if (nrow(stat_subset) == 0) return(plot)
  
  # Get the x-axis levels from the plot (they are treated as discrete)
  x_levels <- plot$data %>%
    distinct(day) %>%
    arrange(day) %>%
    pull(day) %>%
    as.character()
  
  stat_subset <- stat_subset %>%
    mutate(day_char = as.character(day),
           x_pos = match(day_char, x_levels))  # Match day value to its axis position
  
  # Compute y position just above max tick
  y_max_tick <- max(ggplot_build(plot)$layout$panel_params[[1]]$y$get_breaks(), na.rm = TRUE)
  spacing <- y_max_tick * spacing_factor
  
  # Add annotation
  plot +
    geom_text(
      data = stat_subset,
      aes(x = x_pos, y = y_max_tick + spacing, label = symbol),
      inherit.aes = FALSE,
      size = 6,
      fontface = "italic"
    ) +
    labs(caption = paste0(plot$labels$caption, "\n", symbol, " = Friedman significant; pairwise inconclusive"))
}

## ---
## PLOT THEME ----
## ---

theme_pub <- function(base_size = 14, base_family = "helvetica_neue") {
  library(grid)
  library(ggthemes)
  
  theme_foundation(base_size = base_size, base_family = base_family) +
    theme(
      plot.title = element_text(
        face = "bold", size = 20, hjust = 0.5,
        margin = margin(10, 10, 20,  10)
      ),
      text = element_text(),
      panel.background = element_rect(colour = NA, fill = NA), #makes background transparent
      plot.background = element_rect(colour = NA, fill = NA),
      panel.border = element_rect(colour = NA),
      
      # Axis titles and text sizes (your additions)
      axis.title.x = element_text(vjust = -0.2, size = 18),
      axis.title.y = element_text(angle = 90, vjust = 2, size = 18),
      axis.text.x = element_text(size = 18, colour = "black"),
      axis.text.y = element_text(size = 16, colour = "black"),
      
      axis.line.x = element_line(colour = "black"),
      axis.line.y = element_line(colour = "black"),
      axis.ticks = element_line(),
      
      #panel.grid.major = element_line(colour = "#f0f0f0"), #if you want grid lines
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      
      legend.key = element_rect(colour = NA),
      legend.position = "top",
      legend.direction = "horizontal",
      legend.box = "horizontal",
      legend.key.size = unit(0.5, "cm"),
      legend.title = element_text(face = "bold", size = 16),
      legend.text = element_text(size = 16),
      
      #caption settings
      plot.caption = element_text(
        size = 10,
        hjust = 0,             # left-aligned
        vjust = 1.5,           # shift it down
        margin = margin(t = 15, b = 0)
      ),
      plot.caption.position = "plot",

      # for facet labelts 
      strip.background = element_rect(colour = "#f0f0f0", fill = "#f0f0f0"),
      strip.text = element_text(size = 14, face = "bold")
    )
}




## ASSUMPTIONS
check_assumptions_by_group <- function(data,
                                       measurements,
                                       group_col,    # "treatment" or "genotype"
                                       group1,
                                       group2,
                                       within_col,   # "genotype" or "treatment"
                                       id_col = "exp_no") {
  results <- list()
  
  for (measure in measurements$measurement) {
    for (day_val in unique(data$day)) {
      for (within_val in unique(data[[within_col]])) {
        
        filtered <- data %>%
          filter(
            day == day_val,
            !!sym(within_col) == within_val,
            !!sym(group_col) %in% c(group1, group2)
          ) %>%
          select(all_of(c(id_col, group_col)), !!sym(measure)) %>%
          drop_na()
        
        if (nrow(filtered) < 3 || n_distinct(filtered[[group_col]]) < 2) next
        
        group1_vals <- filtered %>% filter(!!sym(group_col) == group1) %>% pull(!!sym(measure))
        group2_vals <- filtered %>% filter(!!sym(group_col) == group2) %>% pull(!!sym(measure))
        
        shapiro1 <- tryCatch(shapiro.test(group1_vals)$p.value, error = function(e) NA)
        shapiro2 <- tryCatch(shapiro.test(group2_vals)$p.value, error = function(e) NA)
        normality <- ifelse(shapiro1 > 0.05 & shapiro2 > 0.05, "Normal", "Not Normal")
        
        levene_p <- tryCatch(
          car::leveneTest(as.formula(paste0(measure, " ~ ", group_col)), data = filtered)[["Pr(>F)"]][1],
          error = function(e) NA
        )
        var_equal <- ifelse(levene_p > 0.05, TRUE, FALSE)
        
        results[[length(results) + 1]] <- tibble(
          day = day_val,
          measurement = measure,
          group_col = group_col,
          group1 = group1,
          group2 = group2,
          stratify_by = within_col,
          stratify_level = within_val,
          shapiro_p_group1 = shapiro1,
          shapiro_p_group2 = shapiro2,
          normality = normality,
          levene_p = levene_p,
          equal_variance = var_equal
        )
      }
    }
  }
  
  bind_rows(results)
}


## ---
## SIGNIFICANCE TESTS 
## ---

## two way comparison
run_pairwise_stats_by_group <- function(data, 
                                        measurements, 
                                        group_col,       # what you are comparing - genotype/treatment
                                        group1, group2,  # what they are e.g., "CTRL", "RA" or "D21", "T21"
                                        within_col,      # within what are you comparing e.g., "genotype" or "treatment"
                                        id_col = "exp_no") {
  results <- list()
  
  for (measure in measurements$measurement) {
    for (day_val in unique(data$day)) {
      for (within_val in unique(data[[within_col]])) {
        
        # Filter the relevant subset
        filtered_data <- data %>%
          filter(
            day == day_val,
            !!sym(within_col) == within_val,
            !!sym(group_col) %in% c(group1, group2)
          ) %>%
          select(all_of(c(id_col, group_col)), !!sym(measure)) %>%
          drop_na()
        
        message(glue::glue("→ Day: {day_val} | {within_col}: {within_val} | N rows: {nrow(filtered_data)} | Unique {group_col}s: {length(unique(filtered_data[[group_col]]))}"))
        
        if (length(unique(filtered_data[[group_col]])) < 2) next
        
        # Convert to wide for paired t-test
        data_wide <- filtered_data %>%
          pivot_wider(names_from = all_of(group_col), values_from = all_of(measure))
        
        # Run tests
        t_test <- tryCatch(
          t.test(as.formula(paste0(measure, " ~ ", group_col)), data = filtered_data, var.equal = FALSE)$p.value,
          error = function(e) NA
        )
        
        wilcox <- tryCatch(
          wilcox.test(as.formula(paste0(measure, " ~ ", group_col)), data = filtered_data)$p.value,
          error = function(e) NA
        )
        
        kruskal <- tryCatch(
          kruskal.test(as.formula(paste0(measure, " ~ ", group_col)), data = filtered_data)$p.value,
          error = function(e) NA
        )
        
        paired_t <- tryCatch({
          if (all(c(group1, group2) %in% names(data_wide))) {
            t.test(data_wide[[group1]], data_wide[[group2]], paired = TRUE)$p.value
          } else {
            NA
          }
        }, error = function(e) NA)
        
        results[[length(results) + 1]] <- tibble(
          measurement = measure,
          day = day_val,
          group_col = group_col,
          group1 = group1,
          group2 = group2,
          stratify_by = within_col,
          stratify_level = within_val,
          t_test_p = t_test,
          wilcox_p = wilcox,
          kruskal_p = kruskal,
          paired_t_p = paired_t
        )
      }
    }
  }
  
  bind_rows(results)
}

run_pairwise_stats_by_group <- function(data, 
                                        measurements, 
                                        group_col,       # "treatment" or "genotype"
                                        group1, group2,  # e.g., "CTRL", "RA" or "D21", "T21"
                                        within_col,      # e.g., "genotype" or "treatment"
                                        id_col = "exp_no") {
  results <- list()
  
  for (measure in measurements$measurement) {
    for (day_val in unique(data$day)) {
      for (within_val in unique(data[[within_col]])) {
        
        # Filter relevant subset
        filtered_data <- data %>%
          filter(
            day == day_val,
            !!sym(within_col) == within_val,
            !!sym(group_col) %in% c(group1, group2)
          ) %>%
          select(all_of(c(id_col, group_col)), !!sym(measure)) %>%
          drop_na()
        
        message(glue::glue("→ Day: {day_val} | {within_col}: {within_val} | N rows: {nrow(filtered_data)} | Unique {group_col}s: {length(unique(filtered_data[[group_col]]))}"))
        
        if (length(unique(filtered_data[[group_col]])) < 2) next
        
        # Convert to wide for paired t-test
        data_wide <- filtered_data %>%
          pivot_wider(names_from = all_of(group_col), values_from = all_of(measure))
        
        # Run tests
        t_test <- tryCatch(
          t.test(as.formula(paste0(measure, " ~ ", group_col)), data = filtered_data, var.equal = FALSE)$p.value,
          error = function(e) NA
        )
        
        wilcox <- tryCatch(
          wilcox.test(as.formula(paste0(measure, " ~ ", group_col)), data = filtered_data, paired = TRUE)$p.value,
          error = function(e) NA
        )
        
        paired_t <- tryCatch({
          if (all(c(group1, group2) %in% names(data_wide))) {
            t.test(data_wide[[group1]], data_wide[[group2]], paired = TRUE)$p.value
          } else {
            NA
          }
        }, error = function(e) NA)
        
        # Store result
        results[[length(results) + 1]] <- tibble(
          measurement = measure,
          day = day_val,
          group_col = group_col,
          group1 = group1,
          group2 = group2,
          stratify_by = within_col,
          stratify_level = within_val,
          t_test_p = t_test,
          wilcox_p = wilcox,
          paired_t_p = paired_t
        )
      }
    }
  }
  
  bind_rows(results)
}

## ---
## PREFERED STAT ANNOTATION ----
## ---
# will default to paired t if reps lessthan/equal to 5 
annotate_preferred_tests <- function(stat_df, assumption_df,
                                     group_col = "group_col",
                                     group1_col = "group1",
                                     group2_col = "group2",
                                     stratify_by_col = "stratify_by",
                                     stratify_level_col = "stratify_level",
                                     default_to_paired_if_small_n = TRUE,
                                     n_replicates = 3,
                                     small_n_threshold = 5) {
  
  # Ensure assumption_df is distinct
  assumption_df <- assumption_df %>%
    select(
      day, measurement,
      !!sym(group_col), !!sym(group1_col), !!sym(group2_col),
      !!sym(stratify_by_col), !!sym(stratify_level_col),
      normality, equal_variance
    ) %>%
    distinct()
  
  stat_df %>%
    left_join(
      assumption_df,
      by = c(
        "day" = "day",
        "measurement" = "measurement",
        group_col = group_col,
        "group1" = group1_col,
        "group2" = group2_col,
        "stratify_by" = stratify_by_col,
        "stratify_level" = stratify_level_col
      )
    ) %>%
    mutate(
      preferred_test = case_when(
        default_to_paired_if_small_n & n_replicates <= small_n_threshold ~ "paired_t",
        normality == "Normal" & equal_variance ~ "paired_t",
        normality != "Normal" | equal_variance == FALSE ~ "wilcox",
        TRUE ~ NA_character_
      ),
      preferred_p = case_when(
        preferred_test == "paired_t" ~ paired_t_p,
        preferred_test == "wilcox" ~ wilcox_p,
        TRUE ~ NA_real_
      ),
      significant = ifelse(!is.na(preferred_p) & preferred_p < 0.05, TRUE, FALSE)
    )
}


## ---
## PLOT STATS ----
## ---

add_stat_annotation_multiple_days <- function(plot, stat_data, measurement, group_col, stratify_by, stratify_level) {
  # Filter to all relevant stats for this measure + group
  all_stats <- stat_data %>%
    filter(
      measurement == !!measurement,
      group_col == !!group_col,
      stratify_by == !!stratify_by,
      stratify_level == !!stratify_level,
      !is.na(preferred_p)
    )
  
  # Check if any relevant stats exist
  if (nrow(all_stats) == 0) {
    cat(glue::glue("❌ No stats found for: {measurement} | {stratify_by}: {stratify_level}"))
    return(plot)
  }
  
  # Now filter to those we want to annotate
  stat_subset <- all_stats %>%
    filter(preferred_p <= 0.07)
  
  cat(glue::glue("📊 Found {nrow(all_stats)} stats for: {measurement} | {stratify_by}: {stratify_level} ({nrow(stat_subset)} with p ≤ 0.07)"))
  
  if (nrow(stat_subset) == 0) return(plot)
  
  # Add labels, bar positions
  stat_subset <- stat_subset %>%
    mutate(
      stars = case_when(
        preferred_p <= 0.001 ~ "***",
        preferred_p <= 0.01  ~ "**",
        preferred_p <= 0.05  ~ "*",
        TRUE                 ~ ""
      ),
      label_pval = paste0("p = ", formatC(preferred_p, format = "f", digits = 3)),
      label_y = 94,
      stars_y = 100,
      bar_y = 98,
      
      x_day = as.numeric(factor(day, levels = sort(unique(plot$data$day)))),
      x_offset = 0.2,
      x_start = x_day - x_offset,
      x_end   = x_day + x_offset
    )
  
  # Add annotations
  plot +
    geom_text(  # Stars
      data = stat_subset,
      aes(x = x_day, y = stars_y, label = stars),
      inherit.aes = FALSE,
      size = 7,
      fontface = "bold"
    ) +
    geom_text(  # P-values (italic)
      data = stat_subset,
      aes(x = x_day, y = label_y, label = label_pval),
      inherit.aes = FALSE,
      size = 2,
      fontface = "italic"
    ) +
    geom_segment(  # Bar (only if star exists)
      data = stat_subset %>% filter(stars != ""),
      aes(x = x_start, xend = x_end, y = bar_y, yend = bar_y),
      inherit.aes = FALSE,
      linewidth = 0.7,
      color = "black"
    )
}



# # ---------  OTHER GGTHEME THEMES STUFF --------
# #https://medium.com/analytics-vidhya/ggplot2-themes-for-publication-ready-plots-including-dark-themes-9cd65cc5a7e3
# 
# scale_fill_Publication <- function(...){
#    library(scales)
#    discrete_scale("fill","Publication",manual_pal(values = c("#386cb0","#f87f01","#7fc97f","#ef3b2c","#feca01","#a6cee3","#fb9a99","#984ea3","#8C591D")), ...)
# 
# }
# 
# scale_colour_Publication <- function(...){
#   library(scales)
#   discrete_scale("colour","Publication",manual_pal(values = c("#386cb0","#f87f01","#7fc97f","#ef3b2c","#feca01","#a6cee3","#fb9a99","#984ea3","#8C591D")), ...)
# 
# }
# 
# ### Dark theme for ggplot plots
# theme_dark_grey <- function(base_size=14, base_family="sans") {
#   library(grid)
#   library(ggthemes)
#   (theme_foundation(base_size=base_size, base_family=base_family)
#     + theme(plot.title = element_text(face = "bold", colour = '#ffffb3',
#                                       size = rel(1.2), hjust = 0.5, margin = margin(0,0,20,0)),
#             text = element_text(),
#             panel.background = element_rect(colour = NA, fill = 'grey20'),
#             plot.background = element_rect(colour = NA, fill = '#262626'),
#             panel.border = element_rect(colour = NA),
#             axis.title = element_text(face = "bold",size = rel(1), colour = 'white'),
#             axis.title.y = element_text(angle=90,vjust =2),
#             axis.title.x = element_text(vjust = -0.2),
#             axis.text = element_text(colour = 'grey70'),
#             axis.line.x = element_line(colour="grey70"),
#             axis.line.y = element_line(colour="grey70"),
#             axis.ticks = element_line(colour="grey70"),
#             panel.grid.major = element_line(colour="#262626"),
#             panel.grid.minor = element_blank(),
#             legend.background = element_rect(fill ='#262626'),
#             legend.text = element_text(color = 'white'),
#             legend.key = element_rect(colour = NA, fill = '#262626'),
#             legend.position = "bottom",
#             legend.direction = "horizontal",
#             legend.box = "vetical",
#             legend.key.size= unit(0.5, "cm"),
#             #legend.margin = unit(0, "cm"),
#             legend.title = element_text(face="italic", colour = 'white'),
#             plot.margin=unit(c(10,5,5,5),"mm"),
#             strip.background=element_rect(colour="#2D3A4C",fill="#2D3A4C"),
#             strip.text = element_text(face="bold", colour = 'white')
#     ))
# }
# 
# scale_fill_Publication_dark <- function(...){
#   library(scales)
#   discrete_scale("fill","Publication",manual_pal(values = c("#fbb4ae","#b3cde3","#ccebc5","#decbe4","#fed9a6","#ffffcc","#e5d8bd","#fddaec","#f2f2f2")), ...)
# 
# }
# 
# scale_colour_Publication_dark <- function(...){
#   library(scales)
#   discrete_scale("colour","Publication",manual_pal(values = c("#fbb4ae","#b3cde3","#ccebc5","#decbe4","#fed9a6","#ffffcc","#e5d8bd","#fddaec","#f2f2f2")), ...)
# 
# }
# 
# theme_transparent <- function(base_size=14, base_family="sans") {
#     library(grid)
#     library(ggthemes)
#    (theme_foundation(base_size=base_size, base_family=base_family)
#       + theme(plot.title = element_text(face = "bold", colour = '#ffffb3',
#                                         size = rel(1.2), hjust = 0.5),
#               text = element_text(),
#               panel.background = element_rect(colour = NA, fill = 'transparent'),
#               plot.background = element_rect(colour = NA, fill = 'transparent'),
#               panel.border = element_rect(colour = NA),
#               axis.title = element_text(face = "bold",size = rel(1), colour = 'white'),
#               axis.title.y = element_text(angle=90,vjust =2),
#               axis.title.x = element_text(vjust = -0.2),
#               axis.text = element_text(colour = 'grey70'),
#               axis.line.x = element_line(colour="grey70"),
#               axis.line.y = element_line(colour="grey70"),
#               axis.ticks = element_line(colour="grey70"),
#               panel.grid.major = element_line(colour="#262626"),
#               panel.grid.minor = element_blank(),
#               legend.background = element_rect(fill = 'transparent'),
#               legend.text = element_text(color = 'white'),
#               legend.key = element_rect(colour = NA, fill = 'grey20'),
#               legend.position = "bottom",
#               legend.direction = "horizontal",
#               legend.box = "vetical",
#               legend.key.size= unit(0.5, "cm"),
#               #legend.margin = unit(0, "cm"),
#               legend.title = element_text(face="italic", colour = 'white'),
#               plot.margin=unit(c(10,5,5,5),"mm"),
#               strip.background=element_rect(colour="#2D3A4C",fill="#2D3A4C"),
#               strip.text = element_text(face="bold", colour = 'white')
#       ))
# }
# 
# theme_dark_blue <- function(base_size=14, base_family="sans") {
#   library(grid)
#   library(ggthemes)
#   (theme_foundation(base_size=base_size, base_family=base_family)
#     + theme(plot.title = element_text(face = "bold", colour = '#ffffb3',
#                                       size = rel(1.2), hjust = 0.5, margin = margin(0,0,20,0)),
#             text = element_text(),
#             panel.background = element_rect(colour = NA, fill = '#282C33'),
#             plot.background = element_rect(colour = NA, fill = '#282C33'),
#             panel.border = element_rect(colour = NA),
#             axis.title = element_text(face = "bold",size = rel(1), colour = 'white'),
#             axis.title.y = element_text(angle=90,vjust =2),
#             axis.title.x = element_text(vjust = -0.2),
#             axis.text = element_text(colour = 'grey70'),
#             axis.line.x = element_line(colour="grey70"),
#             axis.line.y = element_line(colour="grey70"),
#             axis.ticks = element_line(colour="grey70"),
#             panel.grid.major = element_line(colour="#343840"),
#             panel.grid.minor = element_blank(),
#             legend.background = element_rect(fill ='#282C33'),
#             legend.text = element_text(color = 'white'),
#             legend.key = element_rect(colour = NA, fill = '#282C33'),
#             legend.position = "bottom",
#             legend.direction = "horizontal",
#             legend.box = "vetical",
#             legend.key.size= unit(0.5, "cm"),
#             #legend.margin = unit(0, "cm"),
#             legend.title = element_text(face="italic", colour = 'white'),
#             plot.margin=unit(c(10,5,5,5),"mm"),
#             strip.background=element_rect(colour="#2D3A4C",fill="#2D3A4C"),
#             strip.text = element_text(face="bold", colour = 'white')
#     ))
#  }