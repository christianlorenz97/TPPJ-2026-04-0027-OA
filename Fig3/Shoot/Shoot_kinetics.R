# ==========================================================
# Scatterplot (Mean ± SE) + p-value symbols for:
# shoot_side_area, shoot_height, shootA_rootA
# + CSV "Mean ± SE" for each index w/ Tukey (if p<0.05)
# ==========================================================

# Packages
  library(dplyr)
  library(readr)
  library(lubridate)
  library(tidyr)
  library(ggplot2)
  library(agricolae)
  library(purrr)
  library(patchwork)

# File load
file_path <- "PS141_index.txt"

# Helper
se_fun <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) <= 1) return(NA_real_)
  sd(x) / sqrt(length(x))
}

# p -> symbol, only *, **, *** (ns only for non significant)
p_to_symbol <- function(p) {
  dplyr::case_when(
    is.na(p)  ~ "",
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    TRUE      ~ "ns"
  )
}

# Import
df <- read_tsv(
  file_path,
  show_col_types = FALSE,
  locale = locale(decimal_mark = ".", date_names = "it")
) %>%
  mutate(
    date      = dmy(date),
    n_day     = as.integer(n_day),
    treatment = droplevels(as.factor(treatment))
  )

# Index selection
vars_plot <- c("shoot_side_area", "shoot_height", "shootA_rootA")
ylab_map <- c(
  shoot_side_area = "Shoot Side Area [mm²]",
  shoot_height    = "Shoot Height [mm]",
  shootA_rootA    = "Root to Shoot Area Ratio"
)

# calculate Mean ± SE w/ long
long_dat <- df %>%
  select(n_day, treatment, all_of(vars_plot)) %>%
  pivot_longer(cols = all_of(vars_plot), names_to = "metric", values_to = "value")

# (Mean ± SE) per day x treatment
sum_dat <- long_dat %>%
  group_by(metric, treatment, n_day) %>%
  summarise(
    n    = sum(is.finite(value)),
    mean = mean(value, na.rm = TRUE),
    sd   = sd(value,   na.rm = TRUE),
    se   = ifelse(n > 1, sd / sqrt(n), NA_real_),
    .groups = "drop"
  )

# ANOVA one-way per day (value ~ treatment)
aov_dat <- long_dat %>%
  filter(is.finite(value)) %>%
  group_by(metric, n_day) %>%
  summarise(
    p = {
      ok <- n_distinct(treatment) >= 2 && sum(is.finite(value)) >= 3
      if (ok) {
        pval <- tryCatch({
          fit <- stats::aov(value ~ treatment)
          as.numeric(summary(fit)[[1]][["Pr(>F)"]][1])
        }, error = function(e) NA_real_)
        pval
      } else NA_real_
    },
    .groups = "drop"
  ) %>%
  mutate(p_symbol = p_to_symbol(p))

# Tukey (HSD.test) per day x index (if p<0.05)
tukey_letters <- long_dat %>%
  filter(is.finite(value)) %>%
  group_nest(metric, n_day) %>%
  mutate(
    res = map(data, ~{
      dd <- .x
      dd$treatment <- droplevels(dd$treatment)
      if (n_distinct(dd$treatment) < 2 || sum(is.finite(dd$value)) < 3) {
        return(tibble(treatment = levels(dd$treatment), letter = NA_character_))
      }
      p_val <- tryCatch({
        fit0 <- aov(value ~ treatment, data = dd)
        as.numeric(summary(fit0)[[1]][["Pr(>F)"]][1])
      }, error = function(e) NA_real_)
      if (!is.finite(p_val) || p_val >= 0.05) {
        return(tibble(treatment = levels(dd$treatment), letter = NA_character_))
      }
      hs <- tryCatch({
        fit <- aov(value ~ treatment, data = dd)
        agricolae::HSD.test(fit, "treatment", group = TRUE, console = FALSE)
      }, error = function(e) NULL)
      if (is.null(hs) || is.null(hs$groups)) {
        return(tibble(treatment = levels(dd$treatment), letter = NA_character_))
      }
      out <- hs$groups
      out$treatment <- rownames(out)
      tibble(treatment = factor(out$treatment, levels = levels(dd$treatment)),
             letter    = as.character(out$groups))
    })
  ) %>%
  select(-data) %>%
  unnest(res)

# P-value symbol labelling
tops <- sum_dat %>%
  group_by(metric, n_day) %>%
  summarise(
    y_top = max(mean + se, na.rm = TRUE),
    y_min = min(mean - se, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    y_range = ifelse(
      is.finite(y_top - y_min) & (y_top - y_min) > 0,
      y_top - y_min,
      0.05 * pmax(y_top, 1, na.rm = TRUE)
    ),
    y_ann = y_top + 0.05 * y_range
  )

ann <- aov_dat %>%
  left_join(tops, by = c("metric","n_day")) %>%
  filter(is.finite(y_ann))

# Add Tukey letters to sum_dat for CSV
sum_with_letters <- sum_dat %>%
  left_join(tukey_letters, by = c("metric","n_day","treatment"))

# Output directory
dir.create("plots_selected", showWarnings = FALSE)
dir.create("tables_mean_se", showWarnings = FALSE)

# Graph plotting
plot_fun <- function(df_met, ann_met, ylab_txt) {
  brks <- sort(unique(df_met$n_day))
  ggplot(df_met, aes(x = n_day, y = mean, color = treatment, group = treatment)) +
    geom_line(linewidth = 0.7, alpha = 0.9) +
    geom_point(size = 2) +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0) +
    geom_text(
      data = ann_met %>% filter(p_symbol != ""),
      aes(x = n_day, y = y_ann, label = p_symbol),
      inherit.aes = FALSE, vjust = 0, size = 3.8, color = "black"
    ) +
    scale_x_continuous(breaks = brks, expand = expansion(mult = c(0.02, 0.05))) +
    scale_color_manual(
      values = c(
        "WW_MIXR" = "#117733",
        "WW_MIX1" = "#44AA99",
        "WW_NOI"  = "#332288",
        "NW_MIXR" = "#AA4499",
        "NW_MIX1" = "#CC6677",
        "NW_NOI"  = "#882255"
      )
    ) +
    labs(x = "Day", y = ylab_txt, color = "Treatment") +
    theme_bw(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title       = element_blank()
    )
}

# CSV "Mean ± SE" w/ Tukey letters
write_mean_se_csv <- function(df_met_letters, metric_name) {
  days_all <- sort(unique(df_met_letters$n_day))
  df_wide <- df_met_letters %>%
    mutate(
      mean  = round(mean, 3),
      se    = round(se, 3),
      apice = ifelse(is.na(letter) | !nzchar(letter), "", paste0("^", letter)),
      cell  = ifelse(
        is.finite(mean),
        ifelse(is.finite(se), paste0(mean, " ± ", se, apice), as.character(mean)),
        ""
      )
    ) %>%
    select(treatment, n_day, cell) %>%
    mutate(day_col = paste0("Day ", n_day)) %>%
    select(-n_day) %>%
    tidyr::pivot_wider(names_from = day_col, values_from = cell) %>%
    arrange(treatment) %>%
    rename(Treatments = treatment)
  
  # Column order: Treatments, Day 1 ... Day N
  day_cols <- paste0("Day ", days_all)
  df_wide <- df_wide %>% select(any_of(c("Treatments", day_cols)))
  
  out_csv <- file.path("tables_mean_se", paste0("mean_se_", metric_name, ".csv"))
  readr::write_csv(df_wide, out_csv)
  message("CSV created: ", out_csv)
}

plots_list <- list()  # to build 3-panel figure

# Repeat on the indexes
for (met in vars_plot) {
  df_met_plot  <- sum_dat %>% filter(metric == met)
  ann_met_plot <- ann     %>% filter(metric == met)
  df_met_csv   <- sum_with_letters %>% filter(metric == met)
  ylab         <- unname(ylab_map[[met]])
  
  # Plot
  p <- plot_fun(df_met_plot, ann_met_plot, ylab_txt = ylab)
  plots_list[[met]] <- p  #save for combined figure
  print(p)
  
  # Save (PNG + PDF)
  ggsave(
    file.path("plots_selected", paste0("scatter_", met, ".png")),
    plot   = p,
    width  = 7,
    height = 5,
    dpi    = 300
  )
  ggsave(
    file.path("plots_selected", paste0("scatter_", met, ".pdf")),
    plot   = p,
    width  = 7,
    height = 5,
    device = cairo_pdf
  )
  
  # CSV Mean ± SE w/ Tukey
  write_mean_se_csv(df_met_csv, metric_name = met)
}

# Combined figure
if (all(c("shoot_height", "shoot_side_area", "shootA_rootA") %in% names(plots_list))) {
  
  p_combined <- plots_list[["shoot_height"]] +
    plots_list[["shoot_side_area"]] +
    plots_list[["shootA_rootA"]] +
    plot_layout(nrow = 1, guides = "collect") &
    theme(legend.position = "bottom")
  
  ggsave(
    file.path("plots_selected", "scatter_all_indices.png"),
    plot   = p_combined,
    width  = 18,
    height = 6,
    dpi    = 300
  )
  
  ggsave(
    file.path("plots_selected", "scatter_all_indices.pdf"),
    plot   = p_combined,
    width  = 18,
    height = 6,
    device = cairo_pdf
  )
  
  message("Combined figure: plots_selected/scatter_all_indices.(png/pdf)")
}

cat("\nDone! Graphs in 'plots_selected/' and Tables (w/ Tukey, if p<0.05) in 'tables_mean_se/'.\n")

# ==========================================================
# Two-Way ANOVA (regime * inoculant) day-by-day
# for: shoot_side_area, shoot_height, shootA_rootA
# Output: only console printing
# ==========================================================

# Packages
  library(dplyr)
  library(readr)
  library(lubridate)

# File load
file_path_2way <- "PS141_index_2way.txt"

# Import
df2 <- read_tsv(
  file_path_2way,
  show_col_types = FALSE,
  locale = locale(decimal_mark = ".", date_names = "it")
) %>%
  mutate(
    date      = dmy(date),
    n_day     = as.integer(n_day),
    regime    = droplevels(as.factor(regime)),
    inoculant = droplevels(as.factor(inoculant)),
    treatment = droplevels(as.factor(treatment))
  )

# Indexes
metrics_available <- intersect(
  c("shoot_side_area","shoot_height","shootA_rootA"),
  names(df2)
)
if (length(metrics_available) == 0) {
  stop("None of the index is present.")
}

# print ANOVA per day x index
print_aov_2way <- function(day_df, metric, day_label, nice_name) {
  dat <- day_df %>%
    select(regime, inoculant, value = all_of(metric)) %>%
    filter(is.finite(value)) %>%
    droplevels()
  
    ok <- n_distinct(dat$regime) >= 2 &&
    n_distinct(dat$inoculant) >= 2 &&
    nrow(dat) >= 3 &&
    is.finite(var(dat$value)) && var(dat$value) > 0
  
  cat("\n------------------------------------------------------------\n")
  cat(sprintf("[Day %s] 2-way ANOVA (regime * inoculant) — %s\n",
              day_label, nice_name))
  
  if (!ok) {
    cat("No enough data (levels/replicates) or null variance.\n")
    return(invisible(NULL))
  }
  
  fit <- tryCatch(aov(value ~ regime * inoculant, data = dat),
                  error = function(e) NULL)
  if (is.null(fit)) {
    cat("Error in model fit.\n")
    return(invisible(NULL))
  }
  print(summary(fit))
}

# Repeat per day x index
days <- sort(unique(df2$n_day))
nice_names <- c(
  shoot_side_area = "Shoot Side Area",
  shoot_height    = "Shoot Height",
  shootA_rootA    = "Shoot Area / Root Area"
)

for (met in metrics_available) {
  cat("\n============================================================\n")
  cat(sprintf("Index: %s\n", nice_names[[met]]))
  cat("============================================================\n")
  for (d in days) {
    dd <- df2 %>% filter(n_day == d)
    print_aov_2way(dd, metric = met, day_label = d,
                   nice_name = nice_names[[met]])
  }
}

cat("\nDone: printed ANOVA tables in console for each index and day.\n")

# ==========================================================
# One-Way ANOVA within regime (NW e WW), day-by-day
# for: shoot_side_area, shoot_height, shootA_rootA
# Output: only console printing
# ==========================================================

# Packages
  library(dplyr)
  library(readr)
  library(lubridate)
  library(agricolae)

# File load
file_path_1way <- "PS141_index_2way.txt"

# Import
df1 <- read_tsv(
  file_path_1way,
  show_col_types = FALSE,
  locale = locale(decimal_mark = ".", date_names = "it")
) %>%
  mutate(
    date      = dmy(date),
    n_day     = as.integer(n_day),
    regime    = droplevels(as.factor(regime)),
    inoculant = droplevels(as.factor(inoculant))
  )

# Index
metrics_1way <- intersect(
  c("shoot_side_area","shoot_height","shootA_rootA"),
  names(df1)
)
if (length(metrics_1way) == 0) {
  stop("None of the index is present.")
}

nice_names_1way <- c(
  shoot_side_area = "Shoot Side Area",
  shoot_height    = "Shoot Height",
  shootA_rootA    = "Shoot Area / Root Area"
)

# Print One-Way ANOVA (within regime) w/ Tukey (if p<0.05)
run_oneway_in_regime <- function(day_df, metric, day_label, regime_label, pretty_metric) {
  dat <- day_df %>%
    filter(regime == regime_label) %>%
    select(inoculant, value = all_of(metric)) %>%
    filter(is.finite(value)) %>%
    droplevels()
  
  cat("\n------------------------------------------------------------\n")
  cat(sprintf("[Day %s | Regime %s] 1-way ANOVA — %s (factor: inoculant)\n",
              day_label, regime_label, pretty_metric))
  
  ok <- nrow(dat) >= 3 &&
    dplyr::n_distinct(dat$inoculant) >= 2 &&
    (length(unique(dat$value)) > 1)
  
  if (!ok) {
    cat("No enough data (levels/replicates) or null variance.\n")
    return(invisible(NULL))
  }
  
  fit <- tryCatch(aov(value ~ inoculant, data = dat), error = function(e) NULL)
  if (is.null(fit)) {
    cat("Error in model fit.\n")
    return(invisible(NULL))
  }
  
  print(summary(fit))
  
  p_val <- tryCatch(as.numeric(summary(fit)[[1]][["Pr(>F)"]][1]),
                    error = function(e) NA_real_)
  
  # Tukey HSD if significant
  if (is.finite(p_val) && p_val < 0.05) {
    hs <- tryCatch(agricolae::HSD.test(fit, "inoculant", group = TRUE, console = FALSE),
                   error = function(e) NULL)
    if (!is.null(hs) && !is.null(hs$groups)) {

        means_tbl <- dat %>%
        group_by(inoculant) %>%
        summarise(mean = mean(value, na.rm = TRUE), .groups = "drop")
      
      groups_tbl <- hs$groups
      groups_tbl$inoculant <- rownames(groups_tbl)
      
      out_letters <- means_tbl %>%
        left_join(groups_tbl %>% select(inoculant, groups), by = "inoculant") %>%
        arrange(desc(mean))
      
      cat("\nTukey HSD (inoculant) — significance letters:\n")
      print(as.data.frame(out_letters))
    } else {
      cat("\nTukey HSD not available.\n")
    }
  } else {
    cat("\nTukey HSD not ran (ANOVA ns).\n")
  }
}

# Print for each index, regime (NW, WW), for each day
days_1way   <- sort(unique(df1$n_day))
regimes_avl <- intersect(c("NW","WW"), levels(df1$regime))

for (met in metrics_1way) {
  cat("\n============================================================\n")
  cat(sprintf("Index: %s\n", nice_names_1way[[met]]))
  cat("============================================================\n")
  for (reg in regimes_avl) {
    for (d in days_1way) {
      dd <- df1 %>% filter(n_day == d)
      run_oneway_in_regime(dd, metric = met, day_label = d,
                           regime_label = reg,
                           pretty_metric = nice_names_1way[[met]])
    }
  }
}

cat("\nDone: One-Way ANOVA per regime (NW/WW) printed in console for each day and index.\n")
