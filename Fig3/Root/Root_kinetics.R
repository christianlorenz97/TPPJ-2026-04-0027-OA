# ==========================================================
# Scatterplot (Mean ± SE) + p-value symbols for:
# tot_length, tot_area, Area_Efficiency
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

# File load
file_path <- "PS141_index.txt"

# Helper
safe_div <- function(a, b) ifelse(is.finite(a) & is.finite(b) & b > 0, a / b, NA_real_)

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
    treatment = as.factor(treatment),
    sample    = as.factor(sample)
  )

# Index selection
vars_plot <- c("tot_length", "tot_area", "Area_Efficiency")
ylab_map <- c(
  tot_length       = "Root Total Length [mm]",
  tot_area         = "Root Total Area [mm²]",
  Area_Efficiency  = "Root Area Efficiency [mm⁻¹]"
)

# calculate Mean ± SE w/ long
long_dat <- df %>%
  select(n_day, treatment, sample, all_of(vars_plot)) %>%
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
      if (dplyr::n_distinct(dd$treatment) < 2 || sum(is.finite(dd$value)) < 3) {
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
      tibble(
        treatment = factor(out$treatment, levels = levels(dd$treatment)),
        letter    = as.character(out$groups)
      )
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
    y_top  = ifelse(is.finite(y_top), y_top, max(y_min, na.rm = TRUE)),
    y_min  = ifelse(is.finite(y_min), y_min, y_top),
    y_range = ifelse(is.finite(y_top - y_min) & (y_top - y_min) > 0,
                     y_top - y_min, 0.05 * pmax(y_top, 1, na.rm = TRUE)),
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

# Color mapping for treatments
treatment_colors <- c(
  "WW_MIXR" = "#117733",
  "WW_MIX1" = "#44AA99",
  "WW_NOI"  = "#332288",
  "NW_MIXR" = "#AA4499",
  "NW_MIX1" = "#CC6677",
  "NW_NOI"  = "#882255"
)

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
    scale_color_manual(values = treatment_colors, drop = FALSE) +
    labs(x = "Day", y = ylab_txt, color = "Treatment") +
    theme_bw(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_blank()
    )
}

# CSV "Mean ± SE" w/ Tukey letters
write_mean_se_csv <- function(df_met_letters, metric_name) {
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
  
  day_cols <- paste0("Day ", sort(unique(df_met_letters$n_day)))
  df_wide <- df_wide %>%
    select(any_of(c("Treatments", day_cols)))
  
  out_csv <- file.path("tables_mean_se", paste0("mean_se_", metric_name, ".csv"))
  readr::write_csv(df_wide, out_csv)
  message("CSV created: ", out_csv)
}

# Repeat on the indexes
for (met in vars_plot) {
  df_met_plot  <- sum_dat %>% filter(metric == met)
  ann_met_plot <- ann     %>% filter(metric == met)
  df_met_csv   <- sum_with_letters %>% filter(metric == met)
  
  ylab <- unname(ylab_map[[met]])
  
  p <- plot_fun(df_met_plot, ann_met_plot, ylab_txt = ylab)
  print(p)
  
  ggsave(file.path("plots_selected", paste0("scatter_", met, ".png")),
         plot = p, width = 7, height = 5, dpi = 300)
  ggsave(file.path("plots_selected", paste0("scatter_", met, ".pdf")),
         plot = p, width = 7, height = 5, device = cairo_pdf)
  
  write_mean_se_csv(df_met_csv, metric_name = met)
}

# Combined 3-panel figure (row)
suppressPackageStartupMessages({
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    install.packages("patchwork")
  }
  library(patchwork)
})

# Re-create the three individual plots
p_tot_length <- plot_fun(
  df_met   = sum_dat %>% filter(metric == "tot_length"),
  ann_met  = ann      %>% filter(metric == "tot_length"),
  ylab_txt = ylab_map[["tot_length"]]
)

p_tot_area <- plot_fun(
  df_met   = sum_dat %>% filter(metric == "tot_area"),
  ann_met  = ann      %>% filter(metric == "tot_area"),
  ylab_txt = ylab_map[["tot_area"]]
)

p_area_eff <- plot_fun(
  df_met   = sum_dat %>% filter(metric == "Area_Efficiency"),
  ann_met  = ann      %>% filter(metric == "Area_Efficiency"),
  ylab_txt = ylab_map[["Area_Efficiency"]]
)

# Combine them in one row with a shared legend
p_row <- p_tot_length + p_tot_area + p_area_eff +
  plot_layout(nrow = 1, guides = "collect")

# Print
print(p_row)

# Save combined figure
ggsave(
  filename = file.path("plots_selected", "scatter_all_metrics.png"),
  plot     = p_row,
  width    = 18, height = 5, dpi = 300
)

ggsave(
  filename = file.path("plots_selected", "scatter_all_metrics.pdf"),
  plot     = p_row,
  width    = 18, height = 5, device = cairo_pdf
)

cat("\nDone! Graphs in 'plots_selected/' and Tables (w/ Tukey if p<0.05) in 'tables_mean_se/'.\n")

# ==========================================================
# Two-Way ANOVA (regime * inoculant) day-by-day
# for: tot_length, tot_area, area_efficiency
# Output: only console printing
# ==========================================================

# Packages
  library(dplyr)
  library(readr)
  library(lubridate)

# File load
file_path_2way <- "PS141_index_2way.txt"

# Helper
safe_div <- function(a, b) ifelse(is.finite(a) & is.finite(b) & b > 0, a / b, NA_real_)

# Import
df2 <- read_tsv(
  file_path_2way,
  show_col_types = FALSE,
  locale = locale(decimal_mark = ".", date_names = "it")
) %>%
  mutate(
    date          = dmy(date),
    n_day         = as.integer(n_day),
    regime        = as.factor(regime),
    inoculant     = as.factor(inoculant),
    treatment     = as.factor(treatment),
    sample        = as.factor(sample),
    area_efficiency = safe_div(tot_length, tot_area)
  )

# Indeces
vars2 <- c("tot_length", "tot_area", "area_efficiency")
nice  <- c(tot_length = "Total Length",
           tot_area   = "Total Area",
           area_efficiency = "Area Efficiency")

# print ANOVA per day x index
print_aov_2way <- function(dd, metric, day_label) {
  dat <- dd %>%
    select(regime, inoculant, !!sym(metric)) %>%
    rename(value = !!sym(metric)) %>%
    filter(is.finite(value)) %>%
    droplevels()
  
  ok <- n_distinct(dat$regime) >= 2 &&
    n_distinct(dat$inoculant) >= 2 &&
    nrow(dat) >= 3
  
  cat("\n------------------------------------------------------------\n")
  cat(sprintf("[Day %s] 2-way ANOVA (regime * inoculant) — %s\n",
              day_label, nice[[metric]]))
  if (!ok) {
    cat("No enough data to run 2-way ANOVA (levels/replicates).\n")
    return(invisible(NULL))
  }
  
  fit <- tryCatch(aov(value ~ regime * inoculant, data = dat),
                  error = function(e) NULL)
  if (is.null(fit)) {
    cat("Error in model fitting.\n")
    return(invisible(NULL))
  }
  
  print(summary(fit))
}

# Repeat per day x index
days <- sort(unique(df2$n_day))
for (met in vars2) {
  cat("\n============================================================\n")
  cat(sprintf("Index: %s\n", nice[[met]]))
  cat("============================================================\n")
  for (d in days) {
    dd <- df2 %>% filter(n_day == d)
    print_aov_2way(dd, metric = met, day_label = d)
  }
}

cat("\nDone: printed ANOVA tables in console.\n")

# ==========================================================
# One-Way ANOVA within regime (factor: inoculant), day-by-day
# for: tot_length, tot_area, area_efficiency
# Output: only console printing
# ==========================================================

# Packages
  library(dplyr)
  library(readr)
  library(lubridate)
  library(agricolae)

# File load
  file_path_2way_141 <- "PS141_index_2way.txt"

# Import
  df2 <- read_tsv(
    file_path_2way_141,
    show_col_types = FALSE,
    locale = locale(decimal_mark = ".", date_names = "it")
  )

# Preparation
  df2 <- df2 %>%
  mutate(
    n_day     = as.integer(n_day),
    regime    = droplevels(as.factor(regime)),
    inoculant = droplevels(as.factor(inoculant))
  )

# Index
metrics <- intersect(c("tot_length","tot_area","Area_Efficiency"), names(df2))
if (length(metrics) == 0) stop("None of the index is present.")

# Print One-Way ANOVA (within regime) w/ Tukey (if p<0.05)
print_aov_1way <- function(day_df, metric, day_label, reg_label) {
  dat <- day_df %>%
    select(inoculant, value = all_of(metric)) %>%
    filter(is.finite(value)) %>%
    droplevels()
  
  cat("\n------------------------------------------------------------\n")
  cat(sprintf("[Day %s | Regime: %s] 1-way ANOVA — %s (factor: inoculant)\n",
              day_label, reg_label, metric))
  
  ok <- n_distinct(dat$inoculant) >= 2 && nrow(dat) >= 3 &&
    is.finite(var(dat$value)) && var(dat$value) > 0
  if (!ok) {
    cat("No enough data (levels/replicates) or null variance.\n")
    return(invisible(NULL))
  }
  
  fit <- tryCatch(aov(value ~ inoculant, data = dat), error = function(e) NULL)
  if (is.null(fit)) {
    cat("Error in model fitting.\n")
    return(invisible(NULL))
  }
  
  sm <- summary(fit)
  print(sm)
  
  pval <- tryCatch(as.numeric(sm[[1]][["Pr(>F)"]][1]), error = function(e) NA_real_)
  if (is.finite(pval) && pval < 0.05) {
    cat("\nTukey (HSD.test) — significance letters:\n")
    hs <- tryCatch(
      agricolae::HSD.test(fit, "inoculant", group = TRUE, console = FALSE),
      error = function(e) NULL
    )
    if (!is.null(hs) && !is.null(hs$groups)) {
      out <- hs$groups
      out$inoculant <- rownames(out)
      out <- out[, c("inoculant","groups")]
      colnames(out) <- c("inoculant","letter")
      rownames(out) <- NULL
      print(out[order(out$inoculant), ], row.names = FALSE)
    } else {
      cat("Impossible to calculate Tukey letters.\n")
    }
  } else {
    cat(sprintf("\nTukey (HSD.test): ns (p = %s). No letter.\n",
                ifelse(is.finite(pval), signif(pval, 3), "NA")))
  }
}

# Print index × regime (NW, WW) × day
days_avail    <- sort(unique(df2$n_day))
regimes_avail <- intersect(c("NW","WW"), levels(df2$regime))

for (met in metrics) {
  cat("\n============================================================\n")
  cat(sprintf("Index: %s\n", met))
  cat("============================================================\n")
  for (reg in regimes_avail) {
    cat(sprintf("\n######## Regime: %s ########\n", reg))
    for (d in days_avail) {
      dd <- df2 %>% filter(n_day == d, regime == reg)
      print_aov_1way(dd, metric = met, day_label = d, reg_label = reg)
    }
  }
}

cat("\nDone: One-Way ANOVA (inoculant within regime) printed for each day and index.\n")
