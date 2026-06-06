# ==========================================================
# Scatterplot (Mean ± SE) + p-value symbols for:
# shoot_side_area, shoot_height
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

# File load
file_path <- "PS140_index.txt"
file_path_2way <- "PS140_index_2way.txt"

# Fixed treatment order and color palette
treat_levels <- c("NW_MIX1","NW_MIXR","NW_NOI","WW_MIX1","WW_MIXR","WW_NOI")
pal_trt <- c(
  "WW_MIXR" = "#117733",
  "WW_MIX1" = "#44AA99",
  "WW_NOI"  = "#332288",
  "NW_MIXR" = "#AA4499",
  "NW_MIX1" = "#CC6677",
  "NW_NOI"  = "#882255"
)

# Index
vars_plot <- c("shoot_side_area", "shoot_height")

# Labelling Y
ylab_map <- c(
  shoot_side_area = "Shoot Side Area [mm²]",
  shoot_height    = "Shoot Height [mm]"
)

# ==========================================================
# 1) PLOTS + CSV (PS140_index.txt)
# ==========================================================

# Import
df <- read_tsv(
  file_path,
  show_col_types = FALSE,
  locale = locale(decimal_mark = ".", date_names = "it")
) %>%
  mutate(
    date      = dmy(date),
    n_day     = as.integer(n_day),
    treatment = factor(as.character(treatment), levels = treat_levels)
  ) %>%
  filter(!is.na(treatment))

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
      dat <- cur_data()
      ok <- n_distinct(dat$treatment) >= 2 && nrow(dat) >= 3
      if (!ok) NA_real_ else tryCatch({
        fit <- stats::aov(value ~ treatment, data = dat)
        as.numeric(summary(fit)[[1]][["Pr(>F)"]][1])
      }, error = function(e) NA_real_)
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
      dd <- droplevels(.x)
      if (n_distinct(dd$treatment) < 2 || nrow(dd) < 3) {
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
        agricolae::HSD.test(aov(value ~ treatment, data = dd),
                            "treatment", group = TRUE, console = FALSE)
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
    scale_color_manual(values = pal_trt, drop = FALSE) +
    labs(x = "Day", y = ylab_txt, color = "Treatment") +
    theme_bw(base_size = 12) +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_blank())
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
  
  day_cols <- paste0("Day ", days_all)
  df_wide <- df_wide %>% select(any_of(c("Treatments", day_cols)))
  
  out_csv <- file.path("tables_mean_se", paste0("mean_se_", metric_name, ".csv"))
  readr::write_csv(df_wide, out_csv)
  message("CSV created: ", out_csv)
}

# Repeat on the indexes
plots_list <- list()

for (met in vars_plot) {
  df_met_plot  <- sum_dat %>% filter(metric == met)
  ann_met_plot <- ann     %>% filter(metric == met)
  df_met_csv   <- sum_with_letters %>% filter(metric == met)
  
  p <- plot_fun(df_met_plot, ann_met_plot, ylab_txt = unname(ylab_map[[met]]))
  plots_list[[met]] <- p
  print(p)
  
  ggsave(file.path("plots_selected", paste0("scatter_", met, ".png")),
         plot = p, width = 7, height = 5, dpi = 300)
  ggsave(file.path("plots_selected", paste0("scatter_", met, ".pdf")),
         plot = p, width = 7, height = 5, device = cairo_pdf)
  
  write_mean_se_csv(df_met_csv, metric_name = met)
}

# Combined figure
if (all(c("shoot_height", "shoot_side_area") %in% names(plots_list))) {
  
  p_combined <- plots_list[["shoot_height"]] +
    plots_list[["shoot_side_area"]] +
    plot_layout(nrow = 1, guides = "collect") &
    theme(legend.position = "bottom")
  
  ggsave(file.path("plots_selected", "scatter_all_indices.png"),
         plot = p_combined, width = 12, height = 6, dpi = 300)
  
  ggsave(file.path("plots_selected", "scatter_all_indices.pdf"),
         plot = p_combined, width = 12, height = 6, device = cairo_pdf)
  
  message("Saved combined figure in plots_selected/scatter_all_indices.(png/pdf)")
}

cat("\nDone! Graphs in 'plots_selected/' and tables in 'tables_mean_se/'.\n")


# ==========================================================
# Two-Way ANOVA (regime * inoculant) day-by-day (PS140_index_2way.txt)
# ==========================================================

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
    treatment = factor(as.character(treatment), levels = treat_levels)
  ) %>%
  filter(!is.na(treatment))

metrics_available <- intersect(vars_plot, names(df2))
if (length(metrics_available) == 0) stop("None of the index is present.")

print_aov_2way <- function(day_df, metric, day_label) {
  dat <- day_df %>%
    select(regime, inoculant, value = all_of(metric)) %>%
    filter(is.finite(value)) %>%
    droplevels()
  
  ok <- n_distinct(dat$regime) >= 2 &&
    n_distinct(dat$inoculant) >= 2 &&
    nrow(dat) >= 3 &&
    is.finite(var(dat$value)) && var(dat$value) > 0
  
  cat("\n------------------------------------------------------------\n")
  cat(sprintf("[Day %s] 2-way ANOVA (regime * inoculant) — %s\n", day_label, metric))
  
  if (!ok) {
    cat("No enough data (levels/replicates) or null variance.\n")
    return(invisible(NULL))
  }
  
  fit <- tryCatch(aov(value ~ regime * inoculant, data = dat), error = function(e) NULL)
  if (is.null(fit)) {
    cat("Error in fit model.\n")
    return(invisible(NULL))
  }
  print(summary(fit))
}

days <- sort(unique(df2$n_day))

cat("\n==================== 2-WAY ANOVA (regime * inoculant) ====================\n")
for (met in metrics_available) {
  cat("\n============================================================\n")
  cat(sprintf("Index: %s\n", met))
  cat("============================================================\n")
  for (d in days) {
    dd <- df2 %>% filter(n_day == d)
    print_aov_2way(dd, metric = met, day_label = d)
  }
}
cat("\nDone: Two-Way ANOVA tables printed in console.\n")


# ==========================================================
# One-Way ANOVA within regime (NW e WW), day-by-day w/ Tukey
# ==========================================================

run_oneway_in_regime <- function(day_df, metric, day_label, regime_label) {
  dat <- day_df %>%
    filter(regime == regime_label) %>%
    select(inoculant, value = all_of(metric)) %>%
    filter(is.finite(value)) %>%
    droplevels()
  
  cat("\n------------------------------------------------------------\n")
  cat(sprintf("[Day %s | Regime %s] 1-way ANOVA — %s (factor: inoculant)\n",
              day_label, regime_label, metric))
  
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
  
  if (is.finite(p_val) && p_val < 0.05) {
    hs <- tryCatch(agricolae::HSD.test(fit, "inoculant", group = TRUE, console = FALSE),
                   error = function(e) NULL)
    if (!is.null(hs) && !is.null(hs$groups)) {
      out <- hs$groups
      out$inoculant <- rownames(out)
      cat("\nTukey HSD (inoculant) — significance letters:\n")
      print(out %>% select(inoculant, groups))
    } else {
      cat("\nTukey HSD not available.\n")
    }
  } else {
    cat("\nTukey HSD: not ran (ANOVA ns).\n")
  }
}

days_1way   <- sort(unique(df2$n_day))
regimes_avl <- intersect(c("NW","WW"), levels(df2$regime))

cat("\n==================== 1-WAY ANOVA WITHIN REGIME (factor: inoculant) ====================\n")
for (met in metrics_available) {
  cat("\n============================================================\n")
  cat(sprintf("Index: %s\n", met))
  cat("============================================================\n")
  for (reg in regimes_avl) {
    for (d in days_1way) {
      dd <- df2 %>% filter(n_day == d)
      run_oneway_in_regime(dd, metric = met, day_label = d, regime_label = reg)
    }
  }
}

cat("\nDone: One-Way ANOVA per regime (NW/WW) printed in console for each index and day.\n")
