################################################################################
# Barplot (Mean ± SE) + p-value symbols for relative growth rate
# + CSV "Mean ± SE" for each index w/ Tukey (if p<0.05)
################################################################################

# Packages
  library(dplyr)
  library(readr)
  library(tidyr)
  library(ggplot2)
  library(agricolae)
  library(patchwork)
  library(purrr)

# File load + import
file_index      <- "PS140_index.txt"
file_index_2way <- "PS140_index_2way.txt"
dat_index  <- read_tsv(file_index,      show_col_types = FALSE)
dat_index2 <- read_tsv(file_index_2way, show_col_types = FALSE)

# Preparation
treat_order <- c("NW_MIX1","NW_MIXR","NW_NOI",
                 "WW_MIX1","WW_MIXR","WW_NOI")
treat_colors <- c(
  "WW_MIXR" = "#117733",
  "WW_MIX1" = "#44AA99",
  "WW_NOI"  = "#332288",
  "NW_MIXR" = "#AA4499",
  "NW_MIX1" = "#CC6677",
  "NW_NOI"  = "#882255"
)

present_levels <- treat_order[treat_order %in% unique(dat_index2$treatment)]

df2 <- dat_index2 %>%
  mutate(
    n_day     = as.integer(n_day),
    treatment = factor(treatment, levels = present_levels),
    regime    = factor(regime),
    inoculant = factor(inoculant),
    sample    = factor(sample)
  )

# Global max for each variable
global_max <- df2 %>%
  summarise(
    shoot_height_max    = max(shoot_height,    na.rm = TRUE),
    shoot_side_area_max = max(shoot_side_area, na.rm = TRUE)
  )

G_sh_h    <- global_max$shoot_height_max
G_sh_area <- global_max$shoot_side_area_max

# Save max
dir.create("tables_index_RESI_meanse", showWarnings = FALSE)
write_csv(global_max, file.path("tables_index_RESI_meanse", "PS140_RESI_global_max.csv"))

# Function (AUC trapezium method / (global_max * T_tot))
compute_resi_one <- function(df_one, time_col = "n_day", value_col, global_max) {
  dfp <- df_one %>%
    filter(is.finite(.data[[time_col]]),
           is.finite(.data[[value_col]])) %>%
    arrange(.data[[time_col]]) %>%
    distinct(.data[[time_col]], .keep_all = TRUE)
  
  if (nrow(dfp) < 2) return(NA_real_)
  
  t <- dfp[[time_col]]
  x <- dfp[[value_col]]
  
  T_tot <- max(t) - min(t)
  if (!is.finite(T_tot) || T_tot <= 0 ||
      !is.finite(global_max) || global_max <= 0) {
    return(NA_real_)
  }
  
  dt    <- diff(t)
  x_mid <- (x[-1] + x[-length(x)]) / 2
  auc   <- sum(x_mid * dt)
  
  resi <- auc / (global_max * T_tot)
  resi <- pmin(pmax(resi, 0), 1)
  
  as.numeric(resi)
}

# Calculation (group_split + map_dfr)
indices_long <- df2 %>%
  group_by(regime, inoculant, treatment, sample) %>%
  group_split() %>%
  purrr::map_dfr(function(df_one) {
    tibble(
      regime    = df_one$regime[1],
      inoculant = df_one$inoculant[1],
      treatment = df_one$treatment[1],
      sample    = df_one$sample[1],
      
      Shoot_Elongation_Rate =
        compute_resi_one(df_one, value_col = "shoot_height", global_max = G_sh_h),
      
      Shoot_Side_Area_Enlargement_Rate =
        compute_resi_one(df_one, value_col = "shoot_side_area", global_max = G_sh_area)
    )
  }) %>%
  ungroup() %>%
  mutate(
    treatment = factor(treatment, levels = present_levels),
    regime    = factor(regime),
    inoculant = factor(inoculant),
    sample    = factor(sample)
  )

# CSV Mean +- SE
write_csv(indices_long,
          file.path("tables_index_RESI_meanse", "PS140_RESI_all_indices_long.csv"))

# Visualization
p_to_sym <- function(p) {
  dplyr::case_when(
    !is.finite(p) ~ "ns",
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    TRUE      ~ "ns"
  )
}

make_barplot_index <- function(df_idx, var, ylab_txt) {
  df <- df_idx %>%
    filter(is.finite(.data[[var]])) %>%
    droplevels()
  
  trt_levels <- levels(df$treatment)
  last_trt   <- tail(trt_levels, 1)
  
  # ANOVA 1-via (treatment)
  fit <- tryCatch(
    aov(reformulate("treatment", response = var), data = df),
    error = function(e) NULL
  )
  pval <- if (!is.null(fit)) {
    tryCatch(as.numeric(summary(fit)[[1]][["Pr(>F)"]][1]), error = function(e) NA_real_)
  } else NA_real_
  
  p_sym <- p_to_sym(pval)
  
  # Tukey letters (if ANOVA ns -> no letters)
  letters_df <- tibble(treatment = trt_levels, letters = "")
  if (!is.null(fit) && p_sym != "ns") {
    hs <- tryCatch(
      agricolae::HSD.test(fit, "treatment", group = TRUE, console = FALSE),
      error = function(e) NULL
    )
    if (!is.null(hs) && !is.null(hs$groups)) {
      tg <- hs$groups
      tg$treatment <- rownames(tg)
      letters_df <- tg %>%
        select(treatment, letters = groups) %>%
        mutate(treatment = factor(treatment, levels = trt_levels))
    }
  }
  
  # Mean ± SE
  sum_tbl <- df %>%
    group_by(treatment) %>%
    summarise(
      n    = sum(is.finite(.data[[var]])),
      mean = mean(.data[[var]], na.rm = TRUE),
      sd   = sd(.data[[var]],   na.rm = TRUE),
      se   = if_else(n > 1, sd / sqrt(n), NA_real_),
      .groups = "drop"
    ) %>%
    mutate(treatment = factor(treatment, levels = trt_levels)) %>%
    left_join(letters_df, by = "treatment") %>%
    mutate(y_lab = pmin(mean + se + 0.03, 0.95))
  
  y_for_symbol <- 0.97
  
  p <- ggplot(sum_tbl, aes(x = treatment, y = mean, fill = treatment)) +
    geom_col(width = 0.65, alpha = 0.95, color = "black") +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                  width = 0.2, linewidth = 0.5) +
    geom_text(aes(y = y_lab, label = letters),
              vjust = 0, size = 4, color = "black") +
    annotate("text",
             x = last_trt, y = y_for_symbol,
             label = p_sym,
             hjust = 1, vjust = 0,
             size = 5) +
    scale_x_discrete(drop = FALSE, limits = trt_levels) +
    scale_fill_manual(values = treat_colors, drop = FALSE, limits = trt_levels) +
    scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.03))) +
    coord_cartesian(ylim = c(0, 1), clip = "off") +
    labs(x = "Treatment", y = ylab_txt, fill = "Treatment") +
    theme_bw(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position  = "bottom",
      plot.margin      = unit(c(10, 20, 5.5, 5.5), "pt")
    )
  
  list(plot = p, letters = letters_df, pval = pval, summary = sum_tbl)
}

# Barplot
ylabs_resi <- c(
  Shoot_Elongation_Rate            = "Shoot Elongation Rate",
  Shoot_Side_Area_Enlargement_Rate = "Shoot Side Area Enlargement Rate"
)
vars_idx <- names(ylabs_resi)

res_plots <- lapply(vars_idx, function(v) {
  make_barplot_index(indices_long, v, ylabs_resi[[v]])
})
names(res_plots) <- vars_idx

combined_resi <- (res_plots[["Shoot_Elongation_Rate"]]$plot /
                    res_plots[["Shoot_Side_Area_Enlargement_Rate"]]$plot) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

dir.create("plots_index_RESI", showWarnings = FALSE)
ggsave("plots_index_RESI/PS140_RESI_barplots_2indices.tiff",
       plot   = combined_resi,
       width  = 10,
       height = 10,
       dpi    = 300,
       compression = "lzw")

# CSV Mean ± SE (w/ letters) 
write_meanse_tukey_RESI <- function(df_idx, var, letters_tbl) {
  trt_levels <- levels(df_idx$treatment)
  
  sum_tbl <- df_idx %>%
    filter(is.finite(.data[[var]])) %>%
    group_by(treatment) %>%
    summarise(
      n    = sum(is.finite(.data[[var]])),
      mean = mean(.data[[var]], na.rm = TRUE),
      sd   = sd(.data[[var]],   na.rm = TRUE),
      se   = if_else(n > 1, sd / sqrt(n), NA_real_),
      .groups = "drop"
    ) %>%
    mutate(
      mean      = round(mean, 3),
      se        = round(se,   3),
      treatment = factor(treatment, levels = trt_levels)
    ) %>%
    right_join(tibble(treatment = factor(trt_levels, levels = trt_levels)), by = "treatment") %>%
    left_join(letters_tbl, by = "treatment") %>%
    mutate(
      letters = if_else(is.na(letters) | trimws(letters) == "", "", as.character(letters)),
      cell = dplyr::case_when(
        is.finite(se) & letters != "" ~ paste0(mean, " ± ", se, " ", letters),
        is.finite(se)                 ~ paste0(mean, " ± ", se),
        TRUE                          ~ as.character(mean)
      )
    ) %>%
    select(treatment, cell)
  
  out <- tibble(Treatments = as.character(sum_tbl$treatment),
                !!var := sum_tbl$cell) %>%
    bind_rows(tibble(Treatments = "Treatments", !!var := var), .)
  
  write_csv(out, file.path("tables_index_RESI_meanse",
                           paste0("mean_se_Tukey_RESI_", var, "_PS140.csv")))
}

invisible(mapply(function(v, lst) write_meanse_tukey_RESI(indices_long, v, lst$letters),
                 vars_idx, res_plots))

# Two-way ANOVA (regime * inoculant)
for (v in vars_idx) {
  cat("\n============================================================\n")
  cat(sprintf("2-way ANOVA (regime * inoculant) — %s\n", v))
  cat("============================================================\n")
  
  dd <- indices_long %>%
    filter(is.finite(.data[[v]])) %>%
    droplevels()
  
  ok <- dplyr::n_distinct(dd$regime)    >= 2 &&
    dplyr::n_distinct(dd$inoculant) >= 2 &&
    nrow(dd) >= 3 &&
    var(dd[[v]], na.rm = TRUE) > 0
  
  if (!ok) { cat("No enough data or null variance.\n"); next }
  
  fit2 <- tryCatch(
    aov(reformulate(c("regime","inoculant","regime:inoculant"), response = v), data = dd),
    error = function(e) NULL
  )
  if (is.null(fit2)) { cat("Error in model fit.\n"); next }
  
  print(summary(fit2))
}

# One-Way ANOVA within regime (factor: inoculant) w/ Tukey
regimes <- intersect(c("WW","NW"), levels(indices_long$regime))

for (v in vars_idx) {
  cat("\n############################################################\n")
  cat(sprintf("1-way ANOVA within regime — %s (factor: inoculant)\n", v))
  cat("############################################################\n")
  
  for (r in regimes) {
    cat(sprintf("\n-- Regime: %s --\n", r))
    
    dd <- indices_long %>%
      filter(regime == r, is.finite(.data[[v]])) %>%
      droplevels()
    
    ok <- dplyr::n_distinct(dd$inoculant) >= 2 &&
      nrow(dd) >= 3 &&
      var(dd[[v]], na.rm = TRUE) > 0
    
    if (!ok) { cat("No enough data or null variance.\n"); next }
    
    fit <- tryCatch(aov(reformulate("inoculant", response = v), data = dd),
                    error = function(e) NULL)
    if (is.null(fit)) { cat("Error in model fit.\n"); next }
    
    print(summary(fit))
    
    hs <- tryCatch(agricolae::HSD.test(fit, "inoculant", group = TRUE, console = FALSE),
                   error = function(e) NULL)
    
    if (!is.null(hs) && !is.null(hs$groups)) {
      grp <- hs$groups; grp$inoculant <- rownames(grp)
      means <- dd %>% group_by(inoculant) %>% summarise(mean = mean(.data[[v]], na.rm = TRUE), .groups = "drop")
      out <- means %>%
        left_join(grp %>% select(inoculant, groups), by = "inoculant") %>%
        arrange(desc(mean)) %>%
        rename(letter = groups)
      cat("\nTukey HSD — significance letters:\n")
      print(as.data.frame(out), row.names = FALSE)
    } else {
      cat("Tukey HSD not available.\n")
    }
  }
}
