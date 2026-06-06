####################################################################################################
# Barplot Mean±SE + One-Way ANOVA w/ Tukey
####################################################################################################

# Packages
  library(dplyr)
  library(readr)
  library(tidyr)
  library(ggplot2)
  library(agricolae)
  library(patchwork)

# File load
file_destr <- "PS140_destr.txt"

# Import
dat <- read_tsv(file_destr, show_col_types = FALSE)

# Check column
vars_plot <- c(
  "root_FW","shoot_FW","root_DW","shoot_DW",
  "harvested_root_length","harvested_root_area",
  "harvested_root_area_efficiency","harvested_root_to_shoot_area"
)

required_cols <- c("Treatment", vars_plot)
missing_cols <- setdiff(required_cols, names(dat))
if (length(missing_cols) > 0) {
  stop("Lacking column in input file: ", paste(missing_cols, collapse = ", "))
}


# Fixed treatment order
treat_order <- c("NW_MIX1","NW_MIXR","NW_NOI","WW_MIX1","WW_MIXR","WW_NOI")
present_levels <- treat_order[treat_order %in% unique(dat$Treatment)]

dat <- dat %>%
  mutate(Treatment = factor(Treatment, levels = present_levels))

# Treatment color
trt_cols <- c(
  "WW_MIXR" = "#117733",
  "WW_MIX1" = "#44AA99",
  "WW_NOI"  = "#332288",
  "NW_MIXR" = "#AA4499",
  "NW_MIX1" = "#CC6677",
  "NW_NOI"  = "#882255"
)

# Labelling Y
ylabs <- c(
  root_FW = "Root fresh weight [g]",
  shoot_FW = "Shoot fresh weight [g]",
  root_DW = "Root dry weight [g]",
  shoot_DW = "Shoot dry weight [g]",
  harvested_root_length = "Harvested total root length [mm]",
  harvested_root_area = "Harvested total root area [mm²]",
  harvested_root_area_efficiency = "Harvested root area efficiency [mm⁻¹]",
  harvested_root_to_shoot_area = "Harvested root to shoot area"
)

# Limits Y
ylims_pref <- list(
  root_FW = c(0, 1),
  shoot_FW = c(0, 1),
  root_DW = c(0, 0.1),
  shoot_DW = c(0, 0.1),
  harvested_root_length = c(0, 2000),
  harvested_root_area = c(0, 2000),
  harvested_root_area_efficiency = c(0, 2),
  harvested_root_to_shoot_area  = c(0, 5)
)

# p-value -> symbol
p_to_sym <- function(p) {
  dplyr::case_when(
    !is.finite(p) ~ "ns",
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    TRUE      ~ "ns"
  )
}

# robust Y lim
resolve_ylim <- function(df, var, ylim_pref = NULL, pad_frac = 0.12) {
  x <- df[[var]]
  x <- x[is.finite(x)]
  if (length(x) == 0) return(c(0, 1))
  
  dmax <- max(x, na.rm = TRUE)
  dmin <- min(x, na.rm = TRUE)
  
  # default: below 0 if positive variable, otherwise min - pad
  base_low <- if (dmin >= 0) 0 else dmin
  base_high <- dmax
  
  # padding
  rng <- max(1e-9, base_high - base_low)
  auto <- c(base_low, base_high + pad_frac * rng)
  
  if (is.null(ylim_pref) || any(!is.finite(ylim_pref)) || length(ylim_pref) != 2) {
    return(auto)
  }
  
  # if pref cut data, expand (keeping low pref if needed)
  if (ylim_pref[2] < dmax) {
    new_low <- min(ylim_pref[1], auto[1])
    new_high <- max(ylim_pref[2], auto[2])
    message(sprintf("Note: ylim for '%s' expanded from [%g, %g] to [%g, %g] to include data.",
                    var, ylim_pref[1], ylim_pref[2], new_low, new_high))
    return(c(new_low, new_high))
  }
  
  # if pref not cut, use pref
  return(ylim_pref)
}

# BARPLOT Mean ± SE + Tukey + ANOVA
make_barplot_one <- function(df, var, ylab_txt, ylim_pref = NULL) {
  
  df <- df %>% filter(is.finite(.data[[var]])) %>% droplevels()
  trt_levels <- levels(df$Treatment)
  
  # Mean ± SE
  sum_tbl <- df %>%
    group_by(Treatment) %>%
    summarise(
      n = n(),
      mean = mean(.data[[var]], na.rm = TRUE),
      sd   = sd(.data[[var]], na.rm = TRUE),
      se   = if_else(n > 1, sd / sqrt(n), 0),
      .groups = "drop"
    ) %>%
    right_join(tibble(Treatment = factor(trt_levels, levels = trt_levels)),
               by = "Treatment") %>%
    arrange(Treatment)
  
  # ANOVA
  fit  <- tryCatch(aov(reformulate("Treatment", response = var), data = df),
                   error = function(e) NULL)
  pval <- if (!is.null(fit)) summary(fit)[[1]][["Pr(>F)"]][1] else NA_real_
  p_sym <- p_to_sym(pval)
  
  # Tukey letters
  letters_df <- tibble(Treatment = trt_levels, letters = NA_character_)
  if (!is.null(fit)) {
    hs <- tryCatch(HSD.test(fit, "Treatment", group = TRUE, console = FALSE),
                   error = function(e) NULL)
    if (!is.null(hs) && !is.null(hs$groups)) {
      tg <- hs$groups
      tg$Treatment <- rownames(tg)
      letters_df <- tg %>%
        select(Treatment, letters = groups) %>%
        mutate(Treatment = factor(Treatment, levels = trt_levels)) %>%
        arrange(Treatment)
    }
  }
  
  ann_tbl <- sum_tbl %>% left_join(letters_df, by = "Treatment")
  
  # Ylim robust + labels offset
  ylim_vec <- resolve_ylim(df, var, ylim_pref = ylim_pref, pad_frac = 0.12)
  y_rng <- diff(ylim_vec)
  off_letters <- 0.05 * y_rng
  
  ann_tbl <- ann_tbl %>% mutate(y_text = mean + se + off_letters)
  
  p <- ggplot(sum_tbl, aes(Treatment, mean, fill = Treatment)) +
    geom_col(width = 0.7, color = "black", linewidth = 0.35) +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                  width = 0.2, linewidth = 0.35) +
    geom_text(data = ann_tbl,
              aes(y = y_text, label = letters),
              vjust = 0, size = 4, na.rm = TRUE) +
    annotate("text", x = Inf, y = Inf, label = p_sym,
             hjust = 1.05, vjust = 1.3, size = 5) +
    scale_fill_manual(values = trt_cols, limits = trt_levels, drop = FALSE) +
    coord_cartesian(ylim = ylim_vec, clip = "off") +
    labs(x = "Treatment", y = ylab_txt) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      plot.margin = unit(c(8, 22, 8, 8), "pt")
    )
  
  # Summary export
  out_sum <- ann_tbl %>%
    mutate(
      variable = var,
      ylab = ylab_txt,
      anova_p = pval,
      anova_sym = p_sym
    ) %>%
    select(variable, ylab, Treatment, n, mean, se, letters, anova_p, anova_sym)
  
  list(plot = p, summary = out_sum)
}

# Graph plotting and summary
res <- lapply(vars_plot, function(v) {
  make_barplot_one(dat, v, ylabs[[v]], ylims_pref[[v]])
})
names(res) <- vars_plot

plots <- lapply(res, `[[`, "plot")
sum_all <- bind_rows(lapply(res, `[[`, "summary"))

# Combined figure
combined <- (plots$root_FW + plots$shoot_FW) /
  (plots$root_DW + plots$shoot_DW) /
  (plots$harvested_root_length + plots$harvested_root_area) /
  (plots$harvested_root_area_efficiency + plots$harvested_root_to_shoot_area) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

# Output directory
dir.create("plots_destr", showWarnings = FALSE)

ggsave(
  filename = "plots_destr/PS140_barplots_grid.tiff",
  plot = combined,
  width = 12, height = 14,
  dpi = 300, compression = "lzw"
)

write_csv(
  sum_all,
  file = "plots_destr/PS140_mean_se_tukey_anova.csv"
)

message("OK — saved:")
message(" - plots_destr/PS140_barplots_grid.tiff")
message(" - plots_destr/PS140_mean_se_tukey_anova.csv")



####################################################################################################
# CSV Mean ± SE + Tukey (Mean ± SE)
####################################################################################################
dir.create("plots_destr", showWarnings = FALSE)

fmt_num <- function(x, digits = 3) {
  ifelse(is.finite(x), formatC(x, format = "f", digits = digits), NA_character_)
}

for (v in vars_plot) {
  
  tbl_v <- sum_all %>%
    filter(variable == v) %>%
    mutate(
      Mean_fmt = fmt_num(mean, digits = 3),
      SE_fmt   = fmt_num(se,   digits = 3),
      letters  = if_else(is.na(letters) | trimws(letters) == "", "", as.character(letters)),
      `Mean ± SE` = if_else(
        letters == "",
        paste0(Mean_fmt, " ± ", SE_fmt),
        paste0(Mean_fmt, " ± ", SE_fmt, " ", letters)
      )
    ) %>%
    transmute(
      Treatment = as.character(Treatment),
      `Mean ± SE` = `Mean ± SE`
    )
  
  out_name <- file.path("plots_destr", paste0("PS140_mean_se_Tukey_", v, ".csv"))
  write_csv(tbl_v, out_name)
}

message("OK — CSV w/ 'Mean ± SE' and Tukey.")



####################################################################################################
# Two-Way ANOVA (regime * inoculant) for each index, print in console
####################################################################################################

# Packages
  library(dplyr)
  library(readr)

file_destr_2way <- "PS140_destr_2way.txt"

df2 <- read_tsv(file_destr_2way, show_col_types = FALSE) %>%
  mutate(
    regime = as.factor(regime),
    inoculant = as.factor(inoculant)
  )

# Column check
required_cols2 <- c("regime", "inoculant", vars_plot)
missing_cols2 <- setdiff(required_cols2, names(df2))
if (length(missing_cols2) > 0) {
  stop("Lacking columns in PS140_destr_2way.txt: ", paste(missing_cols2, collapse = ", "))
}

vars_two_way <- c("root_FW","shoot_FW","root_DW","shoot_DW",
                  "harvested_root_length","harvested_root_area",
                  "harvested_root_area_efficiency","harvested_root_to_shoot_area")

for (v in vars_two_way) {
  cat("\n============================================================\n")
  cat(sprintf("2-way ANOVA (regime * inoculant) — %s\n", v))
  cat("============================================================\n")
  
  dd <- df2 %>% filter(is.finite(.data[[v]]))
  
  ok <- dplyr::n_distinct(dd$regime) >= 2 &&
    dplyr::n_distinct(dd$inoculant) >= 2 &&
    nrow(dd) >= 3 &&
    stats::var(dd[[v]], na.rm = TRUE) > 0
  
  if (!ok) { cat("No enough data or null variance.\n"); next }
  
  fit <- tryCatch(
    aov(reformulate(c("regime","inoculant","regime:inoculant"), response = v), data = dd),
    error = function(e) NULL
  )
  if (is.null(fit)) { cat("Error in model fit.\n"); next }
  
  print(summary(fit))
}

####################################################################################################
# One-Way ANOVA within regime (NW e WW) + Tukey letters
####################################################################################################

# Packages
  library(dplyr)
  library(agricolae)

vars_one_way <- c("root_FW","shoot_FW","root_DW","shoot_DW",
                  "harvested_root_length","harvested_root_area",
                  "harvested_root_area_efficiency","harvested_root_to_shoot_area")

regimes <- intersect(c("WW","NW"), levels(df2$regime))

for (v in vars_one_way) {
  cat("\n############################################################\n")
  cat(sprintf("1-way ANOVA within regime — %s (factor: inoculant)\n", v))
  cat("############################################################\n")
  
  for (r in regimes) {
    cat(sprintf("\n-- Regime: %s --\n", r))
    
    dd <- df2 %>%
      filter(regime == r, is.finite(.data[[v]])) %>%
      droplevels()
    
    ok <- dplyr::n_distinct(dd$inoculant) >= 2 &&
      nrow(dd) >= 3 &&
      stats::var(dd[[v]], na.rm = TRUE) > 0
    
    if (!ok) { cat("No enough data or null variance.\n"); next }
    
    fit <- tryCatch(aov(reformulate("inoculant", response = v), data = dd),
                    error = function(e) NULL)
    if (is.null(fit)) { cat("Error in model fit.\n"); next }
    
    print(summary(fit))
    
    hs <- tryCatch(agricolae::HSD.test(fit, "inoculant", group = TRUE, console = FALSE),
                   error = function(e) NULL)
    
    if (!is.null(hs) && !is.null(hs$groups)) {
      grp <- hs$groups
      grp$inoculant <- rownames(grp)
      
      means <- dd %>%
        group_by(inoculant) %>%
        summarise(mean = mean(.data[[v]], na.rm = TRUE), .groups = "drop")
      
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
