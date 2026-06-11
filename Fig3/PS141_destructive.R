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
file_destr <- "PS141_destr.txt"

# Import
dat <- read_tsv(file_destr, show_col_types = FALSE)

# Treatment order fixed
treat_order <- c("NW_MIX1","NW_MIXR","NW_NOI","WW_MIX1","WW_MIXR","WW_NOI")
present_levels <- treat_order[treat_order %in% unique(dat$Treatment)]

dat <- dat %>%
  mutate(Treatment = factor(Treatment, levels = present_levels))

# Color treatment
trt_cols <- c(
  "WW_MIXR" = "#117733",
  "WW_MIX1" = "#44AA99",
  "WW_NOI"  = "#332288",
  "NW_MIXR" = "#AA4499",
  "NW_MIX1" = "#CC6677",
  "NW_NOI"  = "#882255"
)

# Labelling Y
vars_plot <- c("root_FW","shoot_FW","root_DW","shoot_DW",
               "harvested_root_length","harvested_root_area")

ylabs <- c(
  root_FW = "Root fresh weight [g]",
  shoot_FW = "Shoot fresh weight [g]",
  root_DW = "Root dry weight [g]",
  shoot_DW = "Shoot dry weight [g]",
  harvested_root_length = "Harvested total root length [mm]",
  harvested_root_area   = "Harvested total root area [mm²]"
)

ylims <- list(
  root_FW = c(0, 1),
  shoot_FW = c(0, 1),
  root_DW = c(0, 0.1),
  shoot_DW = c(0, 0.1),
  harvested_root_length = c(0, 2000),
  harvested_root_area   = c(0, 2000)
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

# BARPLOT Mean ± SE + Tukey + ANOVA
make_barplot_one <- function(df, var, ylab_txt, ylim_vec) {
  
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
               by = "Treatment")
  
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
    if (!is.null(hs$groups)) {
      tg <- hs$groups
      tg$Treatment <- rownames(tg)
      letters_df <- tg %>%
        select(Treatment, letters = groups) %>%
        mutate(Treatment = factor(Treatment, levels = trt_levels))
    }
  }
  
  ann_tbl <- sum_tbl %>% left_join(letters_df, by = "Treatment")
  
  # Offset labels
  y_rng <- diff(ylim_vec)
  off_letters <- 0.05 * y_rng
  
  ann_tbl <- ann_tbl %>%
    mutate(y_text = mean + se + off_letters)
  
  ggplot(sum_tbl, aes(Treatment, mean, fill = Treatment)) +
    geom_col(width = 0.7, color = "black", linewidth = 0.35) +
    geom_errorbar(
      aes(ymin = mean - se, ymax = mean + se),
      width = 0.2, linewidth = 0.35
    ) +
    geom_text(
      data = ann_tbl,
      aes(y = y_text, label = letters),
      vjust = 0, size = 4
    ) +
    annotate(
      "text", x = Inf, y = Inf, label = p_sym,
      hjust = 1.05, vjust = 1.3, size = 5
    ) +
    scale_fill_manual(values = trt_cols, limits = trt_levels, drop = FALSE) +
    scale_x_discrete(drop = TRUE) +
    coord_cartesian(ylim = ylim_vec, clip = "off") +
    labs(x = "Treatment", y = ylab_txt) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      plot.margin = unit(c(8, 22, 8, 8), "pt")
    )
}

# Graph plotting
plots <- lapply(vars_plot, function(v)
  make_barplot_one(dat, v, ylabs[[v]], ylims[[v]])
)
names(plots) <- vars_plot

combined <- (plots$root_FW + plots$shoot_FW) /
  (plots$root_DW + plots$shoot_DW) /
  (plots$harvested_root_length + plots$harvested_root_area) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

dir.create("plots_destr", showWarnings = FALSE)
ggsave("plots_destr/PS141_barplots_grid.tiff",
       plot = combined, width = 12, height = 14,
       dpi = 300, compression = "lzw")



####################################################################################################
# Two-Way ANOVA (regime * inoculant) for each index, print in console
####################################################################################################

# Packages
library(dplyr)
library(readr)

# File load
file_destr_2way <- "PS141_destr_2way.txt"

# Import
df2 <- read_tsv(file_destr_2way, show_col_types = FALSE) %>%
  mutate(regime = as.factor(regime),
         inoculant = as.factor(inoculant))

vars_two_way <- c("root_FW","shoot_FW","root_DW","shoot_DW",
                  "harvested_root_length","harvested_root_area")

for (v in vars_two_way) {
  cat("\n============================================================\n")
  cat(sprintf("2-way ANOVA (regime * inoculant) — %s\n", v))
  cat("============================================================\n")
  dd <- df2 %>% filter(is.finite(.data[[v]]))
  ok <- dplyr::n_distinct(dd$regime) >= 2 && dplyr::n_distinct(dd$inoculant) >= 2 && nrow(dd) >= 3
  if (!ok) { cat("Dati insufficienti per la 2-way ANOVA.\n"); next }
  fit <- tryCatch(aov(reformulate(c("regime","inoculant","regime:inoculant"), response = v), data = dd),
                  error = function(e) NULL)
  if (is.null(fit)) { cat("Errore nel fit del modello.\n"); next }
  print(summary(fit))
}

####################################################################################################
# One-Way ANOVA within regime (NW e WW) + Tukey letters
####################################################################################################

# Packages
  library(dplyr)
  library(readr)
  library(agricolae)

# File load
  file_destr_2way <- "PS141_destr_2way.txt"

# Import
  df2 <- read_tsv(file_destr_2way, show_col_types = FALSE) %>%
    mutate(regime = as.factor(regime), inoculant = as.factor(inoculant))

vars_one_way <- c("root_FW","shoot_FW","root_DW","shoot_DW",
                  "harvested_root_length","harvested_root_area")
regimes <- intersect(c("WW","NW"), levels(df2$regime))

for (v in vars_one_way) {
  cat("\n############################################################\n")
  cat(sprintf("1-way ANOVA within regime — %s (factor: inoculant)\n", v))
  cat("############################################################\n")
  for (r in regimes) {
    cat(sprintf("\n-- Regime: %s --\n", r))
    dd <- df2 %>% filter(regime == r, is.finite(.data[[v]])) %>% droplevels()
    ok <- dplyr::n_distinct(dd$inoculant) >= 2 && nrow(dd) >= 3 && var(dd[[v]]) > 0
    if (!ok) { cat("No enough data or null variance.\n"); next }
    fit <- tryCatch(aov(reformulate("inoculant", response = v), data = dd), error = function(e) NULL)
    if (is.null(fit)) { cat("Error in model fit.\n"); next }
    print(summary(fit))
    # Tukey letters
    hs <- tryCatch(agricolae::HSD.test(fit, "inoculant", group = TRUE, console = FALSE), error = function(e) NULL)
    if (!is.null(hs) && !is.null(hs$groups)) {
      grp <- hs$groups; grp$inoculant <- rownames(grp)
      means <- dd %>% group_by(inoculant) %>% summarise(mean = mean(.data[[v]], na.rm = TRUE), .groups="drop")
      out <- means %>% left_join(grp %>% select(inoculant, groups), by = "inoculant") %>%
        arrange(desc(mean)) %>% rename(letter = groups)
      cat("\nTukey HSD — significance letters:\n")
      print(as.data.frame(out), row.names = FALSE)
    } else {
      cat("Tukey not available.\n")
    }
  }
}
