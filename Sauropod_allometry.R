# ============================================================
# Title : Extreme Value Theory Reveals the Upper Limit of
#         Terrestrial Body Size: Sauropods
# Author: Bastiaan A. Van Velthoven
# File  : allometry_regressions.R
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readxl)
  library(grid)   # for unit()
  library(glue)   # for nice printing
})

DATA_XLSX <- "Data/sauropod_measurements_demic.xlsx"
df_raw    <- read_excel(DATA_XLSX)

# ---------------------------
# Output directory
# ---------------------------
FIG_DIR <- "Figures/Sauropods"
if (!dir.exists(FIG_DIR)) dir.create(FIG_DIR, recursive = TRUE)

theme_science <- theme_minimal(base_family = "Arial", base_size = 12) +
  theme(
    axis.title        = element_text(size = 14, face = "bold"),
    axis.text         = element_text(size = 12),
    legend.title      = element_text(size = 10, face = "bold"),
    legend.text       = element_text(size = 10),
    panel.grid.major  = element_line(color = "gray80", linewidth = 0.3),
    panel.grid.minor  = element_blank(),
    panel.border      = element_rect(color = "black", fill = NA, linewidth = 0.6),
    axis.ticks.length = unit(0.20, "cm"),
    axis.ticks        = element_line(color = "black", linewidth = 0.4),
    plot.margin       = margin(5, 5, 5, 5),
    legend.position   = "right"
  )

# ---------------------------
# Choose scale: log–log 
# ---------------------------
use_log <- TRUE

# ---------------------------
# Column names
# ---------------------------
col_CF <- "femur circ (mm)"
col_CH <- "humerus circ (mm)"
col_LF <- "femur L"
col_LH <- "humerus L"

# ---------------------------
# Extract + numeric coercion + logs
# ---------------------------
df_allom <- df_raw %>%
  transmute(
    specimen = `genus and species`,
    CF = suppressWarnings(as.numeric(.data[[col_CF]])),
    CH = suppressWarnings(as.numeric(.data[[col_CH]])),
    LF = suppressWarnings(as.numeric(.data[[col_LF]])),
    LH = suppressWarnings(as.numeric(.data[[col_LH]]))
  ) %>%
  mutate(
    log_CF = ifelse(is.finite(CF) & CF > 0, log(CF), NA_real_),
    log_CH = ifelse(is.finite(CH) & CH > 0, log(CH), NA_real_),
    log_LF = ifelse(is.finite(LF) & LF > 0, log(LF), NA_real_),
    log_LH = ifelse(is.finite(LH) & LH > 0, log(LH), NA_real_)
  )

`%||%` <- function(a, b) if (!is.null(a)) a else b

fit_pair_lm <- function(dat, x, y, label = NULL) {
  d <- dat %>% filter(is.finite(.data[[x]]), is.finite(.data[[y]]))
  n <- nrow(d)
  if (n < 3) {
    return(tibble(
      model = label %||% paste(y, "~", x),
      n = n,
      intercept = NA_real_, slope = NA_real_,
      se_intercept = NA_real_, se_slope = NA_real_,
      r2 = NA_real_
    ))
  }
  mod <- lm(reformulate(x, response = y), data = d)
  sm  <- summary(mod)
  co  <- coef(sm)
  
  tibble(
    model        = label %||% paste(y, "~", x),
    n            = n,
    intercept    = unname(co[1, "Estimate"]),
    slope        = unname(co[2, "Estimate"]),
    se_intercept = unname(co[1, "Std. Error"]),
    se_slope     = unname(co[2, "Std. Error"]),
    r2           = unname(sm$r.squared)
  )
}

plot_pair <- function(dat, x, y, xlab, ylab, fname) {
  d <- dat %>% filter(is.finite(.data[[x]]), is.finite(.data[[y]]))
  p <- ggplot(d, aes(x = .data[[x]], y = .data[[y]])) +
    geom_point(alpha = 0.75, size = 2) +
    geom_smooth(method = "lm", se = FALSE) +
    labs(x = xlab, y = ylab) +
    theme_science
  ggsave(file.path(FIG_DIR, fname), p, dpi = 600, w = 6.2, h = 4.8, units = "in")
  p
}

# ---------------------------
# Fit regressions (log–log => power laws)
#   - CF vs CH (as before): log(CF) ~ log(CH)
#   - Length vs Circ (requested direction): log(LF) ~ log(CF), log(LH) ~ log(CH)
# ---------------------------
if (isTRUE(use_log)) {
  
  reg_tbl <- bind_rows(
    fit_pair_lm(df_allom, "log_CH", "log_CF", "log(CF) ~ log(CH)"),
    fit_pair_lm(df_allom, "log_CF", "log_LF", "log(LF) ~ log(CF)"),  # <-- swapped
    fit_pair_lm(df_allom, "log_CH", "log_LH", "log(LH) ~ log(CH)")   # <-- swapped
  )
  
  print(reg_tbl)
  write.csv(reg_tbl, file.path(FIG_DIR, "allometry_regressions_R2.csv"), row.names = FALSE)
  
  # --- Plots (requested axes) ---
  plot_pair(df_allom, "log_CH", "log_CF", "log(CH)", "log(CF)", "reg_logCF_logCH.png")
  plot_pair(df_allom, "log_CF", "log_LF", "log(CF)", "log(LF)", "reg_logLF_logCF.png")
  plot_pair(df_allom, "log_CH", "log_LH", "log(CH)", "log(LH)", "reg_logLH_logCH.png")
  
  # --- Print power-law forms on original scale ---
  # If log(Y) = a0 + b log(X), then Y = exp(a0) * X^b
  row_LF_CF <- reg_tbl %>% filter(model == "log(LF) ~ log(CF)")
  row_LH_CH <- reg_tbl %>% filter(model == "log(LH) ~ log(CH)")
  
  if (nrow(row_LF_CF) == 1 && is.finite(row_LF_CF$intercept) && is.finite(row_LF_CF$slope)) {
    a <- exp(row_LF_CF$intercept)
    b <- row_LF_CF$slope
    cat(glue("\nPower-law :  LF = a * CF^b\n"))
    cat(glue("  a = exp({signif(row_LF_CF$intercept, 6)}) = {signif(a, 6)}\n"))
    cat(glue("  b = {signif(b, 6)}\n"))
    cat(glue("  => LF = {signif(a, 6)} * CF^{signif(b, 6)}\n"))
  }
  
  if (nrow(row_LH_CH) == 1 && is.finite(row_LH_CH$intercept) && is.finite(row_LH_CH$slope)) {
    a2 <- exp(row_LH_CH$intercept)
    b2 <- row_LH_CH$slope
    cat(glue("\nPower-law :  LH = a' * CH^b'\n"))
    cat(glue("  a' = exp({signif(row_LH_CH$intercept, 6)}) = {signif(a2, 6)}\n"))
    cat(glue("  b' = {signif(b2, 6)}\n"))
    cat(glue("  => LH = {signif(a2, 6)} * CH^{signif(b2, 6)}\n"))
  }
  
} else {
  
  # If you ever switch off logs, we fit linear relations on the original scale.
  # (No power-law printing in this branch.)
  reg_tbl <- bind_rows(
    fit_pair_lm(df_allom, "CH", "CF", "CF ~ CH"),
    fit_pair_lm(df_allom, "CF", "LF", "LF ~ CF"),  # swapped
    fit_pair_lm(df_allom, "CH", "LH", "LH ~ CH")   # swapped
  )
  
  print(reg_tbl)
  write.csv(reg_tbl, file.path(FIG_DIR, "allometry_regressions_R2.csv"), row.names = FALSE)
  
  plot_pair(df_allom, "CH", "CF", "CH", "CF", "reg_CF_CH.png")
  plot_pair(df_allom, "CF", "LF", "CF", "LF", "reg_LF_CF.png")
  plot_pair(df_allom, "CH", "LH", "CH", "LH", "reg_LH_CH.png")
}

message("Done. Outputs written to: ", normalizePath(FIG_DIR))

# ============================================================
# Alligators: SVL–TL allometry regression (log–log)
# ============================================================

DATA_XLSX_ALLIG <- "Data/alligators_woodward.xlsx"
FIG_DIR_ALLIG   <- "Figures/Alligators"
if (!dir.exists(FIG_DIR_ALLIG)) dir.create(FIG_DIR_ALLIG, recursive = TRUE)

df_allig_raw <- read_excel(DATA_XLSX_ALLIG)

# Keep consistent with your EVT pipeline:
# Deform == 1 or 3: tail broken => TL structurally missing
df_allig_allom <- df_allig_raw %>%
  mutate(
    SVL = suppressWarnings(as.numeric(SVL)),
    TL  = suppressWarnings(as.numeric(TL)),
    TL  = ifelse(Deform %in% c(1, 3), NA_real_, TL)
  ) %>%
  transmute(
    specimen = row_number(),
    SVL, TL
    # optionally restrict, e.g. males only:
    # Sex = Sex
  ) %>%
  mutate(
    log_SVL = ifelse(is.finite(SVL) & SVL > 0, log(SVL), NA_real_),
    log_TL  = ifelse(is.finite(TL)  & TL  > 0, log(TL),  NA_real_)
  )

# If you want to enforce "complete & non-deformed" explicitly:
# df_allig_allom <- df_allig_allom %>% filter(is.finite(log_SVL), is.finite(log_TL))

# Plot helper that lets you choose output dir (so we don't reuse FIG_DIR from sauropods)
plot_pair_to <- function(dat, x, y, xlab, ylab, fig_dir, fname) {
  d <- dat %>% filter(is.finite(.data[[x]]), is.finite(.data[[y]]))
  p <- ggplot(d, aes(x = .data[[x]], y = .data[[y]])) +
    geom_point(alpha = 0.75, size = 2) +
    geom_smooth(method = "lm", se = FALSE) +
    labs(x = xlab, y = ylab) +
    theme_science
  ggsave(file.path(fig_dir, fname), p, dpi = 600, w = 6.2, h = 4.8, units = "in")
  p
}

# Choose scale using your existing flag
if (isTRUE(use_log)) {
  reg_allig_tbl <- bind_rows(
    fit_pair_lm(df_allig_allom, "log_SVL", "log_TL", "log(TL) ~ log(SVL)")
  )
  print(reg_allig_tbl)
  write.csv(reg_allig_tbl, file.path(FIG_DIR_ALLIG, "allometry_regression_logTL_logSVL_R2.csv"),
            row.names = FALSE)
  
  plot_pair_to(df_allig_allom, "log_SVL", "log_TL",
               "log(SVL)", "log(TL)", FIG_DIR_ALLIG, "reg_logTL_logSVL.png")
  
} else {
  reg_allig_tbl <- bind_rows(
    fit_pair_lm(df_allig_allom, "SVL", "TL", "TL ~ SVL")
  )
  print(reg_allig_tbl)
  write.csv(reg_allig_tbl, file.path(FIG_DIR_ALLIG, "allometry_regression_TL_SVL_R2.csv"),
            row.names = FALSE)
  
  plot_pair_to(df_allig_allom, "SVL", "TL",
               "SVL (cm)", "TL (cm)", FIG_DIR_ALLIG, "reg_TL_SVL.png")
}

message("Alligator SVL–TL regression outputs written to: ", normalizePath(FIG_DIR_ALLIG))
