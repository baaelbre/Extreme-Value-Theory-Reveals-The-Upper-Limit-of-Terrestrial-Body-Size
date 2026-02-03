# ============================================================
# Title: Extreme Value Theory Reveals the Upper Limit of 
#        Terrestrial Body Size: Alligators
# Author: Bastiaan A. Van Velthoven
# ============================================================
library(readxl)     # read in excel files
library(dplyr)      # data manipulation
library(tidyverse)  # data manipulation
#library(purrr)
library(ggplot2)    # visualization
library(evd)        # extreme value theory
library(scales)
library(glue)
library(forcats)
library(grid)
library(MASS)
library(copula)

set.seed(42)

# ------------------------------------------------------------
# Load compiled C code (the .dll/.so file acts like a library
# upon using dyn.load)
# This is code from the evd package (Github page cran/evd) that contains the
# negative log likelihood for bivariate logistic dependence
# ------------------------------------------------------------
citation("evd")
# system("R CMD SHLIB bvpot.c") # if you need to compile the C code, with Rtools
# this creates a .dll and a .o file

dyn.load("bvpot.dll")
is.loaded("nllbvclog")

# ---------------------------
# Directories & theming
# ---------------------------
FIG_DIR <- "Figures/Alligators"
if (!dir.exists(FIG_DIR)) dir.create(FIG_DIR, recursive = TRUE)

theme_science <- theme_minimal(base_family = "Arial", base_size = 12) +
  theme(
    axis.title       = element_text(size = 14, face = "bold"),
    axis.text        = element_text(size = 12),
    legend.title     = element_text(size = 10, face = "bold"),
    legend.text      = element_text(size = 10),
    panel.grid.major = element_line(color = "gray80", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.6),
    axis.ticks.length= unit(0.20, "cm"),
    axis.ticks       = element_line(color = "black", linewidth = 0.4),
    plot.margin      = margin(5, 5, 5, 5),
    legend.position  = "right"
  )

# ---------------------------
# Global settings
# ---------------------------
ci_level      <- 0.90     # confidence level
u_lo   <- 0.60     # lower quantile to start scanning thresholds
u_hi   <- 0.99     # upper quantile
u_n      <- 50       # how many to check (mainly for plotting purposes)
min_ex    <- 20       # how many exceedances do you want at least for ML fit

trait_names   <- c("SVL","TL")  # SVL: snout-vent length, TL: total length
q_opt  <- c(SVL = 0.94, TL = 0.94)  # anchor quantiles (6% upper ones)

# Bootstrap settings
B_boot <- 1000L   # parametric bivariate bootstrap replicates
B_KS   <- 1000L   # parametric KS bootstrap per margin

# ============================================================
# Read the data
# ============================================================

DATA_XLSX <- "Data/alligators_woodward.xlsx"
df_raw    <- read_excel(DATA_XLSX)

# Deform == 1 or 3: tail broken ⇒ TL structurally missing (numeric NA); SVL kept.
df <- df_raw %>%
  mutate(
    SVL = as.numeric(SVL),
    TL  = as.numeric(TL),
    TL  = ifelse(Deform %in% c(1, 3), NA_real_, TL)
  ) %>%
  transmute(
    specimen = row_number(),
    SVL, TL
  )

# ============================================================
# Completeness diagnostics
# ============================================================

compl_tbl <- trait_names %>%
  set_names() %>%
  map_df(function(tr) {
    v <- df[[tr]]
    tibble(
      trait        = tr,
      n_total      = length(v),
      n_obs        = sum(!is.na(v)),
      completeness = mean(!is.na(v))
    )
  })
compl_tbl

p_compl <- compl_tbl %>%
  mutate(trait = fct_inorder(trait)) %>%
  ggplot(aes(trait, completeness)) +
  geom_col(fill = "#3B82F6") +
  geom_text(aes(label = percent(completeness, accuracy = 0.1)),
            vjust = -0.2, size = 4) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1.10)) +
  labs(y = "Completeness") +
  theme_science

ggsave(file.path(FIG_DIR, "completeness_SVL_TL.png"), p_compl,
       dpi = 600, w = 6.0, h = 4.2, units = "in")
p_compl # SVL rarely missing, TL in 4.2% of the cases.

# ============================================================
# Log-transform and split patterns (complete / SVL-only / TL-only)
# ============================================================

df_log <- df %>%
  mutate(
    log_SVL = ifelse(is.finite(SVL) & SVL > 0, log(SVL), NA_real_),
    log_TL  = ifelse(is.finite(TL)  & TL  > 0, log(TL),  NA_real_)
  )

# Complete cases: both SVL and TL observed
df_complete <- df_log %>%
  filter(!is.na(log_SVL), !is.na(log_TL))

# SVL only
df_svl_only <- df_log %>%
  filter(!is.na(log_SVL), is.na(log_TL))

# TL only
df_tl_only <- df_log %>%
  filter(is.na(log_SVL), !is.na(log_TL))

cat("\nRow counts by observation pattern:\n")
cat("  complete (SVL & TL) :", nrow(df_complete), "\n")
cat("  SVL only            :", nrow(df_svl_only), "\n")
cat("  TL only             :", nrow(df_tl_only), "\n")

# ------------------------------------------------------------
# Thresholds per trait (tweak q_opt to pick optimal ones)
# ------------------------------------------------------------

u0_by_trait <- setNames(numeric(length(trait_names)), trait_names)
for (tr in trait_names) {
  v_log <- df_log[[paste0("log_", tr)]]
  u0_by_trait[tr] <- as.numeric(quantile(v_log, q_opt[tr], na.rm = TRUE))
}
print(u0_by_trait)       # thresholds on log scale (for the analysis)
print(exp(u0_by_trait))  # thresholds on original scale

# ============================================================
# 4B Pairwise tail regimes & missing/observed/exceed/censored
#     (SVL, TL) — quadrant plot on log-scale
#     Legend uses tuples (E1,E2,O1,O2)
# ============================================================

u1 <- u0_by_trait["SVL"]
u2 <- u0_by_trait["TL"]

x1 <- df_log$log_SVL  # trait 1 = SVL
x2 <- df_log$log_TL   # trait 2 = TL

# Observed / missing indicators
O1 <- is.finite(x1)
O2 <- is.finite(x2)
M1 <- !O1
M2 <- !O2

# Exceedance indicators (only meaningful if observed)
E1 <- O1 & (x1 > u1)
E2 <- O2 & (x2 > u2)

# Convenience: both observed
both_obs <- O1 & O2

# Helper: tuple label in figure order (E1,E2,O1,O2)
tuple <- function(e1, e2, o1, o2) sprintf("(%d,%d,%d,%d)", e1, e2, o1, o2)

make_pair_plot_alligators <- function(t1 = "SVL", t2 = "TL") {
  fn <- file.path(FIG_DIR, glue("pair_scatter_{t1}_{t2}_final.png"))
  
  # Axis ranges based on observed data
  xr <- range(x1[O1], na.rm = TRUE)
  yr <- range(x2[O2], na.rm = TRUE)
  
  pad <- function(r, p = 0.03) {
    w <- diff(r)
    c(r[1] - p * w, r[2] + p * w)
  }
  xr <- pad(xr); yr <- pad(yr)
  x_min <- xr[1]; x_max <- xr[2]
  y_min <- yr[1]; y_max <- yr[2]
  
  # ---------------------------
  # 1) Points: both traits observed (O1=O2=1)
  # ---------------------------
  dd_dot <- tibble(
    specimen = df$specimen[both_obs],
    x        = x1[both_obs],
    y        = x2[both_obs],
    kind     = case_when(
      E1[both_obs] & E2[both_obs]  ~ tuple(1, 1, 1, 1),  # joint exceed
      E1[both_obs] & !E2[both_obs] ~ tuple(1, 0, 1, 1),  # SVL exceed, TL observed subthr
      !E1[both_obs] & E2[both_obs] ~ tuple(0, 1, 1, 1),  # TL exceed, SVL observed subthr
      TRUE                          ~ tuple(0, 0, 1, 1)  # both observed subthr
    )
  )
  
  # ---------------------------
  # 2) Segments: partially observed extremes (one exceed, partner missing)
  # ---------------------------
  
  # Vertical segments: SVL exceed, TL missing => (1,0,1,0)
  dd_v <- tibble(
    specimen = df$specimen[E1 & M2],
    x        = x1[E1 & M2],
    y0       = y_min,
    y1       = y_max,
    kind     = tuple(1, 0, 1, 0)
  )
  
  # Horizontal segments: TL exceed, SVL missing => (0,1,0,1)
  dd_h <- tibble(
    specimen = df$specimen[M1 & E2],
    y        = x2[M1 & E2],
    x0       = x_min,
    x1       = x_max,
    kind     = tuple(0, 1, 0, 1)
  )
  
  # ---------------------------
  # Colour map keyed by (E1,E2,O1,O2)
  # ---------------------------
  keys <- c(
    tuple(1, 1, 1, 1),
    tuple(1, 0, 1, 1),
    tuple(0, 1, 1, 1),
    tuple(0, 0, 1, 1),
    tuple(1, 0, 1, 0),
    tuple(0, 1, 0, 1)
  )
  
  vals <- c(
    "#DC2626",  # (1,1,1,1) joint exceed, both observed
    "#2563EB",  # (1,0,1,1) SVL exceed, TL observed subthr
    "#059669",  # (0,1,1,1) TL exceed, SVL observed subthr
    "grey70",   # (0,0,1,1) both observed subthr
    "#2563EB",  # (1,0,1,0) SVL exceed, TL missing (segment)
    "#059669"   # (0,1,0,1) TL exceed, SVL missing (segment)
  )
  
  col_map <- setNames(vals, keys)
  
  # Legend order
  legend_breaks <- c(
    tuple(0, 0, 1, 1),
    tuple(0, 1, 1, 1),
    tuple(1, 0, 1, 1),
    tuple(1, 1, 1, 1),
    tuple(0, 1, 0, 1),
    tuple(1, 0, 1, 0)
  )
  
  # ---------------------------
  # Background tail regions (based on thresholds u1,u2)
  # ---------------------------
  p <- ggplot() +
    annotate(
      "rect",
      xmin = x_min, xmax = u1,
      ymin = u2,    ymax = y_max,
      fill = col_map[tuple(0, 1, 1, 1)], alpha = 0.04
    ) +
    annotate(
      "rect",
      xmin = u1,    xmax = x_max,
      ymin = y_min, ymax = u2,
      fill = col_map[tuple(1, 0, 1, 1)], alpha = 0.04
    ) +
    annotate(
      "rect",
      xmin = u1,    xmax = x_max,
      ymin = u2,    ymax = y_max,
      fill = col_map[tuple(1, 1, 1, 1)], alpha = 0.07
    ) +
    
    geom_point(
      data = dd_dot,
      aes(x = x, y = y, color = kind),
      size = 2.6, alpha = 0.95
    ) +
    geom_segment(
      data = dd_v,
      aes(x = x, xend = x, y = y0, yend = y1, color = kind),
      linewidth = 0.9, alpha = 0.95
    ) +
    geom_segment(
      data = dd_h,
      aes(x = x0, xend = x1, y = y, yend = y, color = kind),
      linewidth = 0.9, alpha = 0.95
    ) +
    
    geom_vline(xintercept = u1, linetype = "dashed", color = "red") +
    geom_hline(yintercept = u2, linetype = "dashed", color = "red") +
    
    scale_color_manual(
      values = col_map,
      breaks = legend_breaks,
      name   = "(E1,E2,O1,O2)"
    ) +
    coord_cartesian(xlim = xr, ylim = yr, expand = FALSE) +
    labs(
      x = glue("log({t1} [cm])"),
      y = glue("log({t2} [cm])")
    ) +
    theme_science
  
  print(p)
  ggsave(fn, p, dpi = 600, w = 6.8, h = 5.6, units = "in")
  message("Saved: ", normalizePath(fn))
}

# Create the plot
make_pair_plot_alligators("SVL", "TL")

# ============================================================
# 6b. Empirical near-endpoint *unconditional* survival curves
#      using empirical maxima as pseudo-endpoints
#      Ŝ(t) = P(Y > y_max - t)
# ============================================================

# All log-data (non-missing) per margin
y_SVL_all <- df_log$log_SVL[!is.na(df_log$log_SVL)]
y_TL_all  <- df_log$log_TL[ !is.na(df_log$log_TL )]

# Helper: build Ŝ(t) = P(Y > y_max - t) for one trait (on log-scale)
make_endpoint_sf_curve_uncond <- function(y, u, name) {
  y_all <- y[is.finite(y)]
  n_all <- length(y_all)
  
  y_tail <- y_all[y_all > u]
  y_tail <- sort(y_tail)
  n_tail <- length(y_tail)
  
  if (n_tail < 5L) {
    warning(glue("Too few exceedances for {name} (n_tail < 5)."))
    return(tibble())
  }
  
  # pseudo-endpoint = empirical max within tail
  y_max <- max(y_tail)
  
  # indices 1..(n_tail-1) so Ŝ never equals 0 and t_hat > 0
  k     <- seq_len(n_tail - 1L)
  t_hat <- y_max - y_tail[k]          # distance to empirical max (log-scale)
  
  # conditional survival given Y > u:
  #   S_cond(t) = P(Y > y_tail[k] | Y > u) = (n_tail - k)/n_tail
  S_cond <- (n_tail - k) / n_tail
  
  # tail fraction P(Y > u) ≈ n_tail / n_all
  p_tail <- n_tail / n_all
  
  # unconditional survival:
  #   S_uncond(t) = P(Y > u) * S_cond(t)
  S_uncond <- p_tail * S_cond
  
  tibble(
    t_hat = t_hat,
    S_hat = S_uncond,
    Trait = name
  )
}

sf_SVL <- make_endpoint_sf_curve_uncond(y_SVL_all, u1, "SVL")
sf_TL  <- make_endpoint_sf_curve_uncond(y_TL_all,  u2, "TL")

sf_all <- bind_rows(sf_SVL, sf_TL)

trait_cols <- c(
  "SVL" = "#377eb8",  # blue
  "TL"  = "#1b9e77"   # green
)

p_sf_endpoint <- ggplot(sf_all, aes(x = t_hat, y = S_hat, colour = Trait)) +
  geom_point(size = 2, alpha = 0.9) +
  scale_colour_manual(values = trait_cols, name = "Trait") +
  scale_x_log10() +
  scale_y_log10() +
  labs(
    x = expression(hat(t) == y[max] - y),
    y = expression(hat(S)(t) == P(Y > y[max] - t))
  ) +
  theme_science

print(p_sf_endpoint)

ggsave(
  file.path(FIG_DIR, "SVL_TL_empirical_survival_endpoint_uncond.png"),
  p_sf_endpoint,
  dpi   = 600,
  w     = 6.0,
  h     = 4.5,
  units = "in"
)


# ------------------------------------------------------------
# Helper functions for univariate diagnostics
# ------------------------------------------------------------

# generates a list of u_n thresholds that runs from u_lo to u_hi
make_thresholds <- function(y, q_from = u_lo, q_to = u_hi, n = u_n) {
  rng <- quantile(y, c(q_from, q_to), na.rm = TRUE)
  seq(rng[1], rng[2], length.out = n)
}

fit_gpd_at_u <- function(y, u) {
  fit <- evd::fpot(y, threshold = u, model = "gpd")
  par <- fit$estimate
  cov <- fit$var.cov
  list(
    fit      = fit,
    scale    = unname(par["scale"]),
    shape    = unname(par["shape"]),
    cov      = cov,
    n_exceed = sum(y > u) # this is of course the quantile.
  )
}

# the adjusted (threshold-independent) scale is defined as sigma_u - xi(u-u_0)
# which should be constant upon increasing u_0 to u (sigma_u is not), see Coles CH3
adj_scale_fun <- function(scale, xi, u, u0) {
  scale - xi * (u - u0)
}

diagnostic_plots <- function(y, u, sigma_hat, xi_hat, label) {
  exc   <- y[y > u] - u
  n     <- length(exc)
  probs <- ppoints(n)
  
  if (abs(xi_hat) > 1e-10) {
    theo_q <- u + sigma_hat / xi_hat * (probs^(-xi_hat) - 1)
  } else {
    theo_q <- u - sigma_hat * log(probs) # exponential limit.
  }
  
  dfqq <- data.frame(
    Theoretical = rev(theo_q),
    Empirical   = sort(y[y > u])
  )
  pqq <- ggplot(dfqq, aes(Theoretical, Empirical)) +
    geom_point(color = "steelblue") +
    geom_abline(slope = 1, intercept = 0,
                linetype = "dashed", color = "red") +
    labs(title = glue("Q–Q: {label}")) +
    theme_science
  
  F_theo <- if (abs(xi_hat) > 1e-10) {
    1 - (1 + xi_hat * exc / sigma_hat)^(-1 / xi_hat)
  } else {
    1 - exp(-exc / sigma_hat)
  }
  dfpp <- data.frame(
    Theoretical = sort(F_theo),
    Empirical   = (1:n) / n
  )
  ppp <- ggplot(dfpp, aes(Theoretical, Empirical)) +
    geom_point(color = "darkgreen") +
    geom_abline(slope = 1, intercept = 0,
                linetype = "dashed", color = "red") +
    labs(title = glue("P–P: {label}")) +
    theme_science
  
  list(pqq = pqq, ppp = ppp)
}

# ------------------------------------------------------------
# KS goodness-of-fit via parametric bootstrap (univariate!)
# ------------------------------------------------------------

ks_boot <- function(y, u, scale_hat, shape_hat,B = 1000L) {

  y_exceed <- y[y > u]
  ex       <- y_exceed - u
  n_ex     <- length(ex)
  if (n_ex < 20L) {
    warning("Not enough exceedances for KS bootstrap (n_ex < 20) -> NA")
    return(list(D = NA_real_, p = NA_real_, D_boot = numeric(0L), n_exceed = n_ex))
  }
  
  # observed KS statistic under fitted GPD
  Fn_ex      <- ecdf(ex)
  ex_sorted  <- sort(ex)
  Fgpd_obs   <- evd::pgpd(ex_sorted, loc = 0, scale = scale_hat, shape = shape_hat)
  Fn_vals    <- Fn_ex(ex_sorted)
  D_obs      <- max(abs(Fn_vals - Fgpd_obs))
  
  # Parametric bootstrap under fitted GPD with refitting
  D_boot <- numeric(B)
  
  for (b in seq_len(B)) {
    ex_b <- evd::rgpd(n_ex, loc = 0, scale = scale_hat, shape = shape_hat)
    y_b  <- u + ex_b
    
    fit_b <- try(evd::fpot(y_b, threshold = u, model = "gpd"), silent = TRUE)
    if (inherits(fit_b, "try-error")) {
      D_boot[b] <- NA_real_
      next
    }
    
    sc_b <- unname(fit_b$estimate["scale"])
    xi_b <- unname(fit_b$estimate["shape"])
    
    Fn_b       <- ecdf(ex_b)
    ex_b_sort  <- sort(ex_b)
    Fgpd_b     <- evd::pgpd(ex_b_sort, loc = 0, scale = sc_b, shape = xi_b)
    Fn_b_vals  <- Fn_b(ex_b_sort)
    D_boot[b]  <- max(abs(Fn_b_vals - Fgpd_b))
  }
  
  D_boot <- D_boot[is.finite(D_boot)]
  if (!length(D_boot)) {
    warning("All KS bootstrap replicates failed; returning NA.")
    return(list(D = D_obs, p = NA_real_, D_boot = numeric(0L), n_exceed = n_ex))
  }
  
  p_boot <- mean(D_boot >= D_obs)
  
  list(D = D_obs, p = p_boot, D_boot = D_boot, n_exceed = n_ex)
}

# ------------------------------------------------------------
# Diagnostics, KS, and Wald p-values, per trait.
# ------------------------------------------------------------
run_trait <- function(tr_key) {
  v_raw <- df[[tr_key]] |> suppressWarnings(as.numeric())
  y     <- log(v_raw[is.finite(v_raw) & v_raw > 0])
  ylab  <- glue("log {tr_key}")
  stopifnot(length(y) >= 40)
  
  q_anchor <- q_opt[tr_key]
  u_seq    <- make_thresholds(y)
  u0       <- quantile(y, q_anchor, na.rm = TRUE) |> as.numeric()
  
  # GPD fit at u0
  fit0  <- fit_gpd_at_u(y, u0)
  xi0   <- fit0$shape
  sc0   <- fit0$scale
  cov0  <- fit0$cov
  xi_se <- sqrt(cov0[2, 2])
  sigma_se <- sqrt(cov0[1, 1])
  
  # Wald test for xi < 0 (one-sided)
  z_wald <- xi0 / xi_se
  p_wald <- pnorm(z_wald)  # H1: xi < 0
  
  # KS goodness-of-fit
  ks_res <- ks_boot(
    y          = y,
    u          = u0,
    scale_hat  = sc0,
    shape_hat  = xi0,
    B          = B_KS
  )
  
  zcrit <- qnorm(1 - (1 - ci_level) / 2)
  
  df_scan <- tibble(
    u        = u_seq,
    shape    = NA_real_,
    shape_lo = NA_real_,
    shape_hi = NA_real_,
    adj      = NA_real_,
    adj_lo   = NA_real_,
    adj_hi   = NA_real_
  )
  
  for (i in seq_along(u_seq)) {
    u <- u_seq[i]
    ex <- y[y > u] - u
    if (length(ex) < min_ex) next
    
    out <- try(fit_gpd_at_u(y, u), silent = TRUE)
    if (inherits(out, "try-error")) next
    
    se_shape <- sqrt(out$cov[2, 2])
    df_scan$shape[i]    <- out$shape
    df_scan$shape_lo[i] <- out$shape - zcrit * se_shape
    df_scan$shape_hi[i] <- out$shape + zcrit * se_shape
    
    df_scan$adj[i]      <- adj_scale_fun(out$scale, out$shape, u, u0)
    var_adj <- out$cov[1, 1] +
      (u - u0)^2 * out$cov[2, 2] -
      2 * (u - u0) * out$cov[1, 2]
    se_adj  <- sqrt(max(var_adj, 0))
    df_scan$adj_lo[i] <- df_scan$adj[i] - zcrit * se_adj
    df_scan$adj_hi[i] <- df_scan$adj[i] + zcrit * se_adj
  }
  
  p_shape <- ggplot(df_scan, aes(u, shape)) +
    geom_point(color = "blue") +
    geom_errorbar(aes(ymin = shape_lo, ymax = shape_hi),
                  width = 0.03, color = "blue") +
    geom_vline(xintercept = u0, color = "red", linetype = "dashed") +
    geom_hline(yintercept = xi0, color = "red", linetype = "dashed") +
    labs(
      x = glue("Threshold ({ylab})"),
      y = "Shape (xi)") +
    theme_science
  
  p_adj <- ggplot(df_scan, aes(u, adj)) +
    geom_point(color = "blue") +
    geom_errorbar(aes(ymin = adj_lo, ymax = adj_hi),
                  width = 0.03, color = "blue") +
    geom_vline(xintercept = u0, color = "red", linetype = "dashed") +
    geom_hline(yintercept = sc0, color = "red", linetype = "dashed") +
    labs(
      x = glue("Threshold ({ylab})"),
      y = "Adjusted scale") +
    theme_science
  
  ggsave(file.path(FIG_DIR, glue("{tr_key}_log_shape_stability.png")),
         p_shape, dpi = 600, w = 6.2, h = 4.2, units = "in")
  ggsave(file.path(FIG_DIR, glue("{tr_key}_log_adj_scale_stability.png")),
         p_adj,   dpi = 600, w = 6.2, h = 4.2, units = "in")
  
  di <- diagnostic_plots(y, u0, sc0, xi0, ylab)
  ggsave(file.path(FIG_DIR, glue("{tr_key}_log_qq.png")),
         di$pqq, dpi = 600, w = 6.2, h = 4.2, units = "in")
  ggsave(file.path(FIG_DIR, glue("{tr_key}_log_pp.png")),
         di$ppp, dpi = 600, w = 6.2, h = 4.2, units = "in")

  
  cat(glue("\n{tr_key}: Wald test xi<0 at u0={round(u0,3)}: xî={round(xi0,3)}, se={round(xi_se,3)}, z={round(z_wald,3)}, p={signif(p_wald,3)}"))
  cat(glue("\n{tr_key}: KS bootstrap: D={round(ks_res$D,3)}, p_boot={signif(ks_res$p,3)}, n_exceed={ks_res$n_exceed}\n"))
  
  tibble(
    trait           = tr_key,
    u0              = u0,
    xi_hat_diag     = xi0,
    sigma_hat_diag  = sc0,
    xi_se           = xi_se,
    sigma_se        = sigma_se,
    wald_z_xi       = z_wald,
    wald_p_xi_neg   = p_wald,
    ks_D            = ks_res$D,
    ks_p_boot       = ks_res$p,
    ks_n_exceed     = ks_res$n_exceed
  )
}

diag_SVL <- run_trait("SVL")
diag_TL  <- run_trait("TL")
marg_diag_tbl <- bind_rows(diag_SVL, diag_TL)
print(marg_diag_tbl)

# ============================================================
# Bivariate POT with censored likelihood (complete cases)
# ============================================================

# Complete-case matrix for censored bivariate logistic POT (this is handled in bvpot.c)
X_biv <- df_complete %>%
  transmute(
    log_SVL = log_SVL,
    log_TL  = log_TL
  ) %>%
  as.matrix()

n_biv <- nrow(X_biv)
cat("\nBivariate sample size (complete cases SVL & TL):", n_biv, "\n")
u1 <- u0_by_trait["SVL"]
u2 <- u0_by_trait["TL"]
u_vec <- c(u1, u2)

quad_tab <- table(
  SVL_exceed = X_biv[, 1] > u1,
  TL_exceed  = X_biv[, 2] > u2
)
cat("\nQuadrant table (on log scale, complete cases only):\n")
print(quad_tab)

# Missing patterns (not part of the nllbvclog routine!)

# Only cases with an exceedance and the other margin missing:
y1_E_M2 <- df_svl_only$log_SVL[df_svl_only$log_SVL > u1]  # SVL exceeds, TL missing
y2_M1_E <- df_tl_only$log_TL[df_tl_only$log_TL > u2]      # TL exceeds, SVL missing

cat("\nCounts of missing-partner exceedances:\n")
cat("  SVL exceed, TL missing (E1M2):", length(y1_E_M2), "\n")
cat("  TL exceed, SVL missing (M1E2):", length(y2_M1_E), "\n") 
# this is of course quite unlikely, but for sauropods finding a femur or a humerus
# is more or less equally probable ...


make_ll_biv <- function(x, u, cshape = FALSE, cscale = FALSE) {
  # Use sep.bvdata from evd to build (x1, x2, nn, n, thid, lambda, ...)
  sep.bvdata <- getFromNamespace("sep.bvdata", "evd")
  spx        <- sep.bvdata(x = x, method = "cpot", u = u)
  
  function(theta) {
    scale1 <- theta["scale1"]
    shape1 <- theta["shape1"]
    scale2 <- theta["scale2"]
    shape2 <- theta["shape2"]
    dep    <- theta["dep"]
    
    if (isTRUE(cshape)) shape2 <- shape1
    if (isTRUE(cscale)) scale2 <- scale1
    
    # Call compiled C routine: nllbvclog(...)
    out <- .C("nllbvclog",
              data1  = as.double(spx$x1),
              data2  = as.double(spx$x2),
              nn     = as.integer(spx$nn),
              n      = as.integer(spx$n),
              thid   = as.double(spx$thdi),     # NOTE: double* in C
              lambda = as.double(spx$lambda),
              dep    = as.double(dep),
              scale1 = as.double(scale1),
              shape1 = as.double(shape1),
              scale2 = as.double(scale2),
              shape2 = as.double(shape2),
              dns    = as.double(0)
    )
    
    # C routine returns the negative log-likelihood in dns, so to avoid mistakes
    # we negate the output
    -out$dns
  }
}

ll_cens_log <- make_ll_biv(
  x      = X_biv,
  u      = u_vec,
  cshape = FALSE,
  cscale = FALSE
)

# ---------------------------
# Univariate GP contributions from missings
# ---------------------------

lgp_exceed_evd <- function(y, u, scale, shape) {
  ex <- y - u
  evd::dgpd(ex, loc = 0, scale = scale, shape = shape, log = TRUE)
}

make_ll_uni <- function(y1_E_M2, y2_M1_E, u1, u2) {
  function(theta) {
    scale1 <- theta["scale1"]
    shape1 <- theta["shape1"]
    scale2 <- theta["scale2"]
    shape2 <- theta["shape2"]
    
    ll1 <- if (length(y1_E_M2)) {
      sum(lgp_exceed_evd(y1_E_M2, u1, scale1, shape1))
    } else 0
    
    ll2 <- if (length(y2_M1_E)) {
      sum(lgp_exceed_evd(y2_M1_E, u2, scale2, shape2))
    } else 0
    
    ll1 + ll2
  }
}

ll_missing <- make_ll_uni(y1_E_M2, y2_M1_E, u1, u2)

# Full log-likelihood: censored bivariate POT + univariate missings
loglik_full <- function(theta) {
  ll_cens_log(theta) + ll_missing(theta)
}

# Wrapper for optim
negloglik_full <- function(p, common_shape = FALSE) {
  if (!common_shape) {
    # Free shapes: theta = (scale1, shape1, scale2, shape2, dep)
    theta <- c(
      scale1 = p[1],
      shape1 = p[2],
      scale2 = p[3],
      shape2 = p[4],
      dep    = p[5]
    )
  } else {
    # H0: xi1 = xi2 = xi
    theta <- c(
      scale1 = p[1],
      shape1 = p[2],  # xi
      scale2 = p[3],
      shape2 = p[2],  # same xi
      dep    = p[4]
    )
  }
  -loglik_full(theta)
}

# ============================================================
# Maximum likelihood: free shapes vs common shape
# ============================================================
sigma1_hat_uni <- as.numeric(diag_SVL$sigma_hat_diag)
xi1_hat_uni    <- as.numeric(diag_SVL$xi_hat_diag)

sigma2_hat_uni <- as.numeric(diag_TL$sigma_hat_diag)
xi2_hat_uni    <- as.numeric(diag_TL$xi_hat_diag)

# Free tail
# initialization
start_free <- c(
  scale1 = sigma1_hat_uni,
  shape1 = xi1_hat_uni,
  scale2 = sigma2_hat_uni,
  shape2 = xi2_hat_uni,
  dep    = 0.8
)

opt_free <- optim(
  par     = unname(start_free),
  fn      = negloglik_full,
  common_shape = FALSE,
  method  = "Nelder-Mead"
)

par_free <- opt_free$par
names(par_free) <- names(start_free)

ell_hat_free <- -opt_free$value

cat("\nBivariate MLEs (full llh, free shapes):\n")
print(par_free)
cat("Maximized log-likelihood (free) =", ell_hat_free, "\n")

# Common tail (H0: xi1 = xi2 = xi)
# initialization
start_eq <- c(
  scale1 = sigma1_hat_uni,
  xi     = mean(c(xi1_hat_uni, xi2_hat_uni)),
  scale2 = sigma2_hat_uni,
  dep    = 0.8
)

opt_eq <- optim(
  par     = unname(start_eq),
  fn      = negloglik_full,
  common_shape = TRUE,
  method  = "Nelder-Mead"
)

par_eq_phi <- opt_eq$par
names(par_eq_phi) <- names(start_eq)

ell_hat_eq <- -opt_eq$value

# Reconstruct full θ under H0
par_eq <- c(
  scale1 = par_eq_phi["scale1"],
  shape1 = par_eq_phi["xi"],
  scale2 = par_eq_phi["scale2"],
  shape2 = par_eq_phi["xi"],
  dep    = par_eq_phi["dep"]
)

cat("\nCommon tail shape MLEs (full llh, H0: xi1 = xi2):\n")
print(par_eq)
cat("Maximized log-likelihood (H0) =", ell_hat_eq, "\n")

# Likelihood ratio test for H0: xi1 = xi2
LR_shape <- 2 * (ell_hat_free - ell_hat_eq)
df_LR    <- 1
p_LR     <- 1 - pchisq(LR_shape, df = df_LR)

cat("\nLikelihood ratio test for tail equality H0: xi1 = xi2\n")
cat("  LR statistic =", LR_shape, " with df =", df_LR, "\n")
cat("  p-value      =", p_LR, "\n")

# ============================================================
# 9. Endpoint estimation under common tail shape (H0)
# ============================================================
xi_eq_hat  <- unname(par_eq_phi["xi"])      # common tail index (xi)
sigma1_hat <- unname(par_eq_phi["scale1"])  # SVL scale at u1
sigma2_hat <- unname(par_eq_phi["scale2"])  # TL  scale at u2

n_SVL_total <- sum(!is.na(df_log$log_SVL))
n_TL_total  <- sum(!is.na(df_log$log_TL))

n_SVL_exceed_u1 <- sum(df_log$log_SVL > u1, na.rm = TRUE)
n_TL_exceed_u2  <- sum(df_log$log_TL  > u2, na.rm = TRUE)

tail_frac_SVL <- n_SVL_exceed_u1 / n_SVL_total
tail_frac_TL  <- n_TL_exceed_u2  / n_TL_total

cat("\nEmpirical tail fractions:\n")
cat("  P(SVL > u1) ≈", tail_frac_SVL, "(u1 =", u1, ", exp(u1) ≈", exp(u1), "cm)\n")
cat("  P(TL  > u2) ≈", tail_frac_TL,  "(u2 =", u2, ", exp(u2) ≈", exp(u2), "cm)\n")

# ------------------------------------------------------------
# 9.3 Parametric bivariate bootstrap under fitted logistic (Gumbel) dependence
#
#   Complete cases: (log_SVL, log_TL) with dependence via a Gumbel copula,
#   and with each margin generated from a mixture:
#     * below-threshold values are sampled from the empirical sub-threshold data,
#     * above-threshold values are generated from the fitted GPD exceedance law.
# - Missing-pattern cases:
#     * SVL-only sample: simulate SVL marginal only, then keep exceedances (E1,M2)
#     * TL-only sample : simulate TL marginal only, then keep exceedances (M1,E2)
#
# What we refit:
# - For each bootstrap dataset, refit the constrained bivariate POT model
#   (xi shared across margins) to obtain (sigma1, sigma2, xi) and endpoints.
# ------------------------------------------------------------

# Sub-threshold samples for realism (on log scale)
sub_SVL_cc   <- df_complete$log_SVL[df_complete$log_SVL <= u1]
sub_TL_cc    <- df_complete$log_TL[df_complete$log_TL <= u2]
sub_SVL_only <- df_svl_only$log_SVL[df_svl_only$log_SVL <= u1]
sub_TL_only  <- df_tl_only$log_TL[df_tl_only$log_TL <= u2]

n_cc       <- nrow(df_complete)
n_svl_only <- nrow(df_svl_only)
n_tl_only  <- nrow(df_tl_only)

# Logistic EV parameter alpha from constrained bivariate fit; map to Gumbel copula 1/alpha
dep_name <- grep("^dep", names(par_eq), value = TRUE)[1]
alpha_logistic_hat <- unname(par_eq[dep_name])
theta_gumbel_hat   <- 1 / alpha_logistic_hat  # θ ≥ 1

gumbel_cop_hat <- gumbelCopula(param = theta_gumbel_hat, dim = 2)

# Storage for bootstrap draws (constrained refits)
boot_sigma1 <- numeric(B_boot)
boot_sigma2 <- numeric(B_boot)
boot_xi     <- numeric(B_boot)
boot_y1star <- numeric(B_boot)  # log endpoint for SVL*
boot_y2star <- numeric(B_boot)  # log endpoint for TL*
boot_conv   <- logical(B_boot)

set.seed(2026)

for (b in seq_len(B_boot)) {
  
  # ---------- (1) Simulate complete cases via copula + marginal tail mixtures ----------
  U_cc  <- rCopula(n_cc, gumbel_cop_hat)
  Y1_cc <- numeric(n_cc)  # log_SVL
  Y2_cc <- numeric(n_cc)  # log_TL
  
  for (k in seq_len(n_cc)) {
    u1k <- U_cc[k, 1]
    u2k <- U_cc[k, 2]
    
    # SVL margin: bulk vs tail
    if (u1k <= 1 - tail_frac_SVL || length(sub_SVL_cc) == 0L) {
      Y1_cc[k] <- if (length(sub_SVL_cc)) sample(sub_SVL_cc, 1L, replace = TRUE) else u1
    } else {
      v1   <- (u1k - (1 - tail_frac_SVL)) / tail_frac_SVL
      exc1 <- evd::qgpd(v1, loc = 0, scale = sigma1_hat, shape = xi_eq_hat)
      Y1_cc[k] <- u1 + exc1
    }
    
    # TL margin: bulk vs tail
    if (u2k <= 1 - tail_frac_TL || length(sub_TL_cc) == 0L) {
      Y2_cc[k] <- if (length(sub_TL_cc)) sample(sub_TL_cc, 1L, replace = TRUE) else u2
    } else {
      v2   <- (u2k - (1 - tail_frac_TL)) / tail_frac_TL
      exc2 <- evd::qgpd(v2, loc = 0, scale = sigma2_hat, shape = xi_eq_hat)
      Y2_cc[k] <- u2 + exc2
    }
  }
  
  X_biv_boot <- cbind(Y1_cc, Y2_cc)
  
  # ---------- (2) Simulate SVL-only and TL-only marginals, retain only exceedances ----------
  if (n_svl_only > 0L) {
    U1_svl <- runif(n_svl_only)
    Y1_svl <- numeric(n_svl_only)
    for (k in seq_len(n_svl_only)) {
      u1k <- U1_svl[k]
      if (u1k <= 1 - tail_frac_SVL || length(sub_SVL_only) == 0L) {
        Y1_svl[k] <- if (length(sub_SVL_only)) sample(sub_SVL_only, 1L, replace = TRUE) else u1
      } else {
        v1   <- (u1k - (1 - tail_frac_SVL)) / tail_frac_SVL
        exc1 <- evd::qgpd(v1, loc = 0, scale = sigma1_hat, shape = xi_eq_hat)
        Y1_svl[k] <- u1 + exc1
      }
    }
    y1_E_M2_boot <- Y1_svl[Y1_svl > u1]
  } else {
    y1_E_M2_boot <- numeric(0L)
  }
  
  if (n_tl_only > 0L) {
    U2_tl <- runif(n_tl_only)
    Y2_tl <- numeric(n_tl_only)
    for (k in seq_len(n_tl_only)) {
      u2k <- U2_tl[k]
      if (u2k <= 1 - tail_frac_TL || length(sub_TL_only) == 0L) {
        Y2_tl[k] <- if (length(sub_TL_only)) sample(sub_TL_only, 1L, replace = TRUE) else u2
      } else {
        v2   <- (u2k - (1 - tail_frac_TL)) / tail_frac_TL
        exc2 <- evd::qgpd(v2, loc = 0, scale = sigma2_hat, shape = xi_eq_hat)
        Y2_tl[k] <- u2 + exc2
      }
    }
    y2_M1_E_boot <- Y2_tl[Y2_tl > u2]
  } else {
    y2_M1_E_boot <- numeric(0L)
  }
  
  # ---------- (3) Refit constrained model on bootstrap sample ----------
  ll_cens_log_boot <- make_ll_biv(
    x      = X_biv_boot,
    u      = u_vec,
    cshape = FALSE,
    cscale = FALSE
  )
  ll_missing_boot <- make_ll_uni(y1_E_M2_boot, y2_M1_E_boot, u1, u2)
  
  loglik_full_boot <- function(theta) ll_cens_log_boot(theta) + ll_missing_boot(theta)
  
  negloglik_full_boot <- function(p) {
    # p = (scale1, xi, scale2, dep) under constraint xi1=xi2=xi
    theta <- c(
      scale1 = p[1],
      shape1 = p[2],
      scale2 = p[3],
      shape2 = p[2],
      dep    = p[4]
    )
    -loglik_full_boot(theta)
  }
  
  start_boot <- c(
    scale1 = sigma1_hat,
    xi     = xi_eq_hat,
    scale2 = sigma2_hat,
    dep    = alpha_logistic_hat
  )
  
  opt_b <- try(
    optim(
      par    = unname(start_boot),
      fn     = negloglik_full_boot,
      method = "Nelder-Mead"
    ),
    silent = TRUE
  )
  
  if (inherits(opt_b, "try-error") || opt_b$convergence != 0) {
    boot_conv[b] <- FALSE
    next
  }
  
  boot_conv[b] <- TRUE
  phi_b <- opt_b$par
  names(phi_b) <- names(start_boot)
  
  sigma1_b <- phi_b["scale1"]
  xi_b     <- phi_b["xi"]
  sigma2_b <- phi_b["scale2"]
  
  boot_sigma1[b] <- sigma1_b
  boot_sigma2[b] <- sigma2_b
  boot_xi[b]     <- xi_b
  
  # Endpoints on log scale (finite only if xi_b < 0)
  boot_y1star[b] <- u1 - sigma1_b / xi_b
  boot_y2star[b] <- u2 - sigma2_b / xi_b
}

cat("\nParametric bootstrap (common-tail logistic GP):\n")
cat("  Successful replicates:", sum(boot_conv), "out of", B_boot, "\n")
boot_ok <- boot_conv & is.finite(boot_xi) & (boot_xi < 0) &
  is.finite(boot_y1star) & is.finite(boot_y2star)


alpha_half <- (1 - ci_level) / 2

ci_sigma1 <- quantile(boot_sigma1[boot_ok], probs = c(alpha_half, 1 - alpha_half), na.rm = TRUE)
ci_sigma2 <- quantile(boot_sigma2[boot_ok], probs = c(alpha_half, 1 - alpha_half), na.rm = TRUE)
ci_xi     <- quantile(boot_xi[boot_ok],     probs = c(alpha_half, 1 - alpha_half), na.rm = TRUE)

cat("\n", 100 * ci_level, "% bootstrap CIs for GP parameters (log scale):\n", sep = "")
cat("  sigma_SVL:", ci_sigma1[1], "to", ci_sigma1[2], "\n")
cat("  sigma_TL :", ci_sigma2[1], "to", ci_sigma2[2], "\n")
cat("  xi      :", ci_xi[1],     "to", ci_xi[2],     "\n")

# One-sided (upper) bootstrap bounds for log endpoints
ci_y1_boot_log_upper <- quantile(boot_y1star[boot_ok], probs = ci_level, na.rm = TRUE)
ci_y2_boot_log_upper <- quantile(boot_y2star[boot_ok], probs = ci_level, na.rm = TRUE)

# Back-transform to original (cm) scale
ci_y1_boot_orig_upper <- exp(ci_y1_boot_log_upper)
ci_y2_boot_orig_upper <- exp(ci_y2_boot_log_upper)

cat("\nOne-sided ", 100 * ci_level, "% upper bootstrap bounds for endpoints:\n", sep = "")
cat("  SVL* ≤ ", ci_y1_boot_orig_upper, " cm\n", sep = "")
cat("  TL*  ≤ ", ci_y2_boot_orig_upper, " cm\n", sep = "")

svl_boot_cm <- exp(boot_y1star[boot_ok])
tl_boot_cm  <- exp(boot_y2star[boot_ok])

boot_joint <- data.frame(SVL = svl_boot_cm, TL = tl_boot_cm)

svl_window <- as.numeric(quantile(svl_boot_cm, c(0.01, 0.99), na.rm = TRUE))
tl_window  <- as.numeric(quantile(tl_boot_cm,  c(0.01, 0.99), na.rm = TRUE))

inside_win <- with(
  boot_joint,
  SVL >= svl_window[1] & SVL <= svl_window[2] &
    TL  >= tl_window[1]  & TL  <= tl_window[2]
)

boot_joint_kde <- boot_joint[inside_win, , drop = FALSE]


kde_joint <- MASS::kde2d(
  x    = boot_joint_kde$SVL,
  y    = boot_joint_kde$TL,
  n    = 200,
  lims = c(svl_window[1], svl_window[2], tl_window[1], tl_window[2])
)

# Joint KDE mode = argmax_z on the KDE grid
ij_max <- which(kde_joint$z == max(kde_joint$z), arr.ind = TRUE)[1, ]
svl_mode <- kde_joint$x[ij_max[1]]
tl_mode  <- kde_joint$y[ij_max[2]]

cat("\n2D mode (bootstrap) for endpoints:\n")
cat("  SVL*_mode =", round(svl_mode, 1), " cm\n")
cat("  TL*_mode  =", round(tl_mode,  1), " cm\n")

# SVL Stokes
SVL_stokes <- 239
TL_stokes <- 450

# ============================================================
# 10. Exceedance probability for the Stokes alligator
# ============================================================

# --- 10.1 Helper: unconditional exceedance prob under a POT/GPD tail
# If y > u:
#   P(Y > y) = P(Y > u) * P(Y > y | Y > u)
#           ≈ p_u      * S_GPD(y - u; sigma, xi)
pot_exceed_prob <- function(y, u, sigma, xi, p_u) {
  if (!is.finite(y) || !is.finite(u) || !is.finite(sigma) || !is.finite(xi) || !is.finite(p_u))
    return(NA_real_)
  if (y <= u) return(NA_real_)  # outside tail model (by design)
  
  x <- y - u
  if (abs(xi) < 1e-8) {
    # exponential limit
    return(p_u * exp(-x / sigma))
  } else {
    term <- 1 + xi * x / sigma
    if (term <= 0) return(0)    # beyond finite endpoint when xi < 0
    return(p_u * term^(-1 / xi))
  }
}

# --- 10.2 Helper: Gumbel copula C(u,v) (closed form)
gumbel_C <- function(u, v, theta) {
  # theta >= 1
  a <- (-log(u))^theta
  b <- (-log(v))^theta
  exp(-(a + b)^(1 / theta))
}

# --- 10.3 Helper: joint exceedance from marginals + copula
# For U=F1(Y1), V=F2(Y2) with copula C:
#   P(Y1>y1, Y2>y2) = 1 - u - v + C(u,v),  where u=F1(y1), v=F2(y2)
joint_exceed_gumbel <- function(p1, p2, theta) {
  # p1 = P(Y1>y1), p2 = P(Y2>y2)
  u <- 1 - p1
  v <- 1 - p2
  eps <- 1e-12
  u <- min(max(u, eps), 1 - eps)
  v <- min(max(v, eps), 1 - eps)
  Cuv <- gumbel_C(u, v, theta)
  1 - u - v + Cuv
}


log_SVL_stokes <- log(SVL_stokes)
log_TL_stokes  <- log(TL_stokes)

cat("\nStokes alligator (log scale):\n")
cat("  log(SVL) =", log_SVL_stokes, "\n")
cat("  log(TL)  =", log_TL_stokes,  "\n")

# --- Marginal exceedance probs under fitted POT tails (requires Stokes > thresholds)
p_SVL_stokes <- pot_exceed_prob(
  y    = log_SVL_stokes,
  u    = u1,
  sigma= sigma1_hat,
  xi   = xi_eq_hat,
  p_u  = tail_frac_SVL
)

p_TL_stokes <- pot_exceed_prob(
  y    = log_TL_stokes,
  u    = u2,
  sigma= sigma2_hat,
  xi   = xi_eq_hat,
  p_u  = tail_frac_TL
)

cat("\nUnivariate exceedance (plug-in POT):\n")
cat("  P(SVL >", SVL_stokes, "cm) ≈", signif(p_SVL_stokes, 4), "\n")
cat("  P(TL  >", TL_stokes,  "cm) ≈", signif(p_TL_stokes,  4), "\n")

# --- Dependence from fitted logistic model (evd parameter alpha; Gumbel theta = 1/alpha)
alpha_logistic <- alpha_logistic_hat
theta_gumbel   <- 1 / alpha_logistic

cat("\nDependence (from bivariate POT fit):\n")
cat("  alpha (logistic) =", signif(alpha_logistic, 4), "\n")
cat("  theta (Gumbel)   =", signif(theta_gumbel,   4), "\n")

# --- Joint exceedance under fitted dependence (plug-in, closed form)
p_joint <- joint_exceed_gumbel(p_SVL_stokes, p_TL_stokes, theta_gumbel)

cat("\nJoint exceedance under fitted logistic/Gumbel dependence (plug-in):\n")
cat("  P(SVL >", SVL_stokes, ", TL >", TL_stokes, ") ≈", signif(p_joint, 4), "\n", sep = "")


# ============================================================
# 11. Joint bootstrap distribution for (SVL*, TL*)
# ============================================================

# Keep only successful, finite bootstrap draws with xi < 0
boot_ok <- boot_conv & is.finite(boot_xi) & (boot_xi < 0)

if (!any(boot_ok)) {
  stop("No valid bootstrap draws for endpoints (boot_ok is empty).")
}

# Transform log-endpoints to original (cm) scale
svl_boot_cm <- exp(boot_y1star[boot_ok])
tl_boot_cm  <- exp(boot_y2star[boot_ok])

boot_joint <- data.frame(
  SVL = svl_boot_cm,
  TL  = tl_boot_cm
)

q_upper <- ci_level

ub_svl_marg <- as.numeric(quantile(svl_boot_cm, probs = q_upper, na.rm = TRUE))
ub_tl_marg  <- as.numeric(quantile(tl_boot_cm,  probs = q_upper, na.rm = TRUE))

# Reasonable plotting limits
svl_xlim <- quantile(svl_boot_cm, c(0.01, 0.99), na.rm = TRUE)
tl_ylim  <- quantile(tl_boot_cm,  c(0.01, 0.99), na.rm = TRUE)

# 2D KDE on a grid
kde_joint <- MASS::kde2d(
  x = boot_joint$SVL,
  y = boot_joint$TL,
  n = 200,
  lims = c(svl_xlim[1], svl_xlim[2], tl_ylim[1], tl_ylim[2])
)

# helper to evaluate kde2d at arbitrary points (bilinear interpolation)
eval_kde2d <- function(kde, x, y) {
  nx <- length(kde$x)
  ny <- length(kde$y)
  
  ix <- findInterval(x, kde$x, all.inside = TRUE)
  iy <- findInterval(y, kde$y, all.inside = TRUE)
  
  ix[ix >= nx] <- nx - 1
  iy[iy >= ny] <- ny - 1
  
  x1 <- kde$x[ix]
  x2 <- kde$x[ix + 1]
  y1 <- kde$y[iy]
  y2 <- kde$y[iy + 1]
  
  z11 <- kde$z[cbind(ix,     iy)]
  z21 <- kde$z[cbind(ix + 1, iy)]
  z12 <- kde$z[cbind(ix,     iy + 1)]
  z22 <- kde$z[cbind(ix + 1, iy + 1)]
  
  wx <- (x - x1) / (x2 - x1 + 1e-12)
  wy <- (y - y1) / (y2 - y1 + 1e-12)
  
  (1 - wx) * (1 - wy) * z11 +
    wx      * (1 - wy) * z21 +
    (1 - wx) * wy      * z12 +
    wx      * wy       * z22
}



# Data frame for KDE grid (for contour)
kde_df <- with(
  kde_joint,
  expand.grid(SVL = x, TL = y) |>
    transform(z = as.vector(z))
)

# HPD contour threshold from KDE grid
hpd_threshold_from_kde <- function(kde, level = 0.90) {
  z <- as.vector(kde$z)
  dx <- diff(kde$x)[1]
  dy <- diff(kde$y)[1]
  w  <- z * dx * dy
  ord <- order(z, decreasing = TRUE)
  cum <- cumsum(w[ord])
  z[ord][which(cum >= level)[1]]
}

hpd_threshold <- hpd_threshold_from_kde(kde_joint, level = ci_level)

p_joint <- ggplot() +
  # 2D KDE (filled levels)
  geom_density_2d_filled(
    data  = boot_joint,
    aes(x = SVL, y = TL),
    alpha = 0.7
  ) +
  # Bootstrap draws (faint)
  geom_point(
    data  = boot_joint,
    aes(x = SVL, y = TL),
    alpha = 0.25,
    size  = 0.6,
    color = "grey20"
  ) +
  # Point estimate (MAP / plug-in)
  geom_point(
    aes(x = svl_mode, y = tl_mode),
    color = "purple",
    size  = 3
  ) +
  # Stokes alligator (SVL, TL)
  geom_point(
    aes(x = SVL_stokes, y = TL_stokes),
    color = "red",
    size  = 3,
    shape = 4,
    stroke = 1.2
  ) +
  coord_cartesian(xlim = svl_xlim, ylim = tl_ylim) +
  labs(
    x     = "SVL endpoint (cm)",
    y     = "TL endpoint (cm)",
    fill  = "Density level"
  ) +
  theme_science

p_joint

ggsave(
  file.path(FIG_DIR, "SVL_TL_joint_endpoint_bootstrap.png"),
  p_joint,
  dpi   = 600,
  w     = 6.5,
  h     = 4.5,
  units = "in"
)

# ============================================================
# 12. Univariate bootstrap endpoint distributions: SVL* and TL*
#      MAP = mode of 1D KDE within same window
# ============================================================

# Keep only successful, finite bootstrap draws with xi < 0 (should already exist)
boot_ok <- boot_conv & is.finite(boot_xi) & (boot_xi < 0)
if (!any(boot_ok)) stop("No valid bootstrap draws for endpoints (boot_ok is empty).")

# Bootstrap endpoints on original (cm) scale
svl_boot_cm <- exp(boot_y1star[boot_ok])
tl_boot_cm  <- exp(boot_y2star[boot_ok])

# Use same plotting windows as for the 2D KDE
svl_xlim_uni <- svl_window
tl_xlim_uni  <- tl_window

svl_in <- svl_boot_cm[svl_boot_cm >= svl_xlim_uni[1] & svl_boot_cm <= svl_xlim_uni[2]]
tl_in  <- tl_boot_cm[ tl_boot_cm  >= tl_xlim_uni[1]  & tl_boot_cm  <= tl_xlim_uni[2]]

if (length(svl_in) < 30L || length(tl_in) < 30L)
  warning("Few draws inside univariate KDE window; consider widening svl_window/tl_window.")

# 1D KDEs within the plotting window
svl_dens <- density(
  svl_in,
  n    = 512,
  from = svl_xlim_uni[1],
  to   = svl_xlim_uni[2]
)

tl_dens <- density(
  tl_in,
  n    = 512,
  from = tl_xlim_uni[1],
  to   = tl_xlim_uni[2]
)

svl_df      <- data.frame(SVL = svl_dens$x, density = svl_dens$y)
svl_samp_df <- data.frame(SVL = svl_in)

tl_df       <- data.frame(TL = tl_dens$x, density = tl_dens$y)
tl_samp_df  <- data.frame(TL = tl_in)

# Marginal MAPs = modes of 1D KDEs
map_svl_marg <- svl_dens$x[which.max(svl_dens$y)]
map_tl_marg  <- tl_dens$x[which.max(tl_dens$y)]

cat("\nMarginal KDE MAPs for endpoints (within window):\n")
cat("  SVL*_MAP, marginal ≈", round(map_svl_marg, 1), "cm\n")
cat("  TL*_MAP,  marginal ≈", round(map_tl_marg,  1), "cm\n")

# If you didn’t compute marginal one-sided upper bounds yet, uncomment:
# ub_svl_marg <- as.numeric(quantile(svl_boot_cm, probs = ci_level, na.rm = TRUE))
# ub_tl_marg  <- as.numeric(quantile(tl_boot_cm,  probs = ci_level, na.rm = TRUE))

# ---- SVL* bootstrap distribution ----
p_svl <- ggplot() +
  geom_histogram(
    data  = svl_samp_df,
    aes(x = SVL, y = ..density..),
    bins  = 40,
    fill  = "lightblue",
    color = "black",
    alpha = 0.6
  ) +
  geom_line(
    data = svl_df,
    aes(x = SVL, y = density),
    color = "darkblue",
    linewidth = 1.1
  ) +
  geom_vline(
    xintercept = map_svl_marg,
    color = "purple",
    linetype = "dashed",
    linewidth = 1.1
  ) +
  geom_vline(
    xintercept = ub_svl_marg,
    color = "orange",
    linetype = "dotdash",
    linewidth = 1.1
  ) +
  # Optional: mark Stokes on SVL axis
  geom_vline(
    xintercept = SVL_stokes,
    color = "red",
    linetype = "solid",
    linewidth = 0.9,
    alpha = 0.8
  ) +
  labs(x = "SVL endpoint (cm)", y = "Density") +
  coord_cartesian(xlim = svl_xlim_uni) +
  theme_science

# ---- TL* bootstrap distribution ----
p_tl <- ggplot() +
  geom_histogram(
    data  = tl_samp_df,
    aes(x = TL, y = ..density..),
    bins  = 40,
    fill  = "lightgreen",
    color = "black",
    alpha = 0.6
  ) +
  geom_line(
    data = tl_df,
    aes(x = TL, y = density),
    color = "darkgreen",
    linewidth = 1.1
  ) +
  geom_vline(
    xintercept = map_tl_marg,
    color = "purple",
    linetype = "dashed",
    linewidth = 1.1
  ) +
  geom_vline(
    xintercept = ub_tl_marg,
    color = "orange",
    linetype = "dotdash",
    linewidth = 1.1
  ) +
  # Optional: mark Stokes on TL axis
  geom_vline(
    xintercept = TL_stokes,
    color = "red",
    linetype = "solid",
    linewidth = 0.9,
    alpha = 0.8
  ) +
  labs(x = "TL endpoint (cm)", y = "Density") +
  coord_cartesian(xlim = tl_xlim_uni) +
  theme_science

print(p_svl)
print(p_tl)

ggsave(file.path(FIG_DIR, "SVL_endpoint_bootstrap.png"),
       p_svl, dpi = 600, w = 6.5, h = 4.5, units = "in")
ggsave(file.path(FIG_DIR, "TL_endpoint_bootstrap.png"),
       p_tl,  dpi = 600, w = 6.5, h = 4.5, units = "in")



