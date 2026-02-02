#!/usr/bin/env Rscript

# ============================================================
# Bivariate EVT with censored likelihood & logistic dependence
# (SVL, TL) — American alligators
#
# Censored bivariate POT log-likelihood & profile in bvpot.R:
#   - prep_censored_pot()          : censored data pre-processing
#   - gpd_survival_endpoint()      : conditional tail S(y | Y>u) in (xi, y*)
#   - gpd_quantile_endpoint()      : tail quantile in (xi, y*)
#   - nll()                        : core negative log-likelihood (xi, y*)
#   - fit_bvpot()                  : MLE wrapper (free vs common shape)
#   - bvpot_endpoints()            : endpoints y1*, y2* from fit object
#   - profile_ll_joint()           : joint profile at fixed (y1*, y2*)
#   - profile_grid_joint()         : joint profile on a grid (y1*, y2*)
# ============================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(ggplot2)
  library(evd)        # fpot, pgpd, dgpd, qgpd
  library(scales)
  library(glue)
  library(forcats)
  library(grid)
  library(copula)     # Gumbel / logistic EV copula
})

set.seed(42)

# ------------------------------------------------------------
# 0. Load likelihood + fit + profile from bvpot.R
# ------------------------------------------------------------
source("../bvpot.R")

# ---------------------------
# Directories & theming
# ---------------------------
FIG_DIR <- "Figures_gators"
if (!dir.exists(FIG_DIR)) dir.create(FIG_DIR, recursive = TRUE)

theme_science_polished <- theme_minimal(base_family = "Arial", base_size = 12) +
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
ci_level      <- 0.90
u_scan_lo_q   <- 0.60
u_scan_hi_q   <- 0.99
u_scan_n      <- 50
min_ex_mrl    <- 5
min_ex_fit    <- 20

trait_names   <- c("SVL", "TL")
thresh_q_opt  <- c(SVL = 0.8, TL = 0.8)  # anchor quantiles on log scale

# KS bootstrap settings (univariate diagnostics)
B_KS   <- 1000L

# ---------------------------
# Plotting knobs for endpoint profiles
# ---------------------------
# Optional axis limits (in cm) for the joint profile contour
# and the 1D profiles. Keep NULL for automatic ggplot limits.
xlim_prof_joint <- c(220, 320) # e.g. c(220, 320) for SVL*
ylim_prof_joint <- c(380, 560)  # e.g. c(380, 560) for TL*

xlim_prof_y1 <- c(220, 320)     # e.g. c(220, 320) for SVL endpoint 1D profile
xlim_prof_y2 <- c(380, 560)     # e.g. c(380, 560) for TL endpoint 1D profile

# ============================================================
# 1. Ingest and basic cleaning
# ============================================================

DATA_XLSX <- "Data/experimental_alligator_harvest_woodward.xlsx"
df_raw    <- read_excel(DATA_XLSX)

stopifnot(all(c("SVL","TL","Deform") %in% names(df_raw)))

# Deform == 1 or 3: tail broken ⇒ TL structurally missing; SVL kept.
df <- df_raw %>%
  mutate(
    SVL = suppressWarnings(as.numeric(SVL)),
    TL  = suppressWarnings(as.numeric(TL)),
    TL  = ifelse(Deform %in% c(1, 3), NA_real_, TL)
  ) %>%
  transmute(
    specimen = row_number(),
    SVL, TL
  )

# ============================================================
# 2. Completeness diagnostics
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

p_compl <- compl_tbl %>%
  mutate(trait = fct_inorder(trait)) %>%
  ggplot(aes(trait, completeness)) +
  geom_col(fill = "#3B82F6") +
  geom_text(aes(label = percent(completeness, accuracy = 0.1)),
            vjust = -0.2, size = 4) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1.10)) +
  labs(title = "Completeness (SVL, TL)", x = NULL, y = "Completeness") +
  theme_science_polished

p_compl
ggsave(file.path(FIG_DIR, "completeness_SVL_TL.png"), p_compl,
       dpi = 600, w = 6.0, h = 4.2, units = "in")

# ============================================================
# 3. Log-transform and split patterns
# ============================================================

df_log <- df %>%
  mutate(
    log_SVL = ifelse(is.finite(SVL) & SVL > 0, log(SVL), NA_real_),
    log_TL  = ifelse(is.finite(TL)  & TL  > 0, log(TL),  NA_real_)
  )

df_complete <- df_log %>% filter(!is.na(log_SVL), !is.na(log_TL))
df_svl_only <- df_log %>% filter(!is.na(log_SVL),  is.na(log_TL))
df_tl_only  <- df_log %>% filter(is.na(log_SVL),  !is.na(log_TL))

cat("\nRow counts by observation pattern:\n")
cat("  complete (SVL & TL) :", nrow(df_complete), "\n")
cat("  SVL only            :", nrow(df_svl_only), "\n")
cat("  TL only             :", nrow(df_tl_only), "\n")

# ============================================================
# 4. Thresholds on log scale (per trait, using all non-missing)
# ============================================================

u0_by_trait <- setNames(numeric(length(trait_names)), trait_names)
for (tr in trait_names) {
  v_log <- df_log[[paste0("log_", tr)]]
  u0_by_trait[tr] <- as.numeric(quantile(v_log, thresh_q_opt[tr], na.rm = TRUE))
}
print(u0_by_trait)
print(exp(u0_by_trait))

u1    <- u0_by_trait["SVL"]
u2    <- u0_by_trait["TL"]
u_vec <- c(u1, u2)

# ============================================================
# 4B. Pairwise tail regimes & missing/observed/exceed/censored
# ============================================================

x1 <- df_log$log_SVL
x2 <- df_log$log_TL

O1 <- is.finite(x1); M1 <- !O1
O2 <- is.finite(x2); M2 <- !O2

E1 <- O1 & (x1 > u1)
E2 <- O2 & (x2 > u2)
B1 <- O1 & !E1
B2 <- O2 & !E2

pattern_df <- tibble(
  specimen = df$specimen,
  log_SVL  = x1,
  log_TL   = x2,
  O1, O2, M1, M2, E1, E2, B1, B2,
  pattern = case_when(
    E1 & E2           ~ "E1,E2   (both exceed)",
    E1 & B2           ~ "E1,B2   (SVL exceed, TL censored)",
    B1 & E2           ~ "B1,E2   (SVL censored, TL exceed)",
    E1 & M2           ~ "E1,M2   (SVL exceed, TL missing)",
    M1 & E2           ~ "M1,E2   (SVL missing, TL exceed)",
    B1 & B2           ~ "B1,B2   (both censored)",
    B1 & M2           ~ "B1,M2   (SVL censored, TL missing)",
    M1 & B2           ~ "M1,B2   (SVL missing, TL censored)",
    TRUE              ~ "M1,M2   (both missing)"
  )
)

make_pair_plot <- function(t1 = "SVL", t2 = "TL") {
  fn <- file.path(FIG_DIR, glue("pair_scatter_{t1}_{t2}_final.png"))
  
  xr <- range(x1[O1], na.rm = TRUE)
  yr <- range(x2[O2], na.rm = TRUE)
  pad <- function(r, p = 0.03) {
    w <- diff(r); c(r[1] - p * w, r[2] + p * w)
  }
  xr <- pad(xr); yr <- pad(yr)
  
  both_obs <- O1 & O2
  
  dd_dot <- tibble(
    specimen = df$specimen[both_obs],
    x        = x1[both_obs],
    y        = x2[both_obs],
    kind     = case_when(
      E1[both_obs] & E2[both_obs] ~ "(E1,E2)",
      E1[both_obs] & B2[both_obs] ~ "(E1,¬E2)",
      B1[both_obs] & E2[both_obs] ~ "(¬E1,E2)",
      TRUE                        ~ "(¬E1,¬E2)"
    )
  )
  
  dd_v <- tibble(
    specimen = df$specimen[E1 & M2],
    x   = x1[E1 & M2],
    y0  = yr[1],
    y1  = yr[2],
    kind = "(E1,¬O2)"
  )
  
  dd_h <- tibble(
    specimen = df$specimen[M1 & E2],
    y   = x2[M1 & E2],
    x0  = xr[1],
    x1  = xr[2],
    kind = "(¬O1,E2)"
  )
  
  col_map <- c(
    "(E1,E2)"   = "#DC2626",
    "(E1,¬E2)"  = "#2563EB",
    "(¬E1,E2)"  = "#059669",
    "(¬E1,¬E2)" = "grey70",
    "(E1,¬O2)"  = "#2563EB",
    "(¬O1,E2)"  = "#059669"
  )
  
  p <- ggplot() +
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
    scale_color_manual(values = col_map, name = "Regime") +
    coord_cartesian(xlim = xr, ylim = yr, expand = FALSE) +
    labs(
      x = glue("log({t1} [cm]) "),
      y = glue("log({t2} [cm])")
    ) +
    theme_science_polished
  
  print(p)
  ggsave(fn, p, dpi = 600, w = 6.8, h = 5.6, units = "in")
  message("Saved: ", normalizePath(fn))
}

make_pair_plot("SVL", "TL")

# ============================================================
# 5. Univariate diagnostics, KS bootstrap, Wald tests
# ============================================================

make_thresholds <- function(y, q_from = u_scan_lo_q, q_to = u_scan_hi_q, n = u_scan_n) {
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
    n_exceed = sum(y > u)
  )
}

adj_scale_fun <- function(scale, xi, u, u0) {
  scale - xi * (u - u0)
}

mrl_data <- function(y, u_seq, min_ex = 5) {
  sapply(u_seq, function(u) {
    ex <- y[y > u] - u
    if (length(ex) < min_ex) return(NA_real_)
    mean(ex)
  })
}

diagnostic_plots <- function(y, u, sigma_hat, xi_hat, label) {
  exc   <- y[y > u] - u
  n     <- length(exc)
  probs <- ppoints(n)
  
  if (abs(xi_hat) > 1e-10) {
    theo_q <- u + sigma_hat / xi_hat * (probs^(-xi_hat) - 1)
  } else {
    theo_q <- u - sigma_hat * log(probs)
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
    theme_science_polished
  
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
    theme_science_polished
  
  F_hat <- evd::pgpd(exc, loc = 0, scale = sigma_hat, shape = xi_hat)
  pks <- ggplot(data.frame(F_hat = F_hat), aes(F_hat)) +
    geom_histogram(aes(y = ..density..), bins = 20,
                   fill = "skyblue", color = "black", alpha = 0.7) +
    geom_hline(yintercept = 1, color = "red", linetype = "dashed") +
    labs(title = glue("Uniformity (PIT): {label}")) +
    theme_science_polished
  
  list(pqq = pqq, ppp = ppp, pks = pks)
}

ks_gpd_boot_from_fit <- function(y, u, scale_hat, shape_hat,
                                 B = 1000L, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  y_exceed <- y[y > u]
  ex       <- y_exceed - u
  n_ex     <- length(ex)
  if (n_ex < 20L) {
    warning("Not enough exceedances for KS bootstrap (n_ex < 20); returning NA.")
    return(list(D = NA_real_, p = NA_real_, D_boot = numeric(0L), n_exceed = n_ex))
  }
  
  Fn_ex      <- ecdf(ex)
  ex_sorted  <- sort(ex)
  Fgpd_obs   <- evd::pgpd(ex_sorted, loc = 0, scale = scale_hat, shape = shape_hat)
  Fn_vals    <- Fn_ex(ex_sorted)
  D_obs      <- max(abs(Fn_vals - Fgpd_obs))
  
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

run_trait_diag <- function(tr_key) {
  v_raw <- df[[tr_key]] |> suppressWarnings(as.numeric())
  y     <- log(v_raw[is.finite(v_raw) & v_raw > 0])
  ylab  <- glue("log {tr_key}")
  stopifnot(length(y) >= 40)
  
  q_anchor <- thresh_q_opt[tr_key]
  u_seq    <- make_thresholds(y)
  u0       <- quantile(y, q_anchor, na.rm = TRUE) |> as.numeric()
  
  # MRL
  mrl_vals <- mrl_data(y, u_seq, min_ex = min_ex_mrl)
  p_mrl <- ggplot(data.frame(u = u_seq, mrl = mrl_vals), aes(u, mrl)) +
    geom_line() + geom_point() +
    geom_vline(xintercept = u0, linetype = "dashed", color = "red") +
    labs(
      title = glue("MRL plot ({ylab})"),
      x     = glue("Threshold ({ylab})"),
      y     = "Mean excess"
    ) +
    theme_science_polished
  ggsave(file.path(FIG_DIR, glue("{tr_key}_log_mrl.png")),
         p_mrl, dpi = 600, w = 6.2, h = 4.2, units = "in")
  
  # Fit at anchor u0
  fit0  <- fit_gpd_at_u(y, u0)
  xi0   <- fit0$shape
  sc0   <- fit0$scale
  cov0  <- fit0$cov
  xi_se <- sqrt(cov0[2, 2])
  sigma_se <- sqrt(cov0[1, 1])
  
  z_wald <- xi0 / xi_se
  p_wald <- pnorm(z_wald)  # one-sided H1: xi < 0
  
  # KS bootstrap
  seed_ks <- switch(tr_key,
                    "SVL" = 101L,
                    "TL"  = 202L,
                    303L)
  ks_res <- ks_gpd_boot_from_fit(
    y          = y,
    u          = u0,
    scale_hat  = sc0,
    shape_hat  = xi0,
    B          = B_KS,
    seed       = seed_ks
  )
  
  # Threshold scans
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
    if (length(ex) < min_ex_fit) next
    
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
      y = "Shape (xi)"
    ) +
    theme_science_polished
  
  p_adj <- ggplot(df_scan, aes(u, adj)) +
    geom_point(color = "blue") +
    geom_errorbar(aes(ymin = adj_lo, ymax = adj_hi),
                  width = 0.03, color = "blue") +
    geom_vline(xintercept = u0, color = "red", linetype = "dashed") +
    geom_hline(yintercept = sc0, color = "red", linetype = "dashed") +
    labs(
      x = glue("Threshold ({ylab})"),
      y = "Adjusted scale"
    ) +
    theme_science_polished
  
  ggsave(file.path(FIG_DIR, glue("{tr_key}_log_shape_stability.png")),
         p_shape, dpi = 600, w = 6.2, h = 4.2, units = "in")
  ggsave(file.path(FIG_DIR, glue("{tr_key}_log_adj_scale_stability.png")),
         p_adj,   dpi = 600, w = 6.2, h = 4.2, units = "in")
  
  di <- diagnostic_plots(y, u0, sc0, xi0, ylab)
  ggsave(file.path(FIG_DIR, glue("{tr_key}_log_qq.png")),
         di$pqq, dpi = 600, w = 6.2, h = 4.2, units = "in")
  ggsave(file.path(FIG_DIR, glue("{tr_key}_log_pp.png")),
         di$ppp, dpi = 600, w = 6.2, h = 4.2, units = "in")
  ggsave(file.path(FIG_DIR, glue("{tr_key}_log_pit_uniformity.png")),
         di$pks, dpi = 600, w = 6.2, h = 4.2, units = "in")
  
  cat(glue("\n{tr_key}: Wald test xi<0 at u0={round(u0,3)}: ",
           "xî={round(xi0,3)}, se={round(xi_se,3)}, ",
           "z={round(z_wald,3)}, p={signif(p_wald,3)}"))
  cat(glue("\n{tr_key}: KS bootstrap: D={round(ks_res$D,3)}, ",
           "p_boot={signif(ks_res$p,3)}, n_exceed={ks_res$n_exceed}\n"))
  
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

diag_SVL <- run_trait_diag("SVL")
diag_TL  <- run_trait_diag("TL")
marg_diag_tbl <- bind_rows(diag_SVL, diag_TL)
print(marg_diag_tbl)

# ============================================================
# 6. Univariate POT fits using all non-missing data per margin
# ============================================================

y_SVL_all <- df_log$log_SVL[!is.na(df_log$log_SVL)]
y_TL_all  <- df_log$log_TL[!is.na(df_log$log_TL)]

fit_SVL_uni <- evd::fpot(y_SVL_all, threshold = u1, model = "gpd")
fit_TL_uni  <- evd::fpot(y_TL_all,  threshold = u2, model = "gpd")

par_SVL_uni <- fit_SVL_uni$estimate
par_TL_uni  <- fit_TL_uni$estimate

sigma1_hat_uni <- unname(par_SVL_uni["scale"])
xi1_hat_uni    <- unname(par_SVL_uni["shape"])
sigma2_hat_uni <- unname(par_TL_uni["scale"])
xi2_hat_uni    <- unname(par_TL_uni["shape"])

cat("\nUnivariate GP fits (log scale, evd::fpot, ALL marginal data):\n")
cat("  SVL: scale =", sigma1_hat_uni, " shape =", xi1_hat_uni, "\n")
cat("  TL : scale =", sigma2_hat_uni, " shape =", xi2_hat_uni, "\n")

# ============================================================
# 6b. Near-endpoint empirical unconditional survival curves
# ============================================================

make_endpoint_sf_curve_uncond <- function(y, u, name) {
  y_all <- y[is.finite(y)]
  n_all <- length(y_all)
  
  y_tail <- sort(y_all[y_all > u])
  n_tail <- length(y_tail)
  if (n_tail < 5L) {
    warning(glue("Too few exceedances for {name} (n_tail < 5)."))
    return(tibble())
  }
  
  y_max <- max(y_tail)
  k      <- seq_len(n_tail - 1L)
  t_hat  <- y_max - y_tail[k]
  
  S_cond   <- (n_tail - k) / n_tail
  p_tail   <- n_tail / n_all
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

trait_cols <- c("SVL" = "#377eb8", "TL" = "#1b9e77")

p_sf_endpoint <- ggplot(sf_all, aes(x = t_hat, y = S_hat, colour = Trait)) +
  geom_point(size = 2, alpha = 0.9) +
  scale_colour_manual(values = trait_cols, name = "Trait") +
  scale_x_log10() +
  scale_y_log10() +
  labs(
    x = expression(hat(t) == y[max] - y),
    y = expression(hat(S)(t) == P(Y > y[max] - t))
  ) +
  theme_science_polished

p_sf_endpoint

ggsave(
  file.path(FIG_DIR, "SVL_TL_empirical_survival_endpoint_uncond.png"),
  p_sf_endpoint,
  dpi   = 600,
  w     = 6.0,
  h     = 4.5,
  units = "in"
)

# ============================================================
# 7. Bivariate POT with censored likelihood via fit_bvpot()
# ============================================================

X_all <- as.matrix(df_log[, c("log_SVL", "log_TL")])

cat("\nBivariate sample size (complete cases SVL & TL):", nrow(df_complete), "\n")

quad_tab <- table(
  SVL_exceed = df_complete$log_SVL > u1,
  TL_exceed  = df_complete$log_TL > u2
)
cat("\nQuadrant table (on log scale, complete cases only):\n")
print(quad_tab)

# ------------------------------------------------------------
# 8. Maximum likelihood via fit_bvpot(): free vs common shape
# ------------------------------------------------------------

# 8.1 Free shapes (H1: xi1 != xi2)
fit_free <- fit_bvpot(
  X            = X_all,
  u            = u_vec,
  common_shape = FALSE,
  method       = "Nelder-Mead"
)

par_free_ystar <- fit_free$par_ystar      # (xi1, y1*, xi2, y2*, dep)
ell_hat_free   <- fit_free$loglik

cat("\nBivariate MLEs (full likelihood via fit_bvpot, free shapes — (xi,y*)): \n")
print(par_free_ystar)
cat("Maximized log-likelihood (free) =", ell_hat_free, "\n")

# 8.2 Common tail shape (H0: xi1 = xi2 = xi < 0)
fit_eq <- fit_bvpot(
  X            = X_all,
  u            = u_vec,
  common_shape = TRUE,
  method       = "Nelder-Mead"
)

par_eq_ystar  <- fit_eq$par_ystar        # (xi, y1*, y2*, dep)
ell_hat_eq    <- fit_eq$loglik

cat("\nCommon tail shape MLEs (via fit_bvpot, H0: xi1 = xi2, parametrisation (xi, y1*, y2*, dep)):\n")
print(par_eq_ystar)
cat("Maximized log-likelihood (H0) =", ell_hat_eq, "\n")

# LR test H0: xi1 = xi2 vs H1: xi1 != xi2
LR_shape <- 2 * (ell_hat_free - ell_hat_eq)
df_LR    <- 1
p_LR     <- 1 - pchisq(LR_shape, df = df_LR)

cat("\nLikelihood ratio test for tail equality H0: xi1 = xi2\n")
cat("  LR statistic =", LR_shape, " with df =", df_LR, "\n")
cat("  p-value      =", p_LR, "\n")

# ============================================================
# 9. Endpoint estimation under common tail + joint profile LL
# ============================================================

# Common ξ and logistic dep from endpoint parametrisation
xi_eq_hat      <- unname(par_eq_ystar["xi"])
alpha_logistic <- unname(par_eq_ystar["dep"])

# Endpoints (log-scale) from fit_eq$endpoints (Weibull domain)
ystar_hat_log  <- fit_eq$endpoints
ystar1_hat_log <- ystar_hat_log["y1_star"]
ystar2_hat_log <- ystar_hat_log["y2_star"]

cat("\nLog-scale endpoint estimates (common-tail GP, Weibull domain, (xi,y*)):\n")
cat("  xî (common) =", xi_eq_hat, "\n")
cat("  y*_SVL (log) =", ystar1_hat_log, "\n")
cat("  y*_TL  (log) =", ystar2_hat_log, "\n")

ystar1_hat_orig <- exp(ystar1_hat_log)
ystar2_hat_orig <- exp(ystar2_hat_log)

cat("\nPoint estimates for endpoints (original scale):\n")
cat("  SVL* ≈", ystar1_hat_orig, "cm\n")
cat("  TL*  ≈", ystar2_hat_orig, "cm\n")

# Tail fractions (used later for Stokes probabilities / MC)
n_SVL_total <- sum(!is.na(df_log$log_SVL))
n_TL_total  <- sum(!is.na(df_log$log_TL))

n_SVL_exceed_u1 <- sum(df_log$log_SVL > u1, na.rm = TRUE)
n_TL_exceed_u2  <- sum(df_log$log_TL  > u2, na.rm = TRUE)

tail_frac_SVL <- n_SVL_exceed_u1 / n_SVL_total
tail_frac_TL  <- n_TL_exceed_u2  / n_TL_total

cat("\nEmpirical tail fractions:\n")
cat("  P(SVL > u1) ≈", tail_frac_SVL, "(u1 =", u1, ", exp(u1) ≈", exp(u1), "cm)\n")
cat("  P(TL  > u2) ≈", tail_frac_TL,  "(u2 =", u2, ", exp(u2) ≈", exp(u2), "cm)\n")

# Sub-threshold samples for later Monte Carlo (Stokes section)
sub_SVL_cc   <- df_complete$log_SVL[df_complete$log_SVL <= u1]
sub_TL_cc    <- df_complete$log_TL[df_complete$log_TL <= u2]

# -------- Joint profile likelihood for (y1*, y2*) on log scale -----

y1_max_log <- max(df_log$log_SVL, na.rm = TRUE)
y2_max_log <- max(df_log$log_TL, na.rm = TRUE)

# Window around MLE endpoints (on log-scale)
span1 <- 1
span2 <- 1

grid_y1 <- seq(y1_max_log - span1, ystar1_hat_log + span1, length.out = 60)
grid_y2 <- seq(y2_max_log - span2, ystar2_hat_log + span2, length.out = 60)

cat("\nEvaluating joint profile likelihood on a", length(grid_y1),
    "x", length(grid_y2), "grid...\n")

prof_grid <- profile_grid_joint(
  fit        = fit_eq,
  X          = X_all,
  u          = u_vec,
  grid_y1    = grid_y1,
  grid_y2    = grid_y2,
  method     = "Nelder-Mead",
  control    = list(maxit = 200),
  penalty    = 1e6,
  verbose    = TRUE
)

prof_grid <- prof_grid %>%
  filter(is.finite(loglik)) %>%
  mutate(
    loglik_rel   = loglik - max(loglik, na.rm = TRUE),
    SVL_star_cm  = exp(y1_star),
    TL_star_cm   = exp(y2_star)
  )

# 90% joint confidence region for (y1*, y2*) via χ^2_2
lvl_joint_90 <- -0.5 * qchisq(ci_level, df = 2)

region_joint_90 <- prof_grid %>%
  filter(loglik_rel >= lvl_joint_90)

ub_y1_star_log <- max(region_joint_90$y1_star, na.rm = TRUE)
ub_y2_star_log <- max(region_joint_90$y2_star, na.rm = TRUE)

ub_SVL_star_cm <- exp(ub_y1_star_log)
ub_TL_star_cm  <- exp(ub_y2_star_log)

cat("\n", 100 * ci_level,
    "% joint profile upper bounds for endpoints (original scale):\n", sep = "")
cat("  SVL* ≤ ", ub_SVL_star_cm, " cm\n", sep = "")
cat("  TL*  ≤ ", ub_TL_star_cm,  " cm\n", sep = "")

# Joint profile likelihood contour plot (on original scale)
p_prof_joint <- ggplot(
  prof_grid,
  aes(x = SVL_star_cm, y = TL_star_cm, z = loglik_rel)
) +
  geom_contour_filled(bins = 20) +
  geom_contour(breaks = lvl_joint_90, colour = "black", linewidth = 1.0) +
  geom_point(aes(x = ystar1_hat_orig, y = ystar2_hat_orig),
             colour = "white", size = 3) +
  geom_vline(xintercept = ub_SVL_star_cm,
             colour = "orange", linetype = "dashed") +
  geom_hline(yintercept = ub_TL_star_cm,
             colour = "orange", linetype = "dashed") +
  labs(
    x    = "SVL endpoint (cm)",
    y    = "TL endpoint (cm)",
    fill = "Rel. loglik"
  ) +
  theme_science_polished


# Optional axis limits for joint profile plot
if (!is.null(xlim_prof_joint) || !is.null(ylim_prof_joint)) {
  p_prof_joint <- p_prof_joint +
    coord_cartesian(
      xlim   = xlim_prof_joint,
      ylim   = ylim_prof_joint,
      expand = FALSE
    )
}

p_prof_joint

ggsave(
  file.path(FIG_DIR, "SVL_TL_joint_endpoint_profile.png"),
  p_prof_joint,
  dpi   = 600,
  w     = 6.5,
  h     = 4.5,
  units = "in"
)

# (Optional) marginal 1D profile-like upper bounds via projection
#   using χ^2_1(0.90) threshold for 1D profiles:
lvl_1d_90 <- -0.5 * qchisq(ci_level, df = 1)

prof_y1 <- prof_grid %>%
  group_by(y1_star) %>%
  summarise(loglik = max(loglik, na.rm = TRUE), .groups = "drop") %>%
  mutate(loglik_rel = loglik - max(loglik, na.rm = TRUE))

prof_y2 <- prof_grid %>%
  group_by(y2_star) %>%
  summarise(loglik = max(loglik, na.rm = TRUE), .groups = "drop") %>%
  mutate(loglik_rel = loglik - max(loglik, na.rm = TRUE))

ub_y1_star_log_1d <- max(prof_y1$y1_star[prof_y1$loglik_rel >= lvl_1d_90], na.rm = TRUE)
ub_y2_star_log_1d <- max(prof_y2$y2_star[prof_y2$loglik_rel >= lvl_1d_90], na.rm = TRUE)

ub_SVL_star_cm_1d <- exp(ub_y1_star_log_1d)
ub_TL_star_cm_1d  <- exp(ub_y2_star_log_1d)

cat("\nProfile-based 1D upper bounds (projected, χ^2_1) for endpoints:\n")
cat("  SVL* ≤ ", ub_SVL_star_cm_1d, " cm\n", sep = "")
cat("  TL*  ≤ ", ub_TL_star_cm_1d,  " cm\n", sep = "")

# ---- 1D profile likelihood plots (original scale) ------------

# SVL* profile
prof_y1_plot <- prof_y1 %>%
  mutate(SVL_star_cm = exp(y1_star))

p_prof_y1 <- ggplot(prof_y1_plot,
                    aes(x = SVL_star_cm, y = loglik_rel)) +
  geom_line(linewidth = 0.7) +
  geom_hline(yintercept = 0, linetype = "solid", colour = "black") +
  geom_hline(yintercept = lvl_1d_90, linetype = "dashed", colour = "red") +
  geom_vline(xintercept = ystar1_hat_orig,
             linetype = "dotted", colour = "blue") +
  geom_vline(xintercept = ub_SVL_star_cm_1d,
             linetype = "dashed", colour = "orange") +
  labs(
    x = "SVL endpoint (cm)",
    y = "Relative log-likelihood"
  ) +
  theme_science_polished

if (!is.null(xlim_prof_y1)) {
  p_prof_y1 <- p_prof_y1 +
    coord_cartesian(xlim = xlim_prof_y1, expand = FALSE)
}

p_prof_y1
ggsave(
  file.path(FIG_DIR, "SVL_endpoint_profile_1D.png"),
  p_prof_y1,
  dpi   = 600,
  w     = 6.2,
  h     = 4.2,
  units = "in"
)

# TL* profile
prof_y2_plot <- prof_y2 %>%
  mutate(TL_star_cm = exp(y2_star))

p_prof_y2 <- ggplot(prof_y2_plot,
                    aes(x = TL_star_cm, y = loglik_rel)) +
  geom_line(linewidth = 0.7) +
  geom_hline(yintercept = 0, linetype = "solid", colour = "black") +
  geom_hline(yintercept = lvl_1d_90, linetype = "dashed", colour = "red") +
  geom_vline(xintercept = ystar2_hat_orig,
             linetype = "dotted", colour = "blue") +
  geom_vline(xintercept = ub_TL_star_cm_1d,
             linetype = "dashed", colour = "orange") +
  labs(
    x = "TL endpoint (cm)",
    y = "Relative log-likelihood"
  ) +
  theme_science_polished

if (!is.null(xlim_prof_y2)) {
  p_prof_y2 <- p_prof_y2 +
    coord_cartesian(xlim = xlim_prof_y2, expand = FALSE)
}

p_prof_y2
ggsave(
  file.path(FIG_DIR, "TL_endpoint_profile_1D.png"),
  p_prof_y2,
  dpi   = 600,
  w     = 6.2,
  h     = 4.2,
  units = "in"
)

# ============================================================
# 10. Exceedance probability for the Stokes alligator
#     using the endpoint parametrisation S(y) in ξ,y*
# ============================================================

# Unconditional survival in Weibull domain:
#   S(y) = P(Y>y) = P(Y>u) * ((y* - y)/(y* - u))^(-1/xi),  u < y < y*
gpd_survival_endpoint_uncond <- function(y, u, xi, y_star, tail_frac) {
  if (y <= u) return(1)  # below threshold, definitely not "extreme"
  if (xi >= 0) stop("gpd_survival_endpoint_uncond: xi must be < 0.")
  if (y_star <= u) stop("gpd_survival_endpoint_uncond: y_star must be > u.")
  if (y >= y_star) return(0)  # beyond model endpoint
  
  S_cond <- gpd_survival_endpoint(y, u, xi, y_star)  # from bvpot.R
  tail_frac * S_cond
}

SVL_stokes <- 236   # cm
TL_stokes  <- 450   # cm

log_SVL_stokes <- log(SVL_stokes)
log_TL_stokes  <- log(TL_stokes)

cat("\nStokes alligator (log scale):\n")
cat("  log(SVL) =", log_SVL_stokes, "\n")
cat("  log(TL)  =", log_TL_stokes,  "\n")

# Univariate predictive exceedance probabilities using ξ,y* parametrisation
p_SVL_stokes <- gpd_survival_endpoint_uncond(
  y         = log_SVL_stokes,
  u         = u1,
  xi        = xi_eq_hat,
  y_star    = ystar1_hat_log,
  tail_frac = tail_frac_SVL
)

p_TL_stokes <- gpd_survival_endpoint_uncond(
  y         = log_TL_stokes,
  u         = u2,
  xi        = xi_eq_hat,
  y_star    = ystar2_hat_log,
  tail_frac = tail_frac_TL
)

cat("\nUnivariate predictive exceedance probabilities for Stokes (Weibull endpoint model):\n")
cat("  P(SVL >", SVL_stokes, "cm) ≈", p_SVL_stokes,
    " → return period ≈", 1 / p_SVL_stokes, "individuals\n")
cat("  P(TL  >", TL_stokes, "cm) ≈", p_TL_stokes,
    " → return period ≈", 1 / p_TL_stokes,  "individuals\n")

# Empirical P(E | y) where E = {SVL>u1 or TL>u2}
has_SVL <- !is.na(df_log$log_SVL)
has_TL  <- !is.na(df_log$log_TL)

n_ind <- sum(has_SVL | has_TL)

E_emp <- (has_SVL & df_log$log_SVL > u1) |
  (has_TL  & df_log$log_TL  > u2)

n_E <- sum(E_emp, na.rm = TRUE)
tail_frac_E <- n_E / n_ind

cat("\nEmpirical tail fraction for at least one exceedance E:\n")
cat("  P(E | y) ≈", tail_frac_E,
    " (n_E =", n_E, "out of n =", n_ind, ")\n")

# Independence proxy
p_joint_indep <- p_SVL_stokes * p_TL_stokes

cat("\nJoint exceedance under independence approximation:\n")
cat("  P(SVL >", SVL_stokes, "cm, TL >", TL_stokes, "cm) ≈", p_joint_indep,
    " → return period ≈", 1 / p_joint_indep, "individuals\n")

# Logistic tail dependence via Gumbel copula
theta_gumbel   <- 1 / alpha_logistic
cat("\nLogistic EV copula parameter (from common-shape fit):\n")
cat("  alpha (logistic dep) =", alpha_logistic, "\n")
cat("  theta (Gumbel copula) =", theta_gumbel, "\n")

gumbel_cop <- gumbelCopula(param = theta_gumbel, dim = 2)

set.seed(2025)
N_mc <- 2e5L

U_mc  <- rCopula(N_mc, gumbel_cop)
Y1_mc <- numeric(N_mc)
Y2_mc <- numeric(N_mc)

for (k in seq_len(N_mc)) {
  u1k <- U_mc[k, 1]
  u2k <- U_mc[k, 2]
  
  # SVL (log-scale)
  if (u1k <= 1 - tail_frac_SVL || length(sub_SVL_cc) == 0L) {
    Y1_mc[k] <- if (length(sub_SVL_cc)) sample(sub_SVL_cc, 1L, replace = TRUE) else u1
  } else {
    v1   <- (u1k - (1 - tail_frac_SVL)) / tail_frac_SVL  # ∈ (0,1)
    Y1_mc[k] <- gpd_quantile_endpoint(
      v      = v1,
      u      = u1,
      xi     = xi_eq_hat,
      y_star = ystar1_hat_log
    )
  }
  
  # TL (log-scale)
  if (u2k <= 1 - tail_frac_TL || length(sub_TL_cc) == 0L) {
    Y2_mc[k] <- if (length(sub_TL_cc)) sample(sub_TL_cc, 1L, replace = TRUE) else u2
  } else {
    v2   <- (u2k - (1 - tail_frac_TL)) / tail_frac_TL  # ∈ (0,1)
    Y2_mc[k] <- gpd_quantile_endpoint(
      v      = v2,
      u      = u2,
      xi     = xi_eq_hat,
      y_star = ystar2_hat_log
    )
  }
}

E_sim <- (Y1_mc > u1) | (Y2_mc > u2)
A_sim <- (Y1_mc > log_SVL_stokes) & (Y2_mc > log_TL_stokes)

n_E_sim <- sum(E_sim)
if (n_E_sim == 0L) {
  stop("No simulated tail events E_sim; increase N_mc or check model parameters.")
}

p_A_given_E_log <- sum(A_sim & E_sim) / n_E_sim
p_joint_log     <- tail_frac_E * p_A_given_E_log

se_A_given_E <- sqrt(p_A_given_E_log * (1 - p_A_given_E_log) / n_E_sim)
se_joint_log <- tail_frac_E * se_A_given_E

cat("\nLogistic EV tail model (P(A | y) = P(E | y) * P(A | E, y)):\n")
cat("  P(A | E, y) ≈", p_A_given_E_log,
    " (MC s.e. ≈", se_A_given_E, ")\n")
cat("  P(SVL >", SVL_stokes, "cm, TL >", TL_stokes, "cm | y)\n")
cat("    ≈", p_joint_log,
    " (MC s.e. ≈", se_joint_log, ")\n")
cat("    → return period ≈", 1 / p_joint_log, "individuals\n")

p_E_sim <- n_E_sim / N_mc
cat("\nSimulated P(E) (from logistic model) ≈", p_E_sim,
    " vs empirical P(E | y) ≈", tail_frac_E, "\n")

N_year <- 7600  # approx mean annual harvest 2009–2024

RP_years_indep  <- (1 / p_joint_indep) / N_year
RP_years_logist <- (1 / p_joint_log)   / N_year

cat("\nReturn period in years (assuming ≈", N_year, "harvested per year):\n")
cat("  Independence: ≈", RP_years_indep,  "years\n")
cat("  Logistic tail model: ≈", RP_years_logist, "years\n")


# test !!!!!!!

## ------------------------------------------------------------
## 1. Profile of common-shape xi
## ------------------------------------------------------------

# (Re)define, just to be safe
profile_xi_common <- function(xi_val, fit_eq, X, u) {
  # Fix xi and optimise (y1*, y2*, dep)
  nll_wrap <- function(p_rest) {
    y1_star <- p_rest[1]
    y2_star <- p_rest[2]
    dep     <- p_rest[3]
    p <- c(xi = xi_val, y1_star = y1_star, y2_star = y2_star, dep = dep)
    nll(p, X = X, u = u, common_shape = TRUE)
  }
  start_rest <- c(fit_eq$par_ystar["y1_star"],
                  fit_eq$par_ystar["y2_star"],
                  fit_eq$par_ystar["dep"])
  opt <- optim(start_rest, nll_wrap, method = "Nelder-Mead")
  opt$value   # this is the profiled *nll* at xi = xi_val
}

## ------------------------------------------------------------
## 2. Evaluate on a grid of xi values
## ------------------------------------------------------------

# choose a grid around where you expect xi to live
xi_grid <- seq(-0.6, -0.05, length.out = 40)

prof_tbl <- tibble::tibble(
  xi  = xi_grid,
  nll = vapply(
    xi_grid,
    FUN = function(xi) profile_xi_common(xi, fit_eq, X_all, u_vec),
    FUN.VALUE = 0.0
  )
) %>%
  dplyr::filter(is.finite(nll)) %>%
  dplyr::mutate(
    loglik     = -nll,
    loglik_rel = loglik - max(loglik)  # relative profile
  )

print(prof_tbl)

# xi that minimises the profiled nll (i.e. maximises loglik)
xi_prof_hat <- prof_tbl$xi[which.min(prof_tbl$nll)]
cat("\nProfiled common xi-hat ≈", xi_prof_hat, "\n")

# Compare with the MLE from fit_eq
xi_mle <- fit_eq$par_ystar["xi"]
cat("fit_eq$xi =", xi_mle, "\n")

## ------------------------------------------------------------
## 3. Plot the profile
## ------------------------------------------------------------

library(ggplot2)

ggplot(prof_tbl, aes(x = xi, y = loglik_rel)) +
  geom_line() +
  geom_vline(xintercept = xi_mle,
             linetype = "dashed", colour = "red") +
  # optional: add verticals for marginal / free-shape xi's if you have them
  # geom_vline(xintercept = xi1_hat_uni, linetype = "dotted", colour = "blue") +
  # geom_vline(xintercept = xi2_hat_uni, linetype = "dotted", colour = "green") +
  labs(
    x = expression(xi),
    y = "Relative log-likelihood",
    title = "Common-shape profile in xi"
  ) +
  theme_minimal()

