# ============================================================
# Bivariate EVT with censored likelihood & logistic dependence
# (SVL, TL) — American alligators
# Full log-likelihood: censored bivariate POT (bvpot.c:nllbvclog)
# + univariate GP contributions for missing partners (evd::dgpd)
# Parametric bootstrap for (sigma1, sigma2, xi, y1*, y2*)
# + univariate KS bootstrap + Wald tests for xi<0
# ============================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(ggplot2)
  library(evd)        # fpot, pgpd, dgpd, qgpd, sep.bvdata
  library(scales)
  library(glue)
  library(forcats)
  library(grid)
  library(MASS)       # for ginv fallback if Hessian singular
  library(copula)     # Gumbel / logistic EV copula
})

set.seed(42)

# ------------------------------------------------------------
# 0. Load compiled C code (bvpot.c → bvpot.so / bvpot.dll)
# ------------------------------------------------------------
# system("R CMD SHLIB bvpot.c")
if (!is.loaded("nllbvclog")) {
  ## Adjust path / extension as needed:
  ## dyn.load("path/to/bvpot.so") or "bvpot.dll" on Windows
  dyn.load("bvpot.dll")
}

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

trait_names   <- c("SVL","TL")
thresh_q_opt  <- c(SVL = 0.94, TL = 0.94)  # anchor quantiles (log scale) best 0.94!

# Bootstrap settings (bivariate endpoint & KS)
B_boot <- 1000L   # parametric bivariate bootstrap replicates
B_KS   <- 1000L   # parametric KS bootstrap per margin

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
# 3. Log-transform and split patterns (complete / SVL-only / TL-only)
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

# ============================================================
# 4. Thresholds on log scale (per trait, using *all* non-missing)
# ============================================================

u0_by_trait <- setNames(numeric(length(trait_names)), trait_names)
for (tr in trait_names) {
  v_log <- df_log[[paste0("log_", tr)]]
  u0_by_trait[tr] <- as.numeric(quantile(v_log, thresh_q_opt[tr], na.rm = TRUE))
}
print(u0_by_trait)
print(exp(u0_by_trait))  # thresholds on original scale

# ============================================================
# 4B Pairwise tail regimes & missing/observed/exceed/censored
#     (SVL, TL) — quadrant plot on log-scale
# ============================================================

# Use the same thresholds as elsewhere
u1 <- u0_by_trait["SVL"]
u2 <- u0_by_trait["TL"]

# Precompute log vectors
x1 <- df_log$log_SVL
x2 <- df_log$log_TL

# Observed / missing
O1 <- is.finite(x1); M1 <- !O1
O2 <- is.finite(x2); M2 <- !O2

# Exceeding / below threshold (censored) among observed
E1 <- O1 & (x1 > u1)
E2 <- O2 & (x2 > u2)
B1 <- O1 & !E1       # observed but below u1 ⇒ censored at u1
B2 <- O2 & !E2       # observed but below u2 ⇒ censored at u2

# Full pattern per specimen (useful for tables if needed later)
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
  
  # use outer log vectors
  xr <- range(x1[O1], na.rm = TRUE)
  yr <- range(x2[O2], na.rm = TRUE)
  pad <- function(r, p = 0.03) {
    w <- diff(r); c(r[1] - p * w, r[2] + p * w)
  }
  xr <- pad(xr); yr <- pad(yr)
  x_min <- xr[1]; x_max <- xr[2]
  y_min <- yr[1]; y_max <- yr[2]
  
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
    x  = x1[E1 & M2],
    y0 = y_min,
    y1 = y_max,
    kind = "(E1,¬O2)"
  )
  
  dd_h <- tibble(
    specimen = df$specimen[M1 & E2],
    y  = x2[M1 & E2],
    x0 = x_min,
    x1 = x_max,
    kind = "(¬O1,E2)"
  )
  
  col_map <- c(
    "(E1,E2)"   = "#DC2626",  # both exceed
    "(E1,¬E2)"  = "#2563EB",  # SVL exceed, TL censored
    "(¬E1,E2)"  = "#059669",  # SVL censored, TL exceed
    "(¬E1,¬E2)" = "grey70",   # both censored
    "(E1,¬O2)"  = "#2563EB",  # SVL exceed, TL missing
    "(¬O1,E2)"  = "#059669"   # TL exceed, SVL missing
  )
  
  p <- ggplot() +
    # three tail quadrants
    annotate(
      "rect",
      xmin = x_min, xmax = u1,
      ymin = u2,    ymax = y_max,
      fill = col_map["(¬E1,E2)"], alpha = 0.04
    ) +
    annotate(
      "rect",
      xmin = u1,    xmax = x_max,
      ymin = y_min, ymax = u2,
      fill = col_map["(E1,¬E2)"], alpha = 0.04
    ) +
    annotate(
      "rect",
      xmin = u1,    xmax = x_max,
      ymin = u2,    ymax = y_max,
      fill = col_map["(E1,E2)"],  alpha = 0.07
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
    scale_color_manual(values = col_map, name = "Regime") +
    coord_cartesian(xlim = xr, ylim = yr, expand = FALSE) +
    labs(
      x = glue("log({t1} [cm])"),
      y = glue("log({t2} [cm])")
    ) +
    theme_science_polished
  
  print(p)
  ggsave(fn, p, dpi = 600, w = 6.8, h = 5.6, units = "in")
  message("Saved: ", normalizePath(fn))
}

# actually create the plot
make_pair_plot("SVL", "TL")

# ============================================================
# 5. Helper functions for univariate diagnostics (evd only)
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
  
  # PIT using evd::pgpd (exc >= 0, loc = 0)
  F_hat <- evd::pgpd(exc, loc = 0, scale = sigma_hat, shape = xi_hat)
  pks <- ggplot(data.frame(F_hat = F_hat), aes(F_hat)) +
    geom_histogram(aes(y = ..density..), bins = 20,
                   fill = "skyblue", color = "black", alpha = 0.7) +
    geom_hline(yintercept = 1, color = "red", linetype = "dashed") +
    labs(title = glue("Uniformity (PIT): {label}")) +
    theme_science_polished
  
  list(pqq = pqq, ppp = ppp, pks = pks)
}

# ------------------------------------------------------------
# KS goodness-of-fit via parametric bootstrap (univariate GPD)
# ------------------------------------------------------------
ks_gpd_boot_from_fit <- function(y, u, scale_hat, shape_hat,
                                 B = 1000L, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  # Exceedances over u
  y_exceed <- y[y > u]
  ex       <- y_exceed - u
  n_ex     <- length(ex)
  if (n_ex < 20L) {
    warning("Not enough exceedances for KS bootstrap (n_ex < 20); returning NA.")
    return(list(D = NA_real_, p = NA_real_, D_boot = numeric(0L), n_exceed = n_ex))
  }
  
  # Observed KS statistic under fitted GPD
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
# Trait-wise diagnostics, KS, and Wald p-values
# ------------------------------------------------------------
run_trait_diag <- function(tr_key) {
  v_raw <- df[[tr_key]] |> suppressWarnings(as.numeric())
  y     <- log(v_raw[is.finite(v_raw) & v_raw > 0])
  ylab  <- glue("log {tr_key}")
  stopifnot(length(y) >= 40)
  
  q_anchor <- thresh_q_opt[tr_key]
  u_seq    <- make_thresholds(y)
  u0       <- quantile(y, q_anchor, na.rm = TRUE) |> as.numeric()
  
  # Mean residual life
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
  
  # GPD fit at anchor threshold u0
  fit0  <- fit_gpd_at_u(y, u0)
  xi0   <- fit0$shape
  sc0   <- fit0$scale
  cov0  <- fit0$cov
  xi_se <- sqrt(cov0[2, 2])
  sigma_se <- sqrt(cov0[1, 1])
  
  # Wald test for xi < 0 (one-sided)
  z_wald <- xi0 / xi_se
  p_wald <- pnorm(z_wald)  # H1: xi < 0
  
  # KS goodness-of-fit via parametric bootstrap
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
      y = "Shape (xi)") +
    theme_science_polished
  
  p_adj <- ggplot(df_scan, aes(u, adj)) +
    geom_point(color = "blue") +
    geom_errorbar(aes(ymin = adj_lo, ymax = adj_hi),
                  width = 0.03, color = "blue") +
    geom_vline(xintercept = u0, color = "red", linetype = "dashed") +
    geom_hline(yintercept = sc0, color = "red", linetype = "dashed") +
    labs(
      x = glue("Threshold ({ylab})"),
      y = "Adjusted scale") +
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

diag_SVL <- run_trait_diag("SVL")
diag_TL  <- run_trait_diag("TL")
marg_diag_tbl <- bind_rows(diag_SVL, diag_TL)
print(marg_diag_tbl)

# ============================================================
# 6. Univariate POT fits using all non-missing data per margin
# ============================================================

y_SVL_all <- df_log$log_SVL[!is.na(df_log$log_SVL)]
y_TL_all  <- df_log$log_TL[!is.na(df_log$log_TL)]

u1 <- u0_by_trait["SVL"]
u2 <- u0_by_trait["TL"]

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
# 6b. Empirical near-endpoint *unconditional* survival curves
#      using empirical maxima as pseudo-endpoints
#      Ŝ(t) = P(Y > y_max - t)
# ============================================================

# Helper: build Ŝ(t) = P(Y > y_max - t) for one trait (on log-scale)
make_endpoint_sf_curve_uncond <- function(y, u, name) {
  # all log-data (non-missing)
  y_all <- y[is.finite(y)]
  n_all <- length(y_all)
  
  # tail data above threshold
  y_tail <- y_all[y_all > u]
  y_tail <- sort(y_tail)
  
  n_tail <- length(y_tail)
  if (n_tail < 5L) {
    warning(glue("Too few exceedances for {name} (n_tail < 5)."))
    return(tibble())
  }
  
  # pseudo-endpoint = empirical max within tail
  y_max <- max(y_tail)
  
  # indices 1..(n_tail-1) so Ŝ never equals 0
  k      <- seq_len(n_tail - 1L)
  t_hat  <- y_max - y_tail[k]          # distance to empirical max (log-scale)
  
  # conditional survival given Y > u:
  #   S_cond(t) = P(Y > y_tail[k] | Y > u) = (n_tail - k)/n_tail
  S_cond <- (n_tail - k) / n_tail
  
  # tail fraction P(Y > u) ≈ n_tail / n_all
  p_tail <- n_tail / n_all
  
  # unconditional survival:
  #   S_uncond(t) = P(Y > u) * S_cond(t) ≈ (n_tail / n_all) * (n_tail - k)/n_tail
  S_uncond <- p_tail * S_cond
  
  tibble(
    t_hat = t_hat,
    S_hat = S_uncond,
    Trait = name
  )
}

# Build curves for log(SVL) and log(TL)
sf_SVL <- make_endpoint_sf_curve_uncond(y_SVL_all, u1, "SVL")
sf_TL  <- make_endpoint_sf_curve_uncond(y_TL_all,  u2, "TL")

sf_all <- bind_rows(sf_SVL, sf_TL)

# Colours consistent with earlier plots
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
# 7. Bivariate POT with censored likelihood (complete cases)
#    + univariate contributions from missings (E1M2, M1E2)
# ============================================================

# 7.1 Complete-case matrix for censored bivariate logistic POT
X_biv <- df_complete %>%
  transmute(
    log_SVL = log_SVL,
    log_TL  = log_TL
  ) %>%
  as.matrix()

n_biv <- nrow(X_biv)
cat("\nBivariate sample size (complete cases SVL & TL):", n_biv, "\n")

u_vec <- c(u1, u2)

quad_tab <- table(
  SVL_exceed = X_biv[, 1] > u1,
  TL_exceed  = X_biv[, 2] > u2
)
cat("\nQuadrant table (on log scale, complete cases only):\n")
print(quad_tab)

# 7.2 Missing patterns that contribute to the tail likelihood

# Only cases with an exceedance and the other margin truly missing:
y1_E_M2 <- df_svl_only$log_SVL[df_svl_only$log_SVL > u1]  # SVL exceeds, TL missing
y2_M1_E <- df_tl_only$log_TL[df_tl_only$log_TL > u2]      # TL exceeds, SVL missing

cat("\nCounts of missing-partner exceedances:\n")
cat("  SVL exceed, TL missing (E1M2):", length(y1_E_M2), "\n")
cat("  TL exceed, SVL missing (M1E2):", length(y2_M1_E), "\n")

# ---------------------------
# 7.3 Censored logistic log-likelihood using nllbvclog from bvpot.c
# ---------------------------

make_ll_cens_logistic_bvpot <- function(x, u, cshape = FALSE, cscale = FALSE) {
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
    
    # C routine returns the *negative* log-likelihood in dns
    -out$dns
  }
}

ll_cens_log <- make_ll_cens_logistic_bvpot(
  x      = X_biv,
  u      = u_vec,
  cshape = FALSE,
  cscale = FALSE
)

# ---------------------------
# 7.4 Univariate GP contributions from missings (evd::dgpd)
# ---------------------------

lgp_exceed_evd <- function(y, u, scale, shape) {
  ex <- y - u
  evd::dgpd(ex, loc = 0, scale = scale, shape = shape, log = TRUE)
}

make_ll_missing_uni <- function(y1_E_M2, y2_M1_E, u1, u2) {
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

ll_missing <- make_ll_missing_uni(y1_E_M2, y2_M1_E, u1, u2)

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
# 8. Maximum likelihood: free shapes vs common shape
# ============================================================

# 8.1 Free shapes (H1)
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

# 8.2 Common tail shape (H0: xi1 = xi2 = xi)
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
#    + parametric bivariate bootstrap for CIs
# ============================================================

# 9.1 Point estimates of GP parameters and endpoints (log scale)

xi_eq_hat  <- unname(par_eq_phi["xi"])
sigma1_hat <- unname(par_eq_phi["scale1"])
sigma2_hat <- unname(par_eq_phi["scale2"])

ystar1_hat_log <- as.numeric(u1 - sigma1_hat / xi_eq_hat)  # SVL*
ystar2_hat_log <- as.numeric(u2 - sigma2_hat / xi_eq_hat)  # TL*

cat("\nLog-scale endpoint estimates (common-tail GP, H0):\n")
cat("  y*_SVL (log) =", ystar1_hat_log, "\n")
cat("  y*_TL  (log) =", ystar2_hat_log, "\n")

ystar1_hat_orig <- exp(ystar1_hat_log)
ystar2_hat_orig <- exp(ystar2_hat_log)

cat("\nPoint estimates for endpoints (original scale):\n")
cat("  SVL* ≈", ystar1_hat_orig, "\n")
cat("  TL*  ≈", ystar2_hat_orig, "\n")

# 9.2 Empirical tail fractions (used both for bootstrap & Stokes)

n_SVL_total <- sum(!is.na(df_log$log_SVL))
n_TL_total  <- sum(!is.na(df_log$log_TL))

n_SVL_exceed_u1 <- sum(df_log$log_SVL > u1, na.rm = TRUE)
n_TL_exceed_u2  <- sum(df_log$log_TL  > u2, na.rm = TRUE)

tail_frac_SVL <- n_SVL_exceed_u1 / n_SVL_total
tail_frac_TL  <- n_TL_exceed_u2  / n_TL_total

cat("\nEmpirical tail fractions:\n")
cat("  P(SVL > u1) ≈", tail_frac_SVL, "(u1 =", u1, ", exp(u1) ≈", exp(u1), "cm)\n")
cat("  P(TL  > u2) ≈", tail_frac_TL,  "(u2 =", u2, ", exp(u2) ≈", exp(u2), "cm)\n")

# 9.3 Parametric bivariate bootstrap (logistic GP, common tail shape)

# Sub-threshold samples for realism
sub_SVL_cc   <- df_complete$log_SVL[df_complete$log_SVL <= u1]
sub_TL_cc    <- df_complete$log_TL[df_complete$log_TL <= u2]
sub_SVL_only <- df_svl_only$log_SVL[df_svl_only$log_SVL <= u1]
sub_TL_only  <- df_tl_only$log_TL[df_tl_only$log_TL <= u2]

n_cc       <- nrow(df_complete)
n_svl_only <- nrow(df_svl_only)
n_tl_only  <- nrow(df_tl_only)

# Logistic EV parameter α from constrained fit
dep_name <- grep("^dep", names(par_eq), value = TRUE)[1]
alpha_logistic_hat <- unname(par_eq[dep_name])
theta_gumbel_hat   <- 1 / alpha_logistic_hat  # Gumbel copula θ ≥ 1

gumbel_cop_hat <- gumbelCopula(param = theta_gumbel_hat, dim = 2)

# Storage for bootstrap draws
boot_sigma1 <- numeric(B_boot)
boot_sigma2 <- numeric(B_boot)
boot_xi     <- numeric(B_boot)
boot_y1star <- numeric(B_boot)
boot_y2star <- numeric(B_boot)
boot_conv   <- logical(B_boot)

set.seed(2026)

for (b in seq_len(B_boot)) {
  ## (1) Simulate complete cases via copula + marginal tails
  U_cc <- rCopula(n_cc, gumbel_cop_hat)
  Y1_cc <- numeric(n_cc)
  Y2_cc <- numeric(n_cc)
  
  for (k in seq_len(n_cc)) {
    u1k <- U_cc[k, 1]
    u2k <- U_cc[k, 2]
    
    # SVL margin
    if (u1k <= 1 - tail_frac_SVL || length(sub_SVL_cc) == 0L) {
      Y1_cc[k] <- if (length(sub_SVL_cc)) sample(sub_SVL_cc, 1L, replace = TRUE) else u1
    } else {
      v1   <- (u1k - (1 - tail_frac_SVL)) / tail_frac_SVL
      exc1 <- evd::qgpd(v1, loc = 0, scale = sigma1_hat, shape = xi_eq_hat)
      Y1_cc[k] <- u1 + exc1
    }
    
    # TL margin
    if (u2k <= 1 - tail_frac_TL || length(sub_TL_cc) == 0L) {
      Y2_cc[k] <- if (length(sub_TL_cc)) sample(sub_TL_cc, 1L, replace = TRUE) else u2
    } else {
      v2   <- (u2k - (1 - tail_frac_TL)) / tail_frac_TL
      exc2 <- evd::qgpd(v2, loc = 0, scale = sigma2_hat, shape = xi_eq_hat)
      Y2_cc[k] <- u2 + exc2
    }
  }
  
  X_biv_boot <- cbind(Y1_cc, Y2_cc)
  
  ## (2) Simulate SVL-only and TL-only marginals
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
  
  ## (3) Build bootstrap likelihood and refit H0: ξ1 = ξ2
  ll_cens_log_boot <- make_ll_cens_logistic_bvpot(
    x      = X_biv_boot,
    u      = u_vec,
    cshape = FALSE,
    cscale = FALSE
  )
  ll_missing_boot <- make_ll_missing_uni(y1_E_M2_boot, y2_M1_E_boot, u1, u2)
  
  loglik_full_boot <- function(theta) {
    ll_cens_log_boot(theta) + ll_missing_boot(theta)
  }
  negloglik_full_boot <- function(p) {
    theta <- c(
      scale1 = p[1],
      shape1 = p[2],  # common ξ
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
  boot_y1star[b] <- u1 - sigma1_b / xi_b
  boot_y2star[b] <- u2 - sigma2_b / xi_b
}

cat("\nParametric bootstrap (common-tail logistic GP):\n")
cat("  Successful replicates:", sum(boot_conv), "out of", B_boot, "\n")

boot_ok <- boot_conv & is.finite(boot_xi) & (boot_xi < 0)

alpha_half <- (1 - ci_level) / 2

# Bootstrap CIs for parameters (still equal-tail)
ci_sigma1 <- quantile(boot_sigma1[boot_ok], probs = c(alpha_half, 1 - alpha_half), na.rm = TRUE)
ci_sigma2 <- quantile(boot_sigma2[boot_ok], probs = c(alpha_half, 1 - alpha_half), na.rm = TRUE)
ci_xi     <- quantile(boot_xi[boot_ok],     probs = c(alpha_half, 1 - alpha_half), na.rm = TRUE)

cat("\n", 100 * ci_level, "% bootstrap CIs for GP parameters (log scale):\n", sep = "")
cat("  sigma1:", ci_sigma1[1], "to", ci_sigma1[2], "\n")
cat("  sigma2:", ci_sigma2[1], "to", ci_sigma2[2], "\n")
cat("  xi    :", ci_xi[1],     "to", ci_xi[2],     "\n")

# One-sided (upper) bootstrap CIs for log endpoints
ci_y1_boot_log_upper <- quantile(boot_y1star[boot_ok], probs = ci_level, na.rm = TRUE)
ci_y2_boot_log_upper <- quantile(boot_y2star[boot_ok], probs = ci_level, na.rm = TRUE)

cat("\nOne-sided ", 100 * ci_level, "% upper bootstrap bounds for log endpoints:\n", sep = "")
cat("  y*_SVL (log) ≤ ", ci_y1_boot_log_upper, "\n", sep = "")
cat("  y*_TL  (log) ≤ ", ci_y2_boot_log_upper, "\n", sep = "")

# Back-transform to original scale
ci_y1_boot_orig_upper <- exp(ci_y1_boot_log_upper)
ci_y2_boot_orig_upper <- exp(ci_y2_boot_log_upper)

cat("\nOne-sided ", 100 * ci_level, "% upper bootstrap bounds for endpoints (original scale):\n", sep = "")
cat("  SVL* ≤ ", ci_y1_boot_orig_upper, "\n", sep = "")
cat("  TL*  ≤ ", ci_y2_boot_orig_upper, "\n", sep = "")

# ============================================================
# 10. Exceedance probability for the Stokes alligator
#      (SVL = 239 cm, TL = 450 cm)
#      P(A | y) = P(E | y) * P(A | E, y)
# ============================================================

# Helper: GPD survival for y > u on the *log*-scale
gpd_survival <- function(y, u, sigma, xi) {
  x <- y - u
  if (x <= 0) return(1)  # not beyond threshold
  if (abs(xi) < 1e-8) {
    # Gumbel limit
    return(exp(-x / sigma))
  } else {
    term <- 1 + xi * x / sigma
    if (term <= 0) return(0)  # beyond fitted endpoint
    return(term^(-1 / xi))
  }
}

# ------------------------------------------------------------
# 10.1 Stokes alligator on log scale
# ------------------------------------------------------------
SVL_stokes <- 239   # cm
TL_stokes  <- 450   # cm

log_SVL_stokes <- log(SVL_stokes)
log_TL_stokes  <- log(TL_stokes)

cat("\nStokes alligator (log scale):\n")
cat("  log(SVL) =", log_SVL_stokes, "\n")
cat("  log(TL)  =", log_TL_stokes,  "\n")

# ------------------------------------------------------------
# 10.2 Univariate predictive exceedance probabilities
#      P(Y_j > y_j^† | y) = P(Y_j > u_j | y) * P(Y_j > y_j^† | Y_j > u_j, y)
# ------------------------------------------------------------
cond_SVL_stokes <- gpd_survival(
  y     = log_SVL_stokes,
  u     = u1,
  sigma = sigma1_hat,
  xi    = xi_eq_hat
)

cond_TL_stokes <- gpd_survival(
  y     = log_TL_stokes,
  u     = u2,
  sigma = sigma2_hat,
  xi    = xi_eq_hat
)

p_SVL_stokes <- tail_frac_SVL * cond_SVL_stokes
p_TL_stokes  <- tail_frac_TL  * cond_TL_stokes

cat("\nUnivariate predictive exceedance probabilities for Stokes:\n")
cat("  P(SVL > 239 cm) ≈", p_SVL_stokes,
    " → return period ≈", 1 / p_SVL_stokes, "individuals\n")
cat("  P(TL  > 450 cm) ≈", p_TL_stokes,
    " → return period ≈", 1 / p_TL_stokes,  "individuals\n")

# ------------------------------------------------------------
# 10.3 Empirical tail fraction P(E | y)
#      E = {Y1 > u1 or Y2 > u2} based on at least one observed trait
# ------------------------------------------------------------
has_SVL <- !is.na(df_log$log_SVL)
has_TL  <- !is.na(df_log$log_TL)

# Individuals with at least one trait observed
n_ind <- sum(has_SVL | has_TL)

# At least one observed exceedance above (u1, u2)
E_emp <- (has_SVL & df_log$log_SVL > u1) |
  (has_TL  & df_log$log_TL  > u2)

n_E <- sum(E_emp, na.rm = TRUE)

tail_frac_E <- n_E / n_ind   # \hat P(E | y)

cat("\nEmpirical tail fraction for at least one exceedance E:\n")
cat("  P(E | y) ≈", tail_frac_E,
    " (n_E =", n_E, "out of n =", n_ind, ")\n")

# ------------------------------------------------------------
# 10.4 Simple bivariate exceedance under independence
#      (for comparison only)
# ------------------------------------------------------------
p_joint_indep <- p_SVL_stokes * p_TL_stokes

cat("\nJoint exceedance under independence approximation:\n")
cat("  P(SVL > 239 cm, TL > 450 cm) ≈", p_joint_indep,
    " → return period ≈", 1 / p_joint_indep, "individuals\n")

# ------------------------------------------------------------
# 10.5 Logistic EV tail model:
#      P(A | y) = P(E | y) * P(A | E, y)
#      A = {SVL > 239 cm, TL > 450 cm}
#      Approximate P(A | E, y) by Monte Carlo from logistic tail model.
# ------------------------------------------------------------

alpha_logistic <- alpha_logistic_hat
theta_gumbel   <- 1 / alpha_logistic  # Gumbel copula parameter θ ≥ 1

cat("\nLogistic EV copula parameter (from POT fit):\n")
cat("  alpha (logistic dep) =", alpha_logistic, "\n")
cat("  theta (Gumbel copula) =", theta_gumbel, "\n")

# Gumbel copula for extremal dependence
gumbel_cop <- gumbelCopula(param = theta_gumbel, dim = 2)

set.seed(2025)
N_mc <- 2e5L

U_mc <- rCopula(N_mc, gumbel_cop)

# Sub-threshold samples for margins (from complete cases)
# (already defined earlier, but re-used here)
#   sub_SVL_cc, sub_TL_cc, u1, u2, sigma1_hat, sigma2_hat,
#   xi_eq_hat, tail_frac_SVL, tail_frac_TL

Y1_mc <- numeric(N_mc)  # SVL on log-scale
Y2_mc <- numeric(N_mc)  # TL  on log-scale

for (k in seq_len(N_mc)) {
  u1k <- U_mc[k, 1]
  u2k <- U_mc[k, 2]
  
  # ---- SVL margin ----
  if (u1k <= 1 - tail_frac_SVL || length(sub_SVL_cc) == 0L) {
    # sub-threshold draw (empirical)
    Y1_mc[k] <- if (length(sub_SVL_cc)) {
      sample(sub_SVL_cc, 1L, replace = TRUE)
    } else {
      u1
    }
  } else {
    # tail draw via GPD
    v1   <- (u1k - (1 - tail_frac_SVL)) / tail_frac_SVL
    exc1 <- evd::qgpd(v1, loc = 0, scale = sigma1_hat, shape = xi_eq_hat)
    Y1_mc[k] <- u1 + exc1
  }
  
  # ---- TL margin ----
  if (u2k <= 1 - tail_frac_TL || length(sub_TL_cc) == 0L) {
    # sub-threshold draw (empirical)
    Y2_mc[k] <- if (length(sub_TL_cc)) {
      sample(sub_TL_cc, 1L, replace = TRUE)
    } else {
      u2
    }
  } else {
    # tail draw via GPD
    v2   <- (u2k - (1 - tail_frac_TL)) / tail_frac_TL
    exc2 <- evd::qgpd(v2, loc = 0, scale = sigma2_hat, shape = xi_eq_hat)
    Y2_mc[k] <- u2 + exc2
  }
}

# Event E (at least one trait beyond threshold) in simulation
E_sim <- (Y1_mc > u1) | (Y2_mc > u2)

# Event A (Stokes-level exceedance in both traits)
A_sim <- (Y1_mc > log_SVL_stokes) & (Y2_mc > log_TL_stokes)

n_E_sim <- sum(E_sim)
if (n_E_sim == 0L) {
  stop("No simulated tail events E_sim; increase N_mc or check model parameters.")
}

# Conditional probability P(A | E, y) under logistic tail model
p_A_given_E_log <- sum(A_sim & E_sim) / n_E_sim

# Unconditional bivariate exceedance probability via decomposition
p_joint_log <- tail_frac_E * p_A_given_E_log

# Monte Carlo standard error for p_A_given_E (conditional) and p_joint_log
se_A_given_E <- sqrt(p_A_given_E_log * (1 - p_A_given_E_log) / n_E_sim)
se_joint_log <- tail_frac_E * se_A_given_E

cat("\nLogistic EV tail model (decomposition P(A | y) = P(E | y) * P(A | E, y)):\n")
cat("  P(A | E, y) ≈", p_A_given_E_log,
    " (MC s.e. ≈", se_A_given_E, ")\n")
cat("  P(SVL > 239 cm, TL > 450 cm | y)\n")
cat("    ≈", p_joint_log,
    " (MC s.e. ≈", se_joint_log, ")\n")
cat("    → return period ≈", 1 / p_joint_log, "individuals\n")

# Also report simulated P(E) for comparison
p_E_sim <- n_E_sim / N_mc
cat("\nSimulated P(E) (from logistic model) ≈", p_E_sim,
    " vs empirical P(E | y) ≈", tail_frac_E, "\n")

# ------------------------------------------------------------
# 10.6 Return periods in years (given annual harvest)
# ------------------------------------------------------------
N_year <- 7600  # approx mean annual harvest 2009–2024

RP_years_indep  <- (1 / p_joint_indep) / N_year
RP_years_logist <- (1 / p_joint_log)   / N_year

cat("\nReturn period in years (assuming ≈", N_year, "harvested per year):\n")
cat("  Independence: ≈", RP_years_indep,  "years\n")
cat("  Logistic tail model: ≈", RP_years_logist, "years\n")

# ============================================================
# 11. Joint bootstrap distribution for (SVL*, TL*)
#      • 2D KDE + HPD contour
#      • Joint MAP = mode of 2D KDE within window
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

# ------------------------------------------------------------
# 11.1 Marginal one-sided 90% upper bounds (for comparison)
# ------------------------------------------------------------

q_upper <- ci_level

ub_svl_marg <- as.numeric(quantile(svl_boot_cm, probs = q_upper, na.rm = TRUE))
ub_tl_marg  <- as.numeric(quantile(tl_boot_cm,  probs = q_upper, na.rm = TRUE))

cat("\n", 100 * ci_level,
    "% one-sided marginal upper bounds for endpoints (original scale):\n", sep = "")
cat("  SVL* ≤ ", ub_svl_marg, " cm\n", sep = "")
cat("  TL*  ≤ ", ub_tl_marg,  " cm\n", sep = "")

# ------------------------------------------------------------
# 11.2 KDE window & restricted draws
# ------------------------------------------------------------

# Central window for plotting/KDE; tweak quantiles if needed
svl_window <- quantile(svl_boot_cm, c(0.01, 0.99), na.rm = TRUE)
tl_window  <- quantile(tl_boot_cm,  c(0.01, 0.99), na.rm = TRUE)

inside_win <- with(
  boot_joint,
  SVL >= svl_window[1] & SVL <= svl_window[2] &
    TL >= tl_window[1] & TL <= tl_window[2]
)

boot_joint_kde <- boot_joint[inside_win, , drop = FALSE]

if (nrow(boot_joint_kde) < 50L) {
  stop("Too few draws inside KDE window; relax svl_window/tl_window.")
}

# ------------------------------------------------------------
# 11.3 2D KDE on restricted draws + joint MAP
# ------------------------------------------------------------

kde_joint <- MASS::kde2d(
  x    = boot_joint_kde$SVL,
  y    = boot_joint_kde$TL,
  n    = 200,
  lims = c(svl_window[1], svl_window[2], tl_window[1], tl_window[2])
)

## Joint MAP = mode of the 2D KDE within this window
ij_max <- which(kde_joint$z == max(kde_joint$z), arr.ind = TRUE)[1, ]
map_svl_joint <- kde_joint$x[ij_max[1]]
map_tl_joint  <- kde_joint$y[ij_max[2]]

cat("\n2D KDE joint MAP (SVL*, TL*):\n")
cat("  SVL*_MAP ≈", round(map_svl_joint, 2), "cm\n")
cat("  TL*_MAP  ≈", round(map_tl_joint, 2), "cm\n")

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

# KDE estimate at each (restricted) bootstrap draw
dens_hat <- eval_kde2d(kde_joint, boot_joint_kde$SVL, boot_joint_kde$TL)

# Sample-based HPD: keep the 90% highest-density draws
ord      <- order(dens_hat, decreasing = TRUE)
cum_prob <- cumsum(rep(1 / length(dens_hat), length(dens_hat)))
in_hpd   <- logical(length(dens_hat))
in_hpd[ord] <- cum_prob <= ci_level

hpd_threshold <- min(dens_hat[in_hpd])

# Data frame for KDE grid (for contour)
kde_df <- with(
  kde_joint,
  expand.grid(SVL = x, TL = y) |>
    transform(z = as.vector(z))
)

# ------------------------------------------------------------
# 11.4 Plot: 2D KDE + HPD contour + joint MAP + Stokes
# ------------------------------------------------------------

p_joint <- ggplot() +
  # 2D KDE (filled levels) on restricted draws
  geom_density_2d_filled(
    data  = boot_joint_kde,
    aes(x = SVL, y = TL),
    alpha = 0.7
  ) +
  # HPD contour from KDE
  geom_contour(
    data  = kde_df,
    aes(x = SVL, y = TL, z = z),
    breaks = hpd_threshold,
    colour = "black",
    linewidth = 1.0
  ) +
  # Bootstrap draws (faint)
  geom_point(
    data  = boot_joint_kde,
    aes(x = SVL, y = TL),
    alpha = 0.25,
    size  = 0.6,
    color = "grey20"
  ) +
  # Joint MAP from 2D KDE
  geom_point(
    aes(x = map_svl_joint, y = map_tl_joint),
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
  coord_cartesian(xlim = svl_window, ylim = tl_window) +
  labs(
    x     = "SVL endpoint (cm)",
    y     = "TL endpoint (cm)",
    fill  = "Density level"
  ) +
  theme_science_polished

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
#      • MAPs = modes of 1D KDEs within same window
# ============================================================

# One-sided upper marginal bounds (already computed in 11.1)
ub_svl <- ub_svl_marg
ub_tl  <- ub_tl_marg

cat("\nOne-sided ", 100 * ci_level,
    "% upper bootstrap bounds for endpoints (original scale):\n", sep = "")
cat("  SVL* ≤ ", ub_svl, " cm\n", sep = "")
cat("  TL*  ≤ ", ub_tl,  " cm\n", sep = "")

# Stokes values
svl_stokes <- SVL_stokes
tl_stokes  <- TL_stokes

# Restrict bootstrap draws to same window for KDE + histograms
svl_in <- svl_boot_cm[svl_boot_cm >= svl_window[1] & svl_boot_cm <= svl_window[2]]
tl_in  <- tl_boot_cm[ tl_boot_cm  >= tl_window[1]  & tl_boot_cm  <= tl_window[2]]

if (length(svl_in) < 30L || length(tl_in) < 30L) {
  warning("Few draws inside univariate KDE window; consider widening svl_window/tl_window.")
}

# Density estimates for smooth curves, only on chosen range
svl_dens <- density(
  svl_in,
  n    = 512,
  from = svl_window[1],
  to   = svl_window[2]
)

tl_dens <- density(
  tl_in,
  n    = 512,
  from = tl_window[1],
  to   = tl_window[2]
)

svl_df      <- data.frame(SVL = svl_dens$x, density = svl_dens$y)
svl_samp_df <- data.frame(SVL = svl_in)

tl_df      <- data.frame(TL = tl_dens$x, density = tl_dens$y)
tl_samp_df <- data.frame(TL = tl_in)

## Marginal MAPs = modes of 1D KDEs
map_svl_marg <- svl_dens$x[which.max(svl_dens$y)]
map_tl_marg  <- tl_dens$x[which.max(tl_dens$y)]

cat("\nMarginal KDE MAPs for endpoints (within window):\n")
cat("  SVL*_MAP, marginal ≈", round(map_svl_marg, 2), "cm\n")
cat("  TL*_MAP,  marginal ≈", round(map_tl_marg,  2), "cm\n")

# ------------------------------------------------------------
# SVL* bootstrap distribution
# ------------------------------------------------------------
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
    size = 1.2
  ) +
  # marginal MAP
  geom_vline(xintercept = map_svl_marg,
             color = "purple", linetype = "dashed", size = 1.2) +
  # one-sided upper bound
  geom_vline(xintercept = ub_svl,
             color = "orange", linetype = "dotdash", size = 1.2) +
  # Stokes SVL
  geom_vline(xintercept = svl_stokes,
             color = "black", size = 1.0) +
  labs(x = "SVL endpoint (cm)", y = "Density") +
  coord_cartesian(xlim = svl_window) +
  theme_science_polished

# ------------------------------------------------------------
# TL* bootstrap distribution
# ------------------------------------------------------------
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
    size = 1.2
  ) +
  # marginal MAP
  geom_vline(xintercept = map_tl_marg,
             color = "purple", linetype = "dashed", size = 1.2) +
  # one-sided upper bound
  geom_vline(xintercept = ub_tl,
             color = "orange", linetype = "dotdash", size = 1.2) +
  # Stokes TL
  geom_vline(xintercept = tl_stokes,
             color = "black", size = 1.0) +
  labs(x = "TL endpoint (cm)", y = "Density") +
  coord_cartesian(xlim = tl_window) +
  theme_science_polished

p_svl
p_tl

ggsave(file.path(FIG_DIR, "SVL_endpoint_bootstrap.png"),
       p_svl, dpi = 600, w = 6.5, h = 4.5, units = "in")
ggsave(file.path(FIG_DIR, "TL_endpoint_bootstrap.png"),
       p_tl,  dpi = 600, w = 6.5, h = 4.5, units = "in")

