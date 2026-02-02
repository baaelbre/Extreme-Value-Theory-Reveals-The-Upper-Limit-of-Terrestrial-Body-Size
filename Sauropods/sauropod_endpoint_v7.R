# ============================================================
# Bivariate EVT with censored likelihood & logistic dependence
# (CF, CH) — Sauropod circumferences
# Full log-likelihood: censored bivariate POT (bvpot.c:nllbvclog)
# + univariate GP contributions for missing partners (evd::dgpd)
# ============================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(ggplot2)
  library(evd)        # fpot, pgpd, dgpd, sep.bvdata
  library(scales)
  library(glue)
  library(forcats)
  library(grid)
  library(MASS)       # for ginv fallback if Hessian singular
  library(readr)
})

set.seed(42)

# ------------------------------------------------------------
# 0. Load compiled C code (bvpot.c → bvpot.so / bvpot.dll)
# ------------------------------------------------------------
# Compile once from a shell / R console:
#   system("R CMD SHLIB bvpot.c")
if (!is.loaded("nllbvclog")) {
  ## Adjust path / extension as needed:
  ## dyn.load("path/to/bvpot.so") or "bvpot.dll" on Windows
  dyn.load("bvpot.dll")
}

# ---------------------------
# Directories & theming
# ---------------------------
FIG_DIR <- "Figures_sauropods_CF_CH"
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
ci_level      <- 0.9
u_scan_lo_q   <- 0.60
u_scan_hi_q   <- 0.9
u_scan_n      <- 50
min_ex_mrl    <- 5
min_ex_fit    <- 5
topN_show     <- 10

trait_names   <- c("CF", "CH")

# Trait-specific anchor quantiles (log scale)
# (reuse your old choices: 0.79 for CF, 0.82 for CH)
thresh_q_opt  <- c(CF = 0.8, CH = 0.84) # 90% for CH is good

# ============================================================
# 1. Ingest and basic cleaning
# ============================================================

DATA_XLSX <- "Data/sauropod_measurements_demic.xlsx"
df_raw    <- read_excel(DATA_XLSX)

stopifnot(all(c("genus and species",
                "femur circ (mm)",
                "humerus circ (mm)") %in% names(df_raw)))

df <- df_raw %>%
  transmute(
    specimen = `genus and species`,
    CF = suppressWarnings(as.numeric(`femur circ (mm)`)),
    CH = suppressWarnings(as.numeric(`humerus circ (mm)`))
  )

# ============================================================
# 2. Completeness diagnostics (CF, CH)
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

write_csv(compl_tbl, file.path(FIG_DIR, "completeness_CF_CH.csv"))

p_compl <- compl_tbl %>%
  mutate(trait = fct_inorder(trait)) %>%
  ggplot(aes(trait, completeness)) +
  geom_col(fill = "#3B82F6") +
  geom_text(aes(label = percent(completeness, accuracy = 0.1)),
            vjust = -0.2, size = 4) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1.10)) +
  labs(title = "Completeness (CF, CH)", x = NULL, y = "Completeness") +
  theme_science_polished
p_compl
ggsave(file.path(FIG_DIR, "completeness_CF_CH.png"), p_compl,
       dpi = 600, w = 6.0, h = 4.2, units = "in")

# ============================================================
# 3. Largest specimens per trait (Top N, raw scale)
# ============================================================

largest_by_trait <- function(tr, topN = 10) {
  df %>%
    filter(!is.na(.data[[tr]])) %>%
    arrange(desc(.data[[tr]])) %>%
    mutate(rank = row_number()) %>%
    slice_head(n = topN) %>%
    transmute(trait = tr,
              rank,
              value_raw = .data[[tr]],
              specimen)
}

largest_list <- map_df(trait_names, largest_by_trait, topN = topN_show)
largest_list
write_csv(largest_list, file.path(FIG_DIR, "largest_CF_CH_top10.csv"))

# ============================================================
# 4. Log-transform and split patterns (complete / CF-only / CH-only)
# ============================================================

df_log <- df %>%
  mutate(
    log_CF = ifelse(is.finite(CF) & CF > 0, log(CF), NA_real_),
    log_CH = ifelse(is.finite(CH) & CH > 0, log(CH), NA_real_)
  )

# Complete cases: both CF and CH observed
df_complete <- df_log %>%
  filter(!is.na(log_CF), !is.na(log_CH))

# CF only
df_cf_only <- df_log %>%
  filter(!is.na(log_CF), is.na(log_CH))

# CH only
df_ch_only <- df_log %>%
  filter(is.na(log_CF), !is.na(log_CH))

cat("\nRow counts by observation pattern:\n")
cat("  complete (CF & CH) :", nrow(df_complete), "\n")
cat("  CF only            :", nrow(df_cf_only), "\n")
cat("  CH only            :", nrow(df_ch_only), "\n")

# ============================================================
# 5. Thresholds on log scale (per trait, using *all* non-missing)
# ============================================================

u0_by_trait <- setNames(numeric(length(trait_names)), trait_names)
for (tr in trait_names) {
  v_log <- df_log[[paste0("log_", tr)]]
  u0_by_trait[tr] <- as.numeric(
    quantile(v_log, thresh_q_opt[tr], na.rm = TRUE)
  )
}
print(u0_by_trait)
print(exp(u0_by_trait))  # thresholds on original scale

# ============================================================
# 6. Helper functions for univariate diagnostics (evd)
# ============================================================

make_thresholds <- function(y,
                            q_from = u_scan_lo_q,
                            q_to   = u_scan_hi_q,
                            n      = u_scan_n) {
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
    theo_q <- u + sigma_hat/xi_hat * (probs^(-xi_hat) - 1)
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
    1 - (1 + xi_hat * exc/sigma_hat)^(-1/xi_hat)
  } else {
    1 - exp(-exc/sigma_hat)
  }
  dfpp <- data.frame(
    Theoretical = sort(F_theo),
    Empirical   = (1:n)/n
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

# ============================================================
# 7. Univariate POT diagnostics per trait (CF, CH)
# ============================================================

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
  
  fit0  <- fit_gpd_at_u(y, u0)
  xi0   <- fit0$shape
  sc0   <- fit0$scale
  cov0  <- fit0$cov
  z     <- xi0 / sqrt(cov0[2,2])
  p_wald <- pnorm(z)   # P(ξ_hat < 0) under asymptotic N(0, Var)
  
  zcrit <- qnorm(1 - (1 - ci_level)/2)
  
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
    
    se_shape <- sqrt(out$cov[2,2])
    df_scan$shape[i]    <- out$shape
    df_scan$shape_lo[i] <- out$shape - zcrit * se_shape
    df_scan$shape_hi[i] <- out$shape + zcrit * se_shape
    
    df_scan$adj[i]      <- adj_scale_fun(out$scale, out$shape, u, u0)
    var_adj <- out$cov[1,1] +
      (u - u0)^2 * out$cov[2,2] -
      2*(u - u0)*out$cov[1,2]
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
  
  tibble(
    trait          = tr_key,
    u0             = u0,
    xi_hat_diag    = xi0,
    sigma_hat_diag = sc0,
    wald_p_xi_lt_0 = p_wald
  )
}

diag_CF <- run_trait_diag("CF")
diag_CH <- run_trait_diag("CH")
marg_diag_tbl <- bind_rows(diag_CF, diag_CH)
write_csv(marg_diag_tbl, file.path(FIG_DIR, "gpd_diag_CF_CH.csv"))
print(marg_diag_tbl)

# ============================================================
# 8. Univariate POT fits using all non-missing data per margin
# ============================================================

y_CF_all <- df_log$log_CF[!is.na(df_log$log_CF)]
y_CH_all <- df_log$log_CH[!is.na(df_log$log_CH)]

u1 <- u0_by_trait["CF"]
u2 <- u0_by_trait["CH"]

fit_CF_uni <- evd::fpot(y_CF_all, threshold = u1, model = "gpd")
fit_CH_uni <- evd::fpot(y_CH_all, threshold = u2, model = "gpd")

par_CF_uni <- fit_CF_uni$estimate
par_CH_uni <- fit_CH_uni$estimate

sigma1_hat_uni <- unname(par_CF_uni["scale"])
xi1_hat_uni    <- unname(par_CF_uni["shape"])

sigma2_hat_uni <- unname(par_CH_uni["scale"])
xi2_hat_uni    <- unname(par_CH_uni["shape"])

cat("\nUnivariate GP fits (log scale, evd::fpot, ALL marginal data):\n")
cat("  CF: scale =", sigma1_hat_uni, " shape =", xi1_hat_uni, "\n")
cat("  CH: scale =", sigma2_hat_uni, " shape =", xi2_hat_uni, "\n")

# ============================================================
# 9. Bivariate POT with censored likelihood (complete cases)
#    + univariate contributions from missings (E1M2, M1E2)
# ============================================================

# 9.1 Complete-case matrix for censored bivariate logistic POT
X_biv <- df_complete %>%
  transmute(
    log_CF = log_CF,
    log_CH = log_CH
  ) %>%
  as.matrix()

n_biv <- nrow(X_biv)
cat("\nBivariate sample size (complete cases CF & CH):", n_biv, "\n")

u_vec <- c(u1, u2)

quad_tab <- table(
  CF_exceed = X_biv[,1] > u1,
  CH_exceed = X_biv[,2] > u2
)
cat("\nQuadrant table (on log scale, complete cases only):\n")
print(quad_tab)

# 9.2 Missing patterns that contribute to the tail likelihood

# Only cases with an exceedance and the other margin truly missing:
y1_E_M2 <- df_cf_only$log_CF[df_cf_only$log_CF > u1]  # CF exceeds, CH missing
y2_M1_E <- df_ch_only$log_CH[df_ch_only$log_CH > u2]  # CH exceeds, CF missing

cat("\nCounts of missing-partner exceedances:\n")
cat("  CF exceed, CH missing (E1M2):", length(y1_E_M2), "\n")
cat("  CH exceed, CF missing (M1E2):", length(y2_M1_E), "\n")

# ---------------------------
# 9.3 Censored logistic log-likelihood using nllbvclog from bvpot.c
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
              dns    = as.double(0)             # start at 0
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
# 9.4 Univariate GP contributions from missings (evd::dgpd)
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
# 10. Maximum likelihood: free shapes vs common shape
# ============================================================

# 10.1 Free shapes (H1)
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
  method  = "Nelder-Mead",
  hessian = TRUE
)

par_free <- opt_free$par
names(par_free) <- names(start_free)

ell_hat_free <- -opt_free$value

# robust-ish covariance from Hessian
V_free <- tryCatch(
  solve(opt_free$hessian),
  error = function(e) MASS::ginv(opt_free$hessian)
)

cat("\nBivariate MLEs (full llh, free shapes):\n")
print(par_free)
cat("Maximized log-likelihood (free) =", ell_hat_free, "\n")

# 10.2 Common tail shape (H0: xi1 = xi2 = xi)
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
  method  = "Nelder-Mead",
  hessian = TRUE
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

V_eq_phi <- tryCatch(
  solve(opt_eq$hessian),
  error = function(e) MASS::ginv(opt_eq$hessian)
)

cat("\nCommon tail shape MLEs (full llh, H0: xi_CF = xi_CH):\n")
print(par_eq)
cat("Maximized log-likelihood (H0) =", ell_hat_eq, "\n")

# Likelihood ratio test for H0: xi1 = xi2
LR_shape <- 2 * (ell_hat_free - ell_hat_eq)
df_LR    <- 1
p_LR     <- 1 - pchisq(LR_shape, df = df_LR)

cat("\nLikelihood ratio test for tail equality H0: xi_CF = xi_CH\n")
cat("  LR statistic =", LR_shape, " with df =", df_LR, "\n")
cat("  p-value      =", p_LR, "\n")

# ============================================================
# 11. Endpoint estimation under common tail shape (H0)
#     + delta-method standard errors and CIs
# ============================================================

# Extract constrained parameter estimates
xi_eq_hat  <- unname(par_eq_phi["xi"])
sigma1_hat <- unname(par_eq_phi["scale1"])  # CF
sigma2_hat <- unname(par_eq_phi["scale2"])  # CH

if (xi_eq_hat >= 0) {
  warning("Common tail shape xi_eq_hat >= 0: finite endpoints not defined.")
}

# 11.1 Point estimates of endpoints on log scale
# GP endpoint: y* = u - sigma / xi  (for xi < 0)
ystar1_hat_log <- as.numeric(u1 - sigma1_hat / xi_eq_hat)  # CF*
ystar2_hat_log <- as.numeric(u2 - sigma2_hat / xi_eq_hat)  # CH*

cat("\nLog-scale endpoint estimates (common-tail GP, H0):\n")
cat("  y*_CF (log) =", ystar1_hat_log, "\n")
cat("  y*_CH (log) =", ystar2_hat_log, "\n")

# 11.2 Delta-method variances for (y*_CF, y*_CH)
# Work with φ = (scale1, xi, scale2, dep)  = par_eq_phi
param_phi_names <- names(par_eq_phi)
k_phi           <- length(par_eq_phi)

idx_scale1_phi <- which(param_phi_names == "scale1")
idx_xi_phi     <- which(param_phi_names == "xi")
idx_scale2_phi <- which(param_phi_names == "scale2")
# dep index exists but endpoints do not depend on it

# Gradient of y*_CF wrt φ
# y*_CF = u1 - scale1 / xi
g1_phi <- rep(0, k_phi)
g1_phi[idx_scale1_phi] <- -1 / xi_eq_hat
g1_phi[idx_xi_phi]     <-  sigma1_hat / (xi_eq_hat^2)

# Gradient of y*_CH wrt φ
# y*_CH = u2 - scale2 / xi
g2_phi <- rep(0, k_phi)
g2_phi[idx_scale2_phi] <- -1 / xi_eq_hat
g2_phi[idx_xi_phi]     <-  sigma2_hat / (xi_eq_hat^2)

# Delta-method variances
var_y1_log <- as.numeric(t(g1_phi) %*% V_eq_phi %*% g1_phi)
var_y2_log <- as.numeric(t(g2_phi) %*% V_eq_phi %*% g2_phi)

se_y1_log  <- sqrt(var_y1_log)
se_y2_log  <- sqrt(var_y2_log)

# 11.3 Confidence intervals on log scale
zcrit <- qnorm(1 - (1 - ci_level) / 2)

ci_y1_log <- ystar1_hat_log + c(-1, 1) * zcrit * se_y1_log
ci_y2_log <- ystar2_hat_log + c(-1, 1) * zcrit * se_y2_log

cat("\n", 100 * ci_level,
    "% marginal CIs for log endpoints (common-tail, H0):\n", sep = "")
cat("  y*_CF (log) CI: [", ci_y1_log[1], ", ", ci_y1_log[2], "]\n", sep = "")
cat("  y*_CH (log) CI: [", ci_y2_log[1], ", ", ci_y2_log[2], "]\n", sep = "")

# 11.4 Transform to original scale (CF*, CH*)
ystar1_hat_orig <- exp(ystar1_hat_log)
ystar2_hat_orig <- exp(ystar2_hat_log)

ci_y1_orig <- exp(ci_y1_log)
ci_y2_orig <- exp(ci_y2_log)

cat("\n", 100 * ci_level,
    "% marginal CIs for endpoints (original scale, common-tail, H0):\n", sep = "")
cat("  CF* =", ystar1_hat_orig,
    "  CI: [", ci_y1_orig[1], ", ", ci_y1_orig[2], "]\n", sep = "")
cat("  CH* =", ystar2_hat_orig,
    "  CI: [", ci_y2_orig[1], ", ", ci_y2_orig[2], "]\n", sep = "")

# 11.5 Joint CI ellipse in log space
G_eq <- rbind(g1_phi, g2_phi)           # 2 x k_phi
Sigma_y_eq <- G_eq %*% V_eq_phi %*% t(G_eq)
mu_y_eq    <- c(ystar1_hat_log, ystar2_hat_log)

colnames(Sigma_y_eq) <- rownames(Sigma_y_eq) <- c("y_CF_star_log", "y_CH_star_log")

cat("\nJoint covariance matrix for log endpoints (H0):\n")
print(Sigma_y_eq)

c2    <- qchisq(ci_level, df = 2)
theta <- seq(0, 2 * pi, length.out = 400)

C_sy_eq <- chol(Sigma_y_eq)

ellipse_log <- sapply(theta, function(th) {
  v <- c(cos(th), sin(th))
  mu_y_eq + sqrt(c2) * (t(C_sy_eq) %*% v)
})
ellipse_log <- as.data.frame(t(ellipse_log))
colnames(ellipse_log) <- c("y_CF_log", "y_CH_log")

pt_hat_log <- data.frame(
  y_CF_log = ystar1_hat_log,
  y_CH_log = ystar2_hat_log
)

p_joint_log <- ggplot() +
  geom_path(data = ellipse_log,
            aes(x = y_CF_log, y = y_CH_log)) +
  geom_point(data = pt_hat_log,
             aes(x = y_CF_log, y = y_CH_log),
             size = 2, color = "red") +
  labs(
    x = "log endpoint CF*",
    y = "log endpoint CH*",
    title = glue("{round(100 * ci_level)}% joint CI for log endpoints (CF*, CH*)")
  ) +
  theme_science_polished

ggsave(file.path(FIG_DIR, "joint_endpoint_log_ellipse_CF_CH.png"),
       p_joint_log, dpi = 600, w = 5.5, h = 5.0, units = "in")

ellipse_orig <- ellipse_log %>%
  mutate(
    CF_star = exp(y_CF_log),
    CH_star = exp(y_CH_log)
  )

pt_hat_orig <- data.frame(
  CF_star = ystar1_hat_orig,
  CH_star = ystar2_hat_orig
)

p_joint_orig <- ggplot() +
  geom_path(data = ellipse_orig,
            aes(x = CF_star, y = CH_star)) +
  geom_point(data = pt_hat_orig,
             aes(x = CF_star, y = CH_star),
             size = 2, color = "red") +
  labs(
    x = "Endpoint CF*",
    y = "Endpoint CH*",
    title = glue("{round(100 * ci_level)}% joint CI for endpoints (CF*, CH*)")
  ) +
  theme_science_polished

ggsave(file.path(FIG_DIR, "joint_endpoint_orig_ellipse_CF_CH.png"),
       p_joint_orig, dpi = 600, w = 5.5, h = 5.0, units = "in")

# ============================================================
# 12. Transform CF* + CH* endpoint to mass endpoint
#     Campione–Evans: log10 M = -1.104 + 2.75 log10(CF + CH)
# ============================================================

# Point estimate for combined circumference endpoint
CF_star  <- ystar1_hat_orig   # endpoint CF* on original scale (mm)
CH_star  <- ystar2_hat_orig   # endpoint CH* on original scale (mm)
Csum_star <- CF_star + CH_star

# Campione–Evans coefficients
alpha_mass <- -1.104
beta_mass  <-  2.75

# Endpoint mass on log10 and original scale
log10_M_star <- alpha_mass + beta_mass * log10(Csum_star)
M_star       <- 10^log10_M_star

cat("\nEndpoint for combined circumference S* = CF* + CH*:\n")
cat("  S* =", Csum_star, " (same units as CF, CH)\n")

cat("\nMass endpoint implied by Campione–Evans scaling:\n")
cat("  log10 M* =", log10_M_star, "\n")
cat("  M*       =", M_star/1e6, " (same mass units as the regression)\n")


