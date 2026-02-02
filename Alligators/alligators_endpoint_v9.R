# ============================================================
# Bivariate EVT with censored likelihood & logistic dependence
# (SVL, TL) — American alligators
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
})

set.seed(42)

# ------------------------------------------------------------
# 0. Load compiled C code (bvpot.c → bvpot.so / bvpot.dll)
# ------------------------------------------------------------
#  system("R CMD SHLIB bvpot.c")
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
ci_level      <- 0.95
u_scan_lo_q   <- 0.60
u_scan_hi_q   <- 0.99
u_scan_n      <- 50
min_ex_mrl    <- 5
min_ex_fit    <- 20

trait_names   <- c("SVL","TL")
thresh_q_opt  <- c(SVL = 0.94, TL = 0.94)  # anchor quantiles (log scale)

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
  geom_text(aes(label = percent(completeness, accuracy=0.1)),
            vjust = -0.2, size = 4) +
  scale_y_continuous(labels = percent_format(accuracy=1), limits = c(0, 1.10)) +
  labs(title = "Completeness (SVL, TL)", x = NULL, y = "Completeness") +
  theme_science_polished

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
    theo_q <- u + sigma_hat/xi_hat * (probs^(-xi_hat) - 1)
  } else {
    theo_q <- u - sigma_hat * log(probs)
  }
  
  dfqq <- data.frame(
    Theoretical = rev(theo_q),
    Empirical   = sort(y[y>u])
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
    ex <- y[y>u] - u
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
  
  tibble(
    trait          = tr_key,
    u0             = u0,
    xi_hat_diag    = xi0,
    sigma_hat_diag = sc0,
    wald_p_xi_lt_0 = p_wald
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
  SVL_exceed = X_biv[,1] > u1,
  TL_exceed  = X_biv[,2] > u2
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
#    + delta-method standard errors and CIs
# ============================================================

# Extract constrained parameter estimates
xi_eq_hat  <- unname(par_eq_phi["xi"])
sigma1_hat <- unname(par_eq_phi["scale1"])
sigma2_hat <- unname(par_eq_phi["scale2"])

# 9.1 Point estimates of endpoints on log scale
# GP endpoint: y* = u - sigma / xi  (for xi < 0)
ystar1_hat_log <- as.numeric(u1 - sigma1_hat / xi_eq_hat)  # SVL*
ystar2_hat_log <- as.numeric(u2 - sigma2_hat / xi_eq_hat)  # TL*

cat("\nLog-scale endpoint estimates (common-tail GP, H0):\n")
cat("  y*_SVL (log) =", ystar1_hat_log, "\n")
cat("  y*_TL  (log) =", ystar2_hat_log, "\n")

# 9.2 Delta-method variances for (y*_SVL, y*_TL)
# Work with φ = (scale1, xi, scale2, dep)  = par_eq_phi
param_phi_names <- names(par_eq_phi)
k_phi           <- length(par_eq_phi)

idx_scale1_phi <- which(param_phi_names == "scale1")
idx_xi_phi     <- which(param_phi_names == "xi")
idx_scale2_phi <- which(param_phi_names == "scale2")
# dep index exists but endpoints do not depend on it

# Gradient of y*_SVL wrt φ
# y*_SVL = u1 - scale1 / xi
g1_phi <- rep(0, k_phi)
g1_phi[idx_scale1_phi] <- -1 / xi_eq_hat
g1_phi[idx_xi_phi]     <-  sigma1_hat / (xi_eq_hat^2)

# Gradient of y*_TL wrt φ
# y*_TL = u2 - scale2 / xi
g2_phi <- rep(0, k_phi)
g2_phi[idx_scale2_phi] <- -1 / xi_eq_hat
g2_phi[idx_xi_phi]     <-  sigma2_hat / (xi_eq_hat^2)

# Delta-method variances
var_y1_log <- as.numeric(t(g1_phi) %*% V_eq_phi %*% g1_phi)
var_y2_log <- as.numeric(t(g2_phi) %*% V_eq_phi %*% g2_phi)

se_y1_log  <- sqrt(var_y1_log)
se_y2_log  <- sqrt(var_y2_log)

# 9.3 Confidence intervals on log scale
zcrit <- qnorm(1 - (1 - ci_level) / 2)

ci_y1_log <- ystar1_hat_log + c(-1, 1) * zcrit * se_y1_log
ci_y2_log <- ystar2_hat_log + c(-1, 1) * zcrit * se_y2_log

cat("\n", 100 * ci_level,
    "% marginal CIs for log endpoints (common-tail, H0):\n", sep = "")
cat("  y*_SVL (log) CI: [", ci_y1_log[1], ", ", ci_y1_log[2], "]\n", sep = "")
cat("  y*_TL  (log) CI: [", ci_y2_log[1], ", ", ci_y2_log[2], "]\n", sep = "")

# 9.4 Transform to original scale (SVL*, TL*)
ystar1_hat_orig <- exp(ystar1_hat_log)
ystar2_hat_orig <- exp(ystar2_hat_log)

ci_y1_orig <- exp(ci_y1_log)
ci_y2_orig <- exp(ci_y2_log)

cat("\n", 100 * ci_level,
    "% marginal CIs for endpoints (original scale, common-tail, H0):\n", sep = "")
cat("  SVL* =", ystar1_hat_orig,
    "  CI: [", ci_y1_orig[1], ", ", ci_y1_orig[2], "]\n", sep = "")
cat("  TL*  =", ystar2_hat_orig,
    "  CI: [", ci_y2_orig[1], ", ", ci_y2_orig[2], "]\n", sep = "")

# ============================================================
# 10. Exceedance probability for the Stokes alligator
#      (SVL = 239 cm, TL = 450 cm)
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
    return(term^(-1/xi))
  }
}

# ------------------------------------------------------------
# 10.1 Stokes alligator measurements on log scale
# ------------------------------------------------------------

SVL_stokes <- 239   # cm
TL_stokes  <- 450   # cm

log_SVL_stokes <- log(SVL_stokes)
log_TL_stokes  <- log(TL_stokes)

cat("\nStokes alligator (log scale):\n")
cat("  log(SVL) =", log_SVL_stokes, "\n")
cat("  log(TL)  =", log_TL_stokes,  "\n")

# ------------------------------------------------------------
# 10.2 Univariate tail fractions for SVL and TL
# ------------------------------------------------------------

# number of usable observations per trait
n_SVL_total <- sum(!is.na(df_log$log_SVL))
n_TL_total  <- sum(!is.na(df_log$log_TL))

# number of exceedances above the chosen thresholds u1, u2
n_SVL_exceed_u1 <- sum(df_log$log_SVL > u1, na.rm = TRUE)
n_TL_exceed_u2  <- sum(df_log$log_TL  > u2, na.rm = TRUE)

tail_frac_SVL <- n_SVL_exceed_u1 / n_SVL_total
tail_frac_TL  <- n_TL_exceed_u2  / n_TL_total

cat("\nUnivariate tail fractions (empirical):\n")
cat("  P(SVL > u1) ≈", tail_frac_SVL, "(u1 =", u1, ", exp(u1) ≈", exp(u1), "cm)\n")
cat("  P(TL  > u2) ≈", tail_frac_TL,  "(u2 =", u2, ", exp(u2) ≈", exp(u2), "cm)\n")

# ------------------------------------------------------------
# 10.3 Univariate exceedance probabilities for Stokes
#      using the marginal GPD fits (fit_SVL_uni, fit_TL_uni)
# ------------------------------------------------------------

# conditional survival in the GPD tail
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

# predictive exceedance probabilities (univariate)
p_SVL_stokes <- tail_frac_SVL * cond_SVL_stokes
p_TL_stokes  <- tail_frac_TL  * cond_TL_stokes

cat("\nUnivariate predictive exceedance probabilities for Stokes:\n")
cat("  P(SVL > 239 cm) ≈", p_SVL_stokes,
    " → return period ≈", 1 / p_SVL_stokes, "individuals\n")
cat("  P(TL  > 450 cm) ≈", p_TL_stokes,
    " → return period ≈", 1 / p_TL_stokes,  "individuals\n")



# ------------------------------------------------------------
# 10.4 Simple bivariate exceedance probability (independence
#      approximation): P(SVL>239, TL>450) ≈ product of marginals
# ------------------------------------------------------------

p_joint_indep <- p_SVL_stokes * p_TL_stokes

cat("\nJoint exceedance under independence approximation:\n")
cat("  P(SVL > 239 cm, TL > 450 cm) ≈", p_joint_indep,
    " → return period ≈", 1 / p_joint_indep, "individuals\n")

# ------------------------------------------------------------
# 10.5 Logistic EV copula (Gumbel) Monte Carlo for joint tail
#      P(SVL > 239 cm, TL > 450 cm)
# ------------------------------------------------------------
# We use the fitted logistic POT dependence parameter 'dep'
# and interpret it as the logistic EV parameter α.
# The associated Gumbel–Hougaard copula has parameter θ = 1 / α.

alpha_logistic <- unname(par_eq["dep.dep"])

theta_gumbel <- 1 / alpha_logistic  # Gumbel parameter θ ≥ 1

cat("\nLogistic EV copula parameter (from POT fit):\n")
cat("  alpha (logistic dep) =", alpha_logistic, "\n")
cat("  theta (Gumbel copula) =", theta_gumbel, "\n")

# --- (a) CDF thresholds for the Stokes alligator on [0,1] ------------------
# For y > u_j:
#   F_j(y) = 1 - p_j * S_j(y | u_j)
# where p_j is the empirical tail fraction and S_j is the GPD survival.

u1_u <- 1 - tail_frac_SVL * cond_SVL_stokes  # F_SVL(239 cm)
u2_u <- 1 - tail_frac_TL  * cond_TL_stokes   # F_TL (450 cm)

cat("\nCDF thresholds for Stokes (on [0,1] scale):\n")
cat("  u1_u = F_SVL(239 cm) ≈", u1_u, "\n")
cat("  u2_u = F_TL (450 cm) ≈", u2_u, "\n")

# Sanity: thresholds must be in (0, 1)
if (!is.finite(u1_u) || !is.finite(u2_u) ||
    u1_u <= 0 || u1_u >= 1 ||
    u2_u <= 0 || u2_u >= 1) {
  stop("Invalid CDF thresholds u1_u / u2_u. Check tail fractions and GPD fits.")
}

# --- (b) Monte Carlo from the Gumbel (logistic EV) copula ------------------
set.seed(2025)  # separate seed for copula MC
N_mc <- 2e5L    # increase if you want smaller MC error

gumbel_cop <- gumbelCopula(param = theta_gumbel, dim = 2)
U_mc       <- rCopula(N_mc, gumbel_cop)

ind_joint_log <- (U_mc[,1] > u1_u) & (U_mc[,2] > u2_u)
p_joint_log   <- mean(ind_joint_log)
se_joint_log  <- sqrt(p_joint_log * (1 - p_joint_log) / N_mc)

cat("\nLogistic EV copula Monte Carlo estimate (bivariate tail):\n")
cat("  P(SVL > 239 cm, TL > 450 cm)\n")
cat("    ≈", p_joint_log,
    " (MC s.e. ≈", se_joint_log, ")\n")
cat("    → return period ≈", 1 / p_joint_log, "individuals\n")

cat("\nComparison independence vs logistic copula:\n")
cat("  Indep:  P ≈", p_joint_indep,
    " → RP ≈", 1 / p_joint_indep, "\n")
cat("  Logist: P ≈", p_joint_log,
    " → RP ≈", 1 / p_joint_log,  "\n")

N_year <- 7600  # approx mean annual harvest 2009–2024

RP_years_indep  <- (1 / p_joint_indep) / N_year
RP_years_logist <- (1 / p_joint_log)   / N_year

cat("Return period (independence) in years ≈", RP_years_indep,  "\n")
cat("Return period (logistic GP) in years ≈", RP_years_logist, "\n")
