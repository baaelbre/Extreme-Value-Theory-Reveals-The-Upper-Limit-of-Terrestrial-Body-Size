# ============================================================
# Bivariate EVT with censored likelihood & logistic dependence
# (SVL, TL) — American alligators
# ============================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(ggplot2)
  library(extRemes)   # for univariate diag plots (pevd)
  library(evd)        # fpot, fbvpot
  library(scales)
  library(glue)
  library(forcats)
  library(grid)
})

set.seed(42)

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
# 3. Log-transform and thresholds (log scale)
# ============================================================

df_log <- df %>%
  mutate(
    log_SVL = ifelse(is.finite(SVL) & SVL > 0, log(SVL), NA_real_),
    log_TL  = ifelse(is.finite(TL)  & TL  > 0, log(TL),  NA_real_)
  )

u0_by_trait <- setNames(numeric(length(trait_names)), trait_names)
for (tr in trait_names) {
  v_log <- df_log[[paste0("log_", tr)]]
  u0_by_trait[tr] <- as.numeric(quantile(v_log, thresh_q_opt[tr], na.rm = TRUE))
}
print(u0_by_trait)
print(exp(u0_by_trait))  # thresholds on original scale

# ============================================================
# 4. Optional univariate POT diagnostics (using evd::fpot)
# ============================================================

make_thresholds <- function(y, q_from = u_scan_lo_q, q_to = u_scan_hi_q, n = u_scan_n) {
  rng <- quantile(y, c(q_from, q_to), na.rm = TRUE)
  seq(rng[1], rng[2], length.out = n)
}

fit_gpd_at_u <- function(y, u) {
  fit <- fpot(y, threshold = u, model = "gpd")
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
  
  F_hat <- pevd(exc, scale = sigma_hat, shape = xi_hat, type = "GP")
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
# 5. Bivariate POT with censored likelihood (three tail quadrants)
#    Logistic dependence via evd::fbvpot
# ============================================================

# Work only with specimens where both SVL and TL observed
X_biv <- df_log %>%
  filter(!is.na(log_SVL), !is.na(log_TL)) %>%
  transmute(
    log_SVL = log_SVL,
    log_TL  = log_TL
  ) %>%
  as.matrix()

n_biv <- nrow(X_biv)
cat("Bivariate sample size (both SVL & TL observed):", n_biv, "\n")

u1 <- u0_by_trait["SVL"]
u2 <- u0_by_trait["TL"]

u_vec <- c(u1, u2)

# Quadrant counts relative to thresholds (for sanity check)
quad_tab <- table(
  SVL_exceed = X_biv[,1] > u1,
  TL_exceed  = X_biv[,2] > u2
)
cat("\nQuadrant table (on log scale):\n")
print(quad_tab)
# R10, R01, R11 are the three tail quadrants used by the censored likelihood

# ------------------------------------------------------------
# 5.1 Univariate GP fits on log scale (evd::fpot, same thresholds)
# ------------------------------------------------------------

fit_SVL_uni <- fpot(X_biv[,1], threshold = u1, model = "gpd")
fit_TL_uni  <- fpot(X_biv[,2], threshold = u2, model = "gpd")

par_SVL_uni <- fit_SVL_uni$estimate
par_TL_uni  <- fit_TL_uni$estimate

sigma1_hat_uni <- unname(par_SVL_uni["scale"])
xi1_hat_uni    <- unname(par_SVL_uni["shape"])

sigma2_hat_uni <- unname(par_TL_uni["scale"])
xi2_hat_uni    <- unname(par_TL_uni["shape"])

cat("\nUnivariate GP fits (log scale, evd::fpot):\n")
cat("  SVL: scale =", sigma1_hat_uni, " shape =", xi1_hat_uni, "\n")
cat("  TL : scale =", sigma2_hat_uni, " shape =", xi2_hat_uni, "\n")

# ------------------------------------------------------------
# 5.2 Bivariate logistic POT with censored likelihood
#     (uses all three tail quadrants R10, R01, R11)
# ------------------------------------------------------------

fit_bv_free <- evd::fbvpot(
  x          = X_biv,
  threshold  = u_vec,
  model      = "log",       # logistic dependence
  likelihood = "censored",  # uses R10, R01, R11 + count of R00
  std.err    = TRUE,
  method = "Nelder-Mead"
)

cat("\nBivariate logistic POT fit (fbvpot, censored likelihood):\n")
print(fit_bv_free)

est_free <- fit_bv_free$estimate

sigma1_hat_biv <- est_free["scale1"]
xi1_hat_biv    <- est_free["shape1"]
sigma2_hat_biv <- est_free["scale2"]
xi2_hat_biv    <- est_free["shape2"]
dep_hat        <- est_free["dep"]   # logistic dep parameter (1 ≈ independence)

ell_hat_free   <- as.numeric(logLik(fit_bv_free))

cat("\nBivariate MLEs (standard parametrization):\n")
cat("  Margin 1 (SVL, log):  scale1 =", sigma1_hat_biv, "  shape1 =", xi1_hat_biv, "\n")
cat("  Margin 2 (TL,  log):  scale2 =", sigma2_hat_biv, "  shape2 =", xi2_hat_biv, "\n")
cat("  Logistic dep (dep)   =", dep_hat, "\n")
cat("  Maximized log-likelihood (free) =", ell_hat_free, "\n")

# ============================================================
# 6. Likelihood ratio test for tail equality: H0: shape1 = shape2
# ============================================================

# Constrained model with common shape: shape2 = shape1
fit_bv_eq <- fbvpot(
  x          = X_biv,
  threshold  = u_vec,
  model      = "log",
  likelihood = "censored",
  cshape     = TRUE,   # constrain shape2 = shape1
  std.err    = TRUE,
  method = "Nelder-Mead"
)
cat("\nConstrained bivariate fit under H0: shape1 = shape2:\n")
print(fit_bv_eq)

est_eq      <- fit_bv_eq$estimate
xi_eq_hat   <- est_eq["shape1"]     # common shape
sigma1_hat <- est_eq["scale1"]
sigma2_hat <- est_eq["scale2"]
ell_hat_eq  <- as.numeric(logLik(fit_bv_eq))

cat("\nCommon tail shape (H0): xi =", xi_eq_hat, "\n")
cat("Maximized log-likelihood (H0) =", ell_hat_eq, "\n")

# Manual LR statistic
LR_shape <- 2 * (ell_hat_free - ell_hat_eq)
df_LR    <- 1
p_LR     <- 1 - pchisq(LR_shape, df = df_LR)

cat("\nLikelihood ratio test for tail equality H0: xi1 = xi2\n")
cat("  LR statistic =", LR_shape, " with df =", df_LR, "\n")
cat("  p-value      =", p_LR, "\n")

# ============================================================
# 7. Endpoint estimation and joint CI for (SVL*, TL*)
#     based on equal-tail censored bivariate fit (fbvpot, cshape=TRUE)
# ============================================================

# --- 7.1 Extract equal-tail parameters and covariance ------------------
par_eq <- fit_bv_eq$estimate
V_eq   <- fit_bv_eq$var.cov

cat("\nEqual-tail parameter estimates (fbvpot, cshape=TRUE):\n")
print(par_eq)
cat("\nCovariance matrix of equal-tail estimates:\n")
print(V_eq)

# Common shape (xi < 0) and marginal scales
xi_eq_hat   <- unname(par_eq["shape1"])  # common shape
sigma1_hat  <- unname(par_eq["scale1"])  # SVL scale (log scale)
sigma2_hat  <- unname(par_eq["scale2"])  # TL scale (log scale)

# Thresholds on log scale were defined earlier
# u1 <- u0_by_trait["SVL"]; u2 <- u0_by_trait["TL"]

# --- 7.2 Point estimates of endpoints on log scale ---------------------
# GP endpoint for margin j: y*_j = u_j - sigma_j / xi
ystar1_hat <- as.numeric(u1 - sigma1_hat / xi_eq_hat)  # SVL endpoint (log scale)
ystar2_hat <- as.numeric(u2 - sigma2_hat / xi_eq_hat)  # TL  endpoint (log scale)

cat("\nLog-scale endpoint estimates (equal-tail GP):\n")
cat("  y*_SVL (log) =", ystar1_hat, "\n")
cat("  y*_TL  (log) =", ystar2_hat, "\n")

# --- 7.3 Delta-method gradients for (y*_SVL, y*_TL) -------------------
# Parameter ordering must match names(par_eq) and dimnames(V_eq)
param_names <- names(par_eq)
k_par       <- length(par_eq)

idx_scale1  <- which(param_names == "scale1")
idx_scale2  <- which(param_names == "scale2")
idx_shape1  <- which(param_names == "shape1")
# dep parameter present but gradient == 0 w.r.t endpoints

# Gradient for y*_SVL = u1 - scale1 / xi
g1 <- rep(0, k_par)
g1[idx_scale1] <- -1 / xi_eq_hat
g1[idx_shape1] <-  sigma1_hat / (xi_eq_hat^2)

# Gradient for y*_TL = u2 - scale2 / xi
g2 <- rep(0, k_par)
g2[idx_scale2] <- -1 / xi_eq_hat
g2[idx_shape1] <-  sigma2_hat / (xi_eq_hat^2)

# --- 7.4 Univariate delta-method SEs and marginal CIs -----------------
var_y1 <- as.numeric(t(g1) %*% V_eq %*% g1)
var_y2 <- as.numeric(t(g2) %*% V_eq %*% g2)

se_y1  <- sqrt(var_y1)
se_y2  <- sqrt(var_y2)

zcrit <- qnorm(1 - (1 - ci_level) / 2)

ci_y1_log <- ystar1_hat + c(-1, 1) * zcrit * se_y1
ci_y2_log <- ystar2_hat + c(-1, 1) * zcrit * se_y2

cat("\n", 100 * ci_level, "% marginal CIs for log endpoints:\n", sep = "")
cat("  y*_SVL (log) CI: [", ci_y1_log[1], ", ", ci_y1_log[2], "]\n", sep = "")
cat("  y*_TL  (log) CI: [", ci_y2_log[1], ", ", ci_y2_log[2], "]\n", sep = "")

# Transform to original scale
ystar1_hat_orig <- exp(ystar1_hat)
ystar2_hat_orig <- exp(ystar2_hat)

ci_y1_orig <- exp(ci_y1_log)
ci_y2_orig <- exp(ci_y2_log)

cat("\n", 100 * ci_level, "% marginal CIs for endpoints (original scale):\n", sep = "")
cat("  SVL* =", ystar1_hat_orig,
    "  CI: [", ci_y1_orig[1], ", ", ci_y1_orig[2], "]\n", sep = "")
cat("  TL*  =", ystar2_hat_orig,
    "  CI: [", ci_y2_orig[1], ", ", ci_y2_orig[2], "]\n", sep = "")

# --- 7.5 Joint (bivariate) CI on log scale via chi^2_2 ellipse --------
G       <- rbind(g1, g2)                  # 2 x k
Sigma_y <- G %*% V_eq %*% t(G)            # 2 x 2 covariance for (y*_SVL, y*_TL)
mu_y    <- c(ystar1_hat, ystar2_hat)      # mean vector

colnames(Sigma_y) <- rownames(Sigma_y) <- c("y1star_log", "y2star_log")

cat("\nJoint covariance matrix for log endpoints (y*_SVL, y*_TL):\n")
print(Sigma_y)

# Chi-square radius for joint (1 - alpha) CI in 2D
c2    <- qchisq(ci_level, df = 2)
theta <- seq(0, 2 * pi, length.out = 400)

# Cholesky factor (Sigma_y = C^T C)
C_sy  <- chol(Sigma_y)  # upper triangular

ellipse_log <- sapply(theta, function(th) {
  v <- c(cos(th), sin(th))
  # boundary point: mu + sqrt(c2) * C^T v
  mu_y + sqrt(c2) * (t(C_sy) %*% v)
})
ellipse_log <- as.data.frame(t(ellipse_log))
colnames(ellipse_log) <- c("y1_log", "y2_log")

pt_hat_log <- data.frame(
  y1_log = ystar1_hat,
  y2_log = ystar2_hat
)

p_joint_log <- ggplot() +
  geom_path(data = ellipse_log,
            aes(x = y1_log, y = y2_log)) +
  geom_point(data = pt_hat_log,
             aes(x = y1_log, y = y2_log),
             size = 2, color = "red") +
  labs(
    x = "log endpoint SVL*",
    y = "log endpoint TL*",
    title = glue("{round(100 * ci_level)}% joint CI for log endpoints (SVL*, TL*)")
  ) +
  theme_science_polished
p_joint_log
ggsave(file.path(FIG_DIR, "joint_endpoint_log_ellipse.png"),
       p_joint_log, dpi = 600, w = 5.5, h = 5.0, units = "in")

# --- 7.6 Joint region on original scale (SVL*, TL*) -------------------
ellipse_orig <- ellipse_log %>%
  mutate(
    SVL_star = exp(y1_log),
    TL_star  = exp(y2_log)
  )

pt_hat_orig <- data.frame(
  SVL_star = ystar1_hat_orig,
  TL_star  = ystar2_hat_orig
)

p_joint_orig <- ggplot() +
  geom_path(data = ellipse_orig,
            aes(x = SVL_star, y = TL_star)) +
  geom_point(data = pt_hat_orig,
             aes(x = SVL_star, y = TL_star),
             size = 2, color = "red") +
  labs(
    x = "Endpoint SVL*",
    y = "Endpoint TL*",
    title = glue("{round(100 * ci_level)}% joint CI for endpoints (SVL*, TL*)")
  ) +
  theme_science_polished
p_joint_orig
ggsave(file.path(FIG_DIR, "joint_endpoint_orig_ellipse.png"),
       p_joint_orig, dpi = 600, w = 5.5, h = 5.0, units = "in")

# ============================================================
# 8. Endpoint estimation under unconstrained bivariate fit
#    (separate xi1, xi2) + marginal and joint CIs
# ============================================================

# --- 8.1 Extract free-parameter estimates and covariance ---------------
par_free <- fit_bv_free$estimate
V_free   <- fit_bv_free$var.cov

cat("\nUnconstrained bivariate parameter estimates (fbvpot, separate shapes):\n")
print(par_free)
cat("\nCovariance matrix of unconstrained estimates:\n")
print(V_free)

# Marginal shapes and scales from unconstrained fit
sigma1_hat_free <- unname(par_free["scale1"])  # SVL scale
xi1_hat_free    <- unname(par_free["shape1"])  # SVL shape

sigma2_hat_free <- unname(par_free["scale2"])  # TL scale
xi2_hat_free    <- unname(par_free["shape2"])  # TL shape

# --- 8.2 Point estimates of endpoints on log scale (free model) -------
# GP endpoint for margin j: y*_j = u_j - sigma_j / xi_j
ystar1_hat_free <- as.numeric(u1 - sigma1_hat_free / xi1_hat_free)  # SVL endpoint (log)
ystar2_hat_free <- as.numeric(u2 - sigma2_hat_free / xi2_hat_free)  # TL  endpoint (log)

cat("\nLog-scale endpoint estimates (unconstrained GP):\n")
cat("  y*_SVL_free (log) =", ystar1_hat_free, "\n")
cat("  y*_TL_free  (log) =", ystar2_hat_free, "\n")

# --- 8.3 Delta-method gradients for y*_SVL_free, y*_TL_free ----------
param_names_free <- names(par_free)
k_free           <- length(par_free)

idx_scale1_free  <- which(param_names_free == "scale1")
idx_shape1_free  <- which(param_names_free == "shape1")
idx_scale2_free  <- which(param_names_free == "scale2")
idx_shape2_free  <- which(param_names_free == "shape2")
# dep parameter again has zero gradient

# y*_SVL_free = u1 - scale1 / xi1
g1_free <- rep(0, k_free)
g1_free[idx_scale1_free] <- -1 / xi1_hat_free
g1_free[idx_shape1_free] <-  sigma1_hat_free / (xi1_hat_free^2)

# y*_TL_free = u2 - scale2 / xi2
g2_free <- rep(0, k_free)
g2_free[idx_scale2_free] <- -1 / xi2_hat_free
g2_free[idx_shape2_free] <-  sigma2_hat_free / (xi2_hat_free^2)

# --- 8.4 Univariate delta-method SEs and marginal CIs -----------------
var_y1_free <- as.numeric(t(g1_free) %*% V_free %*% g1_free)
var_y2_free <- as.numeric(t(g2_free) %*% V_free %*% g2_free)

se_y1_free  <- sqrt(var_y1_free)
se_y2_free  <- sqrt(var_y2_free)

zcrit <- qnorm(1 - (1 - ci_level) / 2)

ci_y1_log_free <- ystar1_hat_free + c(-1, 1) * zcrit * se_y1_free
ci_y2_log_free <- ystar2_hat_free + c(-1, 1) * zcrit * se_y2_free

cat("\n", 100 * ci_level,
    "% marginal CIs for log endpoints (unconstrained):\n", sep = "")
cat("  y*_SVL_free (log) CI: [",
    ci_y1_log_free[1], ", ", ci_y1_log_free[2], "]\n", sep = "")
cat("  y*_TL_free  (log) CI: [",
    ci_y2_log_free[1], ", ", ci_y2_log_free[2], "]\n", sep = "")

# Transform to original scale
ystar1_hat_free_orig <- exp(ystar1_hat_free)
ystar2_hat_free_orig <- exp(ystar2_hat_free)

ci_y1_orig_free <- exp(ci_y1_log_free)
ci_y2_orig_free <- exp(ci_y2_log_free)

cat("\n", 100 * ci_level,
    "% marginal CIs for endpoints (original scale, unconstrained):\n", sep = "")
cat("  SVL*_free =", ystar1_hat_free_orig,
    "  CI: [", ci_y1_orig_free[1], ", ", ci_y1_orig_free[2], "]\n", sep = "")
cat("  TL*_free  =", ystar2_hat_free_orig,
    "  CI: [", ci_y2_orig_free[1], ", ", ci_y2_orig_free[2], "]\n", sep = "")

# --- 8.5 Joint (bivariate) CI on log scale (unconstrained) ------------
G_free       <- rbind(g1_free, g2_free)          # 2 x k
Sigma_y_free <- G_free %*% V_free %*% t(G_free)  # 2 x 2 covariance
mu_y_free    <- c(ystar1_hat_free, ystar2_hat_free)

colnames(Sigma_y_free) <- rownames(Sigma_y_free) <-
  c("y1star_log_free", "y2star_log_free")

cat("\nJoint covariance matrix for log endpoints (unconstrained):\n")
print(Sigma_y_free)

c2_free <- qchisq(ci_level, df = 2)
theta   <- seq(0, 2 * pi, length.out = 400)

C_sy_free <- chol(Sigma_y_free)  # upper triangular

ellipse_log_free <- sapply(theta, function(th) {
  v <- c(cos(th), sin(th))
  mu_y_free + sqrt(c2_free) * (t(C_sy_free) %*% v)
})
ellipse_log_free <- as.data.frame(t(ellipse_log_free))
colnames(ellipse_log_free) <- c("y1_log", "y2_log")

pt_hat_log_free <- data.frame(
  y1_log = ystar1_hat_free,
  y2_log = ystar2_hat_free
)

p_joint_log_free <- ggplot() +
  geom_path(data = ellipse_log_free,
            aes(x = y1_log, y = y2_log)) +
  geom_point(data = pt_hat_log_free,
             aes(x = y1_log, y = y2_log),
             size = 2, color = "red") +
  labs(
    x = "log endpoint SVL* (free)",
    y = "log endpoint TL* (free)",
    title = glue("{round(100 * ci_level)}% joint CI (log endpoints, unconstrained)")
  ) +
  theme_science_polished

ggsave(file.path(FIG_DIR, "joint_endpoint_log_ellipse_free.png"),
       p_joint_log_free, dpi = 600, w = 5.5, h = 5.0, units = "in")

# --- 8.6 Joint region on original scale (SVL*, TL*, unconstrained) ----
ellipse_orig_free <- ellipse_log_free %>%
  mutate(
    SVL_star_free = exp(y1_log),
    TL_star_free  = exp(y2_log)
  )

pt_hat_orig_free <- data.frame(
  SVL_star_free = ystar1_hat_free_orig,
  TL_star_free  = ystar2_hat_free_orig
)

p_joint_orig_free <- ggplot() +
  geom_path(data = ellipse_orig_free,
            aes(x = SVL_star_free, y = TL_star_free)) +
  geom_point(data = pt_hat_orig_free,
             aes(x = SVL_star_free, y = TL_star_free),
             size = 2, color = "red") +
  labs(
    x = "Endpoint SVL* (free)",
    y = "Endpoint TL* (free)",
    title = glue("{round(100 * ci_level)}% joint CI (endpoints, unconstrained)")
  ) +
  theme_science_polished

ggsave(file.path(FIG_DIR, "joint_endpoint_orig_ellipse_free.png"),
       p_joint_orig_free, dpi = 600, w = 5.5, h = 5.0, units = "in")

