# ============================================================
# Bivariate EVT with censored likelihood & logistic dependence
# (SVL, TL) — American alligators, endpoint parametrization
# ============================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(ggplot2)
  library(extRemes)   # univariate GP fits + diagnostics
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
thresh_q_opt  <- c(SVL = 0.95, TL = 0.95)  # anchor quantiles (log scale)

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
print(exp(u0_by_trait))  # original scale

# ============================================================
# 4. Optional univariate POT diagnostics
# ============================================================

make_thresholds <- function(y, q_from = u_scan_lo_q, q_to = u_scan_hi_q, n = u_scan_n) {
  rng <- quantile(y, c(q_from, q_to), na.rm = TRUE)
  seq(rng[1], rng[2], length.out = n)
}

fit_gpd_at_u <- function(y, u) {
  fit <- fevd(y, type = "GP", threshold = u)
  par <- fit$results$par
  cov <- summary(fit)$cov.theta
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
  p_wald <- pnorm(z)
  
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
# 5. Bivariate censored likelihood — endpoint parametrization
#    theta = (y1_star, xi1, y2_star, xi2, alpha)
# ============================================================

df_biv <- df_log %>%
  filter(!is.na(log_SVL), !is.na(log_TL)) %>%
  transmute(
    specimen,
    x1 = log_SVL,
    x2 = log_TL
  )

n_biv <- nrow(df_biv)
cat("Bivariate sample size (both SVL & TL observed):", n_biv, "\n")

x1_biv <- df_biv$x1
x2_biv <- df_biv$x2

u1 <- u0_by_trait["SVL"]
u2 <- u0_by_trait["TL"]

lambda1 <- mean(x1_biv > u1)
lambda2 <- mean(x2_biv > u2)

cat("Empirical exceedance probs: lambda1 =", lambda1,
    ", lambda2 =", lambda2, "\n")

# bases to ensure y* > max(x, u)
base1 <- max(c(u1, x1_biv), na.rm = TRUE)
base2 <- max(c(u2, x2_biv), na.rm = TRUE)

# ============================================================
# 5.1 Endpoint → σ helper: σ = (u - y*) ξ (ξ<0, y*>u ⇒ σ>0)
# ============================================================

sigma_from_endpoint <- function(y_star, u, xi) {
  (u - y_star) * xi
}

# ============================================================
# 5.2 Logistic stable tail function and derivatives
#      l(v1,v2; alpha) = (v1^{1/alpha} + v2^{1/alpha})^alpha
#      with 0 < alpha ≤ 1 (alpha = 1 ↔ independence)
# ============================================================

logistic_l <- function(v1, v2, alpha) {
  S <- v1^(1/alpha) + v2^(1/alpha)
  S^alpha
}

logistic_l_grad <- function(v1, v2, alpha) {
  S  <- v1^(1/alpha) + v2^(1/alpha)
  l1 <- S^(alpha - 1) * v1^(1/alpha - 1)
  l2 <- S^(alpha - 1) * v2^(1/alpha - 1)
  list(l1 = l1, l2 = l2)
}

logistic_l_hess12 <- function(v1, v2, alpha) {
  S <- v1^(1/alpha) + v2^(1/alpha)
  (alpha - 1) * S^(alpha - 2) * (1/alpha) *
    v2^(1/alpha - 1) * v1^(1/alpha - 1)
}

# ============================================================
# 5.3 Marginal tail transforms v_j(x_j) and dv_j/dx_j via (y*, xi)
#      v_j(x) = λ_j [1 + ξ_j (x - u_j)/σ_j]^{-1/ξ_j}
# ============================================================

v_fun_endpoint <- function(x, u, y_star, xi, lambda) {
  sigma <- sigma_from_endpoint(y_star, u, xi)
  A     <- 1 + xi * (x - u) / sigma
  v     <- lambda * A^(-1/xi)
  v[!is.finite(v) | A <= 0] <- NaN
  v
}

dv_dx_endpoint <- function(x, u, y_star, xi, lambda) {
  sigma <- sigma_from_endpoint(y_star, u, xi)
  A     <- 1 + xi * (x - u) / sigma
  v     <- lambda * A^(-1/xi)
  dv    <- - v / (sigma * A)
  dv[!is.finite(dv) | A <= 0] <- NaN
  dv
}

# ============================================================
# 5.4 Censored likelihood with logistic dependence
# ============================================================

loglik_censored_logistic_endpoint <- function(theta,
                                              x1, x2,
                                              u1, u2,
                                              lambda1, lambda2) {
  
  if (any(!is.finite(theta))) return(-1e10)
  
  y1_star <- theta[1]
  xi1     <- theta[2]
  y2_star <- theta[3]
  xi2     <- theta[4]
  alpha   <- theta[5]
  
  # ξ_j < 0, 0 < α ≤ 1
  if (xi1 >= 0 || xi2 >= 0) return(-1e10)
  if (alpha <= 0 || alpha > 1) return(-1e10)
  
  sigma1 <- sigma_from_endpoint(y1_star, u1, xi1)
  sigma2 <- sigma_from_endpoint(y2_star, u2, xi2)
  if (!is.finite(sigma1) || !is.finite(sigma2) ||
      sigma1 <= 0 || sigma2 <= 0) {
    return(-1e10)
  }
  
  n <- length(x1)
  if (length(x2) != n) return(-1e10)
  
  R00 <- (x1 <= u1 & x2 <= u2)
  R10 <- (x1 >  u1 & x2 <= u2)
  R01 <- (x1 <= u1 & x2 >  u2)
  R11 <- (x1 >  u1 & x2 >  u2)
  
  # v_j at thresholds
  v1_u <- v_fun_endpoint(u1, u1, y1_star, xi1, lambda1)
  v2_u <- v_fun_endpoint(u2, u2, y2_star, xi2, lambda2)
  if (any(!is.finite(c(v1_u, v2_u)))) return(-1e10)
  
  F_u1u2 <- exp(-logistic_l(v1_u, v2_u, alpha))
  
  logL <- rep(NA_real_, n)
  
  # R00: both below thresholds
  if (any(R00)) {
    logL[R00] <- log(F_u1u2)
  }
  
  # R10: x1 > u1, x2 <= u2
  if (any(R10)) {
    x1_10 <- x1[R10]
    v1_10 <- v_fun_endpoint(x1_10, u1, y1_star, xi1, lambda1)
    v2_10 <- v2_u
    if (any(!is.finite(v1_10))) return(-1e10)
    
    grads_10 <- logistic_l_grad(v1_10, v2_10, alpha)
    l1_10    <- grads_10$l1
    F_10     <- exp(-logistic_l(v1_10, v2_10, alpha))
    dv1_10   <- dv_dx_endpoint(x1_10, u1, y1_star, xi1, lambda1)
    
    dens_10 <- - F_10 * l1_10 * dv1_10
    dens_10[!is.finite(dens_10) | dens_10 <= 0] <- .Machine$double.xmin
    
    logL[R10] <- log(dens_10)
  }
  
  # R01: x1 <= u1, x2 > u2
  if (any(R01)) {
    x2_01 <- x2[R01]
    v1_01 <- v1_u
    v2_01 <- v_fun_endpoint(x2_01, u2, y2_star, xi2, lambda2)
    if (any(!is.finite(v2_01))) return(-1e10)
    
    grads_01 <- logistic_l_grad(v1_01, v2_01, alpha)
    l2_01    <- grads_01$l2
    F_01     <- exp(-logistic_l(v1_01, v2_01, alpha))
    dv2_01   <- dv_dx_endpoint(x2_01, u2, y2_star, xi2, lambda2)
    
    dens_01 <- - F_01 * l2_01 * dv2_01
    dens_01[!is.finite(dens_01) | dens_01 <= 0] <- .Machine$double.xmin
    
    logL[R01] <- log(dens_01)
  }
  
  # R11: x1 > u1, x2 > u2
  if (any(R11)) {
    x1_11 <- x1[R11]
    x2_11 <- x2[R11]
    
    v1_11 <- v_fun_endpoint(x1_11, u1, y1_star, xi1, lambda1)
    v2_11 <- v_fun_endpoint(x2_11, u2, y2_star, xi2, lambda2)
    if (any(!is.finite(c(v1_11, v2_11)))) return(-1e10)
    
    grads_11 <- logistic_l_grad(v1_11, v2_11, alpha)
    l1_11    <- grads_11$l1
    l2_11    <- grads_11$l2
    l12_11   <- logistic_l_hess12(v1_11, v2_11, alpha)
    
    F_11   <- exp(-logistic_l(v1_11, v2_11, alpha))
    dv1_11 <- dv_dx_endpoint(x1_11, u1, y1_star, xi1, lambda1)
    dv2_11 <- dv_dx_endpoint(x2_11, u2, y2_star, xi2, lambda2)
    
    dens_11 <- F_11 * dv1_11 * dv2_11 * (l1_11 * l2_11 - l12_11)
    dens_11[!is.finite(dens_11) | dens_11 <= 0] <- .Machine$double.xmin
    
    logL[R11] <- log(dens_11)
  }
  
  if (any(!is.finite(logL))) return(-1e10)
  sum(logL)
}

eval_biv_loglik_endpoint <- function(theta) {
  loglik_censored_logistic_endpoint(
    theta   = theta,
    x1      = x1_biv,
    x2      = x2_biv,
    u1      = u1,
    u2      = u2,
    lambda1 = lambda1,
    lambda2 = lambda2
  )
}

# ============================================================
# 5.5 Reparametrize to unconstrained phi
#     phi = (eta1, zeta1, eta2, zeta2, eta_alpha)
#     y1* = base1 + exp(eta1),  xi1 = -exp(zeta1)
#     y2* = base2 + exp(eta2),  xi2 = -exp(zeta2)
#     alpha = plogis(eta_alpha) ∈ (0,1)
# ============================================================

loglik_from_phi <- function(phi) {
  eta1      <- phi[1]
  zeta1     <- phi[2]
  eta2      <- phi[3]
  zeta2     <- phi[4]
  eta_alpha <- phi[5]
  
  y1_star <- base1 + exp(eta1)
  xi1     <- -exp(zeta1)
  y2_star <- base2 + exp(eta2)
  xi2     <- -exp(zeta2)
  alpha   <- plogis(eta_alpha)  # (0,1)
  
  theta <- c(y1_star, xi1, y2_star, xi2, alpha)
  
  loglik_censored_logistic_endpoint(
    theta   = theta,
    x1      = x1_biv,
    x2      = x2_biv,
    u1      = u1,
    u2      = u2,
    lambda1 = lambda1,
    lambda2 = lambda2
  )
}

neg_loglik_phi <- function(phi) {
  -loglik_from_phi(phi)
}

# ============================================================
# 5.6 Initial values from univariate GP fits and optimization
# ============================================================

fit1_biv <- fevd(x1_biv, type = "GP", threshold = u1)
par1_biv <- fit1_biv$results$par
xi1_hat_uni  <- unname(par1_biv["shape"])
sig1_hat_uni <- unname(par1_biv["scale"])

fit2_biv <- fevd(x2_biv, type = "GP", threshold = u2)
par2_biv <- fit2_biv$results$par
xi2_hat_uni  <- unname(par2_biv["shape"])
sig2_hat_uni <- unname(par2_biv["scale"])

y1_star_0 <- u1 - sig1_hat_uni / xi1_hat_uni
y2_star_0 <- u2 - sig2_hat_uni / xi2_hat_uni
alpha_0   <- 0.7   # moderate tail dependence

phi_start <- c(
  log(y1_star_0 - base1),
  log(-xi1_hat_uni),
  log(y2_star_0 - base2),
  log(-xi2_hat_uni),
  qlogis(alpha_0)     # inverse of plogis
)

opt_phi <- optim(
  par     = phi_start,
  fn      = neg_loglik_phi,
  method  = "BFGS",
  control = list(maxit = 500)
)

phi_hat <- opt_phi$par

eta1_hat      <- phi_hat[1]
zeta1_hat     <- phi_hat[2]
eta2_hat      <- phi_hat[3]
zeta2_hat     <- phi_hat[4]
eta_alpha_hat <- phi_hat[5]

y1_star_hat <- base1 + exp(eta1_hat)
xi1_hat_biv <- -exp(zeta1_hat)
y2_star_hat <- base2 + exp(eta2_hat)
xi2_hat_biv <- -exp(zeta2_hat)
alpha_hat   <- plogis(eta_alpha_hat)

ell_hat <- -opt_phi$value

cat("Bivariate endpoint MLEs (logistic dependence):\n")
cat("  y1* =", y1_star_hat, " xi1 =", xi1_hat_biv, "\n")
cat("  y2* =", y2_star_hat, " xi2 =", xi2_hat_biv, "\n")
cat("  alpha =", alpha_hat, "\n")
cat("Maximized censored log-likelihood =", ell_hat, "\n")


# ============================================================
# 5.7 Likelihood ratio test for tail equality: H0: xi1 = xi2
# ============================================================

# Under H0 we set xi1 = xi2 = xi_common, but keep distinct endpoints y1*, y2*
# and logistic alpha. Reparametrize again to unconstrained φ_eq:
#   φ_eq = (eta1, zeta, eta2, eta_alpha)
#   y1*  = base1 + exp(eta1)
#   y2*  = base2 + exp(eta2)
#   xi   = -exp(zeta)
#   alpha = plogis(eta_alpha) ∈ (0,1)

loglik_from_phi_eq <- function(phi_eq) {
  eta1_eq    <- phi_eq[1]
  zeta_eq    <- phi_eq[2]
  eta2_eq    <- phi_eq[3]
  eta_alpha_eq <- phi_eq[4]
  
  y1_star_eq <- base1 + exp(eta1_eq)
  y2_star_eq <- base2 + exp(eta2_eq)
  xi_eq      <- -exp(zeta_eq)
  alpha_eq   <- plogis(eta_alpha_eq)
  
  theta_eq <- c(y1_star_eq, xi_eq, y2_star_eq, xi_eq, alpha_eq)
  
  loglik_censored_logistic_endpoint(
    theta   = theta_eq,
    x1      = x1_biv,
    x2      = x2_biv,
    u1      = u1,
    u2      = u2,
    lambda1 = lambda1,
    lambda2 = lambda2
  )
}

neg_loglik_phi_eq <- function(phi_eq) {
  -loglik_from_phi_eq(phi_eq)
}

# ---- Starting values for the constrained fit ----
# Use unconstrained MLEs as a sensible anchor:
xi_common_0 <- (xi1_hat_biv + xi2_hat_biv) / 2

phi_eq_start <- c(
  log(y1_star_hat - base1),   # eta1
  log(-xi_common_0),          # zeta (common shape, negative)
  log(y2_star_hat - base2),   # eta2
  qlogis(alpha_hat)           # eta_alpha
)

opt_phi_eq <- optim(
  par     = phi_eq_start,
  fn      = neg_loglik_phi_eq,
  method  = "BFGS",
  control = list(maxit = 500)
)

phi_eq_hat <- opt_phi_eq$par

eta1_eq_hat      <- phi_eq_hat[1]
zeta_eq_hat      <- phi_eq_hat[2]
eta2_eq_hat      <- phi_eq_hat[3]
eta_alpha_eq_hat <- phi_eq_hat[4]

y1_star_eq_hat <- base1 + exp(eta1_eq_hat)
y2_star_eq_hat <- base2 + exp(eta2_eq_hat)
xi_eq_hat      <- -exp(zeta_eq_hat)
alpha_eq_hat   <- plogis(eta_alpha_eq_hat)

ell_hat_eq <- -opt_phi_eq$value

cat("\nConstrained MLEs under H0: xi1 = xi2 = xi\n")
cat("  y1* =", y1_star_eq_hat, " y2* =", y2_star_eq_hat, "\n")
cat("  xi  =", xi_eq_hat, "\n")
cat("  alpha =", alpha_eq_hat, "\n")
cat("Maximized censored log-likelihood (H0) =", ell_hat_eq, "\n")

# ---- Likelihood ratio statistic ----
LR_stat <- 2 * (ell_hat - ell_hat_eq)   # H1 (free) vs H0 (xi1=xi2)
df_LR   <- 1
p_LR    <- 1 - pchisq(LR_stat, df = df_LR)

cat("\nLikelihood ratio test for tail equality H0: xi1 = xi2\n")
cat("  LR statistic =", LR_stat, " with df =", df_LR, "\n")
cat("  p-value      =", p_LR, "\n")

