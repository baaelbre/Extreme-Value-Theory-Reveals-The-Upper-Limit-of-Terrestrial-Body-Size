# ============================================================
# Title : Extreme Value Theory Reveals the Upper Limit of
#         Terrestrial Body Size: Alligators
# Author: Bastiaan A. Van Velthoven
# ============================================================

suppressPackageStartupMessages({
  library(readxl)   # read Excel
  library(dplyr)    # data manipulation
  library(purrr)    # map_df
  library(ggplot2)  # plotting
  library(evd)      # EVT (fpot, pgpd, qgpd, dgpd)
  library(scales)   # percent labels
  library(glue)     # glue strings
  library(forcats)  # factor ordering
  library(grid)     # unit()
  library(MASS)     # kde2d
  library(copula)   # Gumbel (logistic EV) copula
})

set.seed(42)

# ------------------------------------------------------------
# Load compiled C code (from evd's bvpot.c): nllbvclog()
# ------------------------------------------------------------
# citation("evd")
# system("R CMD SHLIB bvpot.c") # (compile if needed; produces bvpot{.dll/.so/.dylib})

dyn.load(paste0("bvpot", .Platform$dynlib.ext))
stopifnot(is.loaded("nllbvclog"))

# ---------------------------
# Directories & theming
# ---------------------------
FIG_DIR <- "Figures/Alligators"
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

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

save_png <- function(p, name, w = 6.2, h = 4.2, dpi = 600) {
  ggsave(file.path(FIG_DIR, name), p, dpi = dpi, width = w, height = h, units = "in")
}

# ---------------------------
# Global settings
# ---------------------------
ci_level <- 0.90    # confidence level
u_lo     <- 0.60    # threshold scan: start quantile
u_hi     <- 0.99    # threshold scan: end quantile
u_n      <- 50L     # number of thresholds to scan
min_ex   <- 20L     # minimum exceedances for stable ML fits

trait_names <- c("SVL", "TL")               # SVL: snout-vent length, TL: total length
q_opt       <- c(SVL = 0.94, TL = 0.94)     # anchor thresholds (upper ~6%)

# Bootstrap settings
B_boot <- 1000L     # parametric bivariate bootstrap replicates
B_KS   <- 1000L     # parametric KS bootstrap per margin

# ============================================================
# 1) Read the data
# ============================================================
DATA_XLSX <- "Data/alligators_woodward.xlsx"
df_raw    <- read_excel(DATA_XLSX)

# Deform == 1 or 3: tail broken ⇒ TL structurally missing; SVL kept
df <- df_raw %>%
  mutate(
    SVL = as.numeric(SVL),
    TL  = as.numeric(TL),
    TL  = ifelse(Deform %in% c(1, 3), NA_real_, TL)
  ) %>%
  transmute(specimen = row_number(), SVL, TL)

# ============================================================
# 2) Completeness diagnostics
# ============================================================
compl_tbl <- map_df(trait_names, function(tr) {
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
  labs(y = "Completeness") +
  theme_science

save_png(p_compl, "completeness_SVL_TL.png", w = 6.0, h = 4.2)
print(p_compl) # SVL rarely missing; TL missing in a small fraction

# ============================================================
# 3) Log-transform + split missingness patterns
# ============================================================
df_log <- df %>%
  mutate(
    log_SVL = ifelse(is.finite(SVL) & SVL > 0, log(SVL), NA_real_),
    log_TL  = ifelse(is.finite(TL)  & TL  > 0, log(TL),  NA_real_)
  )

df_complete <- df_log %>% filter(!is.na(log_SVL), !is.na(log_TL))
df_svl_only <- df_log %>% filter(!is.na(log_SVL),  is.na(log_TL))
df_tl_only  <- df_log %>% filter( is.na(log_SVL), !is.na(log_TL))

cat("\nRow counts by observation pattern:\n")
cat("  complete (SVL & TL):", nrow(df_complete), "\n")
cat("  SVL only           :", nrow(df_svl_only), "\n")
cat("  TL only            :", nrow(df_tl_only), "\n")

# ------------------------------------------------------------
# Thresholds per trait (log scale)
# ------------------------------------------------------------
u0_by_trait <- setNames(numeric(length(trait_names)), trait_names)
for (tr in trait_names) {
  v_log <- df_log[[paste0("log_", tr)]]
  u0_by_trait[tr] <- as.numeric(quantile(v_log, q_opt[tr], na.rm = TRUE))
}

u1    <- unname(u0_by_trait["SVL"])
u2    <- unname(u0_by_trait["TL"])
u_vec <- c(u1, u2)

cat("\nThresholds (log scale):\n"); print(u0_by_trait)
cat("\nThresholds (original scale):\n"); print(exp(u0_by_trait))

# ============================================================
# 4) Univariate diagnostics helpers
# ============================================================

# scan thresholds between quantiles u_lo..u_hi
make_thresholds <- function(y, q_from = u_lo, q_to = u_hi, n = u_n) {
  r <- quantile(y, c(q_from, q_to), na.rm = TRUE)
  seq(r[1], r[2], length.out = n)
}

fit_gpd_at_u <- function(y, u) {
  fit <- evd::fpot(y, threshold = u, model = "gpd")
  par <- fit$estimate
  list(
    scale    = unname(par["scale"]),
    shape    = unname(par["shape"]),
    cov      = fit$var.cov,
    n_exceed = sum(y > u)
  )
}

# adjusted (threshold-invariant) scale: sigma_u - xi(u-u0)  (Coles Ch.3)
adj_scale_fun <- function(scale, xi, u, u0) scale - xi * (u - u0)

diagnostic_plots <- function(y, u, sigma, xi, label) {
  exc   <- y[y > u] - u
  n     <- length(exc)
  probs <- ppoints(n)
  
  theo_q <- if (abs(xi) > 1e-10) {
    u + sigma / xi * (probs^(-xi) - 1)
  } else {
    u - sigma * log(probs) # exponential (xi -> 0) limit
  }
  
  dfqq <- tibble(Theoretical = rev(theo_q), Empirical = sort(y[y > u]))
  pqq  <- ggplot(dfqq, aes(Theoretical, Empirical)) +
    geom_point(color = "steelblue") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    labs(title = glue("Q–Q: {label}")) +
    theme_science
  
  F_theo <- if (abs(xi) > 1e-10) {
    1 - (1 + xi * exc / sigma)^(-1 / xi)
  } else {
    1 - exp(-exc / sigma)
  }
  dfpp <- tibble(Theoretical = sort(F_theo), Empirical = (1:n) / n)
  ppp  <- ggplot(dfpp, aes(Theoretical, Empirical)) +
    geom_point(color = "darkgreen") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    labs(title = glue("P–P: {label}")) +
    theme_science
  
  list(pqq = pqq, ppp = ppp)
}

# ------------------------------------------------------------
# KS goodness-of-fit via parametric bootstrap (univariate GPD)
# ------------------------------------------------------------
ks_boot <- function(y, u, scale_hat, shape_hat, B = 1000L, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  ex <- y[y > u] - u
  n_ex <- length(ex)
  if (n_ex < 20L) {
    warning("Not enough exceedances for KS bootstrap (n_ex < 20) -> NA")
    return(list(D = NA_real_, p = NA_real_, D_boot = numeric(0L), n_exceed = n_ex))
  }
  
  # observed KS statistic under fitted GPD
  Fn   <- ecdf(ex)
  exs  <- sort(ex)
  Fgpd <- evd::pgpd(exs, loc = 0, scale = scale_hat, shape = shape_hat)
  Dobs <- max(abs(Fn(exs) - Fgpd))
  
  # parametric bootstrap with refit
  Dboot <- rep(NA_real_, B)
  for (b in seq_len(B)) {
    ex_b <- evd::rgpd(n_ex, loc = 0, scale = scale_hat, shape = shape_hat)
    y_b  <- u + ex_b
    
    fit_b <- try(evd::fpot(y_b, threshold = u, model = "gpd"), silent = TRUE)
    if (inherits(fit_b, "try-error")) next
    
    sc_b <- unname(fit_b$estimate["scale"])
    xi_b <- unname(fit_b$estimate["shape"])
    
    Fn_b   <- ecdf(ex_b)
    ex_bs  <- sort(ex_b)
    Fgpd_b <- evd::pgpd(ex_bs, loc = 0, scale = sc_b, shape = xi_b)
    Dboot[b] <- max(abs(Fn_b(ex_bs) - Fgpd_b))
  }
  
  Dboot <- Dboot[is.finite(Dboot)]
  if (!length(Dboot)) {
    warning("All KS bootstrap replicates failed; returning NA p-value.")
    return(list(D = Dobs, p = NA_real_, D_boot = numeric(0L), n_exceed = n_ex))
  }
  
  list(D = Dobs, p = mean(Dboot >= Dobs), D_boot = Dboot, n_exceed = n_ex)
}

# ------------------------------------------------------------
# Run univariate diagnostics per margin (xi stability + KS + QQ/PP)
# ------------------------------------------------------------
run_trait <- function(tr_key, u0) {
  y <- df_log[[paste0("log_", tr_key)]]
  y <- y[is.finite(y)]
  stopifnot(length(y) >= 40)
  
  u_seq <- make_thresholds(y)
  fit0  <- fit_gpd_at_u(y, u0)
  
  xi0   <- fit0$shape
  sc0   <- fit0$scale
  cov0  <- fit0$cov
  xi_se <- sqrt(cov0[2, 2])
  sc_se <- sqrt(cov0[1, 1])
  
  # Wald test for xi < 0 (one-sided)
  z_wald <- xi0 / xi_se
  p_wald <- pnorm(z_wald)  # H1: xi < 0
  
  # KS goodness-of-fit (parametric bootstrap)
  ks_res <- ks_boot(y, u0, sc0, xi0, B = B_KS, seed = if (tr_key == "SVL") 101L else 202L)
  
  # threshold stability (xi, adjusted scale)
  zcrit <- qnorm(1 - (1 - ci_level) / 2)
  
  df_scan <- tibble(
    u = u_seq,
    shape = NA_real_, shape_lo = NA_real_, shape_hi = NA_real_,
    adj   = NA_real_, adj_lo   = NA_real_, adj_hi   = NA_real_
  )
  
  for (i in seq_along(u_seq)) {
    u <- u_seq[i]
    if (sum(y > u) < min_ex) next
    
    out <- try(fit_gpd_at_u(y, u), silent = TRUE)
    if (inherits(out, "try-error")) next
    
    se_xi <- sqrt(out$cov[2, 2])
    df_scan$shape[i]    <- out$shape
    df_scan$shape_lo[i] <- out$shape - zcrit * se_xi
    df_scan$shape_hi[i] <- out$shape + zcrit * se_xi
    
    df_scan$adj[i] <- adj_scale_fun(out$scale, out$shape, u, u0)
    
    # delta-method var for adj = scale - (u-u0)*shape
    var_adj <- out$cov[1, 1] + (u - u0)^2 * out$cov[2, 2] - 2 * (u - u0) * out$cov[1, 2]
    se_adj  <- sqrt(max(var_adj, 0))
    df_scan$adj_lo[i] <- df_scan$adj[i] - zcrit * se_adj
    df_scan$adj_hi[i] <- df_scan$adj[i] + zcrit * se_adj
  }
  
  ylab <- glue("log {tr_key}")
  
  p_shape <- ggplot(df_scan, aes(u, shape)) +
    geom_point(color = "blue") +
    geom_errorbar(aes(ymin = shape_lo, ymax = shape_hi), width = 0.03, color = "blue") +
    geom_vline(xintercept = u0, color = "red", linetype = "dashed") +
    geom_hline(yintercept = xi0, color = "red", linetype = "dashed") +
    labs(x = glue("Threshold ({ylab})"), y = "Shape (xi)") +
    theme_science
  
  p_adj <- ggplot(df_scan, aes(u, adj)) +
    geom_point(color = "blue") +
    geom_errorbar(aes(ymin = adj_lo, ymax = adj_hi), width = 0.03, color = "blue") +
    geom_vline(xintercept = u0, color = "red", linetype = "dashed") +
    geom_hline(yintercept = sc0, color = "red", linetype = "dashed") +
    labs(x = glue("Threshold ({ylab})"), y = "Adjusted scale") +
    theme_science
  
  save_png(p_shape, glue("{tr_key}_log_shape_stability.png"))
  save_png(p_adj,   glue("{tr_key}_log_adj_scale_stability.png"))
  
  di <- diagnostic_plots(y, u0, sc0, xi0, ylab)
  save_png(di$pqq, glue("{tr_key}_log_qq.png"))
  save_png(di$ppp, glue("{tr_key}_log_pp.png"))
  
  cat(glue("\n{tr_key}: u0={round(u0,3)} | xî={round(xi0,3)} (se={round(xi_se,3)}) | Wald xi<0: z={round(z_wald,3)}, p={signif(p_wald,3)}"))
  cat(glue("\n{tr_key}: KS bootstrap: D={round(ks_res$D,3)}, p_boot={signif(ks_res$p,3)}, n_exceed={ks_res$n_exceed}\n"))
  
  tibble(
    trait = tr_key,
    u0 = u0,
    sigma_hat = sc0,
    xi_hat    = xi0,
    sigma_se = sc_se,
    xi_se    = xi_se,
    wald_z_xi     = z_wald,
    wald_p_xi_neg = p_wald,
    ks_D       = ks_res$D,
    ks_p_boot  = ks_res$p,
    ks_n_exceed= ks_res$n_exceed
  )
}

diag_SVL <- run_trait("SVL", u1)
diag_TL  <- run_trait("TL",  u2)
marg_diag_tbl <- bind_rows(diag_SVL, diag_TL)
print(marg_diag_tbl)

sigma1_hat_uni <- diag_SVL$sigma_hat
xi1_hat_uni    <- diag_SVL$xi_hat
sigma2_hat_uni <- diag_TL$sigma_hat
xi2_hat_uni    <- diag_TL$xi_hat

# ============================================================
# 5) Bivariate POT (censored likelihood) + missing-partner GP terms
# ============================================================

# Complete-case matrix (log scale)
X_biv <- df_complete %>% transmute(log_SVL, log_TL) %>% as.matrix()
cat("\nBivariate sample size (complete cases):", nrow(X_biv), "\n")

quad_tab <- table(SVL_exceed = X_biv[, 1] > u1, TL_exceed = X_biv[, 2] > u2)
cat("\nQuadrant table (log scale, complete cases only):\n")
print(quad_tab)

# Missing patterns: exceedance with partner missing (not handled in nllbvclog)
y1_E_M2 <- df_svl_only$log_SVL[df_svl_only$log_SVL > u1]  # SVL exceeds, TL missing
y2_M1_E <- df_tl_only$log_TL[df_tl_only$log_TL > u2]      # TL exceeds, SVL missing

cat("\nCounts of missing-partner exceedances:\n")
cat("  SVL exceed, TL missing (E1M2):", length(y1_E_M2), "\n")
cat("  TL exceed, SVL missing (M1E2):", length(y2_M1_E), "\n")

# ------------------------------------------------------------
# Censored bivariate POT log-likelihood (via compiled C)
# ------------------------------------------------------------
sep.bvdata <- getFromNamespace("sep.bvdata", "evd")

make_ll_cens <- function(x, u) {
  spx <- sep.bvdata(x = x, method = "cpot", u = u)
  
  x1     <- as.double(spx$x1)
  x2     <- as.double(spx$x2)
  nn     <- as.integer(spx$nn)
  n      <- as.integer(spx$n)
  thdi   <- as.double(spx$thdi)
  lambda <- as.double(spx$lambda)
  
  function(theta) {
    out <- .C("nllbvclog",
              data1  = x1,
              data2  = x2,
              nn     = nn,
              n      = n,
              thid   = thdi,
              lambda = lambda,
              dep    = as.double(theta["dep"]),
              scale1 = as.double(theta["scale1"]),
              shape1 = as.double(theta["shape1"]),
              scale2 = as.double(theta["scale2"]),
              shape2 = as.double(theta["shape2"]),
              dns    = as.double(0))
    
    # C returns NEGATIVE log-likelihood in dns -> negate to get loglik
    -out$dns
  }
}

ll_cens_log <- make_ll_cens(X_biv, u_vec)

# ------------------------------------------------------------
# Univariate GP contributions from missing partners
# ------------------------------------------------------------
lgp_exceed <- function(y, u, scale, shape) {
  ex <- y - u
  evd::dgpd(ex, loc = 0, scale = scale, shape = shape, log = TRUE)
}

ll_missing <- function(theta) {
  ll1 <- if (length(y1_E_M2)) sum(lgp_exceed(y1_E_M2, u1, theta["scale1"], theta["shape1"])) else 0
  ll2 <- if (length(y2_M1_E)) sum(lgp_exceed(y2_M1_E, u2, theta["scale2"], theta["shape2"])) else 0
  ll1 + ll2
}

# Full log-likelihood: censored bivariate POT + univariate missings
loglik_full <- function(theta) ll_cens_log(theta) + ll_missing(theta)

# Wrapper for optim
negloglik_full <- function(p, common_shape = FALSE) {
  theta <- if (!common_shape) {
    c(scale1 = p[1], shape1 = p[2], scale2 = p[3], shape2 = p[4], dep = p[5])
  } else {
    c(scale1 = p[1], shape1 = p[2], scale2 = p[3], shape2 = p[2], dep = p[4]) # H0: xi1=xi2
  }
  -loglik_full(theta)
}

# ============================================================
# 6) Maximum likelihood: free shapes vs common shape
# ============================================================

# Free tail (xi1, xi2 free)
start_free <- c(scale1 = sigma1_hat_uni, shape1 = xi1_hat_uni,
                scale2 = sigma2_hat_uni, shape2 = xi2_hat_uni,
                dep = 0.8)

opt_free <- optim(unname(start_free), negloglik_full, common_shape = FALSE, method = "Nelder-Mead")
par_free <- setNames(opt_free$par, names(start_free))
ell_hat_free <- -opt_free$value

cat("\nBivariate MLEs (free shapes):\n"); print(par_free)

# Common tail (H0: xi1=xi2=xi)
start_eq <- c(scale1 = sigma1_hat_uni,
              xi     = mean(c(xi1_hat_uni, xi2_hat_uni)),
              scale2 = sigma2_hat_uni,
              dep    = 0.8)

opt_eq <- optim(unname(start_eq), negloglik_full, common_shape = TRUE, method = "Nelder-Mead")
par_eq_phi <- setNames(opt_eq$par, names(start_eq))
ell_hat_eq <- -opt_eq$value

par_eq <- c(scale1 = par_eq_phi["scale1"],
            shape1 = par_eq_phi["xi"],
            scale2 = par_eq_phi["scale2"],
            shape2 = par_eq_phi["xi"],
            dep    = par_eq_phi["dep"])

cat("\nCommon tail shape MLEs (H0: xi1=xi2):\n"); print(par_eq)

# Likelihood ratio test for H0: xi1 = xi2
LR_shape <- 2 * (ell_hat_free - ell_hat_eq)
p_LR     <- 1 - pchisq(LR_shape, df = 1)

cat("\nLikelihood ratio test H0: xi1 = xi2\n")
cat("  LR statistic =", LR_shape, "\n")
cat("  p-value      =", p_LR, "\n")

xi_eq_hat  <- unname(par_eq_phi["xi"])
sigma1_hat <- unname(par_eq_phi["scale1"])
sigma2_hat <- unname(par_eq_phi["scale2"])

# ------------------------------------------------------------
# Empirical tail fractions (used in bootstrap + Stokes)
# ------------------------------------------------------------
n_SVL_total <- sum(!is.na(df_log$log_SVL))
n_TL_total  <- sum(!is.na(df_log$log_TL))
tail_frac_SVL <- sum(df_log$log_SVL > u1, na.rm = TRUE) / n_SVL_total
tail_frac_TL  <- sum(df_log$log_TL  > u2, na.rm = TRUE) / n_TL_total

cat("\nEmpirical tail fractions:\n")
cat("  P(SVL > u1) ≈", tail_frac_SVL, "(u1 =", u1, ", exp(u1) ≈", exp(u1), "cm)\n")
cat("  P(TL  > u2) ≈", tail_frac_TL,  "(u2 =", u2, ", exp(u2) ≈", exp(u2), "cm)\n")

# ------------------------------------------------------------
# Parametric bivariate bootstrap (logistic GP / common tail shape)
# ------------------------------------------------------------
sub_SVL_cc   <- df_complete$log_SVL[df_complete$log_SVL <= u1]
sub_TL_cc    <- df_complete$log_TL [df_complete$log_TL  <= u2]
sub_SVL_only <- df_svl_only$log_SVL[df_svl_only$log_SVL <= u1]
sub_TL_only  <- df_tl_only$log_TL  [df_tl_only$log_TL   <= u2]

n_cc       <- nrow(df_complete)
n_svl_only <- nrow(df_svl_only)
n_tl_only  <- nrow(df_tl_only)

# Logistic EV dependence parameter (alpha) from constrained fit; map to Gumbel theta = 1/alpha
alpha_logistic_hat <- unname(par_eq["dep.dep"])
theta_gumbel_hat   <- 1 / alpha_logistic_hat
gumbel_cop_hat      <- gumbelCopula(param = theta_gumbel_hat, dim = 2)

boot_sigma1 <- boot_sigma2 <- boot_xi <- boot_y1star <- boot_y2star <- rep(NA_real_, B_boot)
boot_conv   <- rep(FALSE, B_boot)

set.seed(2026)

for (b in seq_len(B_boot)) {
  
  # (1) Simulate complete cases via copula + marginal mixture (subthreshold resample / GPD tail)
  U_cc <- rCopula(n_cc, gumbel_cop_hat)
  Y1_cc <- Y2_cc <- numeric(n_cc)
  
  for (k in seq_len(n_cc)) {
    u1k <- U_cc[k, 1]
    u2k <- U_cc[k, 2]
    
    # SVL margin
    if (u1k <= 1 - tail_frac_SVL || !length(sub_SVL_cc)) {
      Y1_cc[k] <- if (length(sub_SVL_cc)) sample(sub_SVL_cc, 1L, TRUE) else u1
    } else {
      v1 <- (u1k - (1 - tail_frac_SVL)) / tail_frac_SVL
      Y1_cc[k] <- u1 + evd::qgpd(v1, loc = 0, scale = sigma1_hat, shape = xi_eq_hat)
    }
    
    # TL margin
    if (u2k <= 1 - tail_frac_TL || !length(sub_TL_cc)) {
      Y2_cc[k] <- if (length(sub_TL_cc)) sample(sub_TL_cc, 1L, TRUE) else u2
    } else {
      v2 <- (u2k - (1 - tail_frac_TL)) / tail_frac_TL
      Y2_cc[k] <- u2 + evd::qgpd(v2, loc = 0, scale = sigma2_hat, shape = xi_eq_hat)
    }
  }
  
  X_biv_boot <- cbind(Y1_cc, Y2_cc)
  
  # (2) Simulate SVL-only and TL-only samples (same marginal mixture)
  if (n_svl_only > 0L) {
    U1 <- runif(n_svl_only)
    Y1 <- numeric(n_svl_only)
    for (k in seq_len(n_svl_only)) {
      if (U1[k] <= 1 - tail_frac_SVL || !length(sub_SVL_only)) {
        Y1[k] <- if (length(sub_SVL_only)) sample(sub_SVL_only, 1L, TRUE) else u1
      } else {
        v1 <- (U1[k] - (1 - tail_frac_SVL)) / tail_frac_SVL
        Y1[k] <- u1 + evd::qgpd(v1, loc = 0, scale = sigma1_hat, shape = xi_eq_hat)
      }
    }
    y1_E_M2_boot <- Y1[Y1 > u1]
  } else {
    y1_E_M2_boot <- numeric(0L)
  }
  
  if (n_tl_only > 0L) {
    U2 <- runif(n_tl_only)
    Y2 <- numeric(n_tl_only)
    for (k in seq_len(n_tl_only)) {
      if (U2[k] <= 1 - tail_frac_TL || !length(sub_TL_only)) {
        Y2[k] <- if (length(sub_TL_only)) sample(sub_TL_only, 1L, TRUE) else u2
      } else {
        v2 <- (U2[k] - (1 - tail_frac_TL)) / tail_frac_TL
        Y2[k] <- u2 + evd::qgpd(v2, loc = 0, scale = sigma2_hat, shape = xi_eq_hat)
      }
    }
    y2_M1_E_boot <- Y2[Y2 > u2]
  } else {
    y2_M1_E_boot <- numeric(0L)
  }
  
  # (3) Refit H0: xi1=xi2 on bootstrap sample
  ll_cens_b <- make_ll_cens(X_biv_boot, u_vec)
  ll_miss_b <- function(theta) {
    ll1 <- if (length(y1_E_M2_boot)) sum(lgp_exceed(y1_E_M2_boot, u1, theta["scale1"], theta["shape1"])) else 0
    ll2 <- if (length(y2_M1_E_boot)) sum(lgp_exceed(y2_M1_E_boot, u2, theta["scale2"], theta["shape2"])) else 0
    ll1 + ll2
  }
  loglik_b <- function(theta) ll_cens_b(theta) + ll_miss_b(theta)
  
  negloglik_b <- function(p) {
    theta <- c(scale1 = p[1], shape1 = p[2], scale2 = p[3], shape2 = p[2], dep = p[4])
    -loglik_b(theta)
  }
  
  start_b <- c(scale1 = sigma1_hat, xi = xi_eq_hat, scale2 = sigma2_hat, dep = alpha_logistic_hat)
  
  opt_b <- try(optim(unname(start_b), negloglik_b, method = "Nelder-Mead"), silent = TRUE)
  if (inherits(opt_b, "try-error") || opt_b$convergence != 0) next
  
  boot_conv[b] <- TRUE
  phi_b <- setNames(opt_b$par, names(start_b))
  
  boot_sigma1[b] <- phi_b["scale1"]
  boot_xi[b]     <- phi_b["xi"]
  boot_sigma2[b] <- phi_b["scale2"]
  boot_y1star[b] <- u1 - boot_sigma1[b] / boot_xi[b]
  boot_y2star[b] <- u2 - boot_sigma2[b] / boot_xi[b]
}

cat("\nParametric bootstrap (common-tail logistic GP):\n")
cat("  Successful replicates:", sum(boot_conv), "out of", B_boot, "\n")

boot_ok <- boot_conv & is.finite(boot_xi) & (boot_xi < 0)
stopifnot(any(boot_ok))

alpha_half <- (1 - ci_level) / 2

# 2-sided (equal-tail) bootstrap CIs for GP parameters
ci_sigma1 <- quantile(boot_sigma1[boot_ok], c(alpha_half, 1 - alpha_half), na.rm = TRUE)
ci_sigma2 <- quantile(boot_sigma2[boot_ok], c(alpha_half, 1 - alpha_half), na.rm = TRUE)
ci_xi     <- quantile(boot_xi[boot_ok],     c(alpha_half, 1 - alpha_half), na.rm = TRUE)

cat("\n", 100 * ci_level, "% bootstrap CIs (log scale):\n", sep = "")
cat("  sigma1:", ci_sigma1[1], "to", ci_sigma1[2], "\n")
cat("  sigma2:", ci_sigma2[1], "to", ci_sigma2[2], "\n")
cat("  xi    :", ci_xi[1],     "to", ci_xi[2],     "\n")

# 1-sided upper bounds for endpoints (log scale) + back-transform
ub_y1_log <- as.numeric(quantile(boot_y1star[boot_ok], probs = ci_level, na.rm = TRUE))
ub_y2_log <- as.numeric(quantile(boot_y2star[boot_ok], probs = ci_level, na.rm = TRUE))

cat("\nOne-sided ", 100 * ci_level, "% upper bounds for log endpoints:\n", sep = "")
cat("  y*_SVL (log) ≤ ", ub_y1_log, "\n", sep = "")
cat("  y*_TL  (log) ≤ ", ub_y2_log, "\n", sep = "")

cat("\nOne-sided ", 100 * ci_level, "% upper bounds (original scale):\n", sep = "")
cat("  SVL* ≤ ", exp(ub_y1_log), "\n", sep = "")
cat("  TL*  ≤ ", exp(ub_y2_log), "\n", sep = "")

# ============================================================
# 8) Exceedance probability for the Stokes alligator
#    (SVL=239 cm, TL=450 cm)
# ============================================================

# GPD survival for y > u on the log-scale (returns conditional survival; if y<=u -> 1)
gpd_survival <- function(y, u, sigma, xi) {
  x <- y - u
  if (x <= 0) return(1)
  if (abs(xi) < 1e-8) return(exp(-x / sigma))
  term <- 1 + xi * x / sigma
  if (term <= 0) return(0) # beyond endpoint
  term^(-1 / xi)
}

SVL_stokes <- 239  # cm
TL_stokes  <- 450  # cm
log_SVL_stokes <- log(SVL_stokes)
log_TL_stokes  <- log(TL_stokes)

cat("\nStokes alligator (log scale):\n")
cat("  log(SVL) =", log_SVL_stokes, "\n")
cat("  log(TL)  =", log_TL_stokes,  "\n")

cond_SVL <- gpd_survival(log_SVL_stokes, u1, sigma1_hat, xi_eq_hat)
cond_TL  <- gpd_survival(log_TL_stokes,  u2, sigma2_hat, xi_eq_hat)

p_SVL <- tail_frac_SVL * cond_SVL
p_TL  <- tail_frac_TL  * cond_TL

cat("\nUnivariate predictive exceedance probabilities for Stokes:\n")
cat("  P(SVL > 239 cm) ≈", p_SVL)
cat("  P(TL  > 450 cm) ≈", p_TL)

# Independence approximation
p_joint_indep <- p_SVL * p_TL
cat("\nJoint exceedance (independence):\n")
cat("  P(SVL>239, TL>450) ≈", p_joint_indep)

# Logistic EV copula Monte Carlo (Gumbel)
theta_gumbel <- 1 / alpha_logistic_hat
gumbel_cop   <- gumbelCopula(param = theta_gumbel, dim = 2)

# transform to copula scale U = F(y): for y>u, F(y)=1-P(Y>y)=1-tail_frac*cond
u1_u <- 1 - tail_frac_SVL * cond_SVL
u2_u <- 1 - tail_frac_TL  * cond_TL

cat("\nCopula thresholds for Stokes (U scale):\n")
cat("  u1_u = F_SVL(239) ≈", u1_u, "\n")
cat("  u2_u = F_TL (450) ≈", u2_u, "\n")

stopifnot(is.finite(u1_u), is.finite(u2_u), u1_u > 0, u1_u < 1, u2_u > 0, u2_u < 1)

set.seed(2025)
N_mc <- 2e5L
U_mc <- rCopula(N_mc, gumbel_cop)

p_joint_log  <- mean((U_mc[, 1] > u1_u) & (U_mc[, 2] > u2_u))
se_joint_log <- sqrt(p_joint_log * (1 - p_joint_log) / N_mc)

cat("\nLogistic EV copula (MC) joint exceedance:\n")
cat("  P(SVL>239, TL>450) ≈", p_joint_log, " (MC s.e. ≈", se_joint_log, ")\n")

cat("\nComparison (independence vs logistic):\n")
cat("  Indep : P ≈", p_joint_indep, " -> RP ≈", 1 / p_joint_indep, "\n")
cat("  Logist: P ≈", p_joint_log,   " -> RP ≈", 1 / p_joint_log,   "\n")

N_year <- 7600  # approx mean annual harvest 2009–2024
cat("\nReturn periods in years:\n")
cat("  Indep : ≈", (1 / p_joint_indep) / N_year, "\n")
cat("  Logist: ≈", (1 / p_joint_log)   / N_year, "\n")

# ============================================================
# 9) Joint bootstrap distribution for (SVL*, TL*) : KDE + HPD contour
#     + JOINT MAP = 2D KDE mode (NOT the MLE plug-in)
# ============================================================

svl_boot_cm <- exp(boot_y1star[boot_ok])
tl_boot_cm  <- exp(boot_y2star[boot_ok])

boot_joint <- tibble(SVL = svl_boot_cm, TL = tl_boot_cm)

# marginal one-sided upper bounds (for comparison)
ub_svl_marg <- as.numeric(quantile(svl_boot_cm, probs = ci_level, na.rm = TRUE))
ub_tl_marg  <- as.numeric(quantile(tl_boot_cm,  probs = ci_level, na.rm = TRUE))

cat("\n", 100 * ci_level, "% one-sided marginal upper bounds (original scale):\n", sep = "")
cat("  SVL* ≤ ", ub_svl_marg, " cm\n", sep = "")
cat("  TL*  ≤ ", ub_tl_marg,  " cm\n", sep = "")

# ---- KDE window (use central mass; adjust if you want) ----
svl_window <- quantile(svl_boot_cm, c(0.01, 0.99), na.rm = TRUE)
tl_window  <- quantile(tl_boot_cm,  c(0.01, 0.99), na.rm = TRUE)

inside_win <- with(boot_joint,
                   SVL >= svl_window[1] & SVL <= svl_window[2] &
                     TL  >= tl_window[1]  & TL  <= tl_window[2])

boot_joint_kde <- boot_joint[inside_win, , drop = FALSE]
if (nrow(boot_joint_kde) < 50L) {
  stop("Too few bootstrap draws inside KDE window; relax window or reduce trimming.")
}

# ---- 2D KDE ----
kde_joint <- MASS::kde2d(
  x    = boot_joint_kde$SVL,
  y    = boot_joint_kde$TL,
  n    = 200,
  lims = c(svl_window[1], svl_window[2], tl_window[1], tl_window[2])
)

# JOINT MAP = mode of the 2D KDE grid
ij_max <- which(kde_joint$z == max(kde_joint$z), arr.ind = TRUE)[1, ]
map_svl_joint <- kde_joint$x[ij_max[1]]
map_tl_joint  <- kde_joint$y[ij_max[2]]

cat("\n2D KDE joint MAP (SVL*, TL*):\n")
cat("  SVL*_MAP ≈", round(map_svl_joint, 2), "cm\n")
cat("  TL*_MAP  ≈", round(map_tl_joint,  2), "cm\n")

# ---- HPD contour threshold using grid mass ----
dx <- kde_joint$x[2] - kde_joint$x[1]
dy <- kde_joint$y[2] - kde_joint$y[1]
zv <- as.vector(kde_joint$z)

ord <- order(zv, decreasing = TRUE)
cum_mass <- cumsum(zv[ord]) * dx * dy
hpd_threshold <- zv[ord][which(cum_mass >= ci_level)[1]]

kde_df <- with(kde_joint, expand.grid(SVL = x, TL = y) |> transform(z = as.vector(z)))

# ---- Plot ----
p_joint <- ggplot() +
  geom_density_2d_filled(data = boot_joint_kde, aes(SVL, TL), alpha = 0.7) +
  geom_contour(data = kde_df, aes(SVL, TL, z = z),
               breaks = hpd_threshold, colour = "black", linewidth = 1.0) +
  geom_point(data = boot_joint_kde, aes(SVL, TL),
             alpha = 0.25, size = 0.6, color = "grey20") +
  # JOINT MAP point (KDE mode)
  geom_point(aes(x = map_svl_joint, y = map_tl_joint), color = "purple", size = 3) +
  # Stokes alligator
  geom_point(aes(x = SVL_stokes, y = TL_stokes), color = "red", size = 3,
             shape = 4, stroke = 1.2) +
  coord_cartesian(xlim = svl_window, ylim = tl_window) +
  labs(x = "SVL endpoint (cm)", y = "TL endpoint (cm)", fill = "Density level") +
  theme_science

print(p_joint)
save_png(p_joint, "SVL_TL_joint_endpoint_bootstrap.png", w = 6.5, h = 4.5)
