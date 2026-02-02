#################################################
# Title: Hierarchical mass allometry (conjugate)
# Description: Mass-from-bones with latent lnM,
# using closed-form Normal–Inverse–Gamma updates.
# Author: Bastiaan Aelbrecht
# Date: 15/08/2025
#################################################
# ---- Packages ----
library(MASS)       # mvrnorm
library(dplyr)
library(tidyr)
library(tibble)
library(purrr)
library(ggplot2)
library(bayesplot)  # optional plotting of draws
library(MASSTIMATE)
color_scheme_set("blue")

# ---- Traits (J = 4) ----
trait_levels <- c("lnCF","lnCH","lnLF","lnLH")

# ---- Prior centres (Alexander 1979), log scale ----
ah_CF <- log(0.43); ah_CH <- log(0.35); ah_LF <- log(5.24); ah_LH <- log(4.24)
bh_CF <- 0.36; bh_CH <- 0.38; bh_LF <- 0.36; bh_LH <- 0.36

alpha_hat <- c(ah_CF, ah_CH, ah_LF, ah_LH)
beta_hat  <- c(bh_CF, bh_CH, bh_LF, bh_LH)

# ---- Conjugate NIG hyperparameters ----
# theta_j | sigùa_j^2 ~ N( mu0_j, sigma_j^2 * V0_j ),  sigma_j^2 ~ IG(a0, b0)
sd_alpha <- 0.5
sd_beta  <- 0.1

V0_list <- list()
mu0_list <- list()

for (j in 1:length(trait_levels)) {
  V0_list[[j]] <- diag(c(sd_alpha^2, sd_beta^2))
  mu0_list[[j]] <- c(alpha_hat[j], beta_hat[j])
}

a0 <- 1.1
b0 <- 0.7
b0/(a0-1)

# inverse gamma sampling
rinvgamma <- function(n, shape, rate) {
  # InvGamma(shape=a, rate=b): density ∝ b^a x^{-(a+1)} exp(-b/x)
  1 / rgamma(n, shape = shape, rate = rate)
}

nig_posterior <- function(y, X, mu0, V0, a0, b0) {
  # Prior: theta|sigma ~ N(mu0, sigma^2 V0), sigma^2 ~ InvGamma(a0,b0)
  # Returns list(mu_n, V_n, a_n, b_n, XTX, XTy)
  XtX <- crossprod(X) # X^T X 2x2 matrix
  XtY <- crossprod(X, y) # X^T Y 2x1 vector ; these are the sufficient statistics for Gaussian regression.
  V0_inv <- solve(V0) # inverse matrix-> prior precision
  Vn_inv <- V0_inv + XtX # posterior precision
  Vn <- solve(Vn_inv) # again invert-> to obtain posterior covariance FOR EVERY UNIT OF SIGMA (because of NIG structure)
  mu_n <- Vn %*% (V0_inv %*% mu0 + XtY) # matrix multiplication
  # interpretation of this precision weighted update: "pseudo data V0_inv mu0 + "score, all scaled by posterior covariance"
  
  a_n <- a0 + length(y) / 2 # degrees of freedom update
  quad <- as.numeric(crossprod(y, y) + crossprod(mu0, V0_inv %*% mu0) - crossprod(mu_n, Vn_inv %*% mu_n))
  b_n <- b0 + 0.5 * quad  # after completing the square
  
  list(mu_n = as.numeric(mu_n), Vn = Vn, a_n = a_n, b_n = b_n,
       XtX = XtX, XtY = XtY, n = length(y))
}

sample_from_nig <- function(n_samp, mu_n, Vn, a_n, b_n) {
  # Draw σ2 ~ InvGamma(a_n,b_n), then θ ~ N(mu_n, σ2 Vn) (aka conditional drawing)
  sig2 <- rinvgamma(n_samp, shape = a_n, rate = b_n)
  # cbind(sig2) turns the vector into a 1-column matrix so apply(..., 1, ...)
  # loops over rows
  # t(.) transposes so rows correspond to samples and columns to parameters (alpha, beta)
  theta <- t(apply(cbind(sig2), 1, function(s2) {
    as.numeric(mvrnorm(1, mu = mu_n, Sigma = s2 * Vn))
  }))
  colnames(theta) <- c("alpha","beta")
  list(theta = theta, sigma2 = sig2)
}

# ---------------------------------------
# Data preparation (extants → long)
# ---------------------------------------
data("extants")

extants_long <- extants %>%
  mutate(
    lnM  = log(as.numeric(BM)),
    lnCF = log(as.numeric(FC)),
    lnCH = log(as.numeric(HC)),
    lnLF = log(as.numeric(Femur.Length)),
    lnLH = log(as.numeric(Humerus.Length))
  ) %>%
  dplyr::select(Species, lnM, lnCF, lnCH, lnLF, lnLH) %>% # apparently MASS also has a "select"
  pivot_longer(cols = starts_with("ln"),
               names_to = "trait", values_to = "val") %>%
  mutate(
    is_mass  = (trait == "lnM"),
    trait_id = match(trait, trait_levels)
  )

# Training rows (require lnM and trait observed & finite)
train_rows <- extants_long %>%
  filter(!is_mass, trait %in% trait_levels) %>%
  group_by(Species) %>%
  mutate(lnM = extants_long$val[match(paste0(Species, "_lnM"),
                                      paste0(extants_long$Species, "_", extants_long$trait))]) %>%
  ungroup() %>%
  filter(is.finite(lnM), is.finite(val))

# Build per-trait datasets
trait_train <- lapply(seq_along(trait_levels), function(j) {
  dfj <- train_rows %>% filter(trait_id == j)
  list(
    y = dfj$val,
    X = cbind(1, dfj$lnM)  # [1, lnM]
  )
})

# ---------------------------------------
# TRAIN: Conjugate posteriors per trait
# ---------------------------------------
post_list <- lapply(seq_along(trait_levels), function(j) {
  nig_posterior(
    y   = trait_train[[j]]$y,
    X   = trait_train[[j]]$X,
    mu0 = mu0_list[[j]],
    V0  = V0_list[[j]],
    a0  = a0, b0 = b0
  )
})
names(post_list) <- trait_levels

# ---------------------------------------
# Posterior draws for plotting / summaries
# ---------------------------------------
S_draws <- 10000
draws_by_trait <- lapply(seq_along(trait_levels), function(j) {
  pj <- post_list[[j]]
  smp <- sample_from_nig(S_draws, mu_n = pj$mu_n, Vn = pj$Vn, a_n = pj$a_n, b_n = pj$b_n)
  tibble::tibble(
    trait = trait_levels[j],
    alpha = smp$theta[, "alpha"],
    beta  = smp$theta[, "beta"],
    sigma = sqrt(smp$sigma2)
  )
})
draws_df <- bind_rows(draws_by_trait)


# Summaries by trait
summ_trait <- draws_df %>%
  pivot_longer(cols = c(alpha, beta, sigma), names_to = "param", values_to = "value") %>%
  group_by(trait, param) %>%
  summarise(
    mean = mean(value),
    sd   = sd(value),
    q025 = quantile(value, 0.025),
    q50  = quantile(value, 0.5),
    q975 = quantile(value, 0.975),
    .groups = "drop"
  )

# Also compute inverse regression parameters: log x = alpha' + beta' log y
# For each draw, calculate alpha' = -alpha/beta, beta' = 1/beta, and log10(exp(alpha'))
inv_draws_df <- draws_df %>%
  mutate(
    alpha_inv = -alpha / beta,
    beta_inv  = 1 / beta,
    log10_exp_alpha_inv = log10(exp(alpha_inv))
  )

summ_trait_inv <- inv_draws_df %>%
  pivot_longer(cols = c(alpha_inv, beta_inv, log10_exp_alpha_inv), names_to = "param", values_to = "value") %>%
  group_by(trait, param) %>%
  summarise(
    mean = mean(value),
    sd   = sd(value),
    q025 = quantile(value, 0.025),
    q50  = quantile(value, 0.5),
    q975 = quantile(value, 0.975),
    .groups = "drop"
  )

print(summ_trait)
print(summ_trait_inv)
print(summ_trait)


# Intercept posterior for all traits
ggplot(draws_df, aes(x = alpha, color = trait, fill = trait)) +
  geom_density(alpha = 0.3) +
  labs(title = "Posterior densities of alpha (intercept, log scale) by trait")

# Slope posterior for all traits
ggplot(draws_df, aes(x = beta, color = trait, fill = trait)) +
  geom_density(alpha = 0.3) +
  labs(title = "Posterior densities of beta (slope) by trait")

# Plot proportionality constant (exp(alpha)) on original scale
ggplot(draws_df, aes(x = exp(alpha), color = trait, fill = trait)) +
  geom_density(alpha = 0.3) +
  labs(title = "Posterior densities of proportionality constant (exp(alpha)) by trait",
       x = "Proportionality constant (original scale)")

# Plot all inverse intercepts (alpha') on one density plot
ggplot(inv_draws_df, aes(x = alpha_inv, color = trait, fill = trait)) +
  geom_density(alpha = 0.3) +
  labs(title = "Posterior densities of alpha' (inverse intercept) by trait",
       x = "alpha' = -alpha/beta")

# Plot all inverse slopes (beta') on one density plot
ggplot(inv_draws_df, aes(x = beta_inv, color = trait, fill = trait)) +
  geom_density(alpha = 0.3) +
  labs(title = "Posterior densities of beta' (inverse slope) by trait",
       x = "beta' = 1/beta")

# Plot posteriors of the residual standard deviation (sigma) by trait
ggplot(draws_df, aes(x = sigma, color = trait, fill = trait)) +
  geom_density(alpha = 0.3) +
  labs(title = "Posterior densities of residual standard deviation (sigma) by trait",
       x = "Residual SD (sigma)")

# goodness of fit plot
# Goodness of fit plots: observed vs fitted (posterior mean) per trait
gof_plots <- lapply(seq_along(trait_levels), function(j) {
  yj <- trait_train[[j]]$y
  Xj <- trait_train[[j]]$X
  pj <- post_list[[j]]
  theta_hat <- pj$mu_n
  muj_hat <- as.vector(Xj %*% theta_hat)
  df_plot <- tibble::tibble(
    observed = yj,
    fitted = muj_hat
  )
  ggplot(df_plot, aes(x = observed, y = fitted)) +
    geom_point(alpha = 0.7) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    labs(
      title = paste("Observed vs Fitted for", trait_levels[j]),
      x = "Observed",
      y = "Fitted"
    ) +
    theme_minimal()
})

# Optionally print or plot all
for (p in gof_plots) print(p)

# ---- GOF diagnostics: R^2, AIC/BIC, WAIC, LOO per trait ----

compute_gof_global <- function() {
  # Stack training data, fitted means and sigmas at posterior mean
  y_all   <- c()
  mu_all  <- c()
  n_total <- 0
  ll_sum  <- 0
  k_total <- 0
  
  # For Bayesian R^2 draws, build per-draw fitted vectors and residual variances
  mu_draw_list <- list()
  sig2_draw_list <- list()
  
  for (j in seq_along(trait_levels)) {
    yj <- trait_train[[j]]$y
    Xj <- trait_train[[j]]$X
    pj <- post_list[[j]]
    
    # Posterior-mean params
    theta_hat  <- pj$mu_n
    sigma2_hat <- pj$b_n / (pj$a_n - 1)
    muj_hat    <- as.vector(Xj %*% theta_hat)
    
    # Accumulate for classical metrics
    y_all  <- c(y_all,  yj)
    mu_all <- c(mu_all, muj_hat)
    n_total <- n_total + length(yj)
    ll_sum  <- ll_sum + sum(dnorm(yj, mean = muj_hat, sd = sqrt(sigma2_hat), log = TRUE))
    k_total <- k_total + 3  # alpha, beta, sigma^2 per trait
    
    # Per-draw fitted means & residual variances
    dr <- draws_df %>% dplyr::filter(trait == trait_levels[j])
    # mu per draw: matrix [n_j x S]
    mu_mat <- sapply(seq_len(nrow(dr)), function(s) {
      as.vector(Xj %*% c(dr$alpha[s], dr$beta[s]))
    })
    # residual variance per draw: vector length S, recycle to n_j
    sig2_vec <- dr$sigma^2
    sig2_mat <- matrix(rep(sig2_vec, each = nrow(mu_mat)),
                       nrow = nrow(mu_mat), ncol = ncol(mu_mat))
    
    mu_draw_list[[j]]   <- mu_mat
    sig2_draw_list[[j]] <- sig2_mat
  }
  
  # Classical global R^2
  R2 <- 1 - sum((y_all - mu_all)^2) / sum((y_all - mean(y_all))^2)
  AIC <- -2 * ll_sum + 2 * k_total
  BIC <- -2 * ll_sum + log(n_total) * k_total
  
  
  # --- Percent prediction error on original (mm) scale ---
  # Use exp(mu) as the median prediction under log-normal errors
  y_all_mm  <- exp(y_all)
  mu_all_mm <- exp(mu_all)
  perc_err  <- 100 * (y_all_mm - mu_all_mm) / y_all_mm
  MAPE <- mean(abs(perc_err))      # mean absolute percent error
  
  tibble::tibble(
    scope        = "global_all_traits",
    n_total      = n_total,
    k_total      = k_total,
    R2   = R2,
    AIC          = AIC,
    BIC          = BIC,
    MAPE         = MAPE
  )
}

gof_global <- compute_gof_global()
print(gof_global)

# ---- Posterior predictive checks ----

# helper: posterior predictive summary for one trait/specimen
ppc_predict_trait <- function(lnM_i, y_obs, post) {
  x_i <- c(1, lnM_i)
  df  <- 2 * post$a_n
  loc <- as.numeric(x_i %*% post$mu_n)
  scale <- sqrt((post$b_n / post$a_n) * (1 + x_i %*% post$Vn %*% x_i))
  
  # draws from Student-t predictive
  yrep <- loc + scale * rt(5000, df = df)
  
  tibble::tibble(
    observed = y_obs,
    mean     = mean(yrep),
    q025     = quantile(yrep, 0.025),
    q50      = quantile(yrep, 0.5),
    q975     = quantile(yrep, 0.975)
  )
}

# apply to all extant specimens and traits
ppc_results <- train_rows %>%
  rowwise() %>%
  mutate(
    ppc = list(ppc_predict_trait(
      lnM_i = lnM,
      y_obs = val,
      post  = post_list[[trait]]
    ))
  ) %>%
  unnest(ppc)

# plot: observed vs predictive mean + intervals
ggplot(ppc_results, aes(x = observed, y = mean, color = trait)) +
  geom_point(alpha = 0.6) +
  geom_errorbar(aes(ymin = q025, ymax = q975), alpha = 0.3) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(
    title = "Posterior predictive checks per trait",
    x = "Observed (log scale)",
    y = "Posterior predictive mean ± 95% CI"
  ) +
  facet_wrap(~trait, scales = "free") +
  theme_minimal()

# ============================================================
# MASS-ONLY PREDICTION (Gaussian trait factors + Laplace in m)
# ============================================================

# NOTE: We replace the Student-t predictive by a Gaussian:
# p(y_j | m, Ω) ≈ N( μ_j(m), v_j(m) ).
# (Old: t_{2 a_{jn}}(loc = x^T μ_n, scale^2 = (b_n/a_n)(1 + x^T Vn x)) )

# ---- trait id sanity helper ----
normalize_trait_ids <- function(x, allowed) {
  x_chr <- trimws(as.character(x))
  bad   <- setdiff(x_chr, allowed)
  if (length(bad)) {
    stop("Unknown trait id(s): ", paste(bad, collapse = ", "),
         "\nAllowed: ", paste(allowed, collapse = ", "))
  }
  x_chr
}

# ---- Gaussian trait factor components (μ_j(m), v_j(m), v'_j(m), v''_j(m)) ----
gaussian_trait_params <- function(post, m) {
  x  <- c(1, m)
  mu <- as.numeric(x %*% post$mu_n)
  
  cfac <- (post$b_n / post$a_n)
  Vn   <- post$Vn
  
  q    <- as.numeric(1 + x %*% Vn %*% x)           # 1 + x^T Vn x
  v    <- cfac * q                                 # v_j(m)
  v1   <- 2 * cfac * (Vn[1,2] + Vn[2,2] * m)       # v'_j(m)
  v2   <- 2 * cfac * Vn[2,2]                       # v''_j(m)
  
  list(mu = mu, v = v, v1 = v1, v2 = v2)
}

# log N(y | μ(m), v(m))
log_norm_factor <- function(y, mu, v) {
  -0.5 * (log(2*pi*v) + (y - mu)^2 / v)
}

# ---- Laplace posterior for m with exact derivatives (your formulas) ----
predict_mass_laplace <- function(y_obs, trait_ids,
                                 prior_mean = 0, prior_sd = 10,
                                 m_init = prior_mean,
                                 maxit = 50, tol = 1e-8, step_shrink = 0.5) {
  allowed <- names(post_list)
  trait_ids <- normalize_trait_ids(trait_ids, allowed)
  if (is.null(names(y_obs))) stop("y_obs must be a named numeric vector.")
  names(y_obs) <- normalize_trait_ids(names(y_obs), allowed)
  stopifnot(all(trait_ids %in% names(y_obs)))
  
  y_vec <- as.numeric(y_obs[trait_ids])
  pjL   <- post_list[trait_ids]
  
  # ℓ(m)
  ell <- function(m) {
    lp <- dnorm(m, mean = prior_mean, sd = prior_sd, log = TRUE)
    for (j in seq_along(pjL)) {
      gp <- gaussian_trait_params(pjL[[j]], m)
      lp <- lp + log_norm_factor(y_vec[j], gp$mu, gp$v)
    }
    lp
  }
  
  # ℓ'(m)
  ell_prime <- function(m) {
    g  <- -(m - prior_mean) / (prior_sd^2)
    for (j in seq_along(pjL)) {
      pj <- pjL[[j]]
      μ2 <- pj$mu_n[2]
      gp <- gaussian_trait_params(pj, m)
      r  <- y_vec[j] - gp$mu
      v  <- gp$v
      v1 <- gp$v1
      g  <- g + (μ2 * r) / v + 0.5 * ((r^2 / (v^2)) - 1 / v) * v1
    }
    g
  }
  
  # ℓ''(m)
  ell_second <- function(m) {
    H  <- -1 / (prior_sd^2)
    for (j in seq_along(pjL)) {
      pj <- pjL[[j]]
      μ2 <- pj$mu_n[2]
      gp <- gaussian_trait_params(pj, m)
      r  <- y_vec[j] - gp$mu
      v  <- gp$v
      v1 <- gp$v1
      v2 <- gp$v2
      H <- H + ( - (μ2^2) / v
                 - (2 * μ2 * r * v1) / (v^2)
                 + 0.5 * ( (r^2 * v2) / (v^2)
                           - (2 * r^2 * (v1^2)) / (v^3)
                           - v2 / v
                           + (v1^2) / (v^2) ) )
    }
    H
  }
  
  # Newton with light backtracking
  m <- m_init
  for (k in 1:maxit) {
    g  <- ell_prime(m)
    H  <- ell_second(m)
    if (!is.finite(H) || H >= 0) {
      m <- m + rnorm(1, 0, 1e-4) # small nudge if curvature unusable
      next
    }
    step  <- -g / H
    m_new <- m + step
    
    tries <- 0
    while (!is.finite(ell(m_new)) || (ell(m_new) < ell(m) && tries < 6)) {
      step  <- step * step_shrink
      m_new <- m + step
      tries <- tries + 1
    }
    if (abs(m_new - m) < tol) { m <- m_new; break }
    m <- m_new
  }
  
  Hhat <- ell_second(m)
  Vhat <- if (is.finite(Hhat) && Hhat < 0) -1 / Hhat else NA_real_
  
  list(
    mean = as.numeric(m),
    var  = as.numeric(Vhat),
    mode = as.numeric(m),
    hess = as.numeric(Hhat),
    logpost_at_mode = as.numeric(ell(m)),
    converged = is.finite(Vhat) && Vhat > 0
  )
}

# ---- helper: density curve for Laplace Normal N( m̂ , V̂ ) on a grid ----
laplace_curve <- function(m_post, grid = NULL, width = 5) {
  stopifnot(is.finite(m_post$var), m_post$var > 0)
  if (is.null(grid)) {
    sd_hat <- sqrt(m_post$var)
    grid <- seq(m_post$mean - width * sd_hat,
                m_post$mean + width * sd_hat, length.out = 400)
  }
  tibble::tibble(
    m = grid,
    dens = dnorm(m, mean = m_post$mean, sd = sqrt(m_post$var))
  )
}

# ============================================================
# Cross-validation with randomly missing traits in TEST data
# ============================================================

set.seed(42)

# --- masking controls ---
R_reps  <- 50       # number of missingness replicates per CV fold
p_miss  <- 0.90     # probability a *trait* is set missing in test rows (independent per trait)

# Optionally set per-trait probabilities instead of a scalar:
# p_miss <- c(lnCF = 0.25, lnCH = 0.25, lnLF = 0.40, lnLH = 0.40)

# --- helper: turn scalar into named per-trait vector ---
mk_p_miss_vec <- function(p) {
  if (length(p) == 1) {
    setNames(rep(p, length(trait_levels)), trait_levels)
  } else {
    pv <- setNames(rep(0, length(trait_levels)), trait_levels)
    pv[names(p)] <- p
    pv
  }
}

# --- helper: mask traits in a single-species LONG data frame (val -> NA with prob p_miss) ---
mask_traits_long <- function(rows_i, p_vec) {
  rows_i %>%
    dplyr::mutate(
      val = ifelse(is.finite(val) & runif(dplyr::n()) < p_vec[trait], NA_real_, val)
    )
}

# -------------------------------------------
# CFH baseline: fold-honest fitting + POINT impute for missing CF/CH
# -------------------------------------------

fit_cfh_impute_models <- function(train_rows_fold) {
  wide_train <- train_rows_fold %>%
    dplyr::select(Species, trait, val, lnM) %>%
    tidyr::pivot_wider(names_from = trait, values_from = val) %>%
    dplyr::filter(is.finite(lnM))
  
  # CFH model needs at least some rows with both CF & CH
  train_cfh <- wide_train %>% dplyr::filter(is.finite(lnCF), is.finite(lnCH))
  if (nrow(train_cfh) < 5) return(NULL)
  
  lm_cfh <- lm(lnM ~ I(log(exp(lnCF) + exp(lnCH))), data = train_cfh)
  sigma_hat <- summary(lm_cfh)$sigma
  
  # Imputation lines (fit where both present)
  cf_from_ch <- try(lm(lnCF ~ lnCH, data = train_cfh), silent = TRUE)
  if (inherits(cf_from_ch, "try-error")) cf_from_ch <- NULL
  
  ch_from_cf <- try(lm(lnCH ~ lnCF, data = train_cfh), silent = TRUE)
  if (inherits(ch_from_cf, "try-error")) ch_from_cf <- NULL
  
  list(lm_cfh = lm_cfh, sigma = sigma_hat,
       cf_from_ch = cf_from_ch, ch_from_cf = ch_from_cf)
}

predict_cfh_with_impute_point <- function(row_one, models) {
  if (is.null(models) || is.null(models$lm_cfh)) return(NULL)
  
  # pull scalars safely (columns exist after pivot_wider, but may be NA)
  a <- suppressWarnings(as.numeric(row_one[["lnCF"]]))
  b <- suppressWarnings(as.numeric(row_one[["lnCH"]]))
  
  # point-impute CF from CH
  if (!is.finite(a) && is.finite(b) && !is.null(models$cf_from_ch)) {
    a <- as.numeric(predict(models$cf_from_ch, newdata = data.frame(lnCH = b)))
  }
  # point-impute CH from CF
  if (!is.finite(b) && is.finite(a) && !is.null(models$ch_from_cf)) {
    b <- as.numeric(predict(models$ch_from_cf, newdata = data.frame(lnCF = a)))
  }
  
  # if either still missing, we cannot make a CFH prediction
  if (!is.finite(a) || !is.finite(b)) return(NULL)
  
  # predict lnM with the fold-honest CFH model; provide lnCF & lnCH and let the formula compute the transform
  pr <- predict(models$lm_cfh,
                newdata = data.frame(lnCF = a, lnCH = b),
                se.fit  = TRUE)
  
  pred_sd <- sqrt(drop(pr$se.fit)^2 + models$sigma^2)
  
  tibble::tibble(pred_lnM = as.numeric(pr$fit),
                 pred_sd  = as.numeric(pred_sd))
}


# -------------------------------------------
# K-fold splits (species-wise)
# -------------------------------------------
K <- 5
species_ids <- unique(train_rows$Species)
fold_id <- sample(rep(1:K, length.out = length(species_ids)))

# -------------------------------------------
# Run CV with missing-at-test masking
# -------------------------------------------
p_vec <- mk_p_miss_vec(p_miss)

lm_results  <- list()
cfh_results <- list()

for (k in 1:K) {
  test_sp   <- species_ids[fold_id == k]
  train_sp  <- setdiff(species_ids, test_sp)
  
  # fold-specific refit (for trait factors + lnM prior)
  fit_k   <- refit_conjugate(train_sp)
  post_k  <- fit_k$post_list
  
  # --- IMPORTANT: make predict_mass_laplace use fold-specific posteriors ---
  post_list_old <- post_list
  post_list     <- post_k
  
  # CFH baseline models fit on this fold
  cfh_models <- fit_cfh_impute_models(train_rows %>% dplyr::filter(Species %in% train_sp))
  
  for (rep in 1:R_reps) {
    # Iterate species in test fold
    for (spi in test_sp) {
      rows_i <- train_rows %>% dplyr::filter(Species == spi)
      
      # Mask traits at random (test-time only)
      rows_i_masked <- mask_traits_long(rows_i, p_vec)
      
      # -------- Latent-mass prediction (use only available traits) --------
      y_obs_i <- setNames(rows_i_masked$val, rows_i_masked$trait)
      y_obs_i <- y_obs_i[names(y_obs_i) %in% trait_levels]
      y_obs_i <- y_obs_i[is.finite(y_obs_i)]
      
      if (length(y_obs_i) > 0) {
        res_lap <- predict_mass_laplace(
          y_obs      = y_obs_i,
          trait_ids  = names(y_obs_i),
          prior_mean = fit_k$lnM_mean,
          prior_sd   = fit_k$lnM_sd
        )
        lm_results[[length(lm_results)+1]] <- tibble::tibble(
          Species   = spi,
          true_lnM  = unique(rows_i$lnM),
          pred_lnM  = res_lap$mean,
          pred_sd   = sqrt(res_lap$var),
          n_traits  = length(y_obs_i),
          fold      = k,
          rep       = rep,
          model     = "Latent-mass"
        )
      }
      
      # -------- CFH baseline with point imputation --------
      wide_i <- rows_i_masked %>%
        dplyr::select(Species, trait, val, lnM) %>%
        tidyr::pivot_wider(names_from = trait, values_from = val)
      
      if (nrow(wide_i) == 1) {
        pr <- predict_cfh_with_impute_point(wide_i, cfh_models)
        if (!is.null(pr)) {
          cfh_results[[length(cfh_results)+1]] <- tibble::tibble(
            Species   = spi,
            true_lnM  = wide_i$lnM,
            pred_lnM  = pr$pred_lnM,
            pred_sd   = pr$pred_sd,
            n_traits  = sum(is.finite(c(wide_i$lnCF, wide_i$lnCH))),  # CF/CH count after masking/impute check
            fold      = k,
            rep       = rep,
            model     = "CF+H (point-impute)"
          )
        }
      }
    } # species
  }   # reps
  
  # restore global post_list
  post_list <- post_list_old
}

cv_latent <- dplyr::bind_rows(lm_results)
cv_cfh    <- dplyr::bind_rows(cfh_results)
cv_all    <- dplyr::bind_rows(cv_latent, cv_cfh)

# -------------------------------------------
# Metrics under missing-at-test masking
# -------------------------------------------
metrics <- cv_all %>%
  dplyr::mutate(
    error     = pred_lnM - true_lnM,
    abs_error = abs(error),
    cover95   = (true_lnM >= pred_lnM - 1.96*pred_sd) &
      (true_lnM <= pred_lnM + 1.96*pred_sd)
  ) %>%
  dplyr::group_by(model) %>%
  dplyr::summarise(
    N       = dplyr::n(),
    RMSE    = sqrt(mean(error^2, na.rm = TRUE)),
    MAE     = mean(abs_error, na.rm = TRUE),
    Cover95 = mean(cover95, na.rm = TRUE),
    .groups = "drop"
  )

print(metrics)

# Optional: stratify by how many traits survived masking (latent only)
by_nt <- cv_latent %>%
  dplyr::mutate(
    error     = pred_lnM - true_lnM,
    abs_error = abs(error),
    cover95   = (true_lnM >= pred_lnM - 1.96*pred_sd) &
      (true_lnM <= pred_lnM + 1.96*pred_sd)
  ) %>%
  dplyr::group_by(n_traits) %>%
  dplyr::summarise(
    n    = dplyr::n(),
    RMSE = sqrt(mean(error^2, na.rm = TRUE)),
    MAE  = mean(abs_error, na.rm = TRUE),
    cov95= mean(cover95, na.rm = TRUE),
    .groups = "drop"
  )
print(by_nt)

# Optional: scatter (pooled reps)
ggplot(cv_all, aes(x = true_lnM, y = pred_lnM, color = model)) +
  geom_point(alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, linetype="dashed") +
  labs(title = sprintf("K-fold CV with random test masking (p=%.2f), R=%d", ifelse(length(p_miss)==1, p_miss, mean(unname(p_miss))), R_reps),
       x = "Observed ln(M)", y = "Predicted ln(M)", color = "Model") +
  theme_minimal()

