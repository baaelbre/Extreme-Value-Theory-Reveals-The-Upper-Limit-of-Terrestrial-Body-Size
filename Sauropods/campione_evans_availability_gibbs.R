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
sd_beta  <- 0.01

V0_list <- list()
mu0_list <- list()

for (j in 1:length(trait_levels)) {
  V0_list[[j]] <- diag(c(sd_alpha^2, sd_beta^2))
  mu0_list[[j]] <- c(alpha_hat[j], beta_hat[j])
}

a0 <- 2.1
b0 <- 1.1
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

#----------------------------------------------------------------
# Cross Validation
# Leave-one-out CV (LOO): each data point is its own fold. Useful when data is small, but computationally heavy.
# “new data”: repeatedly split data into training and test parts:
# train on one, test on the other.
# Every data point is predicted once, always by a model that did not see it during training.
# evaluates how the model generalizes beyond the calibration data.
# uses all available data (each point is in the training set K−1 times, in the test set once).

# Test folds: species left out; you try to predict their masses using only their traits, 
# as if they were “fossil specimens”.
#----------------------------------------------------------------
# Fixes the RNG seed so the fold assignment (and any other randomness) is reproducible.
set.seed(42)
K <- 10 # number of folds (number of subsets used for training)
species_ids <- unique(train_rows$Species) # split at the species level so no species’ data leaks across train/test.
fold_id <- sample(rep(1:K, length.out = length(species_ids))) # randomly assign each species to one of the K folds

# Helper to (re)fit conjugate posteriors on a given training species set
refit_conjugate <- function(train_species) {
  # subset to current training data set
  train_rows_k <- train_rows %>% filter(Species %in% train_species)
  # loop over traits
  trait_train_k <- lapply(seq_along(trait_levels), function(j) {
    # pull rows for trait j in the training dataset
    dfj <- train_rows_k %>% filter(trait_id == j)
    # Builds response y (log‑trait) and design matrix X = [1, lnM] for that trait.
    list(y = dfj$val, X = cbind(1, dfj$lnM))
  })
  # list of, per trait a (y,X)
  
  # Another loop over traits to compute NIG posteriors.
  post_list_k <- lapply(seq_along(trait_levels), function(j) {
    # Runs conjugate update, using prior center mu0_list[[j]], 
    # prior covariance scale V0_list[[j]], and inverse‑gamma hyperparameters a0, b0.
    nig_posterior(
      y   = trait_train_k[[j]]$y,
      X   = trait_train_k[[j]]$X,
      mu0 = mu0_list[[j]],
      V0  = V0_list[[j]],
      a0  = a0, b0 = b0
    )
  })
  
  names(post_list_k) <- trait_levels # name list by trait labels (e.g., lnCF, lnCH, …) for easy indexing.
  list(post_list = post_list_k, # list of fitted posteriors per trait
       lnM_mean = mean(train_rows_k$lnM), # and a weak prior to be used
       lnM_sd   = sd(train_rows_k$lnM)) # if there's no volumetric info available
}

# Pick one fold and one species from its TEST set
k_demo <- 1 # show for one fold.
test_species  <- species_ids[fold_id == k_demo] # species assigned to the test fold
train_species <- setdiff(species_ids, test_species) # all other species form the training data
# basic set operations in R
fit_demo <- refit_conjugate(train_species)

# Choose one species from the TEST set to demo the incremental plot
demo_sp <- test_species[11]
demo_sp # witstaartgnoe

# Grab that species' observed traits (on log scale) and its true lnM
demo_rows <- train_rows %>% filter(Species == demo_sp) # All (long‑format) rows for that species.
lnM_true  <- unique(demo_rows$lnM)
M_true_kg <- exp(lnM_true)/1000 # redelijk lichte gnoe -> misschien nog niet volwassen.

# Build named vector of observed traits for this specimen
y_obs_demo <- setNames(demo_rows$val, demo_rows$trait)
# Keep only traits that are in our model
y_obs_demo <- y_obs_demo[names(y_obs_demo) %in% trait_levels]

# Choose an order in which traits are added (change if you like)
incremental_orders <- list(
  "Prior only" = character(0),
  "1 trait: lnCF" = intersect("lnCF", names(y_obs_demo)),
  "2 traits: lnCF + lnLF" = intersect(c("lnCF","lnLF"), names(y_obs_demo)),
  "3 traits: lnCF + lnLF + lnCH" = intersect(c("lnCF","lnCH","lnLF"), names(y_obs_demo)),
  "4 traits: lnCF + lnCH + lnLF + lnLH" = intersect(trait_levels, names(y_obs_demo))
)
# Drop labels that end up empty because the specimen lacks some traits
incremental_orders <- incremental_orders[sapply(incremental_orders, length) >= 0]

# Compute posteriors via Laplace for each subset
posts_demo <- lapply(names(incremental_orders), function(lbl) {
  ids <- incremental_orders[[lbl]] # Trait IDs to use for this subset.
  if (length(ids) == 0) { # When there are no traits (prior‑only case)…
    list(mean = fit_demo$lnM_mean, var = fit_demo$lnM_sd^2, converged = TRUE) 
    # …return the prior on lnM (Normal with training mean & variance).
  } else { # Otherwise, run the Laplace posterior in lnM given those traits:
    predict_mass_laplace(
      y_obs     = y_obs_demo,
      trait_ids = ids,
      prior_mean = fit_demo$lnM_mean,
      prior_sd   = fit_demo$lnM_sd
    ) # Build the approximate Gaussian posterior for lnM using the Gaussianized trait predictive factors and the fold’s prior.
  }
})
names(posts_demo) <- names(incremental_orders)

# Common grid and curves
all_means <- sapply(posts_demo, `[[`, "mean")
all_sds   <- sqrt(sapply(posts_demo, `[[`, "var"))
grid <- seq(min(all_means - 5*all_sds), max(all_means + 5*all_sds), length.out = 800)

curves_demo <- dplyr::bind_rows(lapply(names(posts_demo), function(lbl) {
  df <- laplace_curve(posts_demo[[lbl]], grid = grid, width = 5)
  df$label <- lbl
  df
}))

# Styling
levels_demo <- names(posts_demo)
curves_demo$label <- factor(curves_demo$label, levels = levels_demo)

linetypes_demo <- setNames(ifelse(levels_demo=="Prior only","dashed","solid"), levels_demo)
colors_demo    <- setNames(RColorBrewer::brewer.pal(n = max(3,length(levels_demo)), "Dark2")[seq_along(levels_demo)], levels_demo)
if ("Prior only" %in% names(colors_demo)) colors_demo["Prior only"] <- "grey40"

# Plot
ggplot(curves_demo, aes(x = m, y = dens, color = label, linetype = label)) +
  geom_line(size = 1.2) +
  geom_vline(xintercept = lnM_true, color = "black", linewidth = 1.0, alpha = 0.85) +
  scale_color_manual(values = colors_demo) +
  scale_linetype_manual(values = linetypes_demo) +
  labs(title = sprintf("%s: (Fold %d test)", demo_sp, k_demo),
       x = "ln(Mass in grams)", y = "Density", color = "Information used", linetype = "Information used") +
  xlim(5,15)+
  theme_minimal(base_size = 13) +
  theme(legend.position = "right")

# Optional numeric summary
summ_demo <- dplyr::bind_rows(lapply(names(posts_demo), function(lbl) {
  mhat <- posts_demo[[lbl]]$mean
  vhat <- posts_demo[[lbl]]$var
  tibble::tibble(
    label    = lbl,
    lnM_mean = mhat,
    lnM_sd   = sqrt(vhat),
    M_median = exp(mhat),
    M_q025   = exp(mhat - 1.96*sqrt(vhat)),
    M_q975   = exp(mhat + 1.96*sqrt(vhat))
  )
}))
print(summ_demo)


#------------------------------------------------------------
# B) K-fold cross-validation (species-wise)
#------------------------------------------------------------
cv_list <- vector("list", K)

for (k in 1:K) {
  test_sp   <- species_ids[fold_id == k]
  train_sp  <- setdiff(species_ids, test_sp)
  fit_k     <- refit_conjugate(train_sp)
  post_k    <- fit_k$post_list
  
  # Build predictions for each held-out species
  preds_k <- lapply(test_sp, function(spi) {
    rows_i <- train_rows %>% filter(Species == spi)
    y_obs_i <- setNames(rows_i$val, rows_i$trait)
    y_obs_i <- y_obs_i[names(y_obs_i) %in% names(post_k)]
    # Skip if no usable traits
    if (length(y_obs_i) == 0) return(NULL)
    res <- predict_mass_laplace(
      y_obs      = y_obs_i,
      trait_ids  = names(y_obs_i),
      prior_mean = fit_k$lnM_mean,
      prior_sd   = fit_k$lnM_sd*10
    )
    tibble::tibble(
      Species = spi,
      true_lnM = unique(rows_i$lnM),
      pred_lnM = res$mean,
      pred_sd  = sqrt(res$var),
      n_traits = length(y_obs_i)
    )
  })
  cv_list[[k]] <- bind_rows(Filter(Negate(is.null), preds_k)) %>% mutate(fold = k)
}

cv_df <- bind_rows(cv_list)

# Metrics
cv_df <- cv_df %>%
  mutate(error = pred_lnM - true_lnM,
         abs_error = abs(error),
         cover95 = (true_lnM >= pred_lnM - 1.96*pred_sd) & (true_lnM <= pred_lnM + 1.96*pred_sd))

# RMSE = how close predicted lnM is to true lnM.
rmse <- sqrt(mean(cv_df$error^2, na.rm=TRUE))
mae  <- mean(cv_df$abs_error, na.rm=TRUE)
# Coverage = do the 90% intervals really cover ~90% of true masses?
coverage95 <- mean(cv_df$cover95, na.rm=TRUE)

cat(sprintf("Latent-mass CV — RMSE: %.4f, MAE: %.4f, 90%% coverage: %.3f\n", rmse, mae, coverage95))

# (Optional) by number of traits actually preserved
cv_by_nt <- cv_df %>%
  group_by(n_traits) %>%
  summarise(
    n = n(),
    RMSE = sqrt(mean(error^2)),
    MAE  = mean(abs_error),
    cov95 = mean(cover95),
    .groups = "drop"
  )
print(cv_by_nt)

# (Optional) diagnostic plot
ggplot(cv_df, aes(x = true_lnM, y = pred_lnM)) +
  geom_point(aes(size = n_traits), alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype="dashed") +
  labs(title = "K-fold CV: observed vs predicted ln(M)",
       x = "Observed ln(M)", y = "Predicted ln(M)", size = "#traits") +
  theme_minimal()

#-------------------------------------------------------------------------------
# compare to CFH prediction
# Comparability: lets you compare different models on the same resampling splits.
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# Baseline model: CF+H regression  (uses explicit dplyr::select to avoid MASS clash)
#-------------------------------------------------------------------------------

cfh_preds <- vector("list", K)

for (k in 1:K) {
  test_sp   <- species_ids[fold_id == k]
  train_sp  <- setdiff(species_ids, test_sp)
  
  # Training fold: require both CF and CH
  wide_train <- train_rows %>%
    dplyr::filter(Species %in% train_sp) %>%
    dplyr::select(Species, trait, val, lnM) %>%
    tidyr::pivot_wider(names_from = trait, values_from = val) %>%
    dplyr::filter(is.finite(lnM), is.finite(lnCF), is.finite(lnCH)) %>%
    dplyr::mutate(CFH = log(exp(lnCF) + exp(lnCH)))   # ln(CF + CH) in mm
  
  if (nrow(wide_train) < 5) {  # skip if too few training points
    cfh_preds[[k]] <- tibble::tibble()
    next
  }
  
  lm_cfh <- lm(lnM ~ CFH + lnLF + lnLH, data = wide_train)
  sigma_hat <- summary(lm_cfh)$sigma  # residual SD for predictive intervals
  
  # Test fold: require both CF and CH
  wide_test <- train_rows %>%
    dplyr::filter(Species %in% test_sp) %>%
    dplyr::select(Species, trait, val, lnM) %>%
    tidyr::pivot_wider(names_from = trait, values_from = val) %>%
    dplyr::filter(is.finite(lnM), is.finite(lnCF), is.finite(lnCH)) %>%
    dplyr::mutate(CFH = log(exp(lnCF) + exp(lnCH)))
  
  if (nrow(wide_test) == 0) { 
    cfh_preds[[k]] <- tibble::tibble()
    next 
  }
  
  pr <- predict(lm_cfh, newdata = wide_test, se.fit = TRUE)
  # predictive SD ≈ sqrt( se(fit)^2 + residual_sigma^2 )
  pred_sd_pred <- sqrt(drop(pr$se.fit)^2 + sigma_hat^2)
  
  cfh_preds[[k]] <- tibble::tibble(
    Species  = wide_test$Species,
    true_lnM = wide_test$lnM,
    pred_lnM = as.numeric(pr$fit),
    pred_sd  = as.numeric(pred_sd_pred),
    fold     = k
  )
}

cfh_df <- dplyr::bind_rows(cfh_preds)

#-------------------------------------------------------------------------------
# Metrics for CFH baseline (95% coverage on Normal approx)
#-------------------------------------------------------------------------------
cfh_df <- cfh_df %>%
  dplyr::mutate(
    error     = pred_lnM - true_lnM,
    abs_error = abs(error),
    cover95   = (true_lnM >= pred_lnM - 1.96*pred_sd) &
      (true_lnM <= pred_lnM + 1.96*pred_sd)
  )

rmse_cfh <- sqrt(mean(cfh_df$error^2, na.rm=TRUE))
mae_cfh  <- mean(cfh_df$abs_error, na.rm=TRUE)
coverage95_cfh <- mean(cfh_df$cover95, na.rm=TRUE)

cat(sprintf("CF+H baseline — RMSE: %.4f, MAE: %.4f, 95%% coverage: %.3f\n",
            rmse_cfh, mae_cfh, coverage95_cfh))

#-------------------------------------------------------------------------------
# Side-by-side comparison (note: your latent-mass block calls this 'coverage90'
# but it’s based on 1.96; align the name or the z-quantile for consistency)
#-------------------------------------------------------------------------------
metrics_compare <- tibble::tibble(
  Model = c("Latent-mass", "CF+H regression"),
  RMSE  = c(rmse, rmse_cfh),
  MAE   = c(mae, mae_cfh),
  Coverage95 = c(coverage95, coverage95_cfh)  
)
print(metrics_compare)

# Optional combined scatter
ggplot() +
  geom_point(data = cv_df,
             aes(x = true_lnM, y = pred_lnM, color = "Latent-mass"), alpha = 0.6) +
  geom_point(data = cfh_df,
             aes(x = true_lnM, y = pred_lnM, color = "CF+H"), alpha = 0.6, shape = 17) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(title = "Cross-validation: latent-mass vs CF+H baseline",
       x = "Observed ln(M)", y = "Predicted ln(M)", color = "Model") +
  theme_minimal()

