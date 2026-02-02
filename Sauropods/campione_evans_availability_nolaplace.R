#################################################
# Title: Hierarchical mass allometry (conjugate)
# Description: Mass-from-bones with latent lnM,
# using closed-form Normal–Inverse–Gamma updates
# + MVN residual correlation + Laplace & Grid.
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
library(bayesplot)
library(MASSTIMATE)
library(mvtnorm)    # MVN density
library(Matrix)     # nearPD
color_scheme_set("blue")

# ---- Controls / toggles ----
prior_scale_demo <- 10   # weak prior sd = fold sd * this factor (demo)
prior_scale_cv   <- 10   # weak prior sd in CV
use_grid_in_cv   <- FALSE  # set TRUE to use grid posterior in CV
grid_width       <- 6
grid_n           <- 2001

# ---- Traits (J = 4) ----
trait_levels <- c("lnCF","lnCH","lnLF","lnLH")

# ---- Prior centres (Alexander 1979), log scale ----
ah_CF <- log(0.43); ah_CH <- log(0.35); ah_LF <- log(5.24); ah_LH <- log(4.24)
bh_CF <- 0.36;     bh_CH <- 0.38;       bh_LF <- 0.36;       bh_LH <- 0.36
alpha_hat <- c(ah_CF, ah_CH, ah_LF, ah_LH)
beta_hat  <- c(bh_CF, bh_CH, bh_LF, bh_LH)

# ---- Conjugate NIG hyperparameters ----
sd_alpha <- 0.5
sd_beta  <- 0.1
V0_list  <- vector("list", length(trait_levels))
mu0_list <- vector("list", length(trait_levels))
for (j in seq_along(trait_levels)) {
  V0_list[[j]]  <- diag(c(sd_alpha^2, sd_beta^2))
  mu0_list[[j]] <- c(alpha_hat[j], beta_hat[j])
}
a0 <- 2.1
b0 <- 1.1

# ---- Inv-Gamma sampling ----
rinvgamma <- function(n, shape, rate) 1 / rgamma(n, shape = shape, rate = rate)

# ---- Conjugate update per trait ----
nig_posterior <- function(y, X, mu0, V0, a0, b0) {
  XtX <- crossprod(X); XtY <- crossprod(X, y)
  V0i <- solve(V0)
  Vni <- V0i + XtX
  Vn  <- solve(Vni)
  mu_n <- Vn %*% (V0i %*% mu0 + XtY)
  a_n <- a0 + length(y)/2
  quad <- as.numeric(crossprod(y, y) + crossprod(mu0, V0i %*% mu0) - crossprod(mu_n, Vni %*% mu_n))
  b_n <- b0 + 0.5 * quad
  list(mu_n = as.numeric(mu_n), Vn = Vn, a_n = a_n, b_n = b_n, n = length(y))
}

# ---- Posterior draws (optional) ----
sample_from_nig <- function(n_samp, mu_n, Vn, a_n, b_n) {
  sig2  <- rinvgamma(n_samp, shape = a_n, rate = b_n)
  theta <- t(apply(cbind(sig2), 1, function(s2) as.numeric(mvrnorm(1, mu = mu_n, Sigma = s2 * Vn))))
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
  dplyr::select(Species, lnM, lnCF, lnCH, lnLF, lnLH) %>%
  pivot_longer(cols = starts_with("ln"),
               names_to = "trait", values_to = "val") %>%
  mutate(
    is_mass  = (trait == "lnM"),
    trait_id = match(trait, trait_levels)
  )

train_rows <- extants_long %>%
  filter(!is_mass, trait %in% trait_levels) %>%
  group_by(Species) %>%
  mutate(lnM = extants_long$val[match(paste0(Species,"_lnM"),
                                      paste0(extants_long$Species,"_",extants_long$trait))]) %>%
  ungroup() %>%
  filter(is.finite(lnM), is.finite(val))

trait_train <- lapply(seq_along(trait_levels), function(j) {
  dfj <- train_rows %>% filter(trait_id == j)
  list(y = dfj$val, X = cbind(1, dfj$lnM))
})

# ---------------------------------------
# TRAIN on full data (for summaries/ppc if needed)
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

# ============================================================
# MVN likelihood helpers (correlated residuals across traits)
# ============================================================
normalize_trait_ids <- function(x, allowed) {
  x_chr <- trimws(as.character(x)); bad <- setdiff(x_chr, allowed)
  if (length(bad))
    stop("Unknown trait id(s): ", paste(bad, collapse = ", "),
         "\nAllowed: ", paste(allowed, collapse = ", "))
  x_chr
}

# Per-trait mean/var at m
gaussian_trait_params <- function(post, m) {
  x  <- c(1, m)
  mu <- as.numeric(x %*% post$mu_n)
  cfac <- (post$b_n / post$a_n); Vn <- post$Vn
  q  <- as.numeric(1 + x %*% Vn %*% x)
  v  <- cfac * q
  list(mu = mu, v = v)
}

# Build MVN mean vector and covariance Σ(m) = D(m) R D(m)
mv_params_at_m <- function(post_list_local, trait_ids, m, R) {
  pjL <- post_list_local[trait_ids]
  mu  <- vapply(pjL, function(pj) gaussian_trait_params(pj, m)$mu, numeric(1))
  v   <- vapply(pjL, function(pj) gaussian_trait_params(pj, m)$v,  numeric(1))
  D   <- diag(sqrt(v))
  Rsub <- R[trait_ids, trait_ids, drop = FALSE]
  S   <- D %*% Rsub %*% D
  list(mu = mu, Sigma = S)
}

logpost_m <- function(m, y_vec, trait_ids, prior_mean, prior_sd, post_list_local, R) {
  lp <- dnorm(m, prior_mean, prior_sd, log = TRUE)
  prm <- mv_params_at_m(post_list_local, trait_ids, m, R)
  lp + mvtnorm::dmvnorm(y_vec, mean = prm$mu, sigma = prm$Sigma, log = TRUE)
}

# ---- numeric derivatives Laplace (robust in 1D) ----
trapz <- function(x, y) sum((x[-1]-x[-length(x)]) * (y[-1]+y[-length(y)])/2)

posterior_m_laplace <- function(y_obs, trait_ids, prior_mean, prior_sd,
                                post_list_local, R,
                                m_init = prior_mean, maxit = 60, tol = 1e-8, step_shrink = 0.5) {
  trait_ids <- normalize_trait_ids(trait_ids, names(post_list_local))
  y_vec <- as.numeric(y_obs[trait_ids])
  
  ell <- function(m) logpost_m(m, y_vec, trait_ids, prior_mean, prior_sd, post_list_local, R)
  step_for <- function(m) { eps <- .Machine$double.eps^(1/3); h <- eps*(abs(m)+1); max(h, 1e-6) }
  ell_prime  <- function(m){ h <- step_for(m); (ell(m+h) - ell(m-h))/(2*h) }
  ell_second <- function(m){ h <- step_for(m); (ell(m+h) - 2*ell(m) + ell(m-h))/(h^2) }
  
  m <- m_init
  for (k in 1:maxit) {
    g <- ell_prime(m); H <- ell_second(m)
    if (!is.finite(H) || H >= 0) { m <- m + rnorm(1,0,1e-4); next }
    step <- -g/H; m_new <- m + step
    tries <- 0
    while (!is.finite(ell(m_new)) || (ell(m_new) < ell(m) && tries < 6)) {
      step <- step * step_shrink; m_new <- m + step; tries <- tries + 1
    }
    if (abs(m_new - m) < tol) { m <- m_new; break }
    m <- m_new
  }
  Hhat <- ell_second(m)
  Vhat <- if (is.finite(Hhat) && Hhat < 0) -1/Hhat else NA_real_
  list(mean = as.numeric(m), var = as.numeric(Vhat), converged = is.finite(Vhat) && Vhat > 0)
}

posterior_m_grid <- function(y_obs, trait_ids, prior_mean, prior_sd,
                             post_list_local, R,
                             grid_center, grid_sd, width = 6, n = 2001) {
  trait_ids <- normalize_trait_ids(trait_ids, names(post_list_local))
  y_vec <- as.numeric(y_obs[trait_ids])
  m_grid <- seq(grid_center - width*grid_sd, grid_center + width*grid_sd, length.out = n)
  logp <- vapply(m_grid, function(m) {
    logpost_m(m, y_vec, trait_ids, prior_mean, prior_sd, post_list_local, R)
  }, numeric(1))
  logp  <- logp - max(logp)
  w     <- exp(logp)
  Z     <- trapz(m_grid, w); w_norm <- w/Z
  mean_m <- trapz(m_grid, m_grid * w_norm)
  var_m  <- trapz(m_grid, (m_grid - mean_m)^2 * w_norm)
  list(mean = mean_m, var = var_m, grid = m_grid, dens = w_norm)
}

# ---------------------------------------
# Fold refit with residual correlation R
# ---------------------------------------
refit_conjugate <- function(train_species) {
  train_rows_k <- train_rows %>% dplyr::filter(Species %in% train_species)
  
  # Per‑trait posteriors
  trait_train_k <- lapply(seq_along(trait_levels), function(j) {
    dfj <- train_rows_k %>% dplyr::filter(trait_id == j)
    list(y = dfj$val, X = cbind(1, dfj$lnM), lnM = dfj$lnM)
  })
  post_list_k <- lapply(seq_along(trait_levels), function(j) {
    nig_posterior(trait_train_k[[j]]$y, trait_train_k[[j]]$X,
                  mu0_list[[j]], V0_list[[j]], a0, b0)
  })
  names(post_list_k) <- trait_levels
  
  # Residual correlation across traits (posterior-mean fits)
  wide_train <- train_rows_k %>%
    dplyr::select(Species, trait, val, lnM) %>%
    tidyr::pivot_wider(names_from = trait, values_from = val)
  
  for (tj in trait_levels) {
    pj <- post_list_k[[tj]]
    wide_train[[paste0(tj,"_fit")]] <- as.numeric(pj$mu_n[1] + pj$mu_n[2] * wide_train$lnM)
    wide_train[[paste0(tj,"_res")]] <- wide_train[[tj]] - wide_train[[paste0(tj,"_fit")]]
  }
  
  res_mat <- wide_train %>% dplyr::select(dplyr::ends_with("_res")) %>% as.matrix()
  ok <- stats::complete.cases(res_mat)
  Rhat <- try(stats::cor(res_mat[ok,, drop = FALSE], use = "pairwise.complete.obs"), silent = TRUE)
  
  if (inherits(Rhat, "try-error") || any(!is.finite(Rhat))) {
    Rhat <- diag(length(trait_levels))
    colnames(Rhat) <- rownames(Rhat) <- paste0(trait_levels,"_res")
    warning("R estimation failed; using identity.")
  }
  # Map to trait order
  colnames(Rhat) <- gsub("_res$","", colnames(Rhat))
  rownames(Rhat) <- gsub("_res$","", rownames(Rhat))
  Rfull <- matrix(0, length(trait_levels), length(trait_levels),
                  dimnames = list(trait_levels, trait_levels))
  common <- intersect(rownames(Rhat), trait_levels)
  Rfull[common, common] <- Rhat[common, common]
  diag(Rfull) <- 1
  R_pd <- as.matrix(nearPD(Rfull, corr = TRUE)$mat)
  
  list(
    post_list = post_list_k,
    lnM_mean  = mean(train_rows_k$lnM),
    lnM_sd    = sd(train_rows_k$lnM),
    R         = R_pd
  )
}

# ============================================================
# DEMO: Laplace vs Grid on one specimen (MVN likelihood)
# ============================================================
set.seed(42)
K <- 10
species_ids <- unique(train_rows$Species)
fold_id <- sample(rep(1:K, length.out = length(species_ids)))

k_demo <- 1
test_species  <- species_ids[fold_id == k_demo]
train_species <- setdiff(species_ids, test_species)
fit_demo <- refit_conjugate(train_species)

demo_sp   <- test_species[1]
demo_rows <- train_rows %>% dplyr::filter(Species == demo_sp)
lnM_true  <- unique(demo_rows$lnM)

y_obs_demo <- setNames(demo_rows$val, demo_rows$trait)
y_obs_demo <- y_obs_demo[names(y_obs_demo) %in% trait_levels]
ids_demo   <- names(y_obs_demo)

prior_mean_demo <- fit_demo$lnM_mean
prior_sd_demo   <- fit_demo$lnM_sd * prior_scale_demo

# Laplace
lap_demo <- posterior_m_laplace(
  y_obs = y_obs_demo, trait_ids = ids_demo,
  prior_mean = prior_mean_demo, prior_sd = prior_sd_demo,
  post_list_local = fit_demo$post_list, R = fit_demo$R
)

# Grid centered at Laplace
grid_demo <- posterior_m_grid(
  y_obs = y_obs_demo, trait_ids = ids_demo,
  prior_mean = prior_mean_demo, prior_sd = prior_sd_demo,
  post_list_local = fit_demo$post_list, R = fit_demo$R,
  grid_center = lap_demo$mean, grid_sd = sqrt(lap_demo$var),
  width = grid_width, n = grid_n
)

print(tibble::tibble(
  method = c("Laplace","Grid"),
  mean   = c(lap_demo$mean, grid_demo$mean),
  sd     = c(sqrt(lap_demo$var), sqrt(grid_demo$var))
))

curve_df <- rbind(
  tibble::tibble(m = grid_demo$grid,
                 dens = dnorm(grid_demo$grid, mean = lap_demo$mean, sd = sqrt(lap_demo$var)),
                 label = "Laplace Normal"),
  tibble::tibble(m = grid_demo$grid,
                 dens = grid_demo$dens,
                 label = "Grid (normalized)")
)
ggplot(curve_df, aes(m, dens, color = label)) +
  geom_line(size = 1.1) +
  geom_vline(xintercept = lnM_true, color = "black", linetype = "dashed") +
  labs(title = sprintf("%s: posterior ln(M) — Laplace vs Grid (MVN)", demo_sp),
       x = "ln(Mass in grams)", y = "Density", color = "Method") +
  theme_minimal()

# ============================================================
# B) K-fold cross-validation (species-wise), MVN likelihood
# ============================================================
cv_list <- vector("list", K)

for (k in 1:K) {
  test_sp   <- species_ids[fold_id == k]
  train_sp  <- setdiff(species_ids, test_sp)
  fit_k     <- refit_conjugate(train_sp)
  
  preds_k <- lapply(test_sp, function(spi) {
    rows_i  <- train_rows %>% dplyr::filter(Species == spi)
    y_obs_i <- setNames(rows_i$val, rows_i$trait)
    y_obs_i <- y_obs_i[names(y_obs_i) %in% names(fit_k$post_list)]
    if (length(y_obs_i) == 0) return(NULL)
    
    prior_mean_k <- fit_k$lnM_mean
    prior_sd_k   <- fit_k$lnM_sd * prior_scale_cv
    
    # Laplace or Grid
    res_lap <- posterior_m_laplace(
      y_obs = y_obs_i, trait_ids = names(y_obs_i),
      prior_mean = prior_mean_k, prior_sd = prior_sd_k,
      post_list_local = fit_k$post_list, R = fit_k$R
    )
    
    if (use_grid_in_cv) {
      res_grid <- posterior_m_grid(
        y_obs = y_obs_i, trait_ids = names(y_obs_i),
        prior_mean = prior_mean_k, prior_sd = prior_sd_k,
        post_list_local = fit_k$post_list, R = fit_k$R,
        grid_center = res_lap$mean, grid_sd = sqrt(res_lap$var),
        width = grid_width, n = grid_n
      )
      use_mean <- res_grid$mean; use_var <- res_grid$var
    } else {
      use_mean <- res_lap$mean; use_var <- res_lap$var
    }
    
    tibble::tibble(
      Species  = spi,
      true_lnM = unique(rows_i$lnM),
      pred_lnM = use_mean,
      pred_sd  = sqrt(use_var),
      n_traits = length(y_obs_i),
      fold     = k
    )
  })
  
  cv_list[[k]] <- dplyr::bind_rows(Filter(Negate(is.null), preds_k))
}

cv_df <- dplyr::bind_rows(cv_list)

# Metrics
cv_df <- cv_df %>%
  mutate(error = pred_lnM - true_lnM,
         abs_error = abs(error),
         cover95 = (true_lnM >= pred_lnM - 1.96*pred_sd) & (true_lnM <= pred_lnM + 1.96*pred_sd))

rmse <- sqrt(mean(cv_df$error^2, na.rm=TRUE))
mae  <- mean(cv_df$abs_error, na.rm=TRUE)
coverage95 <- mean(cv_df$cover95, na.rm=TRUE)

cat(sprintf("Latent-mass (MVN) — RMSE: %.4f, MAE: %.4f, 95%% coverage: %.3f\n", rmse, mae, coverage95))

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

ggplot(cv_df, aes(x = true_lnM, y = pred_lnM)) +
  geom_point(aes(size = n_traits), alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype="dashed") +
  labs(title = "K-fold CV (MVN): observed vs predicted ln(M)",
       x = "Observed ln(M)", y = "Predicted ln(M)", size = "#traits") +
  theme_minimal()

# -----------------------------------------------------------
# CF+H baseline on the same folds
# -----------------------------------------------------------
cfh_preds <- vector("list", K)

for (k in 1:K) {
  test_sp   <- species_ids[fold_id == k]
  train_sp  <- setdiff(species_ids, test_sp)
  
  wide_train <- train_rows %>%
    dplyr::filter(Species %in% train_sp) %>%
    dplyr::select(Species, trait, val, lnM) %>%
    tidyr::pivot_wider(names_from = trait, values_from = val) %>%
    dplyr::filter(is.finite(lnM), is.finite(lnCF), is.finite(lnCH)) %>%
    dplyr::mutate(CFH = log(exp(lnCF) + exp(lnCH)))
  
  if (nrow(wide_train) < 5) { cfh_preds[[k]] <- tibble::tibble(); next }
  
  lm_cfh <- lm(lnM ~ CFH, data = wide_train)
  sigma_hat <- summary(lm_cfh)$sigma
  
  wide_test <- train_rows %>%
    dplyr::filter(Species %in% test_sp) %>%
    dplyr::select(Species, trait, val, lnM) %>%
    tidyr::pivot_wider(names_from = trait, values_from = val) %>%
    dplyr::filter(is.finite(lnM), is.finite(lnCF), is.finite(lnCH)) %>%
    dplyr::mutate(CFH = log(exp(lnCF) + exp(lnCH)))
  
  if (nrow(wide_test) == 0) { cfh_preds[[k]] <- tibble::tibble(); next }
  
  pr <- predict(lm_cfh, newdata = wide_test, se.fit = TRUE)
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

cfh_df <- cfh_df %>%
  mutate(error = pred_lnM - true_lnM,
         abs_error = abs(error),
         cover95 = (true_lnM >= pred_lnM - 1.96*pred_sd) & (true_lnM <= pred_lnM + 1.96*pred_sd))

rmse_cfh <- sqrt(mean(cfh_df$error^2, na.rm=TRUE))
mae_cfh  <- mean(cfh_df$abs_error, na.rm=TRUE)
coverage95_cfh <- mean(cfh_df$cover95, na.rm=TRUE)

cat(sprintf("CF+H baseline — RMSE: %.4f, MAE: %.4f, 95%% coverage: %.3f\n",
            rmse_cfh, mae_cfh, coverage95_cfh))

# Side-by-side
metrics_compare <- tibble::tibble(
  Model = c("Latent-mass (MVN)", "CF+H regression"),
  RMSE  = c(rmse, rmse_cfh),
  MAE   = c(mae, mae_cfh),
  Coverage95 = c(coverage95, coverage95_cfh)
)
print(metrics_compare)

ggplot() +
  geom_point(data = cv_df,
             aes(x = true_lnM, y = pred_lnM, color = "Latent-mass (MVN)"), alpha = 0.6) +
  geom_point(data = cfh_df,
             aes(x = true_lnM, y = pred_lnM, color = "CF+H"), alpha = 0.6, shape = 17) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(title = "Cross-validation: latent-mass (MVN) vs CF+H",
       x = "Observed ln(M)", y = "Predicted ln(M)", color = "Model") +
  theme_minimal()

