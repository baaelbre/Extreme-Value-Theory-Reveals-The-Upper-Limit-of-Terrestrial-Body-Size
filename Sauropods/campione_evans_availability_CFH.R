#################################################
# Title: CFH–centric mass estimation (Bayesian + OLS)
# Description:
#  (1) OLS: lnM ~ lnCFH
#  (2) OLS: lnM ~ lnCFH + lnLF + lnLH
#  (3) Bayes: p(m|CFH) with weak priors when CFH available
#  (4) Missing CF/CH stress test:
#      - Non-Bayes: impute CF or CH via OLS (CF~CH / CH~CF) → build CFH
#      - Bayes: condition on available channels (CF/CH/LF/LH), no point imputation
#  (5) Matched-sample equivalence check:
#      Bayes (closed-form, weak prior, plug-in variance) on {CFH,LF,LH}
#      vs OLS {CFH,LF,LH} with CFH observed (no imputation)
# Author: Bastiaan Aelbrecht (rewrite)
# Date: 24/08/2025
#################################################

# ---- Packages ----
library(MASS)       # mvrnorm
library(dplyr)
library(tidyr)
library(tibble)
library(purrr)
library(ggplot2)
library(MASSTIMATE)

# ============================================================
# 0) Controls
# ============================================================
set.seed(42)
K <- 5                          # species-wise folds
p_drop_circ <- 0.30             # drop exactly one of CF/CH with prob p (0 = disable)
report_models <- c("OLS: lnM ~ lnCFH",
                   "OLS: lnM ~ lnCFH + lnLF + lnLH",
                   "OLS: lnM ~ lnCFH + lnLF + lnLH (obs-only CFH)",
                   "Bayes: CFH available (p(m|CFH))",
                   "Bayes: CFH+LF+LH (weak prior, plug-in, matched)",
                   "Bayes: CFH missing (via CF/CH/LF/LH)")

# ============================================================
# 1) Priors & helpers (NIG calibration per trait)
# ============================================================

# Traits used for calibration (original four)
trait_levels <- c("lnCF","lnCH","lnLF","lnLH")

# Prior centres (Alexander 1979), log scale
ah_CF <- log(0.43); ah_CH <- log(0.35); ah_LF <- log(5.24); ah_LH <- log(4.24)
bh_CF <- 0.36;      bh_CH <- 0.38;      bh_LF <- 0.36;       bh_LH <- 0.36
alpha_hat <- c(lnCF=ah_CF, lnCH=ah_CH, lnLF=ah_LF, lnLH=ah_LH)
beta_hat  <- c(lnCF=bh_CF, lnCH=bh_CH, lnLF=bh_LF, lnLH=bh_LH)

# Conjugate NIG hyperparameters
# theta_j | sigma_j^2 ~ N( mu0_j, sigma_j^2 * V0_j ),  sigma_j^2 ~ IG(a0, b0)
sd_alpha <- 0.8
sd_beta  <- 0.4
V0_list <- lapply(trait_levels, function(j) diag(c(sd_alpha^2, sd_beta^2)))
names(V0_list) <- trait_levels
mu0_list <- list(
  lnCF = c(alpha_hat["lnCF"], beta_hat["lnCF"]),
  lnCH = c(alpha_hat["lnCH"], beta_hat["lnCH"]),
  lnLF = c(alpha_hat["lnLF"], beta_hat["lnLF"]),
  lnLH = c(alpha_hat["lnLH"], beta_hat["lnLH"])
)
a0 <- 1.1
b0 <- 2

# ===== Utilities =====
logsumexp2 <- function(a, b) { m <- pmax(a, b); m + log(exp(a - m) + exp(b - m)) }
rinvgamma <- function(n, shape, rate) 1 / rgamma(n, shape=shape, rate=rate)

nig_posterior <- function(y, X, mu0, V0, a0, b0) {
  XtX <- crossprod(X); XtY <- crossprod(X, y)
  V0_inv <- solve(V0); Vn_inv <- V0_inv + XtX; Vn <- solve(Vn_inv)
  mu_n <- Vn %*% (V0_inv %*% mu0 + XtY)
  a_n <- a0 + length(y)/2
  quad <- as.numeric(crossprod(y, y) + crossprod(mu0, V0_inv %*% mu0) - crossprod(mu_n, Vn_inv %*% mu_n))
  b_n <- b0 + 0.5 * quad
  list(mu_n = as.numeric(mu_n), Vn = Vn, a_n = a_n, b_n = b_n, n = length(y))
}

# Predictive factor p(y_j | m) ≈ N(mu_j(m), v_j(m)) and derivatives (full NIG predictive)
gaussian_trait_params <- function(post, m) {
  x  <- c(1, m); mu <- as.numeric(x %*% post$mu_n)
  cfac <- (post$b_n / post$a_n); Vn <- post$Vn
  q  <- as.numeric(1 + x %*% Vn %*% x)
  v  <- cfac * q
  v1 <- 2 * cfac * (Vn[1,2] + Vn[2,2] * m)
  v2 <- 2 * cfac * Vn[2,2]
  list(mu = mu, v = v, v1 = v1, v2 = v2)
}

# Plug-in (constant) variance for equivalence check: s_j^2 = E[sigma^2 | data]
post_sigma2_plugin <- function(post) {
  if (post$a_n > 1) post$b_n / (post$a_n - 1) else post$b_n / post$a_n
}

# Closed-form posterior under constant variances (weak prior -> ≈ OLS/GLS)
predict_mass_closed_form <- function(y_obs, trait_ids, prior_mean, prior_sd, post_dict) {
  stopifnot(length(trait_ids) > 0)
  if (!all(trait_ids %in% names(post_dict))) {
    stop("Missing trait(s) in post_dict: ", paste(setdiff(trait_ids, names(post_dict)), collapse=", "))
  }
  pjL <- post_dict[trait_ids]
  alpha_j <- vapply(pjL, function(p) p$mu_n[1], numeric(1))
  beta_j  <- vapply(pjL, function(p) p$mu_n[2], numeric(1))
  s2_j    <- vapply(pjL, function(p) post_sigma2_plugin(p), numeric(1))
  y_vec   <- as.numeric(y_obs[trait_ids])
  tau2 <- prior_sd^2
  prec <- 1/tau2 + sum((beta_j^2)/s2_j)
  mean_num <- (prior_mean/tau2) + sum(beta_j*(y_vec - alpha_j)/s2_j)
  mu_post <- mean_num / prec
  var_post <- 1/prec
  list(mean = as.numeric(mu_post), var = as.numeric(var_post))
}

# Full Laplace with safe fallback (kept for CFH-only and missingness cases)
log_norm_factor <- function(y, mu, v) -0.5 * (log(2*pi*v) + (y - mu)^2 / v)
predict_mass_laplace <- function(y_obs, trait_ids, prior_mean, prior_sd, post_dict, search_sd) {
  stopifnot(length(trait_ids) > 0)
  if (!all(trait_ids %in% names(post_dict))) {
    stop("Missing trait(s) in post_dict: ", paste(setdiff(trait_ids, names(post_dict)), collapse=", "))
  }
  y_vec <- as.numeric(y_obs[trait_ids]); pjL <- post_dict[trait_ids]
  ell <- function(m) {
    lp <- dnorm(m, mean=prior_mean, sd=prior_sd, log=TRUE)
    for (j in seq_along(pjL)) {
      gp <- gaussian_trait_params(pjL[[j]], m)
      lp <- lp + log_norm_factor(y_vec[j], gp$mu, gp$v)
    }
    lp
  }
  rng <- c(prior_mean - 6*search_sd, prior_mean + 6*search_sd)
  opt <- optimize(function(mm) -ell(mm), interval=rng)
  m_mode <- opt$minimum
  h <- max(1e-4, 1e-3*search_sd)
  lpp <- (ell(m_mode + h) - 2*ell(m_mode) + ell(m_mode - h)) / (h^2)
  if (is.finite(lpp) && lpp < 0) {
    v_hat <- -1 / lpp
    return(list(mean = as.numeric(m_mode), var = as.numeric(v_hat), converged = TRUE))
  }
  # fallback: constant-variance closed-form using plug-in s2 at the mode
  pjL <- post_dict[trait_ids]  # recompute (scope)
  alpha_j <- vapply(pjL, function(p) p$mu_n[1], numeric(1))
  beta_j  <- vapply(pjL, function(p) p$mu_n[2], numeric(1))
  s2_j    <- vapply(pjL, function(p) post_sigma2_plugin(p), numeric(1))
  y_vec   <- as.numeric(y_obs[trait_ids])
  tau2 <- prior_sd^2
  prec <- 1/tau2 + sum((beta_j^2)/s2_j)
  mean_num <- (prior_mean/tau2) + sum(beta_j*(y_vec - alpha_j)/s2_j)
  mu_post <- mean_num / prec
  var_post <- 1/prec
  list(mean = as.numeric(mu_post), var = as.numeric(var_post), converged = FALSE)
}

# ---- SAFE BUILDERS (avoid NA trait names/ids) ----
compute_lncfh <- function(vals) {
  if (is.finite(vals["lnCF"]) && is.finite(vals["lnCH"])) return(logsumexp2(vals["lnCF"], vals["lnCH"]))
  NA_real_
}

build_y <- function(vals, post_dict, keep) {
  out <- numeric(0)
  if ("lnCFH" %in% keep) {
    lncfh <- compute_lncfh(vals)
    if (is.finite(lncfh)) out["lnCFH"] <- lncfh
  }
  if ("lnCF" %in% keep && is.finite(vals["lnCF"])) out["lnCF"] <- vals["lnCF"]
  if ("lnCH" %in% keep && is.finite(vals["lnCH"])) out["lnCH"] <- vals["lnCH"]
  if ("lnLF" %in% keep && is.finite(vals["lnLF"])) out["lnLF"] <- vals["lnLF"]
  if ("lnLH" %in% keep && is.finite(vals["lnLH"])) out["lnLH"] <- vals["lnLH"]
  trait_ids <- intersect(names(out), names(post_dict))
  list(y = out[trait_ids], trait_ids = trait_ids)
}

# ============================================================
# 2) Data preparation (extants → long & wide)
# ============================================================

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
  pivot_longer(cols = starts_with("ln"), names_to = "trait", values_to = "val")

# Wide table with one row per species (includes lnM and all traits)
wide_all <- extants_long %>%
  dplyr::select(Species, trait, val) %>%
  tidyr::pivot_wider(names_from = trait, values_from = val) %>%
  dplyr::mutate(
    lnCFH = ifelse(is.finite(lnCF) & is.finite(lnCH), log(exp(lnCF) + exp(lnCH)), NA_real_)
  )

# Per-trait training rows with lnM attached
train_rows <- extants_long %>%
  filter(trait %in% trait_levels) %>%
  left_join(wide_all %>% dplyr::select(Species, lnM), by = "Species") %>%
  filter(is.finite(lnM), is.finite(val))

# ============================================================
# 3) Global NIG posteriors for original four traits (reference)
# ============================================================

trait_train <- lapply(trait_levels, function(j) {
  dfj <- train_rows %>% filter(trait == j)
  list(y = dfj$val, X = cbind(1, dfj$lnM))
})
names(trait_train) <- trait_levels

post_list_global <- lapply(trait_levels, function(j) {
  nig_posterior(
    y   = trait_train[[j]]$y,
    X   = trait_train[[j]]$X,
    mu0 = mu0_list[[j]],
    V0  = V0_list[[j]],
    a0  = a0, b0 = b0
  )
})
names(post_list_global) <- trait_levels

# ============================================================
# 4) Fold refit: NIG for original traits + NIG for CFH + OLS modules
# ============================================================

refit_fold <- function(train_species) {
  tr_long <- train_rows %>% filter(Species %in% train_species)
  tr_wide <- wide_all  %>% filter(Species %in% train_species)
  
  # NIG for original 4 traits
  post_list_k <- lapply(trait_levels, function(j) {
    dfj <- tr_long %>% filter(trait == j)
    nig_posterior(
      y   = dfj$val,
      X   = cbind(1, dfj$lnM),
      mu0 = mu0_list[[j]],
      V0  = V0_list[[j]],
      a0  = a0, b0 = b0
    )
  })
  names(post_list_k) <- trait_levels
  
  # CFH NIG posterior (always fit if >=3 rows)
  cfh_rows <- tr_wide %>% filter(is.finite(lnCFH), is.finite(lnM))
  if (nrow(cfh_rows) < 3) stop("Too few CFH rows in training fold.")
  y_cfh <- cfh_rows$lnCFH
  X_cfh <- cbind(1, cfh_rows$lnM)
  mu0_cfh <- c(mean(alpha_hat[c("lnCF","lnCH")]), mean(beta_hat[c("lnCF","lnCH")]))
  V0_cfh  <- diag(c(sd_alpha^2, sd_beta^2))
  post_cfh <- nig_posterior(y_cfh, X_cfh, mu0=mu0_cfh, V0=V0_cfh, a0=a0, b0=b0)
  post_list_k[["lnCFH"]] <- post_cfh
  
  # OLS baselines
  dat_cfh <- tr_wide %>% filter(is.finite(lnM), is.finite(lnCFH))
  lm_cfh  <- if (nrow(dat_cfh) >= 5) lm(lnM ~ lnCFH, data=dat_cfh) else NULL
  
  dat_cfh_len <- tr_wide %>% filter(is.finite(lnM), is.finite(lnCFH),
                                    is.finite(lnLF), is.finite(lnLH))
  lm_cfh_len  <- if (nrow(dat_cfh_len) >= 8) lm(lnM ~ lnCFH + lnLF + lnLH, data=dat_cfh_len) else NULL
  
  # Non-Bayes imputation lines (CF~CH and CH~CF)
  train_cc <- tr_wide %>% filter(is.finite(lnCF), is.finite(lnCH))
  cf_from_ch <- if (nrow(train_cc) >= 5) lm(lnCF ~ lnCH, data=train_cc) else NULL
  ch_from_cf <- if (nrow(train_cc) >= 5) lm(lnCH ~ lnCF, data=train_cc) else NULL
  
  # Prior: weak (won't dominate), search_sd = fold sd for mode bracketing
  fold_mean <- mean(tr_wide$lnM, na.rm=TRUE)
  fold_sd   <- sd(tr_wide$lnM, na.rm=TRUE)
  prior_sd  <- max(5, fold_sd*10)   # even weaker to encourage equivalence
  
  list(
    post_list   = post_list_k,          # includes lnCFH!
    prior_mean  = fold_mean,
    prior_sd    = prior_sd,
    search_sd   = max(1, fold_sd),
    lm_cfh      = lm_cfh,
    lm_cfh_len  = lm_cfh_len,
    imp_cf_from_ch = cf_from_ch,
    imp_ch_from_cf = ch_from_cf
  )
}

# ============================================================
# 5) Test-time masking: drop exactly one of CF/CH with prob p
# ============================================================
mask_one_circ <- function(vals_named, p) {
  CF <- vals_named["lnCF"]; CH <- vals_named["lnCH"]
  if (is.finite(CF) && is.finite(CH) && runif(1) < p) {
    if (runif(1) < 0.5) vals_named["lnCF"] <- NA_real_ else vals_named["lnCH"] <- NA_real_
  }
  vals_named
}

# ============================================================
# 6) CV loop
# ============================================================

species_ids <- unique(wide_all$Species)
fold_id <- sample(rep(1:K, length.out=length(species_ids)))
out_rows <- list()

for (k in seq_len(K)) {
  test_sp   <- species_ids[fold_id == k]
  train_sp  <- setdiff(species_ids, test_sp)
  
  fit_k <- refit_fold(train_sp)
  post_k <- fit_k$post_list
  
  lm_cfh     <- fit_k$lm_cfh
  lm_cfh_len <- fit_k$lm_cfh_len
  lm_cfh_sig <- if (!is.null(lm_cfh)) summary(lm_cfh)$sigma else NA_real_
  lm_cfh_len_sig <- if (!is.null(lm_cfh_len)) summary(lm_cfh_len)$sigma else NA_real_
  
  for (spi in test_sp) {
    v <- wide_all %>% filter(Species == spi) %>% slice(1)
    if (nrow(v) == 0) next
    
    vals0 <- with(v, c(lnM=lnM, lnCF=lnCF, lnCH=lnCH, lnLF=lnLF, lnLH=lnLH))
    vals1 <- mask_one_circ(vals0, p_drop_circ)
    
    # ---------- (A) OLS: lnM ~ lnCFH ----------
    lnCF_ols <- vals1["lnCF"]; lnCH_ols <- vals1["lnCH"]
    if (!is.finite(lnCF_ols) && is.finite(vals1["lnCH"]) && !is.null(fit_k$imp_cf_from_ch))
      lnCF_ols <- as.numeric(predict(fit_k$imp_cf_from_ch, newdata = tibble(lnCH = vals1["lnCH"])))
    if (!is.finite(lnCH_ols) && is.finite(vals1["lnCF"]) && !is.null(fit_k$imp_ch_from_cf))
      lnCH_ols <- as.numeric(predict(fit_k$imp_ch_from_cf, newdata = tibble(lnCF = vals1["lnCF"])))
    lnCFH_ols <- if (is.finite(lnCF_ols) && is.finite(lnCH_ols)) logsumexp2(lnCF_ols, lnCH_ols) else NA_real_
    
    if (!is.null(lm_cfh) && is.finite(lnCFH_ols)) {
      pr <- predict(lm_cfh, newdata = tibble(lnCFH = lnCFH_ols), se.fit = TRUE)
      pred_sd <- sqrt(drop(pr$se.fit)^2 + lm_cfh_sig^2)
      out_rows[[length(out_rows)+1]] <- tibble(
        Species  = spi,
        true_lnM = vals1["lnM"],
        pred_lnM = as.numeric(pr$fit),
        pred_sd  = as.numeric(pred_sd),
        model    = "OLS: lnM ~ lnCFH",
        fold     = k
      )
    }
    
    # ---------- (B) OLS: lnM ~ lnCFH + lnLF + lnLH ----------
    if (!is.null(lm_cfh_len) && is.finite(lnCFH_ols) &&
        is.finite(vals1["lnLF"]) && is.finite(vals1["lnLH"])) {
      newd <- tibble(lnCFH = lnCFH_ols, lnLF = vals1["lnLF"], lnLH = vals1["lnLH"])
      pr2 <- try(predict(lm_cfh_len, newdata = newd, se.fit = TRUE), silent = TRUE)
      if (!inherits(pr2, "try-error")) {
        pred_sd2 <- sqrt(drop(pr2$se.fit)^2 + lm_cfh_len_sig^2)
        out_rows[[length(out_rows)+1]] <- tibble(
          Species  = spi,
          true_lnM = vals1["lnM"],
          pred_lnM = as.numeric(pr2$fit),
          pred_sd  = as.numeric(pred_sd2),
          model    = "OLS: lnM ~ lnCFH + lnLF + lnLH",
          fold     = k
        )
      }
    }
    
    # ---------- (B2) OLS: lnM ~ lnCFH + lnLF + lnLH (obs-only CFH) ----------
    # only if CFH is truly observed (no imputation) and both lengths present
    lnCFH_obs_only <- if (is.finite(vals1["lnCF"]) && is.finite(vals1["lnCH"]))
      logsumexp2(vals1["lnCF"], vals1["lnCH"]) else NA_real_
    if (!is.null(lm_cfh_len) && is.finite(lnCFH_obs_only) &&
        is.finite(vals1["lnLF"]) && is.finite(vals1["lnLH"])) {
      newd2 <- tibble(lnCFH = lnCFH_obs_only, lnLF = vals1["lnLF"], lnLH = vals1["lnLH"])
      pr2o <- try(predict(lm_cfh_len, newdata = newd2, se.fit = TRUE), silent = TRUE)
      if (!inherits(pr2o, "try-error")) {
        pred_sd2o <- sqrt(drop(pr2o$se.fit)^2 + lm_cfh_len_sig^2)
        out_rows[[length(out_rows)+1]] <- tibble(
          Species  = spi,
          true_lnM = vals1["lnM"],
          pred_lnM = as.numeric(pr2o$fit),
          pred_sd  = as.numeric(pred_sd2o),
          model    = "OLS: lnM ~ lnCFH + lnLF + lnLH (obs-only CFH)",
          fold     = k
        )
      }
    }
    
    # ---------- (C) Bayes: p(m | CFH) when CFH available (full NIG predictive) ----------
    by <- build_y(vals1, post_k, keep = c("lnCFH"))
    if (length(by$trait_ids) > 0) {
      res <- predict_mass_laplace(
        y_obs      = by$y,
        trait_ids  = by$trait_ids,
        prior_mean = fit_k$prior_mean,
        prior_sd   = fit_k$prior_sd*1000,
        post_dict  = post_k,
        search_sd  = fit_k$search_sd
      )
      out_rows[[length(out_rows)+1]] <- tibble(
        Species  = spi,
        true_lnM = vals1["lnM"],
        pred_lnM = res$mean,
        pred_sd  = sqrt(res$var),
        model    = "Bayes: CFH available (p(m|CFH))",
        fold     = k
      )
    }
    
    # ---------- (C2) Bayes: CFH+LF+LH (weak prior, plug-in, matched) ----------
    # require CFH actually observed + both lengths present to match OLS obs-only
    by_match <- build_y(vals1, post_k, keep = c("lnCFH","lnLF","lnLH"))
    if (all(c("lnCFH","lnLF","lnLH") %in% by_match$trait_ids)) {
      res_w <- predict_mass_closed_form(
        y_obs      = by_match$y,
        trait_ids  = by_match$trait_ids,
        prior_mean = fit_k$prior_mean,
        prior_sd   = fit_k$prior_sd,  # very weak
        post_dict  = post_k
      )
      out_rows[[length(out_rows)+1]] <- tibble(
        Species  = spi,
        true_lnM = vals1["lnM"],
        pred_lnM = res_w$mean,
        pred_sd  = sqrt(res_w$var),
        model    = "Bayes: CFH+LF+LH (weak prior, plug-in, matched)",
        fold     = k
      )
    }
    
    # ---------- (D) Bayes: CFH not constructible → via CF/CH/LF/LH ----------
    by3 <- build_y(vals1, post_k, keep = c("lnCF","lnCH","lnLF","lnLH"))
    if (length(by3$trait_ids) > 0 && !("lnCFH" %in% by3$trait_ids)) {
      res2 <- predict_mass_laplace(
        y_obs      = by3$y,
        trait_ids  = by3$trait_ids,
        prior_mean = fit_k$prior_mean,
        prior_sd   = fit_k$prior_sd,
        post_dict  = post_k,
        search_sd  = fit_k$search_sd
      )
      out_rows[[length(out_rows)+1]] <- tibble(
        Species  = spi,
        true_lnM = vals1["lnM"],
        pred_lnM = res2$mean,
        pred_sd  = sqrt(res2$var),
        model    = "Bayes: CFH missing (via CF/CH/LF/LH)",
        fold     = k
      )
    }
  }
}

cv_all <- bind_rows(out_rows)

# ============================================================
# 7) Metrics & visuals
# ============================================================

metrics <- cv_all %>%
  mutate(
    error     = pred_lnM - true_lnM,
    abs_error = abs(error),
    cover95   = (true_lnM >= pred_lnM - 1.96*pred_sd) &
      (true_lnM <= pred_lnM + 1.96*pred_sd)
  ) %>%
  group_by(model) %>%
  summarise(
    N      = n(),
    RMSE   = sqrt(mean(error^2, na.rm=TRUE)),
    MAE    = mean(abs_error, na.rm=TRUE),
    Cover95= mean(cover95, na.rm=TRUE),
    .groups="drop"
  ) %>%
  arrange(match(model, report_models), model)

print(metrics)

# Focused comparisons
p1 <- ggplot(cv_all %>% filter(model %in% c("OLS: lnM ~ lnCFH + lnLF + lnLH (obs-only CFH)",
                                            "Bayes: CFH+LF+LH (weak prior, plug-in, matched)")),
             aes(x=true_lnM, y=pred_lnM, color=model)) +
  geom_point(alpha=0.6) +
  geom_abline(slope=1, intercept=0, linetype="dashed") +
  labs(title = "Matched-sample: OLS vs Bayes (weak prior, plug-in) on CFH+LF+LH",
       x = "Observed ln(M)", y = "Predicted ln(M)", color = "Model") +
  theme_minimal()

p2 <- ggplot(cv_all %>% filter(model %in% c("OLS: lnM ~ lnCFH + lnLF + lnLH",
                                            "Bayes: CFH missing (via CF/CH/LF/LH)")),
             aes(x=true_lnM, y=pred_lnM, color=model)) +
  geom_point(alpha=0.6) +
  geom_abline(slope=1, intercept=0, linetype="dashed") +
  labs(title = "When CFH not constructible: OLS(imputed CFH) vs Bayes(no imputation)",
       x = "Observed ln(M)", y = "Predicted ln(M)", color = "Model") +
  theme_minimal()

print(p1); print(p2)
