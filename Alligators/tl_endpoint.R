#### =========================
#### EVT on log Total Length with point imputation from log SVL → log TL
#### Ordered & commented; LOO cross-validation at the end
#### =========================

## ---------------------------
## 0) Libraries & global setup
## ---------------------------
library(ggplot2)
library(extRemes)    # fevd/revd/pevd for GP fits & sims
library(dplyr)
library(readxl)
library(pracma)      # trapz for numerical integration
library(HDInterval)  # hdi: highest posterior density intervals
library(MASS)        # mvrnorm for allometry uncertainty (weight mapping)

## Reproducibility (used by bootstrap/MC draws)
set.seed(42)

## ---------------------------
## 1) Settings & plotting theme
## ---------------------------
ci_level      <- 0.95
n_boot_ks     <- 1000        # parametric bootstrap size for PIT→Uniform KS p-value
n_post_ystar  <- 20000       # posterior draws for TL endpoint y*
thresh_q_opt  <- 0.93        # threshold u0 = q_0.93(log TL)
fig_dir       <- "Figures"   # output directory for figures

# Reference lines (Stokes alligator) for posterior plots
stokes_tl_cm  <- 450   # TL = 450 cm
stokes_w_kg   <- 459   # weight = 459 kg

# Clean, consistent figure theme
theme_science_polished <- theme_minimal(base_family = "Arial", base_size = 12) +
  theme(
    axis.title       = element_text(size = 14, face = "bold"),
    axis.text        = element_text(size = 12),
    panel.grid.major = element_line(color = "gray80", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.6),
    axis.ticks.length= unit(0.20, "cm"),
    axis.ticks       = element_line(color = "black", linewidth = 0.4),
    plot.margin      = margin(5, 5, 5, 5),
    legend.position  = "right"
  )

## ---------------------------
## 2) Data ingest
## ---------------------------
data_path <- "Data/experimental_alligator_harvest_woodward.xlsx"
dat <- read_excel(data_path)

# Expected columns (informational)
svl_col    <- "SVL"     # snout–vent length (cm)
tl_col     <- "TL"      # total length (cm)
deform_col <- "Deform"  # 0=no deformation, 1=broken, 2=missing
weight_col <- "WTkg"    # weight (kg)

# Ensure bare names exist (no-op if columns already present)
dat <- dat %>%
  mutate(SVL, TL, WTkg, Deform)

## ---------------------------
## 3) Helper utilities
## ---------------------------

# Fit a GPD (GP) at threshold u on vector y; return estimates & covariance
fit_gpd_at_u <- function(y, u) {
  fit <- fevd(y, type = "GP", threshold = u)
  par <- fit$results$par
  list(
    fit   = fit,
    scale = unname(par["scale"]),
    shape = unname(par["shape"]),
    cov   = summary(fit)$cov.theta
  )
}

# Mean residual life values on a grid of thresholds
mrl_data <- function(y, u_seq) {
  sapply(u_seq, function(u){
    exc <- y[y > u] - u
    if (length(exc) < 5) return(NA_real_)
    mean(exc)
  })
}

# Build threshold grid from quantiles
make_thresholds <- function(y, q_from = 0.75, q_to = 0.95, n = 50){
  rng <- quantile(y, c(q_from, q_to), na.rm = TRUE)
  seq(rng[1], rng[2], length.out = n)
}

# Adjust scale to common anchor u0 for stability visualization
adj_scale <- function(scale, shape, u, u0) scale - shape * (u - u0)

# Diagnostics: Q–Q, P–P, and PIT-uniformity proxy for fitted GP at u
diagnostic_plots <- function(y, u, scale_hat, shape_hat, label) {
  exc <- y[y > u] - u; n <- length(exc)
  probs  <- ppoints(n)
  theo_q <- if (abs(shape_hat) > 1e-10) {
    u + scale_hat/shape_hat * (probs^(-shape_hat) - 1)
  } else {
    u - scale_hat * log(probs)
  }
  dfqq <- data.frame(Theoretical = rev(theo_q), Empirical = sort(y[y > u]))
  pqq <- ggplot(dfqq, aes(Theoretical, Empirical)) +
    geom_point(color = "steelblue") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    labs(title = paste0("Q–Q (", label, ")"),
         x = "Theoretical", y = "Empirical") +
    theme_science_polished
  
  F_theo <- if (abs(shape_hat) > 1e-10) {
    1 - (1 + shape_hat * exc / scale_hat)^(-1/shape_hat)
  } else {
    1 - exp(-exc / scale_hat)
  }
  dfpp <- data.frame(Theoretical = sort(F_theo), Empirical = (1:n)/n)
  ppp <- ggplot(dfpp, aes(Theoretical, Empirical)) +
    geom_point(color = "darkgreen") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    labs(title = paste0("P–P (", label, ")"),
         x = "Theoretical CDF", y = "Empirical CDF") +
    theme_science_polished
  
  F_hat <- pevd(exc, scale = scale_hat, shape = shape_hat, type="GP")
  pks <- ggplot(data.frame(F_hat = F_hat), aes(F_hat)) +
    geom_histogram(aes(y = ..density..), bins = 20,
                   fill = "skyblue", color = "black", alpha = 0.7) +
    geom_hline(yintercept = 1, color = "red", linetype = "dashed") +
    geom_density(color = "darkblue", linewidth = 1.1, adjust = 1.5) +
    labs(title = paste0("Uniformity: PIT (", label, ")"),
         x = expression(hat(F)(y)), y = "Density") +
    theme_science_polished
  
  list(pqq = pqq, ppp = ppp, pks = pks)
}

# KS bootstrap p-value for PIT ≈ Uniform(0,1)
ks_boot_pvalue <- function(y, u, scale_hat, shape_hat, B = 1000, seed = 42) {
  set.seed(seed)
  exc    <- y[y > u] - u
  n      <- length(exc)
  F_obs  <- pevd(exc, scale = scale_hat, shape = shape_hat, type="GP")
  ks_obs <- suppressWarnings(ks.test(F_obs, "punif")$statistic)
  
  ks_b <- rep(NA_real_, B)
  for (b in seq_len(B)) {
    # simulate exceedances under fitted GP
    sim_exc <- revd(n, scale = scale_hat, shape = shape_hat, type="GP")
    # splice with bulk (not used by fit here but keeps structure similar)
    sim_y <- c(y[y <= u], u + sim_exc)
    fb <- try(fevd(sim_y, type = "GP", threshold = u, time.units = "none"), silent = TRUE)
    if (inherits(fb, "try-error")) next
    par_b <- fb$results$par
    F_b <- pevd(sim_exc, scale = par_b["scale"], shape = par_b["shape"], type = "GP")
    ks_b[b] <- suppressWarnings(ks.test(F_b, "punif")$statistic)
  }
  ks_b <- ks_b[is.finite(ks_b)]
  list(ks_obs = as.numeric(ks_obs), p_value = mean(ks_b >= ks_obs), ks_boot = ks_b)
}

## ---------------------------------------------------------
## 4) Impute missing TL via SVL→TL regression (log–log OLS)
## ---------------------------------------------------------
reg_df <- dat %>%
  filter(Deform == 0L, is.finite(TL), is.finite(SVL)) %>%
  transmute(logTL = log(TL), logSVL = log(SVL))
stopifnot(nrow(reg_df) >= 10)  # sanity check

lm_tl_svl <- lm(logTL ~ logSVL, data = reg_df)

# Point prediction for missing/deformed TL on the log scale
dat <- dat %>%
  mutate(
    logSVL    = ifelse(is.finite(SVL), log(SVL), NA_real_),
    logTL_obs = ifelse(is.finite(TL),  log(TL),  NA_real_),
    logTL_pred= ifelse(is.finite(logSVL),
                       as.numeric(predict(lm_tl_svl, newdata = data.frame(logSVL = logSVL))),
                       NA_real_),
    # Use observed logTL when available; otherwise point-imputed prediction
    logTL_use = ifelse(is.finite(logTL_obs), logTL_obs, logTL_pred)
  )

# Tail variable for POT
y <- dat$logTL_use
y <- y[is.finite(y)]
stopifnot(length(y) >= 50)

## ---------------------------------------------------------
## 5) Threshold selection (MRL; shape & adjusted-scale stability)
## ---------------------------------------------------------
u_seq <- make_thresholds(y, q_from = 0.75, q_to = 0.95, n = 50)
u0    <- as.numeric(quantile(y, thresh_q_opt, na.rm = TRUE))  # chosen anchor

# MRL plot
mrl_vals <- mrl_data(y, u_seq)
p_mrl <- ggplot(data.frame(u = u_seq, mrl = mrl_vals), aes(u, mrl)) +
  geom_line() + geom_point() +
  geom_vline(xintercept = u0, linetype = "dashed", color = "red") +
  labs(title = "MRL plot (log TL with point imputation)",
       x = "Threshold (log TL)", y = "Mean excess") +
  theme_science_polished
p_mrl
ggsave(file.path(fig_dir, "tl_log_mrl.png"), p_mrl, dpi = 600, width = 7, height = 5, units = "in")

# Stability scans
shape <- scale <- shape_lo <- shape_hi <- adj <- adj_lo <- adj_hi <- rep(NA_real_, length(u_seq))
for (i in seq_along(u_seq)) {
  u <- u_seq[i]
  exc <- y[y > u] - u
  if (length(exc) < 20) next
  out <- try(fit_gpd_at_u(y, u), silent = TRUE)
  if (inherits(out, "try-error")) next
  
  scale[i] <- out$scale; shape[i] <- out$shape
  zcrit     <- qnorm(1 - (1 - ci_level)/2)
  se_scale  <- sqrt(out$cov[1,1]); se_shape <- sqrt(out$cov[2,2])
  shape_lo[i] <- shape[i] - zcrit * se_shape
  shape_hi[i] <- shape[i] + zcrit * se_shape
  
  adj[i] <- adj_scale(scale[i], shape[i], u, u0)
  var_adj <- out$cov[1,1] + (u - u0)^2 * out$cov[2,2] - 2*(u - u0)*out$cov[1,2]
  se_adj  <- sqrt(max(var_adj, 0))
  adj_lo[i] <- adj[i] - zcrit * se_adj
  adj_hi[i] <- adj[i] + zcrit * se_adj
}

p_shape <- ggplot(data.frame(u = u_seq, shape, shape_lo, shape_hi), aes(u, shape)) +
  geom_point() +
  geom_errorbar(aes(ymin = shape_lo, ymax = shape_hi), width = 0.03, color = "blue") +
  geom_vline(xintercept = u0, color = "red", linetype = "dashed") +
  labs(title = "Shape stability (log TL)", x = "Threshold (log TL)", y = "Shape (xi)") +
  theme_science_polished
p_shape
ggsave(file.path(fig_dir, "tl_log_shape_stability.png"), p_shape, dpi = 600, width = 7, height = 5, units = "in")

p_adj <- ggplot(data.frame(u = u_seq, adj, adj_lo, adj_hi), aes(u, adj)) +
  geom_point() +
  geom_errorbar(aes(ymin = adj_lo, ymax = adj_hi), width = 0.03, color = "blue") +
  geom_vline(xintercept = u0, color = "red", linetype = "dashed") +
  labs(title = "Adjusted scale stability (log TL)",
       x = "Threshold (log TL)", y = expression(tilde(sigma)(u[0]))) +
  theme_science_polished
p_adj
ggsave(file.path(fig_dir, "tl_log_adj_scale_stability.png"), p_adj, dpi = 600, width = 7, height = 5, units = "in")

## ---------------------------------------------------------
## 6) Final GP fit @ u0 + tail diagnostics
## ---------------------------------------------------------
fit0  <- fit_gpd_at_u(y, u0)
diags <- diagnostic_plots(y, u0, fit0$scale, fit0$shape, "log TL")
ksout <- ks_boot_pvalue(y, u0, fit0$scale, fit0$shape, B = n_boot_ks)

diags$pqq; ggsave(file.path(fig_dir, "tl_log_qq.png"),  diags$pqq, dpi = 600, width = 7, height = 5, units = "in")
diags$ppp; ggsave(file.path(fig_dir, "tl_log_pp.png"),  diags$ppp, dpi = 600, width = 7, height = 5, units = "in")
diags$pks; ggsave(file.path(fig_dir, "tl_log_pit_uniformity.png"), diags$pks, dpi = 600, width = 7, height = 5, units = "in")

message(sprintf("[log TL] u0=%.3f | scale=%.4f shape=%.4f | KS p=%.3f",
                u0, fit0$scale, fit0$shape, ksout$p_value))

## ---------------------------------------------------------
## 7) Endpoint posterior for y* (Weibull tail: xi<0), xi marginalized
## ---------------------------------------------------------
y_above <- y[y > u0]; n_ex <- length(y_above)
if (n_ex < 10) stop("Not enough exceedances for endpoint posterior.")

Ymax <- max(y_above)
eps  <- max(1e-6, 1e-4 * abs(Ymax))
y_star_grid <- seq(Ymax + eps, Ymax + 7, length.out = 2000)   # grid for y*
xi_grid     <- seq(-1.0, -0.02, length.out = 1200)            # xi < 0 (Weibull)

# Log-likelihood ℓ(y*,xi) for exceedances above u0
loglik_endpoint <- function(y_star, xi, u, y_ex) {
  if (xi >= 0 || any(y_star <= y_ex)) return(-Inf)
  n <- length(y_ex)
  term1 <- n * log(y_star - u) / xi
  term2 <- -(1/xi + 1) * sum(log(y_star - y_ex))
  term3 <- -n * log(-xi)
  term1 + term2 + term3
}

# Marginalize xi numerically → log p(y* | data) up to a constant
marg_log_post <- sapply(y_star_grid, function(ys) {
  lv <- sapply(xi_grid, function(xi) loglik_endpoint(ys, xi, u0, y_above))
  lv_max <- max(lv)
  log(trapz(xi_grid, exp(lv - lv_max))) + lv_max
})

post_unnorm <- exp(marg_log_post - max(marg_log_post))
post        <- post_unnorm / trapz(y_star_grid, post_unnorm)

# Posterior draws for y* on log scale, then transform to TL (cm)
ystar_samples <- sample(y_star_grid, size = n_post_ystar, replace = TRUE, prob = post)
ystar_map     <- y_star_grid[which.max(post)]
ystar_hpdi    <- hdi(ystar_samples, credMass = ci_level)

# TL* (cm) summaries
tl_star_cm <- exp(ystar_samples)
tl_map_cm  <- exp(ystar_map)
tl_hpdi_cm <- exp(ystar_hpdi)

# Plot TL* posterior (+ Stokes dashed)
p_tlstar <- ggplot(data.frame(TLstar = tl_star_cm), aes(TLstar)) +
  geom_histogram(aes(y = ..density..), bins = 60,
                 fill = "skyblue", color = "black", alpha = 0.6) +
  geom_density(color = "darkblue", linewidth = 1.2) +
  geom_vline(xintercept = tl_map_cm,  color = "purple",   linetype = "dashed",  linewidth = 1.1) +
  geom_vline(xintercept = tl_hpdi_cm, color = "orange",   linetype = "dotdash", linewidth = 1.0) +
  geom_vline(xintercept = stokes_tl_cm, color = "firebrick", linetype = "dashed", linewidth = 1.0) +
  annotate("text", x = stokes_tl_cm, y = 0.02, vjust = -0.8,
           label = "Stokes (450 cm)", color = "firebrick", angle = 90, size = 4) +
  labs(
       x = "TL* (cm)", y = "Density") +
  xlim(400, 500) +
  theme_science_polished
p_tlstar
ggsave(file.path(fig_dir, "tl_endpoint_posterior.png"),
       p_tlstar, dpi = 600, width = 7, height = 5, units = "in")

cat("\n--- TL endpoint (cm) ---\n",
    "MAP:", round(tl_map_cm, 2), "\n",
    "HPDI:", round(tl_hpdi_cm[1], 2), "-", round(tl_hpdi_cm[2], 2), "\n")

## ---------------------------------------------------------
## 8) Map TL* to Weight endpoint (regression + parameter/residual uncertainty)
## ---------------------------------------------------------
# Fit log-weight ~ log-TL on undeformed specimens with both observed
w_fit_df <- dat %>%
  filter(Deform == 0L, is.finite(WTkg), is.finite(TL)) %>%
  transmute(logW = log(WTkg), logTL = log(TL))

# Fallback if undeformed subset is tiny: allow all finite (keeps structure)
if (nrow(w_fit_df) < 10) {
  w_fit_df <- dat %>%
    filter(is.finite(WTkg), is.finite(TL)) %>%
    transmute(logW = log(WTkg), logTL = log(TL))
}
stopifnot(nrow(w_fit_df) >= 10)

lm_w_tl <- lm(logW ~ logTL, data = w_fit_df)
Vcoef   <- vcov(lm_w_tl)
beta    <- coef(lm_w_tl)                   # (intercept a, slope b)
sigma_e <- summary(lm_w_tl)$sigma          # residual SD on logW

# Individual upper quantile at TL*: W*_q = exp(a + b * y* + z_q * sigma_e)
q_indiv <- 0.95
zq      <- qnorm(q_indiv)

# Draw (a,b) to propagate regression parameter uncertainty
ab_samps <- MASS::mvrnorm(n_post_ystar, mu = beta, Sigma = Vcoef)
a_s <- ab_samps[,1]; b_s <- ab_samps[,2]

logW_star_struct <- a_s + b_s * ystar_samples     # structural mapping at y*
W_star_struct    <- exp(logW_star_struct)

logW_star_q <- logW_star_struct + zq * sigma_e    # add residual quantile
W_star_q    <- exp(logW_star_q)

# Summaries for W*_q
sum_struct <- c(Median = median(W_star_struct),
                HPDI_low = hdi(W_star_struct, credMass = ci_level)[1],
                HPDI_high= hdi(W_star_struct, credMass = ci_level)[2])

sum_q <- c(Median = median(W_star_q),
           HPDI_low = hdi(W_star_q, credMass = ci_level)[1],
           HPDI_high= hdi(W_star_q, credMass = ci_level)[2])

cat("\nTypical individual upper (", q_indiv*100, "%):\n", sep = "")
print(round(sum_q, 2))

# Plot weight endpoint posterior (+ HPDI, median, Stokes dashed)
hpdi_w_q <- hdi(W_star_q, credMass = ci_level)
med_w_q  <- median(W_star_q)

p_w_q <- ggplot(data.frame(W = W_star_q), aes(W)) +
  geom_histogram(aes(y = ..density..), bins = 60,
                 fill = "khaki", color = "black", alpha = 0.6) +
  geom_density(linewidth = 1.2) +
  geom_vline(xintercept = med_w_q,   color = "purple",   linetype = "dashed",  linewidth = 1.1) +  # posterior median
  geom_vline(xintercept = hpdi_w_q,  color = "orange",   linetype = "dotdash", linewidth = 1.0) +  # HPDI bounds
  geom_vline(xintercept = stokes_w_kg, color = "firebrick", linetype = "dashed", linewidth = 1.0) +
  annotate("text", x = stokes_w_kg, y = 0.012, vjust = -0.8,
           label = "Stokes (459 kg)", color = "firebrick", angle = 90, size = 4) +
  labs(title = paste0("Weight endpoint W*_", q_indiv*100, "% (from TL*)"),
       x = paste0("W*_", q_indiv*100, "% (kg)"), y = "Density") +
  xlim(300, 700) +
  theme_science_polished
p_w_q
ggsave(file.path(fig_dir, "weight_endpoint_typical_upper.png"),
       p_w_q, dpi = 600, width = 7, height = 5, units = "in")

# Console recap
cat(sprintf("\n[log TL POT @ u0=%.3f] scale=%.4f, shape=%.4f; KS p=%.3f\n",
            u0, fit0$scale, fit0$shape, ksout$p_value))


