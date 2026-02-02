#### =========================
#### EVT on log Total Length with point imputation from log SVL → log TL
#### Ordered & commented; Woodward prior + LOO cross-validation (fast options)
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

## Ensure figure directory exists
fig_dir <- "Figures"
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

## ---------------------------
## 1) Settings & plotting theme
## ---------------------------
ci_level      <- 0.95
n_boot_ks     <- 1000        # parametric bootstrap size for PIT→Uniform KS p-value
n_post_ystar  <- 20000       # posterior draws for TL endpoint y*
thresh_q_opt  <- 0.93        # threshold u0 = q_0.93(log TL)

# Reference lines (Stokes alligator) for posterior plots
stokes_tl_cm  <- 450   # TL = 450 cm
stokes_w_kg   <- 459   # weight = 459 kg

# Clean, consistent figure theme
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
dat <- dat %>% mutate(SVL, TL, WTkg, Deform)

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
    labs(
         x = "Theoretical quantiles", y = "Empirical quantiles") +
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
    labs(
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
ks_boot_pvalue <- function(y, u, scale_hat, shape_hat, B = 1000, seed = 42, refit = FALSE) {
  set.seed(seed)
  exc    <- y[y > u] - u
  n      <- length(exc)
  F_obs  <- pevd(exc, scale = scale_hat, shape = shape_hat, type = "GP")
  ks_obs <- suppressWarnings(ks.test(F_obs, "punif")$statistic)
  
  ks_b <- rep(NA_real_, B)
  for (b in seq_len(B)) {
    # simulate exceedances under the fitted GP
    sim_exc <- revd(n, scale = scale_hat, shape = shape_hat, type = "GP")
    
    if (refit) {
      # (slower, sometimes unstable) — refit fevd to each bootstrap sample
      fb <- try(fevd(c(y[y <= u], u + sim_exc), type = "GP", threshold = u, time.units = "none"),
                silent = TRUE)
      if (inherits(fb, "try-error")) next
      par_b <- fb$results$par
      F_b   <- pevd(sim_exc, scale = par_b["scale"], shape = par_b["shape"], type = "GP")
    } else {
      # (fast, stable) — evaluate against the original fitted GP
      F_b <- pevd(sim_exc, scale = scale_hat, shape = shape_hat, type = "GP")
    }
    
    ks_b[b] <- suppressWarnings(ks.test(F_b, "punif")$statistic)
  }
  
  ks_b <- ks_b[is.finite(ks_b)]
  list(
    ks_obs  = as.numeric(ks_obs),
    p_value = if (length(ks_b)) mean(ks_b >= ks_obs) else NA_real_,
    ks_boot = ks_b
  )
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

u_seq <- make_thresholds(y, q_from = 0.72, q_to = 0.99, n = 50)
u0    <- as.numeric(quantile(y, thresh_q_opt, na.rm = TRUE))  # chosen anchor
fit0  <- fit_gpd_at_u(y, u0)
shape_hat0 <- fit0$shape
adj_hat0   <- fit0$scale  # adj_scale at u0 equals scale at u0

fit0 <- fit_gpd_at_u(y, u0)
shape_hat0 <- fit0$shape
adj_hat0   <- fit0$scale

shape_se0   <- sqrt(fit0$cov[2,2])
z       <- shape_hat0 / shape_se0
p_wald  <- pnorm(z)  # one-sided for H1: xi<0

p_wald
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
  labs(x = "Threshold (log TL)", y = "Shape parameter") +
  theme_science_polished + 
  geom_hline(yintercept = shape_hat0, color = "red", linetype = "dashed")
p_shape
ggsave(file.path(fig_dir, "tl_log_shape_stability.png"), p_shape, dpi = 600, width = 7, height = 5, units = "in")

p_adj <- ggplot(data.frame(u = u_seq, adj, adj_lo, adj_hi), aes(u, adj)) +
  geom_point() +
  geom_errorbar(aes(ymin = adj_lo, ymax = adj_hi), width = 0.03, color = "blue") +
  geom_vline(xintercept = u0, color = "red", linetype = "dashed") +
  labs(x = "Threshold (log TL)", y = "Adjusted scale parameter") +
  theme_science_polished + 
  geom_hline(yintercept = adj_hat0, color = "red", linetype = "dashed")
p_adj
ggsave(file.path(fig_dir, "tl_log_adj_scale_stability.png"), p_adj, dpi = 600, width = 7, height = 5, units = "in")

## ---------------------------------------------------------
## 6) Tail diagnostics
## ---------------------------------------------------------

diags <- diagnostic_plots(y, u0, fit0$scale, fit0$shape, "log TL")
ksout <- ks_boot_pvalue(y, u0, fit0$scale, fit0$shape, B = n_boot_ks)

diags$pqq; ggsave(file.path(fig_dir, "tl_log_qq.png"),  diags$pqq, dpi = 600, width = 7, height = 5, units = "in")
diags$ppp; ggsave(file.path(fig_dir, "tl_log_pp.png"),  diags$ppp, dpi = 600, width = 7, height = 5, units = "in")
diags$pks; ggsave(file.path(fig_dir, "tl_log_pit_uniformity.png"), diags$pks, dpi = 600, width = 7, height = 5, units = "in")

message(sprintf("[log TL] u0=%.3f | scale=%.4f shape=%.4f | KS p=%.3f",
                u0, fit0$scale, fit0$shape, ksout$p_value))

## ---------------------------------------------------------
## 7) Endpoint posterior with Woodward prior (grid-only; upper bound only)
##     - No sampling anywhere; pure trapezoidal integration on (y*, xi) grids
##     - Outputs: post (density over y*), y_star_grid, xi_grid, ystar_map, ystar_upper
## ---------------------------------------------------------
y_above <- y[y > u0]; n_ex <- length(y_above)
if (n_ex < 10) stop("Not enough exceedances for endpoint posterior.")

Ymax <- max(y_above)
eps  <- max(1e-6, 1e-4 * abs(Ymax))
y_star_grid <- seq(Ymax + eps, Ymax + 10, length.out = 2000)   # grid for y* (log TL*)
xi_grid     <- seq(-1.0, -0.02, length.out = 1000)             # xi < 0 (Weibull)

# Log-likelihood ℓ(y*, xi) for exceedances above u0 in Weibull domain
loglik_endpoint <- function(y_star, xi, u, y_ex) {
  if (xi >= 0 || any(y_star <= y_ex)) return(-Inf)
  n <- length(y_ex)
  term1 <- n * log(y_star - u) / xi
  term2 <- -(1/xi + 1) * sum(log(y_star - y_ex))
  term3 <- -n * log(-xi)
  term1 + term2 + term3
}

# ---- Priors (data-independent) ----
# Woodward prior for TL* centered at 430 cm on the log scale
mu_y <- log(430)
prior_mode <- "woodward_loose"  # {"woodward_loose", "woodward_mod", "woodward_tight"}
sd_y <- switch(prior_mode,
               "woodward_loose" = 0.095,
               "woodward_mod"   = 0.048,
               "woodward_tight" = 0.024, 0.048
)
logprior_y  <- function(y_star) dnorm(y_star, mean = mu_y, sd = sd_y, log = TRUE)

# PC prior on kappa = -xi in (0,1): Exp(lambda_k), truncated (0,1)
lambda_k   <- 3.0
logprior_xi <- function(xi) {
  if (xi >= 0 || xi <= -1) return(-Inf)
  kappa <- -xi
  if (kappa <= 0 || kappa >= 1) return(-Inf)
  log(lambda_k) - lambda_k * kappa
}

# ---- Numeric marginalization over xi for each y* (trapz over xi) ----
marg_log_post <- sapply(y_star_grid, function(ys) {
  lv <- sapply(xi_grid, function(xi) loglik_endpoint(ys, xi, u0, y_above) + logprior_xi(xi))
  lv_max <- max(lv)
  log_int_xi <- log(pracma::trapz(xi_grid, exp(lv - lv_max))) + lv_max
  log_int_xi + logprior_y(ys)
})

# Normalize over y* (trapz over y*)
post_unnorm <- exp(marg_log_post - max(marg_log_post))
Z_y         <- pracma::trapz(y_star_grid, post_unnorm)
if (!is.finite(Z_y) || Z_y <= 0) stop("Normalization failed for p(y*|y).")
post <- post_unnorm / Z_y  # posterior density on the y* grid (log scale)

# ---- Robust grid helpers using pracma::cumtrapz ----
cdf_from_density <- function(x, f) {
  f[!is.finite(f)] <- 0
  F <- pracma::cumtrapz(x, f)          # same length as x; starts at 0
  rng <- max(F, na.rm = TRUE)
  if (!is.finite(rng) || rng <= 0) stop("cdf_from_density: non-positive integral.")
  F / rng
}
# One-sided (upper) quantile from the CDF built on the grid (e.g., p = 0.95)
q_from_cdf <- function(x, f, p) {
  F <- cdf_from_density(x, f)
  as.numeric(approx(x = F, y = x, xout = p, ties = "ordered", rule = 2)$y)
}

# ---- Posterior summaries (log scale) ----
ystar_map   <- y_star_grid[which.max(post)]
ystar_upper <- q_from_cdf(y_star_grid, post, ci_level)   # e.g., 0.95 upper bound

# ---- Transform to TL (cm) for reporting/plotting ----
tl_grid_cm   <- exp(y_star_grid)
dens_tl      <- post / tl_grid_cm      # change-of-variables: f_TL(t) = f_Y(log t) * 1/t
tl_map_cm    <- exp(ystar_map)
tl_upper_cm  <- exp(ystar_upper)

# Optional visualization on TL scale
df_tlpost <- data.frame(TL = tl_grid_cm, dens = dens_tl)
p_tlstar <- ggplot(df_tlpost, aes(TL, dens)) +
  geom_line(linewidth = 1.2, color = "darkblue") +
  geom_vline(xintercept = tl_map_cm,    color = "purple",   linetype = "dashed",  linewidth = 1.1) +
  geom_vline(xintercept = tl_upper_cm,  color = "orange",   linetype = "dotdash", linewidth = 1.0) +
  geom_vline(xintercept = stokes_tl_cm, color = "firebrick", linetype = "dashed", linewidth = 1.0) +
  annotate("text", x = stokes_tl_cm, y = max(dens_tl)*0.08, vjust = -0.8,
           label = "Stokes (450 cm)", color = "firebrick", angle = 90, size = 4) +
  labs(x = "TL* (cm)", y = "Density") +
  xlim(416, 500) +
  theme_science_polished

print(p_tlstar)
ggsave(file.path(fig_dir, "tl_endpoint_posterior.png"),
       p_tlstar, dpi = 600, width = 7, height = 5, units = "in")

cat("\n--- TL endpoint (cm) via numeric integration (no sampling) ---\n",
    "MAP:          ", round(tl_map_cm, 2), "\n",
    sprintf("%d%% upper:  ", round(ci_level*100)), round(tl_upper_cm, 2), "\n")


## ---------------------------------------------------------
## 8) Map TL* to Weight endpoint — grid-only propagation (no sampling)
##     - Build p(m | y, D) = sum_i W_i * N(m ; μ_i, s_i^2) on log scale
##     - Report one-sided upper bounds (no HPDI, no draws)
## ---------------------------------------------------------
## ---------------------------------------------------------
## 8) Mass extrapolation (sauropod-style): sample (a,b), sample y*, add residual quantile
## ---------------------------------------------------------


ystar_draws  <- sample(y_star_grid, size = n_post_ystar, replace = TRUE, prob = post)
n_post_ystar <- 1000000

# 8.1 Fit (or reuse) log-weight ~ log-TL on undeformed specimens
w_fit_df <- dat %>%
  filter(Deform == 0L, is.finite(WTkg), is.finite(TL)) %>%
  transmute(logW = log(WTkg), logTL = log(TL))

if (nrow(w_fit_df) < 10) {
  w_fit_df <- dat %>% filter(is.finite(WTkg), is.finite(TL)) %>%
    transmute(logW = log(WTkg), logTL = log(TL))
}
stopifnot(nrow(w_fit_df) >= 10)

lm_w_tl <- lm(logW ~ logTL, data = w_fit_df)
beta    <- coef(lm_w_tl)           # (a, b)
Vcoef   <- vcov(lm_w_tl)
sigma_e <- summary(lm_w_tl)$sigma  # residual SD on logW

# 8.2 Sample regression parameters exactly like the sauropod code
ab_samps <- MASS::mvrnorm(n = length(ystar_draws), mu = beta, Sigma = Vcoef)
a_s <- ab_samps[, 1]
b_s <- ab_samps[, 2]

# 8.3 Structural mapping to log mass at each y* draw
logM_draws <- a_s + b_s * ystar_draws

# 8.4 Add residual uncertainty for a "typical individual" upper endpoint (q_indiv)
q_indiv <- 0.95
logM_draws <- logM_draws + qnorm(q_indiv) * sigma_e

# 8.5 Transform to kilograms (or tons) and read off summaries
mass_draws_kg <- exp(logM_draws)
mass_up95_kg  <- as.numeric(quantile(mass_draws_kg, probs = ci_level))

# MAP (mode) via kernel density on kg
dens_m <- density(mass_draws_kg, n = 4096, adjust = 1.0)
mass_map_kg <- dens_m$x[which.max(dens_m$y)]

# 8.6 Plot (smooth density with MAP + upper bound)
p_mass <- ggplot(data.frame(W = mass_draws_kg), aes(W)) +
  geom_density(linewidth = 1.2, color = "darkblue") +
  geom_vline(xintercept = mass_map_kg,  color = "purple",  linetype = "dashed",  linewidth = 1) +
  geom_vline(xintercept = mass_up95_kg, color = "orange",  linetype = "dotdash", linewidth = 1) +
  geom_vline(xintercept = stokes_w_kg,  color = "firebrick", linetype = "dashed", linewidth = 1) +
  annotate("text", x = stokes_w_kg, y = max(dens_m$y)*0.08, vjust = -0.8,
           label = "Stokes (459 kg)", color = "firebrick", angle = 90, size = 4) +
  labs(x = "Weight (kg)", y = "Density") +
  theme_science_polished +
  xlim(350, 600)

print(p_mass)
ggsave(file.path(fig_dir, "mass_endpoint.png"),
       p_mass, dpi = 600, width = 7, height = 5, units = "in")

cat(sprintf("\nMass endpoint (typical individual, q=%d%%):\n", round(q_indiv*100)))
cat(sprintf("MAP: %.2f kg | %d%% upper: %.2f kg\n",
            mass_map_kg, round(ci_level*100), mass_up95_kg))

## ---------------------------------------------------------
## 8b) Regression panel: log weight (y) vs log TL (x)
## ---------------------------------------------------------
# Data for scatter (only complete, undeformed weights & TL)
# --- Clean plotting data: drop NA sex to avoid "NA" legend entry
plot_df <- dat %>%
  filter(is.finite(WTkg), is.finite(TL), !is.na(Sex)) %>%   # <— added !is.na(Sex)
  transmute(logTL = log(TL),
            logW  = log(WTkg),
            Sex   = factor(Sex, levels = c("F","M")))

sex_cols <- c("F" = "#ff6fb3", "M" = "#66b3ff")

# Threshold line y-value
y_at_u0 <- as.numeric(predict(lm_w_tl, newdata = data.frame(logTL = u0)))

# Endpoint x (from TL*) and consistent y on the regression line
x_ep    <- ystar_map                # MAP log TL*
x_ep_ub <- ystar_upper              # upper bound for log TL*
y_ep    <- as.numeric(predict(lm_w_tl, data.frame(logTL = x_ep)))    # ON the line
zq      <- qnorm(0.95)
y_ep_ub <- y_ep + zq * sigma_e      # 95% individual at same x

# Stokes
x_stokes <- log(stokes_tl_cm)
y_stokes <- log(stokes_w_kg)

ab <- coef(lm_w_tl)

p_reg <- ggplot(plot_df, aes(x = logTL, y = logW, color = Sex)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_abline(intercept = ab[1], slope = ab[2], linewidth = 1) +
  geom_vline(xintercept = u0, linetype = "dashed", color = "gray30") +
  geom_hline(yintercept = y_at_u0, linetype = "dashed", color = "gray30") +
  # Endpoint dot on the regression line + L-shaped upper bounds
  geom_point(aes(x = x_ep, y = y_ep), color = "darkgreen", size = 3, inherit.aes = FALSE) +
  geom_segment(aes(x = x_ep, y = y_ep, xend = x_ep_ub, yend = y_ep),
               color = "darkgreen", linewidth = 1.1, inherit.aes = FALSE) +
  geom_segment(aes(x = x_ep, y = y_ep, xend = x_ep, yend = y_ep_ub),
               color = "darkgreen", linewidth = 1.1, inherit.aes = FALSE) +
  # Stokes
  geom_point(aes(x = x_stokes, y = y_stokes),
             inherit.aes = FALSE, color = "firebrick", size = 2.8) +
  geom_text(aes(x = x_stokes, y = y_stokes, label = "Stokes"),
            inherit.aes = FALSE, nudge_y = 0.015, color = "firebrick", size = 3.5) +
  scale_color_manual(values = sex_cols,
                     breaks = c("F","M"),           # only show F/M
                     na.translate = FALSE, drop = FALSE) +  # remove NA from legend
  labs(x = "log Total Length", y = "log Weight", color = "Sex") +
  theme_science_polished


# Save
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)
ggsave(file.path(fig_dir, "reg_logW_vs_logTL_with_threshold_endpoint_stokes.png"),
       p_reg, dpi = 600, width = 8, height = 5.2, units = "in")

## ---------------------------------------------------------
## 9) Assessing the exceedance probability of the Stokes alligator (TL = 450 cm)
##     P(Y > y_claim | data) = P(Y > u0 | data) * E_{y*,xi|data}[( (y*-y)/(y*-u0) )^{-1/xi}]
##     Reuses the existing y* and xi grids and trapezoidal weights
## ---------------------------------------------------------

# ---- Helper: trapezoid weights on a grid for a (possibly unnormalized) density vector f(x) ----
trapz_weights <- function(x, f_un) {
  # trapezoid areas for each interval, then convert to point-weights
  # We’ll build "point weights" that sum to 1 and integrate f(x) over x.
  stopifnot(length(x) == length(f_un), length(x) >= 2)
  # interval contributions
  dx <- diff(x)
  mids <- (f_un[-1] + f_un[-length(f_un)]) * dx / 2
  total <- sum(mids)
  if (!is.finite(total) || total <= 0) {
    return(rep(NA_real_, length(x)))
  }
  # Convert interval masses to point masses:
  # interior points get half from left + half from right; endpoints get half of their single adjacent interval
  w <- numeric(length(x))
  w[1]            <- mids[1] / 2
  w[length(x)]    <- mids[length(mids)] / 2
  if (length(x) > 2) {
    w[2:(length(x)-1)] <- (mids[-length(mids)] + mids[-1]) / 2
  }
  w / sum(w)
}

# ---- Objects reused from earlier sections (already defined above this block) ----
# y                : vector of log TL used for POT
# u0               : chosen threshold on log TL
# y_above          : exceedances above u0 on log scale (defined in §7)
# y_star_grid      : grid for y* (log TL*), length I
# xi_grid          : grid for xi (< 0), length J
# loglik_endpoint  : function(y_star, xi, u, y_ex)
# logprior_xi      : PC prior log-density on xi
# post             : normalized posterior density over y* on the grid (vector length I)

stopifnot(exists("y"), exists("u0"), exists("y_above"),
          exists("y_star_grid"), exists("xi_grid"),
          exists("loglik_endpoint"), exists("logprior_xi"),
          exists("post"))

# Sanity: build trapezoidal weights W_i over y* from the (normalized) posterior "post"
W_y <- trapz_weights(y_star_grid, post)
stopifnot(all(is.finite(W_y)), abs(sum(W_y) - 1) < 1e-6)

# Exceedance rate P(Y > u0 | data) estimated empirically
p_exceed_u <- mean(y > u0, na.rm = TRUE)

# Claim: Stokes TL = 450 cm
y_claim_cm  <- stokes_tl_cm    # 450
y_claim_log <- log(y_claim_cm)

if (y_claim_log <= u0) {
  warning("Claim TL is not above the GP threshold; using full formula but note the tail-factor is conditional on Y>u0.")
}

# For each y*_i, compute conditional posterior over xi (via trapezoid-normalization),
# then the conditional survival S_i(y) = E_{xi | y*, data}[ ((y* - y)/(y* - u0))^{-1/xi} ] if y* > y; else 0.
Si_vec <- numeric(length(y_star_grid))

for (i in seq_along(y_star_grid)) {
  ys <- y_star_grid[i]
  if (ys <= y_claim_log) {
    Si_vec[i] <- 0
    next
  }
  # unnormalized conditional over xi at fixed y*:
  lv <- sapply(xi_grid, function(xi) {
    loglik_endpoint(ys, xi, u0, y_above) + logprior_xi(xi)
  })
  # stabilize and exponentiate
  m  <- max(lv)
  fu <- exp(lv - m)           # proportional to p(xi | y*, data)
  V_xi <- trapz_weights(xi_grid, fu)  # conditional weights over xi, sum to 1
  
  # survival term at (ys, xi): ((ys - y)/(ys - u0))^{-1/xi}
  surv_term <- ((ys - y_claim_log) / (ys - u0))^(-1/xi_grid)
  # guard tiny numerical issues (should all be finite for xi<0 and ys>y_claim>u0)
  surv_term[!is.finite(surv_term)] <- 0
  
  Si_vec[i] <- sum(V_xi * surv_term)
}

# Posterior-averaged conditional tail probability P(Y > y | Y > u0, data)
p_tail_cond <- sum(W_y * Si_vec)

# Final estimate: P(Y > y | data)
p_exceed_claim <- p_exceed_u * p_tail_cond

cat(sprintf("\n--- Exceedance probability for Stokes alligator ---\n"))
cat(sprintf("Threshold u0 (log cm)  : %.4f  (≈ %.1f cm)\n", u0, exp(u0)))
cat(sprintf("Claim y (log cm)       : %.4f  (450.0 cm)\n", y_claim_log))
cat(sprintf("P(Y > u0 | data)       : %.4f\n", p_exceed_u)) # 93% quantile.
cat(sprintf("P(Y > y | Y>u0, data)  : %.4f\n", p_tail_cond))

cat(sprintf("P(Y > y | data)        : %.6f\n", p_exceed_claim))

#######################
# CONFIDENCE INTERVALS
#######################
