############################################################
#  Title: POT analysis on log CFH (Normal prior centered at u0+1.05 ≈ 600 t)
#  Author: Bastiaan Aelbrecht (rewrite with prior overlays)
#  Date  : 2025-09-09
############################################################
library(ggplot2)
library(extRemes)      # fevd/revd/pevd
library(dplyr)
library(readxl)
library(pracma)        # trapz
library(HDInterval)    # hdi
library(MASS)          # mvrnorm (for allometry uncertainty if needed)
library(scales)
library(grid)

set.seed(42)

FIG_DIR <- "Figures"
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

## ---------------------------
## 1) Settings
## ---------------------------
ci_level        <- 0.95
n_post_ystar    <- 20000
thresh_q_opt    <- 0.78                    # anchor quantile for u0
xi_prior_mode   <- "improper"              # "improper" (Weibull domain) | "pc"
lambda_k        <- 3.0                     # only used for "pc"

# Prior on y* (log CFH): Normal centered at u0 + delta  (≈600 t anchor)
prior_center_delta <- 0.6                 # center mu_y = u0 + 1.05
prior_sd_y         <- 0.1                 # width on log-CFH scale (tune: 0.15–0.25)

# Plot option: overlay prior densities
plot_prior <- TRUE

## ---------------------------
## 2) Data ingest
## ---------------------------
df <- read_excel("Data/DEMic23_updated_Supplemental_Data_withPubYear.xlsx")
df <- df[, c("hum+fem circ (mm)", "publication year")]
colnames(df) <- c("sum_circ_mm", "year")
df <- na.omit(df)
df$log_circ <- log(as.double(df$sum_circ_mm))
y <- df$log_circ
stopifnot(length(y) >= 80)

## ---------------------------
## 3) Allometry (for mass propagation)
## ---------------------------
coeff_file_cq <- "centered_quadratic_coefficients.rds"
coeffs <- readRDS(coeff_file_cq)
c0 <- coeffs$mean_log_sum_circ
if (is.null(coeffs$V)) {
  coeffs$V <- matrix(c(
    coeffs$alpha_se^2, coeffs$cov_ab, coeffs$cov_ag,
    coeffs$cov_ab,      coeffs$beta_se^2,  coeffs$cov_bg,
    coeffs$cov_ag,      coeffs$cov_bg,     coeffs$gamma_se^2
  ), nrow = 3, byrow = TRUE)
}
map_logM <- function(logC) with(coeffs, alpha + beta*(logC - c0) + gamma*(logC - c0)^2)

## ---------------------------
## 4) Threshold choice: MRL + stability scans (unchanged)
## ---------------------------
make_thresholds <- function(y, q_from=0.75, q_to=0.90, n=50){
  rng <- quantile(y, c(q_from, q_to), na.rm=TRUE)
  seq(rng[1], rng[2], length.out=n)
}
fit_gpd_at_u <- function(y, u) {
  fit <- fevd(y, type="GP", threshold=u)
  par <- fit$results$par
  list(fit=fit, scale=unname(par["scale"]), shape=unname(par["shape"]),
       cov=summary(fit)$cov.theta)
}
adj_scale <- function(scale, shape, u, u0) scale - shape*(u - u0)

u_seq <- make_thresholds(y, q_from=0.40, q_to=0.90, n=50)
u0    <- as.numeric(quantile(y, thresh_q_opt, na.rm=TRUE))

# quick MRL
mrl_data <- function(y, u_seq){
  sapply(u_seq, function(u){
    ex <- y[y>u] - u
    if (length(ex) < 5) return(NA_real_)
    mean(ex)
  })
}
mrl_vals <- mrl_data(y, u_seq)

fit0 <- fit_gpd_at_u(y, u0)
shape_hat0 <- fit0$shape
adj_hat0   <- fit0$scale

zcrit <- qnorm(1 - (1 - ci_level)/2)

shape <- scale <- shape_lo <- shape_hi <- adj <- adj_lo <- adj_hi <- rep(NA_real_, length(u_seq))
for (i in seq_along(u_seq)) {
  u <- u_seq[i]
  ex <- y[y > u] - u
  if (length(ex) < 20) next
  out <- try(fit_gpd_at_u(y, u), silent=TRUE)
  if (inherits(out,"try-error")) next
  scale[i] <- out$scale; shape[i] <- out$shape
  se_scale <- sqrt(out$cov[1,1]); se_shape <- sqrt(out$cov[2,2])
  shape_lo[i] <- shape[i] - zcrit*se_shape
  shape_hi[i] <- shape[i] + zcrit*se_shape
  adj[i] <- adj_scale(scale[i], shape[i], u, u0)
  var_adj <- out$cov[1,1] + (u - u0)^2*out$cov[2,2] - 2*(u - u0)*out$cov[1,2]
  se_adj  <- sqrt(max(var_adj,0))
  adj_lo[i] <- adj[i] - zcrit*se_adj
  adj_hi[i] <- adj[i] + zcrit*se_adj
}

p_mrl <- ggplot(data.frame(u=u_seq, mrl=mrl_vals), aes(u, mrl)) +
  geom_line() + geom_point() +
  geom_vline(xintercept = u0, linetype = "dashed", color = "red") +
  labs(title = "MRL plot (log CFH)", x = "Threshold (log CFH)", y = "Mean excess") +
  theme_science_polished
p_shape <- ggplot(data.frame(u=u_seq, shape, shape_lo, shape_hi), aes(u, shape)) +
  geom_point() +
  geom_errorbar(aes(ymin=shape_lo, ymax=shape_hi), width=0.03, color="blue") +
  geom_vline(xintercept=u0, color="red", linetype="dashed") +
  geom_hline(yintercept=shape_hat0, color="red", linetype="dashed") +
  labs(x="Threshold (log CFH)", y="Shape (xi)") +
  theme_science_polished
p_adj <- ggplot(data.frame(u=u_seq, adj, adj_lo, adj_hi), aes(u, adj)) +
  geom_point() +
  geom_errorbar(aes(ymin=adj_lo, ymax=adj_hi), width=0.03, color="blue") +
  geom_vline(xintercept=u0, color="red", linetype="dashed") +
  geom_hline(yintercept=adj_hat0, color="red", linetype="dashed") +
  labs(x="Threshold (log CFH)", y="Adjusted scale") +
  theme_science_polished

ggsave(file.path(FIG_DIR, "cfh_log_mrl.png"), p_mrl, dpi=600, w=7, h=5, units="in")
ggsave(file.path(FIG_DIR, "cfh_log_shape_stability.png"), p_shape, dpi=600, w=7, h=5, units="in")
ggsave(file.path(FIG_DIR, "cfh_log_adj_scale_stability.png"), p_adj,   dpi=600, w=7, h=5, units="in")

## ---------------------------
## 5) GOF diagnostics (unchanged)
## ---------------------------
diagnostic_plots <- function(y, u, scale_hat, shape_hat, label="log CFH") {
  exc <- y[y>u] - u; n <- length(exc)
  probs <- ppoints(n)
  theo_q <- if (abs(shape_hat) > 1e-10) { u + scale_hat/shape_hat * (probs^(-shape_hat) - 1)
  } else { u - scale_hat * log(probs) }
  dfqq <- data.frame(Theoretical = rev(theo_q), Empirical = sort(y[y>u]))
  pqq <- ggplot(dfqq, aes(Theoretical, Empirical)) +
    geom_point(color="steelblue") +
    geom_abline(slope=1, intercept=0, linetype="dashed", color="red") +
    labs(x="Theoretical quantiles", y="Empirical quantiles") + theme_science_polished
  F_theo <- if (abs(shape_hat)>1e-10) { 1 - (1 + shape_hat*exc/scale_hat)^(-1/shape_hat)
  } else { 1 - exp(-exc/scale_hat) }
  dfpp <- data.frame(Theoretical=sort(F_theo), Empirical=(1:n)/n)
  ppp <- ggplot(dfpp, aes(Theoretical, Empirical)) +
    geom_point(color="darkgreen") +
    geom_abline(slope=1, intercept=0, linetype="dashed", color="red") +
    labs(x="Theoretical CDF", y="Empirical CDF") + theme_science_polished
  F_hat <- pevd(exc, scale=scale_hat, shape=shape_hat, type="GP")
  pks <- ggplot(data.frame(F_hat=F_hat), aes(F_hat)) +
    geom_histogram(aes(y=..density..), bins=20, fill="skyblue", color="black", alpha=0.7) +
    geom_hline(yintercept=1, color="red", linetype="dashed") +
    geom_density(color="darkblue", linewidth = 1.1, adjust=1.5) +
    labs(title=paste0("Uniformity: PIT (", label, ")"),
         x=expression(hat(F)(y)), y="Density") + theme_science_polished
  list(pqq=pqq, ppp=ppp, pks=pks)
}
diags <- diagnostic_plots(y, u0, fit0$scale, fit0$shape, "log CFH")
ggsave(file.path(FIG_DIR, "cfh_log_qq.png"),  diags$pqq, dpi=600, w=7, h=5, units="in")
ggsave(file.path(FIG_DIR, "cfh_log_pp.png"),  diags$ppp, dpi=600, w=7, h=5, units="in")
ggsave(file.path(FIG_DIR, "cfh_log_pit_uniformity.png"), diags$pks, dpi=600, w=7, h=5, units="in")

## ---------------------------
## 6) Endpoint posterior (Normal prior on y* centered at u0+1.05) + prior overlays
## ---------------------------

# Exceedances at the chosen anchor u0
y_above <- y[y > u0]; n_ex <- length(y_above)
stopifnot(n_ex >= 5)

Ymax <- max(y_above)
eps  <- max(1e-10, 1e-6 * abs(Ymax))
L    <- Ymax + eps

# Prior center on y*
U_center <- u0 + prior_center_delta

# Integration grid for y* (cover prior mass + likelihood tail; ensure > L)
right_pad <- 4 * prior_sd_y
y_star_grid <- seq(L, U_center + right_pad, length.out = 1200)
xi_grid     <- seq(-1.0, -0.02, length.out = 1400)

# Priors
logprior_y <- function(y_star) {
  if (y_star <= L) return(-Inf)               # enforce endpoint > max exceedance
  dnorm(y_star, mean = U_center, sd = prior_sd_y, log = TRUE)
}
logprior_xi <- switch(
  xi_prior_mode,
  "improper" = function(xi) if (xi < 0) 0 else -Inf,
  "pc"       = function(xi) if (xi < 0 && xi > -1) log(lambda_k) - lambda_k*(-xi) else -Inf,
  stop("Unknown xi_prior_mode")
)

# Weibull-domain log-likelihood in (y*, ξ)
loglik_endpoint <- function(y_star, xi, u, y_ex) {
  if (xi >= 0 || any(y_star <= y_ex)) return(-Inf)
  n <- length(y_ex)
  n*log(y_star - u)/xi - (1/xi + 1)*sum(log(y_star - y_ex)) - n*log(-xi)
}

# Marginalize ξ
marg_log_post <- vapply(y_star_grid, function(ys){
  lv <- vapply(xi_grid, function(xi){
    loglik_endpoint(ys, xi, u0, y_above) + logprior_xi(xi)
  }, numeric(1))
  m <- max(lv)
  log(pracma::trapz(xi_grid, exp(lv - m))) + m + logprior_y(ys)
}, numeric(1))

post_unnorm <- exp(marg_log_post - max(marg_log_post))
Z_norm      <- pracma::trapz(y_star_grid, post_unnorm)
stopifnot(is.finite(Z_norm) && Z_norm > 0)
post_y      <- post_unnorm / Z_norm

# --- Prior densities for overlays ---
## --- Build FULL prior grids for overlays (no truncation at L) ---
if (plot_prior) {
  k_left  <- 5                    # how many SDs left of the center to show
  k_right <- 5                    # how many SDs right of the center to show
  y_prior_grid <- seq(U_center - k_left*prior_sd_y,
                      U_center + k_right*prior_sd_y,
                      length.out = 2000)
  
  # Full (untruncated) prior on y*
  prior_y_full <- dnorm(y_prior_grid, mean = U_center, sd = prior_sd_y)
  # Normalize for plotting (so it's a proper density over the displayed grid)
  Zp_y <- pracma::trapz(y_prior_grid, prior_y_full)
  prior_y_full <- if (Zp_y > 0) prior_y_full / Zp_y else prior_y_full
  
  # Transform to CFH (mm): x = exp(y), f_C(x) = f_Y(log x) / x
  cfh_prior_grid <- exp(y_prior_grid)
  prior_cfh_full <- prior_y_full / cfh_prior_grid
  Zp_c <- pracma::trapz(cfh_prior_grid, prior_cfh_full)
  prior_cfh_full <- if (Zp_c > 0) prior_cfh_full / Zp_c else prior_cfh_full
}


# Draws & summaries on CFH scale (mm)
ystar_draws <- sample(y_star_grid, size = n_post_ystar, replace = TRUE, prob = post_y)
cfh_star_mm <- exp(ystar_draws)
cfh_map_mm  <- exp(y_star_grid[which.max(post_y)])
cfh_upper_mm <- as.numeric(quantile(cfh_star_mm, probs = ci_level))  # one-sided upper

# --- CFH* plot (posterior + optional prior) ---
# --- CFH* plot (posterior + full prior) ---
df_post_cfh <- data.frame(C = exp(ystar_draws))
p_cfh_star <- ggplot(df_post_cfh, aes(C)) +
  geom_density(color = "darkblue", linewidth = 1.2) +
  { if (plot_prior) geom_line(
    data = data.frame(C = cfh_prior_grid, d = prior_cfh_full),
    aes(x = C, y = d),
    color = "grey40", linetype = "dotted", linewidth = 1.0
  ) else NULL } +
  geom_vline(xintercept = exp(y_star_grid[which.max(post_y)]),
             color = "purple", linetype = "dashed", linewidth = 1.0) +
  geom_vline(xintercept = as.numeric(quantile(exp(ystar_draws), probs = ci_level)),
             color = "orange", linetype = "dotdash", linewidth = 1.0) +
  labs(
    title = sprintf("CFH* posterior (ξ prior: %s; prior μ = u0+%.2f, σ = %.2f)%s",
                    xi_prior_mode, prior_center_delta, prior_sd_y,
                    if (plot_prior) " — grey dotted: prior" else ""),
    x = "CFH* (mm)", y = "Density"
  ) +
  theme_science_polished
p_cfh_star
ggsave(file.path(FIG_DIR, "cfh_endpoint_posterior.png"),
       p_cfh_star, dpi = 600, width = 7, height = 5, units = "in")

# --- y* (log CFH) plot (posterior + full prior) ---
df_post_y <- data.frame(y = ystar_draws)
# expand x-limits to show the full prior (which may extend below L)
xlo <- min(min(y_star_grid), min(y_prior_grid))
xhi <- max(max(y_star_grid), max(y_prior_grid))

p_y_star <- ggplot(df_post_y, aes(y)) +
  geom_density(color = "darkblue", linewidth = 1.2) +
  { if (plot_prior) geom_line(
    data = data.frame(y = y_prior_grid, d = prior_y_full),
    aes(x = y, y = d),
    color = "grey40", linetype = "dotted", linewidth = 1.0
  ) else NULL } +
  geom_vline(xintercept = y_star_grid[which.max(post_y)],
             color = "purple", linetype = "dashed", linewidth = 1.0) +
  geom_vline(xintercept = as.numeric(quantile(ystar_draws, probs = ci_level)),
             color = "orange", linetype = "dotdash", linewidth = 1.0) +
  coord_cartesian(xlim = c(xlo, xhi)) +
  labs(
    title = sprintf("y* = log CFH* posterior (ξ prior: %s; prior μ = u0+%.2f, σ = %.2f)%s",
                    xi_prior_mode, prior_center_delta, prior_sd_y,
                    if (plot_prior) " — grey dotted: prior" else ""),
    x = expression(y^"*  (log~CFH)"), y = "Density"
  ) +
  theme_science_polished
p_y_star
ggsave(file.path(FIG_DIR, "cfh_endpoint_posterior_log_scale.png"),
       p_y_star, dpi = 600, width = 7, height = 5, units = "in")

ggsave(file.path(FIG_DIR, "cfh_endpoint_posterior_log_scale.png"),
       p_y_star, dpi = 600, width = 7, height = 5, units = "in")

cat("\n--- CFH* endpoint (mm), Normal prior centered at u0+",
    prior_center_delta, " ---\n",
    "MAP:", round(cfh_map_mm, 1), "\n",
    "upper limit:", round(cfh_upper_mm, 1), "\n", sep="")

## ---------------------------
## 7) Mass extrapolation (propagate endpoint through allometry)
## ---------------------------
n_post_ystar_big <- 100000
abg_samps <- MASS::mvrnorm(n_post_ystar_big,
                           mu=c(coeffs$alpha, coeffs$beta, coeffs$gamma),
                           Sigma=coeffs$V)
alpha_s <- abg_samps[,1]; beta_s <- abg_samps[,2]; gamma_s <- abg_samps[,3]

# Resample y* draws to match length
ystar_draws_big <- sample(ystar_draws, size = n_post_ystar_big, replace = TRUE)

logM_draws <- alpha_s + beta_s*(ystar_draws_big - coeffs$mean_log_sum_circ) +
  gamma_s*(ystar_draws_big - coeffs$mean_log_sum_circ)^2

# Optional residual uncertainty for predictive upper individual mass endpoint
add_predictive <- TRUE
q_indiv <- 0.95
if (add_predictive) {
  stopifnot(!is.null(coeffs$resid_sd))
  logM_draws <- logM_draws + qnorm(q_indiv) * coeffs$resid_sd
}

mass_draws_t <- exp(logM_draws) / 1e6
dens <- density(mass_draws_t, n = 4096, adjust = 1.0)
df_mass <- data.frame(M = dens$x, d = dens$y)
mass_map_t <- df_mass$M[which.max(df_mass$d)]
mass_up95_t <- as.numeric(quantile(mass_draws_t, probs = ci_level))

p_mass <- ggplot(df_mass, aes(M, d)) +
  geom_line(linewidth = 1.2) +
  geom_vline(xintercept = mass_map_t,  color = "purple", linetype = "dashed",  linewidth = 1) +
  geom_vline(xintercept = mass_up95_t, color = "orange", linetype = "dotdash", linewidth = 1) +
  coord_cartesian(xlim = c(0, 1000)) +
  labs(x = "Mass (tons)", y = "Density") +
  theme_science_polished
p_mass
ggsave(file.path(FIG_DIR, "mass_endpoint_posterior.png"),
       p_mass, dpi = 600, w = 7, h = 5, units = "in")

cat(sprintf("\nMass endpoint (tons)%s:\n",
            if (add_predictive) paste0(" — predictive ", q_indiv*100, "% individual") else ""),
    "MAP:", round(mass_map_t,1), " | upper limit:", round(mass_up95_t,1), "\n")

