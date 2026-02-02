suppressPackageStartupMessages({
  library(rstan)
  library(MASS)      # mvrnorm for allometry uncertainty
  library(ggplot2); library(dplyr); library(readxl); library(tidyr)
  library(scales);   library(grid)
})

rstan_options(auto_write = TRUE)
options(mc.cores = max(1, parallel::detectCores() - 1))
set.seed(42); options(stringsAsFactors = FALSE)

# =============================================================================
# Paths & plotting theme
# =============================================================================
DATA_XLSX  <- "Data/DEMic23_updated_Supplemental_Data_withPubYear.xlsx"
COEFFS_RDS <- "centered_quadratic_coefficients.rds"
# Stan file must be the NEW one with truncated hyperprior on mu_y (mu_y_lo/hi)
STAN_FILE  <- "stan/hier_gp_weibull_sauropod_uniform_2.stan"
FIG_DIR    <- "Figures"
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

# =============================================================================
# Config
# =============================================================================
q_circ   <- 0.78     # threshold quantile for u
ci_level <- 0.95     # one-sided upper quantile for summaries

# Truncated hyperprior bounds for mu_y (log-CFH)
MU_Y_LO <- 7.7
MU_Y_HI <- 8.62

# Hyperprior center & scale for mu_y (before truncation)
MU_Y0    <- 8.0     # ≈ 200 t on mass scale
SD_MU_Y  <- 0.30     # tune (e.g., 0.3–0.6) for more/less pooling

# Prior on xi (hierarchical: xi | mu_xi, tau_xi with mu_xi ~ Normal(mu_xi0, sd_xi0))
mu_xi0 <- -0.25
sd_xi0 <- 0.10

# Truncation support for dataset-level y* (model-enforced)
Y_MIN <- 0.0
Y_MAX <- 8.62         # legacy/user cap (safety)
Y_CAP <- 9.62         # hard cap (~1000 t mass scale)

# SD prior families for tau_y and tau_xi
# 1=half-Cauchy, 2=half-Normal, 3=half-t, 4=PC(Exp on SD), 5=LogNormal(SD),
# 6=Gamma(precision), 7=Inv-Gamma(variance)
PRIOR_TAU_Y  <- 1
PRIOR_TAU_XI <- 1

# Hyperparameters for SD priors (unused ones must still be >0)
hc_scale_y  <- 0.001; hc_scale_xi <- 0.001
hn_scale_y  <- 1.0;  hn_scale_xi <- 1.0
ht_nu_y     <- 3;    ht_scale_y  <- 0.30
ht_nu_xi    <- 5;    ht_scale_xi <- 0.10
pc_lambda_y <- 1.0;  pc_lambda_xi <- 1.0
ln_mu_y     <- log(0.25); ln_sd_y  <- 0.6
ln_mu_xi    <- log(0.08); ln_sd_xi <- 0.6
ga_shape_y  <- 2.5;  ga_rate_y <- 10
ga_shape_xi <- 3.0;  ga_rate_xi <- 120
ig_shape_y  <- 3.0;  ig_scale_y <- 0.05
ig_shape_xi <- 3.0;  ig_scale_xi <- 0.01

# Sampler controls
CHAINS <- 4; ITER <- 10000; WARMUP <- 2000
CONTROL <- list(adapt_delta = 0.995, max_treedepth = 12)

# =============================================================================
# Helpers
# =============================================================================
kde_df <- function(x, from = quantile(x,0.001), to = quantile(x,0.999), n = 4096) {
  d <- density(x, from = from, to = to, n = n, kernel = "gaussian", adjust = 1.0)
  data.frame(t = d$x, d = d$y)
}
mode_from_kde <- function(x) { d <- kde_df(x); d$t[ which.max(d$d) ] }

# Truncated normal sampler
rtruncnorm <- function(n, mu, sd, a = -Inf, b = Inf) {
  alpha <- pnorm((a - mu)/sd)
  beta  <- pnorm((b - mu)/sd)
  u     <- runif(n, alpha, beta)
  qnorm(u) * sd + mu
}
# Half-Cauchy sampler for SDs
rhalfcauchy <- function(n, scale) abs(rcauchy(n, 0, scale))

# Generic SD prior sampler (for prior predictive)
rtau <- function(n, prior_code,
                 hc_scale=1, hn_scale=1,
                 ht_nu=3, ht_scale=1,
                 pc_lambda=1,
                 ln_mu=0, ln_sd=1,
                 ga_shape=2, ga_rate=1,
                 ig_shape=2, ig_scale=1) {
  if (prior_code == 1)      return(rhalfcauchy(n, hc_scale))
  else if (prior_code == 2) return(abs(rnorm(n, 0, hn_scale)))
  else if (prior_code == 3) return(abs(rt(n, df = ht_nu) * ht_scale))
  else if (prior_code == 4) return(rexp(n, rate = pc_lambda))
  else if (prior_code == 5) return(rlnorm(n, ln_mu, ln_sd))
  else if (prior_code == 6) { phi <- rgamma(n, ga_shape, ga_rate); return(1/sqrt(phi)) }
  else if (prior_code == 7) { v <- 1/rgamma(n, ig_shape, ig_scale); return(sqrt(v)) }
  stop("Unknown prior_code")
}

# =============================================================================
# Data ingest (log CFH = log(sum femur + humerus, mm))
# =============================================================================
stopifnot(file.exists(DATA_XLSX), file.exists(COEFFS_RDS), file.exists(STAN_FILE))
df <- read_excel(DATA_XLSX) |>
  dplyr::select(`hum+fem circ (mm)`, `publication year`) |>
  dplyr::rename(sum_circ_mm = `hum+fem circ (mm)`,
                year        = `publication year`) |>
  tidyr::drop_na()

df$log_circ <- log(as.double(df$sum_circ_mm))
y <- df$log_circ
stopifnot(length(y) >= 80)

# Fixed global threshold
u0      <- as.numeric(quantile(y, q_circ, na.rm = TRUE))
y_above <- y[y > u0]
stopifnot(length(y_above) >= 25)
ymax    <- max(y_above)

# Allometry coefficients (for derived mass; no inversion)
coeffs <- readRDS(COEFFS_RDS)
if (!("V" %in% names(coeffs)) &&
    all(c("alpha_se","beta_se","gamma_se","cov_ab","cov_bg","cov_ag") %in% names(coeffs))) {
  coeffs$V <- matrix(c(
    coeffs$alpha_se^2, coeffs$cov_ab,       coeffs$cov_ag,
    coeffs$cov_ab,     coeffs$beta_se^2,    coeffs$cov_bg,
    coeffs$cov_ag,     coeffs$cov_bg,       coeffs$gamma_se^2
  ), 3, 3, byrow = TRUE)
}
stopifnot(all(c("alpha","beta","gamma","mean_log_sum_circ","V","resid_sd") %in% names(coeffs)))
alpha <- as.numeric(coeffs$alpha)
beta  <- as.numeric(coeffs$beta)
gamma <- as.numeric(coeffs$gamma)
c0    <- as.numeric(coeffs$mean_log_sum_circ)
Vabg  <- as.matrix(coeffs$V)
sigma_e <- as.numeric(coeffs$resid_sd)

# =============================================================================
# Build Stan data (Normal(mu_y0, sd_mu_y) with TRUNCATED hyperprior mu_y ∈ [MU_Y_LO, MU_Y_HI])
# =============================================================================
ystar_lb  <- max(u0, ymax) + 1e-6
y_lower   <- max(ystar_lb, Y_MIN)
y_upper   <- min(Y_CAP, Y_MAX)
stopifnot(y_upper > y_lower)

stan_data <- list(
  # data & support
  N    = length(y_above),
  y    = as.vector(y_above),
  u    = u0,
  ymax = ymax,
  y_min = Y_MIN,
  y_max = Y_MAX,
  y_cap = Y_CAP,
  
  # hyperprior on mu_y* (Normal truncated to [MU_Y_LO, MU_Y_HI])
  mu_y0   = MU_Y0,
  sd_mu_y = SD_MU_Y,
  mu_y_lo = MU_Y_LO,
  mu_y_hi = MU_Y_HI,
  
  # xi hyperprior (location) and SD prior choice
  mu_xi0 = mu_xi0,
  sd_xi0 = sd_xi0,
  
  tauy_prior  = PRIOR_TAU_Y,
  tauxi_prior = PRIOR_TAU_XI,
  
  # SD prior hyperparams
  hc_scale_y = hc_scale_y,   hc_scale_xi = hc_scale_xi,
  hn_scale_y = hn_scale_y,   hn_scale_xi = hn_scale_xi,
  ht_nu_y = ht_nu_y,         ht_scale_y = ht_scale_y,
  ht_nu_xi = ht_nu_xi,       ht_scale_xi = ht_scale_xi,
  pc_lambda_y = pc_lambda_y, pc_lambda_xi = pc_lambda_xi,
  ln_mu_y = ln_mu_y, ln_sd_y = ln_sd_y,
  ln_mu_xi = ln_mu_xi, ln_sd_xi = ln_sd_xi,
  ga_shape_y = ga_shape_y, ga_rate_y = ga_rate_y,
  ga_shape_xi = ga_shape_xi, ga_rate_xi = ga_rate_xi,
  ig_shape_y = ig_shape_y, ig_scale_y = ig_scale_y,
  ig_shape_xi = ig_shape_xi, ig_scale_xi = ig_scale_xi,
  
  # allometry (for derived mass output)
  alpha = alpha, beta = beta, gamma = gamma, c0 = c0
)

# =============================================================================
# Fit
# =============================================================================
sm  <- rstan::stan_model(file = STAN_FILE)
fit <- rstan::sampling(
  sm, data = stan_data, seed = 42,
  chains = CHAINS, iter = ITER, warmup = WARMUP, control = CONTROL
)

print(fit,
      pars  = c("ystar","xi","mu_y","tau_y","mu_xi","tau_xi","CFHstar_mm","Mstar_out"),
      probs = c(0.05,0.5,0.95))

post <- as.data.frame(fit)
ystar_draw      <- post$ystar
xi_draw         <- post$xi
CFHstar_draw_mm <- post$CFHstar_mm
Mstar_draw_t    <- post$Mstar_out
mu_y_draw       <- post$mu_y
tau_y_draw      <- post$tau_y

# =============================================================================
# Prior vs Posterior — CFH*
# =============================================================================
set.seed(123)
S <- 50000

# Prior predictive for mu_y (TRUNCATED Normal to [MU_Y_LO, MU_Y_HI])
mu_y_prior  <- rtruncnorm(S, MU_Y0, SD_MU_Y, a = MU_Y_LO, b = MU_Y_HI)
tau_y_prior <- rtau(S, PRIOR_TAU_Y,
                    hc_scale = hc_scale_y,
                    hn_scale = hn_scale_y,
                    ht_nu = ht_nu_y, ht_scale = ht_scale_y,
                    pc_lambda = pc_lambda_y,
                    ln_mu = ln_mu_y, ln_sd = ln_sd_y,
                    ga_shape = ga_shape_y, ga_rate = ga_rate_y,
                    ig_shape = ig_shape_y, ig_scale = ig_scale_y)

ystar_prior <- rtruncnorm(S, mu = mu_y_prior, sd = tau_y_prior, a = y_lower, b = y_upper)
CFHstar_prior_mm <- exp(ystar_prior)

rng_cfh <- range(
  quantile(CFHstar_draw_mm,  c(0.001, 0.999), na.rm = TRUE),
  quantile(CFHstar_prior_mm, c(0.001, 0.999), na.rm = TRUE)
)
dens_post_cfh  <- kde_df(CFHstar_draw_mm,  from = rng_cfh[1], to = rng_cfh[2])
dens_prior_cfh <- kde_df(CFHstar_prior_mm, from = rng_cfh[1], to = rng_cfh[2])

cfh_map_mm <- mode_from_kde(CFHstar_draw_mm)
cfh_up_mm  <- as.numeric(quantile(CFHstar_draw_mm, probs = ci_level))


# ----- Prior draws on LOG-scale, then exponentiate to original (mm)
mu_y_prior  <- rtruncnorm(S, MU_Y0, SD_MU_Y, a = MU_Y_LO, b = MU_Y_HI)
tau_y_prior <- rtau(
  S, PRIOR_TAU_Y,
  hc_scale = hc_scale_y, hn_scale = hn_scale_y,
  ht_nu = ht_nu_y, ht_scale = ht_scale_y,
  pc_lambda = pc_lambda_y,
  ln_mu = ln_mu_y, ln_sd = ln_sd_y,
  ga_shape = ga_shape_y, ga_rate = ga_rate_y,
  ig_shape = ig_shape_y, ig_scale = ig_scale_y
)
ystar_prior        <- rtruncnorm(S, mu = mu_y_prior, sd = tau_y_prior, a = y_lower, b = y_upper)
CFHstar_prior_mm   <- exp(ystar_prior)         # ORIGINAL scale (mm)
# Posterior draws are assumed available as: CFHstar_draw_mm (already on ORIGINAL scale)

# ----- Common x-range for KDEs (robust to tails)
rng_cfh <- range(
  quantile(CFHstar_draw_mm,  c(0.001, 0.999), na.rm = TRUE),
  quantile(CFHstar_prior_mm, c(0.001, 0.999), na.rm = TRUE),
  na.rm = TRUE
)

# ----- Smooth densities on ORIGINAL scale
dens_post_cfh  <- kde_df(CFHstar_draw_mm,  from = rng_cfh[1], to = rng_cfh[2], n = 4096)
dens_prior_cfh <- kde_df(CFHstar_prior_mm, from = rng_cfh[1], to = rng_cfh[2], n = 4096)

# ----- Summaries (ORIGINAL scale)
cfh_map_mm <- mode_from_kde(CFHstar_draw_mm)
cfh_up_mm  <- as.numeric(quantile(CFHstar_draw_mm, probs = ci_level))

# Nicely spaced annotation height
y_annot <- max(dens_post_cfh$d, dens_prior_cfh$d, na.rm = TRUE) * 0.10

# ----- Plot: shaded prior (grey) + posterior (blue), both on ORIGINAL scale (mm)
p_cfh <- ggplot() +
  # Prior: grey shaded area + grey outline
  geom_area(data = dens_prior_cfh, aes(t, d), fill = "grey75", alpha = 0.22) +
  geom_line(data = dens_prior_cfh, aes(t, d), color = "grey50", linewidth = 0.9) +
  # Posterior: blue line
  geom_line(data = dens_post_cfh,  aes(t, d), color = "darkblue", linewidth = 1.25) +
  # Vertical markers
  geom_vline(xintercept = cfh_map_mm, color = "purple",  linetype = "dashed",  linewidth = 1.05) +
  geom_vline(xintercept = cfh_up_mm,  color = "orange",  linetype = "dotdash", linewidth = 1.0) +
  # MAP label
  annotate("text", x = cfh_map_mm, y = y_annot, angle = 90, vjust = -0.5,
           label = sprintf("%.0f mm", cfh_map_mm),
           size = 4.2, fontface = "bold", color = "purple") +
  labs(
    x = "CFH* (mm)", y = "Density",
  ) +
  theme_science_polished

print(p_cfh)

# Optional: save
ggsave(file.path(FIG_DIR, "pop_cfh_endpoint.png"),
       p_cfh, dpi = 600, width = 7, height = 5, units = "in")

cat("\n--- CFH* endpoint (mm) ---\n")
cat(sprintf("Posterior MAP: %.1f | %d%% upper: %.1f mm\n\n",
            cfh_map_mm, round(100*ci_level), cfh_up_mm))

# =============================================================================
# Prior vs Posterior — Mass (tons) via allometry (predictive)
# =============================================================================
set.seed(42)
B    <- nrow(post)
take <- min(30000, B)
idx  <- sample.int(B, size = take, replace = TRUE)

ystar_samp <- ystar_draw[idx]
Z_post     <- ystar_samp - c0
ABG_post   <- MASS::mvrnorm(take, mu = c(alpha,beta,gamma), Sigma = Vabg)
eps_post   <- rnorm(take, 0, sigma_e)
m_log_post <- ABG_post[,1] + ABG_post[,2]*Z_post + ABG_post[,3]*Z_post^2 + eps_post
mass_tons_post <- exp(m_log_post) / 1e6

# Prior mass mixture (consistent with truncated hyperprior on mu_y)
ABG_prior   <- MASS::mvrnorm(S, mu = c(alpha,beta,gamma), Sigma = Vabg)
Z_prior     <- ystar_prior - c0
eps_prior   <- rnorm(S, 0, sigma_e)
m_log_prior <- ABG_prior[,1] + ABG_prior[,2]*Z_prior + ABG_prior[,3]*Z_prior^2 + eps_prior
mass_tons_prior <- exp(m_log_prior) / 1e6

rng_m <- range(
  quantile(mass_tons_post,  c(0.001, 0.999), na.rm = TRUE),
  quantile(mass_tons_prior, c(0.001, 0.999), na.rm = TRUE)
)
dens_post_m  <- kde_df(mass_tons_post,  from = rng_m[1], to = rng_m[2])
dens_prior_m <- kde_df(mass_tons_prior, from = rng_m[1], to = rng_m[2])

mass_map_t  <- mode_from_kde(mass_tons_post)
mass_up95_t <- as.numeric(quantile(mass_tons_post, probs = ci_level))

p_mass <- ggplot() +
  geom_line(data = dens_post_m,  aes(t, d), linewidth = 1.2, color = "darkblue") +
  geom_vline(xintercept = mass_map_t,  color = "purple",  linetype = "dashed",  linewidth = 1) +
  geom_vline(xintercept = mass_up95_t, color = "orange",  linetype = "dotdash", linewidth = 1) +
  labs(x = "Mass (tons)", y = "Density") +
  xlim(0,700) +
  theme_science_polished +
  annotate("text",
           x = mass_map_t, y = max(dens_post_m$d, na.rm = TRUE)*0.08,
           label = sprintf("%.0f t", mass_map_t),
           size = 4.5, fontface = "bold", vjust = -0.5, angle = 90, color = "purple")
print(p_mass)
ggsave(file.path(FIG_DIR, "mass_endpoint_posterior.png"),
       p_mass, dpi = 600, width = 7, height = 5, units = "in")

cat(sprintf("--- Mass endpoint (tons) ---\nPosterior MAP: %.1f | %d%% upper: %.1f\n\n",
            mass_map_t, round(100*ci_level), mass_up95_t))

