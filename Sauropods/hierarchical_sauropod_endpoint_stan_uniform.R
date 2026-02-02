suppressPackageStartupMessages({
  library(rstan)
  library(MASS)      # mvrnorm for allometry uncertainty
  library(ggplot2); library(dplyr); library(readxl)
  library(scales);   library(grid)
})

rstan_options(auto_write = TRUE)
options(mc.cores = max(1, parallel::detectCores() - 1))
set.seed(42); options(stringsAsFactors = FALSE)

# =============================================================================
# Paths & plotting theme
# =============================================================================
DATA_XLSX  <- "Data/DEmic23_updated_Supplemental_Data_withPubYear.xlsx"
COEFFS_RDS <- "centered_quadratic_coefficients.rds"
# Stan file must match the model with Normal(mu_y) hyperprior + y* cap at 8.62
STAN_FILE  <- "stan/hier_gp_weibull_sauropod_uniform.stan"
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

# Prior on xi (hierarchical: xi | mu_xi, tau_xi with mu_xi ~ Normal(mu_xi0, sd_xi0))
mu_xi0 <- -0.25
sd_xi0 <- 0.1

# Truncation support for dataset-level y* (model-enforced)
Y_MIN <- 0.0
Y_MAX <- 8.62         # legacy/user cap (safety)
Y_CAP <- 9.62        # hard prior cap corresponding to ~1000 t on mass scale

# Hyperprior for mu_y*: Normal(8.15, SD_MU_Y)
MU_Y0    <- 8 # 200 tons.
SD_MU_Y  <- 0.2     # tune as needed (e.g., 0.3–0.6)

# SD prior families for tau_y and tau_xi
# 1=half-Cauchy, 2=half-Normal, 3=half-t, 4=PC(Exp on SD), 5=LogNormal(SD),
# 6=Gamma(precision), 7=Inv-Gamma(variance)
PRIOR_TAU_Y  <- 1
PRIOR_TAU_XI <- 1

# Hyperparameters for SD priors (unused ones can be placeholders >0)
hc_scale_y  <- 0.10; hc_scale_xi <- 0.10
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

# =============================================================================
# Data ingest (log CFH = log(sum femur + humerus, mm))
# =============================================================================
stopifnot(file.exists(DATA_XLSX), file.exists(COEFFS_RDS), file.exists(STAN_FILE))
df <- read_excel(DATA_XLSX)
df <- df[, c("hum+fem circ (mm)", "publication year")]
colnames(df) <- c("sum_circ_mm", "year")
df <- na.omit(df)
df$log_circ <- log(as.double(df$sum_circ_mm))
y <- df$log_circ
stopifnot(length(y) >= 80)

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
# Build Stan data for Normal(mu_y) + truncated y* at y_upper=min(Y_CAP, Y_MAX)
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
  y_max = Y_MAX,   # kept for safety inside Stan
  y_cap = Y_CAP,   # hard prior cap (8.62)
  
  # hyperprior on mu_y* (Normal)
  mu_y0   = MU_Y0,
  sd_mu_y = SD_MU_Y,
  
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
  
  # allometry (only for derived mass output)
  alpha = alpha, beta = beta, gamma = gamma, c0 = c0
)

# =============================================================================
# Fit
# =============================================================================
sm  <- rstan::stan_model(file = STAN_FILE)
fit <- rstan::sampling(
  sm, data = stan_data, seed = 42,
  chains = 4, iter = 10000, warmup = 1500,
  control = list(adapt_delta = 0.95, max_treedepth = 12)
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

# Prior predictive for y*:
# mu_y ~ Normal(MU_Y0, SD_MU_Y)
# tau_y ~ chosen SD prior; y* ~ Normal(mu_y, tau_y) T[y_lower, y_upper]
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

mu_y_prior  <- rnorm(S, MU_Y0, SD_MU_Y)
tau_y_prior <- rtau(S, PRIOR_TAU_Y,
                    hc_scale = hc_scale_y,
                    hn_scale = hn_scale_y,
                    ht_nu = ht_nu_y, ht_scale = ht_scale_y,
                    pc_lambda = pc_lambda_y,
                    ln_mu = ln_mu_y, ln_sd = ln_sd_y,
                    ga_shape = ga_shape_y, ga_rate = ga_rate_y,
                    ig_shape = ig_shape_y, ig_scale = ig_scale_y)

ystar_prior <- rtruncnorm(S, mu = mu_y_prior, sd = tau_y_prior, a = y_lower, b = y_upper)

# Hierarchical xi prior
mu_xi_prior  <- rnorm(S, mu_xi0, sd_xi0)
tau_xi_prior <- rtau(S, PRIOR_TAU_XI,
                     hc_scale = hc_scale_xi,
                     hn_scale = hn_scale_xi,
                     ht_nu = ht_nu_xi, ht_scale = ht_scale_xi,
                     pc_lambda = pc_lambda_xi,
                     ln_mu = ln_mu_xi, ln_sd = ln_sd_xi,
                     ga_shape = ga_shape_xi, ga_rate = ga_rate_xi,
                     ig_shape = ig_shape_xi, ig_scale = ig_scale_xi)
xi_prior <- rtruncnorm(S, mu = mu_xi_prior, sd = tau_xi_prior, a = -1, b = 0)

CFHstar_prior_mm <- exp(ystar_prior)

rng_cfh <- range(
  quantile(CFHstar_draw_mm,  c(0.001, 0.999), na.rm = TRUE),
  quantile(CFHstar_prior_mm, c(0.001, 0.999), na.rm = TRUE)
)
dens_post_cfh  <- kde_df(CFHstar_draw_mm,  from = rng_cfh[1], to = rng_cfh[2])
dens_prior_cfh <- kde_df(CFHstar_prior_mm, from = rng_cfh[1], to = rng_cfh[2])

cfh_map_mm <- mode_from_kde(CFHstar_draw_mm)
cfh_up_mm  <- as.numeric(quantile(CFHstar_draw_mm, probs = ci_level))

p_cfh <- ggplot() +
  geom_line(data = dens_post_cfh,  aes(t, d), linewidth = 1.2, color = "darkblue") +
  geom_line(data = dens_prior_cfh, aes(t, d), linewidth = 1.0, linetype = "dashed", color = "gray30") +
  geom_vline(xintercept = cfh_map_mm, color = "purple",  linetype = "dashed",  linewidth = 1.0) +
  geom_vline(xintercept = cfh_up_mm,  color = "orange",  linetype = "dotdash", linewidth = 1.0) +
  labs(
    title = sprintf("CFH* — prior (dashed) vs posterior (solid)\nmu_y ~ N(%.2f, %.2f^2);  y* ∈ [%.3f, %.2f]",
                    MU_Y0, SD_MU_Y, y_lower, y_upper),
    x = "CFH* (mm)", y = "Density"
  ) +
  theme_science_polished
print(p_cfh)
ggsave(file.path(FIG_DIR, "cfh_endpoint_prior_vs_posterior_STAN_normal_mu_cap.png"),
       p_cfh, dpi = 600, width = 7, height = 5, units = "in")

cat("\n--- CFH* endpoint (mm) ---\n")
cat(sprintf("Posterior MAP: %.1f | %d%% upper: %.1f mm\n\n",
            cfh_map_mm, round(100*ci_level), cfh_up_mm))

# =============================================================================
# Prior vs Posterior — Mass (tons) via allometry (no inversion)
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

# Prior mass mixture (consistent with hierarchical y* prior & cap)
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
ggsave(file.path(FIG_DIR, "mass_endpoint_prior_vs_posterior_STAN_normal_mu_cap.png"),
       p_mass, dpi = 600, width = 7, height = 5, units = "in")

cat(sprintf("--- Mass endpoint (tons) ---\nPosterior MAP: %.1f | %d%% upper: %.1f\n\n",
            mass_map_t, round(100*ci_level), mass_up95_t))

