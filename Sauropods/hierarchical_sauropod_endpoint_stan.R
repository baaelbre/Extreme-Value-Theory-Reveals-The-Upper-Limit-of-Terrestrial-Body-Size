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
STAN_FILE  <- "stan/hier_gp_weibull_sauropod.stan"  # <- must be the NEW Stan file with tau_prior_lp()
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

# ---- Choose hyperpriors for tau_y and tau_xi (match Stan):
# 1=half-Cauchy, 2=half-Normal, 3=half-t_nu, 4=PC(Exp on SD),
# 5=LogNormal(on SD), 6=Gamma on precision (phi=1/tau^2), 7=Inv-Gamma on variance (v=tau^2)
PRIOR_TAU_Y  <- 1
PRIOR_TAU_XI <- 1 

# For convenience (labels in prints/plots)
prior_code_name <- function(k) {
  c("half-Cauchy","half-Normal","half-t_nu","PC (Exp on SD)",
    "LogNormal (SD)","Gamma on precision","Inv-Gamma on variance")[k]
}

# =============================================================================
# Helpers
# =============================================================================
kde_df <- function(x, from = quantile(x,0.001), to = quantile(x,0.999), n = 4096) {
  d <- density(x, from = from, to = to, n = n, kernel = "gaussian", adjust = 1.0)
  data.frame(t = d$x, d = d$y)
}
mode_from_kde <- function(x) { d <- kde_df(x); d$t[ which.max(d$d) ] }

# Invert centered quadratic allometry: log M = a + b Z + g Z^2,  Z = y* - c0
inv_mass_to_y <- function(m_tons, alpha, beta, gamma, c0) {
  m_log <- log(m_tons * 1e6) # tons -> grams
  A <- gamma; B <- beta; C <- alpha - m_log
  disc <- B*B - 4*A*C
  if (disc < 0) stop("Allometry inversion: no real root for this mass.")
  z1 <- (-B + sqrt(disc)) / (2*A)
  z2 <- (-B - sqrt(disc)) / (2*A)
  z  <- max(z1, z2)  # choose larger circumference solution
  z + c0             # y* = Z + c0
}

# ---- Calibration utilities for priors on SD tau
pc_lambda_from_Ualpha <- function(U, alpha) { -log(alpha) / U }
half_cauchy_scale_from_Ualpha <- function(U, alpha) {
  # For half-Cauchy(0, s): F(U) = (2/pi) arctan(U/s) = alpha  =>  s = U / tan(alpha*pi/2)
  U / tan(alpha * pi / 2)
}
half_normal_scale_from_Ualpha <- function(U, alpha) {
  # For half-Normal(0, s): F(U) = 2*Phi(U/s) - 1 = alpha  =>  s = U / qnorm((alpha+1)/2)
  U / qnorm((alpha + 1)/2)
}
half_t_scale_from_Ualpha <- function(U, alpha, nu) {
  # For half-t_nu(0, s): F(U) = 2*T_nu(U/s) - 1 = alpha  =>  s = U / qt((alpha+1)/2, df=nu)
  U / qt((alpha + 1)/2, df = nu)
}
lognormal_params_from_quantiles <- function(qlo, qhi, p = c(0.05, 0.95)) {
  # Solve ln-params so that P(tau <= qlo) = p[1], P(tau <= qhi)=p[2]
  z <- qnorm(p)
  sdlog <- (log(qhi) - log(qlo)) / (z[2] - z[1])
  meanlog <- log(qlo) - z[1]*sdlog
  list(meanlog = meanlog, sdlog = sdlog)
}

# ---- Random draws for general tau priors (match Stan choices)
rtau <- function(n, prior_code,
                 # half-Cauchy
                 hc_scale = 1,
                 # half-Normal
                 hn_scale = 1,
                 # half-t
                 ht_nu = 3, ht_scale = 1,
                 # PC prior (Exp on SD)
                 pc_lambda = 1,
                 # LogNormal
                 ln_mu = 0, ln_sd = 1,
                 # Gamma on precision
                 ga_shape = 2, ga_rate = 1,
                 # Inv-Gamma on variance
                 ig_shape = 2, ig_scale = 1) {
  if (prior_code == 1) {
    abs(rcauchy(n, location = 0, scale = hc_scale))
  } else if (prior_code == 2) {
    abs(rnorm(n, 0, hn_scale))
  } else if (prior_code == 3) {
    abs(rt(n, df = ht_nu) * ht_scale)
  } else if (prior_code == 4) {
    rexp(n, rate = pc_lambda)
  } else if (prior_code == 5) {
    rlnorm(n, meanlog = ln_mu, sdlog = ln_sd)
  } else if (prior_code == 6) {
    # tau from phi ~ Gamma => tau = 1/sqrt(phi)
    phi <- rgamma(n, shape = ga_shape, rate = ga_rate)
    1 / sqrt(phi)
  } else if (prior_code == 7) {
    # tau from v ~ Inv-Gamma => v = 1/Gamma(shape, rate=scale), tau = sqrt(v)
    v <- 1 / rgamma(n, shape = ig_shape, rate = ig_scale)
    sqrt(v)
  } else stop("Unknown prior_code for tau")
}

# Truncated normal sampler for xi in (-1,0)
rtruncnorm <- function(n, mu, sd, a = -1, b = 0) {
  alpha <- pnorm((a - mu)/sd)
  beta  <- pnorm((b - mu)/sd)
  u     <- runif(n, alpha, beta)
  qnorm(u) * sd + mu
}

# General prior sampler for (y*, xi) under chosen hyperpriors
sample_prior_yxi_general <- function(S,
                                     mu_y0, sd_y0,  # hyperprior for mu_y
                                     mu_xi0, sd_xi0,  # hyperprior for mu_xi
                                     # tau_y prior choices + params
                                     tauy_prior,
                                     hc_scale_y, hn_scale_y,
                                     ht_nu_y, ht_scale_y,
                                     pc_lambda_y,
                                     ln_mu_y, ln_sd_y,
                                     ga_shape_y, ga_rate_y,
                                     ig_shape_y, ig_scale_y,
                                     # tau_xi prior choices + params
                                     tauxi_prior,
                                     hc_scale_xi, hn_scale_xi,
                                     ht_nu_xi, ht_scale_xi,
                                     pc_lambda_xi,
                                     ln_mu_xi, ln_sd_xi,
                                     ga_shape_xi, ga_rate_xi,
                                     ig_shape_xi, ig_scale_xi) {
  mu_y   <- rnorm(S, mu_y0,  sd_y0)
  tau_y  <- rtau(S, tauy_prior, hc_scale_y, hn_scale_y, ht_nu_y, ht_scale_y,
                 pc_lambda_y, ln_mu_y, ln_sd_y, ga_shape_y, ga_rate_y, ig_shape_y, ig_scale_y)
  mu_xi  <- rnorm(S, mu_xi0, sd_xi0)
  tau_xi <- rtau(S, tauxi_prior, hc_scale_xi, hn_scale_xi, ht_nu_xi, ht_scale_xi,
                 pc_lambda_xi, ln_mu_xi, ln_sd_xi, ga_shape_xi, ga_rate_xi, ig_shape_xi, ig_scale_xi)
  
  ystar <- rnorm(S, mu_y,  tau_y)
  xi    <- rtruncnorm(S, mu_xi, tau_xi, a = -1, b = 0)
  xi    <- pmin(-1e-6, pmax(-0.999, xi))  # match Stan bounds
  list(ystar = ystar, xi = xi)
}

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

# Allometry coefficients (for elicitation + propagation)
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
# Hyper-hyperparameter elicitation from 100–1000 t belief (center at 316 t)
# =============================================================================
y_low  <- 7.72484
y_mid  <- 8.17305
y_high <- 8.62518
width_y <- y_high - y_low

mu_y0 <- y_mid
sd_y0 <- width_y / (2*1.96)     # ~ 95% for mu_y covers [low, high]
mu_xi0 <- -0.25
sd_xi0 <- 0.2

# ---- Calibrate tau priors (examples; tweak as you like) ----------------------
# Use a common interpretability rule: set U and alpha so P(tau > U) = alpha

# y* spread (tau_y) calibration
U_y     <- 0.25   # SD scale you rarely exceed
alpha_y <- 0.05

hc_scale_y <- half_cauchy_scale_from_Ualpha(U = U_y, alpha = 1 - alpha_y) # if using half-Cauchy
hn_scale_y <- half_normal_scale_from_Ualpha(U = U_y, alpha = 1 - alpha_y)
ht_nu_y    <- 3
ht_scale_y <- half_t_scale_from_Ualpha(U = U_y, alpha = 1 - alpha_y, nu = ht_nu_y)
pc_lambda_y<- 5
ln_par_y   <- lognormal_params_from_quantiles(qlo = 0.05, qhi = 0.40, p = c(0.05,0.95))
ln_mu_y    <- ln_par_y$meanlog; ln_sd_y <- ln_par_y$sdlog
ga_shape_y <- 2.5; ga_rate_y <- 10      # implies small tau_y on average
ig_shape_y <- 3.0; ig_scale_y <- 0.05   # places mass near small variance

# xi spread (tau_xi) calibration
U_xi      <- 0.06
alpha_xi  <- 0.05

hc_scale_xi <- half_cauchy_scale_from_Ualpha(U = U_xi, alpha = 1 - alpha_xi)
hn_scale_xi <- half_normal_scale_from_Ualpha(U = U_xi, alpha = 1 - alpha_xi)
ht_nu_xi    <- 5
ht_scale_xi <- half_t_scale_from_Ualpha(U = U_xi, alpha = 1 - alpha_xi, nu = ht_nu_xi)
pc_lambda_xi<- 5
ln_par_xi   <- lognormal_params_from_quantiles(qlo = 0.01, qhi = 0.10, p = c(0.05,0.95))
ln_mu_xi    <- ln_par_xi$meanlog; ln_sd_xi <- ln_par_xi$sdlog
ga_shape_xi <- 3.0; ga_rate_xi <- 120
ig_shape_xi <- 3.0; ig_scale_xi <- 0.01

cat(sprintf("\n[Hyper-hypers]\nmu_y0=%.6f, sd_y0=%.6f\nmu_xi0=%.6f, sd_xi0=%.6f\n", mu_y0, sd_y0, mu_xi0, sd_xi0))
cat(sprintf("[tau priors]  tau_y: %s | tau_xi: %s\n\n",
            prior_code_name(PRIOR_TAU_Y), prior_code_name(PRIOR_TAU_XI)))

# =============================================================================
# Stan data & model
# =============================================================================
stan_data <- list(
  N   = length(y_above),
  y   = as.vector(y_above),
  u   = u0,
  ymax = ymax,
  
  # means
  mu_y0       = mu_y0,
  sd_y0       = sd_y0,
  mu_xi0      = mu_xi0,
  sd_xi0      = sd_xi0,
  
  # prior codes
  tauy_prior  = PRIOR_TAU_Y,
  tauxi_prior = PRIOR_TAU_XI,
  
  # half-Cauchy scales
  hc_scale_y  = hc_scale_y,
  hc_scale_xi = hc_scale_xi,
  
  # half-Normal scales
  hn_scale_y  = hn_scale_y,
  hn_scale_xi = hn_scale_xi,
  
  # half-t: df and scale
  ht_nu_y     = ht_nu_y,
  ht_scale_y  = ht_scale_y,
  ht_nu_xi    = ht_nu_xi,
  ht_scale_xi = ht_scale_xi,
  
  # PC prior (Exp) lambdas
  pc_lambda_y  = pc_lambda_y,
  pc_lambda_xi = pc_lambda_xi,
  
  # LogNormal on SD
  ln_mu_y = ln_mu_y, ln_sd_y = ln_sd_y,
  ln_mu_xi = ln_mu_xi, ln_sd_xi = ln_sd_xi,
  
  # Gamma on precision (shape, rate)
  ga_shape_y = ga_shape_y, ga_rate_y = ga_rate_y,
  ga_shape_xi = ga_shape_xi, ga_rate_xi = ga_rate_xi,
  
  # Inv-Gamma on variance (shape, scale)
  ig_shape_y = ig_shape_y, ig_scale_y = ig_scale_y,
  ig_shape_xi = ig_shape_xi, ig_scale_xi = ig_scale_xi
)

sm  <- rstan::stan_model(file = STAN_FILE)
fit <- rstan::sampling(
  sm, data = stan_data, seed = 42,
  chains = 4, iter = 3000, warmup = 1500,
  control = list(adapt_delta = 0.95, max_treedepth = 12)
)

print(fit,
      pars  = c("ystar","xi","mu_y","tau_y","mu_xi","tau_xi","CFHstar_mm","CFHstar_pop_mm"),
      probs = c(0.05,0.5,0.95))

post <- as.data.frame(fit)
ystar_draw          <- post$ystar
xi_draw             <- post$xi
CFHstar_draw_mm     <- post$CFHstar_mm
CFHstar_pop_mm_draw <- post$CFHstar_pop_mm

# =============================================================================
# Prior vs Posterior — CFH*
# =============================================================================
set.seed(123)
pri <- sample_prior_yxi_general(
  S = 50000,
  mu_y0 = mu_y0, sd_y0 = sd_y0,
  mu_xi0 = mu_xi0, sd_xi0 = sd_xi0,
  # tau_y block
  tauy_prior = PRIOR_TAU_Y,
  hc_scale_y = hc_scale_y, hn_scale_y = hn_scale_y,
  ht_nu_y = ht_nu_y, ht_scale_y = ht_scale_y,
  pc_lambda_y = pc_lambda_y,
  ln_mu_y = ln_mu_y, ln_sd_y = ln_sd_y,
  ga_shape_y = ga_shape_y, ga_rate_y = ga_rate_y,
  ig_shape_y = ig_shape_y, ig_scale_y = ig_scale_y,
  # tau_xi block
  tauxi_prior = PRIOR_TAU_XI,
  hc_scale_xi = hc_scale_xi, hn_scale_xi = hn_scale_xi,
  ht_nu_xi = ht_nu_xi, ht_scale_xi = ht_scale_xi,
  pc_lambda_xi = pc_lambda_xi,
  ln_mu_xi = ln_mu_xi, ln_sd_xi = ln_sd_xi,
  ga_shape_xi = ga_shape_xi, ga_rate_xi = ga_rate_xi,
  ig_shape_xi = ig_shape_xi, ig_scale_xi = ig_scale_xi
)
CFHstar_prior_mm <- exp(pri$ystar)

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
  labs(title = sprintf("CFH* — prior (dashed) vs posterior (solid)\n tau_y: %s | tau_xi: %s",
                       prior_code_name(PRIOR_TAU_Y), prior_code_name(PRIOR_TAU_XI)),
       x = "CFH* (mm)", y = "Density") +
  xlim(80,600)+
  theme_science_polished
print(p_cfh)
ggsave(file.path(FIG_DIR, "cfh_endpoint_prior_vs_posterior_STAN.png"),
       p_cfh, dpi = 600, width = 7, height = 5, units = "in")

cat("\n--- CFH* endpoint (mm) ---\n")
cat(sprintf("Posterior MAP: %.1f | %d%% upper: %.1f mm\n\n",
            cfh_map_mm, round(100*ci_level), cfh_up_mm))

# =============================================================================
# Prior vs Posterior — Mass (tons)
# =============================================================================
# Posterior mass mixture
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

# Prior mass mixture
set.seed(123)
S <- 50000
ABG_prior <- MASS::mvrnorm(S, mu = c(alpha,beta,gamma), Sigma = Vabg)
Z_prior   <- pri$ystar - c0
eps_prior <- rnorm(S, 0, sigma_e)
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
  geom_line(data = dens_prior_m, aes(t, d), linewidth = 1.0, linetype = "dashed", color = "gray30") +
  geom_vline(xintercept = mass_map_t,  color = "purple",  linetype = "dashed",  linewidth = 1) +
  geom_vline(xintercept = mass_up95_t, color = "orange",  linetype = "dotdash", linewidth = 1) +
  labs(x = "Mass (tons)", y = "Density",
       title = sprintf("Maximum mass — prior (dashed) vs posterior (solid)\n tau_y: %s | tau_xi: %s",
                       prior_code_name(PRIOR_TAU_Y), prior_code_name(PRIOR_TAU_XI))) +
  coord_cartesian(xlim = c(0, max(700,
                                  quantile(mass_tons_post,  0.995, na.rm = TRUE),
                                  quantile(mass_tons_prior, 0.995, na.rm = TRUE)))) +
  theme_science_polished +
  annotate("text",
           x = mass_map_t, y = max(dens_post_m$d, na.rm = TRUE)*0.08,
           label = sprintf("%.0f t", mass_map_t),
           size = 4.5, fontface = "bold", vjust = -0.5, angle = 90, color = "purple")
print(p_mass)
ggsave(file.path(FIG_DIR, "mass_endpoint_prior_vs_posterior_STAN.png"),
       p_mass, dpi = 600, width = 7, height = 5, units = "in")

cat(sprintf("--- Mass endpoint (tons) ---\nPosterior MAP: %.1f | %d%% upper: %.1f\n\n",
            mass_map_t, round(100*ci_level), mass_up95_t))

# =============================================================================
# Optional: quick prior mass sanity check
# =============================================================================
cat(sprintf("[PRIOR CHECK] P(M*>1000t)=%.3f, P(M*<100t)=%.3f, median=%.0f t, 95%%=%.0f t\n",
            mean(mass_tons_prior > 1000),
            mean(mass_tons_prior < 100),
            unname(quantile(mass_tons_prior, 0.5)),
            unname(quantile(mass_tons_prior, 0.95))))

