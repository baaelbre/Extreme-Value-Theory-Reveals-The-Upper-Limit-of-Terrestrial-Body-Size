############################################################
#  Title: Hierarchical-style EVT for sauropod (single set)
#         — Stan MCMC endpoint with hierarchical hyperpriors
#         — Optional PC prior on kappa=-xi
#  Author: Bastiaan Aelbrecht (template by ChatGPT)
#  Date  : 2025-09-09
############################################################
library(ggplot2)
library(dplyr)
library(readxl)
library(pracma)   # trapz, cumtrapz
library(rstan)

rstan_options(auto_write = TRUE)
options(mc.cores = max(1, parallel::detectCores() - 1))
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

ci_level <- 0.95

## ---------------------------
## 1) Data ingest (single dataset)
## ---------------------------
df <- read_excel("Data/DEmic23_updated_Supplemental_Data_withPubYear.xlsx")
df <- df[, c("hum+fem circ (mm)", "publication year")]
colnames(df) <- c("sum_circ_mm", "year")
df <- na.omit(df)
df$log_cfh <- log(as.double(df$sum_circ_mm))
y <- df$log_cfh
stopifnot(length(y) >= 80)

## ---------------------------
## 2) Threshold (choose as you did before)
## ---------------------------
thresh_q_opt  <- 0.78
u0   <- as.numeric(quantile(y, thresh_q_opt, na.rm = TRUE))
y_ex <- y[y > u0]
stopifnot(length(y_ex) >= 20)
ymax <- max(y_ex)

## ---------------------------
## 3) Hyperpriors (hierarchical, like gators)
## ---------------------------
# We'll center the location hyperprior for mu_y at the Hokkanen mid within [L_y, U_y]
# where L_y is just above ymax and U_y reflects a biological cap.
L_eps <- max(1e-10, 1e-6 * abs(ymax))
L_y   <- ymax + L_eps
U_y   <- u0 + 1.25       # your previous upper (≈ “600 t” on your mapping)

mu_y0 <- (L_y + U_y)/2   # soft center on log CFH*
sd_y0 <- (U_y - L_y)/4   # broad (≈95% mass in [L_y, U_y] if Normal)
s_y   <- 0.50            # half-Cauchy scale for tau_y (between-dataset spread analog)

# Shape hyperpriors (Weibull domain)
mu_xi0 <- -0.20
sd_xi0 <-  0.20
s_xi   <-  0.50           # half-Cauchy scale for tau_xi

# Optional PC prior on kappa=-xi (toggle). If TRUE, xi uses PC prior instead of hierarchical Normal.
use_pc_xi <- 0L
pc_lambda <- 3.0

## ---------------------------
## 4) Stan model (single-dataset with hyperpriors)
## ---------------------------
stan_code <- "
functions {
  real gp_weibull_lpdf(real y, real u, real xi, real sigma) {
    // y > u, xi<0, sigma>0, support: 0 <= (y-u) <= -sigma/xi
    real x = y - u;
    real z = 1 + xi * x / sigma;
    if (x < 0) return negative_infinity();
    if (z <= 0) return negative_infinity();
    return -log(sigma) - (1/xi + 1) * log(z);
  }
}
data {
  int<lower=1> N;           // exceedances
  vector[N] y;              // on original log CFH scale
  real u;                   // threshold
  real ymax;                // max(y) above u (lower bound for y*)

  // Hyperprior centers/scales (like gators)
  real mu_y0;
  real<lower=0> sd_y0;
  real<lower=0> s_y;        // half-Cauchy scale for tau_y

  real mu_xi0;
  real<lower=0> sd_xi0;
  real<lower=0> s_xi;       // half-Cauchy scale for tau_xi

  // Optional PC prior toggle & rate (overrides hierarchical xi if enabled)
  int<lower=0,upper=1> use_pc_xi;
  real<lower=0> pc_lambda;
}
parameters {
  // Endpoint & shape (group-level, single group)
  real<lower=-1,upper=0> xi;
  real ystar_raw;                 // y* = ymax + exp(raw)

  // Hyperparameters (population-level)
  real mu_y;
  real<lower=0> tau_y;

  real mu_xi;
  real<lower=0> tau_xi;
}
transformed parameters {
  real ystar = ymax + exp(ystar_raw);  // ensures y* > ymax
  real sigma = (u - ystar) * xi;       // >0 because (u - y*)<0 and xi<0
}
model {
  // Hyperpriors
  mu_y  ~ normal(mu_y0,  sd_y0);
  tau_y ~ cauchy(0, s_y);

  if (use_pc_xi == 1) {
    // Use PC prior directly on xi via kappa = -xi (handled below in 'priors on xi')
    // We still define mu_xi, tau_xi but do not use them (weak anchoring avoided).
  } else {
    mu_xi  ~ normal(mu_xi0, sd_xi0);
    tau_xi ~ cauchy(0, s_xi);
  }

  // Priors on y* (hierarchical Normal) + Jacobian of exp()
  ystar ~ normal(mu_y, tau_y);
  target += ystar_raw;

  // Priors on xi
  if (use_pc_xi == 1) {
    // PC prior on kappa = -xi in (0,1): Exp(lambda), truncated to (0,1)
    real kappa = -xi;
    if (kappa <= 0 || kappa >= 1) target += negative_infinity();
    else {
      real logZ = log1m_exp(-pc_lambda); // log(1 - e^{-lambda})
      target += exponential_lpdf(kappa | pc_lambda) - logZ;
    }
  } else {
    xi ~ normal(mu_xi, tau_xi); // bounded to (-1,0) by declaration
  }

  // Likelihood
  for (n in 1:N) target += gp_weibull_lpdf(y[n] | u, xi, sigma);
}
generated quantities {
  real cfh_star_mm = exp(ystar);
  vector[N] pit;
  for (n in 1:N) {
    real x = y[n] - u;
    real F = (xi != 0) ? (1 - pow(1 + xi * x / sigma, -1/xi))
                       : (1 - exp(-x / sigma));
    pit[n] = fmin(fmax(F, 1e-12), 1 - 1e-12);
  }
}
"

## ---------------------------
## 5) Fit
## ---------------------------
stan_dat <- list(
  N = length(y_ex),
  y = y_ex,
  u = u0,
  ymax = ymax,
  mu_y0 = mu_y0,  sd_y0 = sd_y0,  s_y = s_y,
  mu_xi0 = mu_xi0, sd_xi0 = sd_xi0, s_xi = s_xi,
  use_pc_xi = use_pc_xi, pc_lambda = pc_lambda
)

sm  <- stan_model(model_code = stan_code)
fit <- sampling(
  sm, data = stan_dat,
  chains = 4, iter = 4000, warmup = 2000, thin = 1,
  control = list(adapt_delta = 0.96, max_treedepth = 12),
  seed = 42
)

print(fit, pars = c("xi","ystar","sigma","cfh_star_mm","mu_y","tau_y","mu_xi","tau_xi"),
      probs = c(0.1, 0.5, 0.9))

## ---------------------------
## 6) Prior–posterior overlay on CFH*
## ---------------------------
post <- rstan::extract(fit)
ystar_draws   <- post$ystar
cfh_star_draws<- post$cfh_star_mm

# A simple prior overlay: map Normal(mu_y0, sd_y0) to CFH* scale (ignores tau_y mixing for clarity)
y_grid     <- seq(mu_y0 - 4*sd_y0, mu_y0 + 4*sd_y0, length.out = 2000)
prior_y    <- dnorm(y_grid, mean = mu_y0, sd = sd_y0)
cfh_grid   <- exp(y_grid)
prior_cfh  <- prior_y / cfh_grid
prior_cfh  <- prior_cfh / pracma::trapz(cfh_grid, prior_cfh)

dens_post  <- density(cfh_star_draws, n = 2048, adjust = 1)
prior_df   <- data.frame(c = cfh_grid, d = prior_cfh)
post_df    <- data.frame(c = dens_post$x, d = dens_post$y)

p_cfh <- ggplot() +
  geom_area(data = prior_df, aes(c, d), fill = "grey70", alpha = 0.18) +
  geom_line(data = prior_df, aes(c, d), color = "grey50", linewidth = 0.8) +
  geom_line(data = post_df,  aes(c, d), color = "darkblue", linewidth = 1.2) +
  geom_vline(xintercept = median(cfh_star_draws),      color = "purple", linetype = "dashed") +
  geom_vline(xintercept = quantile(cfh_star_draws,.95), color = "orange", linetype = "dotdash") +
  labs(x = "CFH* (mm)", y = "Density",
       caption = sprintf("Grey: prior on exp(μ_y). Blue: posterior. u=%.2f (log CFH).", u0)) +
  theme_science_polished
print(p_cfh)
ggsave(file.path(FIG_DIR, "cfh_endpoint_posterior_stan_hier.png"),
       p_cfh, dpi = 600, width = 7, height = 5, units = "in")

cat(sprintf("\nCFH* (mm): median=%.1f, 90%%=%.1f, 95%%=%.1f\n",
            quantile(cfh_star_draws,.5), quantile(cfh_star_draws,.9), quantile(cfh_star_draws,.95)))
cat(sprintf("xi: median=%.3f  [10%%,90%%]=[%.3f, %.3f]\n",
            median(post$xi), quantile(post$xi,.1), quantile(post$xi,.9)))

## ---------------------------
## 7) Mass extrapolation via numeric integration (S7)
##      — centered quadratic allometry; NO sampling
## ---------------------------
# Load your allometry coefficients (centered quadratic)
coeffs <- readRDS("centered_quadratic_coefficients.rds")
alpha  <- as.numeric(coeffs$alpha)
beta   <- as.numeric(coeffs$beta)
gamma  <- as.numeric(coeffs$gamma)
c0     <- as.numeric(coeffs$mean_log_sum_circ)
Vabg   <- matrix(c(
  coeffs$alpha_se^2, coeffs$cov_ab,      coeffs$cov_ag,
  coeffs$cov_ab,     coeffs$beta_se^2,   coeffs$cov_bg,
  coeffs$cov_ag,     coeffs$cov_bg,      coeffs$gamma_se^2
), nrow = 3, byrow = TRUE)
sigma_e<- as.numeric(coeffs$resid_sd)  # predictive residual SD on log mass

# Build a posterior *density* for y* from Stan draws
dens_y   <- density(ystar_draws, n = 4096, adjust = 1)
y_star_grid <- dens_y$x
post_y      <- pmax(dens_y$y, 0)
# normalize (defensive)
post_y <- post_y / pracma::trapz(y_star_grid, post_y)

# Helper: trapezoid weights that sum to 1 (for mixtures)
trapz_weights <- function(x, f_un) {
  stopifnot(length(x) == length(f_un), length(x) >= 2)
  dx <- diff(x)
  mids <- (f_un[-1] + f_un[-length(f_un)]) * dx / 2
  total <- sum(mids); if (!is.finite(total) || total <= 0) stop("trapz_weights: non-positive mass.")
  w <- numeric(length(x))
  w[1] <- mids[1]/2
  w[length(x)] <- mids[length(mids)]/2
  if (length(x) > 2) w[2:(length(x)-1)] <- (mids[-length(mids)] + mids[-1])/2
  w / sum(w)
}
cdf_from_density <- function(x, f) {
  f[!is.finite(f)] <- 0
  F <- pracma::cumtrapz(x, f)
  F / max(F, na.rm = TRUE)
}

W_y <- trapz_weights(y_star_grid, post_y)

# Predictive mean/var for log mass at each y* (centered quadratic)
z_i    <- y_star_grid - c0
mu_i   <- alpha + beta*z_i + gamma*z_i^2
X_mat  <- cbind(1, z_i, z_i^2)
quad   <- rowSums((X_mat %*% Vabg) * X_mat)
add_predictive <- TRUE
sig2_i <- if (add_predictive) pmax(quad + sigma_e^2, 0) else pmax(quad, 0)

# Build mixture on m = log mass
m_min <- min(mu_i - 6*sqrt(sig2_i))
m_max <- max(mu_i + 6*sqrt(sig2_i))
m_grid <- seq(m_min, m_max, length.out = 4000)

dnorm_vec <- function(x, mean, sd) (1/(sd*sqrt(2*pi))) * exp(-0.5*((x-mean)/sd)^2)
f_m <- rowSums(sapply(seq_along(mu_i), function(i){
  sd_i <- sqrt(sig2_i[i]); if (!is.finite(sd_i) || sd_i <= 0) return(rep(0, length(m_grid)))
  W_y[i] * dnorm_vec(m_grid, mu_i[i], sd_i)
}))
f_m <- f_m / pracma::trapz(m_grid, f_m)

F_m   <- cdf_from_density(m_grid, f_m)
m_map <- m_grid[which.max(f_m)]
m_up95<- as.numeric(approx(x = F_m, y = m_grid, xout = ci_level, ties = "ordered", rule = 2)$y)

# Transform to tons and apply Jacobian
w_t_grid <- exp(m_grid)/1e6
f_wt     <- f_m / (exp(m_grid)/1e6)
f_wt     <- f_wt / pracma::trapz(w_t_grid, f_wt)
F_wt     <- cdf_from_density(w_t_grid, f_wt)

w_map_t  <- w_t_grid[which.max(f_wt)]
w_up95_t <- as.numeric(approx(x = F_wt, y = w_t_grid, xout = ci_level, ties = "ordered", rule = 2)$y)

cat(sprintf("\n--- Mass endpoint via numeric integration (S7) ---\n"))
cat(sprintf("MAP: %.1f tons | 95%% upper: %.1f\n", w_map_t, w_up95_t))

df_mass <- data.frame(M_t = w_t_grid, dens = f_wt)
p_mass  <- ggplot(df_mass, aes(M_t, dens)) +
  geom_line(linewidth = 1.2, color = "darkblue") +
  geom_vline(xintercept = w_map_t,  color = "purple",  linetype = "dashed") +
  geom_vline(xintercept = w_up95_t, color = "orange",  linetype = "dotdash") +
  labs(x = "Mass (tons)", y = "Density") +
  theme_science_polished + xlim(0, 700)
print(p_mass)
ggsave(file.path(FIG_DIR, "mass_endpoint_hier_stan.png"),
       p_mass, dpi = 600, width = 7, height = 5, units = "in")

