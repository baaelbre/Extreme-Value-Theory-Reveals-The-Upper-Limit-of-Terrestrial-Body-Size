############################################################
#  Title: Hierarchical-style EVT for log CFH (single dataset)
#         — finite-endpoint GP with hyperpriors (μ_y, τ_y, μ_ξ, τ_ξ)
#         — posterior for CFH* via quadrature
#         — propagate to mass with centered quadratic allometry
#  Author: Bastiaan Aelbrecht (template by ChatGPT)
#  Date  : 2025-10-05
############################################################
library(ggplot2)
library(dplyr)
library(readxl)
library(pracma)      # trapz, cumtrapz
library(HDInterval)  # hdi (optional)
library(MASS)        # (only for allometry uncertainty if needed)
library(scales)
library(grid)

set.seed(42)
options(stringsAsFactors = FALSE)

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

# ---------------------------
# 1) CONFIG: thresholds, hyperpriors, grids
# ---------------------------
q_circ <- 0.78           # threshold quantile for u
ci_level <- 0.95         # one-sided upper line on densities

# Hyperprior centers (analogous to gator code; adapt as desired)
# We center μ_y a bit above the largest observed y to stay proper but weakly informative.
# τ_y and τ_ξ get half-Cauchy priors; μ_ξ truncated to (-1,0) with a Normal center.
mu_y0   <- NULL  # if NULL we set from data (below) as log(max CFH) + 0.15
sd_y0   <- 0.35
mu_xi0  <- -0.25
sd_xi0  <- 0.20
hc_scale_y  <- 0.6
hc_scale_xi <- 0.5

# Core grid sizes (tune if needed)
Ny_star       <- 600     # y* points (log CFH)
Nxi           <- 220     # xi points in (-1,0)
N_mu_y        <- 161
N_tau_y       <- 41
N_mu_xi       <- 81
N_tau_xi      <- 41

# Bounds for hyper grids (fixed bands; wide but safe)
mu_xi_min  <- -0.95; mu_xi_max  <- -0.02
tau_y_min  <- 0.03;  tau_y_max  <- 1.5
tau_xi_min <- 0.03;  tau_xi_max <- 0.9

# ---------------------------
# 2) Data ingest (log CFH = log(sum of femur + humerus circumferences in mm))
# ---------------------------
df <- read_excel("Data/DEmic23_updated_Supplemental_Data_withPubYear.xlsx")
df <- df[, c("hum+fem circ (mm)", "publication year")]
colnames(df) <- c("sum_circ_mm", "year")
df <- na.omit(df)
df$log_circ <- log(as.double(df$sum_circ_mm))
y <- df$log_circ
stopifnot(length(y) >= 80)
y_max    <- max(y)

# Threshold and exceedances
u0      <- as.numeric(quantile(y, q_circ, na.rm = TRUE))
y_above <- y[y > u0]
stopifnot(length(y_above) >= 25)
ymax    <- max(y_above)

# If μ_y0 not set, choose slightly above the largest observed exceedance
if (is.null(mu_y0)) mu_y0 <- ymax + 0.15

# ---------------------------
# 3) Helpers (densities & numerics)
# ---------------------------
logsumexp <- function(v) { m <- max(v); m + log(sum(exp(v - m))) }

dhalfcauchy <- function(x, scale = 1, log = FALSE) {
  dens <- (2/(pi*scale)) * 1/(1 + (x/scale)^2)
  if (log) return(log(dens)); dens
}

# Truncated normal on (a,b)
dtnorm <- function(x, mean=0, sd=1, a=-Inf, b=Inf, log=FALSE) {
  Z <- pnorm((b-mean)/sd) - pnorm((a-mean)/sd)
  dens <- dnorm((x-mean)/sd)/sd / Z
  if (log) return(log(dens)); dens
}

# Finite-endpoint GP log-pdf for vector y (y <= y*), with sigma=(u - y*)*xi, xi<0
gp_logpdf_vec <- function(y, u, xi, ystar) {
  if (ystar <= u) return(-Inf)
  if (any(y > ystar)) return(-Inf)
  sigma <- (u - ystar) * xi
  if (sigma <= 0) return(-Inf)
  z <- 1 + (y - u)/(u - ystar)  # support check independent of xi sign aside from domain
  if (any(z <= 0)) return(-Inf)
  n <- length(y)
  - n*log(sigma) - (1/xi + 1) * sum(log(z))
}

grid_quantile <- function(x, pdf, p) {
  dx  <- c(0, diff(x))
  cdf <- cumsum(pmax(pdf,0)*dx); cdf <- cdf / max(cdf)
  as.numeric(approx(cdf, x, xout = p, ties = "ordered")$y)
}

prior_cfh_density <- function(mu_log, sd_log, cfh_grid) {
  f <- dnorm(log(cfh_grid), mean = mu_log, sd = sd_log) / cfh_grid
  Z <- trapz(cfh_grid, f)
  if (is.finite(Z) && Z > 0) f/Z else f
}

# ---------------------------
# 4) Build grids (y*, xi, hyperparams)
# ---------------------------
# y* grid (log CFH) — start just above max observation & u0
y_star_min <- max(ymax + 1e-6, u0 + 1e-6)
y_star_max <- mu_y0 + 4*sd_y0 + 0.3
y_star_grid <- seq(y_star_min, y_star_max, length.out = Ny_star)
dy <- diff(y_star_grid)[1]

# xi grid
xi_grid <- seq(-0.95, -0.02, length.out = Nxi)
dxi     <- diff(xi_grid)[1]

# Hyper grids
mu_y_grid   <- seq(mu_y0 - 4*sd_y0, mu_y0 + 4*sd_y0, length.out = N_mu_y)
tau_y_grid  <- exp(seq(log(tau_y_min), log(tau_y_max), length.out = N_tau_y))
mu_xi_grid  <- seq(mu_xi_min, mu_xi_max, length.out = N_mu_xi)
tau_xi_grid <- exp(seq(log(tau_xi_min), log(tau_xi_max), length.out = N_tau_xi))

dmu_y  <- diff(mu_y_grid)[1]
dtau_y <- diff(tau_y_grid)[1]
dmu_xi <- diff(mu_xi_grid)[1]
dtau_xi<- diff(tau_xi_grid)[1]

My <- length(mu_y_grid) * length(tau_y_grid)
Mx <- length(mu_xi_grid) * length(tau_xi_grid)

# ---------------------------
# 5) Likelihood table A = exp(log L(y* , xi)) (Ny* x Nxi)
# ---------------------------
A <- matrix(0.0, nrow = Ny_star, ncol = Nxi)
for (k in 1:Ny_star) {
  ys <- y_star_grid[k]
  logs <- vapply(xi_grid, function(xi) gp_logpdf_vec(y_above, u0, xi, ys), numeric(1))
  A[k, ] <- ifelse(is.finite(logs), exp(logs), 0.0)
}

# ---------------------------
# 6) Weight matrices W_y (Ny* x My) and W_xi (Nxi x Mx)
# ---------------------------
build_Wy <- function(ygrid, mu_grid, tau_grid, dyi) {
  Ny <- length(ygrid); My <- length(mu_grid)*length(tau_grid)
  Wy <- matrix(0.0, nrow = Ny, ncol = My)
  col <- 1L
  for (it in seq_along(tau_grid)) {
    for (im in seq_along(mu_grid)) {
      Wy[, col] <- dnorm(ygrid, mean = mu_grid[im], sd = tau_grid[it]) * dyi
      col <- col + 1L
    }
  }
  Wy
}
W_y <- build_Wy(y_star_grid, mu_y_grid, tau_y_grid, dy)

build_Wxi <- function(xi_grid, mu_grid, tau_grid, dxi) {
  Nxi <- length(xi_grid); Mx <- length(mu_grid)*length(tau_grid)
  Wxi <- matrix(0.0, nrow = Nxi, ncol = Mx)
  col <- 1L
  for (it in seq_along(tau_grid)) {
    for (im in seq_along(mu_grid)) {
      Wxi[, col] <- dtnorm(xi_grid, mean = mu_grid[im], sd = tau_grid[it], a = -1, b = 0) * dxi
      col <- col + 1L
    }
  }
  Wxi
}
W_xi <- build_Wxi(xi_grid, mu_xi_grid, tau_xi_grid, dxi)

# Precompute H = A %*% W_xi (Ny* x Mx): integrates out ξ given (μ_ξ, τ_ξ)
H <- A %*% W_xi

# ---------------------------
# 7) Hyperpriors (with Riemann weights included)
# ---------------------------
prior_y_vec <- as.vector(outer(
  dnorm(mu_y_grid, mean = mu_y0, sd = sd_y0) * dmu_y,
  dhalfcauchy(tau_y_grid, scale = hc_scale_y) * dtau_y,
  FUN = "*"
))
prior_xi_vec <- as.vector(outer(
  dtnorm(mu_xi_grid, mean = mu_xi0, sd = sd_xi0, a = -1, b = 0) * dmu_xi,
  dhalfcauchy(tau_xi_grid, scale = hc_scale_xi) * dtau_xi,
  FUN = "*"
))
eps <- 1e-300
prior_y_vec[prior_y_vec < eps]   <- eps
prior_xi_vec[prior_xi_vec < eps] <- eps

# ---------------------------
# 8) Posterior density of y* (single dataset)
#     dens[k] ∝ sum_m φ_y(y*_k|μ_y,τ_y) * prior_y[m]  ×  sum_n H[k,n] * prior_xi[n]
# ---------------------------
phi_y_all_m <- function(y_star, mu_grid, tau_grid) {
  vec <- numeric(length(mu_grid)*length(tau_grid))
  col <- 1L
  for (it in seq_along(tau_grid)) {
    for (im in seq_along(mu_grid)) {
      vec[col] <- dnorm(y_star, mean = mu_grid[im], sd = tau_grid[it])
      col <- col + 1L
    }
  }
  vec
}

posterior_y_density_single <- function(ygrid, H, prior_y_vec, prior_xi_vec, mu_grid, tau_grid) {
  Ny <- length(ygrid)
  dens <- numeric(Ny)
  # Precompute the xi-mixed likelihood slice: s_k = sum_n H[k,n] * prior_xi[n]
  s_k <- as.numeric(H %*% prior_xi_vec)   # length Ny
  # Then mix in y*-hyper with φ_y(y*_k | μ_y, τ_y) * prior_y
  for (k in 1:Ny) {
    wy <- phi_y_all_m(ygrid[k], mu_grid, tau_grid)
    dens[k] <- s_k[k] * sum(wy * prior_y_vec)
  }
  dens / trapz(ygrid, dens)
}

post_y <- posterior_y_density_single(
  ygrid = y_star_grid, H = H,
  prior_y_vec = prior_y_vec,
  prior_xi_vec = prior_xi_vec,
  mu_grid = mu_y_grid, tau_grid = tau_y_grid
)

# ---------------------------
# 9) CFH* density on mm scale + summaries
# ---------------------------
cfh_grid_mm <- exp(y_star_grid)
pdf_cfh <- post_y / cfh_grid_mm
pdf_cfh <- pdf_cfh / trapz(cfh_grid_mm, pdf_cfh)

cfh_map_mm   <- cfh_grid_mm[ which.max(pdf_cfh) ]
cfh_upper_mm <- grid_quantile(cfh_grid_mm, pdf_cfh, ci_level)

# Wide-grid prior overlay (log-normal from μ_y0, sd_y0)
t_prior_min_mm <- max(1000, exp(y_star_min) * 0.85)
t_prior_max_mm <- exp(y_star_max) * 1.25
t_prior <- seq(t_prior_min_mm, t_prior_max_mm, length.out = 4000)
prior_df <- data.frame(
  t = t_prior,
  d = dlnorm(t_prior, meanlog = mu_y0, sdlog = sd_y0)
)

p_cfh_star <- ggplot() +
  geom_area(data = prior_df, aes(t, d), fill = "grey70", alpha = 0.18) +
  geom_line(data = prior_df, aes(t, d), color = "grey50", linewidth = 0.8, alpha = 0.9) +
  geom_line(data = data.frame(t = cfh_grid_mm, d = pdf_cfh),
            aes(t, d), color = "darkblue", linewidth = 1.2) +
  geom_vline(xintercept = cfh_map_mm,  color = "purple",  linetype = "dashed",  linewidth = 1.0) +
  geom_vline(xintercept = cfh_upper_mm, color = "orange",  linetype = "dotdash", linewidth = 1.0) +
  labs(title = "CFH* posterior (single dataset, hierarchical hyperpriors)",
       x = "CFH* (mm)", y = "Density",
       caption = "Grey: prior on CFH* induced by LogNormal(μ_y0, sd_y0). Blue: posterior via quadrature.") +
  theme_science_polished

print(p_cfh_star)
ggsave(file.path(FIG_DIR, "cfh_endpoint_posterior_hier_single.png"),
       p_cfh_star, dpi = 600, width = 7, height = 5, units = "in")

cat("\n--- CFH* endpoint (mm) — hierarchical mixture (single dataset) ---\n")
cat(sprintf("MAP: %.1f mm | %d%% upper: %.1f mm\n\n",
            cfh_map_mm, round(100*ci_level), cfh_upper_mm))

# ---------------------------------------------------------
# 10) Mass propagation (centered quadratic), as mixture over y*
#      — keeps your variable names so downstream code works
# ---------------------------------------------------------
# Expect an RDS with centered quadratic (or linear) regression parameters.
# If you only have SEs + covariances, we'll assemble V accordingly.
coeff_file_cq <- "centered_quadratic_coefficients.rds"
coeffs <- readRDS(coeff_file_cq)

# Build V from provided SEs and covariances if needed
if (!("V" %in% names(coeffs)) && all(c("alpha_se","beta_se","gamma_se","cov_ab","cov_bg","cov_ag") %in% names(coeffs))) {
  coeffs$V <- matrix(c(
    coeffs$alpha_se^2, coeffs$cov_ab,       coeffs$cov_ag,
    coeffs$cov_ab,     coeffs$beta_se^2,    coeffs$cov_bg,
    coeffs$cov_ag,     coeffs$cov_bg,       coeffs$gamma_se^2
  ), nrow = 3, byrow = TRUE)
}

stopifnot(all(c("alpha","beta","gamma","mean_log_sum_circ","V","resid_sd") %in% names(coeffs)))
alpha  <- as.numeric(coeffs$alpha)
beta   <- as.numeric(coeffs$beta)
gamma  <- as.numeric(coeffs$gamma)
c0     <- as.numeric(coeffs$mean_log_sum_circ)
Vabg   <- as.matrix(coeffs$V)         # 3x3
sigma_e<- as.numeric(coeffs$resid_sd)  # predictive SD on log mass

# Trapezoidal weights on y* posterior (sum to 1)
trapz_weights <- function(x, f_un) {
  stopifnot(length(x) == length(f_un), length(x) >= 2)
  f <- pmax(f_un, 0)
  dx <- diff(x)
  mids <- (f[-1] + f[-length(f)]) * dx / 2
  total <- sum(mids)
  if (!is.finite(total) || total <= 0) stop("trapz_weights: non-positive total mass.")
  w <- numeric(length(x))
  w[1]                 <- mids[1] / 2
  w[length(x)]         <- mids[length(mids)] / 2
  if (length(x) > 2) w[2:(length(x)-1)] <- (mids[-length(mids)] + mids[-1]) / 2
  w / sum(w)
}
cdf_from_density <- function(x, f) {
  f[!is.finite(f)] <- 0
  F <- pracma::cumtrapz(x, f)
  F / max(F, na.rm = TRUE)
}

# Use the already-computed posterior on y*:
# (names kept identical to your older pipeline)
post_y <- post_y
y_star_grid <- y_star_grid

W_y <- trapz_weights(y_star_grid, post_y)

# Predictive Normal for m = log mass at each y*
z_i    <- y_star_grid - c0
mu_i   <- alpha + beta*z_i + gamma*z_i^2
X_mat  <- cbind(1, z_i, z_i^2)
quad   <- rowSums((X_mat %*% Vabg) * X_mat)
sig2_i <- pmax(quad + sigma_e^2, 1e-12)
sd_i   <- sqrt(sig2_i)

# Mixture over m
m_min <- min(mu_i - 6*sd_i)
m_max <- max(mu_i + 6*sd_i)
m_grid <- seq(m_min, m_max, length.out = 4096)

dnorm_vec <- function(x, mean, sd) (1/(sd*sqrt(2*pi))) * exp(-0.5*((x-mean)/sd)^2)

pdf_cols <- sapply(seq_along(mu_i), function(i) W_y[i] * dnorm_vec(m_grid, mu_i[i], sd_i[i]))
f_m <- rowSums(pdf_cols); f_m <- f_m / trapz(m_grid, f_m)
F_m <- cdf_from_density(m_grid, f_m)

# Transform to tons
w_t_grid <- exp(m_grid)/1e6
f_wt     <- f_m / (exp(m_grid)/1e6)
f_wt     <- f_wt / trapz(w_t_grid, f_wt)
F_wt     <- cdf_from_density(w_t_grid, f_wt)

w_map_t  <- w_t_grid[which.max(f_wt)]
w_up95_t <- as.numeric(approx(x = F_wt, y = w_t_grid, xout = ci_level, ties = "ordered", rule = 2)$y)
w_up90_t <- as.numeric(approx(x = F_wt, y = w_t_grid, xout = 0.90,     ties = "ordered", rule = 2)$y)

cat(sprintf("--- Mass endpoint (tons) via hierarchical y* mixture ---\n"))
cat(sprintf("MAP: %.1f | 90%% upper: %.1f | 95%% upper: %.1f\n\n",
            w_map_t, w_up90_t, w_up95_t))

# Plot mass density
df_mass <- data.frame(M_t = w_t_grid, dens = f_wt)
p_mass  <- ggplot(df_mass, aes(M_t, dens)) +
  geom_line(linewidth = 1.2, color = "darkblue") +
  geom_vline(xintercept = w_map_t,  color = "purple",  linetype = "dashed",  linewidth = 1) +
  geom_vline(xintercept = w_up95_t, color = "orange",  linetype = "dotdash", linewidth = 1) +
  labs(x = "Mass (tons)", y = "Density") +
  coord_cartesian(xlim = c(0, 700)) +
  theme_science_polished

# Annotate MAP vertically in same color
p_mass <- p_mass +
  annotate("text",
           x = w_map_t, y = max(f_wt, na.rm = TRUE)*0.08,
           label = sprintf("%.0f tons", w_map_t),
           size = 4.5, fontface = "bold", vjust = -0.5, angle = 90,
           color = "purple")

print(p_mass)
ggsave(file.path(FIG_DIR, "mass_endpoint_hier_single.png"),
       p_mass, dpi = 600, width = 7, height = 5, units = "in")

# Also save CFH* plot as PDF/PNG pair
ggsave(file.path(FIG_DIR, "cfh_endpoint_posterior_hier_single.pdf"),
       p_cfh_star, device = cairo_pdf, width = 7, height = 5, units = "in")
