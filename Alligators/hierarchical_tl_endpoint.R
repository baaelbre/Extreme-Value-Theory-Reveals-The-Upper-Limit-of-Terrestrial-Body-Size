############################################################
# Hierarchical EVT for log TL exceedances (two datasets)
# Author: Bastiaan Aelbrecht
# Date: 5/10/2025
############################################################

## ---------------------------
## 0) Libraries & global setup
## ---------------------------
library(ggplot2)
library(dplyr)
library(readxl)
library(pracma)
library(grid)   # for unit() in theme
options(stringsAsFactors = FALSE)
set.seed(42)

## Figure dir and plot theme
fig_dir <- "Figures"
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

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
## 1) CONFIG: thresholds, priors, grids
## ---------------------------
# Threshold quantiles for each dataset
q_s <- 0.95
q_a <- 0.95

# Stokes reference and CI level for “upper bound” line
stokes_tl_cm <- 450
ci_level     <- 0.95

# Hyperprior centers 
mu_y0  <- log(430);  sd_y0  <- 0.30        # for log TL* population center μ_y
mu_xi0 <- -0.20;     sd_xi0 <- 0.20        # for ξ population center μ_ξ
hc_scale_y  <- 0.5                        # half-Cauchy scale for τ_y
hc_scale_xi <- 0.5                        # half-Cauchy scale for τ_ξ

# Core grid sizes (adjust for speed/accuracy)
Ny_per_group  <- 320     # y* points per dataset (on log TL)
Nxi           <- 200     # xi points in (-1,0)
N_mu_y        <- 201
N_tau_y       <- 41
N_mu_xi       <- 61
N_tau_xi      <- 41

# Bounds for hyper grids
mu_y_min  <- mu_y0 - 4*sd_y0 # mu_y0 - 4*sd_y0
mu_y_max  <- mu_y0 + 4*sd_y0 # mu_y0 + 4*sd_y0
tau_y_min <- 0.01 # 0.03
tau_y_max <- 3 # 1.2
mu_xi_min  <- -0.9 # -0.7
mu_xi_max  <- -0.01 # 0.02
tau_xi_min <- 0.01 # 0.03
tau_xi_max <- 3 # 0.7

## ---------------------------
## 2) Load & prepare both datasets
## ---------------------------
# Sergio (CSV)
sergio_path <- "Data/CaptureData_Gator_allometry_paper_2024.csv"
dat_s <- read.csv(sergio_path) %>%
  mutate(
    SVL  = suppressWarnings(as.numeric(SV.length)),
    TL   = suppressWarnings(as.numeric(Total.Length)),
    Deform = case_when(
      is.na(Deformities_Notes) ~ 0L,
      nchar(trimws(Deformities_Notes)) == 0 ~ 0L,
      grepl("^none$", trimws(tolower(Deformities_Notes))) ~ 0L,
      TRUE ~ 1L
    )
  )
reg_s <- dat_s %>% filter(Deform==0L, is.finite(TL), is.finite(SVL)) %>%
  transmute(logTL = log(TL), logSVL = log(SVL))
stopifnot(nrow(reg_s) >= 10)
lm_s <- lm(logTL ~ logSVL, data = reg_s)
dat_s <- dat_s %>%
  mutate(
    logSVL    = ifelse(is.finite(SVL), log(SVL), NA_real_),
    logTL_obs = ifelse(is.finite(TL),  log(TL),  NA_real_),
    logTL_pred= ifelse(is.finite(logSVL),
                       as.numeric(predict(lm_s, newdata = data.frame(logSVL=logSVL))),
                       NA_real_),
    logTL     = ifelse(is.finite(logTL_obs), logTL_obs, logTL_pred)
  )
y_s <- dat_s$logTL[is.finite(dat_s$logTL)]

# Allan (Excel)
allan_path <- "Data/experimental_alligator_harvest_woodward.xlsx"
dat_a <- read_excel(allan_path) %>%
  mutate(
    SVL   = as.numeric(SVL),
    TL    = as.numeric(TL),
    Deform= as.integer(Deform)
  )
reg_a <- dat_a %>% filter(Deform==0L, is.finite(TL), is.finite(SVL)) %>%
  transmute(logTL = log(TL), logSVL = log(SVL))
stopifnot(nrow(reg_a) >= 10)
lm_a <- lm(logTL ~ logSVL, data = reg_a)
dat_a <- dat_a %>%
  mutate(
    logSVL    = ifelse(is.finite(SVL), log(SVL), NA_real_),
    logTL_obs = ifelse(is.finite(TL),  log(TL),  NA_real_),
    logTL_pred= ifelse(is.finite(logSVL),
                       as.numeric(predict(lm_a, newdata = data.frame(logSVL=logSVL))),
                       NA_real_),
    logTL     = ifelse(is.finite(logTL_obs), logTL_obs, logTL_pred)
  )
y_a <- dat_a$logTL[is.finite(dat_a$logTL)]

stopifnot(length(y_s) >= 50, length(y_a) >= 50)

## ---------------------------
## 3) Thresholds & exceedances
## ---------------------------
u_s <- as.numeric(quantile(y_s, q_s))
u_a <- as.numeric(quantile(y_a, q_a))
ex_s <- y_s[y_s > u_s]
ex_a <- y_a[y_a > u_a]

u_vec    <- c(u_s, u_a)
ymax_vec <- c(max(ex_s), max(ex_a))
y_list   <- list(ex_s, ex_a)

## ---------------------------
## 4) Helpers (densities & numerics)
## ---------------------------
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
  z <- 1 + (y - u)/(u - ystar)  # support check; independent of xi
  if (any(z <= 0)) return(-Inf)
  n <- length(y)
  - n*log(sigma) - (1/xi + 1) * sum(log(z))
}

grid_quantile <- function(x, pdf, p) {
  dx  <- c(0, diff(x))
  cdf <- cumsum(pmax(pdf,0)*dx); cdf <- cdf / max(cdf)
  as.numeric(approx(cdf, x, xout = p, ties = "ordered")$y)
}

prior_tl_density <- function(mu_log, sd_log, tl_grid) {
  f <- dnorm(log(tl_grid), mean = mu_log, sd = sd_log) / tl_grid
  Z <- trapz(tl_grid, f)
  if (is.finite(Z) && Z > 0) f/Z else f
}

annot_y <- function(d) max(d, na.rm = TRUE) * 0.08  # vertical label position

## ---------------------------
## 5) Build grids (y*, xi, hyperparams)
## ---------------------------
# y* per group (on log-TL scale)
make_ygrid <- function(ymax, u, mu0, sd0, n=Ny_per_group) {
  lo <- max(ymax + 1e-6, u + 1e-6)
  hi <- max(lo + 1e-3, mu0 + 4*sd0) + 0.3
  seq(lo, hi, length.out = n)
}
ygrid_list <- list(
  make_ygrid(ymax_vec[1], u_s, mu_y0, sd_y0, Ny_per_group),
  make_ygrid(ymax_vec[2], u_a, mu_y0, sd_y0, Ny_per_group)
)
dy1 <- diff(ygrid_list[[1]])[1]
dy2 <- diff(ygrid_list[[2]])[1]

# xi grid on (-1,0)
xi_grid <- seq(-0.95, -0.02, length.out = Nxi)
dxi     <- diff(xi_grid)[1]

# Hyper grids
mu_y_grid   <- seq(mu_y_min,  mu_y_max,  length.out = N_mu_y)
tau_y_grid  <- exp(seq(log(tau_y_min), log(tau_y_max), length.out = N_tau_y))
mu_xi_grid  <- seq(mu_xi_min, mu_xi_max, length.out = N_mu_xi)
tau_xi_grid <- exp(seq(log(tau_xi_min), log(tau_xi_max), length.out = N_tau_xi))

dmu_y  <- diff(mu_y_grid)[1]
dtau_y <- diff(tau_y_grid)[1]   # integrate on tau-scale directly
dmu_xi  <- diff(mu_xi_grid)[1]
dtau_xi <- diff(tau_xi_grid)[1]

My <- N_mu_y * N_tau_y
Mx <- N_mu_xi * N_tau_xi

## ---------------------------
## 6) Precompute per-group likelihood tables A_j = exp(log L_j(y*,xi))
## ---------------------------
A_list <- vector("list", 2)
for (j in 1:2) {
  yj  <- y_list[[j]]
  uj  <- u_vec[j]
  ygj <- ygrid_list[[j]]
  Ny  <- length(ygj)
  M   <- matrix(0.0, nrow = Ny, ncol = Nxi)
  for (k in 1:Ny) {
    ystar_k <- ygj[k]
    logs <- vapply(xi_grid, function(xi) gp_logpdf_vec(yj, uj, xi, ystar_k), numeric(1))
    M[k, ] <- ifelse(is.finite(logs), exp(logs), 0.0)
  }
  A_list[[j]] <- M   # Ny x Nxi
}
A1 <- A_list[[1]]
A2 <- A_list[[2]]

## ---------------------------
## 7) Weight matrices W_y (Ny x My) and W_xi (Nxi x Mx)
##     columns = densities * Δ (Riemann weights)
## ---------------------------
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
W_y1 <- build_Wy(ygrid_list[[1]], mu_y_grid, tau_y_grid, dy1)
W_y2 <- build_Wy(ygrid_list[[2]], mu_y_grid, tau_y_grid, dy2)

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


## ---------------------------
## 8) Group integrals on hyper-grids:
##     M_j = t(W_yj) %*% A_j %*% W_xi    (My x Mx)
## ---------------------------
M1 <- t(W_y1) %*% A1 %*% W_xi
M2 <- t(W_y2) %*% A2 %*% W_xi
eps <- 1e-300
M1[M1 < eps] <- eps
M2[M2 < eps] <- eps

## ---------------------------
## 9) Hyperpriors (with Riemann weights included)
## ---------------------------
prior_y_vec <- as.vector(outer(
  dnorm(mu_y_grid, mean = mu_y0, sd = sd_y0) * dmu_y,
  dhalfcauchy(tau_y_grid, scale = hc_scale_y) * dtau_y,
  FUN = "*"
))
prior_xi_vec <- as.vector(outer(
  dnorm(mu_xi_grid, mean = mu_xi0, sd = sd_xi0) * dmu_xi,
  dhalfcauchy(tau_xi_grid, scale = hc_scale_xi) * dtau_xi,
  FUN = "*"
))
prior_y_vec[prior_y_vec < eps]   <- eps
prior_xi_vec[prior_xi_vec < eps] <- eps

log_prior_y  <- log(prior_y_vec)
log_prior_xi <- log(prior_xi_vec)

## ---------------------------
## 10) Hyperposterior surface Post(m,n) ∝ prior_y[m] prior_xi[n] M1[m,n] M2[m,n]
## ---------------------------
logPost <- (matrix(log_prior_y,  nrow = My, ncol = Mx, byrow = FALSE)) +
  (matrix(log_prior_xi, nrow = My, ncol = Mx, byrow = TRUE))  +
  log(M1) + log(M2)

norm_const <- logsumexp(as.vector(logPost))
Post <- exp(logPost - norm_const)  # My x Mx, sums~1 under rectangle rule in priors

## ---------------------------
## 11) Population TL* ≈ exp(μ_y)
## ---------------------------
# Marginal over (τ_y, μ_ξ, τ_ξ)
Post_array <- array(Post, dim = c(N_mu_y, N_tau_y, Mx))
post_mu_y_mass <- apply(Post_array, 1, sum)   # length N_mu_y
post_mu_y_mass <- post_mu_y_mass / sum(post_mu_y_mass)

tl_grid_pop <- exp(mu_y_grid)
pdf_tl_pop_unnorm <- post_mu_y_mass / tl_grid_pop
pdf_tl_pop <- pdf_tl_pop_unnorm / trapz(tl_grid_pop, pdf_tl_pop_unnorm)

# Population MAP & 95% upper
pop_map <- tl_grid_pop[ which.max(pdf_tl_pop) ]
pop_up  <- grid_quantile(tl_grid_pop, pdf_tl_pop, ci_level)

## ---------------------------
## 12) Posterior of y*_j for each dataset (mix over hyperposterior and ξ)
## ---------------------------
# H_j = A_j %*% W_xi   (Ny_j x Mx)
H1 <- A1 %*% W_xi
H2 <- A2 %*% W_xi

# R_j = (prior_y ⊗ prior_xi) ∘ M_other
R1 <- (matrix(prior_y_vec,  nrow = My, ncol = Mx, byrow = FALSE)) *
  (matrix(prior_xi_vec, nrow = My, ncol = Mx, byrow = TRUE))  * M2
R2 <- (matrix(prior_y_vec,  nrow = My, ncol = Mx, byrow = FALSE)) *
  (matrix(prior_xi_vec, nrow = My, ncol = Mx, byrow = TRUE))  * M1

# Numerically stabilize (scales cancel in normalization)
scale_R1 <- max(R1); if (!is.finite(scale_R1) || scale_R1 <= 0) scale_R1 <- 1
scale_R2 <- max(R2); if (!is.finite(scale_R2) || scale_R2 <= 0) scale_R2 <- 1
R1s <- R1 / scale_R1
R2s <- R2 / scale_R2

# φ_y(y* | all (μ_y, τ_y)) as length-My vector
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

posterior_y_density <- function(ygrid, Hj, Rjs, mu_grid, tau_grid) {
  Ny <- length(ygrid)
  dens <- numeric(Ny)
  for (k in 1:Ny) {
    z_k  <- Hj[k, ]             # length Mx
    inn  <- Rjs %*% z_k         # length My
    wy   <- phi_y_all_m(ygrid[k], mu_grid, tau_grid)  # length My
    dens[k] <- sum(wy * inn)
  }
  dens / trapz(ygrid, dens)
}

pdf_y1 <- posterior_y_density(ygrid_list[[1]], H1, R1s, mu_y_grid, tau_y_grid)
pdf_y2 <- posterior_y_density(ygrid_list[[2]], H2, R2s, mu_y_grid, tau_y_grid)

# Transform to TL scale (change-of-variables)
to_tl_density <- function(ygrid, pdf_y) {
  tgrid <- exp(ygrid)
  pdf_t <- pdf_y / tgrid
  pdf_t <- pdf_t / trapz(tgrid, pdf_t)
  list(tgrid = tgrid, pdf_t = pdf_t)
}
dens_s <- to_tl_density(ygrid_list[[1]], pdf_y1)  # Sergio
dens_a <- to_tl_density(ygrid_list[[2]], pdf_y2)  # Allan

# MAP and 95% upper per dataset
sergio_map <- dens_s$tgrid[ which.max(dens_s$pdf_t) ]
sergio_up  <- grid_quantile(dens_s$tgrid, dens_s$pdf_t, ci_level)
allan_map  <- dens_a$tgrid[ which.max(dens_a$pdf_t) ]
allan_up   <- grid_quantile(dens_a$tgrid, dens_a$pdf_t, ci_level)

# =========================
# Prior-only wide grid (no change to posterior grids)
# =========================

# prior density on posterior plots
t_prior_min_cm <- 300
t_prior_max_cm <- 800
n_prior_points <- 4000      # resolution for prior only
keep_axis_to_posterior <- FALSE  # keep x-limits to the posterior’s TL grid

# Exact prior on TL* using log-normal density (no grid renormalization artifacts)
t_prior <- seq(t_prior_min_cm, t_prior_max_cm, length.out = n_prior_points)
prior_df <- data.frame(
  t = t_prior,
  d = dlnorm(t_prior, meanlog = mu_y0, sdlog = sd_y0)
)

# 3) Helper: make single plot where PRIOR uses its own (bigger) grid
make_single_plot_prior_df <- function(tgrid, post, prior_df, title_lab, map_val, up_val,
                                      map_label = NULL, clamp_axis = TRUE) {
  if (is.null(map_label)) map_label <- sprintf("%.0f cm", map_val)
  d_post <- data.frame(t = tgrid, post = post)
  
  p <- ggplot() +
    # prior on its own wide grid
    geom_area(data = prior_df, aes(t, d), fill = "grey70", alpha = 0.18) +
    geom_line(data = prior_df, aes(t, d), color = "grey50", linewidth = 0.8, alpha = 0.9) +
    # posterior on its native grid
    geom_line(data = d_post, aes(t, post), color = "darkblue", linewidth = 1.2) +
    # MAP / upper / Stokes
    geom_vline(xintercept = map_val, color = "purple",  linetype = "dashed",  linewidth = 1.05) +
    geom_vline(xintercept = up_val,  color = "orange",  linetype = "dotdash", linewidth = 1.0) +
    geom_vline(xintercept = stokes_tl_cm, color = "firebrick", linetype = "dashed", linewidth = 1.0) +
    annotate("text", x = map_val,
             y = max(post, na.rm = TRUE) * 0.08, angle = 90, vjust = -0.4,
             label = map_label, size = 4.2, fontface = "bold", color = "purple") +
    labs(x = "TL* (cm)", y = "Density", title = title_lab,
         caption = "Grey: prior on TL* (log-normal). Blue: posterior (quadrature).") +
    theme_science_polished+
    xlim(300,700)
  
  p
}

# 4) Rebuild the four figures using the wide-grid prior
# Sergio
p_sergio <- make_single_plot_prior_df(
  tgrid = dens_s$tgrid, post = dens_s$pdf_t, prior_df = prior_df,
  title_lab = "TL* — Sergio dataset",
  map_val = dens_s$tgrid[which.max(dens_s$pdf_t)],
  up_val  = { grid_quantile(dens_s$tgrid, dens_s$pdf_t, ci_level) },
  map_label = sprintf("%.0f cm",
                      dens_s$tgrid[which.max(dens_s$pdf_t)]),
  clamp_axis = keep_axis_to_posterior
)
p_sergio
ggsave(file.path(fig_dir, "hier_full_Sergio_TLstar.png"),
       p_sergio, dpi = 600, width = 7, height = 5, units = "in")

# Allan
p_allan <- make_single_plot_prior_df(
  tgrid = dens_a$tgrid, post = dens_a$pdf_t, prior_df = prior_df,
  title_lab = "TL* — Allan dataset",
  map_val = dens_a$tgrid[which.max(dens_a$pdf_t)],
  up_val  = { grid_quantile(dens_a$tgrid, dens_a$pdf_t, ci_level) },
  map_label = sprintf("%.0f cm",
                      dens_a$tgrid[which.max(dens_a$pdf_t)]),
  clamp_axis = keep_axis_to_posterior
)
p_allan
ggsave(file.path(fig_dir, "hier_full_Allan_TLstar.png"),
       p_allan, dpi = 600, width = 7, height = 5, units = "in")

# Overlay: prior stays wide; posteriors stay on their own grids (no resampling)
p_both <- ggplot() +
  geom_area(data = prior_df, aes(t, d), fill = "grey70", alpha = 0.18) +
  geom_line(data = prior_df, aes(t, d), color = "grey50", linewidth = 0.8, alpha = 0.9) +
  geom_line(data = data.frame(t = dens_s$tgrid, d = dens_s$pdf_t),
            aes(t, d, color = "Sergio"), linewidth = 1.25) +
  geom_line(data = data.frame(t = dens_a$tgrid, d = dens_a$pdf_t),
            aes(t, d, color = "Allan"), linewidth = 1.25) +
  scale_color_manual(values = c("Sergio" = "darkblue", "Allan" = "#2a9df4")) +
  geom_vline(xintercept = stokes_tl_cm, color = "firebrick", linetype = "dashed", linewidth = 1.0) +
  geom_vline(xintercept = dens_s$tgrid[which.max(dens_s$pdf_t)], color = "purple", linetype = "dashed") +
  geom_vline(xintercept = dens_a$tgrid[which.max(dens_a$pdf_t)], color = "purple", linetype = "dashed") +
  geom_vline(xintercept = grid_quantile(dens_s$tgrid, dens_s$pdf_t, ci_level), color = "orange", linetype = "dotdash") +
  geom_vline(xintercept = grid_quantile(dens_a$tgrid, dens_a$pdf_t, ci_level), color = "orange", linetype = "dotdash") +
  labs(x = "TL* (cm)", y = "Density", color = NULL,
       caption = "Grey: prior (wide grid); Blue tones: posteriors (original grids).") +
  theme_science_polished +
  theme(legend.position = "top") +xlim(300,700)

p_both
ggsave(file.path(fig_dir, "hier_full_Overlay_TLstar.png"),
       p_both, dpi = 600, width = 7.5, height = 5.2, units = "in")

# Population: prior wide, posterior unchanged
prior_pop_df <- prior_df  # same prior on TL* scale
p_pop <- ggplot() +
  geom_area(data = prior_pop_df, aes(t, d), fill = "grey70", alpha = 0.18) +
  geom_line(data = prior_pop_df, aes(t, d), color = "grey50", linewidth = 0.8, alpha = 0.9) +
  geom_line(data = data.frame(t = tl_grid_pop, d = pdf_tl_pop),
            aes(t, d), color = "darkblue", linewidth = 1.2) +
  geom_vline(xintercept = tl_grid_pop[which.max(pdf_tl_pop)], color = "purple", linetype = "dashed", linewidth = 1.05) +
  geom_vline(xintercept = grid_quantile(tl_grid_pop, pdf_tl_pop, ci_level), color = "orange", linetype = "dotdash", linewidth = 1.0) +
  geom_vline(xintercept = stokes_tl_cm, color = "firebrick", linetype = "dashed", linewidth = 1.0) +
  labs(x = "TL* (cm)", y = "Density") +
  theme_science_polished + xlim(300,700)

p_pop
## ---------------------------
## 15) Console summary
## ---------------------------
cat("\n--- Hierarchical TL* (quadrature) summaries ---\n")
cat(sprintf("Sergio:    MAP = %.1f cm | %d%% upper = %.1f cm\n",
            sergio_map, round(100*ci_level), sergio_up))
cat(sprintf("Allan:     MAP = %.1f cm | %d%% upper = %.1f cm\n",
            allan_map,  round(100*ci_level), allan_up))
cat(sprintf("Population MAP = %.1f cm | %d%% upper = %.1f cm\n\n",
            pop_map,    round(100*ci_level), pop_up))

## =========================================================
## 16) Mass propagation using POPULATION TL* posterior
##     - Allometry (a,b,Vbeta,sigma_e) from Sergio data
##     - Mixture weights from population p(TL* | all data)
## =========================================================

# Stokes reference (kg)
if (!exists("stokes_w_kg")) stokes_w_kg <- 459

# Helper: trapezoidal weights that sum to 1
trapz_weights <- function(x, f_un) {
  stopifnot(length(x) == length(f_un), length(x) >= 2)
  f <- pmax(f_un, 0)
  dx <- diff(x)
  mids <- (f[-1] + f[-length(f)]) * dx / 2
  total <- sum(mids)
  if (!is.finite(total) || total <= 0) return(rep(NA_real_, length(x)))
  w <- numeric(length(x))
  w[1]                 <- mids[1] / 2
  w[length(x)]         <- mids[length(mids)] / 2
  if (length(x) > 2) w[2:(length(x)-1)] <- (mids[-length(mids)] + mids[-1]) / 2
  w / sum(w)
}

# --- Population TL* posterior grid (TL scale) and weights ---
# Use the already-computed population posterior on TL*: (tl_grid_pop, pdf_tl_pop)
# Convert to log-TL for the allometric predictor: y*_i = log TL*_i
y_star_log_pop <- log(tl_grid_pop)
w_pop          <- trapz_weights(tl_grid_pop, pdf_tl_pop)  # mixture weights
stopifnot(all(is.finite(w_pop)), abs(sum(w_pop) - 1) < 1e-6)

# --- Allometric calibration on Sergio: logW ~ logTL ---
# (reuse Sergio table "dat_s"; weight column is "Weight" in that CSV)
w_fit_df_s <- dat_s %>%
  mutate(WTkg = suppressWarnings(as.numeric(Weight))) %>%
  filter(is.finite(WTkg), is.finite(TL)) %>%
  transmute(logW = log(WTkg), logTL = log(TL))
stopifnot(nrow(w_fit_df_s) >= 10)

lm_w_tl_s <- lm(logW ~ logTL, data = w_fit_df_s)
beta_s    <- coef(lm_w_tl_s)         # (a, b)
Vbeta_s   <- vcov(lm_w_tl_s)         # 2x2
sigma_e_s <- summary(lm_w_tl_s)$sigma

# --- Per-grid predictive Normal for m = log Weight ---
mu_i   <- beta_s[1] + beta_s[2] * y_star_log_pop
x1     <- rep(1, length(y_star_log_pop))
# predictive variance = sigma_e^2 + x' V x
sigma2_i <- sigma_e_s^2 +
  (x1 * Vbeta_s[1,1] * x1) +
  (2 * x1 * Vbeta_s[1,2] * y_star_log_pop) +
  (y_star_log_pop * Vbeta_s[2,2] * y_star_log_pop)
sd_i <- sqrt(pmax(sigma2_i, 1e-12))

# --- Mixture-of-Gaussians for m (log kg) ---
m_min <- min(mu_i - 6*sd_i)
m_max <- max(mu_i + 6*sd_i)
m_grid <- seq(m_min, m_max, length.out = 4096)

# pdf_m(m) = sum_i w_i N(m; mu_i, sd_i^2),  cdf_m(m) = sum_i w_i Phi((m-mu_i)/sd_i)
pdf_cols <- sapply(seq_along(mu_i), function(i) w_pop[i] * dnorm(m_grid, mean = mu_i[i], sd = sd_i[i]))
pdf_m <- rowSums(pdf_cols)
cdf_cols <- sapply(seq_along(mu_i), function(i) w_pop[i] * pnorm((m_grid - mu_i[i]) / sd_i[i]))
cdf_m <- rowSums(cdf_cols)
pdf_m <- pmax(pdf_m, 0)
cdf_m <- pmin(pmax(cdf_m, 0), 1)

# Transform to kg scale
mass_grid_kg <- exp(m_grid)
pdf_mass_kg  <- pdf_m / mass_grid_kg      # change of variables
cdf_mass_kg  <- cdf_m

# Summaries (MAP and 95% upper)
grid_quantile_vec <- function(x, cdf, p) {
  as.numeric(approx(x = cdf, y = x, xout = p, ties = "ordered", rule = 2)$y)
}
mass_map_kg  <- mass_grid_kg[ which.max(pdf_mass_kg) ]
mass_up95_kg <- grid_quantile_vec(mass_grid_kg, cdf_mass_kg, ci_level)

# --- Plot: density on kg with annotations ---
dens_m_y <- max(pdf_mass_kg, na.rm = TRUE)
p_mass_pop <- ggplot(data.frame(W = mass_grid_kg, d = pdf_mass_kg), aes(W, d)) +
  geom_line(color = "darkblue", linewidth = 1.2) +
  geom_vline(xintercept = mass_map_kg,  color = "purple",  linetype = "dashed",  linewidth = 1.0) +
  geom_vline(xintercept = mass_up95_kg, color = "orange",  linetype = "dotdash", linewidth = 1.0) +
  geom_vline(xintercept = stokes_w_kg,  color = "firebrick", linetype = "dashed", linewidth = 1.0) +
  annotate("text", x = mass_map_kg,  y = dens_m_y*0.10, label = sprintf("%.0f kg", mass_map_kg),
           angle = 90, vjust = -0.5, size = 4.2, fontface = "bold", color = "purple") +
  annotate("text", x = stokes_w_kg,  y = dens_m_y*0.12, label = "Stokes (459 kg)",
           angle = 90, vjust = -0.6, size = 3.6, color = "firebrick") +
  labs(x = "Weight* (kg)", y = "Density",
       title = "Maximum mass — population TL* posterior + Sergio allometry",
       caption = "Mixture-of-Gaussians over log mass; weights from population p(TL*|data).") +
  coord_cartesian(xlim = c(200, 800)) +
  theme_science_polished

print(p_mass_pop)
ggsave(file.path(fig_dir, "mass_endpoint_population_mixture.png"),
       p_mass_pop, dpi = 600, width = 7, height = 5, units = "in")

cat("\nMass endpoint (Population TL* + Sergio allometry):\n")
cat(sprintf("MAP: %.2f kg | %d%% upper: %.2f kg\n",
            mass_map_kg, round(ci_level*100), mass_up95_kg))

# --- Optional: regression panel (Sergio scatter) with POPULATION endpoint markers ---
# x-markers from population TL* posterior (μ_y = log TL*)
cdf_mu <- cumsum(post_mu_y_mass); cdf_mu <- cdf_mu / max(cdf_mu)
x_ep_map   <- mu_y_grid[ which.max(post_mu_y_mass) ]                            # log TL* (pop MAP)
x_ep_upper <- as.numeric(approx(x = cdf_mu, y = mu_y_grid, xout = ci_level)$y)   # log TL* (pop 95% up)
y_ep_map   <- log(mass_map_kg)
y_ep_upper <- log(mass_up95_kg)

plot_df_s <- dat_s %>%
  mutate(WTkg = suppressWarnings(as.numeric(Weight)),
         Sex  = toupper(trimws(Sex))) %>%
  filter(is.finite(WTkg), is.finite(TL), !is.na(Sex)) %>%
  transmute(logTL = log(TL), logW = log(WTkg), Sex = factor(Sex, levels = c("F","M")))
ab_s <- coef(lm_w_tl_s)
y_at_u_s <- as.numeric(predict(lm_w_tl_s, newdata = data.frame(logTL = u_s)))
x_stokes_log <- log(stokes_tl_cm); y_stokes_log <- log(stokes_w_kg)
sex_cols <- c("F" = "#ff6fb3", "M" = "#66b3ff")

p_reg_pop <- ggplot(plot_df_s, aes(x = logTL, y = logW, color = Sex)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_abline(intercept = ab_s[1], slope = ab_s[2], linewidth = 1) +
  geom_vline(xintercept = u_s, linetype = "dashed", color = "gray30") +
  geom_hline(yintercept = y_at_u_s, linetype = "dashed", color = "gray30") +
  geom_point(aes(x = x_ep_map, y = y_ep_map), color = "darkgreen", size = 3, inherit.aes = FALSE) +
  geom_segment(aes(x = x_ep_map, y = y_ep_map, xend = x_ep_upper, yend = y_ep_map),
               color = "darkgreen", linewidth = 1.1, inherit.aes = FALSE) +
  geom_segment(aes(x = x_ep_map, y = y_ep_map, xend = x_ep_map, yend = y_ep_upper),
               color = "darkgreen", linewidth = 1.1, inherit.aes = FALSE) +
  geom_point(aes(x = x_stokes_log, y = y_stokes_log),
             inherit.aes = FALSE, color = "firebrick", size = 2.8) +
  geom_text(aes(x = x_stokes_log, y = y_stokes_log, label = "Stokes"),
            inherit.aes = FALSE, nudge_y = 0.015, color = "firebrick", size = 3.5) +
  scale_color_manual(values = sex_cols, breaks = c("F","M"), na.translate = FALSE, drop = FALSE) +
  labs(x = "log Total Length", y = "log Weight", color = "Sex") +
  theme_science_polished

print(p_reg_pop)
ggsave(file.path(fig_dir, "regression_panel_population_markers.png"),
       p_reg_pop, dpi = 600, width = 8, height = 5.2, units = "in")

## ---------------------------------------------------------
## 16b) Male-only regression panel (like old code)
##      - points: Sergio males only
##      - line: male-only OLS (logW ~ logTL)
##      - crosshair: male threshold at q=0.95 of male logTL
##      - endpoint markers: POPULATION TL* (x) + POPULATION mass (y)
## ---------------------------------------------------------

# Male subset from Sergio, with weight
male_df_s <- dat_s %>%
  mutate(
    WTkg = suppressWarnings(as.numeric(Weight)),
    Sex  = toupper(trimws(Sex))
  ) %>%
  filter(Sex == "M", is.finite(WTkg), is.finite(TL)) %>%
  transmute(logTL = log(TL), logW = log(WTkg))

stopifnot(nrow(male_df_s) >= 10)

# Male-only OLS for display
lm_w_tl_m <- lm(logW ~ logTL, data = male_df_s)
ab_m      <- coef(lm_w_tl_m)

# Male threshold (q=0.95 on male logTL, like old crosshairs)
u0_M <- as.numeric(quantile(male_df_s$logTL, 0.95, na.rm = TRUE))
y_at_u0_M <- as.numeric(predict(lm_w_tl_m, newdata = data.frame(logTL = u0_M)))

# Population TL* markers on x (from μ_y grid)
cdf_mu <- cumsum(post_mu_y_mass); cdf_mu <- cdf_mu / max(cdf_mu)
x_ep_map   <- mu_y_grid[ which.max(post_mu_y_mass) ]                            # log TL* (pop MAP)
x_ep_upper <- as.numeric(approx(x = cdf_mu, y = mu_y_grid, xout = ci_level)$y)   # log TL* (pop 95% up)

# Population mass markers on y (from Section 16 mixture over TL* population)
y_ep_map   <- log(mass_map_kg)
y_ep_upper <- log(mass_up95_kg)

# Stokes in log space
x_stokes_log <- log(stokes_tl_cm)
y_stokes_log <- log(stokes_w_kg)

# Aesthetics
male_col     <- "#66b3ff"
cross_col    <- "gray30"
endpoint_col <- "darkgreen"

p_male_only <- ggplot(male_df_s, aes(x = logTL, y = logW)) +
  # points + male-only regression
  geom_point(color = male_col, alpha = 0.75, size = 2) +
  geom_abline(intercept = ab_m[1], slope = ab_m[2], linewidth = 1) +
  # threshold crosshair
  geom_vline(xintercept = u0_M, linetype = "dashed", color = cross_col) +
  geom_hline(yintercept = y_at_u0_M, linetype = "dashed", color = cross_col) +
  # endpoint (population TL* on x; population mass on y)
  geom_point(aes(x = x_ep_map, y = y_ep_map),
             inherit.aes = FALSE, color = endpoint_col, size = 3) +
  geom_segment(aes(x = x_ep_map, y = y_ep_map, xend = x_ep_upper, yend = y_ep_map),
               inherit.aes = FALSE, color = endpoint_col, linewidth = 1.1) +
  geom_segment(aes(x = x_ep_map, y = y_ep_map, xend = x_ep_map, yend = y_ep_upper),
               inherit.aes = FALSE, color = endpoint_col, linewidth = 1.1) +
  # Stokes
  geom_point(aes(x = x_stokes_log, y = y_stokes_log),
             inherit.aes = FALSE, color = "firebrick", size = 2.8) +
  geom_text(aes(x = x_stokes_log, y = y_stokes_log, label = "Stokes"),
            inherit.aes = FALSE, nudge_y = 0.015, color = "firebrick", size = 3.5) +
  labs(x = "log(Total Length [cm])", y = "log(Mass [kg])") +
  theme_science_polished

print(p_male_only)
ggsave(file.path(fig_dir, "regression_panel_male_only_hierarchical.png"),
       p_male_only, dpi = 600, width = 8, height = 5.2, units = "in")

## =========================================================
## Deterministic tail probabilities
## =========================================================

# Ensure density (not area-weighted) bases are available
dens_y1 <- W_y1 / dy1                 # Ny1 x My
dens_y2 <- W_y2 / dy2                 # Ny2 x My
dens_x  <- W_xi / dxi                 # Nxi x Mx

# Weibull survival on the (y*, xi) grid for a given target
# Vectorized Weibull survival over the (y*, xi) grid
S_weibull_mat <- function(y_log, u_log, ystar_grid, xi_grid) {
  Ny  <- length(ystar_grid)
  Nxi <- length(xi_grid)
  
  YS  <- matrix(ystar_grid, nrow = Ny, ncol = Nxi, byrow = FALSE)  # y* by column
  XI  <- matrix(xi_grid,   nrow = Ny, ncol = Nxi, byrow = TRUE)    # xi by row
  
  # Ratio = ((y* - y) / (y* - u)); valid only when y < y*
  ratio <- (YS - y_log) / (YS - u_log)
  
  # Start with zeros; fill only where y < y* and ratio>0
  S <- matrix(0.0, nrow = Ny, ncol = Nxi)
  ok <- (YS > y_log) & is.finite(ratio) & (ratio > 0)
  
  # Weibull-domain survival S_W = ratio^(-1/xi)
  S[ok] <- ratio[ok]^(-1 / XI[ok])
  
  # Clean any numerical junk
  S[!is.finite(S)] <- 0.0
  S
}


# Weighted quantiles on a discrete grid
weighted_quantile <- function(values, weights, probs = c(0.05, 0.5, 0.95)) {
  o   <- order(values)
  v   <- values[o]; w <- pmax(weights[o], 0)
  W   <- sum(w); if (!is.finite(W) || W <= 0) return(rep(NA_real_, length(probs)))
  cw  <- cumsum(w) / W
  as.numeric(approx(x = cw, y = v, xout = probs, ties = "ordered", rule = 2)$y)
}

# -------- Dataset-specific (Sergio / Allan) ----------
# Build joint posterior over (y*, xi) ∝ [Σ_{m,n} R_js[m,n] φ_y_m(y*) φ_x_n(ξ)] × A_j(y*, ξ)
# Then normalize over grid cells (include Δy*Δξ), and integrate S_W against it.

tail_eval_dataset_quad <- function(target_cm,
                                   A_j, dens_y_j, Rjs,        # Ny×Nxi, Ny×My, My×Mx
                                   ygrid, xi_grid, u_log,     # grids & dataset threshold
                                   dy, dxi, Nu, N,            # cell areas & tail fraction pieces
                                   qlo = 0.05, qhi = 0.95) {
  # Mixture kernel over (y*, xi) from hyperposterior (no sampling)
  K_mid  <- dens_y_j %*% Rjs          # Ny x Mx
  K      <- K_mid %*% t(dens_x)       # Ny x Nxi
  Joint  <- pmax(K * A_j, 0)          # Ny x Nxi (unnormalized)
  
  # Normalize to a discrete joint over the grid
  Joint  <- Joint * dy * dxi
  Z      <- sum(Joint)
  if (!is.finite(Z) || Z <= 0) stop("Joint normalization failed (dataset case).")
  W_joint <- Joint / Z
  
  eval_one <- function(t_cm) {
    ylog <- log(t_cm)
    S    <- S_weibull_mat(ylog, u_log, ygrid, xi_grid)
    p_c  <- sum(W_joint * S)                        # conditional tail probability
    # credible band from the posterior distribution of S over the grid
    qs   <- weighted_quantile(as.vector(S), as.vector(W_joint), c(qlo, qhi))
    data.frame(
      target_cm = t_cm,
      p_cond    = pmin(pmax(p_c, 0), 1),
      p_cond_lo = pmin(pmax(qs[1], 0), 1),
      p_cond_hi = pmin(pmax(qs[2], 0), 1),
      p_full    = (Nu/N) * p_c,
      p_full_lo = (Nu/N) * qs[1],
      p_full_hi = (Nu/N) * qs[2]
    )
  }
  do.call(rbind, lapply(as.numeric(target_cm), eval_one))
}

tail_eval_sergio_quad <- function(target_cm) {
  Nu <- length(ex_s); N <- length(y_s)
  tail_eval_dataset_quad(
    target_cm,
    A_j = A1, dens_y_j = dens_y1, Rjs = R1s,
    ygrid = ygrid_list[[1]], xi_grid = xi_grid, u_log = u_s,
    dy = dy1, dxi = dxi, Nu = Nu, N = N
  )
}

tail_eval_allan_quad <- function(target_cm) {
  Nu <- length(ex_a); N <- length(y_a)
  tail_eval_dataset_quad(
    target_cm,
    A_j = A2, dens_y_j = dens_y2, Rjs = R2s,
    ygrid = ygrid_list[[2]], xi_grid = xi_grid, u_log = u_a,
    dy = dy2, dxi = dxi, Nu = Nu, N = N
  )
}

# -------- Population (hyperposterior predictive) ----------
# For a new draw from the population: (y*, ξ) ~ Σ_{m,n} Post[m,n] φ_y_m ⊗ φ_x_n.
# No dataset likelihood multiplier here (predictive at the population level).
tail_eval_population_quad <- function(target_cm, q_pop = max(q_s, q_a),
                                      Ny_pop = 320, qlo = 0.05, qhi = 0.95) {
  y_pool <- c(y_s, y_a); N <- length(y_pool)
  u_pop  <- as.numeric(quantile(y_pool, q_pop))
  Nu     <- sum(y_pool > u_pop)
  
  # Population y* grid + densities
  ymax_pool <- max(c(ex_s, ex_a))
  ygrid_pop <- make_ygrid(ymax = ymax_pool, u = u_pop, mu0 = mu_y0, sd0 = sd_y0, n = Ny_pop)
  dy_pop    <- diff(ygrid_pop)[1]
  dens_y_pop <- build_Wy(ygrid_pop, mu_y_grid, tau_y_grid, dyi = dy_pop) / dy_pop  # Ny_pop x My
  
  # Mixture kernel over (y*, ξ): dens_y_pop %*% Post %*% t(dens_x)
  K_mid  <- dens_y_pop %*% Post         # Ny_pop x Mx
  K      <- K_mid %*% t(dens_x)         # Ny_pop x Nxi
  Joint  <- pmax(K, 0)
  Joint  <- Joint * dy_pop * dxi
  Z      <- sum(Joint)
  if (!is.finite(Z) || Z <= 0) stop("Joint normalization failed (population case).")
  W_joint <- Joint / Z
  
  eval_one <- function(t_cm) {
    ylog <- log(t_cm)
    S    <- S_weibull_mat(ylog, u_pop, ygrid_pop, xi_grid)
    p_c  <- sum(W_joint * S)
    qs   <- weighted_quantile(as.vector(S), as.vector(W_joint), c(qlo, qhi))
    data.frame(
      target_cm = t_cm,
      p_cond    = pmin(pmax(p_c, 0), 1),
      p_cond_lo = pmin(pmax(qs[1], 0), 1),
      p_cond_hi = pmin(pmax(qs[2], 0), 1),
      p_full    = (Nu/N) * p_c,
      p_full_lo = (Nu/N) * qs[1],
      p_full_hi = (Nu/N) * qs[2]
    )
  }
  do.call(rbind, lapply(as.numeric(target_cm), eval_one))
}

## -------------------------
## Example usage
## -------------------------
targets <- c(400, 425, 450, 475)

res_sergio_q <- tail_eval_sergio_quad(targets)
res_allan_q  <- tail_eval_allan_quad(targets)
res_pop_q    <- tail_eval_population_quad(targets)

print(res_sergio_q)
print(res_allan_q)
print(res_pop_q)

# Quick lines for the Stokes length (450 cm)
stokes <- 450
cat(sprintf("\nSergio (quad):  p_cond=%.3e [%.3e, %.3e] ; p_full=%.3e\n",
            subset(res_sergio_q, target_cm==stokes)$p_cond,
            subset(res_sergio_q, target_cm==stokes)$p_cond_lo,
            subset(res_sergio_q, target_cm==stokes)$p_cond_hi,
            subset(res_sergio_q, target_cm==stokes)$p_full))
cat(sprintf("Allan  (quad):  p_cond=%.3e [%.3e, %.3e] ; p_full=%.3e\n",
            subset(res_allan_q, target_cm==stokes)$p_cond,
            subset(res_allan_q, target_cm==stokes)$p_cond_lo,
            subset(res_allan_q, target_cm==stokes)$p_cond_hi,
            subset(res_allan_q, target_cm==stokes)$p_full))
cat(sprintf("Popul. (quad):  p_cond=%.3e [%.3e, %.3e] ; p_full=%.3e\n\n",
            subset(res_pop_q, target_cm==stokes)$p_cond,
            subset(res_pop_q, target_cm==stokes)$p_cond_lo,
            subset(res_pop_q, target_cm==stokes)$p_cond_hi,
            subset(res_pop_q, target_cm==stokes)$p_full))

