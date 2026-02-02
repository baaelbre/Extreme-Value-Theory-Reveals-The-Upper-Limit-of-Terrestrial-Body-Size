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
q_s <- 0.97
q_a <- 0.95

# Stokes reference and CI level for “upper bound” line
stokes_tl_cm <- 450
ci_level     <- 0.95

# Hyperprior centers 
mu_y0  <- log(420);  sd_y0  <- 0.20        # for log TL* population center μ_y
mu_xi0 <- -0.25;     sd_xi0 <- 0.10        # for ξ population center μ_ξ
hc_scale_y  <- 0.1                        # half-Cauchy scale for τ_y best: 0.1
hc_scale_xi <- 0.1                        # half-Cauchy scale for τ_ξ best: 0.1

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
## 4) Stan setup (compile once)
## ---------------------------
suppressPackageStartupMessages({
  library(rstan)
  library(MASS)      # mvrnorm for allometry uncertainty
})

rstan_options(auto_write = TRUE)
options(mc.cores = max(1, parallel::detectCores() - 1))

stan_file <- "stan/hier_gp_weibull.stan"


## ---------------------------
## 5) Build Stan data
## ---------------------------
N1 <- length(ex_s); y1 <- as.vector(ex_s); u1 <- u_s; ymax1 <- max(ex_s)
N2 <- length(ex_a); y2 <- as.vector(ex_a); u2 <- u_a; ymax2 <- max(ex_a)

stan_data <- list(
  N1 = N1, y1 = y1, u1 = u1, ymax1 = ymax1,
  N2 = N2, y2 = y2, u2 = u2, ymax2 = ymax2,
  mu_y0 = mu_y0,    sd_y0 = sd_y0,
  mu_xi0 = mu_xi0,  sd_xi0 = sd_xi0,
  hc_scale_y = hc_scale_y,
  hc_scale_xi = hc_scale_xi
)

## ---------------------------
## 6) Sample
## ---------------------------
fit <- rstan::stan(
  file   = stan_file,
  data   = stan_data,
  seed   = 42,
  chains = 4,
  iter   = 10000,   # adjust if you want more draws
  warmup = 2000,
  control = list(adapt_delta = 0.98, max_treedepth = 10)
)

#fit <- rstan::stan(
#file   = stan_file,
#data   = stan_data,
#seed   = 42,
#chains = 4,
#iter   = 10000,   # adjust if you want more draws
#warmup = 2000,
#control = list(adapt_delta = 0.95, max_treedepth = 12)
#) works extremely well./

print(fit, pars = c("ystar1","ystar2","xi1","xi2","mu_y","tau_y","mu_xi","tau_xi",
                    "TLstar1_cm","TLstar2_cm","TLstar_pop_cm"),
      probs = c(0.05,0.5,0.95))

## =========================================================
## 6b) Diagnostics: R-hat, ESS, trace plots, NUTS checks
## =========================================================
suppressPackageStartupMessages({
  library(bayesplot)
  library(posterior)
})
color_scheme_set("brightblue")

# Parameters to diagnose (add/remove as needed)
params_main <- c("ystar1","ystar2","xi1","xi2","mu_y","tau_y","mu_xi","tau_xi")
params_derived <- c("TLstar1_cm","TLstar2_cm","TLstar_pop_cm")

# ---------- Summaries: R-hat & ESS ----------
drws <- as_draws_array(fit)                   # posterior::as_draws_array
summ_tbl <- summarise_draws(
  drws,
  mean, sd, ~quantile(.x, probs = c(0.05,0.5,0.95)),
  rhat, ess_bulk, ess_tail
)
diag_tbl <- subset(summ_tbl, variable %in% c(params_main, params_derived))
print(diag_tbl, n = nrow(diag_tbl))

# Save as CSV for the supplement
diag_csv_path <- file.path(fig_dir, "stan_diagnostics_rhat_ess.csv")
write.csv(as.data.frame(diag_tbl), diag_csv_path, row.names = FALSE)
cat(sprintf("Wrote diagnostics table to: %s\n", diag_csv_path))

# ---------- Trace plots ----------
mcmc_arr <- rstan::As.mcmc.list(fit)          # for bayesplot convenience
posterior_array <- rstan::extract(fit, permuted = FALSE)  #  iterations x chains x params

# Trace: core hierarchical parameters
p_trace_core <- mcmc_trace(posterior_array, pars = params_main, facet_args = list(ncol = 2))
ggsave(file.path(fig_dir, "trace_core_params.png"), p_trace_core, dpi = 300, width = 12, height = 8, units = "in")

# Trace: derived endpoints (optional)
p_trace_der <- mcmc_trace(posterior_array, pars = params_derived, facet_args = list(ncol = 2))
ggsave(file.path(fig_dir, "trace_derived_params.png"), p_trace_der, dpi = 300, width = 12, height = 6, units = "in")

# ---------- Rank histograms (uniform ranks indicate good mixing) ----------
p_rank <- mcmc_rank_hist(drws, pars = params_main)
p_rank
ggsave(file.path(fig_dir, "rank_hist_core.png"), p_rank, dpi = 300, width = 10, height = 6, units = "in")

# ---------- Autocorrelation (short-lag) ----------
p_acf <- mcmc_acf(posterior_array, pars = params_main, lags = 30)
p_acf
ggsave(file.path(fig_dir, "acf_core_params.png"), p_acf, dpi = 300, width = 12, height = 8, units = "in")

# ---------- Density overlays by chain ----------
p_dens <- mcmc_dens_overlay(posterior_array, pars = params_main)
p_dens
ggsave(file.path(fig_dir, "dens_overlay_core.png"), p_dens, dpi = 300, width = 12, height = 8, units = "in")

# ---------- NUTS diagnostics: divergences, treedepth, energy ----------
samp <- rstan::get_sampler_params(fit, inc_warmup = FALSE)

n_div <- sum(vapply(samp, function(x) sum(x[, "divergent__"]), numeric(1)))
n_treedepth_sat <- sum(vapply(samp, function(x) sum(x[, "treedepth__"] >= 12), numeric(1)))  # match your max_treedepth
bfmi <- vapply(samp, function(x) {
  # Bayesian Fraction of Missing Information (energy diagnostic)
  numer <- sum(diff(x[, "energy__"])^2)
  denom <- sum((x[, "energy__"] - mean(x[, "energy__"]))^2)
  if (denom > 0) numer / denom else NA_real_
}, numeric(1))
cat(sprintf("Divergences: %d\n", n_div))
cat(sprintf("Transitions hitting max treedepth: %d\n", n_treedepth_sat))
cat(sprintf("BFMI (per chain): %s ; min BFMI = %.3f\n", paste(sprintf("%.3f", bfmi), collapse = ", "), min(bfmi, na.rm = TRUE)))

# Plots for NUTS diagnostics
np <- nuts_params(fit)
p_div <- mcmc_parcoord(
  as.matrix(fit, pars = params_main),
  np = np,
  pars = params_main
) + ggtitle("Parallel coordinates (red = divergent transitions)")
p_div
ggsave(file.path(fig_dir, "nuts_parcoord_divergences.png"), p_div, dpi = 300, width = 10, height = 6, units = "in")

# Energy diagnostics (should look roughly Gaussian; low BFMI is problematic)
p_energy <- mcmc_nuts_energy(np)
p_energy
ggsave(file.path(fig_dir, "nuts_energy.png"), p_energy, dpi = 300, width = 8, height = 5, units = "in")

# ---------- Quick sanity checks & guidance ----------
bad_rhat <- subset(diag_tbl, rhat >= 1.01)
low_ess  <- subset(diag_tbl, ess_bulk < 1000 | ess_tail < 1000)

if (nrow(bad_rhat) > 0) {
  message("Parameters with R-hat >= 1.01:\n",
          paste(bad_rhat$variable, collapse = ", "),
          "\nConsider: non-centered parameterization, higher adapt_delta (e.g., 0.98), more warmup, or stronger priors.")
}
if (nrow(low_ess) > 0) {
  message("Parameters with low ESS (< 1000):\n",
          paste(low_ess$variable, collapse = ", "),
          "\nConsider: reparameterization, longer chains, or thinning AFTER improving geometry (not before).")
}
if (n_div > 0) {
  message(sprintf("Found %d divergent transition(s). Consider raising adapt_delta to 0.98–0.99, non-centering, or tightening priors.", n_div))
}
if (n_treedepth_sat > 0) {
  message(sprintf("Found %d transitions hitting max treedepth. Consider increasing max_treedepth (e.g., 13) and/or improving geometry.", n_treedepth_sat))
}
if (min(bfmi, na.rm = TRUE) < 0.2) {
  message("Low BFMI detected (< 0.2). Energy exploration may be poor; consider reparameterization or stronger priors.")
}


post <- as.data.frame(fit)

## Extract draws
ystar1_draw <- post$ystar1      # log TL*
ystar2_draw <- post$ystar2
xi1_draw    <- post$xi1
xi2_draw    <- post$xi2
mu_y_draw   <- post$mu_y
tau_y_draw  <- post$tau_y
mu_xi_draw  <- post$mu_xi
tau_xi_draw <- post$tau_xi

TLstar1_cm_draw <- exp(ystar1_draw)
TLstar2_cm_draw <- exp(ystar2_draw)
TLstar_pop_cm_draw <- exp(mu_y_draw)

## ---------------------------
## 7) Helper density & MAP
## ---------------------------
kde_df <- function(x, from=300, to=800, n=5012) {
  d <- density(x, from = from, to = to, n = n, kernel = "gaussian", adjust = 1.0)
  data.frame(t = d$x, d = d$y)
}
mode_from_kde <- function(x, from=300, to=800) {
  d <- kde_df(x, from, to)
  d$t[ which.max(d$d) ]
}

## ---------------------------
## 8) Posterior densities (Sergio/Allan/Population TL*)
## ---------------------------
dens_s <- kde_df(TLstar1_cm_draw, from = 300, to = 800)
dens_a <- kde_df(TLstar2_cm_draw, from = 300, to = 800)
dens_p <- kde_df(TLstar_pop_cm_draw, from = 300, to = 800)

sergio_map <- mode_from_kde(TLstar1_cm_draw, 300, 800)
allan_map  <- mode_from_kde(TLstar2_cm_draw, 300, 800)
pop_map    <- mode_from_kde(TLstar_pop_cm_draw, 300, 800)

# Means & medians on the TL (cm) scale
sergio_mean   <- mean(TLstar1_cm_draw)
sergio_median <- median(TLstar1_cm_draw)

allan_mean    <- mean(TLstar2_cm_draw)
allan_median  <- median(TLstar2_cm_draw)

pop_mean      <- mean(TLstar_pop_cm_draw)
pop_median    <- median(TLstar_pop_cm_draw)


sergio_up  <- as.numeric(quantile(TLstar1_cm_draw, probs = ci_level))
allan_up   <- as.numeric(quantile(TLstar2_cm_draw, probs = ci_level))
pop_up     <- as.numeric(quantile(TLstar_pop_cm_draw, probs = ci_level))

## Prior overlay for TL*pop (as before)
t_prior_min_cm <- 300
t_prior_max_cm <- 800
n_prior_points <- 4000
t_prior <- seq(t_prior_min_cm, t_prior_max_cm, length.out = n_prior_points)
prior_df <- data.frame(
  t = t_prior,
  d = dlnorm(t_prior, meanlog = mu_y0, sdlog = sd_y0)
)

make_single_plot_prior_df <- function(kde_df_obj, prior_df, title_lab, map_val, up_val, map_label=NULL) {
  if (is.null(map_label)) map_label <- sprintf("%.0f cm", map_val)
  p <- ggplot() +
    geom_area(data = prior_df, aes(t, d), fill = "grey70", alpha = 0.18) +
    geom_line(data = prior_df, aes(t, d), color = "grey50", linewidth = 0.8, alpha = 0.9) +
    geom_line(data = kde_df_obj, aes(t, d), color = "darkblue", linewidth = 1.2) +
    geom_vline(xintercept = map_val, color = "purple",  linetype = "dashed",  linewidth = 1.05) +
    geom_vline(xintercept = up_val,  color = "orange",  linetype = "dotdash", linewidth = 1.0) +
    geom_vline(xintercept = stokes_tl_cm, color = "firebrick", linetype = "dashed", linewidth = 1.0) +
    annotate("text", x = map_val,
             y = max(kde_df_obj$d, na.rm = TRUE) * 0.08, angle = 90, vjust = -0.4,
             label = map_label, size = 4.2, fontface = "bold", color = "purple") +
    labs(x = "TL* (cm)", y = "Density", title = title_lab,
         caption = "Grey: prior on TL* (log-normal). Blue: posterior (Stan).") +
    theme_science_polished + xlim(300,700)
  p
}

p_sergio <- make_single_plot_prior_df(dens_s, prior_df,
                                      "TL* — Sergio dataset (Stan)", sergio_map, sergio_up,
                                      sprintf("%.0f cm", sergio_map))
p_sergio
ggsave(file.path(fig_dir, "hier_stan_Sergio_TLstar.png"),
       p_sergio, dpi = 600, width = 7, height = 5, units = "in")

p_allan <- make_single_plot_prior_df(dens_a, prior_df,
                                     "TL* — Allan dataset (Stan)", allan_map, allan_up,
                                     sprintf("%.0f cm", allan_map))
p_allan
ggsave(file.path(fig_dir, "hier_stan_Allan_TLstar.png"),
       p_allan, dpi = 600, width = 7, height = 5, units = "in")

p_pop <- ggplot() +
  geom_line(data = dens_p, aes(t, d), color = "darkblue", linewidth = 1.2) +
  geom_vline(xintercept = pop_map, color = "purple", linetype = "dashed", linewidth = 1.05) +
  annotate("text", x = pop_map,
           y = 0.002, angle = 90, vjust = -0.4,
           label = "452 cm", size = 4.2, fontface = "bold", color = "purple")+
  geom_vline(xintercept = pop_up,  color = "orange", linetype = "dotdash", linewidth = 1.0) +
  geom_vline(xintercept = stokes_tl_cm, color = "firebrick", linetype = "dashed", linewidth = 1.0) +
  labs(x = "TL* (cm)", y = "Density") +
  theme_science_polished + xlim(300,700)
p_pop

post_both_df <- dplyr::bind_rows(
  dplyr::mutate(dens_s, dataset = "Sergio"),
  dplyr::mutate(dens_a, dataset = "Allan")
)
p_both <- ggplot() +
  geom_area(data = prior_df, aes(t, d), fill = "grey75", alpha = 0.22, show.legend = FALSE) +
  geom_line(data = prior_df, aes(t, d), color = "grey50", linewidth = 0.9, show.legend = FALSE) +
  geom_line(data = post_both_df, aes(t, d, color = dataset), linewidth = 1.25, show.legend = FALSE) +
  scale_color_manual(values = c("Sergio"="#1f77b4","Allan"="#2ca02c"), guide = "none") +
  coord_cartesian(xlim = c(300, 700)) +
  labs(x="TL* (cm)", y="Density") +
  theme_science_polished +
  theme(legend.position = "none")

p_both
ggsave(file.path(fig_dir, "sergio_allan_tl_endpoint.png"),
       p_both, dpi = 600, width = 7.0, height = 5.0, units = "in")

ggsave(file.path(fig_dir, "pop_tl_endpoint.png"),
       p_pop, dpi = 600, width = 7, height = 5, units = "in")

cat("\n--- Hierarchical TL* (Stan) summaries ---\n")
cat(sprintf("Sergio:     MAP = %.1f cm | median = %.1f cm | mean = %.1f cm | %d%% upper = %.1f cm\n",
            sergio_map, sergio_median, sergio_mean, round(100*ci_level), sergio_up))
cat(sprintf("Allan:      MAP = %.1f cm | median = %.1f cm | mean = %.1f cm | %d%% upper = %.1f cm\n",
            allan_map,  allan_median,  allan_mean,  round(100*ci_level), allan_up))
cat(sprintf("Population: MAP = %.1f cm | median = %.1f cm | mean = %.1f cm | %d%% upper = %.1f cm\n\n",
            pop_map,    pop_median,    pop_mean,    round(100*ci_level), pop_up))

## =========================================================
## 9) Mass propagation from TL* (Sergio, Allan, Population)
## =========================================================
# Goal: map endpoint draws on TL (cm) to mass (kg) via log-linear allometry.
# We fit pooled CSV allometry (Sergio dataset), then push the uncertainty in
# (a, b, epsilon) through the TL* draws from Stan.

# Reference (for vertical line in plots)
if (!exists("stokes_w_kg")) stokes_w_kg <- 459

# 9.1 Fit pooled CSV allometry: log(M) = a + b * log(TL)
w_fit_df_s <- dat_s %>%
  mutate(WTkg = suppressWarnings(as.numeric(Weight))) %>%
  filter(is.finite(WTkg), is.finite(TL)) %>%
  transmute(logW = log(WTkg), logTL = log(TL))
stopifnot(nrow(w_fit_df_s) >= 10)

lm_w_tl_s <- lm(logW ~ logTL, data = w_fit_df_s)
beta_s    <- coef(lm_w_tl_s)         # (a, b)
Vbeta_s   <- vcov(lm_w_tl_s)         # 2x2 covariance for (a, b)
sigma_e_s <- summary(lm_w_tl_s)$sigma # residual SD on log scale

# 9.2 Helper to propagate TL* draws to mass draws (kg)
#     For each TL* draw, sample (a, b) ~ N(beta_s, Vbeta_s) and epsilon ~ N(0, sigma_e)
mass_from_tlstar <- function(tlstar_cm_draw, beta_mu, beta_V, sigma_e,
                             ndraw = 20000, seed = 42) {
  set.seed(seed)
  idx <- sample.int(length(tlstar_cm_draw), size = min(ndraw, length(tlstar_cm_draw)), replace = TRUE)
  tl_log <- log(tlstar_cm_draw[idx])
  betas  <- MASS::mvrnorm(n = length(idx), mu = beta_mu, Sigma = beta_V)  # draws of (a, b)
  eps    <- rnorm(length(idx), mean = 0, sd = sigma_e)                    # residuals
  logM   <- betas[,1] + betas[,2] * tl_log + eps
  exp(logM)
}

# 9.3 Build mass posteriors for each target (Sergio, Allan, Population)
#     - Sergio & Allan use dataset-specific TL* draws (ystar1/ystar2 on TL scale)
#     - Population uses mu_y (log TL*) draws, which represent the population TL* parameter
mass_sergio <- mass_from_tlstar(TLstar1_cm_draw, beta_s, Vbeta_s, sigma_e_s)
mass_allan  <- mass_from_tlstar(TLstar2_cm_draw, beta_s, Vbeta_s, sigma_e_s)

set.seed(42)
nd   <- min(20000, length(mu_y_draw))
idp  <- sample.int(length(mu_y_draw), size = nd, replace = TRUE)
betp <- MASS::mvrnorm(n = nd, mu = beta_s, Sigma = Vbeta_s)
epsp <- rnorm(nd, 0, sigma_e_s)
logM_pop <- betp[,1] + betp[,2] * mu_y_draw[idp] + epsp   # NOTE: mu_y is already log(TL*)
mass_pop <- exp(logM_pop)

# 9.4 Summaries: MAP (via KDE mode), mean, median, and 95% upper
#     We'll re-use the earlier mode_from_kde() but define a mass-specific KDE wrapper here.
kde_df_mass   <- function(x, from = 100, to = 800, n = 4096) {
  d <- density(x, from = from, to = to, n = n, kernel = "gaussian", adjust = 1.0)
  data.frame(t = d$x, d = d$y)
}
mass_mode_kg  <- function(x, from = 100, to = 800) mode_from_kde(x, from, to)

mass_mode <- function(x, from = 50, to = 1200) mode_from_kde(x, from, to)
summ_mass <- function(x, name, ci = 0.95) {
  data.frame(
    which   = name,
    MAP     = mass_mode(x),
    mean    = mean(x),
    median  = median(x),
    up95    = as.numeric(quantile(x, probs = ci))
  )
}

sum_sergio <- summ_mass(mass_sergio, "Sergio")
sum_allan  <- summ_mass(mass_allan,  "Allan")
sum_pop    <- summ_mass(mass_pop,    "Population")

print(rbind(sum_sergio, sum_allan, sum_pop), row.names = FALSE)

cat("\n--- Mass endpoints (predictive; TL* mapped through pooled CSV allometry) ---\n")
cat(sprintf("Sergio:     MAP = %.1f kg | median = %.1f kg | mean = %.1f kg | %d%% upper = %.1f kg\n",
            sum_sergio$MAP, sum_sergio$median, sum_sergio$mean, round(100*ci_level), sum_sergio$up95))
cat(sprintf("Allan:      MAP = %.1f kg | median = %.1f kg | mean = %.1f kg | %d%% upper = %.1f kg\n",
            sum_allan$MAP,  sum_allan$median,  sum_allan$mean,  round(100*ci_level), sum_allan$up95))
cat(sprintf("Population: MAP = %.1f kg | median = %.1f kg | mean = %.1f kg | %d%% upper = %.1f kg\n\n",
            sum_pop$MAP,    sum_pop$median,    sum_pop$mean,    round(100*ci_level), sum_pop$up95))

## =========================================================
## 9b) Mass plots (Sergio, Allan, Population) + combined
##      — with 95% upper bounds, no legends, Stokes line only
## =========================================================

# Densities for plotting
kde_df_mass   <- function(x, from = 100, to = 800, n = 4096) {
  d <- density(x, from = from, to = to, n = n, kernel = "gaussian", adjust = 1.0)
  data.frame(t = d$x, d = d$y)
}
dens_m_sergio <- kde_df_mass(mass_sergio, from = 100, to = 800)
dens_m_allan  <- kde_df_mass(mass_allan,  from = 100, to = 800)
dens_m_pop    <- kde_df_mass(mass_pop,    from = 100, to = 800)

# MAPs (for vertical lines)
mass_mode_kg  <- function(x, from = 100, to = 800) mode_from_kde(x, from, to)
mass_map_sergio <- mass_mode_kg(mass_sergio, 100, 800)
mass_map_allan  <- mass_mode_kg(mass_allan,  100, 800)
mass_map_pop    <- mass_mode_kg(mass_pop,    100, 800)

# 95% uppers
mass_up95_sergio <- as.numeric(quantile(mass_sergio, probs = ci_level))
mass_up95_allan  <- as.numeric(quantile(mass_allan,  probs = ci_level))
mass_up95_pop    <- as.numeric(quantile(mass_pop,    probs = ci_level))

# Colors
col_sergio <- "#1f77b4"
col_allan  <- "#2ca02c"
col_pop    <- "darkblue"

# -------- Separate plots --------
p_mass_sergio <- ggplot(dens_m_sergio, aes(t, d)) +
  geom_line(linewidth = 1.2, color = col_sergio) +
  geom_vline(xintercept = mass_map_sergio,  color = "purple",   linetype = "dashed",  linewidth = 1.0) +
  geom_vline(xintercept = mass_up95_sergio, color = "orange",   linetype = "dotdash", linewidth = 1.0) +
  geom_vline(xintercept = stokes_w_kg,      color = "firebrick",linetype = "dashed",  linewidth = 1.0) +
  annotate("text", x = mass_map_sergio, y = max(dens_m_sergio$d)*0.10,
           label = sprintf("%.0f kg", mass_map_sergio),
           angle = 90, vjust = -0.5, size = 4.0, fontface = "bold", color = "purple") +
  labs(x = "Weight* (kg)", y = "Density", title = "Mass endpoint — Sergio (predictive at TL*)") +
  coord_cartesian(xlim = c(100, 800)) +
  theme_science_polished +
  theme(legend.position = "none")
p_mass_sergio
ggsave(file.path(fig_dir, "mass_endpoint_sergio.png"), p_mass_sergio, dpi = 600, width = 7, height = 5)

p_mass_allan <- ggplot(dens_m_allan, aes(t, d)) +
  geom_line(linewidth = 1.2, color = col_allan) +
  geom_vline(xintercept = mass_map_allan,   color = "purple",   linetype = "dashed",  linewidth = 1.0) +
  geom_vline(xintercept = mass_up95_allan,  color = "orange",   linetype = "dotdash", linewidth = 1.0) +
  geom_vline(xintercept = stokes_w_kg,      color = "firebrick",linetype = "dashed",  linewidth = 1.0) +
  annotate("text", x = mass_map_allan, y = max(dens_m_allan$d)*0.10,
           label = sprintf("%.0f kg", mass_map_allan),
           angle = 90, vjust = -0.5, size = 4.0, fontface = "bold", color = "purple") +
  labs(x = "Weight* (kg)", y = "Density", title = "Mass endpoint — Allan (predictive at TL*)") +
  coord_cartesian(xlim = c(100, 800)) +
  theme_science_polished +
  theme(legend.position = "none")
p_mass_allan
ggsave(file.path(fig_dir, "mass_endpoint_allan.png"), p_mass_allan, dpi = 600, width = 7, height = 5)

p_mass_pop <- ggplot(dens_m_pop, aes(t, d)) +
  geom_line(linewidth = 1.2, color = col_pop) +
  geom_vline(xintercept = mass_map_pop,     color = "purple",   linetype = "dashed",  linewidth = 1.0) +
  geom_vline(xintercept = mass_up95_pop,    color = "orange",   linetype = "dotdash", linewidth = 1.0) +
  geom_vline(xintercept = stokes_w_kg,      color = "firebrick",linetype = "dashed",  linewidth = 1.0) +
  annotate("text", x = mass_map_pop, y = max(dens_m_pop$d)*0.10,
           label = sprintf("%.0f kg", mass_map_pop),
           angle = 90, vjust = -0.5, size = 4.0, fontface = "bold", color = "purple") +
  labs(x = "Weight* (kg)", y = "Density") +
  coord_cartesian(xlim = c(100, 800)) +
  theme_science_polished +
  theme(legend.position = "none")
p_mass_pop
ggsave(file.path(fig_dir, "mass_endpoint_population.png"), p_mass_pop, dpi = 600, width = 7, height = 5)

# -------- Combined overlay: Sergio vs Allan --------
both_mass_df <- bind_rows(
  mutate(dens_m_sergio, set = "Sergio"),
  mutate(dens_m_allan,  set = "Allan")
)

p_mass_both <- ggplot(both_mass_df, aes(t, d, color = set)) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = c("Sergio" = col_sergio, "Allan" = col_allan)) +
  labs(x = "Weight* (kg)", y = "Density", color = NULL) +
  coord_cartesian(xlim = c(100, 800)) +
  theme_science_polished +
  theme(legend.position = "none")
p_mass_both
ggsave(file.path(fig_dir, "mass_endpoint_sergio_allan.png"),
       p_mass_both, dpi = 600, width = 8, height = 5.5)

## =========================================================
## Log–log regression plot with thresholds & endpoints
## =========================================================

## =========================================================
## Log–log regression plot (MAP point + orthogonal lines to uppers)
## =========================================================

a_hat <- as.numeric(beta_s[1]); b_hat <- as.numeric(beta_s[2])

# Threshold on log TL (Sergio) and mapped log W (deterministic)
logTL_thresh <- u_s
logW_thresh  <- a_hat + b_hat * logTL_thresh

# Endpoint (Sergio dataset): MAP & 95% upper on TL, mapped deterministically to W
logTL_map   <- log(sergio_map)
logTL_up95  <- log(sergio_up)
logW_map    <- a_hat + b_hat * logTL_map
logW_up95   <- a_hat + b_hat * logTL_up95

# Stokes alligator (brown dot)
stokes_logTL <- log(stokes_tl_cm)
stokes_logW  <- log(stokes_w_kg)

# Plot window
rng_x <- range(w_fit_df_s$logTL, na.rm = TRUE)
rng_y <- range(w_fit_df_s$logW,  na.rm = TRUE)
pad_x <- diff(rng_x) * 0.05; pad_y <- diff(rng_y) * 0.05
x_min <- min(rng_x[1], logTL_thresh, logTL_map, logTL_up95) - pad_x
x_max <- max(rng_x[2], logTL_thresh, logTL_map, logTL_up95) + pad_x
y_min <- min(rng_y[1], logW_thresh,  logW_map,  logW_up95,  stokes_logW) - pad_y
y_max <- max(rng_y[2], logW_thresh,  logW_map,  logW_up95,  stokes_logW) + pad_y

p_loglog <- ggplot(w_fit_df_s, aes(x = logTL, y = logW)) +
  # data used in W~TL regression
  geom_point(alpha = 0.55, size = 1.8, color = "grey25") +
  # fitted regression line
  geom_abline(intercept = a_hat, slope = b_hat, linewidth = 1.0, color = "steelblue") +
  
  # thresholds (red dashed)
  geom_vline(xintercept = logTL_thresh, color = "red3", linetype = "dashed", linewidth = 0.9) +
  geom_hline(yintercept = logW_thresh,  color = "red3", linetype = "dashed", linewidth = 0.9) +
  
  # MAP as a purple point
  geom_point(aes(x = logTL_map, y = logW_map), inherit.aes = FALSE,
             color = "purple", size = 3.2) +
  
  # Orthogonal segments from MAP to the 95% uppers (orange dot-dash)
  geom_segment(aes(x = logTL_map, y = logW_map,
                   xend = logTL_up95, yend = logW_map),
               inherit.aes = FALSE, color = "orange", linetype = "dotdash", linewidth = 1.0) +
  geom_segment(aes(x = logTL_map, y = logW_map,
                   xend = logTL_map, yend = logW_up95),
               inherit.aes = FALSE, color = "orange", linetype = "dotdash", linewidth = 1.0) +
  
  
  # Stokes alligator (brown)
  geom_point(aes(x = stokes_logTL, y = stokes_logW),
             inherit.aes = FALSE, color = "#8B4513", size = 3.2) +

  coord_cartesian(xlim = c(x_min, x_max), ylim = c(y_min, y_max)) +
  labs(
    x = "log TL (cm)", y = "log Weight (kg)",
  ) +
  theme_science_polished +
  theme(legend.position = "none")

p_loglog
ggsave(file.path(fig_dir, "mass_loglog_regression_plot.png"),
       p_loglog, dpi = 600, width = 7.5, height = 5.5, units = "in")

ggsave(file.path(fig_dir, "mass_loglog_regression_plot.png"),
       p_loglog, dpi = 600, width = 7.5, height = 5.5, units = "in")


## =========================================================
## 10) Tail probabilities (Stan draws) — dataset & population
## =========================================================

# Weighted survival in Weibull domain, using log-scale formula:
# S_W(y | u, y*, xi) = ((y* - y)/(y* - u))^(-1/xi) for y < y*, else 0.
S_weibull_draws <- function(target_cm, u_log, ystar_draw, xi_draw) {
  y_log <- log(target_cm)
  ratio <- (ystar_draw - y_log) / (ystar_draw - u_log)
  out   <- ifelse(ystar_draw > y_log & is.finite(ratio) & ratio > 0,
                  ratio ^ (-1/xi_draw), 0.0)
  pmin(pmax(out, 0), 1)
}

tail_eval_dataset_stan <- function(target_cm, u_log, ystar_draw, xi_draw, Nu, N,
                                   qlo = 0.05, qhi = 0.95) {
  SW <- S_weibull_draws(target_cm, u_log, ystar_draw, xi_draw)
  data.frame(
    target_cm = target_cm,
    p_cond    = mean(SW),
    p_cond_lo = quantile(SW, probs = qlo),
    p_cond_hi = quantile(SW, probs = qhi),
    p_full    = (Nu/N) * mean(SW),
    p_full_lo = (Nu/N) * quantile(SW, probs = qlo),
    p_full_hi = (Nu/N) * quantile(SW, probs = qhi)
  )
}

# Population predictive: draw (y*, xi) at population from (mu_y, tau_y, mu_xi, tau_xi)
rtnorm1 <- function(mu, sd, a, b) {
  u <- runif(1, pnorm(a, mu, sd), pnorm(b, mu, sd))
  qnorm(u, mu, sd)
}
tail_eval_population_stan <- function(target_cm, q_pop = max(q_s, q_a),
                                      qlo = 0.05, qhi = 0.95) {
  y_pool <- c(y_s, y_a); N <- length(y_pool)
  u_pop  <- as.numeric(quantile(y_pool, q_pop))
  Nu     <- sum(y_pool > u_pop)
  
  # sample predictive (y*, xi) from hyper draws
  B <- length(mu_y_draw)
  idx <- sample.int(B, size = min(20000, B), replace = TRUE)
  ystar_pop <- rnorm(length(idx), mean = mu_y_draw[idx], sd = pmax(tau_y_draw[idx], 1e-8))
  # xi is truncated Normal on (-1,0)
  xi_pop <- vapply(seq_along(idx),
                   function(i) rtnorm1(mu_xi_draw[idx[i]], pmax(tau_xi_draw[idx[i]], 1e-8), -1, 0),
                   numeric(1))
  SW <- S_weibull_draws(target_cm, u_log = u_pop, ystar_draw = ystar_pop, xi_draw = xi_pop)
  
  data.frame(
    target_cm = target_cm,
    p_cond    = mean(SW),
    p_cond_lo = quantile(SW, probs = qlo),
    p_cond_hi = quantile(SW, probs = qhi),
    p_full    = (Nu/N) * mean(SW),
    p_full_lo = (Nu/N) * quantile(SW, probs = qlo),
    p_full_hi = (Nu/N) * quantile(SW, probs = qhi)
  )
}

## Example usage (same targets as before)
targets <- c(400, 425, 450, 475)
res_sergio <- do.call(rbind, lapply(targets, function(t)
  tail_eval_dataset_stan(t, u_log = u_s, ystar_draw = ystar1_draw, xi_draw = xi1_draw,
                         Nu = length(ex_s), N = length(y_s))))
res_allan  <- do.call(rbind, lapply(targets, function(t)
  tail_eval_dataset_stan(t, u_log = u_a, ystar_draw = ystar2_draw, xi_draw = xi2_draw,
                         Nu = length(ex_a), N = length(y_a))))
res_pop    <- do.call(rbind, lapply(targets, function(t)
  tail_eval_population_stan(t)))

print(res_sergio)
print(res_allan)
print(res_pop)

stokes <- 450
cat(sprintf("\nSergio (Stan):  p_cond=%.3e [%.3e, %.3e] ; p_full=%.3e\n",
            subset(res_sergio, target_cm==stokes)$p_cond,
            subset(res_sergio, target_cm==stokes)$p_cond_lo,
            subset(res_sergio, target_cm==stokes)$p_cond_hi,
            subset(res_sergio, target_cm==stokes)$p_full))
cat(sprintf("Allan  (Stan):  p_cond=%.3e [%.3e, %.3e] ; p_full=%.3e\n",
            subset(res_allan, target_cm==stokes)$p_cond,
            subset(res_allan, target_cm==stokes)$p_cond_lo,
            subset(res_allan, target_cm==stokes)$p_cond_hi,
            subset(res_allan, target_cm==stokes)$p_full))
cat(sprintf("Popul. (Stan):  p_cond=%.3e [%.3e, %.3e] ; p_full=%.3e\n\n",
            subset(res_pop, target_cm==stokes)$p_cond,
            subset(res_pop, target_cm==stokes)$p_cond_lo,
            subset(res_pop, target_cm==stokes)$p_cond_hi,
            subset(res_pop, target_cm==stokes)$p_full))

