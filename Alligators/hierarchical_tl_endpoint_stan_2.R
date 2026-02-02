############################################################
# Hierarchical EVT (two datasets) with PC prior on the gap
# δ = y* - y_max, with λ ~ Gamma(a_lambda, b_lambda)
# Author: Bastiaan Aelbrecht
# Date  : 2025-10-05
############################################################

## ---------------------------
## 0) Libraries & global setup
## ---------------------------
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readxl)
  library(pracma)
  library(grid)   # for unit() in theme
  library(rstan)
  library(MASS)   # mvrnorm for allometry uncertainty
})

options(stringsAsFactors = FALSE)
set.seed(42)

# Stan parallelization
rstan_options(auto_write = TRUE)
options(mc.cores = max(1, parallel::detectCores() - 1))

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
## 1) CONFIG: thresholds, PC prior hyperparams
## ---------------------------
# Threshold quantiles for each dataset
q_s <- 0.95
q_a <- 0.95

# Reference and CI level
stokes_tl_cm <- 450
ci_level     <- 0.95

# PC prior calibration for log-gap δ = y* - y_max:
# Choose a relative TL headroom r so U_log = log(1+r), and a prior tail prob alpha:
pc_rel_gap <- 0.1   # 2% above y_max (TL scale)
alpha_pc   <- 0.1   # P(δ > U_log) under mean λ
cv_lambda  <- 1.0    # coeff. of variation for λ hyperprior (Gamma)

U_log <- log(1 + pc_rel_gap)

pc_gamma_from_U_alpha <- function(U_log, alpha, cv = 1.0) {
  lam_bar <- -log(alpha) / max(U_log, 1e-12)  # E[λ]
  a <- 1/(cv^2)
  b <- a / lam_bar
  list(a_lambda = a, b_lambda = b)
}
pc_hyp   <- pc_gamma_from_U_alpha(U_log, alpha_pc, cv_lambda)
a_lambda <- pc_hyp$a_lambda
b_lambda <- pc_hyp$b_lambda

# Shape hierarchy for ξ
mu_xi0      <- -0.20
sd_xi0      <-  0.1
hc_scale_xi <-  0.5

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

u_vec      <- c(u_s, u_a)
ymax_vec   <- c(max(ex_s), max(ex_a))        # dataset log-maxima
y_list     <- list(ex_s, ex_a)
y_max_pool <- max(c(ex_s, ex_a))             # pooled log-maximum

## ---------------------------
## 4) STAN (PC gap prior + Gamma hyperprior)
## ---------------------------
stan_file <- "stan/hier_gp_weibull_pc_gap.stan"

## ---------------------------
## 5) Build Stan data
## ---------------------------
N1 <- length(ex_s); y1 <- as.vector(ex_s); u1 <- u_s; ymax1 <- max(ex_s)
N2 <- length(ex_a); y2 <- as.vector(ex_a); u2 <- u_a; ymax2 <- max(ex_a)

stan_data <- list(
  N1 = N1, y1 = y1, u1 = u1, ymax1 = ymax1,
  N2 = N2, y2 = y2, u2 = u2, ymax2 = ymax2,
  a_lambda = a_lambda, b_lambda = b_lambda,
  mu_xi0 = mu_xi0, sd_xi0 = sd_xi0, hc_scale_xi = hc_scale_xi,
  y_max_pool = y_max_pool
)

## ---------------------------
## 6) Sample
## ---------------------------
fit <- rstan::stan(
  file   = stan_file,
  data   = stan_data,
  seed   = 42,
  chains = 4,
  iter   = 3000,
  warmup = 1500,
  control = list(adapt_delta = 0.95, max_treedepth = 12)
)

print(fit, pars = c("delta1","delta2","lambda","mu_xi","tau_xi","xi1","xi2",
                    "ystar1","ystar2","TLstar1_cm","TLstar2_cm","TLstar_pop_cm"),
      probs = c(0.05,0.5,0.95))

post <- as.data.frame(fit)

## ---------------------------
## 7) Extract draws
## ---------------------------
ystar1_draw <- post$ystar1
ystar2_draw <- post$ystar2
xi1_draw    <- post$xi1
xi2_draw    <- post$xi2
mu_xi_draw  <- post$mu_xi
tau_xi_draw <- post$tau_xi
lambda_draw <- post$lambda
delta1_draw <- post$delta1
delta2_draw <- post$delta2

TLstar1_cm_draw    <- post$TLstar1_cm
TLstar2_cm_draw    <- post$TLstar2_cm
TLstar_pop_cm_draw <- post$TLstar_pop_cm  # predictive population TL*

## ---------------------------
## 8) Helper density & MAP
## ---------------------------
kde_df <- function(x, from=300, to=800, n=4096) {
  d <- density(x, from = from, to = to, n = n, kernel = "gaussian", adjust = 1.0)
  data.frame(t = d$x, d = d$y)
}
mode_from_kde <- function(x, from=300, to=800) {
  d <- kde_df(x, from, to)
  d$t[ which.max(d$d) ]
}

## ---------------------------
## 9) Posterior densities (Sergio/Allan/Population TL*)
## ---------------------------
dens_s <- kde_df(TLstar1_cm_draw, from = 300, to = 800)
dens_a <- kde_df(TLstar2_cm_draw, from = 300, to = 800)
dens_p <- kde_df(TLstar_pop_cm_draw, from = 300, to = 800)

sergio_map <- mode_from_kde(TLstar1_cm_draw, 300, 800)
allan_map  <- mode_from_kde(TLstar2_cm_draw, 300, 800)
pop_map    <- mode_from_kde(TLstar_pop_cm_draw, 300, 800)

sergio_up  <- as.numeric(quantile(TLstar1_cm_draw, probs = ci_level))
allan_up   <- as.numeric(quantile(TLstar2_cm_draw, probs = ci_level))
pop_up     <- as.numeric(quantile(TLstar_pop_cm_draw, probs = ci_level))

## ---------------------------
## 10) Prior overlays (PC + Gamma -> Lomax on log-gap)
## ---------------------------
prior_pc_tl <- function(t, y_max, a, b) {
  delta <- pmax(log(t) - y_max, 0)
  dens  <- (a * b^a) / ((b + delta)^(a+1)) * (1 / t)
  dens[ t < exp(y_max) ] <- 0
  dens
}

t_prior_min_cm <- 300
t_prior_max_cm <- 800
n_prior_points <- 4000
t_prior <- seq(t_prior_min_cm, t_prior_max_cm, length.out = n_prior_points)

prior_df_s <- data.frame(t = t_prior,
                         d = prior_pc_tl(t_prior, y_max = ymax1, a = a_lambda, b = b_lambda))
prior_df_a <- data.frame(t = t_prior,
                         d = prior_pc_tl(t_prior, y_max = ymax2, a = a_lambda, b = b_lambda))
prior_df_p <- data.frame(t = t_prior,
                         d = prior_pc_tl(t_prior, y_max = y_max_pool, a = a_lambda, b = b_lambda))

make_single_plot_prior_df <- function(kde_df_obj, prior_df, title_lab, map_val, up_val, map_label=NULL) {
  if (is.null(map_label)) map_label <- sprintf("%.0f cm", map_val)
  ggplot() +
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
         caption = "Grey: PC(δ) Gamma hyperprior (Lomax on log-gap). Blue: posterior (Stan).") +
    theme_science_polished + xlim(300,700)
}

p_sergio <- make_single_plot_prior_df(dens_s, prior_df_s,
                                      "TL* — Sergio dataset (Stan, PC gap prior)",
                                      sergio_map, sergio_up,
                                      sprintf("%.0f cm", sergio_map))
ggsave(file.path(fig_dir, "pcgap_stan_Sergio_TLstar.png"),
       p_sergio, dpi = 600, width = 7, height = 5, units = "in")

p_allan <- make_single_plot_prior_df(dens_a, prior_df_a,
                                     "TL* — Allan dataset (Stan, PC gap prior)",
                                     allan_map, allan_up,
                                     sprintf("%.0f cm", allan_map))
ggsave(file.path(fig_dir, "pcgap_stan_Allan_TLstar.png"),
       p_allan, dpi = 600, width = 7, height = 5, units = "in")

p_pop <- ggplot() +
  geom_area(data = prior_df_p, aes(t, d), fill = "grey70", alpha = 0.18) +
  geom_line(data = prior_df_p, aes(t, d), color = "grey50", linewidth = 0.8, alpha = 0.9) +
  geom_line(data = dens_p, aes(t, d), color = "darkblue", linewidth = 1.2) +
  geom_vline(xintercept = pop_map, color = "purple", linetype = "dashed", linewidth = 1.05) +
  geom_vline(xintercept = pop_up,  color = "orange", linetype = "dotdash", linewidth = 1.0) +
  geom_vline(xintercept = stokes_tl_cm, color = "firebrick", linetype = "dashed", linewidth = 1.0) +
  labs(x = "TL* (cm)", y = "Density", title = "Population TL* (Stan, PC gap prior)",
       caption = "Grey: PC(δ) Gamma hyperprior prior (anchored at pooled y_max). Blue: posterior (Stan).") +
  theme_science_polished + xlim(300,700)

ggsave(file.path(fig_dir, "pcgap_stan_Population_TLstar.png"),
       p_pop, dpi = 600, width = 7, height = 5, units = "in")

cat("\n--- Hierarchical TL* (Stan, PC gap prior) summaries ---\n")
cat(sprintf("Sergio:    MAP = %.1f cm | %d%% upper = %.1f cm\n",
            sergio_map, round(100*ci_level), sergio_up))
cat(sprintf("Allan:     MAP = %.1f cm | %d%% upper = %.1f cm\n",
            allan_map,  round(100*ci_level), allan_up))
cat(sprintf("Population MAP = %.1f cm | %d%% upper = %.1f cm\n\n",
            pop_map,    round(100*ci_level), pop_up))

## =========================================================
## 11) Mass propagation with population TL* (Stan PC draws)
## =========================================================
if (!exists("stokes_w_kg")) stokes_w_kg <- 459

w_fit_df_s <- dat_s %>%
  mutate(WTkg = suppressWarnings(as.numeric(Weight))) %>%
  filter(is.finite(WTkg), is.finite(TL)) %>%
  transmute(logW = log(WTkg), logTL = log(TL))
stopifnot(nrow(w_fit_df_s) >= 10)

lm_w_tl_s <- lm(logW ~ logTL, data = w_fit_df_s)
beta_s    <- coef(lm_w_tl_s)         # (a, b)
Vbeta_s   <- vcov(lm_w_tl_s)         # 2x2
sigma_e_s <- summary(lm_w_tl_s)$sigma

set.seed(42)
B  <- length(TLstar_pop_cm_draw)
idx <- sample.int(B, size = min(20000, B), replace = TRUE)

mu_y_samp <- log(TLstar_pop_cm_draw[idx])  # population endpoint (log TL*) draws
betas     <- MASS::mvrnorm(length(idx), mu = beta_s, Sigma = Vbeta_s)
eps       <- rnorm(length(idx), mean = 0, sd = sigma_e_s)

m_log <- betas[,1] + betas[,2] * mu_y_samp + eps
mass_kg_draw <- exp(m_log)

mass_map_kg  <- mode_from_kde(mass_kg_draw, 200, 800)
mass_up95_kg <- as.numeric(quantile(mass_kg_draw, probs = ci_level))

df_mass <- kde_df(mass_kg_draw, from = 200, to = 800)
p_mass_pop <- ggplot(df_mass, aes(t, d)) +
  geom_line(linewidth = 1.2, color = "darkblue") +
  geom_vline(xintercept = mass_map_kg,  color = "purple",  linetype = "dashed",  linewidth = 1.0) +
  geom_vline(xintercept = mass_up95_kg, color = "orange",  linetype = "dotdash", linewidth = 1.0) +
  geom_vline(xintercept = stokes_w_kg,  color = "firebrick", linetype = "dashed", linewidth = 1.0) +
  annotate("text", x = mass_map_kg, y = max(df_mass$d)*0.10,
           label = sprintf("%.0f kg", mass_map_kg),
           angle = 90, vjust = -0.5, size = 4.2, fontface = "bold", color = "purple") +
  labs(x = "Weight* (kg)", y = "Density",
       title = "Maximum mass — population TL* (Stan, PC prior) + Sergio allometry",
       caption = "Monte Carlo over (TL*pop, β, ε) using Stan posterior draws.") +
  coord_cartesian(xlim = c(200, 800)) +
  theme_science_polished

ggsave(file.path(fig_dir, "pcgap_mass_endpoint_population_stan.png"),
       p_mass_pop, dpi = 600, width = 7, height = 5, units = "in")

cat("\nMass endpoint (Population TL* via Stan PC + Sergio allometry):\n")
cat(sprintf("MAP: %.2f kg | %d%% upper: %.2f kg\n\n",
            mass_map_kg, round(ci_level*100), mass_up95_kg))

## =========================================================
## 12) Tail probabilities (Stan PC draws) — dataset & population
## =========================================================
# Weibull survival in log-scale:
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

# Truncated Normal RNG helper for xi predictive
rtnorm1 <- function(mu, sd, a, b) {
  pa <- pnorm(a, mu, sd); pb <- pnorm(b, mu, sd)
  u  <- runif(1, pa, pb)
  qnorm(u, mu, sd)
}

tail_eval_population_stan <- function(target_cm, q_pop = max(q_s, q_a),
                                      qlo = 0.05, qhi = 0.95) {
  y_pool <- c(y_s, y_a); N <- length(y_pool)
  u_pop  <- as.numeric(quantile(y_pool, q_pop))
  Nu     <- sum(y_pool > u_pop)
  
  # Predictive draws: y* from TLstar_pop_cm_draw; xi from hyper (mu_xi, tau_xi)
  B  <- length(TLstar_pop_cm_draw)
  id <- sample.int(B, size = min(20000, B), replace = TRUE)
  ystar_pop <- log(TLstar_pop_cm_draw[id])
  
  xi_pop <- vapply(seq_along(id), function(i) {
    m <- mu_xi_draw[id[i]]; s <- max(tau_xi_draw[id[i]], 1e-8)
    rtnorm1(m, s, -1, 0)
  }, numeric(1))
  
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

## Example usage
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
cat(sprintf("\nSergio (Stan PC):  p_cond=%.3e [%.3e, %.3e] ; p_full=%.3e\n",
            subset(res_sergio, target_cm==stokes)$p_cond,
            subset(res_sergio, target_cm==stokes)$p_cond_lo,
            subset(res_sergio, target_cm==stokes)$p_cond_hi,
            subset(res_sergio, target_cm==stokes)$p_full))
cat(sprintf("Allan  (Stan PC):  p_cond=%.3e [%.3e, %.3e] ; p_full=%.3e\n",
            subset(res_allan, target_cm==stokes)$p_cond,
            subset(res_allan, target_cm==stokes)$p_cond_lo,
            subset(res_allan, target_cm==stokes)$p_cond_hi,
            subset(res_allan, target_cm==stokes)$p_full))
cat(sprintf("Popul. (Stan PC):  p_cond=%.3e [%.3e, %.3e] ; p_full=%.3e\n\n",
            subset(res_pop, target_cm==stokes)$p_cond,
            subset(res_pop, target_cm==stokes)$p_cond_lo,
            subset(res_pop, target_cm==stokes)$p_cond_hi,
            subset(res_pop, target_cm==stokes)$p_full))
