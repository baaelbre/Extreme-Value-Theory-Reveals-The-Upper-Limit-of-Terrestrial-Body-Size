#### =========================
#### EVT on log Total Length (FEMALES ONLY)
#### with point imputation from log SVL → log TL
#### Updated for CaptureData_Gator_allometry_paper_2024.csv
#### Woodward prior + LOO cross-validation (fast options)
#### =========================

## ---------------------------
## 0) Libraries & global setup
## ---------------------------
library(ggplot2)
library(extRemes)    # fevd/revd/pevd for GP fits & sims
library(dplyr)
library(readxl)      # in case you load other files elsewhere
library(pracma)      # trapz for numerical integration
library(HDInterval)  # hdi: highest posterior density intervals
library(MASS)        # mvrnorm for allometry uncertainty (weight mapping)

set.seed(42)

fig_dir <- "Figures"
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

## ---------------------------
## 1) Settings & plotting theme
## ---------------------------
ci_level      <- 0.95
n_boot_ks     <- 1000        # parametric bootstrap size for PIT→Uniform KS p-value
n_post_ystar  <- 50000       # posterior draws for TL endpoint y*
thresh_q_opt  <- 0.935        # threshold u0 = q_0.95(log TL)

# Reference lines (Stokes alligator) for posterior plots
stokes_tl_cm  <- 450   # TL = 450 cm
stokes_w_kg   <- 459   # weight = 459 kg

# Lake Apopka (female benchmark)
apopka_tl_cm  <- 322
apopka_w_kg   <- 170

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
## 2) Data ingest (filter FEMALES)
## ---------------------------
data_path <- "Data/CaptureData_Gator_allometry_paper_2024.csv"
dat_raw <- read.csv(data_path, stringsAsFactors = FALSE)

dat <- dat_raw %>%
  mutate(
    SVL   = suppressWarnings(as.numeric(SV.length)),
    TL    = suppressWarnings(as.numeric(Total.Length)),
    WTkg  = suppressWarnings(as.numeric(Weight)),
    Sex   = toupper(trimws(Sex)),
    Deform = case_when(
      is.na(Deformities_Notes) ~ 0L,
      nchar(trimws(Deformities_Notes)) == 0 ~ 0L,
      grepl("^none$", trimws(tolower(Deformities_Notes))) ~ 0L,
      TRUE ~ 1L
    )
  ) %>%
  # FEMALES ONLY from here on
  filter(Sex == "F") %>%
  mutate(
    SVL  = ifelse(is.finite(SVL)  & SVL  > 0, SVL,  NA_real_),
    TL   = ifelse(is.finite(TL)   & TL   > 0, TL,   NA_real_),
    WTkg = ifelse(is.finite(WTkg) & WTkg > 0, WTkg, NA_real_)
  )

stopifnot(nrow(dat) >= 50)  # sanity: enough females

## ---------------------------
## 3) Helper utilities
## ---------------------------
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

mrl_data <- function(y, u_seq) {
  sapply(u_seq, function(u){
    exc <- y[y > u] - u
    if (length(exc) < 5) return(NA_real_)
    mean(exc)
  })
}

make_thresholds <- function(y, q_from = 0.8, q_to = 0.999, n = 50){
  rng <- quantile(y, c(q_from, q_to), na.rm = TRUE)
  seq(rng[1], rng[2], length.out = n)
}

adj_scale <- function(scale, shape, u, u0) scale - shape * (u - u0)

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
    labs(x = "Theoretical quantiles", y = "Empirical quantiles") +
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
    labs(x = "Theoretical CDF", y = "Empirical CDF") +
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

ks_boot_pvalue <- function(y, u, scale_hat, shape_hat, B = 1000, seed = 42, refit = FALSE) {
  set.seed(seed)
  exc    <- y[y > u] - u
  n      <- length(exc)
  F_obs  <- pevd(exc, scale = scale_hat, shape = shape_hat, type = "GP")
  ks_obs <- suppressWarnings(ks.test(F_obs, "punif")$statistic)
  
  ks_b <- rep(NA_real_, B)
  for (b in seq_len(B)) {
    sim_exc <- revd(n, scale = scale_hat, shape = shape_hat, type = "GP")
    if (refit) {
      fb <- try(fevd(c(y[y <= u], u + sim_exc), type = "GP", threshold = u, time.units = "none"),
                silent = TRUE)
      if (inherits(fb, "try-error")) next
      par_b <- fb$results$par
      F_b   <- pevd(sim_exc, scale = par_b["scale"], shape = par_b["shape"], type = "GP")
    } else {
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
## 4) Impute missing TL via FEMALE SVL→TL regression (log–log OLS)
## ---------------------------------------------------------
reg_df <- dat %>%
  filter(Deform == 0L, is.finite(TL), is.finite(SVL)) %>%
  transmute(logTL = log(TL), logSVL = log(SVL))
stopifnot(nrow(reg_df) >= 10)

lm_tl_svl <- lm(logTL ~ logSVL, data = reg_df)

dat <- dat %>%
  mutate(
    logSVL    = ifelse(is.finite(SVL), log(SVL), NA_real_),
    logTL_obs = ifelse(is.finite(TL),  log(TL),  NA_real_),
    logTL_pred= ifelse(is.finite(logSVL),
                       as.numeric(predict(lm_tl_svl, newdata = data.frame(logSVL = logSVL))),
                       NA_real_),
    logTL_use = ifelse(is.finite(logTL_obs), logTL_obs, logTL_pred)
  )

# Tail variable for POT (FEMALES)
y <- dat$logTL_use
y <- y[is.finite(y)]
stopifnot(length(y) >= 50)

## ---------------------------------------------------------
## 5) Threshold selection (MRL; shape & adjusted-scale stability)
## ---------------------------------------------------------
u_seq <- make_thresholds(y, q_from = 0.7, q_to = 0.97, n = 50)
u0    <- as.numeric(quantile(y, thresh_q_opt, na.rm = TRUE))
fit0  <- fit_gpd_at_u(y, u0)
shape_hat0 <- fit0$shape
adj_hat0   <- fit0$scale

shape_se0 <- sqrt(fit0$cov[2,2])
z         <- shape_hat0 / shape_se0
p_wald    <- pnorm(z)  # one-sided for H1: xi < 0
p_wald

# MRL
mrl_vals <- mrl_data(y, u_seq)
p_mrl <- ggplot(data.frame(u = u_seq, mrl = mrl_vals), aes(u, mrl)) +
  geom_line() + geom_point() +
  geom_vline(xintercept = u0, linetype = "dashed", color = "red") +
  labs(
       x = "Threshold (log TL)", y = "Mean excess") +
  theme_science_polished
ggsave(file.path(fig_dir, "f_tl_mrl_sergio.png"), p_mrl, dpi = 600, width = 7, height = 5, units = "in")

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
  labs(
       x = "Threshold (log TL)", y = "Shape parameter") +
  theme_science_polished +
  geom_hline(yintercept = shape_hat0, color = "red", linetype = "dashed")
p_shape
ggsave(file.path(fig_dir, "f_shape_stability_sergio.png"), p_shape, dpi = 600, width = 7, height = 5, units = "in")

p_adj <- ggplot(data.frame(u = u_seq, adj, adj_lo, adj_hi), aes(u, adj)) +
  geom_point() +
  geom_errorbar(aes(ymin = adj_lo, ymax = adj_hi), width = 0.03, color = "blue") +
  geom_vline(xintercept = u0, color = "red", linetype = "dashed") +
  labs(
       x = "Threshold (log TL)", y = "Adjusted scale parameter") +
  theme_science_polished +
  geom_hline(yintercept = adj_hat0, color = "red", linetype = "dashed")
p_adj
ggsave(file.path(fig_dir, "f_adj_scale_stability_sergio.png"), p_adj, dpi = 600, width = 7, height = 5, units = "in")

## ---------------------------------------------------------
## 6) Tail diagnostics
## ---------------------------------------------------------
diags <- diagnostic_plots(y, u0, fit0$scale, fit0$shape, "log TL")
ksout <- ks_boot_pvalue(y, u0, fit0$scale, fit0$shape, B = n_boot_ks)

ggsave(file.path(fig_dir, "f_qq_sergio.png"),  diags$pqq, dpi = 600, width = 7, height = 5, units = "in")
ggsave(file.path(fig_dir, "f_qq_sergio.png"),  diags$ppp, dpi = 600, width = 7, height = 5, units = "in")
ggsave(file.path(fig_dir, "f_pit_sergio.png"), diags$pks, dpi = 600, width = 7, height = 5, units = "in")

message(sprintf("[FEMALES log TL] u0=%.3f | scale=%.4f shape=%.4f | KS p=%.3f",
                u0, fit0$scale, fit0$shape, ksout$p_value))

## ---------------------------------------------------------
## 7) Endpoint posterior (Woodward prior; grid-only integration)
## ---------------------------------------------------------
y_above <- y[y > u0]; n_ex <- length(y_above)
if (n_ex < 10) stop("Not enough female exceedances for endpoint posterior.")

Ymax <- max(y_above)
eps  <- max(1e-6, 1e-4 * abs(Ymax))
y_star_grid <- seq(Ymax + eps, Ymax + 10, length.out = 2000)  # grid for y* (log TL*)
xi_grid     <- seq(-1.0, -0.02, length.out = 1000)            # xi < 0 (Weibull)

loglik_endpoint <- function(y_star, xi, u, y_ex) {
  if (xi >= 0 || any(y_star <= y_ex)) return(-Inf)
  n <- length(y_ex)
  term1 <- n * log(y_star - u) / xi
  term2 <- -(1/xi + 1) * sum(log(y_star - y_ex))
  term3 <- -n * log(-xi)
  term1 + term2 + term3
}

# Woodward prior for TL* centered at 330 cm on the log scale (species-wide reference)
mu_y <- log(330)
prior_mode <- "woodward_loose"  # loose|mod|tight
sd_y <- switch(prior_mode,
               "woodward_loose" = 0.095,
               "woodward_mod"   = 0.048,
               "woodward_tight" = 0.024, 0.048)
logprior_y  <- function(y_star) dnorm(y_star, mean = mu_y, sd = sd_y, log = TRUE)

# PC prior on kappa = -xi in (0,1): Exp(lambda_k), truncated (0,1)
lambda_k   <- 3.0
logprior_xi <- function(xi) {
  if (xi >= 0 || xi <= -1) return(-Inf)
  kappa <- -xi
  if (kappa <= 0 || kappa >= 1) return(-Inf)
  log(lambda_k) - lambda_k * kappa
}

marg_log_post <- sapply(y_star_grid, function(ys) {
  lv <- sapply(xi_grid, function(xi) loglik_endpoint(ys, xi, u0, y_above) + logprior_xi(xi))
  lv_max <- max(lv)
  log_int_xi <- log(pracma::trapz(xi_grid, exp(lv - lv_max))) + lv_max
  log_int_xi + logprior_y(ys)
})

post_unnorm <- exp(marg_log_post - max(marg_log_post))
Z_y         <- pracma::trapz(y_star_grid, post_unnorm)
if (!is.finite(Z_y) || Z_y <= 0) stop("Normalization failed for p(y*|y).")
post <- post_unnorm / Z_y

cdf_from_density <- function(x, f) {
  f[!is.finite(f)] <- 0
  F <- pracma::cumtrapz(x, f)
  rng <- max(F, na.rm = TRUE)
  if (!is.finite(rng) || rng <= 0) stop("cdf_from_density: non-positive integral.")
  F / rng
}
q_from_cdf <- function(x, f, p) {
  F <- cdf_from_density(x, f)
  as.numeric(approx(x = F, y = x, xout = p, ties = "ordered", rule = 2)$y)
}

ystar_map   <- y_star_grid[which.max(post)]
ystar_upper <- q_from_cdf(y_star_grid, post, ci_level)

tl_grid_cm   <- exp(y_star_grid)
dens_tl      <- post / tl_grid_cm
tl_map_cm    <- exp(ystar_map)
tl_upper_cm  <- exp(ystar_upper)

df_tlpost <- data.frame(TL = tl_grid_cm, dens = dens_tl)
p_tlstar <- ggplot(df_tlpost, aes(TL, dens)) +
  geom_line(linewidth = 1.2, color = "darkblue") +
  geom_vline(xintercept = tl_map_cm,    color = "purple",   linetype = "dashed",  linewidth = 1.1) +
  geom_vline(xintercept = tl_upper_cm,  color = "orange",   linetype = "dotdash", linewidth = 1.0) +
  geom_vline(xintercept = stokes_tl_cm, color = "firebrick", linetype = "dashed", linewidth = 1.0) +
  annotate("text", x = stokes_tl_cm, y = max(dens_tl)*0.08, vjust = -0.8,
           label = "Stokes (450 cm)", color = "firebrick", angle = 90, size = 4) +
  labs(x = "TL* (cm)", y = "Density") +
  theme_science_polished
ggsave(file.path(fig_dir, "f_tl_endpoint_posterior_sergio.png"),
       p_tlstar, dpi = 600, width = 7, height = 5, units = "in")

cat("\n--- FEMALE TL endpoint (cm) via numeric integration ---\n",
    "MAP:          ", round(tl_map_cm, 2), "\n",
    sprintf("%d%% upper:  ", round(ci_level*100)), round(tl_upper_cm, 2), "\n")

## ---------------------------------------------------------
## 8) Mass extrapolation (female a,b; female y*; one-sided predictive uplift)
## ---------------------------------------------------------
ystar_draws <- sample(y_star_grid, size = 5*n_post_ystar, replace = TRUE, prob = post)

w_fit_df <- dat %>%
  filter(Deform == 0L, is.finite(WTkg), is.finite(TL)) %>%
  transmute(logW = log(WTkg), logTL = log(TL))
if (nrow(w_fit_df) < 10) {
  w_fit_df <- dat %>% filter(is.finite(WTkg), is.finite(TL)) %>%
    transmute(logW = log(WTkg), logTL = log(TL))
}
stopifnot(nrow(w_fit_df) >= 10)

lm_w_tl <- lm(logW ~ logTL, data = w_fit_df)  # FEMALE regression
beta    <- coef(lm_w_tl)
Vcoef   <- vcov(lm_w_tl)
sigma_e <- summary(lm_w_tl)$sigma

ab_samps <- MASS::mvrnorm(n = length(ystar_draws), mu = beta, Sigma = Vcoef)
a_s <- ab_samps[, 1]
b_s <- ab_samps[, 2]

logM_draws <- a_s + b_s * ystar_draws
q_indiv <- 0.95
logM_draws <- logM_draws + qnorm(q_indiv) * sigma_e

mass_draws_kg <- exp(logM_draws)
mass_up95_kg  <- as.numeric(quantile(mass_draws_kg, probs = ci_level))
dens_m <- density(mass_draws_kg, n = 4096, adjust = 1.0)
mass_map_kg <- dens_m$x[which.max(dens_m$y)]

p_mass <- ggplot(data.frame(W = mass_draws_kg), aes(W)) +
  geom_density(linewidth = 1.2, color = "darkblue") +
  geom_vline(xintercept = mass_map_kg,  color = "purple",  linetype = "dashed",  linewidth = 1) +
  geom_vline(xintercept = mass_up95_kg, color = "orange",  linetype = "dotdash", linewidth = 1) +
  geom_vline(xintercept = stokes_w_kg,  color = "firebrick", linetype = "dashed", linewidth = 1) +
  annotate("text", x = stokes_w_kg, y = max(dens_m$y)*0.08, vjust = -0.8,
           label = "Stokes (459 kg)", color = "firebrick", angle = 90, size = 4) +
  labs(x = "Weight* (kg)", y = "Density") +
  theme_science_polished
ggsave(file.path(fig_dir, "mass_endpoint_females.png"),
       p_mass, dpi = 600, width = 7, height = 5, units = "in")

cat(sprintf("\nFEMALE mass endpoint (typical individual, q=%d%%):\n", round(q_indiv*100)))
cat(sprintf("MAP: %.2f kg | %d%% upper: %.2f kg\n",
            mass_map_kg, round(ci_level*100), mass_up95_kg))

## ---------------------------------------------------------
## 8b) Regression panel (FEMALES ONLY) + Lake Apopka annotation
## ---------------------------------------------------------
plot_df_f <- dat %>%
  filter(is.finite(WTkg), is.finite(TL)) %>%
  transmute(logTL = log(TL),
            logW  = log(WTkg))

y_at_u0 <- as.numeric(predict(lm_w_tl, newdata = data.frame(logTL = u0)))
x_ep    <- ystar_map
x_ep_ub <- ystar_upper
y_ep    <- log(mass_map_kg)
y_ep_ub <- log(mass_up95_kg)

x_stokes  <- log(stokes_tl_cm)
y_stokes  <- log(stokes_w_kg)
x_apopka  <- log(apopka_tl_cm)
y_apopka  <- log(apopka_w_kg)
ab        <- coef(lm_w_tl)

p_reg_f <- ggplot(plot_df_f, aes(x = logTL, y = logW)) +
  geom_point(alpha = 0.85, size = 2.2, color = "#ff6fb3") +            # pink female dots
  geom_abline(intercept = ab[1], slope = ab[2], linewidth = 1) +       # female regression
  geom_vline(xintercept = u0, linetype = "dashed", color = "gray30") + # threshold (x)
  geom_hline(yintercept = y_at_u0, linetype = "dashed", color = "gray30") + # threshold (y)
  geom_point(aes(x = x_ep, y = y_ep), color = "darkgreen", size = 3, inherit.aes = FALSE) +
  geom_segment(aes(x = x_ep, y = y_ep, xend = x_ep_ub, yend = y_ep),
               color = "darkgreen", linewidth = 1.1, inherit.aes = FALSE) +
  geom_segment(aes(x = x_ep, y = y_ep, xend = x_ep, yend = y_ep_ub),
               color = "darkgreen", linewidth = 1.1, inherit.aes = FALSE) +
  geom_point(aes(x = x_stokes, y = y_stokes),
             inherit.aes = FALSE, color = "firebrick", size = 2.8) +
  geom_text(aes(x = x_stokes, y = y_stokes, label = "Stokes"),
            inherit.aes = FALSE, nudge_y = 0.015, color = "firebrick", size = 3.5) +
  # Lake Apopka label (filled diamond to stand out)
  geom_point(aes(x = x_apopka, y = y_apopka),
             inherit.aes = FALSE, shape = 23, size = 3, fill = "black", color = "white", stroke = 0.4) +
  geom_text(aes(x = x_apopka, y = y_apopka, label = "Lake Apopka"),
            inherit.aes = FALSE, nudge_y = 0.03, size = 3.4) +
  labs(title = "log Weight vs log TL (FEMALES)",
       x = "log Total Length", y = "log Weight") +
  theme_science_polished +
  theme(legend.position = "none")

ggsave(file.path(fig_dir, "reg_logW_vs_logTL_females_apopka.png"),
       p_reg_f, dpi = 600, width = 8, height = 5.2, units = "in")

## ---------------------------------------------------------
## 9) Exceedance probability of Stokes TL (reuses female grids)
## ---------------------------------------------------------
trapz_weights <- function(x, f_un) {
  stopifnot(length(x) == length(f_un), length(x) >= 2)
  dx <- diff(x)
  mids <- (f_un[-1] + f_un[-length(f_un)]) * dx / 2
  total <- sum(mids)
  if (!is.finite(total) || total <= 0) return(rep(NA_real_, length(x)))
  w <- numeric(length(x))
  w[1]                 <- mids[1] / 2
  w[length(x)]         <- mids[length(mids)] / 2
  if (length(x) > 2) w[2:(length(x)-1)] <- (mids[-length(mids)] + mids[-1]) / 2
  w / sum(w)
}

W_y <- trapz_weights(y_star_grid, post)
stopifnot(all(is.finite(W_y)), abs(sum(W_y) - 1) < 1e-6)

p_exceed_u   <- mean(y > u0, na.rm = TRUE)
y_claim_cm   <- stokes_tl_cm
y_claim_log  <- log(y_claim_cm)

Si_vec <- numeric(length(y_star_grid))
for (i in seq_along(y_star_grid)) {
  ys <- y_star_grid[i]
  if (ys <= y_claim_log) { Si_vec[i] <- 0; next }
  lv <- sapply(xi_grid, function(xi) {
    loglik_endpoint(ys, xi, u0, y_above) + logprior_xi(xi)
  })
  m  <- max(lv)
  fu <- exp(lv - m)
  V_xi <- trapz_weights(xi_grid, fu)
  surv_term <- ((ys - y_claim_log) / (ys - u0))^(-1/xi_grid)
  surv_term[!is.finite(surv_term)] <- 0
  Si_vec[i] <- sum(V_xi * surv_term)
}
p_tail_cond     <- sum(W_y * Si_vec)
p_exceed_claim  <- p_exceed_u * p_tail_cond

cat(sprintf("\n--- Exceedance probability for Stokes TL (FEMALES) ---\n"))
cat(sprintf("Threshold u0 (log cm)  : %.4f  (≈ %.1f cm)\n", u0, exp(u0)))
cat(sprintf("Claim y (log cm)       : %.4f  (450.0 cm)\n", y_claim_log))
cat(sprintf("P(Y > u0 | data)       : %.4f\n", p_exceed_u))
cat(sprintf("P(Y > y | Y>u0, data)  : %.4f\n", p_tail_cond))
cat(sprintf("P(Y > y | data)        : %.6f\n", p_exceed_claim))
