#### =========================
#### EVT on log Total Length (FEMALES) with point imputation from log SVL → log TL
#### Ordered & commented; Woodward prior + LOO cross-validation (fast options)
#### =========================

## ---------------------------
## 0) Libraries & global setup
## ---------------------------
library(ggplot2)
library(extRemes)    # fevd/revd/pevd for GP fits & sims
library(dplyr)
library(readxl)
library(pracma)      # trapz / cumtrapz (numeric integration)
library(HDInterval)  # hdi (used only for quick intervals on samples if needed)

set.seed(42)

## Ensure figure directory exists
fig_dir <- "Figures"
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

## ---------------------------
## 1) Settings & theme
## ---------------------------
ci_level      <- 0.95
n_boot_ks     <- 1000
thresh_q_opt  <- 0.98       # u0 = q0.95(log TL) — anchor in stability region

# Lake Apopka reference (female)
apopka_tl_cm  <- 322
apopka_w_kg   <- 170

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
## 2) Data ingest
## ---------------------------
data_path <- "Data/experimental_alligator_harvest_woodward.xlsx"
dat0 <- read_excel(data_path)

# Keep females only as analysis base; keep original for fallbacks if needed
datF <- dat0 %>% mutate(Sex = as.character(Sex)) %>% filter(Sex %in% c("F","f","Female","female","Fem","FEMALE"))
datF$Sex <- "F"

# Ensure expected columns exist (no-op if already present)
datF <- datF %>% mutate(SVL, TL, WTkg, Deform)

## ---------------------------
## 3) Helpers
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

make_thresholds <- function(y, q_from = 0.8, q_to = 0.99, n = 50){
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
  list(ks_obs = as.numeric(ks_obs), p_value = if (length(ks_b)) mean(ks_b >= ks_obs) else NA_real_)
}

# Trapezoid point-weights that sum to 1 for (x, f(x))
trapz_weights <- function(x, f_un) {
  stopifnot(length(x) == length(f_un), length(x) >= 2)
  dx <- diff(x)
  mids <- (f_un[-1] + f_un[-length(f_un)]) * dx / 2
  total <- sum(mids); if (!is.finite(total) || total <= 0) return(rep(NA_real_, length(x)))
  w <- numeric(length(x))
  w[1] <- mids[1] / 2
  w[length(x)] <- mids[length(mids)] / 2
  if (length(x) > 2) w[2:(length(x)-1)] <- (mids[-length(mids)] + mids[-1]) / 2
  w / sum(w)
}

cdf_from_density <- function(x, f) {
  f[!is.finite(f)] <- 0
  F <- pracma::cumtrapz(x, f)
  F / max(F, na.rm = TRUE)
}
q_from_cdf <- function(x, f, p) {
  F <- cdf_from_density(x, f)
  as.numeric(approx(x = F, y = x, xout = p, ties = "ordered", rule = 2)$y)
}

## ---------------------------------------------------------
## 4) Impute missing TL (SVL→TL, log–log OLS) — FEMALES
## ---------------------------------------------------------
reg_df_F <- datF %>%
  filter(Deform == 0L, is.finite(TL), is.finite(SVL)) %>%
  transmute(logTL = log(TL), logSVL = log(SVL))

# If too few female pairs, fallback to all undeformed pairs
if (nrow(reg_df_F) < 10) {
  reg_df_F <- dat0 %>%
    filter(Deform == 0L, is.finite(TL), is.finite(SVL)) %>%
    transmute(logTL = log(TL), logSVL = log(SVL))
}
stopifnot(nrow(reg_df_F) >= 10)

lm_tl_svl <- lm(logTL ~ logSVL, data = reg_df_F)

datF <- datF %>%
  mutate(
    logSVL    = ifelse(is.finite(SVL), log(SVL), NA_real_),
    logTL_obs = ifelse(is.finite(TL),  log(TL),  NA_real_),
    logTL_pred= ifelse(is.finite(logSVL),
                       as.numeric(predict(lm_tl_svl, newdata = data.frame(logSVL = logSVL))),
                       NA_real_),
    logTL_use = ifelse(is.finite(logTL_obs), logTL_obs, logTL_pred)
  )

y <- datF$logTL_use
y <- y[is.finite(y)]
stopifnot(length(y) >= 40)  # female sample can be smaller than males

## ---------------------------------------------------------
## 5) Threshold selection (MRL; shape & adjusted-scale stability)
## ---------------------------------------------------------
u_seq <- make_thresholds(y, q_from = 0.6, q_to = 0.99, n = 50)
u0    <- as.numeric(quantile(y, thresh_q_opt, na.rm = TRUE))  # chosen anchor
u0 <- log(253)

# ---- 5a) Mean Residual Life (MRL) plot ----
mrl_vals <- mrl_data(y, u_seq)

p_mrl <- ggplot(data.frame(u = u_seq, mrl = mrl_vals), aes(u, mrl)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2) +
  geom_vline(xintercept = u0, linetype = "dashed", color = "red") +
  labs(title = "MRL plot (log TL; females, point imputation)",
       x = "Threshold u (log TL)", y = "Mean excess E[Y-u | Y>u]") +
  theme_science_polished
p_mrl
ggsave(file.path(fig_dir, "F_tl_log_mrl.png"),
       p_mrl, dpi = 600, width = 7, height = 5, units = "in")

# ---- 5b) Shape (xi) stability with Wald CIs at each u ----
shape <- scale <- rep(NA_real_, length(u_seq))
shape_lo <- shape_hi <- rep(NA_real_, length(u_seq))
adj <- adj_lo <- adj_hi <- rep(NA_real_, length(u_seq))

# Reference fit at u0 (adj at u0 equals scale at u0)
fit0 <- fit_gpd_at_u(y, u0)
adj_hat0 <- fit0$scale
shape_hat0 <- fit0$shape

for (i in seq_along(u_seq)) {
  u <- u_seq[i]
  exc <- y[y > u] - u
  if (length(exc) < 15) next  # fewer females → allow slightly lower cutoff
  
  out <- try(fit_gpd_at_u(y, u), silent = TRUE)
  if (inherits(out, "try-error")) next
  
  scale[i] <- out$scale
  shape[i] <- out$shape
  
  zcrit    <- qnorm(1 - (1 - ci_level)/2)
  se_scale <- sqrt(out$cov[1,1])
  se_shape <- sqrt(out$cov[2,2])
  
  shape_lo[i] <- shape[i] - zcrit * se_shape
  shape_hi[i] <- shape[i] + zcrit * se_shape
  
  # Adjust scale to common anchor u0
  adj[i] <- adj_scale(scale[i], shape[i], u, u0)
  var_adj <- out$cov[1,1] + (u - u0)^2 * out$cov[2,2] - 2*(u - u0)*out$cov[1,2]
  se_adj  <- sqrt(max(var_adj, 0))
  adj_lo[i] <- adj[i] - zcrit * se_adj
  adj_hi[i] <- adj[i] + zcrit * se_adj
}

# Shape stability plot
df_shape <- data.frame(u = u_seq, shape, shape_lo, shape_hi)
p_shape <- ggplot(df_shape, aes(u, shape)) +
  geom_point(na.rm = TRUE) +
  geom_errorbar(aes(ymin = shape_lo, ymax = shape_hi), width = 0.03,
                color = "steelblue", na.rm = TRUE) +
  geom_vline(xintercept = u0, color = "red", linetype = "dashed") +
  geom_hline(yintercept = shape_hat0, color = "red", linetype = "dotted") +
  labs(
       x = "Threshold u (log TL)", y = "Shape") +
  theme_science_polished
p_shape
ggsave(file.path(fig_dir, "F_tl_log_shape_stability.png"),
       p_shape, dpi = 600, width = 7, height = 5, units = "in")

# Adjusted-scale stability plot (σ̃ at common anchor u0)
df_adj <- data.frame(u = u_seq, adj, adj_lo, adj_hi)
p_adj <- ggplot(df_adj, aes(u, adj)) +
  geom_point(na.rm = TRUE) +
  geom_errorbar(aes(ymin = adj_lo, ymax = adj_hi), width = 0.03,
                color = "steelblue", na.rm = TRUE) +
  geom_vline(xintercept = u0, color = "red", linetype = "dashed") +
  geom_hline(yintercept = adj_hat0, color = "red", linetype = "dotted") +
  labs(
       x = "Threshold u (log TL)", y = "Adjusted scale") +
  theme_science_polished
ggsave(file.path(fig_dir, "F_tl_log_adj_scale_stability.png"),
       p_adj, dpi = 600, width = 7, height = 5, units = "in")

## ---------------------------------------------------------
## 6) Tail diagnostics
## ---------------------------------------------------------
diags <- diagnostic_plots(y, u0, fit0$scale, fit0$shape, "log TL (F)")
ggsave(file.path(fig_dir, "F_tl_log_qq.png"), diags$pqq, dpi = 600, width = 7, height = 5, units = "in")
ggsave(file.path(fig_dir, "F_tl_log_pp.png"), diags$ppp, dpi = 600, width = 7, height = 5, units = "in")
ggsave(file.path(fig_dir, "F_tl_log_pit.png"), diags$pks, dpi = 600, width = 7, height = 5, units = "in")
ksout <- ks_boot_pvalue(y, u0, fit0$scale, fit0$shape, B = n_boot_ks)
shape_se0   <- sqrt(fit0$cov[2,2])
z       <- shape_hat0 / shape_se0
p_wald  <- pnorm(z)  # one-sided for H1: xi<0

message(sprintf("[F log TL] u0=%.3f | scale=%.4f shape=%.4f | KS p=%.3f | Wald p(one-sided xi<0)=%.3g",
                u0, fit0$scale, fit0$shape, ksout$p_value, p_wald))

## ---------------------------------------------------------
## 7) Endpoint posterior with Woodward prior (grid-only)
## ---------------------------------------------------------
y_above <- y[y > u0]; n_ex <- length(y_above)
if (n_ex < 10) stop("Not enough exceedances for endpoint posterior (F).")

Ymax <- max(y_above)
eps  <- max(1e-6, 1e-4 * abs(Ymax))
y_star_grid <- seq(Ymax + eps, Ymax + 10, length.out = 2000)    # grid for y* (log TL*)
xi_grid     <- seq(-1.0, -0.02, length.out = 1000)              # xi < 0

loglik_endpoint <- function(y_star, xi, u, y_ex) {
  if (xi >= 0 || any(y_star <= y_ex)) return(-Inf)
  n <- length(y_ex)
  n * log(y_star - u) / xi - (1/xi + 1) * sum(log(y_star - y_ex)) - n * log(-xi)
}

# Priors (kept as in the male script unless you want a female-centered prior)
mu_y <- log(332)        # prior center (can be lowered for females if desired)
prior_mode <- "woodward_mod"
sd_y <- switch(prior_mode,
               "woodward_loose" = 0.096,
               "woodward_mod"   = 0.048,
               "woodward_tight" = 0.024, 0.048)
logprior_y  <- function(y_star) dnorm(y_star, mean = mu_y, sd = sd_y, log = TRUE)

logprior_xi <- function(xi) {
  if (xi >= 0 || xi <= -1) return(-Inf)
  0  # improper flat on (-1,0)
}

# Marginalize xi for each y* (trapz over xi), then multiply by prior on y*
marg_log_post <- sapply(y_star_grid, function(ys) {
  lv <- sapply(xi_grid, function(xi) loglik_endpoint(ys, xi, u0, y_above) + logprior_xi(xi))
  lv_max <- max(lv)
  log_int_xi <- log(pracma::trapz(xi_grid, exp(lv - lv_max))) + lv_max
  log_int_xi + logprior_y(ys)
})

# Normalize over y*
post_unnorm <- exp(marg_log_post - max(marg_log_post))
Z_y         <- pracma::trapz(y_star_grid, post_unnorm)
if (!is.finite(Z_y) || Z_y <= 0) stop("Normalization failed for p_F(y*|y).")
post <- post_unnorm / Z_y

# Posterior summaries for y*
ystar_map   <- y_star_grid[which.max(post)]
ystar_upper <- q_from_cdf(y_star_grid, post, ci_level)

tl_map_cm   <- exp(ystar_map)
tl_upper_cm <- exp(ystar_upper)

cat("\n--- FEMALE TL* (cm) via numeric integration ---\n",
    "MAP:    ", round(tl_map_cm, 1), "\n",
    sprintf("%d%% upper:", round(100*ci_level)), round(tl_upper_cm, 1), "\n")

## ---------------------------------------------------------
## 7b) Posterior plot for TL* (cm scale) — FEMALES
## ---------------------------------------------------------
# Change-of-variables: f_TL(t) = f_Y(log t) * 1/t
tl_grid_cm <- exp(y_star_grid)
dens_tl    <- post / tl_grid_cm
dens_tl    <- dens_tl / pracma::trapz(tl_grid_cm, dens_tl)

df_tlpost <- data.frame(TL = tl_grid_cm, dens = dens_tl)

p_tlstar <- ggplot(df_tlpost, aes(TL, dens)) +
  geom_line(linewidth = 1.2, color = "darkblue") +
  geom_vline(xintercept = tl_map_cm,
             color = "purple",   linetype = "dashed",  linewidth = 1.1) +
  geom_vline(xintercept = tl_upper_cm,
             color = "orange",   linetype = "dotdash", linewidth = 1.0) +
  geom_vline(xintercept = apopka_tl_cm,
             color = "firebrick", linetype = "dashed", linewidth = 1.0) +
  annotate("text", x = apopka_tl_cm, y = max(dens_tl)*0.2, vjust = -0.8,
           label = "Lake Apopka (322 cm)", color = "firebrick", angle = 90, size = 4) +
  labs(x = "Total length TL* (cm)", y = "Density") +
  xlim(297,400)+
  theme_science_polished
p_tlstar
ggsave(file.path(fig_dir, "F_tl_endpoint.png"),
       p_tlstar, dpi = 600, width = 7, height = 5, units = "in")

cat("\n--- FEMALE TL endpoint (cm) ---\n",
    "MAP:    ", round(tl_map_cm, 1), "\n",
    sprintf("%d%% upper:", round(100*ci_level)), round(tl_upper_cm, 1), "\n")

## =========================================================
## 8) Mass propagation via (S7): numeric mixture (NO sampling) — FEMALES
## =========================================================
# Fit allometric regression: logW ~ logTL on female undeformed with both observed
w_fit_df_F <- datF %>%
  filter(Deform == 0L, is.finite(WTkg), is.finite(TL)) %>%
  transmute(logW = log(WTkg), logTL = log(TL))

# Fallback: any females with both; final fallback: all individuals
if (nrow(w_fit_df_F) < 10) {
  w_fit_df_F <- datF %>% filter(is.finite(WTkg), is.finite(TL)) %>%
    transmute(logW = log(WTkg), logTL = log(TL))
}
if (nrow(w_fit_df_F) < 10) {
  w_fit_df_F <- dat0 %>% filter(is.finite(WTkg), is.finite(TL)) %>%
    transmute(logW = log(WTkg), logTL = log(TL))
}
stopifnot(nrow(w_fit_df_F) >= 10)

lm_w_tl <- lm(logW ~ logTL, data = w_fit_df_F)
beta    <- coef(lm_w_tl)           # (a, b)
Vcoef   <- vcov(lm_w_tl)
sigma_e <- summary(lm_w_tl)$sigma  # residual SD on logW (individual predictive)

# Trapezoidal weights on y* posterior
W_y <- trapz_weights(y_star_grid, post)
stopifnot(all(is.finite(W_y)), abs(sum(W_y) - 1) < 1e-6)

# For each y*_i: predictive mean & variance on log mass
mu_i   <- beta[1] + beta[2] * y_star_grid
x1     <- rep(1, length(y_star_grid))
quad   <- Vcoef[1,1]*x1^2 + 2*Vcoef[1,2]*x1*y_star_grid + Vcoef[2,2]*y_star_grid^2
sig2_i <- pmax(quad + sigma_e^2, 0)

# Mixture on m = log mass
m_min <- min(mu_i - 6*sqrt(sig2_i))
m_max <- max(mu_i + 6*sqrt(sig2_i))
m_grid <- seq(m_min, m_max, length.out = 4000)

dnorm_vec <- function(x, mean, sd) (1/(sd*sqrt(2*pi))) * exp(-0.5*((x-mean)/sd)^2)

f_m <- rowSums( sapply(seq_along(mu_i), function(i) {
  W_y[i] * dnorm_vec(m_grid, mu_i[i], sqrt(sig2_i[i]))
}) )

f_m <- f_m / pracma::trapz(m_grid, f_m)

# Summaries on log scale
F_m <- cdf_from_density(m_grid, f_m)
m_map <- m_grid[which.max(f_m)]
m_up  <- as.numeric(approx(x = F_m, y = m_grid, xout = ci_level, ties = "ordered", rule = 2)$y)

# Transform to kg (Jacobian)
w_grid <- exp(m_grid)
f_w    <- f_m / w_grid
f_w    <- f_w / pracma::trapz(w_grid, f_w)
F_w    <- cdf_from_density(w_grid, f_w)

w_map  <- w_grid[which.max(f_w)]
w_up   <- as.numeric(approx(x = F_w, y = w_grid, xout = ci_level, ties = "ordered", rule = 2)$y)

cat(sprintf("\n--- FEMALE Weight endpoint (typical individual) via (S7) ---\n"))
cat(sprintf("MAP: %.1f kg | %d%% upper: %.1f kg\n", w_map, round(100*ci_level), w_up))

# Plot weight density with MAP + upper + Lake Apopka
df_w <- data.frame(w = w_grid, dens = f_w)
p_w <- ggplot(df_w, aes(w, dens)) +
  geom_line(linewidth = 1.2, color = "darkblue") +
  geom_vline(xintercept = w_map, color = "purple", linetype = "dashed", linewidth = 1) +
  geom_vline(xintercept = w_up,  color = "orange", linetype = "dotdash", linewidth = 1) +
  geom_vline(xintercept = apopka_w_kg, color = "firebrick", linetype = "dashed", linewidth = 1) +
  annotate("text", x = apopka_w_kg, y = max(df_w$dens)*0.2, vjust = -0.8,
           label = "Lake Apopka (170 kg)", color = "firebrick", angle = 90, size = 4) +
  labs(x = "Weight (kg)", y = "Density") +
  theme_science_polished + 
  xlim(50,250)
p_w
ggsave(file.path(fig_dir, "F_mass_endpoint.png"),
       p_w, dpi = 600, width = 7, height = 5, units = "in")

## ---------------------------------------------------------
## 8b) Regression panel: log weight vs log TL (+ endpoint & Lake Apopka)
## ---------------------------------------------------------
plot_df <- datF %>%
  filter(is.finite(WTkg), is.finite(TL)) %>%
  transmute(logTL = log(TL), logW = log(WTkg))

ab <- coef(lm_w_tl)
y_at_u0 <- as.numeric(predict(lm_w_tl, newdata = data.frame(logTL = u0)))

# Endpoint point ON the line at x=MAP(y*)
x_ep    <- ystar_map
x_ep_ub <- ystar_upper
y_ep    <- log(w_map)
y_ep_ub <- log(w_up)

x_apopka <- log(apopka_tl_cm)
y_apopka <- log(apopka_w_kg)

p_reg <- ggplot(plot_df, aes(x = logTL, y = logW)) +
  geom_point(alpha = 0.7, size = 2, color = "#ff6fb3") +
  geom_abline(intercept = ab[1], slope = ab[2], linewidth = 1) +
  geom_vline(xintercept = u0, linetype = "dashed", color = "gray30") +
  geom_hline(yintercept = y_at_u0, linetype = "dashed", color = "gray30") +
  geom_point(aes(x = x_ep, y = y_ep), color = "darkgreen", size = 3, inherit.aes = FALSE) +
  geom_segment(aes(x = x_ep, y = y_ep, xend = x_ep_ub, yend = y_ep),
               color = "darkgreen", linewidth = 1.1, inherit.aes = FALSE) +
  geom_segment(aes(x = x_ep, y = y_ep, xend = x_ep, yend = y_ep_ub),
               color = "darkgreen", linewidth = 1.1, inherit.aes = FALSE) +
  geom_point(aes(x = x_apopka, y = y_apopka),
             inherit.aes = FALSE, color = "firebrick", size = 2.8) +
  labs(x = "log(Total Length [cm])", y = "log(Mass [kg])") +
  theme_science_polished
p_reg
ggsave(file.path(fig_dir, "F_regression_panel.png"),
       p_reg, dpi = 600, width = 8, height = 5.2, units = "in")

## ---------------------------------------------------------
## 9) Exceedance probability for Lake Apopka TL (as in S8–S9)
## ---------------------------------------------------------
p_exceed_u <- mean(y > u0, na.rm = TRUE)           # Pr(Y>u0)
y_claim_log <- log(apopka_tl_cm)

Si_vec <- numeric(length(y_star_grid))
for (i in seq_along(y_star_grid)) {
  ys <- y_star_grid[i]
  if (ys <= y_claim_log) { Si_vec[i] <- 0; next }
  lv <- sapply(xi_grid, function(xi) loglik_endpoint(ys, xi, u0, y_above) + logprior_xi(xi))
  m  <- max(lv); fu <- exp(lv - m)
  V_xi <- trapz_weights(xi_grid, fu)                 # p(xi | ys, data) up to norm
  surv_term <- ((ys - y_claim_log)/(ys - u0))^(-1/xi_grid)  # GP survival above u0
  surv_term[!is.finite(surv_term)] <- 0
  Si_vec[i] <- sum(V_xi * surv_term)
}

p_tail_cond <- sum(trapz_weights(y_star_grid, post) * Si_vec)  # E_{y*|data}[ E_{xi|y*,data}[...] ]
p_exceed_claim <- p_exceed_u * p_tail_cond
cat(sprintf("\n[FEMALES] P(Y > %d cm | data) ≈ %.2e\n", apopka_tl_cm, p_exceed_claim))

