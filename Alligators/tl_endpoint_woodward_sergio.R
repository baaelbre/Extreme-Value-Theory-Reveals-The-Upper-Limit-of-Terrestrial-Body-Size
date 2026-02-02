#### =========================
#### EVT on log Total Length with point imputation from log SVL → log TL
#### Analysis on WOODWARD; mass mapping from renewed CSV
#### Woodward prior + LOO cross-validation (fast options)
#### =========================

## ---------------------------
## 0) Libraries & global setup
## ---------------------------
library(ggplot2)
library(extRemes)    # fevd/revd/pevd for GP fits & sims
library(dplyr)
library(readxl)
library(pracma)      # trapz for numerical integration
library(HDInterval)  # hdi
library(MASS)        # mvrnorm

set.seed(42)

fig_dir <- "Figures"
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

## ---------------------------
## 1) Settings & theme
## ---------------------------
ci_level      <- 0.95
n_boot_ks     <- 1000
n_post_ystar  <- 20000
thresh_q_opt_allan  <- 0.95  # u0_allan = q_0.94(log TL) — tune if needed
thresh_q_opt_sergio <- 0.97

stokes_tl_cm  <- 450
stokes_w_kg   <- 459

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
## 2a) Renewed allometry CSV — used for TL→Weight regression ONLY
csv_path <- "Data/CaptureData_Gator_allometry_paper_2024.csv"
dat_raw  <- read.csv(csv_path, stringsAsFactors = FALSE)

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
  mutate(
    SVL  = ifelse(is.finite(SVL)  & SVL  > 0, SVL,  NA_real_),
    TL   = ifelse(is.finite(TL)   & TL   > 0, TL,   NA_real_),
    WTkg = ifelse(is.finite(WTkg) & WTkg > 0, WTkg, NA_real_)
  )

## 2b) WOODWARD dataset — used for EVT (TL) analysis
wood_path <- "Data/experimental_alligator_harvest_woodward.xlsx"
wood <- readxl::read_excel(wood_path) %>%
  dplyr::mutate(SVL, TL, WTkg, Deform) %>%
  dplyr::mutate(
    TL    = ifelse(is.finite(TL)   & TL   > 0, TL,   NA_real_),
    WTkg  = ifelse(is.finite(WTkg) & WTkg > 0, WTkg, NA_real_),
    Deform = dplyr::case_when(is.na(Deform) ~ 0L, TRUE ~ as.integer(Deform))
  )

## ---------------------------
## 2c) Print & count TL ≥ 350 cm in BOTH datasets
## ---------------------------
cut_cm <- 327
cat("\n========== TL ≥ 350 cm (WOODWARD) ==========\n")
wood_350 <- wood %>% filter(is.finite(TL), TL >= cut_cm)
print(wood_350)
cat(sprintf("Count (WOODWARD, TL ≥ %d cm): %d\n", cut_cm, nrow(wood_350)))

cat("\n========== TL ≥ 350 cm (SERGIO) ==========\n")
dat_350 <- dat %>% filter(is.finite(TL), TL >= cut_cm)
print(dat_350)
cat(sprintf("Count (CSV, TL ≥ %d cm): %d\n\n", cut_cm, nrow(dat_350)))

cut_kg <- 177
cat("\n========== W ≥ 300 kg (WOODWARD) ==========\n")
wood_300kg <- wood %>% filter(is.finite(WTkg), WTkg >= cut_kg)
print(wood_300kg)
cat(sprintf("Count (WOODWARD, W ≥ %d kg): %d\n", cut_kg, nrow(wood_300kg)))

cat("\n========== TL ≥ 350 cm (SERGIO) ==========\n")
dat_300kg <- dat %>% filter(is.finite(Weight), Weight >= cut_kg)
print(dat_300kg)
cat(sprintf("Count (CSV, W ≥ %d kg): %d\n\n", cut_kg, nrow(dat_300kg)))

# (Optional) write to disk for inspection
# write.csv(wood_350, file.path(fig_dir, "woodward_TL_ge_350.csv"), row.names = FALSE)
# write.csv(dat_350,  file.path(fig_dir, "csv_TL_ge_350.csv"),      row.names = FALSE)

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

mrl_data <- function(y, u_seq_allan) {
  sapply(u_seq_allan, function(u){
    exc <- y[y > u] - u
    if (length(exc) < 5) return(NA_real_)
    mean(exc)
  })
}

make_thresholds <- function(y, q_from = 0.8, q_to = 0.999, n = 50){
  rng <- quantile(y, c(q_from, q_to), na.rm = TRUE)
  seq(rng[1], rng[2], length.out = n)
}

adj_scale <- function(scale, shape, u, u0_allan) scale - shape * (u - u0_allan)

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
  
  exc <- y[y > u] - u
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
## 4) Impute missing TL via SVL→TL (WOODWARD, used for EVT)
## ---------------------------------------------------------
reg_df_wood <- wood %>%
  filter(Deform == 0L, is.finite(TL), is.finite(SVL)) %>%
  transmute(logTL = log(TL), logSVL = log(SVL))
stopifnot(nrow(reg_df_wood) >= 10)

lm_tl_svl_wood <- lm(logTL ~ logSVL, data = reg_df_wood)

wood <- wood %>%
  mutate(
    logSVL    = ifelse(is.finite(SVL), log(SVL), NA_real_),
    logTL_obs = ifelse(is.finite(TL),  log(TL),  NA_real_),
    logTL_pred= ifelse(is.finite(logSVL),
                       as.numeric(predict(lm_tl_svl_wood, newdata = data.frame(logSVL = logSVL))),
                       NA_real_),
    logTL_use = ifelse(is.finite(logTL_obs), logTL_obs, logTL_pred)
  )

# Tail variable for POT — WOODWARD
y_allan <- wood$logTL_use
y_allan <- y_allan[is.finite(y_allan)]
stopifnot(length(y_allan) >= 50)

# SERGIO
y_sergio <- log(dat$TL)
y_sergio <- y_sergio[is.finite(y_sergio)]
stopifnot(length(y_sergio) >= 50)
## ---------------------------------------------------------
## 5) Threshold selection
## ---------------------------------------------------------
# Woodward
u_seq_allan <- make_thresholds(y_allan, q_from = 0.7, q_to = 0.99, n = 50)
u0_allan    <- as.numeric(quantile(y_allan, thresh_q_opt_allan, na.rm = TRUE))
fit0_allan  <- fit_gpd_at_u(y_allan, u0_allan)
shape0_allan <- fit0_allan$shape
adj0_allan   <- fit0_allan$scale

shape0_se_allan <- sqrt(fit0_allan$cov[2,2])
z_allan         <- shape0_allan / shape0_se_allan
p_wald_allan    <- pnorm(z_allan)  # one-sided H1: xi<0
p_wald_allan

# Balaguera-Reina
u_seq_sergio <- make_thresholds(y_sergio, q_from = 0.7, q_to = 0.99, n = 50)
u0_sergio    <- as.numeric(quantile(y_sergio, thresh_q_opt_sergio, na.rm = TRUE))
fit0_sergio  <- fit_gpd_at_u(y_sergio, u0_sergio)
shape0_sergio <- fit0_sergio$shape
adj0_sergio   <- fit0_sergio$scale

shape0_se_sergio <- sqrt(fit0_sergio$cov[2,2])
z_sergio         <- shape0_sergio / shape0_se_sergio
p_wald_sergio    <- pnorm(z_sergio)  # one-sided H1: xi<0
p_wald_sergio

# Woodward
mrl_vals_allan <- mrl_data(y_allan, u_seq_allan)
p_mrl_allan <- ggplot(data.frame(u = u_seq_allan, mrl = mrl_vals_allan), aes(u, mrl)) +
  geom_line() + geom_point() +
  geom_vline(xintercept = u0_allan, linetype = "dashed", color = "red") +
  labs(title = "MRL plot (log TL, WOODWARD)",
       x = "Threshold (log TL)", y = "Mean excess") +
  theme_science_polished
p_mrl_allan
ggsave(file.path(fig_dir, "wood_tl_log_mrl.png"), p_mrl_allan, dpi = 600, width = 7, height = 5)

shape_allan <- scale_allan <- shape_lo_allan <- shape_hi_allan <- adj_allan <- adj_lo_allan <- adj_hi_allan <- rep(NA_real_, length(u_seq_allan))
for (i in seq_along(u_seq_allan)) {
  u_allan <- u_seq_allan[i]
  exc_allan <- y_allan[y_allan > u_allan] - u_allan
  if (length(exc_allan) < 20) next
  out_allan <- try(fit_gpd_at_u(y_allan, u_allan), silent = TRUE)
  if (inherits(out_allan, "try-error")) next
  scale_allan[i] <- out_allan$scale; shape_allan[i] <- out_allan$shape
  zcrit     <- qnorm(1 - (1 - ci_level)/2)
  se_scale_allan  <- sqrt(out_allan$cov[1,1]); se_shape_allan <- sqrt(out_allan$cov[2,2])
  shape_lo_allan[i] <- shape_allan[i] - zcrit * se_shape_allan
  shape_hi_allan[i] <- shape_allan[i] + zcrit * se_shape_allan
  adj_allan[i] <- adj_scale(scale_allan[i], shape_allan[i], u_allan, u0_allan)
  var_adj_allan <- out_allan$cov[1,1] + (u_allan - u0_allan)^2 * out_allan$cov[2,2] - 2*(u_allan - u0_allan)*out_allan$cov[1,2]
  se_adj_allan  <- sqrt(max(var_adj_allan, 0))
  adj_lo_allan[i] <- adj_allan[i] - zcrit * se_adj_allan
  adj_hi_allan[i] <- adj_allan[i] + zcrit * se_adj_allan
}

p_shape_allan <- ggplot(data.frame(u = u_seq_allan, shape_allan, shape_lo_allan, shape_hi_allan), aes(u, shape_allan)) +
  geom_point() +
  geom_errorbar(aes(ymin = shape_lo_allan, ymax = shape_hi_allan), width = 0.03, color = "blue") +
  geom_vline(xintercept = u0_allan, color = "red", linetype = "dashed") +
  labs(x = "Threshold (log TL)", y = "Shape (xi)") +
  theme_science_polished + 
  geom_hline(yintercept = shape0_allan, color = "red", linetype = "dashed")
ggsave(file.path(fig_dir, "shape_stability_allan.png"), p_shape_allan, dpi = 600, width = 7, height = 5)
p_shape_allan

p_adj_allan <- ggplot(data.frame(u = u_seq_allan, adj_allan, adj_lo_allan, adj_hi_allan), aes(u, adj_allan)) +
  geom_point() +
  geom_errorbar(aes(ymin = adj_lo_allan, ymax = adj_hi_allan), width = 0.03, color = "blue") +
  geom_vline(xintercept = u0_allan, color = "red", linetype = "dashed") +
  labs(x = "Threshold (log TL)", y = "Adjusted scale") +
  theme_science_polished + 
  geom_hline(yintercept = adj0_allan, color = "red", linetype = "dashed")
p_adj_allan
ggsave(file.path(fig_dir, "adj_scale_stability_allan.png"), p_adj_allan, dpi = 600, width = 7, height = 5)


# Balaguera-Reina
mrl_vals_sergio <- mrl_data(y_sergio, u_seq_sergio)
p_mrl_sergio <- ggplot(data.frame(u = u_seq_sergio, mrl = mrl_vals_sergio), aes(u, mrl)) +
  geom_line() + geom_point() +
  geom_vline(xintercept = u0_sergio, linetype = "dashed", color = "red") +
  labs(title = "MRL plot (log TL, WOODWARD)",
       x = "Threshold (log TL)", y = "Mean excess") +
  theme_science_polished
p_mrl_sergio
ggsave(file.path(fig_dir, "wood_tl_log_mrl.png"), p_mrl_sergio, dpi = 600, width = 7, height = 5)

shape_sergio <- scale_sergio <- shape_lo_sergio <- shape_hi_sergio <- adj_sergio <- adj_lo_sergio <- adj_hi_sergio <- rep(NA_real_, length(u_seq_sergio))
for (i in seq_along(u_seq_sergio)) {
  u_sergio <- u_seq_sergio[i]
  exc_sergio <- y_sergio[y_sergio > u_sergio] - u_sergio
  if (length(exc_sergio) < 20) next
  out_sergio <- try(fit_gpd_at_u(y_sergio, u_sergio), silent = TRUE)
  if (inherits(out_sergio, "try-error")) next
  scale_sergio[i] <- out_sergio$scale; shape_sergio[i] <- out_sergio$shape
  zcrit     <- qnorm(1 - (1 - ci_level)/2)
  se_scale_sergio  <- sqrt(out_sergio$cov[1,1]); se_shape_sergio <- sqrt(out_sergio$cov[2,2])
  shape_lo_sergio[i] <- shape_sergio[i] - zcrit * se_shape_sergio
  shape_hi_sergio[i] <- shape_sergio[i] + zcrit * se_shape_sergio
  adj_sergio[i] <- adj_scale(scale_sergio[i], shape_sergio[i], u_sergio, u0_sergio)
  var_adj_sergio <- out_sergio$cov[1,1] + (u_sergio - u0_sergio)^2 * out_sergio$cov[2,2] - 2*(u_sergio - u0_sergio)*out_sergio$cov[1,2]
  se_adj_sergio  <- sqrt(max(var_adj_sergio, 0))
  adj_lo_sergio[i] <- adj_sergio[i] - zcrit * se_adj_sergio
  adj_hi_sergio[i] <- adj_sergio[i] + zcrit * se_adj_sergio
}

p_shape_sergio <- ggplot(data.frame(u = u_seq_sergio, shape_sergio, shape_lo_sergio, shape_hi_sergio), aes(u, shape_sergio)) +
  geom_point() +
  geom_errorbar(aes(ymin = shape_lo_sergio, ymax = shape_hi_sergio), width = 0.03, color = "blue") +
  geom_vline(xintercept = u0_sergio, color = "red", linetype = "dashed") +
  labs(x = "Threshold (log TL)", y = "Shape (xi)") +
  theme_science_polished + 
  geom_hline(yintercept = shape0_sergio, color = "red", linetype = "dashed")
p_shape_sergio
ggsave(file.path(fig_dir, "shape_stability_sergio.png"), p_shape_sergio, dpi = 600, width = 7, height = 5)


p_adj_sergio <- ggplot(data.frame(u = u_seq_sergio, adj_sergio, adj_lo_sergio, adj_hi_sergio), aes(u, adj_sergio)) +
  geom_point() +
  geom_errorbar(aes(ymin = adj_lo_sergio, ymax = adj_hi_sergio), width = 0.03, color = "blue") +
  geom_vline(xintercept = u0_sergio, color = "red", linetype = "dashed") +
  labs(x = "Threshold (log TL)", y = "Adjusted scale") +
  theme_science_polished + 
  geom_hline(yintercept = adj0_sergio, color = "red", linetype = "dashed")
p_adj_sergio
ggsave(file.path(fig_dir, "adj_scale_stability_sergio.png"), p_adj_sergio, dpi = 600, width = 7, height = 5)



## ---------------------------------------------------------
## 6) Tail diagnostics
## ---------------------------------------------------------
# Woodward
diags_allan <- diagnostic_plots(y_allan, u0_allan, fit0_allan$scale, fit0_allan$shape, "log TL (WOODWARD)")
ksout_allan <- ks_boot_pvalue(y_allan, u0_allan, fit0_allan$scale, fit0_allan$shape, B = n_boot_ks)

ggsave(file.path(fig_dir, "allan_tl_log_qq.png"),  diags_allan$pqq, dpi = 600, width = 7, height = 5)
ggsave(file.path(fig_dir, "allan_tl_log_pp.png"),  diags_allan$ppp, dpi = 600, width = 7, height = 5)
ggsave(file.path(fig_dir, "allan_tl_log_pit_uniformity.png"), diags_allan$pks, dpi = 600, width = 7, height = 5)

message(sprintf("[WOODWARD log TL] u0_allan=%.3f | scale=%.4f shape=%.4f | KS p=%.3f",
                u0_allan, fit0_allan$scale, fit0_allan$shape, ksout_allan$p_value))

# Balaguera-Reina
diags_sergio <- diagnostic_plots(y_sergio, u0_sergio, fit0_sergio$scale, fit0_sergio$shape, "log TL (WOODWARD)")
ksout_sergio <- ks_boot_pvalue(y_sergio, u0_sergio, fit0_sergio$scale, fit0_sergio$shape, B = n_boot_ks)

ggsave(file.path(fig_dir, "sergio_tl_log_qq.png"),  diags_sergio$pqq, dpi = 600, width = 7, height = 5)
ggsave(file.path(fig_dir, "sergio_tl_log_pp.png"),  diags_sergio$ppp, dpi = 600, width = 7, height = 5)
ggsave(file.path(fig_dir, "sergio_tl_log_pit_uniformity.png"), diags_sergio$pks, dpi = 600, width = 7, height = 5)

message(sprintf("[BR log TL] u0_sergio=%.3f | scale=%.4f shape=%.4f | KS p=%.3f",
                u0_sergio, fit0_sergio$scale, fit0_sergio$shape, ksout_sergio$p_value))

diags_sergio

## ---------------------------------------------------------
## 7) Endpoint posterior (WOODWARD)
## ---------------------------------------------------------
y_allan_above <- y_allan[y_allan > u0_allan]; n_ex <- length(y_allan_above)
if (n_ex < 10) stop("Not enough exceedances for endpoint posterior.")

y_allanmax <- max(y_allan_above)
eps  <- max(1e-6, 1e-4 * abs(y_allanmax))
y_allan_star_grid <- seq(y_allanmax + eps, y_allanmax + 10, length.out = 2000)   # y_allan* grid (log TL)
xi_grid     <- seq(-1.0, -0.02, length.out = 1000)             # xi < 0

loglik_endpoint <- function(y_star, xi, u, y_ex) {
  if (xi >= 0 || any(y_star <= y_ex)) return(-Inf)
  n <- length(y_ex)
  term1 <- n * log(y_star - u) / xi
  term2 <- -(1/xi + 1) * sum(log(y_star - y_ex))
  term3 <- -n * log(-xi)
  term1 + term2 + term3
}

# Woodward prior on TL* (center 430 cm on log scale)
mu_y_allan <- log(430)
prior_mode <- "woodward_loose"  # "woodward_mod" or "woodward_tight"
sd_y_allan <- switch(prior_mode,
               "woodward_loose" = 0.095,
               "woodward_mod"   = 0.048,
               "woodward_tight" = 0.024, 0.048)
logprior_y_allan  <- function(y_allan_star) dnorm(y_allan_star, mean = mu_y_allan, sd = sd_y_allan, log = TRUE)

# PC prior on kappa=-xi ~ Exp(lambda), xi in (-1,0)
lambda_k   <- 3.0
logprior_xi <- function(xi) {
  if (xi >= 0 || xi <= -1) return(-Inf)
  kappa <- -xi
  if (kappa <= 0 || kappa >= 1) return(-Inf)
  log(lambda_k) - lambda_k * kappa
}

marg_log_post <- sapply(y_allan_star_grid, function(y_allans) {
  lv <- sapply(xi_grid, function(xi) loglik_endpoint(y_allans, xi, u0_allan, y_allan_above) + logprior_xi(xi))
  lv_max <- max(lv)
  log_int_xi <- log(pracma::trapz(xi_grid, exp(lv - lv_max))) + lv_max
  log_int_xi + logprior_y_allan(y_allans)
})

post_unnorm <- exp(marg_log_post - max(marg_log_post))
Z_y_allan         <- pracma::trapz(y_allan_star_grid, post_unnorm)
if (!is.finite(Z_y_allan) || Z_y_allan <= 0) stop("Normalization failed for p(y_allan*|y_allan).")
post <- post_unnorm / Z_y_allan

cdf_from_density <- function(x, f) {
  f[!is.finite(f)] <- 0
  F <- pracma::cumtrapz(x, f)
  rng <- max(F, na.rm = TRUE)
  if (!is.finite(rng) || rng <= 0) stop("cdf_from_density: non-positive integral.")
  F / rng
}
q_from_cdf <- function(x, f, p) {
  F <- cdf_from_density(x, f)
  as.numeric(approx(x = F, y = x, xout = p, ties = "ordered", rule = 2)$y_allan)
}

y_allanstar_map   <- y_allan_star_grid[which.max(post)]
y_allanstar_upper <- q_from_cdf(y_allan_star_grid, post, ci_level)

tl_grid_cm   <- exp(y_allan_star_grid)
dens_tl      <- post / tl_grid_cm
tl_map_cm    <- exp(y_allanstar_map)
tl_upper_cm  <- exp(y_allanstar_upper)

df_tlpost <- data.frame(TL = tl_grid_cm, dens = dens_tl)
p_tlstar <- ggplot(df_tlpost, aes(TL, dens)) +
  geom_line(linewidth = 1.2, color = "darkblue") +
  geom_vline(xintercept = tl_map_cm,    color = "purple",   linetype = "dashed",  linewidth = 1.1) +
  geom_vline(xintercept = tl_upper_cm,  color = "orange",   linetype = "dotdash", linewidth = 1.0) +
  geom_vline(xintercept = stokes_tl_cm, color = "firebrick", linetype = "dashed", linewidth = 1.0) +
  annotate("text", x = stokes_tl_cm, y = max(dens_tl)*0.08, vjust = -0.8,
           label = "Stokes (450 cm)", color = "firebrick", angle = 90, size = 4) +
  labs(x = "TL* (cm)", y = "Density") +
  xlim(350, 550) +
  theme_science_polished
ggsave(file.path(fig_dir, "wood_tl_endpoint_posterior.png"),
       p_tlstar, dpi = 600, width = 7, height = 5)

cat("\n--- TL* (cm) via numeric integration on WOODWARD ---\n",
    "MAP:          ", round(tl_map_cm, 2), "\n",
    sprintf("%d%% upper:  ", round(ci_level*100)), round(tl_upper_cm, 2), "\n")

## ---------------------------------------------------------
## 8) Mass extrapolation — regression fitted on RENEWED CSV
## ---------------------------------------------------------
y_allanstar_draws <- sample(y_allan_star_grid, size = n_post_y_allanstar, replace = TRUE, prob = post)

w_fit_df <- dat %>%
  filter(Deform == 0L, is.finite(WTkg), is.finite(TL)) %>%
  transmute(logW = log(WTkg), logTL = log(TL))
if (nrow(w_fit_df) < 10) {
  w_fit_df <- dat %>% filter(is.finite(WTkg), is.finite(TL)) %>%
    transmute(logW = log(WTkg), logTL = log(TL))
}
stopifnot(nrow(w_fit_df) >= 10)

lm_w_tl <- lm(logW ~ logTL, data = w_fit_df)
beta    <- coef(lm_w_tl)
Vcoef   <- vcov(lm_w_tl)
sigma_e <- summary(lm_w_tl)$sigma

ab_samps <- MASS::mvrnorm(n = length(y_allanstar_draws), mu = beta, Sigma = Vcoef)
a_s <- ab_samps[, 1]; b_s <- ab_samps[, 2]
logM_draws <- a_s + b_s * y_allanstar_draws

q_indiv <- 0.95
logM_draws <- logM_draws + qnorm(q_indiv) * sigma_e

mass_draws_kg <- exp(logM_draws)
mass_up95_kg  <- as.numeric(quantile(mass_draws_kg, probs = ci_level))
dens_m <- density(mass_draws_kg, n = 4096, adjust = 1.0)
mass_map_kg <- dens_m$x[which.max(dens_m$y_allan)]

p_mass <- ggplot(data.frame(W = mass_draws_kg), aes(W)) +
  geom_density(linewidth = 1.2, color = "darkblue") +
  geom_vline(xintercept = mass_map_kg,  color = "purple",  linetype = "dashed",  linewidth = 1) +
  geom_vline(xintercept = mass_up95_kg, color = "orange",  linetype = "dotdash", linewidth = 1) +
  geom_vline(xintercept = stokes_w_kg,  color = "firebrick", linetype = "dashed", linewidth = 1) +
  annotate("text", x = stokes_w_kg, y = max(dens_m$y_allan)*0.08, vjust = -0.8,
           label = "Stokes (459 kg)", color = "firebrick", angle = 90, size = 4) +
  labs(x = "Weight (kg)", y = "Density") +
  theme_science_polished
p_mass
ggsave(file.path(fig_dir, "mass_endpoint_from_csv_regression.png"),
       p_mass, dpi = 600, width = 7, height = 5)

cat(sprintf("\nMass endpoint (typical individual, q=%d%%) using CSV regression:\n", round(q_indiv*100)))
cat(sprintf("MAP: %.2f kg | %d%% upper: %.2f kg\n",
            mass_map_kg, round(ci_level*100), mass_up95_kg))

## ---------------------------------------------------------
## 8b) Regression panel — scatter from CSV, endpoint from WOODWARD TL*
## ---------------------------------------------------------
plot_df <- dat %>%
  filter(is.finite(WTkg), is.finite(TL)) %>%
  transmute(logTL = log(TL),
            logW  = log(WTkg),
            Sex   = factor(ifelse(is.na(Sex), "NA", Sex), levels = c("F","M","NA")))
sex_cols <- c("F" = "#ff6fb3", "M" = "#66b3ff", "NA" = "gray70")

y_at_u0_allan <- as.numeric(predict(lm_w_tl, newdata = data.frame(logTL = u0_allan)))
x_ep    <- y_allanstar_map
x_ep_ub <- y_allanstar_upper
y_allan_ep    <- as.numeric(predict(lm_w_tl, data.frame(logTL = x_ep)))
y_allan_ep_ub <- y_allan_ep + qnorm(0.95) * sigma_e

x_stokes <- log(stokes_tl_cm)
y_stokes <- log(stokes_w_kg)
ab <- coef(lm_w_tl)

p_reg <- ggplot(plot_df, aes(x = logTL, y = logW, color = Sex)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_abline(intercept = ab[1], slope = ab[2], linewidth = 1) +
  geom_vline(xintercept = u0_allan, linetype = "dashed", color = "gray30") +
  geom_hline(yintercept = y_at_u0_allan, linetype = "dashed", color = "gray30") +
  geom_point(aes(x = x_ep, y = y_allan_ep), inherit.aes = FALSE, color = "darkgreen", size = 3) +
  geom_segment(aes(x = x_ep, y = y_allan_ep, xend = x_ep_ub, yend = y_allan_ep),
               inherit.aes = FALSE, color = "darkgreen", linewidth = 1.1) +
  geom_segment(aes(x = x_ep, y = y_allan_ep, xend = x_ep, yend = y_allan_ep_ub),
               inherit.aes = FALSE, color = "darkgreen", linewidth = 1.1) +
  geom_point(aes(x = x_stokes, y = y_stokes),
             inherit.aes = FALSE, color = "firebrick", size = 2.8) +
  geom_text(aes(x = x_stokes, y = y_stokes, label = "Stokes"),
            inherit.aes = FALSE, nudge_y = 0.015, color = "firebrick", size = 3.5) +
  scale_color_manual(values = sex_cols, breaks = c("F","M","NA"), drop = FALSE) +
  labs(x = "log Total Length (CSV)", y = "log Weight (CSV)", color = "Sex") +
  theme_science_polished
ggsave(file.path(fig_dir, "reg_logW_vs_logTL_threshold_endpoint_stokes_CSV.png"),
       p_reg, dpi = 600, width = 8, height = 5.2)

## ---------------------------------------------------------
## 9) Stokes exceedance using WOODWARD posterior
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

stopifnot(exists("y_allan"), exists("u0_allan"), exists("y_allan_above"),
          exists("y_allan_star_grid"), exists("xi_grid"),
          exists("loglik_endpoint"), exists("logprior_xi"),
          exists("post"))

W_y_allan <- trapz_weights(y_allan_star_grid, post)
stopifnot(all(is.finite(W_y_allan)), abs(sum(W_y_allan) - 1) < 1e-6)

p_exceed_u <- mean(y_allan > u0_allan, na.rm = TRUE)
y_claim_cm  <- stokes_tl_cm
y_claim_log <- log(y_claim_cm)

Si_vec <- numeric(length(y_allan_star_grid))
for (i in seq_along(y_allan_star_grid)) {
  y_allans <- y_allan_star_grid[i]
  if (y_allans <= y_claim_log) { Si_vec[i] <- 0; next }
  lv <- sapply(xi_grid, function(xi) {
    loglik_endpoint(y_allans, xi, u0_allan, y_allan_above) + logprior_xi(xi)
  })
  m  <- max(lv)
  fu <- exp(lv - m)
  V_xi <- trapz_weights(xi_grid, fu)
  surv_term <- ((y_allans - y_claim_log) / (y_allans - u0_allan))^(-1/xi_grid)
  surv_term[!is.finite(surv_term)] <- 0
  Si_vec[i] <- sum(V_xi * surv_term)
}
p_tail_cond <- sum(W_y_allan * Si_vec)
p_exceed_claim <- p_exceed_u * p_tail_cond

cat(sprintf("\n--- Exceedance probability for Stokes alligator (WOODWARD posterior) ---\n"))
cat(sprintf("Threshold u0_allan (log cm)  : %.4f  (≈ %.1f cm)\n", u0_allan, exp(u0_allan)))
cat(sprintf("Claim y (log cm)       : %.4f  (450.0 cm)\n", y_claim_log))
cat(sprintf("P(Y > u0_allan | data)       : %.4f\n", p_exceed_u))
cat(sprintf("P(Y > y | Y>u0_allan, data)  : %.4f\n", p_tail_cond))
cat(sprintf("P(Y > y | data)        : %.6f\n", p_exceed_claim))

