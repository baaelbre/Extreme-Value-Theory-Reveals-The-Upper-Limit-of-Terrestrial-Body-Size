#### =========================
#### EVT on log Total Length (FEMALES) with point imputation from log SVL → log TL
#### Analysis on WOODWARD (F); mass mapping from renewed CSV (F)
#### Female prior + fast diagnostics
#### =========================

## ---------------------------
## 0) Libraries & global setup
## ---------------------------
suppressPackageStartupMessages({
  library(ggplot2)
  library(extRemes)    # fevd/revd/pevd for GP fits & sims
  library(dplyr)
  library(readxl)
  library(pracma)      # trapz for numerical integration
  library(HDInterval)  # hdi
  library(MASS)        # mvrnorm
  library(grid)
})

set.seed(42)

fig_dir <- "Figures"
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

## ---------------------------
## 1) Settings & theme
## ---------------------------
ci_level             <- 0.95
n_boot_ks            <- 1000
n_post_ystar         <- 20000
thresh_q_opt_allan   <- 0.90
thresh_q_opt_sergio  <- 0.985

# Female-relevant reference values (set to NA_real_ to hide on plots)
ref_tl_cm_female <- 360
ref_w_kg_female  <- NA_real_

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
## 2) Data ingest (FEMALES ONLY)
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
  filter(Sex == "F") %>%   # FEMALES ONLY
  mutate(
    SVL  = ifelse(is.finite(SVL)  & SVL  > 0, SVL,  NA_real_),
    TL   = ifelse(is.finite(TL)   & TL   > 0, TL,   NA_real_),
    WTkg = ifelse(is.finite(WTkg) & WTkg > 0, WTkg, NA_real_)
  )

## 2b) WOODWARD dataset — used for EVT (TL) analysis
wood_path <- "Data/experimental_alligator_harvest_woodward.xlsx"
wood_raw  <- readxl::read_excel(wood_path)

# Normalize Sex in Woodward, filter females
wood <- wood_raw %>%
  mutate(
    Sex   = tolower(trimws(as.character(Sex))),
    Sex   = case_when(
      Sex %in% c("f","female","fem","vrouw","f?","female?") ~ "f",
      TRUE ~ Sex
    ),
    SVL   = suppressWarnings(as.numeric(SVL)),
    TL    = suppressWarnings(as.numeric(TL)),
    WTkg  = suppressWarnings(as.numeric(if ("WTkg" %in% names(.)) WTkg else NA_real_)),
    Deform = case_when(is.na(Deform) ~ 0L, TRUE ~ as.integer(Deform))
  ) %>%
  filter(Sex == "f") %>%   # FEMALES ONLY
  mutate(
    TL   = ifelse(is.finite(TL)   & TL   > 0, TL,   NA_real_),
    WTkg = ifelse(is.finite(WTkg) & WTkg > 0, WTkg, NA_real_),
    SVL  = ifelse(is.finite(SVL)  & SVL  > 0, SVL,  NA_real_)
  )

## ---------------------------
## 2c) Print & count TL/W cutoffs (FEMALES)
## ---------------------------
cut_cm <- 327
cut_kg <- 177

cat("\n========== FEMALES: TL ≥ ", cut_cm, " cm (WOODWARD) ==========\n", sep = "")
print(wood %>% filter(is.finite(TL), TL >= cut_cm))
cat(sprintf("Count (WOODWARD-F, TL ≥ %d cm): %d\n", cut_cm, sum(is.finite(wood$TL) & wood$TL >= cut_cm)))

cat("\n========== FEMALES: TL ≥ ", cut_cm, " cm (CSV) ==========\n", sep = "")
print(dat %>% filter(is.finite(TL), TL >= cut_cm))
cat(sprintf("Count (CSV-F, TL ≥ %d cm): %d\n\n", cut_cm, sum(is.finite(dat$TL) & dat$TL >= cut_cm)))

cat("\n========== FEMALES: W ≥ ", cut_kg, " kg (WOODWARD) ==========\n", sep = "")
print(wood %>% filter(is.finite(WTkg), WTkg >= cut_kg))
cat(sprintf("Count (WOODWARD-F, W ≥ %d kg): %d\n", cut_kg, sum(is.finite(wood$WTkg) & wood$WTkg >= cut_kg)))

cat("\n========== FEMALES: W ≥ ", cut_kg, " kg (CSV) ==========\n", sep = "")
print(dat %>% filter(is.finite(WTkg), WTkg >= cut_kg))
cat(sprintf("Count (CSV-F, W ≥ %d kg): %d\n\n", cut_kg, sum(is.finite(dat$WTkg) & dat$WTkg >= cut_kg)))

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
  
  F_hat <- pevd(exc, scale = scale_hat, shape = shape_hat, type="GP")
  dfpp <- data.frame(Theoretical = sort(F_hat), Empirical = (1:n)/n)
  ppp <- ggplot(dfpp, aes(Theoretical, Empirical)) +
    geom_point(color = "darkgreen") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    labs(x = "Theoretical CDF", y = "Empirical CDF") +
    theme_science_polished
  
  pks <- ggplot(data.frame(F_hat = F_hat), aes(F_hat)) +
    geom_histogram(aes(y = ..density..), bins = 20,
                   fill = "skyblue", color = "black", alpha = 0.7) +
    geom_hline(yintercept = 1, color = "red", linetype = "dashed") +
    geom_density(color = "darkblue", linewidth = 1.1, adjust = 1.5) +
    labs(title = paste0("PIT Uniformity — ", label), x = expression(hat(F)(y)), y = "Density") +
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
  list(ks_obs = as.numeric(ks_obs),
       p_value = if (length(ks_b)) mean(ks_b >= ks_obs) else NA_real_,
       ks_boot = ks_b)
}

## ---------------------------------------------------------
## 4) Impute missing TL via SVL→TL (WOODWARD-F, used for EVT)
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

# Tail variables (log TL)
y_allan  <- wood$logTL_use; y_allan  <- y_allan[is.finite(y_allan)]
y_sergio <- log(dat$TL);     y_sergio <- y_sergio[is.finite(y_sergio)]
stopifnot(length(y_allan)  >= 50, length(y_sergio) >= 50)

## ---------------------------------------------------------
## 5) Threshold selection
## ---------------------------------------------------------
# Woodward (F)
u_seq_allan <- make_thresholds(y_allan, q_from = 0.5, q_to = 0.97, n = 50)
u0_allan    <- as.numeric(quantile(y_allan, thresh_q_opt_allan, na.rm = TRUE))
fit0_allan  <- fit_gpd_at_u(y_allan, u0_allan)
shape0_allan <- fit0_allan$shape
adj0_allan   <- fit0_allan$scale
shape0_se_allan <- sqrt(fit0_allan$cov[2,2])
z_allan         <- shape0_allan / shape0_se_allan
p_wald_allan    <- pnorm(z_allan)  # one-sided H1: xi<0
p_wald_allan

# CSV (F)
u_seq_sergio <- make_thresholds(y_sergio, q_from = 0.7, q_to = 0.99, n = 50)
u0_sergio    <- as.numeric(quantile(y_sergio, thresh_q_opt_sergio, na.rm = TRUE))
fit0_sergio  <- fit_gpd_at_u(y_sergio, u0_sergio)
shape0_sergio <- fit0_sergio$shape
adj0_sergio   <- fit0_sergio$scale
shape0_se_sergio <- sqrt(fit0_sergio$cov[2,2])
z_sergio         <- shape0_sergio / shape0_se_sergio
p_wald_sergio    <- pnorm(z_sergio)  # one-sided H1: xi<0
p_wald_sergio

# MRL & stability (Woodward-F)
mrl_vals_allan <- mrl_data(y_allan, u_seq_allan)
p_mrl_allan <- ggplot(data.frame(u = u_seq_allan, mrl = mrl_vals_allan), aes(u, mrl)) +
  geom_line() + geom_point() +
  geom_vline(xintercept = u0_allan, linetype = "dashed", color = "red") +
  labs(
       x = "Threshold (log TL)", y = "Mean excess") +
  theme_science_polished
p_mrl_allan
ggsave(file.path(fig_dir, "F_wood_tl_log_mrl.png"), p_mrl_allan, dpi = 600, width = 7, height = 5)

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
  theme_science_polished
p_shape_allan
ggsave(file.path(fig_dir, "f_shape_stability_allan.png"), p_shape_allan, dpi = 600, width = 7, height = 5)

p_adj_allan <- ggplot(data.frame(u = u_seq_allan, adj_allan, adj_lo_allan, adj_hi_allan), aes(u, adj_allan)) +
  geom_point() +
  geom_errorbar(aes(ymin = adj_lo_allan, ymax = adj_hi_allan), width = 0.03, color = "blue") +
  geom_vline(xintercept = u0_allan, color = "red", linetype = "dashed") +
  labs(x = "Threshold (log TL)", y = "Adjusted scale") +
  theme_science_polished
p_adj_allan
ggsave(file.path(fig_dir, "f_adj_scale_stability_allan.png"), p_adj_allan, dpi = 600, width = 7, height = 5)

# MRL & stability (CSV-F) — optional diagnostics
mrl_vals_sergio <- mrl_data(y_sergio, u_seq_sergio)
p_mrl_sergio <- ggplot(data.frame(u = u_seq_sergio, mrl = mrl_vals_sergio), aes(u, mrl)) +
  geom_line() + geom_point() +
  geom_vline(xintercept = u0_sergio, linetype = "dashed", color = "red") +
  labs(
       x = "Threshold (log TL)", y = "Mean excess") +
  theme_science_polished
p_mrl_sergio
ggsave(file.path(fig_dir, "f_sergio_tl_log_mrl.png"), p_mrl_sergio, dpi = 600, width = 7, height = 5)

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
  theme_science_polished
p_shape_sergio
ggsave(file.path(fig_dir, "f_shape_stability_sergio.png"), p_shape_sergio, dpi = 600, width = 7, height = 5)

p_adj_sergio <- ggplot(data.frame(u = u_seq_sergio, adj_sergio, adj_lo_sergio, adj_hi_sergio), aes(u, adj_sergio)) +
  geom_point() +
  geom_errorbar(aes(ymin = adj_lo_sergio, ymax = adj_hi_sergio), width = 0.03, color = "blue") +
  geom_vline(xintercept = u0_sergio, color = "red", linetype = "dashed") +
  labs(x = "Threshold (log TL)", y = "Adjusted scale") +
  theme_science_polished
p_adj_sergio
ggsave(file.path(fig_dir, "f_adj_scale_stability_sergio.png"), p_adj_sergio, dpi = 600, width = 7, height = 5)

## ---------------------------------------------------------
## 6) Tail diagnostics
## ---------------------------------------------------------
diags_allan <- diagnostic_plots(y_allan, u0_allan, fit0_allan$scale, fit0_allan$shape, "log TL (WOODWARD — F)")
ksout_allan <- ks_boot_pvalue(y_allan, u0_allan, fit0_allan$scale, fit0_allan$shape, B = n_boot_ks)
ggsave(file.path(fig_dir, "f_qq_allan.png"),  diags_allan$pqq, dpi = 600, width = 7, height = 5)
ggsave(file.path(fig_dir, "f_pp_allan.png"),  diags_allan$ppp, dpi = 600, width = 7, height = 5)
ggsave(file.path(fig_dir, "F_allan_tl_log_pit_uniformity.png"), diags_allan$pks, dpi = 600, width = 7, height = 5)
message(sprintf("[WOODWARD-F log TL] u0=%.3f | scale=%.4f shape=%.4f | KS p=%.3f",
                u0_allan, fit0_allan$scale, fit0_allan$shape, ksout_allan$p_value))

diags_sergio <- diagnostic_plots(y_sergio, u0_sergio, fit0_sergio$scale, fit0_sergio$shape, "log TL (CSV — F)")
ksout_sergio <- ks_boot_pvalue(y_sergio, u0_sergio, fit0_sergio$scale, fit0_sergio$shape, B = n_boot_ks)
ggsave(file.path(fig_dir, "f_qq_sergio.png"),  diags_sergio$pqq, dpi = 600, width = 7, height = 5)
ggsave(file.path(fig_dir, "f_pp_sergio.png"),  diags_sergio$ppp, dpi = 600, width = 7, height = 5)
ggsave(file.path(fig_dir, "F_sergio_tl_log_pit_uniformity.png"), diags_sergio$pks, dpi = 600, width = 7, height = 5)
message(sprintf("[CSV-F log TL] u0=%.3f | scale=%.4f shape=%.4f | KS p=%.3f",
                u0_sergio, fit0_sergio$scale, fit0_sergio$shape, ksout_sergio$p_value))

## ---------------------------------------------------------
## 7) Endpoint posterior (WOODWARD-F)
## ---------------------------------------------------------
y_allan_above <- y_allan[y_allan > u0_allan]; n_ex <- length(y_allan_above)
if (n_ex < 10) stop("Not enough exceedances for endpoint posterior (F).")

ymax <- max(y_allan_above)
eps  <- max(1e-6, 1e-4 * abs(ymax))
y_star_grid <- seq(ymax + eps, ymax + 10, length.out = 2000)  # log TL*
xi_grid     <- seq(-1.0, -0.02, length.out = 1000)            # xi < 0

loglik_endpoint <- function(y_star, xi, u, y_ex) {
  if (xi >= 0 || any(y_star <= y_ex)) return(-Inf)
  n <- length(y_ex)
  n * log(y_star - u) / xi - (1/xi + 1) * sum(log(y_star - y_ex)) - n * log(-xi)
}

# Female prior for TL*: center around ~360 cm (log scale)
mu_y_allan <- log(360)
prior_mode <- "mod"  # "loose" / "mod" / "tight"
sd_y_allan <- switch(prior_mode,
                     "loose" = 0.095,
                     "mod"   = 0.048,
                     "tight" = 0.024, 0.048)
logprior_y <- function(y_star) dnorm(y_star, mean = mu_y_allan, sd = sd_y_allan, log = TRUE)

# PC prior on kappa=-xi ~ Exp(lambda), xi in (-1,0)
lambda_k   <- 3.0
logprior_xi <- function(xi) {
  if (xi >= 0 || xi <= -1) return(-Inf)
  kappa <- -xi
  if (kappa <= 0 || kappa >= 1) return(-Inf)
  log(lambda_k) - lambda_k * kappa
}

marg_log_post <- sapply(y_star_grid, function(ys) {
  lv <- sapply(xi_grid, function(xi) loglik_endpoint(ys, xi, u0_allan, y_allan_above) + logprior_xi(xi))
  lv_max <- max(lv)
  log_int_xi <- log(pracma::trapz(xi_grid, exp(lv - lv_max))) + lv_max
  log_int_xi + logprior_y(ys)
})

post_unnorm <- exp(marg_log_post - max(marg_log_post))
Z_norm      <- pracma::trapz(y_star_grid, post_unnorm)
if (!is.finite(Z_norm) || Z_norm <= 0) stop("Normalization failed for p(y*|data) (F).")
post <- post_unnorm / Z_norm

cdf_from_density <- function(x, f) {
  f[!is.finite(f)] <- 0
  F <- pracma::cumtrapz(x, f)
  F / max(F, na.rm = TRUE)
}
q_from_cdf <- function(x, f, p) {
  F <- cdf_from_density(x, f)
  as.numeric(approx(x = F, y = x, xout = p, ties = "ordered", rule = 2)$y)
}

y_star_map   <- y_star_grid[which.max(post)]
y_star_upper <- q_from_cdf(y_star_grid, post, ci_level)

tl_grid_cm   <- exp(y_star_grid)
dens_tl      <- post / tl_grid_cm
tl_map_cm    <- exp(y_star_map)
tl_upper_cm  <- exp(y_star_upper)

df_tlpost <- data.frame(TL = tl_grid_cm, dens = dens_tl)
p_tlstar <- ggplot(df_tlpost, aes(TL, dens)) +
  geom_line(linewidth = 1.2, color = "darkblue") +
  geom_vline(xintercept = tl_map_cm,    color = "purple",   linetype = "dashed",  linewidth = 1.1) +
  geom_vline(xintercept = tl_upper_cm,  color = "orange",   linetype = "dotdash", linewidth = 1.0) +
  { if (is.finite(ref_tl_cm_female)) geom_vline(xintercept = ref_tl_cm_female, color = "firebrick", linetype = "dashed", linewidth = 1.0) } +
  labs(x = "TL* (cm)", y = "Density", title = "Endpoint posterior — WOODWARD (F)") +
  coord_cartesian(xlim = c(280, 500)) +
  theme_science_polished
ggsave(file.path(fig_dir, "F_wood_tl_endpoint_posterior.png"),
       p_tlstar, dpi = 600, width = 7, height = 5)

cat("\n--- FEMALES: TL* (cm) via numeric integration on WOODWARD ---\n",
    "MAP:          ", round(tl_map_cm, 2), "\n",
    sprintf("%d%% upper:  ", round(ci_level*100)), round(tl_upper_cm, 2), "\n", sep = "")

## ---------------------------------------------------------
## 8) Mass extrapolation — regression fitted on CSV (F)
## ---------------------------------------------------------
y_star_draws <- sample(y_star_grid, size = n_post_ystar, replace = TRUE, prob = post)

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

ab_samps <- MASS::mvrnorm(n = length(y_star_draws), mu = beta, Sigma = Vcoef)
a_s <- ab_samps[, 1]; b_s <- ab_samps[, 2]
logM_draws <- a_s + b_s * y_star_draws

q_indiv <- 0.95
logM_draws <- logM_draws + qnorm(q_indiv) * sigma_e

mass_draws_kg <- exp(logM_draws)
mass_up95_kg  <- as.numeric(quantile(mass_draws_kg, probs = ci_level))
dens_m <- density(mass_draws_kg, n = 4096, adjust = 1.0)
mass_map_kg <- dens_m$x[ which.max(dens_m$y) ]

p_mass <- ggplot(data.frame(W = mass_draws_kg), aes(W)) +
  geom_density(linewidth = 1.2, color = "darkblue") +
  geom_vline(xintercept = mass_map_kg,  color = "purple",  linetype = "dashed",  linewidth = 1) +
  geom_vline(xintercept = mass_up95_kg, color = "orange",  linetype = "dotdash", linewidth = 1) +
  { if (is.finite(ref_w_kg_female)) geom_vline(xintercept = ref_w_kg_female, color = "firebrick", linetype = "dashed", linewidth = 1) } +
  labs(x = "Weight (kg)", y = "Density", title = "Mass endpoint from CSV regression (F)") +
  theme_science_polished
ggsave(file.path(fig_dir, "F_mass_endpoint_from_csv_regression.png"),
       p_mass, dpi = 600, width = 7, height = 5)

cat(sprintf("\nFEMALES — Mass endpoint (typical individual, q=%d%%) using CSV regression:\n", round(q_indiv*100)))
cat(sprintf("MAP: %.1f kg | %d%% upper: %.1f kg\n", mass_map_kg, round(ci_level*100), mass_up95_kg))

## ---------------------------------------------------------
## 8b) Regression panel — scatter from CSV (F), endpoint from WOODWARD TL*
## ---------------------------------------------------------
plot_df <- dat %>% filter(is.finite(WTkg), is.finite(TL)) %>% transmute(logTL = log(TL), logW = log(WTkg))
y_at_u0 <- as.numeric(predict(lm_w_tl, newdata = data.frame(logTL = u0_allan)))
x_ep    <- y_star_map
x_ep_ub <- y_star_upper
y_ep    <- as.numeric(predict(lm_w_tl, data.frame(logTL = x_ep)))
y_ep_ub <- y_ep + qnorm(0.95) * sigma_e

p_reg <- ggplot(plot_df, aes(x = logTL, y = logW)) +
  geom_point(alpha = 0.7, size = 2, color = "#ff6fb3") +
  geom_abline(intercept = beta[1], slope = beta[2], linewidth = 1) +
  geom_vline(xintercept = u0_allan, linetype = "dashed", color = "gray30") +
  geom_hline(yintercept = y_at_u0,   linetype = "dashed", color = "gray30") +
  geom_point(aes(x = x_ep, y = y_ep), inherit.aes = FALSE, color = "darkgreen", size = 3) +
  geom_segment(aes(x = x_ep, y = y_ep, xend = x_ep_ub, yend = y_ep),
               inherit.aes = FALSE, color = "darkgreen", linewidth = 1.1) +
  geom_segment(aes(x = x_ep, y = y_ep, xend = x_ep, yend = y_ep_ub),
               inherit.aes = FALSE, color = "darkgreen", linewidth = 1.1) +
  labs(x = "log Total Length (CSV — F)", y = "log Weight (CSV — F)",
       title = "Female CSV: logW vs logTL with threshold & endpoint") +
  theme_science_polished
ggsave(file.path(fig_dir, "F_reg_logW_vs_logTL_threshold_endpoint_CSV.png"),
       p_reg, dpi = 600, width = 8, height = 5.2)

## ---------------------------------------------------------
## 9) Exceedance at a female TL claim using WOODWARD-F posterior
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
p_exceed_u <- mean(y_allan > u0_allan, na.rm = TRUE)

# Choose claim (female-relevant). Default 360 cm; change as needed.
y_claim_cm  <- if (is.finite(ref_tl_cm_female)) ref_tl_cm_female else 360
y_claim_log <- log(y_claim_cm)

Si_vec <- numeric(length(y_star_grid))
for (i in seq_along(y_star_grid)) {
  ys <- y_star_grid[i]
  if (ys <= y_claim_log) { Si_vec[i] <- 0; next }
  lv <- sapply(xi_grid, function(xi) loglik_endpoint(ys, xi, u0_allan, y_allan_above) + logprior_xi(xi))
  m  <- max(lv)
  fu <- exp(lv - m)
  V_xi <- trapz_weights(xi_grid, fu)
  surv_term <- ((ys - y_claim_log) / (ys - u0_allan))^(-1/xi_grid)
  surv_term[!is.finite(surv_term)] <- 0
  Si_vec[i] <- sum(V_xi * surv_term)
}
p_tail_cond    <- sum(W_y * Si_vec)
p_exceed_claim <- p_exceed_u * p_tail_cond

cat(sprintf("\n--- FEMALES: Exceedance probability for TL = %.0f cm (WOODWARD posterior) ---\n", y_claim_cm))
cat(sprintf("Threshold u0 (log cm)  : %.4f  (≈ %.1f cm)\n", u0_allan, exp(u0_allan)))
cat(sprintf("Claim y (log cm)       : %.4f  (%.1f cm)\n", y_claim_log, y_claim_cm))
cat(sprintf("P(Y > u0 | data)       : %.4f\n", p_exceed_u))
cat(sprintf("P(Y > y | Y>u0, data)  : %.4f\n", p_tail_cond))
cat(sprintf("P(Y > y | data)        : %.6f\n", p_exceed_claim))
