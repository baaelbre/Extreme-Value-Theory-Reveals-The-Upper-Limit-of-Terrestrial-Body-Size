#### =========================
#### Alligator EVA (SVL, TL, Weight) — POT/GPD
#### =========================

### LIBRARIES ###
library(ggplot2)
library(extRemes)   # fevd (GPD)
library(dplyr)
library(readxl)
library(pracma)     # trapz
library(HDInterval) # hdi

### SETTINGS ###
ci_level   <- 0.95
n_boot_ks  <- 1000
fig_dir    <- "Figures"     


theme_science_polished <- theme_minimal(base_family = "Arial", base_size = 12) +
  theme(
    axis.title = element_text(size = 14, face = "bold"),
    axis.text  = element_text(size = 12),
    legend.title = element_text(size = 10, face = "bold"),
    legend.text  = element_text(size = 10),
    legend.key.size = unit(0.5, "cm"),
    panel.grid.major = element_line(color = "gray80", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    axis.ticks.length = unit(0.20, "cm"),
    axis.ticks = element_line(color = "black", linewidth = 0.4),
    plot.margin = margin(5, 5, 5, 5),
    legend.position = "right"
  )

### LOAD DATA ###
data_path <- "Data/experimental_alligator_harvest_woodward.xlsx"
dat <- read_excel(data_path)

# Column names (adapt if your sheet uses different names)
svl_col    <- "SVL"    # cm
tl_col     <- "TL"     # cm
deform_col <- "Deform" # 0/1; 1 = deformed tail
weight_col <- "WTkg"

# Basic cleaning (keep finite)
dat <- dat %>%
  mutate(
    SVL = as.numeric(.data[[svl_col]]),
    TL  = suppressWarnings(as.numeric(.data[[tl_col]])),
    WKG = suppressWarnings(as.numeric(.data[[weight_col]])),
    Deform = if (deform_col %in% names(.)) as.integer(.data[[deform_col]]) else 0L
  )

# --- Small utilities ---

# Robust GPD CDF / RNG
pgpd_ <- function(x, sigma, xi) {
  ifelse(abs(xi) > 1e-10, 1 - (1 + xi * x / sigma)^(-1/xi), 1 - exp(-x / sigma))
}
rgpd_ <- function(n, sigma, xi) {
  u <- runif(n)
  if (abs(xi) > 1e-10) sigma * ((1 - u)^(-xi) - 1) / xi else -sigma * log(1 - u)
}

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

# Mean Residual Life across u
mrl_data <- function(y, u_seq) {
  sapply(u_seq, function(u) {
    exc <- y[y > u] - u
    if (length(exc) < 5) return(NA_real_)
    mean(exc)
  })
}

make_thresholds <- function(y, q_from = 0.75, q_to = 0.95, n = 50) {
  rng <- quantile(y, c(q_from, q_to), na.rm = TRUE)
  seq(rng[1], rng[2], length.out = n)
}

adj_scale <- function(scale, shape, u, u0) scale - shape * (u - u0)

diagnostic_plots <- function(y, u, scale_hat, shape_hat, var_label) {
  exc <- y[y > u] - u
  n   <- length(exc)
  # Q–Q
  probs <- ppoints(n)
  theo_q <- if (abs(shape_hat) > 1e-10) {
    u + scale_hat/shape_hat * ((probs^(-shape_hat)) - 1)
  } else {
    u - scale_hat * log(probs)
  }
  dfqq <- data.frame(Theoretical = rev(theo_q), Empirical = sort(y[y > u]))
  pqq <- ggplot(dfqq, aes(Theoretical, Empirical)) +
    geom_point(color = "steelblue") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    labs(title = paste0("Q–Q (", var_label, ")"), x = "Theoretical", y = "Empirical") +
    theme_science_polished
  
  # P–P
  F_theo <- if (abs(shape_hat) > 1e-10) {
    1 - (1 + shape_hat * exc / scale_hat)^(-1/shape_hat)
  } else {
    1 - exp(-exc / scale_hat)
  }
  dfpp <- data.frame(Theoretical = sort(F_theo), Empirical = (1:n)/n)
  ppp <- ggplot(dfpp, aes(Theoretical, Empirical)) +
    geom_point(color = "darkgreen") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    labs(title = paste0("P–P (", var_label, ")"), x = "Theoretical CDF", y = "Empirical CDF") +
    theme_science_polished
  
  # PIT uniformity
  F_hat <- pgpd_(exc, sigma = scale_hat, xi = shape_hat)
  pks <- ggplot(data.frame(F_hat = F_hat), aes(F_hat)) +
    geom_histogram(aes(y = ..density..), bins = 20, fill = "skyblue", color = "black", alpha = 0.7) +
    geom_hline(yintercept = 1, color = "red", linetype = "dashed") +
    geom_density(color = "darkblue", linewidth = 1.1, adjust = 1.5) +
    labs(title = paste0("Uniformity: PIT (", var_label, ")"),
         x = expression(hat(F)(y)), y = "Density") +
    theme_science_polished
  
  list(pqq = pqq, ppp = ppp, pks = pks, F_hat = F_hat)
}

ks_boot_pvalue <- function(y, u, scale_hat, shape_hat, B = 1000, seed = 42) {
  set.seed(seed)
  exc <- y[y > u] - u
  n   <- length(exc)
  F_obs <- pgpd_(exc, sigma = scale_hat, xi = shape_hat)
  ks_obs <- suppressWarnings(ks.test(F_obs, "punif")$statistic)
  
  ks_b <- rep(NA_real_, B)
  for (b in seq_len(B)) {
    sim_exc <- rgpd_(n, sigma = scale_hat, xi = shape_hat)
    sim_y   <- c(y[y <= u], u + sim_exc) # keep subthresholds; simulate exceedances
    fb <- try(fevd(sim_y, type = "GP", threshold = u, time.units = "none"), silent = TRUE)
    if (inherits(fb, "try-error")) next
    par_b <- fb$results$par
    F_b <- pgpd_(sim_exc, sigma = par_b["scale"], xi = par_b["shape"])
    ks_b[b] <- suppressWarnings(ks.test(F_b, "punif")$statistic)
  }
  ks_b <- ks_b[is.finite(ks_b)]
  pval <- mean(ks_b >= ks_obs)
  list(ks_obs = as.numeric(ks_obs), p_value = pval, ks_boot = ks_b)
}

# Main pipeline for one variable
analyze_tail <- function(vec_raw, var_label,
                         log_transform = TRUE,
                         q_from = 0.75, q_to = 0.95, n_u = 50,
                         u_anchor_quantile = 0.85,
                         save_prefix = NULL,
                         out_dir = fig_dir) {
  
  vec_raw <- vec_raw[is.finite(vec_raw)]
  if (length(vec_raw) < 50) stop(paste0("Not enough data for ", var_label))
  y <- if (log_transform) log(vec_raw) else vec_raw
  
  # Threshold grid + anchor
  u_seq <- make_thresholds(y, q_from = q_from, q_to = q_to, n = n_u)
  u0    <- as.numeric(quantile(y, u_anchor_quantile, na.rm = TRUE))
  
  # MRL
  mrl_vals <- mrl_data(y, u_seq)
  p_mrl <- ggplot(data.frame(u = u_seq, mrl = mrl_vals), aes(u, mrl)) +
    geom_line() + geom_point() +
    geom_vline(xintercept = u0, linetype = "dashed", color = "red") +
    labs(title = paste0("MRL plot (", var_label, ")"),
         x = if (log_transform) paste0("Threshold (log ", var_label, ")") else paste0("Threshold (", var_label, ")"),
         y = "Mean excess") +
    theme_science_polished
  
  # Shape & adjusted scale stability
  shape <- scale <- shape_lo <- shape_hi <- adj <- adj_lo <- adj_hi <- rep(NA_real_, length(u_seq))
  for (i in seq_along(u_seq)) {
    u <- u_seq[i]
    exc <- y[y > u] - u
    if (length(exc) < 20) next
    out <- try(fit_gpd_at_u(y, u), silent = TRUE)
    if (inherits(out, "try-error")) next
    
    scale[i] <- out$scale
    shape[i] <- out$shape
    
    # Wald CI from covariance
    zcrit <- qnorm(1 - (1 - ci_level)/2)
    se_scale <- sqrt(out$cov[1,1])
    se_shape <- sqrt(out$cov[2,2])
    shape_lo[i] <- shape[i] - zcrit * se_shape
    shape_hi[i] <- shape[i] + zcrit * se_shape
    
    # adjusted scale @ u0 with delta method
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
    labs(title = paste0("Shape stability (", var_label, ")"),
         x = if (log_transform) paste0("Threshold (log ", var_label, ")") else paste0("Threshold (", var_label, ")"),
         y = "Shape (xi)") +
    theme_science_polished
  
  p_adj <- ggplot(data.frame(u = u_seq, adj, adj_lo, adj_hi), aes(u, adj)) +
    geom_point() +
    geom_errorbar(aes(ymin = adj_lo, ymax = adj_hi), width = 0.03, color = "blue") +
    geom_vline(xintercept = u0, color = "red", linetype = "dashed") +
    labs(title = paste0("Adjusted scale stability (", var_label, ")"),
         x = if (log_transform) paste0("Threshold (log ", var_label, ")") else paste0("Threshold (", var_label, ")"),
         y = expression(tilde(sigma)(u[0]))) +
    theme_science_polished
  
  # Final fit at u0 + diagnostics
  fit0 <- fit_gpd_at_u(y, u0)
  diags <- diagnostic_plots(y, u0, fit0$scale, fit0$shape, var_label)
  ksout <- ks_boot_pvalue(y, u0, fit0$scale, fit0$shape, B = n_boot_ks)
  
  message(sprintf("[%s] u0=%.3f | scale=%.4f shape=%.4f | KS p=%.3f",
                  var_label, u0, fit0$scale, fit0$shape, ksout$p_value))
  
  # Save PNGs
  prefix <- if (is.null(save_prefix)) gsub("\\s+", "_", tolower(var_label)) else save_prefix
  save_png <- function(p, name) {
    ggsave(file.path(out_dir, paste0(prefix, "_", name, ".png")), p,
           dpi = 600, width = 7, height = 5, units = "in")
  }
  save_png(p_mrl,  "mrl")
  save_png(p_shape,"shape_stability")
  save_png(p_adj,  "adj_scale_stability")
  save_png(diags$pqq, "qq")
  save_png(diags$ppp, "pp")
  save_png(diags$pks, "pit_uniformity")
  
  # Return compact summary
  list(
    var_label = var_label,
    u0 = u0,
    par = c(scale = fit0$scale, shape = fit0$shape),
    ks = ksout
  )
}

#### =========================
#### RUN THE THREE ANALYSES
#### =========================

# Flags to omit deformed tails in TL and Weight
omit_deformed_TL     <- TRUE
omit_deformed_weight <- TRUE

# 1) SVL (cm) — we analyze raw SVL, log-transform inside
svl_vec <- dat$SVL
res_SVL <- analyze_tail(
  vec_raw = svl_vec,
  var_label = "SVL (cm)",
  log_transform = TRUE,
  q_from = 0.75, q_to = 0.95, n_u = 50,
  u_anchor_quantile = 0.92,
  save_prefix = "svl",
  out_dir = fig_dir
)

# 2) TL (cm) — optionally omit deformed tails
tl_vec <- if (omit_deformed_TL && ("Deform" %in% names(dat))) {
  dat %>% filter(Deform == 0L) %>% pull(TL)
} else {
  dat$TL
}
res_TL <- analyze_tail(
  vec_raw = tl_vec,
  var_label = "TL (cm)",
  log_transform = TRUE,
  q_from = 0.75, q_to = 0.95, n_u = 50,
  u_anchor_quantile = 0.90,
  save_prefix = if (omit_deformed_TL) "tl_nodeform" else "tl_all",
  out_dir = fig_dir
)

# 3) Weight (kg) — optionally omit deformed tails (same Deform flag)
w_vec <- if (omit_deformed_weight && ("Deform" %in% names(dat))) {
  dat %>% filter(Deform == 0L) %>% pull(WKG)
} else {
  dat$WKG
}
res_W <- analyze_tail(
  vec_raw = w_vec,
  var_label = "Weight (kg)",
  log_transform = TRUE,      # many weights are skewed; keep log to stabilize tail
  q_from = 0.75, q_to = 0.99, n_u = 50,
  u_anchor_quantile = 0.95,
  save_prefix = if (omit_deformed_weight) "weight_nodeform" else "weight_all",
  out_dir = fig_dir
)

#### =========================
#### PRINT SUMMARY
#### =========================
cat("\n===== EVA summary (parameters at chosen u0) =====\n")
for (res in list(res_SVL, res_TL, res_W)) {
  cat(sprintf("%-12s | u0(log) = %.3f | scale = %.4f | shape = %.4f | KS p = %.3f\n",
              res$var_label, res$u0, res$par["scale"], res$par["shape"], res$ks$p_value))
}

# ML estimates for endpoint: 231 cm SVL, 425 cm TL, 

