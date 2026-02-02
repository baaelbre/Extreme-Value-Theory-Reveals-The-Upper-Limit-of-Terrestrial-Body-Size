############################################################
#  Title: POT analysis on log CFH       
#  Author: Bastiaan Aelbrecht                               #
#  Date  : 2025‑08‑08                                        #
############################################################

############################
# 1.        SETUP          #
############################
library(ggplot2)
library(pracma)
library(HDInterval)
library(readxl)
library(dplyr)
library(knitr)
library(scales)
theme_science <- theme_minimal(base_family = "Arial", base_size = 12) +
  theme(
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 16),
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 10),
    legend.key.size = unit(0.5, "cm"),
    panel.grid.major = element_line(color = "gray80", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    axis.ticks.length = unit(0.20, "cm"),
    axis.ticks = element_line(color = "black", linewidth = 0.4),
    plot.margin = margin(5, 5, 5, 5),
    legend.position = "right"
  )

# ---- 1.3  Helper: trapezoidal integration wrapper ----
trapz_safe <- function(x, y) pracma::trapz(x, y)

# Directory for figure output
FIG_DIR <- "Figures_sauropods"
if (!dir.exists(FIG_DIR)) dir.create(FIG_DIR)

############################
# 2.  DATA LOADING         #
############################

load_cfh_data <- function(path) {
  read_excel(path) %>%
    select(`hum+fem circ (mm)`, `publication year`) %>%
    setNames(c("sum_circ_mm", "year")) %>%
    na.omit() %>%
    mutate(log_circ = log(as.double(sum_circ_mm)))
}

df_cfh <- load_cfh_data("Data/DEmic23_updated_Supplemental_Data_withPubYear.xlsx")
log_circ  <- df_cfh$log_circ
x_max     <- max(log_circ)

############################
# 3.  EXPLORATORY PLOTS    #
############################

plot_hist_density <- function(df) {
  ggplot(df, aes(x = log_circ)) +
    geom_histogram(aes(y = ..density..), binwidth = 0.05,
                   fill = "skyblue", color = "black", alpha = 0.7) +
    geom_density(color = "darkblue", size = 1.2, adjust = 1.2) +
    scale_x_continuous(
      name = expression(ln(C[F+H])~"[ln(mm)]"),
      sec.axis = sec_axis(~exp(.), name = expression(C[F+H]~"[mm]"))
    ) +
    labs(y = "Density") +
    theme_science
}

plot_qq <- function(df) {
  ggplot(df, aes(sample = log_circ)) +
    stat_qq(color = "darkred", size = 2) +
    stat_qq_line(color = "black", linetype = "dashed", size = 1) +
    labs(x = "Theoretical Quantiles", y = "Sample Quantiles") +
    theme_science
}

p_hist <- plot_hist_density(df_cfh)
p_qq   <- plot_qq(df_cfh)

p_hist
p_qq

############################
# 4.  THRESHOLD STABILITY  #
############################

# -- 4.1  Estimate parameters & 95% Wald CIs for a grid of thresholds --
find_threshold <- function(x, u0 = 6.31, n_grid = 50, alpha = 0.05) {
  z_crit <- qnorm(1 - alpha/2)
  thresholds <- seq(u0, 7.44, length.out = n_grid)
  
  results <- purrr::map(thresholds, function(u) {
    fit <- fevd(x, threshold = u, type = "GP")
    pars <- fit$results$par
    
    # Wald SEs via observed information (extRemes stores them in $se)
    se  <- summary(fit)$se
    cov <- tryCatch(summary(fit)$cov.theta, error = function(e) diag(2))
    
    # Parameter estimates
    sigma_hat <- pars["scale"]
    xi_hat    <- pars["shape"]
    
    # Wald CI for sigma & xi
    sigma_ci <- sigma_hat + c(-1, 1) * z_crit * se["scale"]
    xi_ci    <- xi_hat    + c(-1, 1) * z_crit * se["shape"]
    
    # Adjusted scale σ*(u) = σ(u) − ξ(u)(u − u0)
    sigma_adj_hat <- sigma_hat - xi_hat * (u - u0)
    
    # Delta‑method variance for σ* : Var(σ) + (u−u0)^2 Var(ξ) − 2(u−u0)Cov(σ,ξ)
    var_adj <- cov[1,1] + (u - u0)^2 * cov[2,2] - 2*(u - u0) * cov[1,2]
    sd_adj  <- sqrt(var_adj)
    adj_ci  <- sigma_adj_hat + c(-1, 1) * z_crit * sd_adj
    
    tibble::tibble(
      threshold = u,
      sigma     = sigma_hat,
      sigma_l   = sigma_ci[1],
      sigma_u   = sigma_ci[2],
      xi        = xi_hat,
      xi_l      = xi_ci[1],
      xi_u      = xi_ci[2],
      sigma_adj = sigma_adj_hat,
      sigma_adj_l = adj_ci[1],
      sigma_adj_u = adj_ci[2]
    )
  }) %>%
    dplyr::bind_rows()
  
  results
}

th_tbl <- find_threshold(log_circ)
threshold_opt <- th_tbl$threshold[22]  # chosen manually (index 22 ≈ 6.77)

# -- 4.2  Plot helpers with 95% error bars --
plot_param_vs_threshold <- function(tbl, y, y_low, y_up, ylab, ref_idx) {
  ggplot(tbl, aes(x = threshold, y = .data[[y]])) +
    geom_point() +
    geom_errorbar(aes(ymin = .data[[y_low]], ymax = .data[[y_up]]), width = 0.08, colour = "blue") +
    geom_vline(xintercept = tbl$threshold[ref_idx], colour = "red", linetype = "dashed") +
    labs(x = expression(bold("Threshold (log " * C[F+H] * ")")), y = ylab) +
    theme_science
}

p_param_adj_scale <- plot_param_vs_threshold(th_tbl, "sigma_adj", "sigma_adj_l", "sigma_adj_u",
                                             "Adjusted scale", 22)

p_param_shape     <- plot_param_vs_threshold(th_tbl, "xi", "xi_l", "xi_u",
                                             "Shape (xi)", 22)
p_param_adj_scale
p_param_shape
############################
# 5.  GPD FIT & GOF       #
############################

fit_gpd_model <- function(x, u) fevd(x, threshold = u, type = "GP")
fit_gpd       <- fit_gpd_model(log_circ, threshold_opt)
xi_hat        <- fit_gpd$results$par["shape"]
sigma_hat     <- fit_gpd$results$par["scale"]

make_gpd_qq_pp <- function(x, u, xi, sigma) {
  excess  <- x[x > u] - u
  n       <- length(excess)
  probs   <- ppoints(n)
  theo_q  <- if (abs(xi) > 1e-6) u + sigma/xi * (probs^(-xi) - 1) else u - sigma * log(probs)
  qq_df   <- data.frame(Theoretical = rev(theo_q), Empirical = sort(x[x > u]))
  pp_df   <- data.frame(Theoretical = sort(if (abs(xi) > 1e-6)
    1 - (1 + xi * excess / sigma)^(-1/xi)
    else 1 - exp(-excess / sigma)),
    Empirical = (1:n)/n)
  list(qq_df = qq_df, pp_df = pp_df)
}

gof_df <- make_gpd_qq_pp(log_circ, threshold_opt, xi_hat, sigma_hat)

p_qq_gpd <- ggplot(gof_df$qq_df, aes(Theoretical, Empirical)) +
  geom_point(color = "steelblue") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Q–Q plot: GPD", x = "Theoretical", y = "Empirical") +
  theme_science

p_pp_gpd <- ggplot(gof_df$pp_df, aes(Theoretical, Empirical)) +
  geom_point(color = "darkgreen") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(title = "P–P plot: GPD", x = "Theoretical CDF", y = "Empirical CDF") +
  theme_science

p_qq_gpd
p_pp_gpd

############################
# 6.  ENDPOINT POSTERIOR   #
############################

# Joint RW–MH in (phi, rho) where
#   phi = log(y* - y_min)   ∈ ℝ
#   rho = log(-xi)          ∈ ℝ
posterior_endpoint <- function(x, u,
                               n_iter = 60000,
                               burn   = 10000,
                               thin   = 20,
                               prop_sd = c(0.15, 0.10)) {
  
  excess <- x[x > u] - u
  n      <- length(excess)
  y_min  <- max(excess) + 1e-6        # smallest legal (y* - u)
  
  # ----- log-posterior in transformed coords ---------------------
  log_post <- function(phi, rho) {
    y_star <- u + y_min + exp(phi)    # back-transform
    xi     <- -exp(rho)               # negative
    
    # log-likelihood (GPD, eq. (1) in notes)
    ll <- n * log(y_star - u) / xi -
      (1 / xi + 1) * sum(log(y_star - u - excess)) -
      n * log(-xi)
    
    # Flat priors  ⇒  Jacobian term φ+ρ
    ll + phi + rho
  }
  
  # ----- Metropolis sampler --------------------------------------
  phi_cur <- log(0.5)                 # y* starts 0.5 above y_min
  rho_cur <- log(0.1)                 # xi starts at −0.1
  lp_cur  <- log_post(phi_cur, rho_cur)
  
  chain <- matrix(NA_real_, n_iter, 2)
  acc   <- 0
  
  for (t in seq_len(n_iter)) {
    prop    <- rnorm(2, c(phi_cur, rho_cur), prop_sd)
    lp_prop <- log_post(prop[1], prop[2])
    
    if (log(runif(1)) < lp_prop - lp_cur) {
      phi_cur <- prop[1]; rho_cur <- prop[2]; lp_cur <- lp_prop
      acc <- acc + 1
    }
    chain[t, ] <- c(phi_cur, rho_cur)
  }
  message("Acceptance rate: ", round(acc / n_iter, 3))
  
  # ----- Burn-in, thinning, back-transform -----------------------
  keep   <- seq(burn + 1, n_iter, by = thin)
  phi    <- chain[keep, 1]
  rho    <- chain[keep, 2]
  
  y_samp <- u + y_min + exp(phi)      # draws of y*
  xi_samp <- -exp(rho)                # (optional) draws of xi
  
  # Density for plotting
  kd <- density(y_samp, n = 1024)
  list(
    y       = kd$x,
    density = kd$y / trapz_safe(kd$x, kd$y),
    samples = y_samp,
    xi      = xi_samp                 # returned for diagnostics
  )
}

# ---- 6.3  Run sampler & summaries -------------------------------
set.seed(42)
post_ep <- posterior_endpoint(log_circ, threshold_opt, prop_sd=c(0.4,0.5))

hpdi_ep <- HDInterval::hdi(post_ep$samples, credMass = 0.95)
map_ep  <- post_ep$y[which.max(post_ep$density)]

############################
# 7.  MASS EXTRAPOLATION   #
############################

library(mvtnorm)

mass_extrapolation_mc <- function(y_star_samples, coeff_file,
                                  n_draws = 5e4,   # total MC size
                                  dens_n  = 1024,  # grid for KDE
                                  max_mass_t = 1000) {
  
  ## --- 7.1  Load regression coefficients ------------------------
  coeffs <- readRDS(coeff_file)
  mu_alpha <- coeffs$alpha
  mu_beta  <- coeffs$beta
  mu_gamma <- coeffs$gamma
  
  Sigma <- matrix(c(coeffs$alpha_se^2,  coeffs$cov_ab,       coeffs$cov_ag,
                    coeffs$cov_ab,      coeffs$beta_se^2,    coeffs$cov_bg,
                    coeffs$cov_ag,      coeffs$cov_bg,       coeffs$gamma_se^2),
                  3, 3, byrow = TRUE)
  
  resid_sd <- coeffs$resid_sd
  xbar     <- coeffs$mean_log_sum_circ
  
  ## --- 7.2  Monte-Carlo sampling --------------------------------
  y_star_mc <- sample(y_star_samples, n_draws, replace = TRUE)
  
  theta_mc  <- mvtnorm::rmvnorm(n_draws,
                                mean = c(mu_alpha, mu_beta, mu_gamma),
                                sigma = Sigma)
  alpha_mc <- theta_mc[, 1]
  beta_mc  <- theta_mc[, 2]
  gamma_mc <- theta_mc[, 3]
  
  centered_x <- y_star_mc - xbar
  mu_logM    <- alpha_mc +
    beta_mc  * centered_x +
    gamma_mc * centered_x^2
  
  log_mass_mc <- rnorm(n_draws, mean = mu_logM, sd = resid_sd)
  
  ## --- 7.3  Convert to tonnes & density estimate ----------------
  mass_t_mc <- exp(log_mass_mc) / 1e6          # kg → t
  kd <- density(mass_t_mc,
                n   = dens_n,
                from = 0,
                to   = max_mass_t)
  
  ## --- 7.4  summaries -------------------------------------------
  map_mass <- kd$x[which.max(kd$y)]
  ci95     <- HDInterval::hdi(mass_t_mc, credMass = 0.95)   # <-- HPDI
  
  list(
    df        = data.frame(mass_tons = kd$x,
                           density   = kd$y / trapz_safe(kd$x, kd$y)),
    map       = map_mass,
    ci95      = ci95,
    samples   = mass_t_mc               # optional: further analyses
  )
}


# ---------- Run Monte-Carlo propagation --------------------------
set.seed(123)
mass_out <- mass_extrapolation_mc(post_ep$samples,
                                  "centered_quadratic_coefficients.rds")

mass_df  <- mass_out$df
map_mass <- mass_out$map
ci95     <- mass_out$ci95

## ---------- Plot -------------------------------------------------
p_mass_endpoint <- ggplot(mass_df, aes(mass_tons, density)) +
  geom_line(color = "darkblue", linewidth = 1.2) +
  geom_area(data = subset(mass_df,
                          mass_tons >= ci95[1] & mass_tons <= ci95[2]),
            fill = "skyblue", alpha = 0.4) +
  geom_vline(xintercept = map_mass, colour = "purple",
             linetype = "dashed", linewidth = 1) +
  geom_vline(xintercept = ci95, colour = "black",
             linetype = "dashed") +
  labs(x = "Mass endpoint (tons)",
       y = "Posterior density") +
  theme_science
p_mass_endpoint


############################
# 8.  EXPORT PLOTS         #
############################

plot_list <- list(p_hist, p_qq, p_param_adj_scale, p_param_shape,
                  p_qq_gpd, p_pp_gpd, p_endpoint_mm, p_mass_endpoint)
file_names <- c("hist_density", "qq", "adj_scale", "shape_param",
                "gpd_qq", "gpd_pp", "endpoint_mm", "mass_endpoint")

for (i in seq_along(plot_list)) {
  ggsave(
    filename = file.path(FIG_DIR, paste0(file_names[i], ".png")),
    plot = plot_list[[i]],
    device = "png", dpi = 600,
    width = 7, height = 5, units = "in"
  )
  ggsave(
    filename = file.path(FIG_DIR, paste0(file_names[i], ".pdf")),
    plot = plot_list[[i]],
    device = cairo_pdf, 
    width = 7, height = 5, units = "in"
  )}

