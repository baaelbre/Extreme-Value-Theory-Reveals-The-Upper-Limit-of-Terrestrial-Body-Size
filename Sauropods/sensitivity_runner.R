# =========================
# Sensitivity runner
# =========================

# --- helpers ---
lambda0_from_k_alpha <- function(k, alpha) -log(alpha)/log(k)

kde_df <- function(x, from = quantile(x,0.001), to = quantile(x,0.999), n = 4096) {
  d <- density(x, from = from, to = to, n = n, kernel = "gaussian", adjust = 1.0)
  data.frame(t = d$x, d = d$y)
}
map_from_draws <- function(x) { d <- density(x, n = 4096); d$x[which.max(d$y)] }

overlap_trapz <- function(x, f, g) {
  stopifnot(length(x)==length(f), length(f)==length(g))
  pracma::trapz(x, pmin(f,g))
}

rtrunc_norm <- function(n, mu, sd, lo, hi) {
  out <- numeric(0)
  while (length(out) < n) {
    cand <- rnorm(n, mu, sd)
    cand <- cand[cand > lo & cand < hi]
    out  <- c(out, cand)
  }
  out[1:n]
}

# --- compile once ---
sm <- rstan::stan_model(stan_file)

# --- base stan data (from your script) ---
base_sd <- list(
  N  = length(y_above),
  y  = as.vector(y_above),
  u  = u0,
  ymax = ymax
)

# --- one fit with given hyper-hyperparameters ---
fit_once <- function(k_factor, alpha_pc, a_lambda,
                     mu_xi0, sd_xi0, hc_scale_xi,
                     iter = 2000, warmup = 1000, chains = 2, seed = 42,
                     compute_mass = FALSE, coeffs = NULL, include_eps = FALSE) {
  
  lambda0  <- lambda0_from_k_alpha(k_factor, alpha_pc)
  b_lambda <- a_lambda / lambda0   # E[lambda]=lambda0
  
  sd <- c(base_sd,
          list(mu_xi0 = mu_xi0, sd_xi0 = sd_xi0, hc_scale_xi = hc_scale_xi,
               a_lambda = a_lambda, b_lambda = b_lambda))
  
  fit <- rstan::sampling(sm, data = sd, seed = seed, chains = chains,
                         iter = iter, warmup = warmup,
                         control = list(adapt_delta = 0.95, max_treedepth = 12),
                         refresh = 0)
  
  post <- as.data.frame(fit)
  ystar_draw      <- post$ystar
  CFHstar_draw_mm <- exp(ystar_draw)
  xi_draw         <- post$xi
  delta_med       <- median(post$delta)
  lambda_med      <- median(post$lambda_gap)
  
  # CFH* summaries
  cfh_map  <- map_from_draws(CFHstar_draw_mm)
  cfh_up95 <- as.numeric(quantile(CFHstar_draw_mm, 0.95))
  
  # --- prior vs posterior overlap for CFH* (quantifies prior influence) ---
  R <- 80000
  lambda_prior <- rgamma(R, shape = a_lambda, rate = b_lambda)
  delta_prior  <- rexp(R, rate = lambda_prior)
  ystar_prior  <- ymax + delta_prior
  CFH_prior_mm <- exp(ystar_prior)
  # (xi prior not used in overlap for CFH*, but available if you want)
  
  dens_prior <- kde_df(CFH_prior_mm,
                       from = max(1000, exp(u0)*0.95),
                       to   = max(CFH_prior_mm, CFHstar_draw_mm)*1.25)
  dens_post  <- kde_df(CFHstar_draw_mm,
                       from = min(dens_prior$t), to = max(dens_prior$t))
  xg     <- dens_post$t
  f_prior<- approx(dens_prior$t, dens_prior$d, xout = xg, rule = 2)$y
  f_post <- dens_post$d
  ovl_cfh <- overlap_trapz(xg, f_prior, f_post)  # in [0,1], higher = more prior-like
  
  # Optional mass propagation
  mass_map <- NA_real_; mass_up95 <- NA_real_
  if (isTRUE(compute_mass) && !is.null(coeffs)) {
    stopifnot(all(c("alpha","beta","gamma","mean_log_sum_circ","V","resid_sd") %in% names(coeffs)))
    alpha <- as.numeric(coeffs$alpha)
    beta  <- as.numeric(coeffs$beta)
    gamma <- as.numeric(coeffs$gamma)
    c0    <- as.numeric(coeffs$mean_log_sum_circ)
    Vabg  <- as.matrix(coeffs$V)
    sigma_e <- as.numeric(coeffs$resid_sd)
    
    B    <- nrow(post); take <- min(20000, B)
    idx  <- sample.int(B, size = take, replace = TRUE)
    ystar_samp <- ystar_draw[idx]
    Z <- ystar_samp - c0
    ABG <- MASS::mvrnorm(take, mu = c(alpha,beta,gamma), Sigma = Vabg)
    if (include_eps) {
      eps <- rnorm(take, 0, sigma_e)
      m_log <- ABG[,1] + ABG[,2]*Z + ABG[,3]*Z^2 + eps
    } else {
      m_log <- ABG[,1] + ABG[,2]*Z + ABG[,3]*Z^2
    }
    mass_tons <- exp(m_log)/1e6
    mass_map  <- map_from_draws(mass_tons)
    mass_up95 <- as.numeric(quantile(mass_tons, 0.95))
  }
  
  data.frame(
    k_factor = k_factor, alpha_pc = alpha_pc, a_lambda = a_lambda,
    mu_xi0 = mu_xi0, sd_xi0 = sd_xi0, hc_scale_xi = hc_scale_xi,
    lambda0 = lambda0,
    cfh_map_mm = cfh_map, cfh_up95_mm = cfh_up95,
    xi_mean = mean(xi_draw), xi_q05 = quantile(xi_draw, 0.05),
    xi_q95 = quantile(xi_draw, 0.95),
    delta_med = delta_med, lambda_med = lambda_med,
    overlap_CFHstar = ovl_cfh,
    mass_map_t = mass_map, mass_up95_t = mass_up95
  )
}

# --- grid to scan (tune freely) ---
grid <- expand.grid(
  k_factor = c(1.2, 1.5, 2, 3),
  alpha_pc = c(0.05, 0.10, 0.20),
  a_lambda = c(1.0, 2.0, 5.0),
  mu_xi0   = c(-0.30, -0.25),
  sd_xi0   = c(0.05, 0.10, 0.20),
  hc_scale_xi = c(0.05, 0.10)
)

# Optional: load allometry coefficients if you want mass outputs
# coeffs <- readRDS("centered_quadratic_coefficients.rds")

# --- run (small defaults to keep it quick; increase iter if needed) ---
sens_res <- dplyr::bind_rows(lapply(seq_len(nrow(grid)), function(i) {
  g <- grid[i,]
  fit_once(
    k_factor = g$k_factor, alpha_pc = g$alpha_pc, a_lambda = g$a_lambda,
    mu_xi0 = g$mu_xi0, sd_xi0 = g$sd_xi0, hc_scale_xi = g$hc_scale_xi,
    iter = 1500, warmup = 750, chains = 2, seed = 42,
    compute_mass = FALSE, coeffs = NULL, include_eps = FALSE
  )
}))

# --- quick looks ---
print(dplyr::arrange(sens_res, desc(cfh_up95_mm))[1:10, ])

# CFH* 95% upper sensitivity vs k (faceted by alpha, linetype = a_lambda)
ggplot(sens_res,
       aes(k_factor, cfh_up95_mm, color = factor(sd_xi0),
           linetype = factor(a_lambda), group = interaction(sd_xi0, a_lambda))) +
  geom_line() + geom_point() +
  facet_wrap(~ alpha_pc, labeller = label_bquote(alpha==.(alpha_pc))) +
  labs(x = "k (factor beyond max with prob α)", y = "CFH* 95% upper (mm)",
       color = "sd_xi0", linetype = "a_lambda") +
  theme_science_polished
ggsave(file.path(FIG_DIR, "sens_CFHstar_95upper_vs_k.png"), dpi=600, w=8, h=5)

# Prior influence metric: lower overlap => data dominate; higher => prior dominates
ggplot(sens_res, aes(lambda0, overlap_CFHstar, color=factor(a_lambda))) +
  geom_point() + geom_line() +
  labs(x = expression(lambda[0]==-log(alpha)/log(k)),
       y = "Overlap(prior, posterior) on CFH*",
       color = "a_lambda") +
  theme_science_polished
ggsave(file.path(FIG_DIR, "sens_overlap_vs_lambda0.png"), dpi=600, w=7, h=5)
