# ============================================================
# Bivariate EVT with censored likelihood & logistic dependence
# (SVL, TL) — American alligators
# Full log-likelihood: censored bivariate POT (bvpot.c:nllbvclog)
# + univariate GP contributions for missing partners (evd::dgpd)
# Parametric bootstrap for (sigma1, sigma2, xi, y1*, y2*)
# + univariate KS bootstrap + Wald tests for xi<0
# ============================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(ggplot2)
  library(evd)        # fpot, pgpd, dgpd, qgpd, sep.bvdata
  library(scales)
  library(glue)
  library(forcats)
  library(grid)
  library(MASS)       # for ginv fallback if Hessian singular
  library(copula)     # Gumbel / logistic EV copula
})

set.seed(42)

# ------------------------------------------------------------
# 0. Load compiled C code (bvpot.c → bvpot.so / bvpot.dll)
# ------------------------------------------------------------
# system("R CMD SHLIB bvpot.c")
if (!is.loaded("nllbvclog")) {
  ## Adjust path / extension as needed:
  ## dyn.load("path/to/bvpot.so") or "bvpot.dll" on Windows
  dyn.load("bvpot.dll")
}

# ---------------------------
# Directories & theming
# ---------------------------
FIG_DIR <- "Figures_gators"
if (!dir.exists(FIG_DIR)) dir.create(FIG_DIR, recursive = TRUE)

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

# ---------------------------
# Global settings
# ---------------------------
ci_level      <- 0.90
u_scan_lo_q   <- 0.60
u_scan_hi_q   <- 0.99
u_scan_n      <- 50
min_ex_mrl    <- 5
min_ex_fit    <- 20

trait_names   <- c("SVL","TL")
thresh_q_opt  <- c(SVL = 0.94, TL = 0.94)  # anchor quantiles (log scale)

# Bootstrap settings (bivariate endpoint & KS)
B_boot <- 1000L   # parametric bivariate bootstrap replicates
B_KS   <- 1000L   # parametric KS bootstrap per margin

# ============================================================
# 1. Ingest and basic cleaning
# ============================================================

DATA_XLSX <- "Data/experimental_alligator_harvest_woodward.xlsx"
df_raw    <- read_excel(DATA_XLSX)

stopifnot(all(c("SVL","TL","Deform") %in% names(df_raw)))

# Deform == 1 or 3: tail broken ⇒ TL structurally missing; SVL kept.
df <- df_raw %>%
  mutate(
    SVL = suppressWarnings(as.numeric(SVL)),
    TL  = suppressWarnings(as.numeric(TL)),
    TL  = ifelse(Deform %in% c(1, 3), NA_real_, TL)
  ) %>%
  transmute(
    specimen = row_number(),
    SVL, TL
  )

# ============================================================
# 2. Completeness diagnostics
# ============================================================

compl_tbl <- trait_names %>%
  set_names() %>%
  map_df(function(tr) {
    v <- df[[tr]]
    tibble(
      trait        = tr,
      n_total      = length(v),
      n_obs        = sum(!is.na(v)),
      completeness = mean(!is.na(v))
    )
  })

p_compl <- compl_tbl %>%
  mutate(trait = fct_inorder(trait)) %>%
  ggplot(aes(trait, completeness)) +
  geom_col(fill = "#3B82F6") +
  geom_text(aes(label = percent(completeness, accuracy = 0.1)),
            vjust = -0.2, size = 4) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1.10)) +
  labs(title = "Completeness (SVL, TL)", x = NULL, y = "Completeness") +
  theme_science_polished

ggsave(file.path(FIG_DIR, "completeness_SVL_TL.png"), p_compl,
       dpi = 600, w = 6.0, h = 4.2, units = "in")

# ============================================================
# 3. Log-transform and split patterns (complete / SVL-only / TL-only)
# ============================================================

df_log <- df %>%
  mutate(
    log_SVL = ifelse(is.finite(SVL) & SVL > 0, log(SVL), NA_real_),
    log_TL  = ifelse(is.finite(TL)  & TL  > 0, log(TL),  NA_real_)
  )

# Complete cases: both SVL and TL observed
df_complete <- df_log %>%
  filter(!is.na(log_SVL), !is.na(log_TL))

# SVL only
df_svl_only <- df_log %>%
  filter(!is.na(log_SVL), is.na(log_TL))

# TL only
df_tl_only <- df_log %>%
  filter(is.na(log_SVL), !is.na(log_TL))

cat("\nRow counts by observation pattern:\n")
cat("  complete (SVL & TL) :", nrow(df_complete), "\n")
cat("  SVL only            :", nrow(df_svl_only), "\n")
cat("  TL only             :", nrow(df_tl_only), "\n")

# ============================================================
# 4. Thresholds on log scale (per trait, using *all* non-missing)
# ============================================================

u0_by_trait <- setNames(numeric(length(trait_names)), trait_names)
for (tr in trait_names) {
  v_log <- df_log[[paste0("log_", tr)]]
  u0_by_trait[tr] <- as.numeric(quantile(v_log, thresh_q_opt[tr], na.rm = TRUE))
}
print(u0_by_trait)
print(exp(u0_by_trait))  # thresholds on original scale

# ============================================================
# 5. Helper functions for univariate diagnostics (evd only)
# ============================================================

make_thresholds <- function(y, q_from = u_scan_lo_q, q_to = u_scan_hi_q, n = u_scan_n) {
  rng <- quantile(y, c(q_from, q_to), na.rm = TRUE)
  seq(rng[1], rng[2], length.out = n)
}

fit_gpd_at_u <- function(y, u) {
  fit <- evd::fpot(y, threshold = u, model = "gpd")
  par <- fit$estimate
  cov <- fit$var.cov
  list(
    fit      = fit,
    scale    = unname(par["scale"]),
    shape    = unname(par["shape"]),
    cov      = cov,
    n_exceed = sum(y > u)
  )
}

adj_scale_fun <- function(scale, xi, u, u0) {
  scale - xi * (u - u0)
}

mrl_data <- function(y, u_seq, min_ex = 5) {
  sapply(u_seq, function(u) {
    ex <- y[y > u] - u
    if (length(ex) < min_ex) return(NA_real_)
    mean(ex)
  })
}

diagnostic_plots <- function(y, u, sigma_hat, xi_hat, label) {
  exc   <- y[y > u] - u
  n     <- length(exc)
  probs <- ppoints(n)
  
  if (abs(xi_hat) > 1e-10) {
    theo_q <- u + sigma_hat / xi_hat * (probs^(-xi_hat) - 1)
  } else {
    theo_q <- u - sigma_hat * log(probs)
  }
  
  dfqq <- data.frame(
    Theoretical = rev(theo_q),
    Empirical   = sort(y[y > u])
  )
  pqq <- ggplot(dfqq, aes(Theoretical, Empirical)) +
    geom_point(color = "steelblue") +
    geom_abline(slope = 1, intercept = 0,
                linetype = "dashed", color = "red") +
    labs(title = glue("Q–Q: {label}")) +
    theme_science_polished
  
  F_theo <- if (abs(xi_hat) > 1e-10) {
    1 - (1 + xi_hat * exc / sigma_hat)^(-1 / xi_hat)
  } else {
    1 - exp(-exc / sigma_hat)
  }
  dfpp <- data.frame(
    Theoretical = sort(F_theo),
    Empirical   = (1:n) / n
  )
  ppp <- ggplot(dfpp, aes(Theoretical, Empirical)) +
    geom_point(color = "darkgreen") +
    geom_abline(slope = 1, intercept = 0,
                linetype = "dashed", color = "red") +
    labs(title = glue("P–P: {label}")) +
    theme_science_polished
  
  # PIT using evd::pgpd (exc >= 0, loc = 0)
  F_hat <- evd::pgpd(exc, loc = 0, scale = sigma_hat, shape = xi_hat)
  pks <- ggplot(data.frame(F_hat = F_hat), aes(F_hat)) +
    geom_histogram(aes(y = ..density..), bins = 20,
                   fill = "skyblue", color = "black", alpha = 0.7) +
    geom_hline(yintercept = 1, color = "red", linetype = "dashed") +
    labs(title = glue("Uniformity (PIT): {label}")) +
    theme_science_polished
  
  list(pqq = pqq, ppp = ppp, pks = pks)
}

# ------------------------------------------------------------
# KS goodness-of-fit via parametric bootstrap (univariate GPD)
# ------------------------------------------------------------
ks_gpd_boot_from_fit <- function(y, u, scale_hat, shape_hat,
                                 B = 1000L, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  # Exceedances over u
  y_exceed <- y[y > u]
  ex       <- y_exceed - u
  n_ex     <- length(ex)
  if (n_ex < 20L) {
    warning("Not enough exceedances for KS bootstrap (n_ex < 20); returning NA.")
    return(list(D = NA_real_, p = NA_real_, D_boot = numeric(0L), n_exceed = n_ex))
  }
  
  # Observed KS statistic under fitted GPD
  Fn_ex      <- ecdf(ex)
  ex_sorted  <- sort(ex)
  Fgpd_obs   <- evd::pgpd(ex_sorted, loc = 0, scale = scale_hat, shape = shape_hat)
  Fn_vals    <- Fn_ex(ex_sorted)
  D_obs      <- max(abs(Fn_vals - Fgpd_obs))
  
  # Parametric bootstrap under fitted GPD with refitting
  D_boot <- numeric(B)
  
  for (b in seq_len(B)) {
    ex_b <- evd::rgpd(n_ex, loc = 0, scale = scale_hat, shape = shape_hat)
    y_b  <- u + ex_b
    
    fit_b <- try(evd::fpot(y_b, threshold = u, model = "gpd"), silent = TRUE)
    if (inherits(fit_b, "try-error")) {
      D_boot[b] <- NA_real_
      next
    }
    
    sc_b <- unname(fit_b$estimate["scale"])
    xi_b <- unname(fit_b$estimate["shape"])
    
    Fn_b       <- ecdf(ex_b)
    ex_b_sort  <- sort(ex_b)
    Fgpd_b     <- evd::pgpd(ex_b_sort, loc = 0, scale = sc_b, shape = xi_b)
    Fn_b_vals  <- Fn_b(ex_b_sort)
    D_boot[b]  <- max(abs(Fn_b_vals - Fgpd_b))
  }
  
  D_boot <- D_boot[is.finite(D_boot)]
  if (!length(D_boot)) {
    warning("All KS bootstrap replicates failed; returning NA.")
    return(list(D = D_obs, p = NA_real_, D_boot = numeric(0L), n_exceed = n_ex))
  }
  
  p_boot <- mean(D_boot >= D_obs)
  
  list(D = D_obs, p = p_boot, D_boot = D_boot, n_exceed = n_ex)
}

# ------------------------------------------------------------
# Trait-wise diagnostics, KS, and Wald p-values
# ------------------------------------------------------------
run_trait_diag <- function(tr_key) {
  v_raw <- df[[tr_key]] |> suppressWarnings(as.numeric())
  y     <- log(v_raw[is.finite(v_raw) & v_raw > 0])
  ylab  <- glue("log {tr_key}")
  stopifnot(length(y) >= 40)
  
  q_anchor <- thresh_q_opt[tr_key]
  u_seq    <- make_thresholds(y)
  u0       <- quantile(y, q_anchor, na.rm = TRUE) |> as.numeric()
  
  # Mean residual life
  mrl_vals <- mrl_data(y, u_seq, min_ex = min_ex_mrl)
  p_mrl <- ggplot(data.frame(u = u_seq, mrl = mrl_vals), aes(u, mrl)) +
    geom_line() + geom_point() +
    geom_vline(xintercept = u0, linetype = "dashed", color = "red") +
    labs(
      title = glue("MRL plot ({ylab})"),
      x     = glue("Threshold ({ylab})"),
      y     = "Mean excess"
    ) +
    theme_science_polished
  ggsave(file.path(FIG_DIR, glue("{tr_key}_log_mrl.png")),
         p_mrl, dpi = 600, w = 6.2, h = 4.2, units = "in")
  
  # GPD fit at anchor threshold u0
  fit0  <- fit_gpd_at_u(y, u0)
  xi0   <- fit0$shape
  sc0   <- fit0$scale
  cov0  <- fit0$cov
  xi_se <- sqrt(cov0[2, 2])
  sigma_se <- sqrt(cov0[1, 1])
  
  # Wald test for xi < 0 (one-sided)
  z_wald <- xi0 / xi_se
  p_wald <- pnorm(z_wald)  # H1: xi < 0
  
  # KS goodness-of-fit via parametric bootstrap
  seed_ks <- switch(tr_key,
                    "SVL" = 101L,
                    "TL"  = 202L,
                    303L)
  ks_res <- ks_gpd_boot_from_fit(
    y          = y,
    u          = u0,
    scale_hat  = sc0,
    shape_hat  = xi0,
    B          = B_KS,
    seed       = seed_ks
  )
  
  zcrit <- qnorm(1 - (1 - ci_level) / 2)
  
  df_scan <- tibble(
    u        = u_seq,
    shape    = NA_real_,
    shape_lo = NA_real_,
    shape_hi = NA_real_,
    adj      = NA_real_,
    adj_lo   = NA_real_,
    adj_hi   = NA_real_
  )
  
  for (i in seq_along(u_seq)) {
    u <- u_seq[i]
    ex <- y[y > u] - u
    if (length(ex) < min_ex_fit) next
    
    out <- try(fit_gpd_at_u(y, u), silent = TRUE)
    if (inherits(out, "try-error")) next
    
    se_shape <- sqrt(out$cov[2, 2])
    df_scan$shape[i]    <- out$shape
    df_scan$shape_lo[i] <- out$shape - zcrit * se_shape
    df_scan$shape_hi[i] <- out$shape + zcrit * se_shape
    
    df_scan$adj[i]      <- adj_scale_fun(out$scale, out$shape, u, u0)
    var_adj <- out$cov[1, 1] +
      (u - u0)^2 * out$cov[2, 2] -
      2 * (u - u0) * out$cov[1, 2]
    se_adj  <- sqrt(max(var_adj, 0))
    df_scan$adj_lo[i] <- df_scan$adj[i] - zcrit * se_adj
    df_scan$adj_hi[i] <- df_scan$adj[i] + zcrit * se_adj
  }
  
  p_shape <- ggplot(df_scan, aes(u, shape)) +
    geom_point(color = "blue") +
    geom_errorbar(aes(ymin = shape_lo, ymax = shape_hi),
                  width = 0.03, color = "blue") +
    geom_vline(xintercept = u0, color = "red", linetype = "dashed") +
    geom_hline(yintercept = xi0, color = "red", linetype = "dashed") +
    labs(
      x = glue("Threshold ({ylab})"),
      y = "Shape (xi)") +
    theme_science_polished
  
  p_adj <- ggplot(df_scan, aes(u, adj)) +
    geom_point(color = "blue") +
    geom_errorbar(aes(ymin = adj_lo, ymax = adj_hi),
                  width = 0.03, color = "blue") +
    geom_vline(xintercept = u0, color = "red", linetype = "dashed") +
    geom_hline(yintercept = sc0, color = "red", linetype = "dashed") +
    labs(
      x = glue("Threshold ({ylab})"),
      y = "Adjusted scale") +
    theme_science_polished
  
  ggsave(file.path(FIG_DIR, glue("{tr_key}_log_shape_stability.png")),
         p_shape, dpi = 600, w = 6.2, h = 4.2, units = "in")
  ggsave(file.path(FIG_DIR, glue("{tr_key}_log_adj_scale_stability.png")),
         p_adj,   dpi = 600, w = 6.2, h = 4.2, units = "in")
  
  di <- diagnostic_plots(y, u0, sc0, xi0, ylab)
  ggsave(file.path(FIG_DIR, glue("{tr_key}_log_qq.png")),
         di$pqq, dpi = 600, w = 6.2, h = 4.2, units = "in")
  ggsave(file.path(FIG_DIR, glue("{tr_key}_log_pp.png")),
         di$ppp, dpi = 600, w = 6.2, h = 4.2, units = "in")
  ggsave(file.path(FIG_DIR, glue("{tr_key}_log_pit_uniformity.png")),
         di$pks, dpi = 600, w = 6.2, h = 4.2, units = "in")
  
  cat(glue("\n{tr_key}: Wald test xi<0 at u0={round(u0,3)}: xî={round(xi0,3)}, se={round(xi_se,3)}, z={round(z_wald,3)}, p={signif(p_wald,3)}"))
  cat(glue("\n{tr_key}: KS bootstrap: D={round(ks_res$D,3)}, p_boot={signif(ks_res$p,3)}, n_exceed={ks_res$n_exceed}\n"))
  
  tibble(
    trait           = tr_key,
    u0              = u0,
    xi_hat_diag     = xi0,
    sigma_hat_diag  = sc0,
    xi_se           = xi_se,
    sigma_se        = sigma_se,
    wald_z_xi       = z_wald,
    wald_p_xi_neg   = p_wald,
    ks_D            = ks_res$D,
    ks_p_boot       = ks_res$p,
    ks_n_exceed     = ks_res$n_exceed
  )
}

diag_SVL <- run_trait_diag("SVL")
diag_TL  <- run_trait_diag("TL")
marg_diag_tbl <- bind_rows(diag_SVL, diag_TL)
print(marg_diag_tbl)

# ============================================================
# 6. Univariate POT fits using all non-missing data per margin
# ============================================================

y_SVL_all <- df_log$log_SVL[!is.na(df_log$log_SVL)]
y_TL_all  <- df_log$log_TL[!is.na(df_log$log_TL)]

u1 <- u0_by_trait["SVL"]
u2 <- u0_by_trait["TL"]

fit_SVL_uni <- evd::fpot(y_SVL_all, threshold = u1, model = "gpd")
fit_TL_uni  <- evd::fpot(y_TL_all,  threshold = u2, model = "gpd")

par_SVL_uni <- fit_SVL_uni$estimate
par_TL_uni  <- fit_TL_uni$estimate

sigma1_hat_uni <- unname(par_SVL_uni["scale"])
xi1_hat_uni    <- unname(par_SVL_uni["shape"])

sigma2_hat_uni <- unname(par_TL_uni["scale"])
xi2_hat_uni    <- unname(par_TL_uni["shape"])

cat("\nUnivariate GP fits (log scale, evd::fpot, ALL marginal data):\n")
cat("  SVL: scale =", sigma1_hat_uni, " shape =", xi1_hat_uni, "\n")
cat("  TL : scale =", sigma2_hat_uni, " shape =", xi2_hat_uni, "\n")

# ============================================================
# 7. Bivariate POT with censored likelihood (complete cases)
#    + univariate contributions from missings (E1M2, M1E2)
# ============================================================

# 7.1 Complete-case matrix for censored bivariate logistic POT
X_biv <- df_complete %>%
  transmute(
    log_SVL = log_SVL,
    log_TL  = log_TL
  ) %>%
  as.matrix()

n_biv <- nrow(X_biv)
cat("\nBivariate sample size (complete cases SVL & TL):", n_biv, "\n")

u_vec <- c(u1, u2)

quad_tab <- table(
  SVL_exceed = X_biv[, 1] > u1,
  TL_exceed  = X_biv[, 2] > u2
)
cat("\nQuadrant table (on log scale, complete cases only):\n")
print(quad_tab)

# 7.2 Missing patterns that contribute to the tail likelihood

# Only cases with an exceedance and the other margin truly missing:
y1_E_M2 <- df_svl_only$log_SVL[df_svl_only$log_SVL > u1]  # SVL exceeds, TL missing
y2_M1_E <- df_tl_only$log_TL[df_tl_only$log_TL > u2]      # TL exceeds, SVL missing

cat("\nCounts of missing-partner exceedances:\n")
cat("  SVL exceed, TL missing (E1M2):", length(y1_E_M2), "\n")
cat("  TL exceed, SVL missing (M1E2):", length(y2_M1_E), "\n")

# ---------------------------
# 7.3 Censored logistic log-likelihood using nllbvclog from bvpot.c
# ---------------------------

make_ll_cens_logistic_bvpot <- function(x, u, cshape = FALSE, cscale = FALSE) {
  # Use sep.bvdata from evd to build (x1, x2, nn, n, thid, lambda, ...)
  sep.bvdata <- getFromNamespace("sep.bvdata", "evd")
  spx        <- sep.bvdata(x = x, method = "cpot", u = u)
  
  function(theta) {
    scale1 <- theta["scale1"]
    shape1 <- theta["shape1"]
    scale2 <- theta["scale2"]
    shape2 <- theta["shape2"]
    dep    <- theta["dep"]
    
    if (isTRUE(cshape)) shape2 <- shape1
    if (isTRUE(cscale)) scale2 <- scale1
    
    # Call compiled C routine: nllbvclog(...)
    out <- .C("nllbvclog",
              data1  = as.double(spx$x1),
              data2  = as.double(spx$x2),
              nn     = as.integer(spx$nn),
              n      = as.integer(spx$n),
              thid   = as.double(spx$thdi),     # NOTE: double* in C
              lambda = as.double(spx$lambda),
              dep    = as.double(dep),
              scale1 = as.double(scale1),
              shape1 = as.double(shape1),
              scale2 = as.double(scale2),
              shape2 = as.double(shape2),
              dns    = as.double(0)
    )
    
    # C routine returns the *negative* log-likelihood in dns
    -out$dns
  }
}

ll_cens_log <- make_ll_cens_logistic_bvpot(
  x      = X_biv,
  u      = u_vec,
  cshape = FALSE,
  cscale = FALSE
)

# ---------------------------
# 7.4 Univariate GP contributions from missings (evd::dgpd)
# ---------------------------

lgp_exceed_evd <- function(y, u, scale, shape) {
  ex <- y - u
  evd::dgpd(ex, loc = 0, scale = scale, shape = shape, log = TRUE)
}

make_ll_missing_uni <- function(y1_E_M2, y2_M1_E, u1, u2) {
  function(theta) {
    scale1 <- theta["scale1"]
    shape1 <- theta["shape1"]
    scale2 <- theta["scale2"]
    shape2 <- theta["shape2"]
    
    ll1 <- if (length(y1_E_M2)) {
      sum(lgp_exceed_evd(y1_E_M2, u1, scale1, shape1))
    } else 0
    
    ll2 <- if (length(y2_M1_E)) {
      sum(lgp_exceed_evd(y2_M1_E, u2, scale2, shape2))
    } else 0
    
    ll1 + ll2
  }
}

ll_missing <- make_ll_missing_uni(y1_E_M2, y2_M1_E, u1, u2)

# Full log-likelihood: censored bivariate POT + univariate missings
loglik_full <- function(theta) {
  ll_cens_log(theta) + ll_missing(theta)
}

# Wrapper for optim
negloglik_full <- function(p, common_shape = FALSE) {
  if (!common_shape) {
    # Free shapes: theta = (scale1, shape1, scale2, shape2, dep)
    theta <- c(
      scale1 = p[1],
      shape1 = p[2],
      scale2 = p[3],
      shape2 = p[4],
      dep    = p[5]
    )
  } else {
    # H0: xi1 = xi2 = xi
    theta <- c(
      scale1 = p[1],
      shape1 = p[2],  # xi
      scale2 = p[3],
      shape2 = p[2],  # same xi
      dep    = p[4]
    )
  }
  -loglik_full(theta)
}

# ============================================================
# 8. Maximum likelihood: free shapes vs common shape
# ============================================================

# 8.1 Free shapes (H1)
start_free <- c(
  scale1 = sigma1_hat_uni,
  shape1 = xi1_hat_uni,
  scale2 = sigma2_hat_uni,
  shape2 = xi2_hat_uni,
  dep    = 0.8
)

opt_free <- optim(
  par     = unname(start_free),
  fn      = negloglik_full,
  common_shape = FALSE,
  method  = "Nelder-Mead"
)

par_free <- opt_free$par
names(par_free) <- names(start_free)

ell_hat_free <- -opt_free$value

cat("\nBivariate MLEs (full llh, free shapes):\n")
print(par_free)
cat("Maximized log-likelihood (free) =", ell_hat_free, "\n")

# 8.2 Common tail shape (H0: xi1 = xi2 = xi)
start_eq <- c(
  scale1 = sigma1_hat_uni,
  xi     = mean(c(xi1_hat_uni, xi2_hat_uni)),
  scale2 = sigma2_hat_uni,
  dep    = 0.8
)

opt_eq <- optim(
  par     = unname(start_eq),
  fn      = negloglik_full,
  common_shape = TRUE,
  method  = "Nelder-Mead"
)

par_eq_phi <- opt_eq$par
names(par_eq_phi) <- names(start_eq)

ell_hat_eq <- -opt_eq$value

# Reconstruct full θ under H0
par_eq <- c(
  scale1 = par_eq_phi["scale1"],
  shape1 = par_eq_phi["xi"],
  scale2 = par_eq_phi["scale2"],
  shape2 = par_eq_phi["xi"],
  dep    = par_eq_phi["dep"]
)

cat("\nCommon tail shape MLEs (full llh, H0: xi1 = xi2):\n")
print(par_eq)
cat("Maximized log-likelihood (H0) =", ell_hat_eq, "\n")

# Likelihood ratio test for H0: xi1 = xi2
LR_shape <- 2 * (ell_hat_free - ell_hat_eq)
df_LR    <- 1
p_LR     <- 1 - pchisq(LR_shape, df = df_LR)

cat("\nLikelihood ratio test for tail equality H0: xi1 = xi2\n")
cat("  LR statistic =", LR_shape, " with df =", df_LR, "\n")
cat("  p-value      =", p_LR, "\n")


# ============================================================
# 9. Endpoint estimation under common tail shape (H0)
#    + parametric bivariate bootstrap for CIs
# ============================================================

# 9.1 Point estimates of GP parameters and endpoints (log scale)

xi_eq_hat  <- unname(par_eq_phi["xi"])
sigma1_hat <- unname(par_eq_phi["scale1"])
sigma2_hat <- unname(par_eq_phi["scale2"])

ystar1_hat_log <- as.numeric(u1 - sigma1_hat / xi_eq_hat)  # SVL*
ystar2_hat_log <- as.numeric(u2 - sigma2_hat / xi_eq_hat)  # TL*

cat("\nLog-scale endpoint estimates (common-tail GP, H0):\n")
cat("  y*_SVL (log) =", ystar1_hat_log, "\n")
cat("  y*_TL  (log) =", ystar2_hat_log, "\n")

ystar1_hat_orig <- exp(ystar1_hat_log)
ystar2_hat_orig <- exp(ystar2_hat_log)

cat("\nPoint estimates for endpoints (original scale):\n")
cat("  SVL* ≈", ystar1_hat_orig, "\n")
cat("  TL*  ≈", ystar2_hat_orig, "\n")

# 9.2 Empirical tail fractions (used both for bootstrap & Stokes)

n_SVL_total <- sum(!is.na(df_log$log_SVL))
n_TL_total  <- sum(!is.na(df_log$log_TL))

n_SVL_exceed_u1 <- sum(df_log$log_SVL > u1, na.rm = TRUE)
n_TL_exceed_u2  <- sum(df_log$log_TL  > u2, na.rm = TRUE)

tail_frac_SVL <- n_SVL_exceed_u1 / n_SVL_total
tail_frac_TL  <- n_TL_exceed_u2  / n_TL_total

cat("\nEmpirical tail fractions:\n")
cat("  P(SVL > u1) ≈", tail_frac_SVL, "(u1 =", u1, ", exp(u1) ≈", exp(u1), "cm)\n")
cat("  P(TL  > u2) ≈", tail_frac_TL,  "(u2 =", u2, ", exp(u2) ≈", exp(u2), "cm)\n")

# 9.3 Parametric bivariate bootstrap (logistic GP, common tail shape)

# Sub-threshold samples for realism
sub_SVL_cc   <- df_complete$log_SVL[df_complete$log_SVL <= u1]
sub_TL_cc    <- df_complete$log_TL[df_complete$log_TL <= u2]
sub_SVL_only <- df_svl_only$log_SVL[df_svl_only$log_SVL <= u1]
sub_TL_only  <- df_tl_only$log_TL[df_tl_only$log_TL <= u2]

n_cc       <- nrow(df_complete)
n_svl_only <- nrow(df_svl_only)
n_tl_only  <- nrow(df_tl_only)

# Logistic EV parameter α from constrained fit
dep_name <- grep("^dep", names(par_eq), value = TRUE)[1]
alpha_logistic_hat <- unname(par_eq[dep_name])
theta_gumbel_hat   <- 1 / alpha_logistic_hat  # Gumbel copula θ ≥ 1

gumbel_cop_hat <- gumbelCopula(param = theta_gumbel_hat, dim = 2)

# Storage for bootstrap draws
boot_sigma1 <- numeric(B_boot)
boot_sigma2 <- numeric(B_boot)
boot_xi     <- numeric(B_boot)
boot_y1star <- numeric(B_boot)
boot_y2star <- numeric(B_boot)
boot_conv   <- logical(B_boot)

set.seed(2026)

for (b in seq_len(B_boot)) {
  ## (1) Simulate complete cases via copula + marginal tails
  U_cc <- rCopula(n_cc, gumbel_cop_hat)
  Y1_cc <- numeric(n_cc)
  Y2_cc <- numeric(n_cc)
  
  for (k in seq_len(n_cc)) {
    u1k <- U_cc[k, 1]
    u2k <- U_cc[k, 2]
    
    # SVL margin
    if (u1k <= 1 - tail_frac_SVL || length(sub_SVL_cc) == 0L) {
      Y1_cc[k] <- if (length(sub_SVL_cc)) sample(sub_SVL_cc, 1L, replace = TRUE) else u1
    } else {
      v1   <- (u1k - (1 - tail_frac_SVL)) / tail_frac_SVL
      exc1 <- evd::qgpd(v1, loc = 0, scale = sigma1_hat, shape = xi_eq_hat)
      Y1_cc[k] <- u1 + exc1
    }
    
    # TL margin
    if (u2k <= 1 - tail_frac_TL || length(sub_TL_cc) == 0L) {
      Y2_cc[k] <- if (length(sub_TL_cc)) sample(sub_TL_cc, 1L, replace = TRUE) else u2
    } else {
      v2   <- (u2k - (1 - tail_frac_TL)) / tail_frac_TL
      exc2 <- evd::qgpd(v2, loc = 0, scale = sigma2_hat, shape = xi_eq_hat)
      Y2_cc[k] <- u2 + exc2
    }
  }
  
  X_biv_boot <- cbind(Y1_cc, Y2_cc)
  
  ## (2) Simulate SVL-only and TL-only marginals
  if (n_svl_only > 0L) {
    U1_svl <- runif(n_svl_only)
    Y1_svl <- numeric(n_svl_only)
    for (k in seq_len(n_svl_only)) {
      u1k <- U1_svl[k]
      if (u1k <= 1 - tail_frac_SVL || length(sub_SVL_only) == 0L) {
        Y1_svl[k] <- if (length(sub_SVL_only)) sample(sub_SVL_only, 1L, replace = TRUE) else u1
      } else {
        v1   <- (u1k - (1 - tail_frac_SVL)) / tail_frac_SVL
        exc1 <- evd::qgpd(v1, loc = 0, scale = sigma1_hat, shape = xi_eq_hat)
        Y1_svl[k] <- u1 + exc1
      }
    }
    y1_E_M2_boot <- Y1_svl[Y1_svl > u1]
  } else {
    y1_E_M2_boot <- numeric(0L)
  }
  
  if (n_tl_only > 0L) {
    U2_tl <- runif(n_tl_only)
    Y2_tl <- numeric(n_tl_only)
    for (k in seq_len(n_tl_only)) {
      u2k <- U2_tl[k]
      if (u2k <= 1 - tail_frac_TL || length(sub_TL_only) == 0L) {
        Y2_tl[k] <- if (length(sub_TL_only)) sample(sub_TL_only, 1L, replace = TRUE) else u2
      } else {
        v2   <- (u2k - (1 - tail_frac_TL)) / tail_frac_TL
        exc2 <- evd::qgpd(v2, loc = 0, scale = sigma2_hat, shape = xi_eq_hat)
        Y2_tl[k] <- u2 + exc2
      }
    }
    y2_M1_E_boot <- Y2_tl[Y2_tl > u2]
  } else {
    y2_M1_E_boot <- numeric(0L)
  }
  
  ## (3) Build bootstrap likelihood and refit H0: ξ1 = ξ2
  ll_cens_log_boot <- make_ll_cens_logistic_bvpot(
    x      = X_biv_boot,
    u      = u_vec,
    cshape = FALSE,
    cscale = FALSE
  )
  ll_missing_boot <- make_ll_missing_uni(y1_E_M2_boot, y2_M1_E_boot, u1, u2)
  
  loglik_full_boot <- function(theta) {
    ll_cens_log_boot(theta) + ll_missing_boot(theta)
  }
  negloglik_full_boot <- function(p) {
    theta <- c(
      scale1 = p[1],
      shape1 = p[2],  # common ξ
      scale2 = p[3],
      shape2 = p[2],
      dep    = p[4]
    )
    -loglik_full_boot(theta)
  }
  
  start_boot <- c(
    scale1 = sigma1_hat,
    xi     = xi_eq_hat,
    scale2 = sigma2_hat,
    dep    = alpha_logistic_hat
  )
  
  opt_b <- try(
    optim(
      par    = unname(start_boot),
      fn     = negloglik_full_boot,
      method = "Nelder-Mead"
    ),
    silent = TRUE
  )
  
  if (inherits(opt_b, "try-error") || opt_b$convergence != 0) {
    boot_conv[b] <- FALSE
    next
  }
  
  boot_conv[b] <- TRUE
  phi_b <- opt_b$par
  names(phi_b) <- names(start_boot)
  
  sigma1_b <- phi_b["scale1"]
  xi_b     <- phi_b["xi"]
  sigma2_b <- phi_b["scale2"]
  
  boot_sigma1[b] <- sigma1_b
  boot_sigma2[b] <- sigma2_b
  boot_xi[b]     <- xi_b
  boot_y1star[b] <- u1 - sigma1_b / xi_b
  boot_y2star[b] <- u2 - sigma2_b / xi_b
}

cat("\nParametric bootstrap (common-tail logistic GP):\n")
cat("  Successful replicates:", sum(boot_conv), "out of", B_boot, "\n")

boot_ok <- boot_conv & is.finite(boot_xi) & (boot_xi < 0)

alpha_half <- (1 - ci_level) / 2

# Bootstrap CIs for parameters (still equal-tail)
ci_sigma1 <- quantile(boot_sigma1[boot_ok], probs = c(alpha_half, 1 - alpha_half), na.rm = TRUE)
ci_sigma2 <- quantile(boot_sigma2[boot_ok], probs = c(alpha_half, 1 - alpha_half), na.rm = TRUE)
ci_xi     <- quantile(boot_xi[boot_ok],     probs = c(alpha_half, 1 - alpha_half), na.rm = TRUE)

cat("\n", 100 * ci_level, "% bootstrap CIs for GP parameters (log scale):\n", sep = "")
cat("  sigma1:", ci_sigma1[1], "to", ci_sigma1[2], "\n")
cat("  sigma2:", ci_sigma2[1], "to", ci_sigma2[2], "\n")
cat("  xi    :", ci_xi[1],     "to", ci_xi[2],     "\n")

# One-sided (upper) bootstrap CIs for log endpoints
ci_y1_boot_log_upper <- quantile(boot_y1star[boot_ok], probs = ci_level, na.rm = TRUE)
ci_y2_boot_log_upper <- quantile(boot_y2star[boot_ok], probs = ci_level, na.rm = TRUE)

cat("\nOne-sided ", 100 * ci_level, "% upper bootstrap bounds for log endpoints:\n", sep = "")
cat("  y*_SVL (log) ≤ ", ci_y1_boot_log_upper, "\n", sep = "")
cat("  y*_TL  (log) ≤ ", ci_y2_boot_log_upper, "\n", sep = "")

# Back-transform to original scale
ci_y1_boot_orig_upper <- exp(ci_y1_boot_log_upper)
ci_y2_boot_orig_upper <- exp(ci_y2_boot_log_upper)

cat("\nOne-sided ", 100 * ci_level, "% upper bootstrap bounds for endpoints (original scale):\n", sep = "")
cat("  SVL* ≤ ", ci_y1_boot_orig_upper, "\n", sep = "")
cat("  TL*  ≤ ", ci_y2_boot_orig_upper, "\n", sep = "")

# ============================================================
# 10. Exceedance probability for the Stokes alligator
#      (SVL = 239 cm, TL = 450 cm)
# ============================================================

# Helper: GPD survival for y > u on the *log*-scale
gpd_survival <- function(y, u, sigma, xi) {
  x <- y - u
  if (x <= 0) return(1)  # not beyond threshold
  if (abs(xi) < 1e-8) {
    # Gumbel limit
    return(exp(-x / sigma))
  } else {
    term <- 1 + xi * x / sigma
    if (term <= 0) return(0)  # beyond fitted endpoint
    return(term^(-1 / xi))
  }
}

# 10.1 Stokes alligator measurements on log scale

SVL_stokes <- 239   # cm
TL_stokes  <- 450   # cm

log_SVL_stokes <- log(SVL_stokes)
log_TL_stokes  <- log(TL_stokes)

cat("\nStokes alligator (log scale):\n")
cat("  log(SVL) =", log_SVL_stokes, "\n")
cat("  log(TL)  =", log_TL_stokes,  "\n")

# 10.2 Univariate exceedance probabilities for Stokes

cond_SVL_stokes <- gpd_survival(
  y     = log_SVL_stokes,
  u     = u1,
  sigma = sigma1_hat,
  xi    = xi_eq_hat
)

cond_TL_stokes <- gpd_survival(
  y     = log_TL_stokes,
  u     = u2,
  sigma = sigma2_hat,
  xi    = xi_eq_hat
)

p_SVL_stokes <- tail_frac_SVL * cond_SVL_stokes
p_TL_stokes  <- tail_frac_TL  * cond_TL_stokes

cat("\nUnivariate predictive exceedance probabilities for Stokes:\n")
cat("  P(SVL > 239 cm) ≈", p_SVL_stokes,
    " → return period ≈", 1 / p_SVL_stokes, "individuals\n")
cat("  P(TL  > 450 cm) ≈", p_TL_stokes,
    " → return period ≈", 1 / p_TL_stokes,  "individuals\n")

# 10.3 Simple bivariate exceedance probability (independence)

p_joint_indep <- p_SVL_stokes * p_TL_stokes

cat("\nJoint exceedance under independence approximation:\n")
cat("  P(SVL > 239 cm, TL > 450 cm) ≈", p_joint_indep,
    " → return period ≈", 1 / p_joint_indep, "individuals\n")

# 10.4 Logistic EV copula (Gumbel) Monte Carlo for joint tail

alpha_logistic <- alpha_logistic_hat
theta_gumbel   <- 1 / alpha_logistic  # Gumbel parameter θ ≥ 1

cat("\nLogistic EV copula parameter (from POT fit):\n")
cat("  alpha (logistic dep) =", alpha_logistic, "\n")
cat("  theta (Gumbel copula) =", theta_gumbel, "\n")

# CDF thresholds for the Stokes alligator on [0,1]
u1_u <- 1 - tail_frac_SVL * cond_SVL_stokes  # F_SVL(239 cm)
u2_u <- 1 - tail_frac_TL  * cond_TL_stokes   # F_TL (450 cm)

cat("\nCDF thresholds for Stokes (on [0,1] scale):\n")
cat("  u1_u = F_SVL(239 cm) ≈", u1_u, "\n")
cat("  u2_u = F_TL (450 cm) ≈", u2_u, "\n")

if (!is.finite(u1_u) || !is.finite(u2_u) ||
    u1_u <= 0 || u1_u >= 1 ||
    u2_u <= 0 || u2_u >= 1) {
  stop("Invalid CDF thresholds u1_u / u2_u. Check tail fractions and GPD fits.")
}

set.seed(2025)
N_mc <- 2e5L

gumbel_cop <- gumbelCopula(param = theta_gumbel, dim = 2)
U_mc       <- rCopula(N_mc, gumbel_cop)

ind_joint_log <- (U_mc[, 1] > u1_u) & (U_mc[, 2] > u2_u)
p_joint_log   <- mean(ind_joint_log)
se_joint_log  <- sqrt(p_joint_log * (1 - p_joint_log) / N_mc)

cat("\nLogistic EV copula Monte Carlo estimate (bivariate tail):\n")
cat("  P(SVL > 239 cm, TL > 450 cm)\n")
cat("    ≈", p_joint_log,
    " (MC s.e. ≈", se_joint_log, ")\n")
cat("    → return period ≈", 1 / p_joint_log, "individuals\n")

cat("\nComparison independence vs logistic copula:\n")
cat("  Indep:  P ≈", p_joint_indep,
    " → RP ≈", 1 / p_joint_indep, "\n")
cat("  Logist: P ≈", p_joint_log,
    " → RP ≈", 1 / p_joint_log,  "\n")

N_year <- 7600  # approx mean annual harvest 2009–2024

RP_years_indep  <- (1 / p_joint_indep) / N_year
RP_years_logist <- (1 / p_joint_log)   / N_year

cat("Return period (independence) in years ≈", RP_years_indep,  "\n")
cat("Return period (logistic GP) in years ≈", RP_years_logist, "\n")

# ============================================================
# 11. Joint bootstrap distribution for (SVL*, TL*)
#      • 2D KDE + HPD contour
# ============================================================

# Keep only successful, finite bootstrap draws with xi < 0
boot_ok <- boot_conv & is.finite(boot_xi) & (boot_xi < 0)

if (!any(boot_ok)) {
  stop("No valid bootstrap draws for endpoints (boot_ok is empty).")
}

# Transform log-endpoints to original (cm) scale
svl_boot_cm <- exp(boot_y1star[boot_ok])
tl_boot_cm  <- exp(boot_y2star[boot_ok])

boot_joint <- data.frame(
  SVL = svl_boot_cm,
  TL  = tl_boot_cm
)

# ------------------------------------------------------------
# 11.1 Marginal one-sided 90% upper bounds (for comparison)
# ------------------------------------------------------------

q_upper <- ci_level

ub_svl_marg <- as.numeric(quantile(svl_boot_cm, probs = q_upper, na.rm = TRUE))
ub_tl_marg  <- as.numeric(quantile(tl_boot_cm,  probs = q_upper, na.rm = TRUE))

cat("\n", 100 * ci_level,
    "% one-sided marginal upper bounds for endpoints (original scale):\n", sep = "")
cat("  SVL* ≤ ", ub_svl_marg, " cm\n", sep = "")
cat("  TL*  ≤ ", ub_tl_marg,  " cm\n", sep = "")


# ------------------------------------------------------------
# 11.3 KDE-based HPD region (for plotting)
# ------------------------------------------------------------

# Reasonable plotting limits
svl_xlim <- quantile(svl_boot_cm, c(0.01, 0.99), na.rm = TRUE)
tl_ylim  <- quantile(tl_boot_cm,  c(0.01, 0.99), na.rm = TRUE)

# 2D KDE on a grid
kde_joint <- MASS::kde2d(
  x = boot_joint$SVL,
  y = boot_joint$TL,
  n = 200,
  lims = c(svl_xlim[1], svl_xlim[2], tl_ylim[1], tl_ylim[2])
)

# helper to evaluate kde2d at arbitrary points (bilinear interpolation)
eval_kde2d <- function(kde, x, y) {
  nx <- length(kde$x)
  ny <- length(kde$y)
  
  ix <- findInterval(x, kde$x, all.inside = TRUE)
  iy <- findInterval(y, kde$y, all.inside = TRUE)
  
  ix[ix >= nx] <- nx - 1
  iy[iy >= ny] <- ny - 1
  
  x1 <- kde$x[ix]
  x2 <- kde$x[ix + 1]
  y1 <- kde$y[iy]
  y2 <- kde$y[iy + 1]
  
  z11 <- kde$z[cbind(ix,     iy)]
  z21 <- kde$z[cbind(ix + 1, iy)]
  z12 <- kde$z[cbind(ix,     iy + 1)]
  z22 <- kde$z[cbind(ix + 1, iy + 1)]
  
  wx <- (x - x1) / (x2 - x1 + 1e-12)
  wy <- (y - y1) / (y2 - y1 + 1e-12)
  
  (1 - wx) * (1 - wy) * z11 +
    wx      * (1 - wy) * z21 +
    (1 - wx) * wy      * z12 +
    wx      * wy       * z22
}



# Data frame for KDE grid (for contour)
kde_df <- with(
  kde_joint,
  expand.grid(SVL = x, TL = y) |>
    transform(z = as.vector(z))
)

# ------------------------------------------------------------
# 11.4 Plot: 2D KDE + HPD contour + joint bounds
# ------------------------------------------------------------

p_joint <- ggplot() +
  # 2D KDE (filled levels)
  geom_density_2d_filled(
    data  = boot_joint,
    aes(x = SVL, y = TL),
    alpha = 0.7
  ) +
  # HPD contour from KDE
  geom_contour(
    data  = kde_df,
    aes(x = SVL, y = TL, z = z),
    breaks = hpd_threshold,
    colour = "black",
    linewidth = 1.0
  ) +
  # Bootstrap draws (faint)
  geom_point(
    data  = boot_joint,
    aes(x = SVL, y = TL),
    alpha = 0.25,
    size  = 0.6,
    color = "grey20"
  ) +
  # Point estimate (MAP / plug-in)
  geom_point(
    aes(x = map_svl, y = map_tl),
    color = "purple",
    size  = 3
  ) +
  # Stokes alligator (SVL, TL)
  geom_point(
    aes(x = SVL_stokes, y = TL_stokes),
    color = "red",
    size  = 3,
    shape = 4,
    stroke = 1.2
  ) +
  coord_cartesian(xlim = svl_xlim, ylim = tl_ylim) +
  labs(
    x     = "SVL endpoint (cm)",
    y     = "TL endpoint (cm)",
    fill  = "Density level"
  ) +
  theme_science_polished

p_joint

ggsave(
  file.path(FIG_DIR, "SVL_TL_joint_endpoint_bootstrap.png"),
  p_joint,
  dpi   = 600,
  w     = 6.5,
  h     = 4.5,
  units = "in"
)
