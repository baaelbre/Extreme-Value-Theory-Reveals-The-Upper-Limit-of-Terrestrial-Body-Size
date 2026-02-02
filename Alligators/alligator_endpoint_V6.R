# ============================================================
# Bivariate EVT (SVL, TL) — American alligators
# Endpoint-parameterised marginals + Kiriliouk-style preprocessing
# ============================================================

library(readxl)
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(extRemes)
library(scales)
library(glue)
library(forcats)
library(grid)
library(readr)

set.seed(42)

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
ci_level      <- 0.95
u_scan_lo_q   <- 0.75
u_scan_hi_q   <- 0.95
u_scan_n      <- 50
min_ex_mrl    <- 5
min_ex_fit    <- 20
topN_show     <- 10

trait_names   <- c("SVL","TL")
thresh_q_opt  <- c(SVL = 0.95, TL = 0.95)  # anchor quantiles (log scale)

# ============================================================
# 1. Ingest and basic cleaning
# ============================================================

DATA_XLSX <- "Data/experimental_alligator_harvest_woodward.xlsx"
df_raw    <- read_excel(DATA_XLSX)

stopifnot(all(c("SVL","TL","Deform") %in% names(df_raw)))

# Deform == 1 or 3: tail broken ⇒ TL missing; SVL kept.
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
  geom_text(aes(label = percent(completeness, accuracy=0.1)),
            vjust = -0.2, size = 4) +
  scale_y_continuous(labels = percent_format(accuracy=1), limits = c(0, 1.10)) +
  labs(title = "Completeness (SVL, TL)", x = NULL, y = "Completeness") +
  theme_science_polished

ggsave(file.path(FIG_DIR, "completeness_SVL_TL.png"), p_compl,
       dpi = 600, w = 6.0, h = 4.2, units = "in")

# ============================================================
# 3. Log-transform and thresholds (log scale)
# ============================================================

df_log <- df %>%
  mutate(
    log_SVL = ifelse(is.finite(SVL) & SVL > 0, log(SVL), NA_real_),
    log_TL  = ifelse(is.finite(TL)  & TL  > 0, log(TL),  NA_real_)
  )

u0_by_trait <- setNames(numeric(length(trait_names)), trait_names)
for (tr in trait_names) {
  v_log <- df_log[[paste0("log_", tr)]]
  u0_by_trait[tr] <- as.numeric(quantile(v_log, thresh_q_opt[tr], na.rm = TRUE))
}
print(u0_by_trait)
print(exp(u0_by_trait))  # original scale

# ============================================================
# 4. Univariate POT (log scale, endpoint parametrization)
# ============================================================

make_thresholds <- function(y, q_from = u_scan_lo_q, q_to = u_scan_hi_q, n = u_scan_n) {
  rng <- quantile(y, c(q_from, q_to), na.rm = TRUE)
  seq(rng[1], rng[2], length.out = n)
}

fit_gpd_at_u <- function(y, u) {
  fit <- fevd(y, type = "GP", threshold = u)
  par <- fit$results$par         # location (fixed = u), scale, shape
  # cov.theta is 2x2 for (scale, shape) when threshold fixed
  cov <- summary(fit)$cov.theta
  list(
    fit       = fit,
    scale     = unname(par["scale"]),
    shape     = unname(par["shape"]),
    cov       = cov,
    n_exceed  = sum(y > u)
  )
}

adj_scale_fun <- function(scale, shape, u, u0) {
  # σ(u0) = σ(u) - ξ (u - u0)
  scale - shape * (u - u0)
}

mrl_data <- function(y, u_seq, min_ex = 5) {
  sapply(u_seq, function(u) {
    ex <- y[y > u] - u
    if (length(ex) < min_ex) return(NA_real_)
    mean(ex)
  })
}

diagnostic_plots <- function(y, u, scale_hat, xi_hat, label) {
  exc <- y[y>u] - u
  n   <- length(exc)
  probs <- ppoints(n)
  
  # theoretical quantiles
  if (abs(xi_hat) > 1e-10) {
    theo_q <- u + scale_hat/xi_hat * (probs^(-xi_hat) - 1)
  } else {
    theo_q <- u - scale_hat * log(probs)
  }
  
  dfqq <- data.frame(
    Theoretical = rev(theo_q),
    Empirical   = sort(y[y>u])
  )
  pqq <- ggplot(dfqq, aes(Theoretical, Empirical)) +
    geom_point(color = "steelblue") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    labs(title = glue("Q–Q: {label}")) +
    theme_science_polished
  
  # P–P plot
  F_theo <- if (abs(xi_hat) > 1e-10) {
    1 - (1 + xi_hat * exc/scale_hat)^(-1/xi_hat)
  } else {
    1 - exp(-exc/scale_hat)
  }
  dfpp <- data.frame(
    Theoretical = sort(F_theo),
    Empirical   = (1:n)/n
  )
  ppp <- ggplot(dfpp, aes(Theoretical, Empirical)) +
    geom_point(color = "darkgreen") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    labs(title = glue("P–P: {label}")) +
    theme_science_polished
  
  # PIT histogram
  F_hat <- pevd(exc, scale = scale_hat, shape = xi_hat, type = "GP")
  pks <- ggplot(data.frame(F_hat = F_hat), aes(F_hat)) +
    geom_histogram(aes(y = ..density..), bins = 20,
                   fill = "skyblue", color = "black", alpha = 0.7) +
    geom_hline(yintercept = 1, color = "red", linetype = "dashed") +
    labs(title = glue("Uniformity (PIT): {label}")) +
    theme_science_polished
  
  list(pqq = pqq, ppp = ppp, pks = pks)
}

run_trait <- function(tr_key) {
  v_raw <- df[[tr_key]] |> suppressWarnings(as.numeric())
  y     <- log(v_raw[is.finite(v_raw) & v_raw > 0])
  ylab  <- glue("log {tr_key}")
  stopifnot(length(y) >= 40)
  
  q_anchor <- thresh_q_opt[tr_key]
  u_seq    <- make_thresholds(y)
  u0       <- quantile(y, q_anchor, na.rm = TRUE) |> as.numeric()
  
  # Mean Residual Life
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
  
  # Anchor GPD fit at u0
  fit0  <- fit_gpd_at_u(y, u0)
  xi0   <- fit0$shape
  sc0   <- fit0$scale
  cov0  <- fit0$cov
  z     <- xi0 / sqrt(cov0[2,2])     # Wald H0: xi = 0 vs H1: xi < 0
  p_wald <- pnorm(z)
  
  # Threshold stability for xi and adjusted scale at u0
  zcrit  <- qnorm(1 - (1 - ci_level)/2)
  
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
    ex <- y[y>u] - u
    if (length(ex) < min_ex_fit) next
    
    out <- try(fit_gpd_at_u(y, u), silent = TRUE)
    if (inherits(out, "try-error")) next
    
    se_shape <- sqrt(out$cov[2,2])
    df_scan$shape[i]    <- out$shape
    df_scan$shape_lo[i] <- out$shape - zcrit * se_shape
    df_scan$shape_hi[i] <- out$shape + zcrit * se_shape
    
    df_scan$adj[i]      <- adj_scale_fun(out$scale, out$shape, u, u0)
    var_adj <- out$cov[1,1] +
      (u - u0)^2 * out$cov[2,2] -
      2*(u - u0) * out$cov[1,2]
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
      y = "Shape (xi)",
      title = glue("Shape stability ({ylab})")
    ) +
    theme_science_polished
  
  p_adj <- ggplot(df_scan, aes(u, adj)) +
    geom_point(color = "blue") +
    geom_errorbar(aes(ymin = adj_lo, ymax = adj_hi),
                  width = 0.03, color = "blue") +
    geom_vline(xintercept = u0, color = "red", linetype = "dashed") +
    geom_hline(yintercept = sc0, color = "red", linetype = "dashed") +
    labs(
      x = glue("Threshold ({ylab})"),
      y = "Adjusted scale",
      title = glue("Adjusted-scale stability ({ylab})")
    ) +
    theme_science_polished
  
  ggsave(file.path(FIG_DIR, glue("{tr_key}_log_shape_stability.png")),
         p_shape, dpi = 600, w = 6.2, h = 4.2, units = "in")
  ggsave(file.path(FIG_DIR, glue("{tr_key}_log_adj_scale_stability.png")),
         p_adj,   dpi = 600, w = 6.2, h = 4.2, units = "in")
  
  # Diagnostics
  di <- diagnostic_plots(y, u0, sc0, xi0, ylab)
  ggsave(file.path(FIG_DIR, glue("{tr_key}_log_qq.png")),
         di$pqq, dpi = 600, w = 6.2, h = 4.2, units = "in")
  ggsave(file.path(FIG_DIR, glue("{tr_key}_log_pp.png")),
         di$ppp, dpi = 600, w = 6.2, h = 4.2, units = "in")
  ggsave(file.path(FIG_DIR, glue("{tr_key}_log_pit_uniformity.png")),
         di$pks, dpi = 600, w = 6.2, h = 4.2, units = "in")
  
  # Weibull endpoint on log-scale: y* = u0 - σ/xi (xi<0)
  y_star_hat <- if (xi0 < 0) u0 - sc0/xi0 else NA_real_
  
  tibble(
    trait          = tr_key,
    u0_quantile    = q_anchor,
    u0             = u0,
    xi_hat         = xi0,
    sigma_hat      = sc0,
    y_star_hat_log = y_star_hat,
    wald_p_xi_lt_0 = p_wald,
    n_exceed       = fit0$n_exceed,
    n_total_used   = length(y)
  )
}

svl_tbl <- run_trait("SVL")
tl_tbl  <- run_trait("TL")

marginals_tbl <- bind_rows(svl_tbl, tl_tbl)
write_csv(marginals_tbl, file.path(FIG_DIR, "pot_summary_SVL_TL_log.csv"))
print(marginals_tbl)

# ============================================================
# 5. Pairwise tail regimes & missing/observed/exceed/censored
# ============================================================

# Precompute log vectors and thresholds
x1 <- df_log$log_SVL
x2 <- df_log$log_TL
u1 <- u0_by_trait["SVL"]
u2 <- u0_by_trait["TL"]

# Observed / missing
O1 <- is.finite(x1); M1 <- !O1
O2 <- is.finite(x2); M2 <- !O2

# Exceeding / below threshold (censored) among observed
E1 <- O1 & (x1 > u1)
E2 <- O2 & (x2 > u2)
B1 <- O1 & !E1       # observed but below u1 ⇒ censored at u1
B2 <- O2 & !E2       # observed but below u2 ⇒ censored at u2

# Keep a table with the full pattern per specimen (for later use)
pattern_df <- tibble(
  specimen = df$specimen,
  log_SVL  = x1,
  log_TL   = x2,
  O1, O2, M1, M2, E1, E2, B1, B2,
  pattern = case_when(
    E1 & E2           ~ "E1,E2   (both exceed)",
    E1 & B2           ~ "E1,B2   (SVL exceed, TL censored)",
    B1 & E2           ~ "B1,E2   (SVL censored, TL exceed)",
    E1 & M2           ~ "E1,M2   (SVL exceed, TL missing)",
    M1 & E2           ~ "M1,E2   (SVL missing, TL exceed)",
    B1 & B2           ~ "B1,B2   (both censored)",
    B1 & M2           ~ "B1,M2   (SVL censored, TL missing)",
    M1 & B2           ~ "M1,B2   (SVL missing, TL censored)",
    TRUE              ~ "M1,M2   (both missing)"
  )
)

# Quick scatter illustration: dots = both observed, lines = one missing
make_pair_plot <- function(t1, t2) {
  fn <- file.path(FIG_DIR, glue("pair_scatter_{t1}_{t2}_final.png"))
  
  xr <- range(x1[O1], na.rm = TRUE)
  yr <- range(x2[O2], na.rm = TRUE)
  pad <- function(r, p = 0.03) {
    w <- diff(r); c(r[1] - p * w, r[2] + p * w)
  }
  xr <- pad(xr); yr <- pad(yr)
  
  both_obs <- O1 & O2
  
  dd_dot <- tibble(
    specimen = df$specimen[both_obs],
    x        = x1[both_obs],
    y        = x2[both_obs],
    kind     = case_when(
      E1[both_obs] & E2[both_obs] ~ "(E1,E2) •",
      E1[both_obs] & B2[both_obs] ~ "(E1,¬E2) •",
      B1[both_obs] & E2[both_obs] ~ "(¬E1,E2) •",
      TRUE                        ~ "none •"
    )
  )
  
  dd_v <- tibble(
    specimen = df$specimen[E1 & M2],
    x  = x1[E1 & M2],
    y0 = yr[1],
    y1 = yr[2],
    kind = "(E1,M2) |"
  )
  
  dd_h <- tibble(
    specimen = df$specimen[M1 & E2],
    y  = x2[M1 & E2],
    x0 = xr[1],
    x1 = xr[2],
    kind = "(M1,E2) —"
  )
  
  col_map <- c(
    "(E1,E2) •" = "#DC2626",
    "(E1,¬E2) •"= "#2563EB",
    "(¬E1,E2) •"= "#059669",
    "none •"    = "grey70",
    "(E1,M2) |" = "#2563EB",
    "(M1,E2) —" = "#059669"
  )
  
  p <- ggplot() +
    geom_point(data = dd_dot, aes(x = x, y = y, color = kind),
               size = 2.6, alpha = 0.95) +
    geom_segment(data = dd_v,
                 aes(x = x, xend = x, y = y0, yend = y1, color = kind),
                 linewidth = 0.9, alpha = 0.95) +
    geom_segment(data = dd_h,
                 aes(x = x0, xend = x1, y = y, yend = y, color = kind),
                 linewidth = 0.9, alpha = 0.95) +
    geom_vline(xintercept = u1, linetype = "dashed", color = "red") +
    geom_hline(yintercept = u2, linetype = "dashed", color = "red") +
    scale_color_manual(values = col_map, name = "Regime") +
    coord_cartesian(xlim = xr, ylim = yr, expand = FALSE) +
    labs(
      x = glue("{t1} (log)"),
      y = glue("{t2} (log)")
    ) +
    theme_science_polished
  print(p)
  
  ggsave(fn, p, dpi = 600, w = 6.8, h = 5.6, units = "in")
  message("Saved: ", normalizePath(fn))
}

make_pair_plot("SVL", "TL")

# ============================================================
# 6. Empirical χ(q) for SVL–TL (using both observed)
# ============================================================

emp_chi <- function(Ymat, q) {
  U <- apply(Ymat, 2, function(v) {
    r <- rank(v, na.last = "keep", ties.method = "average")
    r / sum(is.finite(v))
  })
  ok_row <- apply(U, 1, function(z) all(!is.na(z)))
  if (!any(ok_row)) return(NA_real_)
  
  ind <- rep(TRUE, nrow(U))
  for (j in seq_len(ncol(U))) {
    ind <- ind & (U[, j] > q)
  }
  sum(ind & ok_row) / (sum(ok_row) * (1 - q))
}

log_mat <- as.matrix(df_log %>% transmute(SVL = log_SVL, TL = log_TL))
q_grid  <- seq(0.80, 0.99, by = 0.01)

chi_pairs <- tibble(
  pair = "SVL-TL",
  q    = rep(q_grid, each = 1),
  chi  = map_dbl(q_grid, ~ emp_chi(log_mat, .x))
)

pick_qstar <- function(dfq) {
  chi99 <- dfq$chi[dfq$q == max(dfq$q, na.rm = TRUE)]
  cand  <- dfq %>% filter(q >= 0.90, is.finite(chi)) %>% arrange(q)
  if (!length(chi99) || !nrow(cand)) return(0.95)
  hit <- cand %>% filter(abs(chi - chi99) <= 0.05 * abs(chi99))
  if (nrow(hit)) hit$q[1] else 0.95
}

q_star_pairs <- chi_pairs %>%
  group_by(pair) %>%
  summarize(q_star = pick_qstar(cur_data_all()), .groups = "drop")

chi_pairs_hat <- chi_pairs %>%
  inner_join(q_star_pairs, by = "pair") %>%
  filter(q == q_star) %>%
  transmute(scope = pair, q_star, chi_hat = chi)

write_csv(chi_pairs_hat, file.path(FIG_DIR, "kiriliouk_chi_dependence_SVL_TL.csv"))

p_pairs <- chi_pairs %>%
  inner_join(q_star_pairs, by = "pair") %>%
  ggplot(aes(q, chi)) +
  geom_line() + geom_point() +
  geom_vline(aes(xintercept = q_star),
             linetype = "dashed", color = "red") +
  facet_wrap(~pair, ncol = 1, scales = "free_y") +
  labs(
    title = expression(paste("Empirical pairwise ", chi[i*j], "(q) — plateau & ", q["*"])),
    y     = expression(chi),
    x     = "q"
  ) +
  theme_science_polished

ggsave(file.path(FIG_DIR, "chi_pairs_SVL_TL_plateau.png"),
       p_pairs, dpi = 600, w = 6.2, h = 4.6, units = "in")

# ------------------------------------------------------------
# Compact marginal outputs
# ------------------------------------------------------------
write_csv(
  marginals_tbl %>%
    select(trait, u0, xi_hat, sigma_hat, y_star_hat_log,
           wald_p_xi_lt_0, n_exceed, n_total_used),
  file.path(FIG_DIR, "marginals_POT_SVL_TL_log.csv")
)

message("Marginal POT done. Outputs in: ", normalizePath(FIG_DIR))

# ============================================================
# 7. Endpoint-based unit-Pareto standardisation (for spectral fit)
# ============================================================

# Helper: build σ(u) from (y*, xi) in Weibull domain
sigma_from_endpoint <- function(u, y_star, xi) {
  # y* = u - σ/xi ⇒ σ = (u - y*) * xi  (xi < 0, u - y* < 0 ⇒ σ > 0)
  (u - y_star) * xi
}

gp_tail_cdf_endpoint <- function(y, u, y_star, xi) {
  # Conditional GP CDF on excesses z = y - u > 0, in endpoint parametrisation
  z <- y - u
  out <- rep(NA_real_, length(z))
  
  finite <- is.finite(z) & is.finite(u) & is.finite(y_star) & is.finite(xi)
  posz   <- z > 0
  ok     <- finite & posz & (xi != 0)
  
  ok[is.na(ok)] <- FALSE
  if (any(ok)) {
    sig <- sigma_from_endpoint(u, y_star, xi)
    supp <- 1 + xi * z[ok] / sig > 0
    ok2  <- ok
    ok2[ok] <- supp
    if (any(ok2)) {
      sig2 <- sigma_from_endpoint(u, y_star, xi)
      out[ok2] <- 1 - (1 + xi * z[ok2] / sig2)^(-1/xi)
    }
  }
  out
}

to_unit_pareto_endpoint <- function(y, u, y_star, xi) {
  # T >= 1 for exceedances; NA otherwise
  Fz <- gp_tail_cdf_endpoint(y, u, y_star, xi)
  T  <- rep(NA_real_, length(y))
  ok <- is.finite(Fz)
  ok[is.na(ok)] <- FALSE
  if (any(ok)) {
    T[ok] <- 1 / pmax(1 - Fz[ok], .Machine$double.eps)
  }
  T
}

# Marginal endpoint/xi estimates at u0 (log scale)
xi_hat_uni  <- c(
  SVL = marginals_tbl$xi_hat[marginals_tbl$trait == "SVL"],
  TL  = marginals_tbl$xi_hat[marginals_tbl$trait == "TL"]
)
ystar_hat   <- c(
  SVL = marginals_tbl$y_star_hat_log[marginals_tbl$trait == "SVL"],
  TL  = marginals_tbl$y_star_hat_log[marginals_tbl$trait == "TL"]
)

# Exceedance indicators at thresholds u0 (log scale)
exceed_SVL <- O1 & (x1 > u1)
exceed_TL  <- O2 & (x2 > u2)
both_ex    <- exceed_SVL & exceed_TL  # *only* these contribute to W for spectral fit

# Endpoint-based unit Pareto standardization
T_SVL <- to_unit_pareto_endpoint(x1, u1, ystar_hat["SVL"], xi_hat_uni["SVL"])
T_TL  <- to_unit_pareto_endpoint(x2, u2, ystar_hat["TL"],  xi_hat_uni["TL"])

W <- T_SVL[both_ex] / (T_SVL[both_ex] + T_TL[both_ex])
W <- W[is.finite(W) & W > 0 & W < 1]
nW <- length(W)
message("Both-exceedances used for spectral fit (E1,E2, observed): n = ", nW)

# ============================================================
# 8. Symmetric logistic spectral model on W (angles)
# ============================================================

# NOTE: this is the same functional form you used before; here we only
# rewrap it. We are effectively fitting the parametric spectral density
# h_alpha(w) on (0,1) using only E1,E2 (no missing, no censored margin).

h_logistic <- function(w, alpha) {
  w <- pmin(pmax(w, 1e-12), 1 - 1e-12)
  a <- alpha
  u <- w^(1/a)
  v <- (1 - w)^(1/a)
  S <- u + v
  # Standard symmetric logistic spectral density
  num <- (1/a) * (w^((1 - a)/a) + (1 - w)^((1 - a)/a)) * S^(a - 2)
  den <- w * (1 - w)
  out <- num / den
  out[!is.finite(out)] <- .Machine$double.eps
  out
}

negloglik_alpha <- function(eta, w) {
  alpha <- 1 / (1 + exp(-eta))  # map R -> (0,1)
  -sum(log(h_logistic(w, alpha)))
}

eta0      <- qlogis(0.7)  # moderate dep start
fit_alpha <- optim(eta0, negloglik_alpha, w = W, method = "BFGS")
alpha_hat <- 1 / (1 + exp(-fit_alpha$par))
logLik_spec <- -fit_alpha$value

# Stable tail dependence function A(1/2) and χ for logistic
A_half <- 2^(alpha_hat - 1)    # A(1/2)
chi_hat <- 2 - 2^(alpha_hat)   # χ = 2 - 2^α

spec_tbl <- tibble(
  model          = "symmetric_logistic",
  alpha          = as.numeric(alpha_hat),
  A_half         = as.numeric(A_half),
  chi_hat        = as.numeric(chi_hat),
  n_both_exceed  = nW,
  logLik_spectral= as.numeric(logLik_spec)
)
write_csv(spec_tbl, file.path(FIG_DIR, "spectral_fit_logistic_SVL_TL.csv"))
print(spec_tbl)

# Quick plot of W with fitted h_alpha
wgrid <- seq(1e-3, 1 - 1e-3, length.out = 400)
df_w  <- tibble(w = W)

p_w <- ggplot(df_w, aes(w)) +
  geom_histogram(aes(y = ..density..), bins = 30,
                 fill = "#93C5FD", color = "black", alpha = 0.8) +
  geom_line(
    data = tibble(w = wgrid, d = h_logistic(wgrid, alpha_hat)),
    aes(w, d),
    linewidth = 1
  ) +
  labs(
    title = glue("Spectral angles W (E1,E2, both observed) with logistic h(w), alpha={round(alpha_hat,3)}"),
    x     = "W = T_SVL / (T_SVL + T_TL)",
    y     = "Density"
  ) +
  theme_science_polished
p_w

ggsave(file.path(FIG_DIR, "spectral_angles_logistic_fit.png"),
       p_w, dpi = 600, w = 6.6, h = 4.6, units = "in")

# ============================================================
# 9. Common-tail test H0: xi_SVL = xi_TL (log-scale POT)
# ============================================================

gp_loglik <- function(y, u, sigma, xi){
  # drop missing values
  y <- y[is.finite(y)]
  
  # basic parameter checks
  if (!is.finite(sigma) || !is.finite(xi) || sigma <= 0) return(-Inf)
  
  z <- y[y > u] - u
  
  # no exceedances → undefined loglik for our purposes
  if (length(z) == 0) return(-Inf)
  
  # support condition: 1 + xi z / sigma > 0, done NA-safely
  supp <- 1 + xi * z / sigma
  if (any(supp <= 0, na.rm = TRUE)) return(-Inf)
  
  -length(z) * log(sigma) - (1/xi + 1) * sum(log(supp), na.rm = TRUE)
}


# Unconstrained: treat (sigma_j, xi_j) separate, use POT estimates
sig_hat_uni <- c(
  SVL = marginals_tbl$sigma_hat[marginals_tbl$trait == "SVL"],
  TL  = marginals_tbl$sigma_hat[marginals_tbl$trait == "TL"]
)

L_uc <- 0
for (tr in trait_names) {
  v_log <- df_log[[paste0("log_", tr)]]
  u     <- u0_by_trait[tr]
  s     <- sig_hat_uni[tr]
  x     <- xi_hat_uni[tr]
  L_uc  <- L_uc + gp_loglik(v_log, u, s, x)
}

# Constrained: xi common, reoptimise (sigma_SVL, sigma_TL, xi)
par_start_c <- c(
  log(sig_hat_uni["SVL"]),
  log(sig_hat_uni["TL"]),
  mean(xi_hat_uni)
)

LL_c <- function(p) {
  s1 <- exp(p[1])
  s2 <- exp(p[2])
  x  <- p[3]
  
  L1 <- gp_loglik(df_log$log_SVL, u0_by_trait["SVL"], s1, x)
  L2 <- gp_loglik(df_log$log_TL,  u0_by_trait["TL"],  s2, x)
  
  Lsum <- L1 + L2
  
  # If any piece is non-finite, return a large penalty (we minimize LL_c)
  if (!is.finite(L1) || !is.finite(L2) || !is.finite(Lsum)) {
    return(1e10)   # big positive value; optimizer will move away
  }
  
  -(Lsum)          # Negative total loglik (for minimization)
}


opt_c <- optim(
  par_start_c, LL_c,
  lower  = c(log(1e-6), log(1e-6), -5),
  upper  = c(log(1e+6), log(1e+6), 1 - 1e-6),
  method = "L-BFGS-B"
)

L_c  <- -opt_c$value
lrt  <- 2 * (L_uc - L_c)
pval <- pchisq(lrt, df = 1, lower.tail = FALSE)

common_tail_tbl <- tibble(
  test          = "H0: xi_SVL = xi_TL (log scale GP)",
  L_uc          = L_uc,
  L_c           = L_c,
  LRT           = lrt,
  df            = 1,
  pvalue        = pval,
  xi_common_mle = opt_c$par[3],
  sigma_SVL_mle = exp(opt_c$par[1]),
  sigma_TL_mle  = exp(opt_c$par[2])
)
write_csv(common_tail_tbl, file.path(FIG_DIR, "lrt_common_tail_SVL_TL.csv"))
print(common_tail_tbl)

message("All done. See ", normalizePath(FIG_DIR))
