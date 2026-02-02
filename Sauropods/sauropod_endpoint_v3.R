# =========================
# Multivariate EVT (per-trait) — Sauropod measurements
# Data columns (exact): 'genus and species', 'humerus L', 'femur L',
#                       'humerus circ (mm)', 'femur circ (mm)'
# Traits analyzed on log scale: LF, LH, CF, CH
# =========================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(extRemes)   # fevd, pevd
  library(purrr)
  library(scales)
  library(glue)
  library(forcats)
  library(grid)
})

set.seed(42)

# ---------------------------
# Directories & theming
# ---------------------------
FIG_DIR <- "Figures_v6"
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
# Settings
# ---------------------------
ci_level      <- 0.95
u_scan_lo_q   <- 0.60     # threshold scan lower quantile
u_scan_hi_q   <- 0.90     # threshold scan upper quantile
u_scan_n      <- 50
min_ex_mrl    <- 5
min_ex_fit    <- 10
topN_show     <- 10

# >>> Trait-specific anchor quantiles (edit as desired)
thresh_q_opt <- c(
  LF = 0.78,   # log femur length
  LH = 0.78,   # log humerus length
  CF = 0.78,   # log femur circumference
  CH = 0.78    # log humerus circumference
)

# ---------------------------
# Ingest (exact column names)
# ---------------------------
DATA_XLSX <- "Data/sauropod_measurements_demic.xlsx"
df_raw <- read_excel(DATA_XLSX)

# Standardize names we need
stopifnot(all(c("genus and species",
                "humerus L", "femur L",
                "humerus circ (mm)", "femur circ (mm)") %in% names(df_raw)))

df <- df_raw %>%
  transmute(
    specimen = `genus and species`,
    LF  = suppressWarnings(as.numeric(`femur L`)),
    CF  = suppressWarnings(as.numeric(`femur circ (mm)`)),
    CH  = suppressWarnings(as.numeric(`humerus circ (mm)`)),
    LH  = suppressWarnings(as.numeric(`humerus L`))
  )

trait_names <- c("LF","CF","CH","LH")

# ---------------------------
# Completeness per trait
# ---------------------------
compl_tbl <- trait_names %>%
  set_names() %>%
  map_df(function(tr) {
    v <- df[[tr]]
    data.frame(
      trait = tr,
      n_total = length(v),
      n_obs   = sum(!is.na(v)),
      completeness = mean(!is.na(v))
    )
  })

print(compl_tbl)

p_compl <- compl_tbl %>%
  mutate(trait = fct_inorder(trait)) %>%
  ggplot(aes(trait, completeness)) +
  geom_col(fill="#3B82F6") +
  geom_text(aes(label=scales::percent(completeness, accuracy=0.1)),
            vjust=-0.2, size=4) +
  scale_y_continuous(labels=percent_format(accuracy=1), limits=c(0,1.10)) +
  labs(title="Completeness per trait", x=NULL, y="Completeness") +
  theme_science_polished
p_compl

ggsave(file.path(FIG_DIR, "completeness_per_trait.png"),
       p_compl, dpi=600, w=6.5, h=4.5, units="in")
write.csv(compl_tbl, file.path(FIG_DIR, "completeness_per_trait.csv"), row.names = FALSE)

# ---------------------------
# Largest specimens per trait (Top N)
# ---------------------------
largest_by_trait <- function(tr, topN=10) {
  df %>%
    filter(!is.na(.data[[tr]])) %>%
    arrange(desc(.data[[tr]])) %>%
    mutate(rank = row_number()) %>%
    slice_head(n = topN) %>%
    transmute(trait = tr,
              rank,
              value_raw = .data[[tr]],
              specimen)
}

largest_list <- map_df(trait_names, largest_by_trait, topN = topN_show)
print(largest_list)
write.csv(largest_list,
          file.path(FIG_DIR, "largest_specimens_per_trait_top10.csv"),
          row.names = FALSE)

plot_topN <- function(tr) {
  dd <- largest_by_trait(tr, topN = topN_show)
  ggplot(dd, aes(x = reorder(specimen, value_raw), y = value_raw)) +
    geom_col(fill="#10B981") +
    coord_flip() +
    labs(title = glue("Top {topN_show} largest: {tr}"),
         x = "Specimen", y = glue("{tr} (mm)")) +
    theme_science_polished
}
walk(trait_names, ~ ggsave(file.path(FIG_DIR, glue("topN_{.x}.png")),
                           plot_topN(.x), dpi=600, w=6.5, h=5, units="in"))
# ---------------------------
# EVT helpers
# ---------------------------
make_thresholds <- function(y, q_from=u_scan_lo_q, q_to=u_scan_hi_q, n=u_scan_n){
  rng <- quantile(y, c(q_from, q_to), na.rm=TRUE)
  seq(rng[1], rng[2], length.out=n)
}
fit_gpd_at_u <- function(y, u) {
  fit <- fevd(y, type = "GP", threshold = u)
  par <- fit$results$par
  list(fit=fit,
       scale=unname(par["scale"]),
       shape=unname(par["shape"]),
       cov=summary(fit)$cov.theta)
}
adj_scale_fun <- function(scale, shape, u, u0) scale - shape*(u - u0)

mrl_data <- function(y, u_seq, min_ex=5){
  sapply(u_seq, function(u){
    ex <- y[y>u] - u
    if (length(ex) < min_ex) return(NA_real_)
    mean(ex)
  })
}

diagnostic_plots <- function(y, u, scale_hat, shape_hat, label){
  exc <- y[y>u] - u; n <- length(exc)
  probs <- ppoints(n)
  # Theoretical quantiles on original (not excess) scale
  theo_q <- if (abs(shape_hat) > 1e-10) {
    u + scale_hat/shape_hat * (probs^(-shape_hat) - 1)
  } else {
    u - scale_hat * log(probs)
  }
  dfqq <- data.frame(Theoretical = rev(theo_q), Empirical = sort(y[y>u]))
  pqq <- ggplot(dfqq, aes(Theoretical, Empirical)) +
    geom_point(color="steelblue") +
    geom_abline(slope=1, intercept=0, linetype="dashed", color="red") +
    labs(title = glue("Q–Q: {label}"),
         x="Theoretical quantiles", y="Empirical quantiles") +
    theme_science_polished
  pqq
  
  F_theo <- if (abs(shape_hat)>1e-10) {
    1 - (1 + shape_hat*exc/scale_hat)^(-1/shape_hat)
  } else {
    1 - exp(-exc/scale_hat)
  }
  dfpp <- data.frame(Theoretical=sort(F_theo), Empirical=(1:n)/n)
  ppp <- ggplot(dfpp, aes(Theoretical, Empirical)) +
    geom_point(color="darkgreen") +
    geom_abline(slope=1, intercept=0, linetype="dashed", color="red") +
    labs(title = glue("P–P: {label}"),
         x="Theoretical CDF", y="Empirical CDF") +
    theme_science_polished
  ppp
  
  F_hat <- pevd(exc, scale=scale_hat, shape=shape_hat, type="GP")
  pks <- ggplot(data.frame(F_hat=F_hat), aes(F_hat)) +
    geom_histogram(aes(y=..density..), bins=20, fill="skyblue", color="black", alpha=0.7) +
    geom_hline(yintercept=1, color="red", linetype="dashed") +
    geom_density(color="darkblue", linewidth=1.1, adjust=1.5) +
    labs(title=glue("Uniformity (PIT): {label}"),
         x=expression(hat(F)(y)), y="Density") +
    theme_science_polished
  
  list(pqq=pqq, ppp=ppp, pks=pks)
}

# ---------------------------
# Master routine per trait (uses trait-specific u0 quantile)
# ---------------------------
run_trait <- function(tr_key) {
  v_raw <- df[[tr_key]]
  v_raw <- suppressWarnings(as.numeric(v_raw))
  nn <- sum(!is.na(v_raw))
  if (nn < 80) {
    message(glue("Trait {tr_key}: only {nn} non-missing; proceeding but stability may be weak."))
  }
  # Log scale (drop nonpositive)
  y <- log(v_raw[is.finite(v_raw) & v_raw > 0])
  ylab <- glue("log {tr_key}")
  stopifnot(length(y) >= 40)
  
  # Trait-specific anchor quantile
  q_anchor <- if (!is.null(thresh_q_opt[tr_key])) as.numeric(thresh_q_opt[tr_key]) else 0.78
  
  # Thresholds
  u_seq <- make_thresholds(y, q_from=u_scan_lo_q, q_to=u_scan_hi_q, n=u_scan_n)
  u0    <- as.numeric(quantile(y, q_anchor, na.rm=TRUE))
  
  # MRL
  mrl_vals <- mrl_data(y, u_seq, min_ex=min_ex_mrl)
  p_mrl <- ggplot(data.frame(u=u_seq, mrl=mrl_vals), aes(u, mrl)) +
    geom_line() + geom_point() +
    geom_vline(xintercept = u0, linetype = "dashed", color = "red") +
    labs(title = glue("MRL plot ({ylab})"),
         x = glue("Threshold ({ylab})"), y = "Mean excess") +
    theme_science_polished
  ggsave(file.path(FIG_DIR, glue("{tr_key}_log_mrl.png")),
         p_mrl, dpi=600, w=7, h=5, units="in")
  
  # Anchor fit at u0
  fit0 <- fit_gpd_at_u(y, u0)
  shape_hat0 <- fit0$shape
  adj_hat0   <- fit0$scale
  
  shape_se0  <- sqrt(fit0$cov[2,2])
  z          <- shape_hat0 / shape_se0
  p_wald     <- pnorm(z)              # one-sided H1: xi < 0
  
  # Stability scans
  zcrit <- qnorm(1 - (1 - ci_level)/2)
  shape <- scale <- shape_lo <- shape_hi <- adj <- adj_lo <- adj_hi <- rep(NA_real_, length(u_seq))
  
  for (i in seq_along(u_seq)) {
    u <- u_seq[i]
    ex <- y[y > u] - u
    if (length(ex) < min_ex_fit) next
    out <- try(fit_gpd_at_u(y, u), silent=TRUE)
    if (inherits(out,"try-error")) next
    scale[i] <- out$scale; shape[i] <- out$shape
    se_scale <- sqrt(out$cov[1,1]); se_shape <- sqrt(out$cov[2,2])
    shape_lo[i] <- shape[i] - zcrit*se_shape
    shape_hi[i] <- shape[i] + zcrit*se_shape
    adj[i] <- adj_scale_fun(scale[i], shape[i], u, u0)
    var_adj <- out$cov[1,1] + (u - u0)^2*out$cov[2,2] - 2*(u - u0)*out$cov[1,2]
    se_adj  <- sqrt(max(var_adj,0))
    adj_lo[i] <- adj[i] - zcrit*se_adj
    adj_hi[i] <- adj[i] + zcrit*se_adj
  }
  
  p_shape <- ggplot(data.frame(u=u_seq, shape, shape_lo, shape_hi), aes(u, shape)) +
    geom_point() +
    geom_errorbar(aes(ymin=shape_lo, ymax=shape_hi), width=0.03, color="blue") +
    geom_vline(xintercept=u0, color="red", linetype="dashed") +
    geom_hline(yintercept=shape_hat0, color="red", linetype="dashed") +
    labs(x=glue("Threshold ({ylab})"), y="Shape (xi)",
         title = glue("Shape stability ({ylab})")) +
    theme_science_polished
  p_shape
  
  p_adj <- ggplot(data.frame(u=u_seq, adj, adj_lo, adj_hi), aes(u, adj)) +
    geom_point() +
    geom_errorbar(aes(ymin=adj_lo, ymax=adj_hi), width=0.03, color="blue") +
    geom_vline(xintercept=u0, color="red", linetype="dashed") +
    geom_hline(yintercept=adj_hat0, color="red", linetype="dashed") +
    labs(x=glue("Threshold ({ylab})"), y="Adjusted scale",
         title = glue("Adjusted-scale stability ({ylab})")) +
    theme_science_polished
  p_adj
  
  ggsave(file.path(FIG_DIR, glue("{tr_key}_log_shape_stability.png")),
         p_shape, dpi=600, w=7, h=5, units="in")
  ggsave(file.path(FIG_DIR, glue("{tr_key}_log_adj_scale_stability.png")),
         p_adj,   dpi=600, w=7, h=5, units="in")
  
  # Diagnostics at u0
  diags <- diagnostic_plots(y, u0, fit0$scale, fit0$shape, ylab)
  ggsave(file.path(FIG_DIR, glue("{tr_key}_log_qq.png")),  diags$pqq, dpi=600, w=7, h=5, units="in")
  ggsave(file.path(FIG_DIR, glue("{tr_key}_log_pp.png")),  diags$ppp, dpi=600, w=7, h=5, units="in")
  ggsave(file.path(FIG_DIR, glue("{tr_key}_log_pit_uniformity.png")), diags$pks, dpi=600, w=7, h=5, units="in")
  
  # Console summary + plug-in endpoint at u0
  y_star_hat <- u0 - fit0$scale/fit0$shape
  summary_row <- tibble(
    trait = tr_key,
    u0_quantile = q_anchor,
    u0 = u0,
    xi_hat = shape_hat0,
    scale_hat = fit0$scale,
    y_star_hat = y_star_hat,
    wald_p_xi_lt_0 = p_wald,
    n_exceed = sum(y > u0),
    n_total_used = length(y)
  )
  print(summary_row)
  
  list(summary = summary_row)
}

# ---------------------------
# Run all four traits
# ---------------------------
results <- map(trait_names, run_trait)
summary_all <- bind_rows(map(results, "summary"))
write.csv(summary_all, file.path(FIG_DIR, "gpd_anchor_fits_summary.csv"), row.names = FALSE)

message("Done. Outputs written to: ", normalizePath(FIG_DIR))

# =========================
# Pairwise tail regimes (dots if both observed; lines if one missing & observed one exceeds)
# Assumes df_log, trait_names, u0_by_trait, FIG_DIR, theme_science_polished exist
# =========================

pair_levels <- list(
  c("LF","LH"), c("LF","CF"), c("LF","CH"),
  c("LH","CF"), c("LH","CH"), c("CF","CH")
)

library(dplyr)
library(ggplot2)
library(glue)
library(tidyr)

make_pair_plot <- function(t1, t2) {
  # ---- file target defined first & FIG_DIR ensured ----
  if (!dir.exists(FIG_DIR)) dir.create(FIG_DIR, recursive = TRUE)
  fn <- file.path(FIG_DIR, glue::glue("pair_scatter_{t1}_{t2}_final.png"))
  
  x1 <- df_log[[paste0("log_", t1)]]
  x2 <- df_log[[paste0("log_", t2)]]
  u1 <- u0_by_trait[t1]
  u2 <- u0_by_trait[t2]
  
  O1 <- is.finite(x1);  O2 <- is.finite(x2)
  E1 <- O1 & (x1 > u1); E2 <- O2 & (x2 > u2)
  M1 <- !O1;            M2 <- !O2
  
  # ranges from observed only (fallback if all NA)
  xr <- range(x1[O1], na.rm = TRUE); if (!all(is.finite(xr))) xr <- c(-1, 1)
  yr <- range(x2[O2], na.rm = TRUE); if (!all(is.finite(yr))) yr <- c(-1, 1)
  pad <- function(r, p = 0.03) { w <- diff(r); c(r[1] - p*w, r[2] + p*w) }
  xr <- pad(xr); yr <- pad(yr)
  
  # dots: both observed only
  both <- O1 & O2
  dd_dot <- tibble::tibble(
    specimen = df$specimen[both],
    x = x1[both],
    y = x2[both],
    kind = dplyr::case_when(
      (x1[both] > u1) & (x2[both] > u2) ~ "(E1,E2) •",
      (x1[both] > u1) & (x2[both] <= u2) ~ "(E1,¬E2) •",
      (x1[both] <= u1) & (x2[both] > u2) ~ "(¬E1,E2) •",
      TRUE ~ "none •"
    )
  )
  
  # lines: exactly one observed, and that one exceeds
  dd_v <- tibble::tibble(
    specimen = df$specimen[E1 & M2],
    x = x1[E1 & M2],
    y0 = yr[1], y1 = yr[2],
    kind = "(E1,M2) |"
  )
  dd_h <- tibble::tibble(
    specimen = df$specimen[M1 & E2],
    y = x2[M1 & E2],
    x0 = xr[1], x1 = xr[2],
    kind = "(M1,E2) —"
  )
  
  # counts for subtitle
  cnt_dot <- dd_dot %>% dplyr::count(kind) %>%
    tidyr::complete(kind, fill = list(n = 0))
  n_e1e2  <- cnt_dot$n[cnt_dot$kind == "(E1,E2) •"];  if (length(n_e1e2)==0) n_e1e2 <- 0
  n_e1ne2 <- cnt_dot$n[cnt_dot$kind == "(E1,¬E2) •"]; if (length(n_e1ne2)==0) n_e1ne2 <- 0
  n_ne1e2 <- cnt_dot$n[cnt_dot$kind == "(¬E1,E2) •"]; if (length(n_ne1e2)==0) n_ne1e2 <- 0
  subt <- glue::glue(
    "dots  (E1,E2): {n_e1e2}, (E1,¬E2): {n_e1ne2}, (¬E1,E2): {n_ne1e2};  ",
    "lines (E1,M2): {nrow(dd_v)}, (M1,E2): {nrow(dd_h)}"
  )
  
  col_map <- c(
    "(E1,E2) •" = "#DC2626",  # red dot
    "(E1,¬E2) •" = "#2563EB", # blue dot
    "(¬E1,E2) •" = "#059669", # green dot
    "none •"     = "grey70",
    "(E1,M2) |"  = "#2563EB", # blue line
    "(M1,E2) —"  = "#059669"  # green line
  )
  
  p <- ggplot2::ggplot() +
    ggplot2::geom_point(data = dd_dot, ggplot2::aes(x = x, y = y, color = kind),
                        size = 2.6, alpha = 0.95, show.legend = TRUE) +
    ggplot2::geom_segment(data = dd_v,
                          ggplot2::aes(x = x, xend = x, y = y0, yend = y1, color = kind),
                          linewidth = 0.9, alpha = 0.95) +
    ggplot2::geom_segment(data = dd_h,
                          ggplot2::aes(x = x0, xend = x1, y = y, yend = y, color = kind),
                          linewidth = 0.9, alpha = 0.95) +
    ggplot2::geom_vline(xintercept = u1, linetype = "dashed", color = "red") +
    ggplot2::geom_hline(yintercept = u2, linetype = "dashed", color = "red") +
    ggplot2::scale_color_manual(
      values = col_map,
      breaks = c("(E1,E2) •","(E1,¬E2) •","(¬E1,E2) •","none •","(E1,M2) |","(M1,E2) —"),
      labels = c(
        "(E1,E2)  red dot",
        "(E1,¬E2) blue dot",
        "(¬E1,E2) green dot",
        "none     grey dot",
        "(E1,M2)  blue line",
        "(M1,E2)  green line"
      ),
      name = "Regime"
    ) +
    ggplot2::coord_cartesian(xlim = xr, ylim = yr, expand = FALSE) +
    ggplot2::labs(
      title = glue::glue("{t1} vs {t2} (log scale)"),
      subtitle = subt,
      x = glue::glue("{t1} (log)"),
      y = glue::glue("{t2} (log)")
    ) +
    theme_science_polished +
    ggplot2::theme(legend.position = "right")
  
  tryCatch(
    {
      ggplot2::ggsave(fn, p, dpi = 600, w = 6.8, h = 5.6, units = "in")
      message("Saved: ", normalizePath(fn))
    },
    error = function(e) {
      message("ggsave failed for ", fn, " — ", conditionMessage(e))
    }
  )
  invisible(p)
}

invisible(lapply(pair_levels, function(pr) make_pair_plot(pr[1], pr[2])))

# ============================================================
# Kiriliouk-style POT summary: marginals (xi, sigma, y*) + χ dependence
# ============================================================

suppressPackageStartupMessages({
  library(dplyr); library(purrr); library(tidyr)
  library(glue);  library(ggplot2); library(extRemes)
  library(readr); library(forcats); library(scales)
})

# -- 0) Ensure log data + trait thresholds exist (matches your earlier convention)
if (!exists("df_log")) {
  df_log <- df %>%
    mutate(across(all_of(trait_names), ~ suppressWarnings(as.numeric(.x)))) %>%
    mutate(across(all_of(trait_names),
                  ~ ifelse(is.finite(.) & . > 0, log(.), NA_real_),
                  .names = "log_{col}"))
}
if (!exists("u0_by_trait")) {
  u0_by_trait <- setNames(numeric(length(trait_names)), trait_names)
  for (tr in trait_names) {
    v <- df_log[[paste0("log_", tr)]]
    q_anchor <- if (!is.null(thresh_q_opt[tr])) as.numeric(thresh_q_opt[tr]) else 0.78
    u0_by_trait[tr] <- as.numeric(quantile(v, q_anchor, na.rm = TRUE))
  }
}

# -- 1) Marginal GPD fits at trait-specific u_j on log-scale (Sec. 2; threshold-stability)
fit_one_gpd <- function(tr){
  y <- df_log[[paste0("log_", tr)]]
  u <- u0_by_trait[tr]
  yy <- y[is.finite(y)]
  fit <- fevd(yy, threshold = u, type = "GP")  # MLE
  p  <- fit$results$par
  xi <- unname(p["shape"])
  sc <- unname(p["scale"])
  ystar <- if (xi < 0) u - sc/xi else NA_real_        # finite right-endpoint on log scale
  tibble(trait = tr,
         u_log = u,
         xi_hat = xi,
         sigma_hat = sc,
         y_star_hat_log = ystar,
         n_exceed = sum(yy > u),
         n_total  = length(yy))
}

marginals_tbl <- map_dfr(trait_names, fit_one_gpd)

# Save + quick print
write_csv(marginals_tbl, file.path(FIG_DIR, "kiriliouk_pot_marginals_log.csv"))
print(marginals_tbl)

# -- 2) Dependence via χ (Kiriliouk Sec. 4.3): empirical χ_{J}(q) plateaus above a joint threshold
# χ_{1:d}(q) = P(F1(Y1)>q, ..., Fd(Yd)>q) / (1-q).
emp_chi <- function(Ymat, q){
  # Ymat: columns are variables on original scale of interest; we use log traits here
  U <- apply(Ymat, 2, function(v){
    v_fin <- v[is.finite(v)]
    # empirical CDF with ranks; map back keeping NA where missing
    r <- rank(v, ties.method = "average", na.last = "keep")
    p <- r / sum(is.finite(v))
    p
  })
  # Indicator that all components exceed their marginal q-quantile
  ind <- rep(TRUE, nrow(U))
  for (j in seq_len(ncol(U))) ind <- ind & (U[, j] > q)
  num <- sum(ind, na.rm = TRUE)
  den <- sum(apply(U, 1, function(z) all(!is.na(z)))) * (1 - q)
  if (den == 0) return(NA_real_)
  num / den
}

# Helper to pick a *joint* dependence threshold vector from a scalar q* (Sec. 4.3)
# We’ll scan q-grid, show plateau, and report χ at the first stable high-q.
log_mat <- as.matrix(df_log %>% transmute(across(all_of(paste0("log_", trait_names)))))
colnames(log_mat) <- trait_names

q_grid <- seq(0.80, 0.99, by = 0.01)

# Full 4D chi(q)
chi_full <- tibble(q = q_grid,
                   chi = map_dbl(q_grid, ~ emp_chi(log_mat, .x)))

# Pairwise chi_ij(q)
pairs <- combn(trait_names, 2, simplify = FALSE)
chi_pairs <- map_dfr(pairs, function(pr){
  sub <- log_mat[, pr, drop = FALSE]
  tibble(pair = glue("{pr[1]}-{pr[2]}"),
         q = q_grid,
         chi = map_dbl(q_grid, ~ emp_chi(sub, .x)))
})

# Pick reporting level q* as the smallest q >= 0.90 where χ(q) is within 5% of its value at 0.99
pick_qstar <- function(dfq){
  chi99 <- dfq$chi[dfq$q == max(dfq$q, na.rm=TRUE)]
  cand <- dfq %>% filter(q >= 0.90, is.finite(chi)) %>% arrange(q)
  if (length(chi99) == 0 || nrow(cand) == 0) return(0.95)
  hit <- cand %>% filter(abs(chi - chi99) <= 0.05*abs(chi99))
  if (nrow(hit)) hit$q[1] else 0.95
}

q_star_full  <- pick_qstar(chi_full)
q_star_pairs <- chi_pairs %>% group_by(pair) %>% summarize(q_star = pick_qstar(cur_data_all()))

# Report χ at q*
chi_full_hat  <- chi_full %>% filter(q == q_star_full)  %>% transmute(scope = "LF-CF-CH-LH", q_star = q, chi_hat = chi)
chi_pairs_hat <- chi_pairs %>%
  inner_join(q_star_pairs, by = "pair") %>%
  filter(q == q_star) %>%
  transmute(scope = pair, q_star, chi_hat = chi)

dep_tbl <- bind_rows(chi_full_hat, chi_pairs_hat)
write_csv(dep_tbl, file.path(FIG_DIR, "kiriliouk_chi_dependence.csv"))

# Quick plots (optional): χ(q) curves with selected q*
p_full <- ggplot(chi_full, aes(q, chi)) +
  geom_line() + geom_point() +
  geom_vline(xintercept = q_star_full, linetype = "dashed", color = "red") +
  labs(title = expression(paste("Empirical ", chi["LF,CF,CH,LH"], "(q) — plateau & q*")),
       y = expression(chi), x = "q") +
  theme_science_polished
p_full
ggsave(file.path(FIG_DIR, "chi_full_plateau.png"), p_full, dpi=600, w=6.5, h=4.5, units="in")

p_pairs <- chi_pairs %>%
  inner_join(q_star_pairs, by = "pair") %>%
  ggplot(aes(q, chi)) +
  geom_line() + geom_point() +
  geom_vline(aes(xintercept = q_star), linetype = "dashed", color = "red") +
  facet_wrap(~pair, ncol = 3, scales = "free_y") +
  labs(title = expression(paste("Empirical pairwise ", chi[i*j], "(q) — plateaus & q*")),
       y = expression(chi), x = "q") +
  theme_science_polished
ggsave(file.path(FIG_DIR, "chi_pairs_plateau.png"), p_pairs, dpi=600, w=8, h=5, units="in")

# -- 3) Final compact table (marginals + dependence)
final_summary <- list(
  marginals_log = marginals_tbl %>%
    mutate(y_star_note = ifelse(is.na(y_star_hat_log), "infinite (xi>=0)", "finite")) %>%
    select(trait, u_log, xi_hat, sigma_hat, y_star_hat_log, y_star_note, n_exceed, n_total),
  dependence_chi = dep_tbl
)

print(final_summary$marginals_log)
print(final_summary$dependence_chi)

