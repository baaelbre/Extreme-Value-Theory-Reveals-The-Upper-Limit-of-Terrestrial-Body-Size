# =========================
# Multivariate EVT (LF, CF, CH) — Sauropod measurements
# =========================

suppressPackageStartupMessages({
  library(readxl); library(dplyr); library(tidyr); library(purrr)
  library(ggplot2); library(extRemes); library(scales); library(glue)
  library(forcats); library(grid); library(readr)
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
u_scan_lo_q   <- 0.60
u_scan_hi_q   <- 0.95
u_scan_n      <- 50
min_ex_mrl    <- 5
min_ex_fit    <- 5
topN_show     <- 10

# Trait-specific anchor quantiles (log scale)
thresh_q_opt <- c(LF = 0.79, CF = 0.79, CH = 0.70)
trait_names  <- c("LF", "CF", "CH")

# ---------------------------
# Ingest (exact column names)
# ---------------------------
DATA_XLSX <- "Data/sauropod_measurements_demic.xlsx"
df_raw <- read_excel(DATA_XLSX)

stopifnot(all(c("genus and species", "femur L", "femur circ (mm)", "humerus circ (mm)")
              %in% names(df_raw)))

df <- df_raw %>%
  transmute(
    specimen = `genus and species`,
    LF = suppressWarnings(as.numeric(`femur L`)),
    CF = suppressWarnings(as.numeric(`femur circ (mm)`)),
    CH = suppressWarnings(as.numeric(`humerus circ (mm)`))
  )

# ---------------------------
# Completeness (LF, CF, CH)
# ---------------------------
compl_tbl <- trait_names %>%
  set_names() %>%
  map_df(function(tr) {
    v <- df[[tr]]
    tibble(trait = tr,
           n_total = length(v),
           n_obs   = sum(!is.na(v)),
           completeness = mean(!is.na(v)))
  })

write_csv(compl_tbl, file.path(FIG_DIR, "completeness_LF_CF_CH.csv"))

p_compl <- compl_tbl %>%
  mutate(trait = fct_inorder(trait)) %>%
  ggplot(aes(trait, completeness)) +
  geom_col(fill="#3B82F6") +
  geom_text(aes(label=scales::percent(completeness, accuracy=0.1)),
            vjust=-0.2, size=4) +
  scale_y_continuous(labels=percent_format(accuracy=1), limits=c(0,1.10)) +
  labs(title="Completeness (LF, CF, CH)", x=NULL, y="Completeness") +
  theme_science_polished
p_compl
ggsave(file.path(FIG_DIR, "completeness_LF_CF_CH.png"), p_compl,
       dpi=600, w=6.0, h=4.2, units="in")

# ---------------------------
# Largest specimens per trait (Top N)
# ---------------------------
largest_by_trait <- function(tr, topN=10) {
  df %>%
    filter(!is.na(.data[[tr]])) %>%
    arrange(desc(.data[[tr]])) %>%
    mutate(rank = row_number()) %>%
    slice_head(n = topN) %>%
    transmute(trait = tr, rank, value_raw = .data[[tr]], specimen)
}
largest_list <- map_df(trait_names, largest_by_trait, topN = topN_show)
write_csv(largest_list, file.path(FIG_DIR, "largest_LF_CF_CH_top10.csv"))

# ---------------------------
# Log transform + thresholds u0 (log scale)
# ---------------------------
df_log <- df %>%
  mutate(across(all_of(trait_names), ~ suppressWarnings(as.numeric(.x)))) %>%
  mutate(across(all_of(trait_names),
                ~ ifelse(is.finite(.) & . > 0, log(.), NA_real_),
                .names = "log_{col}"))

u0_by_trait <- setNames(numeric(length(trait_names)), trait_names)
for (tr in trait_names) {
  v <- df_log[[paste0("log_", tr)]]
  u0_by_trait[tr] <- as.numeric(quantile(v, as.numeric(thresh_q_opt[tr]), na.rm = TRUE))
}
print(u0_by_trait)

# ---------------------------
# Univariate POT helpers (log scale)
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
  exc <- y[y>u] - u; n <- length(exc); probs <- ppoints(n)
  theo_q <- if (abs(shape_hat) > 1e-10) u + scale_hat/shape_hat * (probs^(-shape_hat) - 1)
  else u - scale_hat * log(probs)
  dfqq <- data.frame(Theoretical = rev(theo_q), Empirical = sort(y[y>u]))
  pqq <- ggplot(dfqq, aes(Theoretical, Empirical)) +
    geom_point(color="steelblue") +
    geom_abline(slope=1, intercept=0, linetype="dashed", color="red") +
    labs(title = glue("Q–Q: {label}")) + theme_science_polished
  F_theo <- if (abs(shape_hat)>1e-10) 1 - (1 + shape_hat*exc/scale_hat)^(-1/shape_hat) else 1 - exp(-exc/scale_hat)
  dfpp <- data.frame(Theoretical=sort(F_theo), Empirical=(1:n)/n)
  ppp <- ggplot(dfpp, aes(Theoretical, Empirical)) +
    geom_point(color="darkgreen") +
    geom_abline(slope=1, intercept=0, linetype="dashed", color="red") +
    labs(title = glue("P–P: {label}")) + theme_science_polished
  F_hat <- pevd(exc, scale=scale_hat, shape=shape_hat, type="GP")
  pks <- ggplot(data.frame(F_hat=F_hat), aes(F_hat)) +
    geom_histogram(aes(y=..density..), bins=20, fill="skyblue", color="black", alpha=0.7) +
    geom_hline(yintercept=1, color="red", linetype="dashed") +
    labs(title=glue("Uniformity (PIT): {label}")) + theme_science_polished
  list(pqq=pqq, ppp=ppp, pks=pks)
}

# ---------------------------
# Run POT per trait (LF, CF, CH)
# ---------------------------
run_trait <- function(tr_key) {
  v_raw <- df[[tr_key]] |> suppressWarnings(as.numeric())
  y <- log(v_raw[is.finite(v_raw) & v_raw > 0]); ylab <- glue("log {tr_key}")
  stopifnot(length(y) >= 40)
  q_anchor <- as.numeric(thresh_q_opt[tr_key])
  u_seq <- make_thresholds(y); u0 <- as.numeric(quantile(y, q_anchor, na.rm=TRUE))
  
  # MRL
  mrl_vals <- mrl_data(y, u_seq, min_ex=min_ex_mrl)
  p_mrl <- ggplot(data.frame(u=u_seq, mrl=mrl_vals), aes(u, mrl)) +
    geom_line() + geom_point() +
    geom_vline(xintercept = u0, linetype = "dashed", color = "red") +
    labs(title = glue("MRL plot ({ylab})"), x = glue("Threshold ({ylab})"), y = "Mean excess") +
    theme_science_polished
  ggsave(file.path(FIG_DIR, glue("{tr_key}_log_mrl.png")), p_mrl, dpi=600, w=6.2, h=4.2, units="in")
  
  # Anchor fit at u0
  fit0 <- fit_gpd_at_u(y, u0); xi0 <- fit0$shape; sc0 <- fit0$scale
  z <- xi0 / sqrt(fit0$cov[2,2]); p_wald <- pnorm(z)  # one-sided H1: xi < 0
  
  # Stability scans
  zcrit <- qnorm(1 - (1 - ci_level)/2)
  df_scan <- tibble(u=u_seq, shape=NA_real_, shape_lo=NA_real_, shape_hi=NA_real_,
                    adj=NA_real_, adj_lo=NA_real_, adj_hi=NA_real_)
  for (i in seq_along(u_seq)) {
    u <- u_seq[i]; ex <- y[y>u]-u; if (length(ex) < min_ex_fit) next
    out <- try(fit_gpd_at_u(y, u), silent=TRUE); if (inherits(out,"try-error")) next
    se_s <- sqrt(out$cov[2,2])
    df_scan$shape[i]    <- out$shape
    df_scan$shape_lo[i] <- out$shape - zcrit*se_s
    df_scan$shape_hi[i] <- out$shape + zcrit*se_s
    df_scan$adj[i]      <- adj_scale_fun(out$scale, out$shape, u, u0)
    var_adj <- out$cov[1,1] + (u-u0)^2*out$cov[2,2] - 2*(u-u0)*out$cov[1,2]
    se_adj  <- sqrt(max(var_adj,0))
    df_scan$adj_lo[i] <- df_scan$adj[i] - zcrit*se_adj
    df_scan$adj_hi[i] <- df_scan$adj[i] + zcrit*se_adj
  }
  
  p_shape <- ggplot(df_scan, aes(u, shape)) +
    geom_point(color='blue') +
    geom_errorbar(aes(ymin=shape_lo, ymax=shape_hi), width=0.03, color="blue") +
    geom_vline(xintercept=u0, color="red", linetype="dashed") +
    geom_hline(yintercept=xi0, color="red", linetype="dashed") +
    labs(x=glue("Threshold ({ylab})"), y="Shape (xi)",
         title = glue("Shape stability ({ylab})")) + theme_science_polished
  p_adj <- ggplot(df_scan, aes(u, adj)) +
    geom_point(color='blue') +
    geom_errorbar(aes(ymin=adj_lo, ymax=adj_hi), width=0.03, color="blue") +
    geom_vline(xintercept=u0, color="red", linetype="dashed") +
    geom_hline(yintercept=sc0, color="red", linetype="dashed") +
    labs(x=glue("Threshold ({ylab})"), y="Adjusted scale",
         title = glue("Adjusted-scale stability ({ylab})")) + theme_science_polished
  p_shape
  p_adj
  ggsave(file.path(FIG_DIR, glue("{tr_key}_log_shape_stability.png")), p_shape, dpi=600, w=6.2, h=4.2, units="in")
  ggsave(file.path(FIG_DIR, glue("{tr_key}_log_adj_scale_stability.png")), p_adj, dpi=600, w=6.2, h=4.2, units="in")
  
  # Diagnostics + summary
  di <- diagnostic_plots(y, u0, sc0, xi0, ylab)
  ggsave(file.path(FIG_DIR, glue("{tr_key}_log_qq.png")),  di$pqq, dpi=600, w=6.2, h=4.2, units="in")
  ggsave(file.path(FIG_DIR, glue("{tr_key}_log_pp.png")),  di$ppp, dpi=600, w=6.2, h=4.2, units="in")
  ggsave(file.path(FIG_DIR, glue("{tr_key}_log_pit_uniformity.png")), di$pks, dpi=600, w=6.2, h=4.2, units="in")
  
  y_star_hat <- if (xi0 < 0) u0 - sc0/xi0 else NA_real_
  tibble(trait = tr_key, u0_quantile = q_anchor, u0 = u0,
         xi_hat = xi0, sigma_hat = sc0, y_star_hat_log = y_star_hat,
         wald_p_xi_lt_0 = p_wald, n_exceed = sum(y > u0), n_total_used = length(y))
}

marginals_tbl <- map_dfr(trait_names, run_trait)
write_csv(marginals_tbl, file.path(FIG_DIR, "gpd_anchor_fits_summary_LF_CF_CH.csv"))
print(marginals_tbl)

# =========================
# Pairwise tail regimes — dots (both observed), lines (one missing & observed exceeds)
# =========================
pair_levels <- list(c("LF","CF"), c("LF","CH"), c("CF","CH"))

make_pair_plot <- function(t1, t2) {
  if (!dir.exists(FIG_DIR)) dir.create(FIG_DIR, recursive = TRUE)
  fn <- file.path(FIG_DIR, glue("pair_scatter_{t1}_{t2}_final.png"))
  
  x1 <- df_log[[paste0("log_", t1)]]
  x2 <- df_log[[paste0("log_", t2)]]
  u1 <- u0_by_trait[t1]; u2 <- u0_by_trait[t2]
  
  O1 <- is.finite(x1);  O2 <- is.finite(x2)
  E1 <- O1 & (x1 > u1); E2 <- O2 & (x2 > u2)
  M1 <- !O1;            M2 <- !O2
  
  xr <- range(x1[O1], na.rm=TRUE); yr <- range(x2[O2], na.rm=TRUE)
  pad <- function(r, p=0.03){ w <- diff(r); c(r[1]-p*w, r[2]+p*w) }
  xr <- pad(xr); yr <- pad(yr)
  
  both <- O1 & O2
  dd_dot <- tibble(
    specimen = df$specimen[both],
    x = x1[both], y = x2[both],
    kind = case_when(
      (x1[both] > u1) & (x2[both] > u2)  ~ "(E1,E2) •",
      (x1[both] > u1) & (x2[both] <= u2) ~ "(E1,¬E2) •",
      (x1[both] <= u1) & (x2[both] > u2) ~ "(¬E1,E2) •",
      TRUE ~ "none •"
    )
  )
  
  dd_v <- tibble( specimen = df$specimen[E1 & M2],
                  x = x1[E1 & M2], y0 = yr[1], y1 = yr[2], kind = "(E1,M2) |" )
  dd_h <- tibble( specimen = df$specimen[M1 & E2],
                  y = x2[M1 & E2], x0 = xr[1], x1 = xr[2], kind = "(M1,E2) —" )
  
  col_map <- c("(E1,E2) •"="#DC2626", "(E1,¬E2) •"="#2563EB",
               "(¬E1,E2) •"="#059669", "none •"="grey70",
               "(E1,M2) |"="#2563EB", "(M1,E2) —"="#059669")
  
  p <- ggplot() +
    geom_point(data=dd_dot, aes(x=x, y=y, color=kind), size=2.6, alpha=0.95) +
    geom_segment(data=dd_v, aes(x=x, xend=x, y=y0, yend=y1, color=kind),
                 linewidth=0.9, alpha=0.95) +
    geom_segment(data=dd_h, aes(x=x0, xend=x1, y=y, yend=y, color=kind),
                 linewidth=0.9, alpha=0.95) +
    geom_vline(xintercept=u1, linetype="dashed", color="red") +
    geom_hline(yintercept=u2, linetype="dashed", color="red") +
    scale_color_manual(values=col_map, name="Regime") +
    coord_cartesian(xlim=xr, ylim=yr, expand=FALSE) +
    labs(title=glue("{t1} vs {t2} (log scale)"),
         x=glue("{t1} (log)"), y=glue("{t2} (log)")) +
    theme_science_polished
  p
  
  ggsave(fn, p, dpi=600, w=6.8, h=5.6, units="in")
  message("Saved: ", normalizePath(fn))
  invisible(p)
}
invisible(lapply(pair_levels, function(pr) make_pair_plot(pr[1], pr[2])))

# ============================================================
# Kiriliouk-style POT dependence: empirical χ(q) (3D + pairs)
# ============================================================
emp_chi <- function(Ymat, q){
  U <- apply(Ymat, 2, function(v){
    r <- rank(v, na.last="keep", ties.method="average")
    r / sum(is.finite(v))
  })
  ok <- apply(U, 1, function(z) all(!is.na(z)))
  if (!any(ok)) return(NA_real_)
  ind <- rep(TRUE, nrow(U))
  for (j in seq_len(ncol(U))) ind <- ind & (U[,j] > q)
  sum(ind & ok) / (sum(ok) * (1 - q))
}

log_mat <- as.matrix(df_log %>% transmute(LF = log_LF, CF = log_CF, CH = log_CH))
q_grid  <- seq(0.80, 0.99, by = 0.01)

chi_full <- tibble(q = q_grid, chi = map_dbl(q_grid, ~ emp_chi(log_mat, .x)))

pairs <- list(c("LF","CF"), c("LF","CH"), c("CF","CH"))
chi_pairs <- map_dfr(pairs, function(pr){
  sub <- log_mat[, pr, drop=FALSE]
  tibble(pair = glue("{pr[1]}-{pr[2]}"),
         q = q_grid,
         chi = map_dbl(q_grid, ~ emp_chi(sub, .x)))
})

pick_qstar <- function(dfq){
  chi99 <- dfq$chi[dfq$q == max(dfq$q, na.rm=TRUE)]
  cand  <- dfq %>% filter(q >= 0.90, is.finite(chi)) %>% arrange(q)
  if (!length(chi99) || !nrow(cand)) return(0.95)
  hit <- cand %>% filter(abs(chi - chi99) <= 0.05*abs(chi99))
  if (nrow(hit)) hit$q[1] else 0.95
}

q_star_full  <- pick_qstar(chi_full)
q_star_pairs <- chi_pairs %>% group_by(pair) %>% summarize(q_star = pick_qstar(cur_data_all()))

chi_full_hat  <- chi_full %>% filter(q == q_star_full) %>%
  transmute(scope = "LF-CF-CH", q_star = q, chi_hat = chi)
chi_pairs_hat <- chi_pairs %>%
  inner_join(q_star_pairs, by="pair") %>%
  filter(q == q_star) %>%
  transmute(scope = pair, q_star, chi_hat = chi)

dep_tbl <- bind_rows(chi_full_hat, chi_pairs_hat)
write_csv(dep_tbl, file.path(FIG_DIR, "kiriliouk_chi_dependence_LF_CF_CH.csv"))

p_full <- ggplot(chi_full, aes(q, chi)) +
  geom_line() + geom_point() +
  geom_vline(xintercept = q_star_full, linetype = "dashed", color = "red") +
  labs(title = expression(paste("Empirical ", chi["LF,CF,CH"], "(q) — plateau & q*")),
       y = expression(chi), x = "q") +
  theme_science_polished
ggsave(file.path(FIG_DIR, "chi_full_LF_CF_CH_plateau.png"),
       p_full, dpi=600, w=6.0, h=4.2, units="in")

p_pairs <- chi_pairs %>%
  inner_join(q_star_pairs, by="pair") %>%
  ggplot(aes(q, chi)) +
  geom_line() + geom_point() +
  geom_vline(aes(xintercept = q_star), linetype="dashed", color="red") +
  facet_wrap(~pair, ncol = 3, scales="free_y") +
  labs(title = expression(paste("Empirical pairwise ", chi[i*j], "(q) — plateaus & q*")),
       y = expression(chi), x = "q") +
  theme_science_polished
p_pairs
ggsave(file.path(FIG_DIR, "chi_pairs_LF_CF_CH_plateau.png"),
       p_pairs, dpi=600, w=8.0, h=4.6, units="in")

# ---------------------------
# Final compact outputs
# ---------------------------
write_csv(marginals_tbl %>%
            select(trait, u0, xi_hat, sigma_hat, y_star_hat_log,
                   wald_p_xi_lt_0, n_exceed, n_total_used),
          file.path(FIG_DIR, "marginals_POT_LF_CF_CH_log.csv"))

message("Done. Outputs in: ", normalizePath(FIG_DIR))

