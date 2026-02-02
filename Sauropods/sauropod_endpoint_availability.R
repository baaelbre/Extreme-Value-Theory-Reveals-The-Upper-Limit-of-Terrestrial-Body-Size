############################################################
#  Title: Multi-dimension POT analysis (ln-scale)
#  Author: Bastiaan Aelbrecht
#  Date  : 2025‑08‑08
############################################################

############################
# 1.  SETUP & UTILITIES    #
############################
library(ggplot2)
library(extRemes)
library(purrr)
library(dplyr)
library(tibble)
library(stringr)
library(readxl)
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

FIG_DIR <- "Figures_sauropods"
if (!dir.exists(FIG_DIR)) dir.create(FIG_DIR, recursive = TRUE)

# Safe filename
safe_name <- function(x) {
  x %>%
    tolower() %>%
    str_replace_all("[^a-z0-9]+", "_") %>%
    str_replace("^_|_$", "")
}

############################
# 2.  DATA LOADING         #
############################
df_all <- read.csv("Data/DEmic23_updated_Supplemental_Data_withPubYear_withAvailability.csv")

# Candidate dimensions (will silently skip those not present)
dim_cols <- c(
  "hum+fem circ (mm)",
  "humerus circ (mm)", "femur circ (mm)",
  "tibia circumference (mm)",
  "humerus AP", "humerus ML", "humerus L",
  "femur AP", "femur ML", "femur L"
)

# Pretty labels for x-axes
pretty_label <- function(col) {
  switch(col,
         "hum+fem circ (mm)"     = expression(ln(C[F+H])~"[ln(mm)]"),
         "humerus circ (mm)"     = expression(ln(C[H])~"[ln(mm)]"),
         "femur circ (mm)"       = expression(ln(C[F])~"[ln(mm)]"),
         "tibia circumference (mm)" = expression(ln(C[Tibia])~"[ln(mm)]"),
         "humerus AP"            = expression(ln(AP[H])~"[ln(mm)]"),
         "humerus ML"            = expression(ln(ML[H])~"[ln(mm)]"),
         "humerus L"             = expression(ln(L[H])~"[ln(mm)]"),
         "femur AP"              = expression(ln(AP[F])~"[ln(mm)]"),
         "femur ML"              = expression(ln(ML[F])~"[ln(mm)]"),
         "femur L"               = expression(ln(L[F])~"[ln(mm)]"),
         # fallback
         bquote(ln(.(col))~"[ln(mm)]")
  )
}

############################
# 3.  GENERIC PLOTTING     #
############################
plot_hist_density <- function(x_log, xlab, title_suffix) {
  df <- tibble(x = x_log)
  ggplot(df, aes(x = x)) +
    geom_histogram(aes(y = ..density..), binwidth = 0.05,
                   fill = "skyblue", color = "black", alpha = 0.7) +
    geom_density(size = 1.2, adjust = 1.2) +
    scale_x_continuous(name = xlab,
                       sec.axis = sec_axis(~exp(.), name = "Original scale [mm]")) +
    labs(y = "Density", title = paste("Histogram & Density —", title_suffix)) +
    theme_science
}

plot_qq <- function(x_log, title_suffix) {
  df <- tibble(x = x_log)
  ggplot(df, aes(sample = x)) +
    stat_qq(size = 2) +
    stat_qq_line(linetype = "dashed", size = 1) +
    labs(title = paste("Q–Q (Normal) —", title_suffix),
         x = "Theoretical Quantiles", y = "Sample Quantiles") +
    theme_science
}

############################
# 4.  THRESHOLD STABILITY  #
############################
find_threshold <- function(x_log, n_grid = 50, alpha = 0.05, q_start = 0.85, q_end = 0.98) {
  stopifnot(is.numeric(x_log), length(x_log) > 10)
  x_log <- sort(x_log)
  
  u0 <- as.numeric(quantile(x_log, q_start, na.rm = TRUE))
  umax <- as.numeric(quantile(x_log, q_end,   na.rm = TRUE))
  if (!is.finite(u0) || !is.finite(umax) || u0 >= umax) {
    # Fallback: use mid-to-high region
    u0 <- x_log[round(0.8 * length(x_log))]
    umax <- x_log[round(0.98 * length(x_log))]
  }
  thresholds <- unique(seq(u0, umax, length.out = n_grid))
  z_crit <- qnorm(1 - alpha/2)
  
  rows <- map(thresholds, function(u) {
    xs <- x_log[x_log > u]
    if (length(xs) < 20) return(NULL)  # avoid too-few exceedances
    fit <- fevd(x_log, threshold = u, type = "GP")
    pars <- fit$results$par
    se   <- tryCatch(summary(fit)$se,        error = function(e) c(scale = NA, shape = NA))
    cov  <- tryCatch(summary(fit)$cov.theta, error = function(e) diag(2))
    
    sigma_hat <- pars["scale"]
    xi_hat    <- pars["shape"]
    sigma_ci  <- sigma_hat + c(-1, 1) * z_crit * se["scale"]
    xi_ci     <- xi_hat    + c(-1, 1) * z_crit * se["shape"]
    
    # Adjusted scale σ*(u) = σ(u) − ξ(u)(u − u0)
    sigma_adj_hat <- sigma_hat - xi_hat * (u - u0)
    var_adj <- cov[1,1] + (u - u0)^2 * cov[2,2] - 2*(u - u0) * cov[1,2]
    sd_adj  <- if (is.finite(var_adj) && var_adj >= 0) sqrt(var_adj) else NA_real_
    adj_ci  <- sigma_adj_hat + c(-1, 1) * z_crit * sd_adj
    
    tibble(
      threshold = u,
      n_exc = length(xs),
      sigma = sigma_hat,
      sigma_l = sigma_ci[1],
      sigma_u = sigma_ci[2],
      xi = xi_hat,
      xi_l = xi_ci[1],
      xi_u = xi_ci[2],
      sigma_adj = sigma_adj_hat,
      sigma_adj_l = adj_ci[1],
      sigma_adj_u = adj_ci[2]
    )
  }) %>% compact() %>% bind_rows()
  
  rows
}

choose_threshold <- function(th_tbl) {
  # Heuristic: pick threshold where xi is "most stable" (closest to median)
  if (nrow(th_tbl) == 0) return(NA_real_)
  med_xi <- median(th_tbl$xi, na.rm = TRUE)
  idx <- which.min(abs(th_tbl$xi - med_xi))
  th_tbl$threshold[idx]
}

plot_param_vs_threshold <- function(tbl, y, y_low, y_up, ylab, u_star, title_suffix) {
  ggplot(tbl, aes(x = threshold, y = .data[[y]])) +
    geom_point() +
    geom_errorbar(aes(ymin = .data[[y_low]], ymax = .data[[y_up]]),
                  width = 0.08, colour = "blue") +
    geom_vline(xintercept = u_star, colour = "red", linetype = "dashed") +
    labs(
      x = "Threshold (ln scale)",
      y = ylab,
      title = paste("Threshold Stability —", title_suffix)
    ) +
    theme_science
}

############################
# 5.  GPD FIT & GOF         #
############################
gpd_gof <- function(x_log, u, xi, sigma) {
  excess <- x_log[x_log > u] - u
  n <- length(excess)
  if (n < 20) return(NULL)
  
  probs  <- ppoints(n)
  theo_q <- if (abs(xi) > 1e-6) u + sigma/xi * (probs^(-xi) - 1) else u - sigma * log(probs)
  qq_df  <- tibble(Theoretical = rev(theo_q), Empirical = sort(x_log[x_log > u]))
  
  pp_theo <- if (abs(xi) > 1e-6) 1 - (1 + xi * excess / sigma)^(-1/xi) else 1 - exp(-excess / sigma)
  pp_df   <- tibble(Theoretical = sort(pp_theo), Empirical = (1:n)/n)
  
  list(qq_df = qq_df, pp_df = pp_df)
}

plot_gof <- function(gof_df, title_suffix) {
  p1 <- ggplot(gof_df$qq_df, aes(Theoretical, Empirical)) +
    geom_point(color = "steelblue") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    labs(title = paste("Q–Q plot: GPD —", title_suffix),
         x = "Theoretical", y = "Empirical") +
    theme_science
  
  p2 <- ggplot(gof_df$pp_df, aes(Theoretical, Empirical)) +
    geom_point(color = "darkgreen") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    labs(title = paste("P–P plot: GPD —", title_suffix),
         x = "Theoretical CDF", y = "Empirical CDF") +
    theme_science
  
  list(qq = p1, pp = p2)
}

############################
# 6.  RUN FOR ALL COLS     #
############################
analyze_one_dimension <- function(df, col) {
  if (!col %in% names(df)) return(NULL)
  x_raw <- suppressWarnings(as.numeric(df[[col]]))
  x_raw <- x_raw[is.finite(x_raw) & !is.na(x_raw)]
  x_raw <- x_raw[x_raw > 0]  # drop nonpositive
  if (length(x_raw) < 50) return(NULL)  # too few values for stable POT
  
  x_log <- log(x_raw)
  label_expr <- pretty_label(col)
  tag <- safe_name(col)
  
  # (A) Exploratory plots
  p_hist <- plot_hist_density(x_log, label_expr, col)
  p_qqn  <- plot_qq(x_log, col)
  p_hist
  p_qqn
  
  ggsave(file.path(FIG_DIR, paste0(tag, "_hist_density.png")), p_hist, width = 6.5, height = 4.5, dpi = 300)
  ggsave(file.path(FIG_DIR, paste0(tag, "_qq_normal.png")),   p_qqn,  width = 6.5, height = 4.5, dpi = 300)
  
  # (B) Threshold stability
  th_tbl <- find_threshold(x_log, n_grid = 50, q_start = 0.85, q_end = 0.98)
  if (nrow(th_tbl) == 0) return(NULL)
  u_star <- choose_threshold(th_tbl)
  
  p_sigma_adj <- plot_param_vs_threshold(th_tbl, "sigma_adj", "sigma_adj_l", "sigma_adj_u",
                                         "Adjusted scale", u_star, col)
  p_xi        <- plot_param_vs_threshold(th_tbl, "xi", "xi_l", "xi_u",
                                         "Shape (xi)", u_star, col)
  
  ggsave(file.path(FIG_DIR, paste0(tag, "_thresh_sigma_adj.png")), p_sigma_adj, width = 6.5, height = 4.5, dpi = 300)
  ggsave(file.path(FIG_DIR, paste0(tag, "_thresh_xi.png")),        p_xi,        width = 6.5, height = 4.5, dpi = 300)
  
  # (C) GPD fit & GOF at chosen threshold
  fit <- fevd(x_log, threshold = u_star, type = "GP")
  xi_hat    <- fit$results$par["shape"]
  sigma_hat <- fit$results$par["scale"]
  
  gof_df <- gpd_gof(x_log, u_star, xi_hat, sigma_hat)
  if (!is.null(gof_df)) {
    gof_plots <- plot_gof(gof_df, col)
    ggsave(file.path(FIG_DIR, paste0(tag, "_qq_gpd.png")), gof_plots$qq, width = 6.5, height = 4.5, dpi = 300)
    ggsave(file.path(FIG_DIR, paste0(tag, "_pp_gpd.png")), gof_plots$pp, width = 6.5, height = 4.5, dpi = 300)
  }
  
  tibble(
    dimension = col,
    n_total = length(x_raw),
    u_star = u_star,
    n_exceed = sum(x_log > u_star),
    xi_hat = as.numeric(xi_hat),
    sigma_hat = as.numeric(sigma_hat)
  )
}

present_cols <- intersect(dim_cols, names(df_all))
results <- map(present_cols, ~ analyze_one_dimension(df_all, .x)) %>% compact() %>% bind_rows()

# Save a CSV summary of fitted params per dimension
write.csv(results, file.path(FIG_DIR, "pot_summary_by_dimension.csv"), row.names = FALSE)

# Print summary to console
print(results)

