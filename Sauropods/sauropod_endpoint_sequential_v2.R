############################################################
#  Title: Sequential POT analysis on log CFH (Uniform y*, improper ξ)
#  Author: Bastiaan Aelbrecht (sequential, solid lines, ≥10 ex, overlay all exceedances)
#  Date  : 2025-09-09
#  Notes : Matches the single-shot endpoint pipeline:
#          - u0 fixed (global quantile).
#          - y* ~ Uniform(L, U) with L = max exceedance + eps, U = u0 + 1.05
#          - ξ prior: improper 1{ξ<0}
#          - Integrate ξ by trapz → p(y* | data ≤ t)
#          - Mass mapping via centered quadratic; add one-sided predictive uplift
#          - Start at ≥ 10 exceedances
#          - Plot: Top-10 with two-sided PIs; Extremosaurus MAP (purple, solid)
#                  & upper 95% (orange, solid); all exceedances as light-gray points
############################################################

## ---------------------------
## 0) Libraries & setup
## ---------------------------
library(ggplot2); library(dplyr); library(readxl)
library(pracma);   library(HDInterval); library(MASS)
library(scales);   library(grid); library(tidyr)

set.seed(42)

FIG_DIR <- "Figures"
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

## ---------------------------
## 1) Settings
## ---------------------------
ci_upper        <- 0.95            # one-sided upper (Extremosaurus)
n_post_ystar    <- 30000
thresh_q_opt    <- 0.78
min_exceedances <- 10              # start analysis at ≥ 10 exceedances
xi_grid         <- seq(-1.0, -0.02, length.out = 1200)

# predictive uplift on mass (individual-at-endpoint)
add_predictive  <- TRUE
q_indiv         <- ci_upper

## ---------------------------
## 2) Data ingest
## ---------------------------
dat_path <- "Data/DEMic23_updated_Supplemental_Data_withPubYear.xlsx"
df <- read_excel(dat_path) |>
  dplyr::select(`hum+fem circ (mm)`, `publication year`) |>
  dplyr::rename(sum_circ_mm = `hum+fem circ (mm)`,
                year        = `publication year`) |>
  tidyr::drop_na()
df$log_circ <- log(as.double(df$sum_circ_mm))
stopifnot(nrow(df) >= 80)

y_full <- df$log_circ
u0     <- as.numeric(quantile(y_full, thresh_q_opt, na.rm = TRUE))
years_eval <- sort(unique(df$year))

## ---------------------------
## 3) Allometry (centered quadratic)
## ---------------------------
coeff_file_cq <- "centered_quadratic_coefficients.rds"
coeffs <- readRDS(coeff_file_cq)
if (is.null(coeffs$V)) {
  coeffs$V <- matrix(c(
    coeffs$alpha_se^2, coeffs$cov_ab,      coeffs$cov_ag,
    coeffs$cov_ab,      coeffs$beta_se^2,  coeffs$cov_bg,
    coeffs$cov_ag,      coeffs$cov_bg,     coeffs$gamma_se^2
  ), nrow = 3, byrow = TRUE)
}
map_logM_cq <- function(logC) {
  a <- coeffs$alpha; b <- coeffs$beta; g <- coeffs$gamma; c0 <- coeffs$mean_log_sum_circ
  a + b*(logC - c0) + g*(logC - c0)^2
}

## ---------------------------
## 4) Priors & likelihood  (Uniform y*, improper ξ)
## ---------------------------
logprior_xi <- function(xi) if (xi < 0) 0 else -Inf  # improper ξ prior

# Weibull-domain endpoint parameterization
loglik_endpoint <- function(y_star, xi, u, y_ex) {
  if (xi >= 0 || any(y_star <= y_ex)) return(-Inf)
  n <- length(y_ex)
  n*log(y_star - u)/xi - (1/xi + 1)*sum(log(y_star - y_ex)) - n*log(-xi)
}

## ---------------------------
## 5) Sequential loop (year ≤ t) — start at ≥ 10 exceedances
## ---------------------------
results <- tibble()

for (t in years_eval) {
  print(t)
  y_sub   <- df$log_circ[df$year <= t]
  y_above <- y_sub[y_sub > u0]
  n_ex    <- length(y_above)
  if (n_ex < min_exceedances) next
  
  Ymax_t <- max(y_above)
  eps    <- max(1e-10, 1e-6 * abs(Ymax_t))
  L      <- Ymax_t + eps
  U      <- u0 + 1.05
  
  # Uniform prior on y* in (L,U)
  logprior_y <- function(y_star) if (y_star > L && y_star < U) -log(U - L) else -Inf
  
  # grid for y*
  y_star_grid <- seq(L, U, length.out = 5000)
  
  # integrate out ξ
  marg_log_post <- vapply(y_star_grid, function(ys){
    lv <- vapply(xi_grid, function(xi){
      loglik_endpoint(ys, xi, u0, y_above) + logprior_xi(xi)
    }, numeric(1))
    m <- max(lv)
    log(pracma::trapz(xi_grid, exp(lv - m))) + m + logprior_y(ys)
  }, numeric(1))
  
  post_unn <- exp(marg_log_post - max(marg_log_post))
  Z_norm   <- pracma::trapz(y_star_grid, post_unn)
  if (!is.finite(Z_norm) || Z_norm <= 0) next
  post_y   <- post_unn / Z_norm
  
  # draws for y*
  ystar_draws <- sample(y_star_grid, size = n_post_ystar, replace = TRUE, prob = post_y)
  
  # propagate to MASS with parameter uncertainty + one-sided predictive uplift
  abg <- MASS::mvrnorm(n_post_ystar, mu = c(coeffs$alpha, coeffs$beta, coeffs$gamma), Sigma = coeffs$V)
  a <- abg[,1]; b <- abg[,2]; g <- abg[,3]; x <- (ystar_draws - coeffs$mean_log_sum_circ)
  logM_draws <- a + b*x + g*x^2
  if (add_predictive) {
    stopifnot(!is.null(coeffs$resid_sd))
    logM_draws <- logM_draws + qnorm(q_indiv) * coeffs$resid_sd
  }
  mass_draws_t <- exp(logM_draws) / 1e6
  
  # MASS MAP from the predictive draws (matches single-shot)
  kd <- density(mass_draws_t, n = 4096, adjust = 1.0)
  mass_map_t <- kd$x[which.max(kd$y)]
  
  # one-sided upper bound
  mass_upper95_t <- as.numeric(quantile(mass_draws_t, probs = ci_upper))
  mass_upper90_t  <- as.numeric(quantile(mass_draws_t, probs = 0.90))
  
  results <- bind_rows(
    results,
    tibble(
      year              = t,
      n_ex              = n_ex,
      u0                = u0,
      map_mass_tons     = mass_map_t,     
      mass_upper95_tons = mass_upper95_t,
      mass_upper90_tons = mass_upper90_t
    )
  )
}

stopifnot(nrow(results) > 0)
saveRDS(results, "results_sequential.RDS")
results <- readRDS('results_sequential.RDS')

## ---------------------------
## 6) Build overlays
## ---------------------------

# Top-10 (fixed data you provided)
top10_data <- data.frame(
  species = c("Ruyangosaurus","Turiasaurus","Yunmenglong","Australotitan",
              "Patagotitan","Dreadnoughtus","Notocolossus","Chucarosaurus",
              "Brachiosaurus","Argentinosaurus"),
  year = c(2017,2006,2013,2021,2017,2014,2016,2024,1903,1993),
  sum_circ_mm = c(1639,1672,1687,1692,1694,1695,1704,1827,1870,2016),
  mass_quad   = c(43.4,45.7,46.7,47.1,47.2,47.3,48.0,57.5,61.0,74.1),
  pi_lower_quad = c(19.8,20.8,21.2,21.4,21.4,21.5,21.7,25.8,27.3,32.7),
  pi_upper_quad = c(95.0,100.0,103.0,104.0,104.0,104.0,106.0,128.0,137.0,168.0)
)
top10_data$type <- "Specimen"

# Mass estimates for ALL exceedances (y > u0) using the same quadratic mapping (no PI bars; points only)
df$mass_quad_all <- exp(map_logM_cq(df$log_circ)) / 1e6
exceed_pts <- df |> dplyr::filter(log_circ > u0) |> dplyr::select(year, mass_quad_all)

# palette for named specimens
species_colors <- c(
  "Ruyangosaurus"="#a6cee3","Turiasaurus"="#1f78b4","Yunmenglong"="#b2df8a",
  "Australotitan"="#33a02c","Patagotitan"="#fb9a99","Dreadnoughtus"="#e31a1c",
  "Notocolossus"="#fdbf6f","Chucarosaurus"="#ff7f00","Brachiosaurus"="#cab2d6",
  "Argentinosaurus"="#6a3d9a","Extremosaurus"="firebrick"
)

## ---------------------------
## 7) Plot — Top-10 (two-sided) + all exceedances (gray points) + Extremosaurus (solid lines)
## ---------------------------
## ---- add after species_colors and before plotting ----------------------

# 1) Create a data frame for the MAP line with a 'species' label
ext_map <- results |>
  dplyr::transmute(year, mass = map_mass_tons, species = "Extremosaurus")

# 2) Extend the color palette with a legend entry for Extremosaurus (MAP)
species_colors2 <- c(
  species_colors,
  "Extremosaurus (MAP)" = "firebrick"  # matches the MAP line color
)

# 3) Desired legend order (top-10 + Extremosaurus)
legend_order <- c(
  "Ruyangosaurus","Turiasaurus","Yunmenglong","Australotitan",
  "Patagotitan","Dreadnoughtus","Notocolossus","Chucarosaurus",
  "Brachiosaurus","Argentinosaurus","Extremosaurus"
)

## ---- re-build p_mass with legend --------------------------------------

p_mass <- ggplot() +
  # All exceedances (no legend)
  geom_point(
    data = exceed_pts,
    aes(x = year, y = mass_quad_all),
    color = "gray70", alpha = 0.5, size = 1.8, shape = 16, show.legend = FALSE
  ) +
  # Top-10: keep color legend, but avoid duplicate guides from errorbars/points
  geom_errorbar(
    data = top10_data,
    aes(x = year, ymin = pi_lower_quad, ymax = pi_upper_quad, color = species),
    width = 0.5, linewidth = 1, show.legend = FALSE
  ) +
  geom_line(
    data = top10_data,
    aes(x = year, y = mass_quad, color = species),
    linewidth = 1.1, alpha = 0.9, show.legend = TRUE
  ) +
  geom_point(
    data = top10_data,
    aes(x = year, y = mass_quad, color = species),
    size = 2.6, shape = 16, show.legend = FALSE
  ) +
  # Extremosaurus MAP as a colored line with its own legend entry
  geom_line(
    data = ext_map,
    aes(x = year, y = mass, color = species),
    linewidth = 1.3, show.legend = TRUE
  ) +
  labs(x = "Publication Year", y = "Estimated Mass [tons]", color = "Specimens") +
  scale_color_manual(values = species_colors2, limits = legend_order, breaks = legend_order) +
  scale_x_continuous(limits = c(1892, max(c(results$year, top10_data$year, df$year), na.rm = TRUE))) +
  ylim(0, 250) +
  theme_science_polished +
  guides(color = guide_legend(override.aes = list(linewidth = 1.3, alpha = 1)))

p_mass
ggsave(file.path(FIG_DIR, "mass_sequential_map.png"),
       p_mass, dpi = 600, w = 7, h = 5, units = "in")


p_mass_supplement <- ggplot() +
  # All exceedances as light-gray points (faint)
  geom_point(
    data = exceed_pts,
    aes(x = year, y = mass_quad_all),
    color = "gray70", alpha = 0.5, size = 1.8, shape = 16, show.legend = FALSE
  ) +
 
  # Extremosaurus: MAP (purple, solid) & upper 95% (orange, solid)
  geom_line(
    data = results,
    aes(x = year, y = map_mass_tons),
    color = "purple", linewidth = 1.3
  ) +
  geom_line(
    data = results,
    aes(x = year, y = mass_upper90_tons),
    color = "firebrick", linewidth = 1.3
  )+
  geom_line(
    data = results,
    aes(x = year, y = mass_upper95_tons),
    color = "orange", linewidth = 1.3
  )+
  
  labs(x = "Publication Year", y = "Estimated Mass [tons]") +
  scale_color_manual(values = species_colors) +
  scale_x_continuous(limits = c(1892, max(c(results$year, top10_data$year, df$year), na.rm = TRUE))) +
  ylim(0,800) +
  theme_science_polished

p_mass_supplement

ggsave(file.path(FIG_DIR, "mass_sequential_supplement.png"),
       p_mass_supplement, dpi = 600, w = 7, h = 5, units = "in")


## ---------------------------
## 7b) Same plot, log-scale Y with custom ticks
## ---------------------------
p_mass_log <- ggplot() +
  # All exceedances (gray) first so colored points sit on top
  geom_point(
    data = exceed_pts,
    aes(x = year, y = mass_quad_all),
    color = "gray70", alpha = 0.5, size = 1.8, shape = 16, show.legend = FALSE
  ) +
  geom_errorbar(
    data = top10_data,
    aes(x = year, ymin = pi_lower_quad, ymax = pi_upper_quad, color = species),
    width = 0.5, linewidth = 1, show.legend = FALSE
  ) +
  geom_line(
    data = top10_data,
    aes(x = year, y = mass_quad, color = species),
    linewidth = 1.1, alpha = 0.9, show.legend = FALSE
  ) +
  geom_point(
    data = top10_data,
    aes(x = year, y = mass_quad, color = species),
    size = 2.6, shape = 16, show.legend = FALSE
  ) +
  geom_line(
    data = results,
    aes(x = year, y = map_mass_tons),
    color = "purple", linewidth = 1.3
  ) +
  geom_line(
    data = results,
    aes(x = year, y = mass_upper95_tons),
    color = "orange", linewidth = 1.3
  ) +
  labs(x = "Publication Year", y = "Estimated Mass [tons] (log scale)") +
  scale_color_manual(values = species_colors) +
  scale_x_continuous(limits = c(1892, max(c(results$year, top10_data$year, df$year), na.rm = TRUE))) +
  scale_y_log10(
    limits = c(10, max(c(top10_data$pi_upper_quad, results$mass_upper95_tons, exceed_pts$mass_quad_all), na.rm = TRUE) * 1.1),
    breaks = c(10, 20, 50, 100, 200, 500, 1000),
    labels = scales::comma
  ) +
  theme_science_polished

p_mass_log
ggsave(file.path(FIG_DIR, "mass_sequential_uniform_improper_one_sided_log.png"),
       p_mass_log, dpi = 600, w = 7, h = 5, units = "in")
ggsave(file.path(FIG_DIR, "mass_sequential_uniform_improper_one_sided_log.pdf"),
       p_mass_log, device = cairo_pdf, w = 7, h = 5, units = "in")

