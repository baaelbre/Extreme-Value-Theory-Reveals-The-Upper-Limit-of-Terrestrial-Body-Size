suppressPackageStartupMessages({
  library(rstan)
  library(MASS)
  library(ggplot2); library(dplyr); library(readxl); library(tidyr)
  library(scales);  library(grid)
})

rstan_options(auto_write = TRUE)
options(mc.cores = max(1, parallel::detectCores() - 1))
set.seed(42); options(stringsAsFactors = FALSE)

# =============================================================================
# Paths & theme
# =============================================================================
DATA_XLSX  <- "Data/DEMic23_updated_Supplemental_Data_withPubYear.xlsx"
COEFFS_RDS <- "centered_quadratic_coefficients.rds"
# >>> Collapsed Stan model (no hyperpriors)
STAN_FILE  <- "stan/hier_gp_weibull_sauropod_3.stan"
FIG_DIR    <- "Figures"
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

# =============================================================================
# Config (matches collapsed model)
# =============================================================================
ci_upper        <- 0.95
ci_upper_alt    <- 0.90
min_exceedances <- 10
q_circ          <- 0.78        # fixed global threshold quantile

# Direct priors (collapsed): y* ~ N(MU_Y0, SD_MU_Y) (trunc), xi ~ N(mu_xi0, sd_xi0) (trunc)
MU_Y0   <- 8.0
SD_MU_Y <- 0.20
mu_xi0  <- -0.25
sd_xi0  <- 0.10

# y* bounds (Stan enforces via y_min/y_max/y_cap)
Y_MIN <- 0.0
Y_MAX <- 8.62     # legacy/user cap (safety)
Y_CAP <- 9.62     # hard cap (~1000 t mass scale)

# Stan sampler controls (tweak for speed/accuracy)
CHAINS <- 4; ITER <- 10000; WARMUP <- 2000
CONTROL <- list(adapt_delta = 0.95, max_treedepth = 12)

# draws for mass propagation per year
N_PROP <- 30000

# =============================================================================
# Helpers
# =============================================================================
kde_df <- function(x, from = quantile(x,0.001), to = quantile(x,0.999), n = 4096) {
  d <- density(x, from = from, to = to, n = n, kernel = "gaussian", adjust = 1.0)
  data.frame(t = d$x, d = d$y)
}
mode_from_kde <- function(x) { d <- kde_df(x); d$t[ which.max(d$d) ] }

# =============================================================================
# Data ingest (log CFH)
# =============================================================================
stopifnot(file.exists(DATA_XLSX), file.exists(COEFFS_RDS), file.exists(STAN_FILE))
df <- read_excel(DATA_XLSX) |>
  dplyr::select(`hum+fem circ (mm)`, `publication year`) |>
  dplyr::rename(sum_circ_mm = `hum+fem circ (mm)`,
                year        = `publication year`) |>
  tidyr::drop_na()
df$log_circ <- log(as.double(df$sum_circ_mm))
stopifnot(nrow(df) >= 80)

y_all      <- df$log_circ
u0         <- as.numeric(quantile(y_all, q_circ, na.rm = TRUE))
years_eval <- sort(unique(df$year))

# Allometry coefficients (for derived mass)
coeffs  <- readRDS(COEFFS_RDS)
stopifnot(all(c("alpha","beta","gamma","mean_log_sum_circ","resid_sd") %in% names(coeffs)))
if (is.null(coeffs$V) && all(c("alpha_se","beta_se","gamma_se","cov_ab","cov_bg","cov_ag") %in% names(coeffs))) {
  coeffs$V <- matrix(c(
    coeffs$alpha_se^2, coeffs$cov_ab,       coeffs$cov_ag,
    coeffs$cov_ab,     coeffs$beta_se^2,    coeffs$cov_bg,
    coeffs$cov_ag,     coeffs$cov_bg,       coeffs$gamma_se^2
  ), 3, 3, byrow = TRUE)
}
alpha  <- as.numeric(coeffs$alpha)
beta   <- as.numeric(coeffs$beta)
gamma  <- as.numeric(coeffs$gamma)
c0     <- as.numeric(coeffs$mean_log_sum_circ)
Vabg   <- as.matrix(coeffs$V)
sigma_e<- as.numeric(coeffs$resid_sd)

# Precompile Stan
sm <- rstan::stan_model(file = STAN_FILE)

# Helper to build Stan data for a given subset (collapsed priors only)
build_stan_data <- function(y_sub, u_fixed) {
  y_above <- y_sub[y_sub > u_fixed]
  N <- length(y_above)
  ymax <- if (N > 0) max(y_above) else -Inf
  
  ystar_lb <- max(u_fixed, ymax) + 1e-6
  y_lower  <- max(ystar_lb, Y_MIN)
  y_upper  <- min(Y_CAP, Y_MAX)
  stopifnot(is.finite(y_lower), is.finite(y_upper), y_upper > y_lower)
  
  list(
    N = N,
    y = as.vector(y_above),
    u = u_fixed,
    ymax = ymax,
    y_min = Y_MIN, y_max = Y_MAX, y_cap = Y_CAP,
    mu_y0 = MU_Y0, sd_mu_y = SD_MU_Y,
    mu_xi0 = mu_xi0, sd_xi0 = sd_xi0,
    alpha = alpha, beta = beta, gamma = gamma, c0 = c0
  )
}

# Warm-start helper: use previous posterior means (only ystar, xi remain)
make_inits <- function(prev) {
  if (is.null(prev)) return("random")
  function() list(ystar = prev$ystar, xi = prev$xi)
}

# =============================================================================
# Sequential loop over years
# =============================================================================
results <- dplyr::tibble()
prev_summ <- NULL

for (t in years_eval) {
  y_sub <- df$log_circ[df$year <= t]
  y_abv <- y_sub[y_sub > u0]
  n_ex  <- length(y_abv)
  if (n_ex < min_exceedances) next
  
  stan_data <- build_stan_data(y_sub, u_fixed = u0)
  
  fit_t <- rstan::sampling(
    sm, data = stan_data, seed = 42,
    chains = CHAINS, iter = ITER, warmup = WARMUP,
    init = make_inits(prev_summ), control = CONTROL, refresh = 0
  )
  
  post <- as.data.frame(fit_t)
  # cache posterior means for warm-start of next year
  prev_summ <- list(
    ystar = mean(post$ystar),
    xi    = mean(post$xi)
  )
  
  # Draws for y* (from Stan posterior)
  ystar_draw <- post$ystar
  
  # Propagate to mass with ABG + residual uncertainty (predictive uplift)
  take <- min(N_PROP, length(ystar_draw))
  idx  <- sample.int(length(ystar_draw), size = take, replace = TRUE)
  Z    <- ystar_draw[idx] - c0
  ABG  <- MASS::mvrnorm(take, mu = c(alpha,beta,gamma), Sigma = Vabg)
  eps  <- rnorm(take, 0, sigma_e)
  mlog <- ABG[,1] + ABG[,2]*Z + ABG[,3]*Z^2 + eps
  mass_draws <- exp(mlog) / 1e6
  
  # Summaries: MAP via KDE, and one-sided uppers
  kd <- density(mass_draws, n = 4096, adjust = 1.0)
  mass_map   <- kd$x[which.max(kd$y)]
  mass_u95   <- as.numeric(quantile(mass_draws, probs = ci_upper))
  mass_u90   <- as.numeric(quantile(mass_draws, probs = ci_upper_alt))
  
  results <- dplyr::bind_rows(
    results,
    tibble::tibble(
      year = t, n_ex = n_ex, u0 = u0,
      map_mass_tons = mass_map,
      mass_upper95_tons = mass_u95,
      mass_upper90_tons = mass_u90
    )
  )
  
  message(sprintf("Year %d | n_ex=%d | MAP=%.1f t | 95%% upper=%.1f t",
                  t, n_ex, mass_map, mass_u95))
}

stopifnot(nrow(results) > 0)
saveRDS(results, "results_sequential_stan_3.RDS")
results <- readRDS("results_sequential_stan_3.RDS")

# =============================================================================
# Overlays (top-10 + all exceedances) — unchanged
# =============================================================================
map_logM_cq <- function(logC) {
  a <- alpha; b <- beta; g <- gamma; c0v <- c0
  a + b*(logC - c0v) + g*(logC - c0v)^2
}
df$mass_quad_all <- exp(map_logM_cq(df$log_circ)) / 1e6
exceed_pts <- df |> dplyr::filter(log_circ > u0) |> dplyr::select(year, mass_quad_all)

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
species_colors <- c(
  "Ruyangosaurus"="#a6cee3","Turiasaurus"="#1f78b4","Yunmenglong"="#b2df8a",
  "Australotitan"="#33a02c","Patagotitan"="#fb9a99","Dreadnoughtus"="#e31a1c",
  "Notocolossus"="#fdbf6f","Chucarosaurus"="#ff7f00","Brachiosaurus"="#cab2d6",
  "Argentinosaurus"="#6a3d9a","Extremosaurus"="firebrick"
)

ext_map <- results |> dplyr::transmute(year, mass = map_mass_tons, species = "Extremosaurus")
species_colors2 <- c(species_colors, "Extremosaurus (MAP)" = "firebrick")
legend_order <- c(
  "Ruyangosaurus","Turiasaurus","Yunmenglong","Australotitan",
  "Patagotitan","Dreadnoughtus","Notocolossus","Chucarosaurus",
  "Brachiosaurus","Argentinosaurus","Extremosaurus"
)

p_mass <- ggplot() +
  geom_point(data = exceed_pts, aes(x = year, y = mass_quad_all),
             color = "gray70", alpha = 0.5, size = 1.8, shape = 16, show.legend = FALSE) +
  geom_errorbar(data = top10_data,
                aes(x = year, ymin = pi_lower_quad, ymax = pi_upper_quad, color = species),
                width = 0.5, linewidth = 1, show.legend = FALSE) +
  geom_line(data = top10_data, aes(x = year, y = mass_quad, color = species),
            linewidth = 1.1, alpha = 0.9, show.legend = TRUE) +
  geom_point(data = top10_data, aes(x = year, y = mass_quad, color = species),
             size = 2.6, shape = 16, show.legend = FALSE) +
  geom_line(data = ext_map, aes(x = year, y = mass, color = species),
            linewidth = 1.3, show.legend = TRUE) +
  labs(x = "Publication Year", y = "Estimated Mass [tons]", color = "Specimens") +
  scale_color_manual(values = species_colors2, limits = legend_order, breaks = legend_order) +
  scale_x_continuous(limits = c(1892, max(c(results$year, top10_data$year, df$year), na.rm = TRUE))) +
  ylim(0, 200) +
  theme_science_polished +
  guides(color = guide_legend(override.aes = list(linewidth = 1.3, alpha = 1)))
p_mass
ggsave(file.path(FIG_DIR, "mass_sequential.png"),
       p_mass, dpi = 600, w = 7, h = 5, units = "in")

p_mass_supplement <- ggplot() +
  geom_point(data = exceed_pts, aes(x = year, y = mass_quad_all),
             color = "gray70", alpha = 0.5, size = 1.8, shape = 16, show.legend = FALSE) +
  geom_line(data = results, aes(x = year, y = map_mass_tons),
            color = "purple", linewidth = 1.3) +
  geom_line(data = results, aes(x = year, y = mass_upper90_tons),
            color = "firebrick", linewidth = 1.3) +
  geom_line(data = results, aes(x = year, y = mass_upper95_tons),
            color = "orange", linewidth = 1.3) +
  labs(x = "Publication Year", y = "Estimated Mass [tons]") +
  scale_x_continuous(limits = c(1892, max(c(results$year, top10_data$year, df$year), na.rm = TRUE))) +
  ylim(0, 800) +
  theme_science_polished
p_mass_supplement
ggsave(file.path(FIG_DIR, "mass_sequential_supplement_STAN_collapsed.png"),
       p_mass_supplement, dpi = 600, w = 7, h = 5, units = "in")
# --- Legend + colors for log plot: top10 + Extremosaurus ---
legend_order_log   <- c(legend_order_top10, "Extremosaurus")
species_colors_log <- c(species_colors_top10, "Extremosaurus" = "firebrick")

# Ensure ext_map has the right species label
ext_map <- results |>
  dplyr::transmute(year, mass = map_mass_tons, species = "Extremosaurus")

# --- Log (natural) scale plot with Extremosaurus in legend ---
p_mass_log <- ggplot() +
  # background points (all exceedances; no legend)
  geom_point(data = exceed_pts, aes(x = year, y = mass_quad_all),
             color = "gray70", alpha = 0.5, size = 1.8, shape = 16, show.legend = FALSE) +
  # top-10 bands + lines + points
  geom_errorbar(data = top10_data,
                aes(x = year, ymin = pi_lower_quad, ymax = pi_upper_quad, color = species),
                width = 0.5, linewidth = 1, show.legend = FALSE) +
  geom_line(data = top10_data,
            aes(x = year, y = mass_quad, color = species),
            linewidth = 1.1, alpha = 0.9, show.legend = TRUE) +
  geom_point(data = top10_data,
             aes(x = year, y = mass_quad, color = species),
             size = 2.6, shape = 16, show.legend = FALSE) +
  # Extremosaurus MAP (goes in legend)
  geom_line(data = ext_map,
            aes(x = year, y = mass, color = species),
            linewidth = 1.3, show.legend = TRUE) +
  labs(x = "Publication Year", y = "Log Mass [grams]", color = "Specimens") +
  scale_x_continuous(limits = c(1892, max(c(results$year, top10_data$year, df$year), na.rm = TRUE))) +
  scale_y_continuous() +  # natural log
  scale_color_manual(values = species_colors_log,
                     limits = legend_order_log,
                     breaks = legend_order_log) +
  theme_science_polished +
  guides(color = guide_legend(override.aes = list(linewidth = 1.3, alpha = 1)))

p_mass_log
ggsave(file.path(FIG_DIR, "mass_sequential_log.png"),
       p_mass_log, dpi = 600, w = 7, h = 5, units = "in")


p_mass_log
ggsave(file.path(FIG_DIR, "mass_sequential_uniform_STAN_log_collapsed.png"),
       p_mass_log, dpi = 600, w = 7, h = 5, units = "in")
ggsave(file.path(FIG_DIR, "mass_sequential_uniform_STAN_log_collapsed.pdf"),
       p_mass_log, device = cairo_pdf, w = 7, h = 5, units = "in")

ggsave(file.path(FIG_DIR, "mass_sequential_uniform_STAN_log_collapsed.png"),
       p_mass_log, dpi = 600, w = 7, h = 5, units = "in")
ggsave(file.path(FIG_DIR, "mass_sequential_uniform_STAN_log_collapsed.pdf"),
       p_mass_log, device = cairo_pdf, w = 7, h = 5, units = "in")

