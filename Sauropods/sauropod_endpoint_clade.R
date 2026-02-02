library(ggplot2)
library(pracma)       # For trapz
library(HDInterval)   # For HPDI
library(readxl)
library(dplyr)
library(knitr)
library(scales)

theme_science_polished <- theme_minimal(base_family = "Arial", base_size = 12) +
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

# Load and preprocess the dataset
df_raw <- read_excel("Data/DEmic23_updated_Supplemental_Data_withPubYear.xlsx")
taxon_summary <- df_raw %>% count(taxon) %>% arrange(desc(n))
kable(taxon_summary, caption = "Number of specimens per taxon (clade abbreviation)")


# --- Clade-wise endpoint estimation ---

taxon_list <- unique(df_raw$taxon)
clade_results <- list()

for (clade in taxon_list) {
  df_clade <- df_raw %>% filter(taxon == clade) %>% select(`hum+fem circ (mm)`) %>% na.omit()
  df_clade$`hum+fem circ (mm)` <- as.numeric(df_clade$`hum+fem circ (mm)`)
  if (nrow(df_clade) < 5) next
  
  y_excess <- log(df_clade$`hum+fem circ (mm)`)
  y_excess <- y_excess[y_excess > log(873)]
  if (length(y_excess) < 5) next
  
  Y_max <- max(y_excess)
  eps <- max(1e-6, 1e-4 * Y_max)
  y_grid <- seq(Y_max + eps, Y_max + 1, length.out = 1000)
  xi_grid <- seq(-0.6, -1e-3, length.out = 1000)
  
  loglik <- function(y_star, xi) {
    if (xi >= 0 || any(y_star <= y_excess)) return(-Inf)
    term1 <- length(y_excess) * log(y_star - log(873)) / xi
    term2 <- -(1 / xi + 1) * sum(log(y_star - y_excess))
    term3 <- -length(y_excess) * log(-xi)
    return(term1 + term2 + term3)
  }
  
  marginal_log_post <- sapply(y_grid, function(y_star) {
    loglik_vals <- sapply(xi_grid, function(xi) loglik(y_star, xi))
    integrand <- exp(loglik_vals - max(loglik_vals))
    log(trapz(xi_grid, integrand)) + max(loglik_vals)
  })
  
  prior <- ifelse(y_grid > Y_max, 1, 0)
  posterior_unnorm <- exp(marginal_log_post - max(marginal_log_post)) * prior
  posterior <- posterior_unnorm / trapz(y_grid, posterior_unnorm)
  
  samples <- sample(y_grid, size = 10000, replace = TRUE, prob = posterior)
  map <- y_grid[which.max(posterior)]
  hpd <- hdi(samples, credMass = 0.95)
  
  clade_results[[clade]] <- list(
    clade = clade,
    map = exp(map),
    hpd_lower = exp(hpd[1]),
    hpd_upper = exp(hpd[2]),
    n = length(y_excess),
    samples = exp(samples)
  )
}

# Convert list to data frame
clade_summary <- do.call(rbind, lapply(clade_results, function(x) {
  data.frame(
    clade = x$clade,
    n = x$n,
    map = x$map,
    hpd_lower = x$hpd_lower,
    hpd_upper = x$hpd_upper
  )
}))

print(kable(clade_summary, digits = 1, caption = "Endpoint estimates (in mm) per clade"))

# Plot posterior distributions per clade
library(tidyr)
library(purrr)

posterior_df <- map_dfr(clade_results, function(x) {
  data.frame(value = x$samples, clade = x$clade)
})

p_clade_posteriors <- ggplot(posterior_df, aes(x = value, fill = clade)) +
  geom_density(alpha = 0.4) +
  labs(x = expression(C[F+H]~"endpoint [mm]"), y = "Density", title = "Posterior distributions per clade") +
  theme_science_polished

print(p_clade_posteriors)

# --- Main dataset for sequential endpoint analysis ---
df <- df_raw[, c("hum+fem circ (mm)", "publication year")]
colnames(df) <- c("sum_circ_mm", "year")
df <- na.omit(df)
df$log_circ <- log(as.double(df$sum_circ_mm))

# Threshold and evaluation years
threshold <- log(873)
years_eval <- sort(unique(df$year[df$year >= 1903]))

# Load regression model coefficients
quad_coeffs <- readRDS("centered_quadratic_coefficients.rds")
lin_coeffs <- readRDS("linear_model_coefficients.rds")

# Sampling settings
n_samples <- 10000
ci_level <- 0.95
lower_prob <- (1 - ci_level) / 2
upper_prob <- 1 - lower_prob

# Sample regression coefficients
sampled_alpha <- rnorm(n_samples, mean = quad_coeffs$alpha, sd = quad_coeffs$alpha_se)
sampled_beta  <- rnorm(n_samples, mean = quad_coeffs$beta,  sd = quad_coeffs$beta_se)
sampled_gamma <- rnorm(n_samples, mean = quad_coeffs$gamma, sd = quad_coeffs$gamma_se)
sampled_alpha_lin <- rnorm(n_samples, mean = lin_coeffs$alpha, sd = lin_coeffs$alpha_se)
sampled_beta_lin <- rnorm(n_samples, mean = lin_coeffs$beta, sd = lin_coeffs$beta_se)

# Storage
results <- data.frame()

# Main loop
for (t in years_eval) {
  y_sub <- df$log_circ[df$year <= t]
  y_excess <- y_sub[y_sub > threshold]
  
  if (length(y_excess) < 5) next
  
  n_excess <- length(y_excess)
  Y_max <- max(y_excess)
  eps <- max(1e-6, 1e-4 * Y_max)
  y_grid <- seq(Y_max + eps, Y_max + 1, length.out = 1000)
  xi_grid <- seq(-0.6, -1e-3, length.out = 1000)
  
  loglik <- function(y_star, xi) {
    if (xi >= 0 || any(y_star <= y_excess)) return(-Inf)
    term1 <- n_excess * log(y_star - threshold) / xi
    term2 <- -(1 / xi + 1) * sum(log(y_star - y_excess))
    term3 <- -n_excess * log(-xi)
    return(term1 + term2 + term3)
  }
  
  marginal_log_post <- sapply(y_grid, function(y_star) {
    loglik_vals <- sapply(xi_grid, function(xi) loglik(y_star, xi))
    integrand <- exp(loglik_vals - max(loglik_vals))
    log(trapz(xi_grid, integrand)) + max(loglik_vals)
  })
  
  prior <- ifelse(y_grid > Y_max, 1, 0)
  posterior_unnorm <- exp(marginal_log_post - max(marginal_log_post)) * prior
  posterior <- posterior_unnorm / trapz(y_grid, posterior_unnorm)
  
  samples <- sample(y_grid, size = n_samples, replace = TRUE, prob = posterior)
  map <- y_grid[which.max(posterior)]
  hpd <- hdi(samples, credMass = ci_level)
  
  # ---- Quadratic Model Inference ----
  centered_samples <- samples - quad_coeffs$mean_log_sum_circ
  sampled_log_mass_quad <- rnorm(n_samples,
                                 mean = sampled_alpha + sampled_beta * centered_samples + sampled_gamma * centered_samples^2,
                                 sd = quad_coeffs$resid_sd
  )
  sampled_mass_quad_tons <- exp(sampled_log_mass_quad) / 1e6
  
  mass_median_quad <- median(sampled_mass_quad_tons)
  mass_hpdi_quad <- hdi(sampled_mass_quad_tons, credMass = ci_level)
  
  map_centered <- map - quad_coeffs$mean_log_sum_circ
  map_log_mass_quad <- quad_coeffs$alpha + quad_coeffs$beta * map_centered + quad_coeffs$gamma * map_centered^2
  map_mass_quad <- exp(map_log_mass_quad) / 1e6
  
  # ---- Linear Model Inference ----
  sampled_log_mass_lin <- rnorm(n_samples,
                                mean = sampled_alpha_lin + sampled_beta_lin * samples,
                                sd = lin_coeffs$resid_sd
  )
  sampled_mass_lin_tons <- exp(sampled_log_mass_lin) / 1e6
  
  mass_median_lin <- median(sampled_mass_lin_tons)
  mass_hpdi_lin <- hdi(sampled_mass_lin_tons, credMass = ci_level)
  
  map_log_mass_lin <- lin_coeffs$alpha + lin_coeffs$beta * map
  map_mass_lin <- exp(map_log_mass_lin) / 1e6
  
  # Store results
  results <- rbind(results, data.frame(
    year = t,
    map = map,
    hpd_lower = hpd[1],
    hpd_upper = hpd[2],
    n = length(y_excess),
    
    # Quadratic estimates
    map_mass_quad = map_mass_quad,
    mass_median_quad = mass_median_quad,
    mass_hpdi_lower_quad = mass_hpdi_quad[1],
    mass_hpdi_upper_quad = mass_hpdi_quad[2],
    
    # Linear estimates
    map_mass_lin = map_mass_lin,
    mass_median_lin = mass_median_lin,
    mass_hpdi_lower_lin = mass_hpdi_lin[1],
    mass_hpdi_upper_lin = mass_hpdi_lin[2]
  ))
}
