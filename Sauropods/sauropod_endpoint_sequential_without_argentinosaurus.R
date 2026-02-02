library(ggplot2)
library(pracma)       # For trapz
library(HDInterval)   # For HPDI
library(readxl)
library(dplyr)
library(knitr)
library(scales)
theme_science_polished <- theme_minimal(base_family = "Arial", base_size = 12) +
  theme(
    axis.title = element_text(size = 16, face = "bold"),
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
df <- read_excel("Data/DEmic23_updated_Supplemental_Data_withPubYear.xlsx")
df <- df[, c("hum+fem circ (mm)", "publication year")]
colnames(df) <- c("sum_circ_mm", "year")
df <- na.omit(df)
df$log_circ <- log(as.double(df$sum_circ_mm))

# Exclude Argentinosaurus (2016 mm, year 1993)
df <- df[!(df$sum_circ_mm == 2016.222 & df$year == 1993), ]

# Threshold and evaluation years
threshold <- 2.94*log(10)
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
    samples=samples,
    sampled_mass_lin_tons=sampled_mass_lin_tons,
    sampled_mass_quad_tons=sampled_mass_quad_tons,
    
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


saveRDS(results, 'results_without_argentinosaurus.RDS')
results <- readRDS('results_without_argentinosaurus.RDS')
samples <- results$samples
sampled_mass_lin_tons <- results$sampled_mass_lin_tons
sampled_mass_quad_tons <- results$sampled_mass_quad_tons
# Remove the heavy samples columns
results$samples <- NULL
results$sampled_mass_lin_tons <- NULL
results$sampled_mass_quad_tons <- NULL
# Plot 1: Endpoint estimate over time
exceedances <- df[df$log_circ > threshold, c("year", "log_circ")]

p_CFH <- ggplot(results, aes(x = year)) +
  geom_ribbon(aes(ymin = exp(hpd_lower), ymax = exp(hpd_upper)), fill = "lightblue", alpha = 0.5) +
  geom_line(aes(y = exp(map)), color = "black", size = 1.2) +
  geom_point(aes(y = exp(map)), color = "blue", size = 2) +
  geom_point(data = exceedances, aes(x = year, y = exp(log_circ)), shape = 4, color = "red", size = 2) +
  labs(
       x = "Publication Year",
       y = expression(C[F+H]~"[mm]")) +
  theme_science_polished
p_CFH

p_log_CFH <- ggplot(results, aes(x = year)) +
  geom_ribbon(aes(ymin = hpd_lower, ymax = hpd_upper), fill = "lightblue", alpha = 0.5) +
  geom_line(aes(y = map), color = "black", size = 1.2) +
  geom_point(aes(y = map), color = "blue", size = 2) +
  geom_point(data = exceedances, aes(x = year, y = log_circ), shape = 4, color = "red", size = 2) +
  labs(
    x = "Publication Year",
    y = expression(log(C[F+H])~"[log(mm)]")
  ) +
  scale_x_continuous(limits = c(1870, max(results$year))) +
  theme_science_polished
p_log_CFH 

# Plot 2: Threshold exceedances
p_exceedances <- ggplot(results, aes(x = year, y = n)) +
  geom_col(fill = "gray40") +
  labs(title = "Number of Threshold-Exceeding Specimens Over Time",
       x = "Year", y = "Count > u") +
  theme_science_polished
p_exceedances
# Data for top 5 sauropods with quadratic mass estimates
top10_data <- data.frame(
  species = c("Ruyangosaurus", "Turiasaurus", "Yunmenglong", "Australotitan", "Patagotitan",
              "Dreadnoughtus", "Notocolossus", "Chucarosaurus", "Brachiosaurus", "Argentinosaurus"),
  year = c(2017, 2006, 2013, 2021, 2017, 2014, 2016, 2024, 1903, 1993),
  sum_circ_mm = c(1639, 1672, 1687, 1692, 1694, 1695, 1704, 1827, 1870, 2016),
  mass_quad = c(43, 46, 47, 47, 47, 47, 48, 57, 61, 74),
  pi_lower_quad = c(20, 20, 19, 19, 18, 20, 20, 24, 25, 30),
  pi_upper_quad = c(104, 110, 112, 113, 114, 114, 116, 140, 150, 185)
)

# Add labels
results$label <- "Extremosaurus"
top10_data$label <- top10_data$species
top10_data$type <- "Specimen"


extremo_data <- results %>%
  dplyr::mutate(
    type = "Extremosaurus",
    mass_quad = map_mass_quad,
    pi_lower_quad = mass_hpdi_lower_quad,
    pi_upper_quad = mass_hpdi_upper_quad,
    species = "Extremosaurus"
  ) %>%
  dplyr::select(year, mass_quad, pi_lower_quad, pi_upper_quad, species, type)

overlay_data <- top10_data %>%
  dplyr::select(year, mass_quad, pi_lower_quad, pi_upper_quad, species, type)

plot_data <- bind_rows(overlay_data, extremo_data)

plot_data$species <- factor(plot_data$species, levels = c(
  "Ruyangosaurus", "Turiasaurus", "Yunmenglong", "Australotitan", "Patagotitan",
  "Dreadnoughtus", "Notocolossus", "Chucarosaurus", "Brachiosaurus", "Argentinosaurus", "Extremosaurus"
))


species_colors <- c(
  "Ruyangosaurus" = "#1b9e77",
  "Turiasaurus" = "#d95f02",
  "Yunmenglong" = "#7570b3",
  "Australotitan" = "#e7298a",
  "Patagotitan" = "#66a61e",
  "Dreadnoughtus" = "#e6ab02",
  "Notocolossus" = "#a6761d",
  "Chucarosaurus" = "#666666",
  "Brachiosaurus" = "#a6cee3",
  "Argentinosaurus" = "#fb9a99",
  "Extremosaurus" = "firebrick"
)

p_mass <- ggplot() +
  geom_ribbon(data = results, aes(x = year, ymin = mass_hpdi_lower_quad, ymax = mass_hpdi_upper_quad), 
              fill = "mistyrose", alpha = 0.5) +
  geom_line(data = plot_data, aes(x = year, y = mass_quad, color = species), linewidth = 1.2) +
  geom_point(data = plot_data, aes(x = year, y = mass_quad, color = species), 
             size = 2.5, shape = 16, show.legend = FALSE) +
  geom_errorbar(data = overlay_data, 
                aes(x = year, ymin = pi_lower_quad, ymax = pi_upper_quad, color = species), 
                width = 0.5, linewidth = 1, show.legend=FALSE) +
  labs(
    x = "Publication Year",
    y = "Estimated Mass [tons]",
    color = "Specimen"
  ) +
  scale_color_manual(values = species_colors) +
  scale_x_continuous(limits = c(1870, max(results$year))) +
  theme_science_polished
p_mass

# now with x axis removed (Mike asked this)

p_mass_no_title_text_ticks <- p_mass +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
p_mass_no_title_text_ticks

p_mass_no_title <- p_mass +
  theme(axis.title.x = element_blank())
p_mass_no_title

p_mass_no_text <- p_mass +
  theme(axis.text.x = element_blank())
p_mass_no_text

p_mass_no_ticks <- p_mass +
  theme(axis.ticks.x = element_blank())
p_mass_no_ticks

p_mass_no_title_text <- p_mass +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank())
p_mass_no_title_text

p_mass_no_title_ticks <- p_mass +
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank())
p_mass_no_title_ticks

p_mass_no_text_ticks <- p_mass +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
p_mass_no_text_ticks
# Full plot with Extremosaurus estimate (ribbon and line)
p_log_mass <- ggplot() +
  geom_ribbon(data = results, aes(x = year, ymin = log(mass_hpdi_lower_quad), ymax = log(mass_hpdi_upper_quad)), 
              fill = "mistyrose", alpha = 0.5) +
  geom_line(data = plot_data, aes(x = year, y = log(mass_quad), color = species), linewidth = 1.2) +
  geom_point(data = plot_data, aes(x = year, y = log(mass_quad), color = species), 
             size = 2.5, shape = 16, show.legend = FALSE) +
  geom_errorbar(data = overlay_data, 
                aes(x = year, ymin = log(pi_lower_quad), ymax = log(pi_upper_quad), color = species), 
                width = 0.5, linewidth = 1, show.legend=FALSE) +
  labs(
    x = "Publication Year",
    y = "Estimated log(Mass) [tons]",
    color = "Specimen"
  ) +
  scale_color_manual(values = species_colors) +
  scale_x_continuous(limits = c(1870, max(results$year))) +
  theme_science_polished

## --- remove Extremosaurus before plotting ---------------------------------
overlay_data_no_extremo <- overlay_data %>% 
  filter(species != "Extremosaurus")

species_colors_no_extremo <- species_colors[names(species_colors) != "Extremosaurus"]
species_levels_no_extremo <- levels(plot_data$species)[levels(plot_data$species) != "Extremosaurus"]

## --- log-mass plot without Extremosaurus ----------------------------------
p_log_mass_no_extremo <- ggplot() +
  geom_point(data = overlay_data_no_extremo,
             aes(x = year, y = log(mass_quad), color = species),
             size = 2.5, shape = 16, show.legend = TRUE) +
  geom_errorbar(data = overlay_data_no_extremo,
                aes(x = year, ymin = log(pi_lower_quad), ymax = log(pi_upper_quad),
                    color = species),
                width = 0.5, linewidth = 1, show.legend = TRUE) +
  labs(x = "Publication Year",
       y = "Estimated log(Mass) [tons]",
       color = "Specimen") +
  scale_color_manual(values = species_colors_no_extremo,
                     breaks = species_levels_no_extremo,   
                     drop   = TRUE)               +  

  scale_x_continuous(limits = c(1870, max(results$year))) +
  theme_science_polished

p_log_mass_no_extremo

# now with x axis removed (Mike asked this)

p_log_mass_no_title_text_ticks <- p_log_mass +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
p_log_mass_no_title_text_ticks

p_log_mass_no_title <- p_log_mass +
  theme(axis.title.x = element_blank())
p_log_mass_no_title

p_log_mass_no_text <- p_log_mass +
  theme(axis.text.x = element_blank())
p_log_mass_no_text

p_log_mass_no_ticks <- p_log_mass +
  theme(axis.ticks.x = element_blank())
p_log_mass_no_ticks

p_log_mass_no_title_text <- p_log_mass +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank())
p_log_mass_no_title_text

p_log_mass_no_title_ticks <- p_log_mass +
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank())
p_log_mass_no_title_ticks

p_log_mass_no_text_ticks <- p_log_mass +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
p_log_mass_no_text_ticks

# Now for no_extremosaurus plots
p_log_mass_no_extremo_no_title_text_ticks <- p_log_mass_no_extremo +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
p_log_mass_no_extremo_no_title_text_ticks

p_log_mass_no_extremo_no_title <- p_log_mass_no_extremo +
  theme(axis.title.x = element_blank())
p_log_mass_no_extremo_no_title

p_log_mass_no_extremo_no_text <- p_log_mass_no_extremo +
  theme(axis.text.x = element_blank())
p_log_mass_no_extremo_no_text

p_log_mass_no_extremo_no_ticks <- p_log_mass_no_extremo +
  theme(axis.ticks.x = element_blank())
p_log_mass_no_extremo_no_ticks

p_log_mass_no_extremo_no_title_text <- p_log_mass_no_extremo +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank())
p_log_mass_no_extremo_no_title_text

p_log_mass_no_extremo_no_title_ticks <- p_log_mass_no_extremo +
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank())
p_log_mass_no_extremo_no_title_ticks

p_log_mass_no_extremo_no_text_ticks <- p_log_mass_no_extremo +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
p_log_mass_no_extremo_no_text_ticks

ggplot(results, aes(x = year)) +
  geom_line(aes(y = map_mass_quad, color = "Quadratic"), linewidth = 1.2) +
  geom_line(aes(y = map_mass_lin, color = "Linear"), linewidth = 1.2, linetype = "dashed") +
  geom_ribbon(aes(ymin = mass_hpdi_lower_quad, ymax = mass_hpdi_upper_quad, fill = "Quadratic"), alpha = 0.2) +
  geom_ribbon(aes(ymin = mass_hpdi_lower_lin, ymax = mass_hpdi_upper_lin, fill = "Linear"), alpha = 0.2) +
  labs(x = "Year", y = "Estimated Mass [tons]", color = "Model", fill = "Model") +
  theme_science_polished

####################
# Plot Posteriors: 
####################
map_cf <- tail(results$map, 1)
hpdi_cf_lower <- tail(results$hpd_lower, 1)
hpdi_cf_upper <- tail(results$hpd_upper, 1)
hpdi_cf <- c(hpdi_cf_lower, hpdi_cf_upper)

map_lin <- tail(results$map_mass_lin,1)
hpdi_lin_lower <- tail(results$mass_hpdi_lower_lin, 1)
hpdi_lin_upper <- tail(results$mass_hpdi_upper_lin, 1)
hpdi_lin <- c(hpdi_lin_lower, hpdi_lin_upper)

map_quad <- tail(results$map_mass_quad, 1)
hpdi_quad_lower <- tail(results$mass_hpdi_lower_quad, 1)
hpdi_quad_upper <- tail(results$mass_hpdi_upper_quad, 1)
hpdi_quad <- c(hpdi_quad_lower, hpdi_quad_upper)


# Transform posterior samples to original scale
samples_exp <- exp(samples)
map_cf_exp <- exp(map_cf)
hpdi_cf_exp <- exp(hpdi_cf)

# Plot on original scale
p_CFH_posterior <- ggplot(data.frame(cfh = samples_exp), aes(x = cfh)) +
  geom_histogram(aes(y = ..density..), bins = 40, fill = "skyblue", color = "black", alpha = 0.6) +
  geom_density(color = "darkblue", size = 1.2) +
  geom_vline(xintercept = map_cf_exp, linetype = "dashed", color = "purple", size = 1.2) +
  geom_vline(xintercept = hpdi_cf_exp, linetype = "dotdash", color = "orange", size = 1.2) +
  labs(x = expression(C[F+H]~"[mm] endpoint"), y = "Density") +
  scale_x_continuous(limits = c(1800, 2700)) +  # adjust as needed
  theme_science_polished
p_CFH_posterior

# ---- Plot posterior of linear mass ----
p_mass_lin_posterior <- ggplot(data.frame(mass = sampled_mass_lin_tons), aes(x = mass)) +
  geom_histogram(aes(y = ..density..), bins = 40, fill = "lightgreen", color = "black", alpha=.6) +
  geom_density(color = "darkgreen", size = 1.2) +
  geom_vline(xintercept = map_lin, linetype = "dashed", color = "purple", size=1.2) +
  geom_vline(xintercept = hpdi_lin, linetype = "dotdash", color = "orange", size=1.2) +
  labs(x = "Mass [tons]", y = "Density") +
  theme_science_polished + 
  scale_x_continuous(limits=c(0,300))
p_mass_lin_posterior

# ---- Plot posterior of quadratic mass ----
p_mass_quad_posterior <- ggplot(data.frame(mass = sampled_mass_quad_tons), aes(x = mass)) +
  geom_histogram(aes(y = ..density..), bins = 50, fill = "lightgreen", color = "black", alpha=.6) +
  geom_density(color = "darkgreen", size = 1.2) +
  geom_vline(xintercept = map_quad, linetype = "dashed", color = "purple", size=1.2) +
  geom_vline(xintercept = hpdi_quad, linetype = "dotdash", color = "orange", size=1.2) +
  labs(x = "Mass [tons]", y = "Density") +
  theme_science_polished + 
  scale_x_continuous(limits=c(0,300))
p_mass_quad_posterior

plot_list <- list(
  p_CFH,
  p_log_CFH,
  p_exceedances,
  p_mass,
  p_mass_no_title_text_ticks,
  p_mass_no_title,
  p_mass_no_text,
  p_mass_no_ticks,
  p_mass_no_title_text,
  p_mass_no_title_ticks,
  p_mass_no_text_ticks,
  p_log_mass,
  p_log_mass_no_title_text_ticks,
  p_log_mass_no_title,
  p_log_mass_no_text,
  p_log_mass_no_ticks,
  p_log_mass_no_title_text,
  p_log_mass_no_title_ticks,
  p_log_mass_no_text_ticks,
  p_log_mass_no_extremo,
  p_log_mass_no_extremo_no_title_text_ticks,
  p_log_mass_no_extremo_no_title,
  p_log_mass_no_extremo_no_text,
  p_log_mass_no_extremo_no_ticks,
  p_log_mass_no_extremo_no_title_text,
  p_log_mass_no_extremo_no_title_ticks,
  p_log_mass_no_extremo_no_text_ticks,
  p_CFH_posterior,
  p_mass_lin_posterior,
  p_mass_quad_posterior
)

file_names <- c(
  "sequential_endpoint_without_Argentinosaurus",
  "sequential_log_endpoint_without_Argentinosaurus",
  "exceedances_over_time_without_Argentinosaurus",
  "sequential_mass_endpoint_without_Argentinosaurus",
  "sequential_mass_endpoint_without_Argentinosaurus_no_title_text_ticks",
  "sequential_mass_endpoint_without_Argentinosaurus_no_title",
  "sequential_mass_endpoint_without_Argentinosaurus_no_text",
  "sequential_mass_endpoint_without_Argentinosaurus_no_ticks",
  "sequential_mass_endpoint_without_Argentinosaurus_no_title_text",
  "sequential_mass_endpoint_without_Argentinosaurus_no_title_ticks",
  "sequential_mass_endpoint_without_Argentinosaurus_no_text_ticks",
  "sequential_log_mass_endpoint_without_Argentinosaurus",
  "sequential_log_mass_endpoint_without_Argentinosaurus_no_title_text_ticks",
  "sequential_log_mass_endpoint_without_Argentinosaurus_no_title",
  "sequential_log_mass_endpoint_without_Argentinosaurus_no_text",
  "sequential_log_mass_endpoint_without_Argentinosaurus_no_ticks",
  "sequential_log_mass_endpoint_without_Argentinosaurus_no_title_text",
  "sequential_log_mass_endpoint_without_Argentinosaurus_no_title_ticks",
  "sequential_log_mass_endpoint_without_Argentinosaurus_no_text_ticks",
  "sequential_log_mass_endpoint_without_Argentinosaurus_no_extremo",
  "sequential_log_mass_endpoint_without_Argentinosaurus_no_extremo_no_title_text_ticks",
  "sequential_log_mass_endpoint_without_Argentinosaurus_no_extremo_no_title",
  "sequential_log_mass_endpoint_without_Argentinosaurus_no_extremo_no_text",
  "sequential_log_mass_endpoint_without_Argentinosaurus_no_extremo_no_ticks",
  "sequential_log_mass_endpoint_without_Argentinosaurus_no_extremo_no_title_text",
  "sequential_log_mass_endpoint_without_Argentinosaurus_no_extremo_no_title_ticks",
  "sequential_log_mass_endpoint_without_Argentinosaurus_no_extremo_no_text_ticks",
  "posterior_logCFH_without_Argentinosaurus",
  "posterior_mass_lin_without_Argentinosaurus",
  "posterior_mass_quad_without_Argentinosaurus"
)

# Output directory
output_dir <- "Figures_sauropods"
if (!dir.exists(output_dir)) dir.create(output_dir)


# Save all plots as PNG and PDF (via cairo_pdf)
for (i in seq_along(plot_list)) {
  ggsave(
    filename = file.path(output_dir, paste0(file_names[i], ".png")),
    plot = plot_list[[i]],
    device = "png", dpi = 600,
    width = 7, height = 5, units = "in"
  )
  ggsave(
    filename = file.path(output_dir, paste0(file_names[i], ".pdf")),
    plot = plot_list[[i]],
    device = cairo_pdf, 
    width = 7, height = 5, units = "in"
  )
}






