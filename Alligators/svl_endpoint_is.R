library(ggplot2)
library(extRemes)
library(ismev)
library(numDeriv)
library(dplyr)
library(MASS)
library(readxl)
library(pracma)
library(HDInterval)
library(grid)
library(gridExtra)
library(coda)
library(truncnorm)

# Ensure the output directory exists
output_dir <- "Figures_alligators"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Helper function to save plot in TIFF and PNG formats
save_plot <- function(plot, filename, width = 7, height = 5, dpi = 600) {
  ggsave(file.path(output_dir, paste0(filename, ".tiff")), plot = plot,
         dpi = dpi, width = width, height = height, units = "in", compression = "lzw")
  ggsave(file.path(output_dir, paste0(filename, ".png")), plot = plot,
         dpi = dpi, width = width, height = height, units = "in")
}


##############################
### LOAD DATA ###
##############################

df <- read_excel("Data/experimental_alligator_harvest_woodward.xlsx")

# Keep males only, remove NAs
df_males <- df %>% filter(Sex == "M", !is.na(SVL), !is.na(TL))

log_SVL <- log(df_males$SVL)
log_TL <- log(df_males$TL)

##############################
### POT ANALYSIS ON LOG SVL ###
##############################

u_0 <- 4.9
thresholds <- seq(u_0, 5.3, length.out = 50)

scale_params <- shape_params <- shape_params_lower <- shape_params_upper <- numeric(length(thresholds))
scale_params_lower <- scale_params_upper <- adjusted_scale_params <- adjusted_scale_params_lower <- adjusted_scale_params_upper <- numeric(length(thresholds))

for (i in seq_along(thresholds)) {
  fit <- fevd(log_SVL, threshold = thresholds[i], type = "GP")
  scale_params[i] <- fit$results$par["scale"]
  shape_params[i] <- fit$results$par["shape"]
  adjusted_scale_params[i] <- scale_params[i] - shape_params[i] * (thresholds[i] - u_0)
  
  tryCatch({
    ci <- ci(fit, type = 'parameter', alpha = 0.05)
    cov <- summary(fit)$cov.theta
    shape_params_lower[i] <- ci[2,1]
    shape_params_upper[i] <- ci[2,3]
    scale_params_lower[i] <- ci[1,1]
    scale_params_upper[i] <- ci[1,3]
    var_adj <- cov[1,1] + cov[2,2] * (thresholds[i] - u_0)^2
    adjusted_scale_params_lower[i] <- adjusted_scale_params[i] - qnorm(0.975) * sqrt(var_adj)
    adjusted_scale_params_upper[i] <- adjusted_scale_params[i] + qnorm(0.9) * sqrt(var_adj)
  }, error = function(e) {})
}

number <- 39
threshold_opt <- thresholds[number]
threshold_opt <- log(187)
#threshold_opt <- log(187)

##############################
### PLOT PARAMETER STABILITY ###
##############################

# Polished theme for Science-style figures
theme_science_polished <- theme_minimal(base_family = "Arial", base_size = 11) +
  theme(
    axis.title = element_text(size = 11, face = "bold"),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 10),
    legend.key.size = unit(0.5, "cm"),
    panel.grid.major = element_line(color = "gray80", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    axis.ticks.length = unit(0.15, "cm"),
    axis.ticks = element_line(color = "black", linewidth = 0.3),
    plot.margin = margin(5, 5, 5, 5),
    legend.position = "right"
  )

p_scl <- ggplot(data.frame(thresholds, adjusted_scale_params, adjusted_scale_params_lower, adjusted_scale_params_upper),
                aes(x = thresholds, y = adjusted_scale_params)) +
  geom_point(size = 1.5) +
  geom_errorbar(aes(ymin = adjusted_scale_params_lower, ymax = adjusted_scale_params_upper), 
                width = 0.1, color = "blue") +
  geom_vline(xintercept = threshold_opt, color = "red", linetype = "dashed") +
  geom_hline(yintercept = adjusted_scale_params[number], color = "red", linetype = "dashed") +
  labs(x = "Threshold (log SVL)", y = "Adjusted scale parameter") +
  theme_science_polished
p_scl

save_plot(p_scl, "fig_alligator_adjusted_scale")


p_shape <- ggplot(data.frame(thresholds, shape_params, shape_params_lower, shape_params_upper),
                  aes(x = thresholds, y = shape_params)) +
  geom_point(size = 1.5) +
  geom_errorbar(aes(ymin = shape_params_lower, ymax = shape_params_upper), 
                width = 0.1, color = "blue") +
  geom_vline(xintercept = threshold_opt, color = "red", linetype = "dashed") +
  geom_hline(yintercept = shape_params[number], color = "red", linetype = "dashed") +
  labs(x = "Threshold (log SVL)", y = "Shape parameter") +
  theme_science_polished
p_shape
save_plot(p_shape, "fig_alligator_shape_param")

##############################
### ONE-SIDED HYPOTHESIS TEST FOR ξ < 0 ###
##############################
fit_test <- fevd(log_SVL, threshold = threshold_opt, type = "GP")
xi_hat <- fit_test$results$par["shape"]
sigma_hat <- fit_test$results$par["scale"]
xi_se <- summary(fit_test)$se["shape"]

# Wald z-test statistic for H0: ξ >= 0 vs H1: ξ < 0
z_value <- (xi_hat - 0) / xi_se
p_value <- pnorm(z_value)  # lower-tail test (ξ < 0)

cat("\n--- One-sided test for ξ < 0 ---\n")
cat("Estimated ξ:", round(xi_hat, 4), "\n")
cat("Standard Error:", round(xi_se, 4), "\n")
cat("Z statistic:", round(z_value, 4), "\n")
cat("One-sided p-value:", signif(p_value, 4), "\n")

##############################
# PP Plot
##############################
y_excess <- log_SVL[log_SVL > threshold_opt]
excesses <- y_excess - threshold_opt
n <- length(excesses)
# P-P plot
F_empirical <- (1:n) / n
if (abs(xi_hat) > 1e-6) {
  F_theoretical <- 1 - (1 + xi_hat * excesses / sigma_hat)^(-1/xi_hat)
} else {
  F_theoretical <- 1 - exp(-excesses / sigma_hat)
}

pp_df <- data.frame(
  Theoretical = sort(F_theoretical),
  Empirical = F_empirical
)

p8 <- ggplot(pp_df, aes(x = Theoretical, y = Empirical)) +
  geom_point(color = "darkgreen") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  labs(x = "Theoretical CDF",
       y = "Empirical CDF") +
  theme_science_polished
p8
save_plot(p8, "fig_alligator_pp_plot")


##############################
# Kolmogorov Smirnov
##############################
F_hat <- pgpd(excesses, loc = 0, scale = sigma_hat, shape = xi_hat)
ks_obs <- ks.test(F_hat, "punif")$statistic

# Parametric bootstrap
B <- 1000
ks_boot <- numeric(B)

set.seed(0)
for (b in 1:B) {
  # Generate synthetic data from the fitted GPD
  sim_data <- rgpd(n, loc = 0, scale = sigma_hat, shape = xi_hat)
  
  # Refit GPD
  fit_b <- tryCatch(gpd.fit(sim_data + threshold_opt, threshold_opt, show = FALSE), error = function(e) NULL)
  if (!is.null(fit_b)) {
    xi_b <- fit_b$mle[2]
    sigma_b <- fit_b$mle[1]
    # Transform back to uniform
    F_b <- pgpd(sim_data, loc = 0, scale = sigma_b, shape = xi_b)
    ks_boot[b] <- ks.test(F_b, "punif")$statistic
  } else {
    ks_boot[b] <- NA  # mark NA in case of fitting error
  }
}

# Remove failed bootstrap iterations
ks_boot <- ks_boot[!is.na(ks_boot)]

# Compute p-value
p_value <- mean(ks_boot >= ks_obs)

# Output
cat("Observed KS statistic:", ks_obs, "\n")
cat("Bootstrap p-value:", p_value, "\n")

# Set up grid
y_data <- log_SVL[log_SVL > threshold_opt]
n <- length(y_data)
ymax <- max(y_data)
u <- threshold_opt

# Settings
n_y_star <- 300  # Grid points for y_star
n_samples <- 5000  # Samples for xi


# Proposal: truncated normal for xi < 0
r_proposal <- function(n) rtruncnorm(n, a = -Inf, b = 0, mean = -0.3, sd = 0.2)
d_proposal <- function(x) dtruncnorm(x, a = -Inf, b = 0, mean = -0.3, sd = 0.2)

# Grids
y_star_grid <- seq(ymax, ymax + 1, length.out = n_y_star)
xi_samples <- r_proposal(n_samples)

# Initialize posterior density
posterior_y_star <- numeric(length(y_star_grid))

for (i in seq_along(y_star_grid)) {
  y_star <- y_star_grid[i]
  
  if (y_star <= ymax) next  # skip illegal values
  
  log_lik <- numeric(n_samples)
  
  for (j in seq_along(xi_samples)) {
    xi <- xi_samples[j]
    if (xi >= 0) next
    
    term1 <- n * log(y_star - u) / xi
    term2 <- -(1/xi + 1) * sum(log(y_star - y_data))
    term3 <- -n * log(-xi)
    log_lik[j] <- term1 + term2 + term3
  }
  
  # Compute log-weighted likelihood and importance weights
  weights <- exp(log_lik - log(d_proposal(xi_samples)))
  posterior_y_star[i] <- mean(weights, na.rm = TRUE)
}

posterior_y_star <- posterior_y_star / trapz(y_star_grid, posterior_y_star)
samples <- sample(y_star_grid, size = n_samples, replace = TRUE, prob = posterior_y_star)
map <- y_grid[which.max(posterior_y_star)]
hpd <- hdi(samples, credMass = ci_level)
print(exp(map))
print(exp(hpd))

plot(exp(y_star_grid), posterior_y_star, type = "l",
     xlab = "SVL endpoint (cm)", ylab = "Posterior density", col = "darkgreen", lwd = 2)
