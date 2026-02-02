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

library(cmdstanr)
library(posterior)
library(bayesplot)

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


fit_svl_gpd_stan <- function(log_SVL, threshold) {
  # Compute excesses
  y_excess <- log_SVL[log_SVL > threshold]
  stan_data <- list(
    N = length(y_excess),
    y = y_excess,
    thresh = threshold,
    y_max = max(log_SVL)
  )
  
  # Compile Stan model (compile once, reuse)
  stan_model <- cmdstan_model("gpd_endpoint.stan")
  
  # Run HMC sampling
  fit <- stan_model$sample(
    data = stan_data,
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 1000,
    iter_sampling = 2000,
    adapt_delta = 0.99,  # important for endpoint models
    seed = 42
  )
  
  return(fit)
}

# Call the wrapper
fit_stan <- fit_svl_gpd_stan(log_SVL, threshold_opt)

# Extract samples
posterior_df <- fit_stan$draws(format = "data.frame")

# Convert to posterior::draws_df
posterior_draws <- as_draws_df(posterior_df)

# Summarize
posterior::summarise_draws(posterior_draws[, c("y_star", "xi")])

# Trace and density
mcmc_trace(posterior_draws, pars = c("y_star", "xi"))
mcmc_pairs(posterior_draws, pars = c("y_star", "xi"))
