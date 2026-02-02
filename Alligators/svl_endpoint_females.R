library(ggplot2)
library(extRemes)
library(ismev)
library(numDeriv)
library(dplyr)
library(MASS)
library(readxl)
library(pracma)
library(HDInterval)
# Polished theme for Science-style figures
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

##############################
### LOAD DATA ###
##############################
df <- read_excel("Data/experimental_alligator_harvest_woodward.xlsx")
# Count the number of males and females
table_sex <- table(df$Sex)
print(table_sex)

# Calculate fractions
fraction_male <- table_sex["M"] / sum(table_sex)
fraction_female <- table_sex["F"] / sum(table_sex)

cat("\n--- Sex distribution ---\n")
cat("Fraction male:", round(fraction_male, 3), "\n")
cat("Fraction female:", round(fraction_female, 3), "\n")

df_females <- df %>% filter(Sex == "F", !is.na(SVL), !is.na(TL))

# Log-transform the measurements
log_SVL <- log(df_females$SVL)
log_TL <- log(df_females$TL)


##############################
### POT ANALYSIS ON LOG SVL ###
##############################

u_0 <- 4.9
thresholds <- seq(u_0-0.25, u_0+0.04, length.out = 50)

scale_params <- shape_params <- shape_params_lower <- shape_params_upper <- numeric(length(thresholds))
scale_params_lower <- scale_params_upper <- adjusted_scale_params <- adjusted_scale_params_lower <- adjusted_scale_params_upper <- numeric(length(thresholds))

for (i in seq_along(thresholds)) {
  fit <- fevd(log_SVL, threshold = thresholds[i], type = "GP")
  scale_params[i] <- fit$results$par["scale"]
  shape_params[i] <- fit$results$par["shape"]
  adjusted_scale_params[i] <- scale_params[i] - shape_params[i] * (thresholds[i] - u_0)
  
  tryCatch({
    ci <- ci(fit, type = 'parameter', alpha = 0.1)
    cov <- summary(fit)$cov.theta
    shape_params_lower[i] <- ci[2,1]
    shape_params_upper[i] <- ci[2,3]
    scale_params_lower[i] <- ci[1,1]
    scale_params_upper[i] <- ci[1,3]
    var_adj <- cov[1,1] + cov[2,2] * (thresholds[i] - u_0)^2
    adjusted_scale_params_lower[i] <- adjusted_scale_params[i] - qnorm(0.95) * sqrt(var_adj)
    adjusted_scale_params_upper[i] <- adjusted_scale_params[i] + qnorm(0.95) * sqrt(var_adj)
  }, error = function(e) {})
}

number <- 40
threshold_opt <- thresholds[number]
#threshold_opt <- log(139)

##############################
### PLOT PARAMETER STABILITY ###
##############################

# Adjusted Scale Parameter vs. Threshold
p_param_adj_scale <- ggplot(data.frame(thresholds, adjusted_scale_params, adjusted_scale_params_lower, adjusted_scale_params_upper),
       aes(x = thresholds, y = adjusted_scale_params)) +
  geom_point() +
  geom_errorbar(aes(ymin = adjusted_scale_params_lower, ymax = adjusted_scale_params_upper), 
                width = 0.1, color = "blue") +
  geom_vline(xintercept = threshold_opt, color = "red", linetype = "dashed") +
  geom_hline(yintercept = adjusted_scale_params[number], color = "red", linetype = "dashed") +
  labs(
    x = "Threshold (log SVL)", 
    y = "Adjusted scale parameter") +
  theme_science_polished
p_param_adj_scale

ggsave("Figures_alligators/fig_param_adj_scale_female.tiff", plot = p_param_adj_scale, dpi = 600,
       width = 7, height = 5, units = "in", compression = "lzw")
ggsave("Figures_alligators/fig_param_adj_scale_female.png", plot = p_param_adj_scale, dpi = 600,
       width = 7, height = 5, units = "in")

# Shape Parameter vs. Threshold
p_param_shape <- ggplot(data.frame(thresholds, shape_params, shape_params_lower, shape_params_upper),
       aes(x = thresholds, y = shape_params)) +
  geom_point() +
  geom_errorbar(aes(ymin = shape_params_lower, ymax = shape_params_upper), 
                width = 0.1, color = "blue") +
  geom_vline(xintercept = threshold_opt, color = "red", linetype = "dashed") +
  geom_hline(yintercept = shape_params[number], color = "red", linetype = "dashed") +
  labs(
    x = "Threshold (log SVL)", 
    y = "Shape parameter") +
  theme_science_polished
p_param_shape

ggsave("Figures_alligators/fig_param_shape_female.tiff", plot = p_param_adj_scale, dpi = 600,
       width = 7, height = 5, units = "in", compression = "lzw")
ggsave("Figures_alligators/fig_param_shape_female.png", plot = p_param_adj_scale, dpi = 600,
       width = 7, height = 5, units = "in")

fit_test <- fevd(log_SVL, threshold = threshold_opt, type = "GP")
xi_hat <- fit_test$results$par["shape"]
xi_se <- summary(fit_test)$se["shape"]
sigma_hat <- fit_test$results$par["scale"]

# Wald z-test statistic for H0: ξ >= 0 vs H1: ξ < 0
z_value <- (xi_hat - 0) / xi_se
p_value <- pnorm(z_value)  # lower-tail test (ξ < 0)

cat("\n--- One-sided test for ξ < 0 ---\n")
cat("Estimated ξ:", round(xi_hat, 4), "\n")
cat("Standard Error:", round(xi_se, 4), "\n")
cat("Z statistic:", round(z_value, 4), "\n")
cat("One-sided p-value:", signif(p_value, 4), "\n")

##############################
# Kolmogorov Smirnov
##############################
F_hat <- pgpd(excesses, loc = 0, scale = sigma_hat, shape = xi_hat)

# Step 2: Visual check – histogram with reference uniform line
p9 <- ggplot(data.frame(F_hat = F_hat), aes(x = F_hat)) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = "skyblue", color = "black", alpha = 0.7) +
  geom_hline(yintercept = 1, color = "red", linetype = "dashed") +
  geom_density(color = "darkblue", size = 1.2, adjust = 1.5) +
  labs(title = "Uniformity Check of Transformed Exceedances",
       x = expression(hat(F)(y)),
       y = "Density") +
  theme_science_polished
p9

ggsave("Figures_alligators/fig_KS_female.tiff", plot = p9, dpi = 600,
       width = 7, height = 5, units = "in", compression = "lzw")
ggsave("Figures_alligators/fig_KS_female.png", plot = p9, dpi = 600,
       width = 7, height = 5, units = "in")

ks_obs <- ks.test(F_hat, "punif")$statistic

# Parametric bootstrap
B <- 10000
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

##############################################
### Q-Q and P-P PLOTS FOR FEMALE ALLIGATORS ###
##############################################

# Fit GPD to female log SVL above threshold
gpd_fit_female <- fevd(log_SVL, threshold = threshold_opt, type = "GP")
scale_hat_f <- gpd_fit_female$results$par["scale"]
shape_hat_f <- gpd_fit_female$results$par["shape"]

# Exceedances
y_excess_f <- log_SVL[log_SVL > threshold_opt]
excesses_f <- y_excess_f - threshold_opt
n_f <- length(excesses_f)

# Probabilities and theoretical quantiles for Q-Q plot
probs_f <- ppoints(n_f)

if (abs(shape_hat_f) > 1e-6) {
  theo_q_f <- threshold_opt + scale_hat_f / shape_hat_f * ((probs_f^(-shape_hat_f)) - 1)
} else {
  theo_q_f <- threshold_opt - scale_hat_f * log(probs_f)
}

qq_df_f <- data.frame(
  Theoretical = rev(theo_q_f),
  Empirical = sort(y_excess_f + threshold_opt)
)

p_qq_female <- ggplot(qq_df_f, aes(x = Theoretical, y = Empirical)) +
  geom_point(color = "steelblue") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  labs(title = "Q-Q Plot for GPD Fit (Females)",
       x = "Theoretical Quantiles",
       y = "Empirical Quantiles") +
  theme_science_polished
p_qq_female
save_plot(p_qq_female, "fig_diagnostics_qq_alligator_female")


# Theoretical and empirical CDF values for P-P plot
F_empirical_f <- (1:n_f) / n_f
if (abs(shape_hat_f) > 1e-6) {
  F_theoretical_f <- 1 - (1 + shape_hat_f * excesses_f / scale_hat_f)^(-1 / shape_hat_f)
} else {
  F_theoretical_f <- 1 - exp(-excesses_f / scale_hat_f)
}

pp_df_f <- data.frame(
  Theoretical = sort(F_theoretical_f),
  Empirical = F_empirical_f
)

p_pp_female <- ggplot(pp_df_f, aes(x = Theoretical, y = Empirical)) +
  geom_point(color = "darkgreen") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  labs(
       x = "Theoretical CDF",
       y = "Empirical CDF") +
  theme_science_polished
p_pp_female
save_plot(p_pp_female, "fig_diagnostics_pp_alligator_female")


##############################
### POSTERIOR SVL ###
##############################
set.seed(0)
run_mh_sampler_female <- function(
    log_data,
    threshold,
    n_iter = 50000,
    init = c(y_star = max(log_data) + 0.05, xi = -0.2),
    proposal_sd = c(0.01, 0.02),
    burn_in = 10000,
    thin = 20
) {
  y_data <- log_data[log_data > threshold]
  n <- length(y_data)
  ymax <- max(y_data)
  
  # Log-likelihood in (y*, xi) parametrization
  loglik_joint <- function(y_star, xi) {
    if (xi >= 0 || y_star <= ymax) return(-Inf)
    term1 <- n * log(y_star - threshold) / xi
    term2 <- -(1 / xi + 1) * sum(log(y_star - y_data))
    term3 <- -n * log(-xi)
    return(term1 + term2 + term3)
  }
  
  samples <- matrix(NA, nrow = n_iter, ncol = 2)
  colnames(samples) <- c("y_star", "xi")
  samples[1, ] <- init
  acc <- 0
  
  for (i in 2:n_iter) {
    current <- samples[i - 1, ]
    proposal <- rnorm(2, mean = current, sd = proposal_sd)
    y_star_prop <- proposal[1]
    xi_prop <- proposal[2]
    
    log_post_prop <- loglik_joint(y_star_prop, xi_prop)
    log_post_curr <- loglik_joint(current[1], current[2])
    
    log_alpha <- log_post_prop - log_post_curr
    if (log(runif(1)) < log_alpha) {
      samples[i, ] <- proposal
      acc <- acc + 1
    } else {
      samples[i, ] <- current
    }
  }
  
  cat("Acceptance rate:", round(acc / n_iter, 3), "\n")
  
  samples_mcmc <- as.data.frame(samples)
  samples_mcmc <- samples_mcmc[(burn_in + 1):n_iter, ]
  samples_mcmc <- samples_mcmc[seq(1, nrow(samples_mcmc), by = thin), ]
  return(samples_mcmc)
}
mh_samples_female <- run_mh_sampler_female(log_SVL, threshold_opt)
# Use samples from the Metropolis-Hastings sampler
y_star_post <- mh_samples_female$y_star

# Posterior density estimate (KDE)
posterior_kde <- density(y_star_post)

# MAP estimate
map_bayes <- posterior_kde$x[which.max(posterior_kde$y)]

# HPDI using 95% credible mass
hpd_bayes <- hdi(y_star_post, credMass = 0.95)

# Optional: equal-tailed CI for comparison
ci_bayes <- quantile(y_star_post, probs = c(0.025, 0.975))

# Exponentiate for original scale (SVL)
map_bayes_exp <- exp(map_bayes)
hpd_bayes_exp <- exp(hpd_bayes)
ci_bayes_exp  <- exp(ci_bayes)

# Print summary
cat("\n--- Posterior from MH samples ---\n")
cat("MAP (log SVL):", round(map_bayes, 4), "\n")
cat("95% HPDI:", round(hpd_bayes[1], 4), "-", round(hpd_bayes[2], 4), "\n")
cat("Exponentiated MAP:", round(map_bayes_exp, 2), "cm\n")
cat("95% HPDI (cm):", round(hpd_bayes_exp[1], 2), "-", round(hpd_bayes_exp[2], 2), "cm\n")

##############################
### PLOT POSTERIOR IN ORIGINAL SVL SCALE ###
##############################
posterior_kde_svl_exp <- density(exp(mh_samples_female$y_star))
map_bayes_exp_female <- posterior_kde_svl_exp$x[which.max(posterior_kde_svl_exp$y)]
hpd_bayes_exp_female <- hdi(exp(mh_samples_female$y_star), credMass = 0.95)

# Constants
svl_record_female <- 185  # Record female alligator SVL in cm
max_SVL_empirical_female <- max(df_females$SVL)

# Plot
p_svl_female <- ggplot(data.frame(samples = exp(mh_samples_female$y_star)), aes(x = samples)) +
  geom_histogram(aes(y = ..density..),
                 fill = "lightgreen", color = "black", bins = 40, alpha = 0.4) +
  geom_line(data = data.frame(x = posterior_kde_svl_exp$x, y = posterior_kde_svl_exp$y),
            aes(x = x, y = y), color = "darkgreen", size = 1.2) +
  geom_vline(xintercept = map_bayes_exp_female, color = "purple", linetype = "dashed", linewidth = 0.8) +
  geom_vline(xintercept = hpd_bayes_exp_female, color = "orange", linetype = "dotdash", linewidth = 0.8) +
  geom_vline(xintercept = max_SVL_empirical_female, color = "gray", linetype = "dotted") +
  geom_vline(xintercept = svl_record_female, color = "black", linetype = "solid", linewidth = 1.0) +
  annotate("text", x = svl_record_female - 2, y = max(posterior_kde_svl_exp$y) * 0.5,
           label = "Record female (185 cm)", color = "black", angle = 90, vjust = -0.5) +
  labs(
    x = "SVL endpoint (cm)",
    y = "Posterior Density"
  ) +
  theme_minimal()

# Optional: display
p_svl_female

# Optional: save high-resolution TIFF and PNG
ggsave("Figures_alligators/fig_female_svl_endpoint_MH.tiff", plot = p_svl_female, dpi = 600,
       width = 7, height = 5, units = "in", compression = "lzw")
ggsave("Figures_alligators/fig_female_svl_endpoint_MH.png", plot = p_svl_female, dpi = 600,
       width = 7, height = 5, units = "in")

##############################
### log(TL) ~ log(SVL) REGRESSION ###
##############################

df_clean <- df %>%
  filter(Deform == 0, !is.na(SVL), !is.na(TL)) %>%
  mutate(log_SVL = log(SVL), log_TL = log(TL))

fit_log <- lm(log_TL ~ log_SVL, data = df_clean)
summary(fit_log)

##############################
### PROPAGATE UNCERTAINTY: log(SVL) posterior -> log(TL) ###
##############################

# Regression coefficients and their standard errors
coef_log <- coef(fit_log)
intercept_est <- coef_log[1]
slope_est <- coef_log[2]

se <- summary(fit_log)$coefficients[,2]
intercept_se <- se[1]
slope_se <- se[2]

# Sample regression parameters
n_samples <- length(y_star_post)
sampled_intercepts <- rnorm(n_samples, mean = intercept_est, sd = intercept_se)
sampled_slopes <- rnorm(n_samples, mean = slope_est, sd = slope_se)

# Estimate residual standard deviation from the regression
residual_sd <- summary(fit_log)$sigma

# Add residual variability to log(TL)
log_TL_samples <- sampled_intercepts + sampled_slopes * y_star_post +
  rnorm(n_samples, mean = 0, sd = residual_sd)

# Transform to original scale
TL_samples <- exp(log_TL_samples)

TL_median <- median(TL_samples)
TL_CI <- quantile(TL_samples, probs = c(0.255, 0.975))
TL_HPDI <- hdi(TL_samples, credMass = 0.95)
posterior_kde_tl <- density(TL_samples)
TL_MAP <- posterior_kde_tl$x[which.max(posterior_kde_tl$y)]  # MAP estimate

# Reprint including MAP
cat("\n--- TL endpoint propagated uncertainty ---\n")
cat("MAP TL:", round(TL_MAP, 2), "cm\n")
cat("95% Equal-tailed CI:", round(TL_CI[1], 2), "-", round(TL_CI[2], 2), "cm\n")
cat("95% HPDI:", round(TL_HPDI[1], 2), "-", round(TL_HPDI[2], 2), "cm\n")


cat("\n--- TL endpoint propagated uncertainty ---\n")
cat("Median TL:", round(TL_median, 2), "cm\n")
cat("95% Equal-tailed CI:", round(TL_CI[1], 2), "-", round(TL_CI[2], 2), "cm\n")
cat("95% HPDI:", round(TL_HPDI[1], 2), "-", round(TL_HPDI[2], 2), "cm\n")

##############################
### EXTRAPOLATE MAX log(TL) ###
##############################

log_svl_max <- map_bayes
log_TL_estimated <- coef_log[2] * log_svl_max + coef_log[1]
TL_max_estimated <- exp(log_TL_estimated)

max_TL_empirical <- max(df_females$TL)
max_SVL_empirical <- max(df_females$SVL)

cat("\n--- Maximum TL extrapolated from EVT SVL ---\n")
cat("MAP SVL estimate:", round(exp(log_svl_max), 2), "cm\n")
cat("Estimated TL:", round(TL_max_estimated, 2), "cm\n")
cat("Empirical max TL:", max_TL_empirical, "cm\n")

##############################
### PLOT POSTERIOR FOR TL ENDPOINT ###
##############################

p_tl_female <- ggplot(data.frame(TL = TL_samples), aes(x = TL)) +
  geom_histogram(aes(y = ..density..), bins = 40, fill = "lightblue", color = "black", alpha = 0.6) +
  geom_line(data = data.frame(x = posterior_kde_tl$x, y = posterior_kde_tl$y),
            aes(x = x, y = y), color = "darkblue", size = 1.2) +
  geom_vline(xintercept = TL_MAP, color = "purple", linetype = "dashed", linewidth = 1.2) +
  geom_vline(xintercept = TL_HPDI, color = "orange", linetype = "dotdash", linewidth = 1.2) +
  geom_vline(xintercept = max_TL_empirical, color = "gray", linetype = "dashed", linewidth = .8) +
  geom_vline(xintercept = 322, color = "black", linetype = "solid", linewidth = 1.2) +
  annotate("text", x = 322 - 8, y = max(posterior_kde_tl$y) * 0.4,
           label = "Record female (322 cm)", color = "black", angle = 90, vjust = -0.5) +
  labs(
    x = "TL endpoint (cm)",
    y = "Posterior Density"
  ) +
  theme_science_polished

# Show or save
p_tl_female
ggsave("Figures_alligators/fig_female_tl_endpoint_MH.tiff", plot = p_tl_female, dpi = 600,
       width = 7, height = 5, units = "in", compression = "lzw")
ggsave("Figures_alligators/fig_female_tl_endpoint_MH.png", plot = p_tl_female, dpi = 600,
       width = 7, height = 5, units = "in")

##############################
### PLOT log(TL) ~ log(SVL) ###
##############################

log_SVL_seq <- seq(min(df_clean$log_SVL), max(df_clean$log_SVL) + 0.3, length.out = 100)
predicted_log_TL <- predict(fit_log, newdata = data.frame(log_SVL = log_SVL_seq))

line_data <- data.frame(log_SVL = log_SVL_seq, log_TL = predicted_log_TL)
line_data$linetype <- ifelse(log_SVL_seq <= max(df_clean$log_SVL), "solid", "dashed")

ggplot(df_females %>% mutate(log_SVL = log(SVL), log_TL = log(TL)), aes(x = log_SVL, y = log_TL)) +
  geom_point(alpha = 0.5, color = "gray") +
  geom_line(data = subset(line_data, linetype == "solid"), aes(x = log_SVL, y = log_TL), color = "blue") +
  geom_line(data = subset(line_data, linetype == "dashed"), aes(x = log_SVL, y = log_TL), color = "blue", linetype = "dashed") +
  geom_errorbar(aes(x = map_bayes, ymin = log(TL_HPDI[1]), ymax = log(TL_HPDI[2])), width = 0.02, color = "darkgreen") +
  geom_errorbarh(aes(y = log_TL_estimated, xmin = hpd_bayes[1], xmax = hpd_bayes[2]), height = 0.02, color = "darkgreen") +
  geom_point(aes(x = map_bayes, y = log_TL_estimated), color = "darkgreen", size = 3) +
  geom_point(aes(x = log(max_SVL_empirical), y = log(max_TL_empirical)), color = "red", size = 3) +
  geom_point(aes(x = log(185.1), y = log(322)), color = "black", size = 3) +
  geom_vline(xintercept = threshold_opt, linetype = "dashed", color = "black") +
  annotate("text", x = map_bayes + 0.09, y = log_TL_estimated - 0.06, label = "EVT extrapolation", color = "darkgreen") +
  annotate("text", x = log(max_SVL_empirical) - 0.1, y = log(max_TL_empirical) + 0.03, label = "Max observed", color = "red") +
  annotate("text", x = log(185.1) - 0.05, y = log(322) + 0.04, label = "Record female", color = "black") +
  annotate("text", x = threshold_opt, y = mean(df_clean$log_TL), label = "threshold", color = "black") +
  labs(
    title = "log(TL) vs log(SVL) with EVT extrapolation",
    x = "log Snout-Vent Length",
    y = "log Total Length"
  ) + coord_fixed()




