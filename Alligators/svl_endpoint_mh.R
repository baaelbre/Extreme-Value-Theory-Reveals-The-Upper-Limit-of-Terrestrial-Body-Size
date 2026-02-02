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

##############################
### POSTERIOR SVL
##############################
set.seed(0)
run_mh_sampler_svl <- function(
    log_data,
    threshold,
    n_iter = 50000,
    init = c(y_star = max(log_data) + 0.1, psi = log(0.3)),  # xi = -exp(psi) = -0.3
    proposal_sd = c(0.02, 0.02),
    burn_frac = 0.1,
    thin = 1
) {
  y_data <- log_data[log_data > threshold]
  n <- length(y_data)
  ymax <- max(y_data)
  burn_in <- floor(burn_frac * n_iter)
  
  loglik_joint <- function(y_star, psi) {
    xi <- -exp(psi)  # ensures xi < 0
    if (y_star <= ymax) return(-Inf)
    term1 <- n * log(y_star - threshold) / xi
    term2 <- -(1 / xi + 1) * sum(log(y_star - y_data))
    term3 <- -n * log(-xi)
    return(term1 + term2 + term3)
  }
  
  samples <- matrix(NA, nrow = n_iter, ncol = 2)
  colnames(samples) <- c("y_star", "psi")
  samples[1, ] <- init
  acc <- 0
  
  for (i in 2:n_iter) {
    current <- samples[i - 1, ]
    proposal <- rnorm(2, mean = current, sd = proposal_sd)
    log_post_curr <- loglik_joint(current[1], current[2])
    log_post_prop <- loglik_joint(proposal[1], proposal[2])
    
    log_alpha <- log_post_prop - log_post_curr
    if (log(runif(1)) < log_alpha) {
      samples[i, ] <- proposal
      acc <- acc + 1
    } else {
      samples[i, ] <- current
    }
  }
  
  cat("Acceptance rate:", round(acc / n_iter, 3), "\n")
  
  samples_df <- as.data.frame(samples)
  samples_df <- samples_df[(burn_in + 1):n_iter, ]
  samples_df <- samples_df[seq(1, nrow(samples_df), by = thin), ]
  
  # Transform back to xi for posterior summaries
  samples_df$xi <- -exp(samples_df$psi)
  
  return(samples_df)
}


mh_samples <- run_mh_sampler_svl(
  log_SVL, threshold_opt,
  n_iter = 1e6,
  proposal_sd = c(0.15, 0.01),  # you can tune this
  burn_frac = 0.1,
  thin = 1
)

mcmc_reparam <- as.mcmc(mh_samples[, c("y_star", "xi")])
effectiveSize(mcmc_reparam)
acf(mcmc_reparam[, "y_star"])


# Trace plots
plot(mcmc_samples)

# Or manually for each parameter:
par(mfrow = c(2,2))
plot(mcmc_samples[, "y_star"], type = "l", main = "Trace plot: y_star")
plot(density(mcmc_samples[, "y_star"]), main = "Density: y_star")
par(mfrow = c(1,1))

# Extract posterior samples and summaries
y_star_post_svl <- mh_samples_svl$y_star
posterior_kde_svl <- density(y_star_post_svl)
map_bayes_svl <- posterior_kde_svl$x[which.max(posterior_kde_svl$y)]
ci_bayes_svl <- quantile(y_star_post_svl, probs = c(0.025, 0.975))
hpd_bayes_svl <- hdi(y_star_post_svl, credMass = 0.95)

# Convert to SVL (cm) scale
map_bayes_exp <- exp(map_bayes_svl)
ci_bayes_exp <- exp(ci_bayes_svl)
hpd_bayes_exp <- exp(hpd_bayes_svl)

cat("\n--- MH Posterior Inference for log(SVL) ---\n")
cat("MAP (log SVL):", round(map_bayes_svl, 4), "\n")
cat("95% CI:", round(ci_bayes_svl[1], 4), "-", round(ci_bayes_svl[2], 4), "\n")
cat("95% HPDI:", round(hpd_bayes_svl[1], 4), "-", round(hpd_bayes_svl[2], 4), "\n")

cat("\n--- SVL Posterior Inference (original scale, cm) ---\n")
cat("MAP:", round(map_bayes_exp, 2), "cm\n")
cat("95% CI:", round(ci_bayes_exp[1], 2), "-", round(ci_bayes_exp[2], 2), "cm\n")
cat("95% HPDI:", round(hpd_bayes_exp[1], 2), "-", round(hpd_bayes_exp[2], 2), "cm\n")

##############################
### PLOT POSTERIOR IN ORIGINAL SVL SCALE ###
##############################

max_SVL_empirical <- exp(max(log_SVL))
svl_stokes <- 236  # cm

posterior_kde_svl_exp <- density(exp(y_star_post_svl))
map_bayes_exp <- posterior_kde_svl_exp$x[which.max(posterior_kde_svl_exp$y)]

p_svl <- ggplot(data.frame(samples = exp(y_star_post_svl)), aes(x = samples)) +
  geom_histogram(aes(y = ..density..),
                 fill = "lightgreen", color = "black", bins = 40, alpha = 0.4) +
  geom_line(data = data.frame(x = posterior_kde_svl_exp$x, y = posterior_kde_svl_exp$y), 
            aes(x = x, y = y), color = "darkgreen", size = 1.2) +
  geom_vline(xintercept = map_bayes_exp, color = "purple", linetype = "dashed", linewidth = 0.8) +
  geom_vline(xintercept = hpd_bayes_exp, color = "orange", linetype = "dotdash", linewidth = 0.8) +
  geom_vline(xintercept = max_SVL_empirical, color = "gray", linetype = "dotted") +
  geom_vline(xintercept = svl_stokes, color = "black", linetype = "solid", linewidth = 1.0) +
  #annotate("text", x = svl_stokes - 2, y = max(posterior_kde_svl_exp$y) * 0.5, 
  #         label = "Stokes alligator (236 cm)", color = "black", angle = 90, vjust = -0.5) +
  labs(x = "SVL endpoint (cm)", y = "Posterior density") +
  theme_science_polished
p_svl

save_plot(p_svl, "fig_alligator_svl_posterior_MH")



##############################
### log(TL) ~ log(SVL) REGRESSION ###
##############################
# Filter only individuals with no deformities and valid SVL + TL
df_clean <- df %>%
  filter(Deform == 0, !is.na(SVL), !is.na(TL)) %>%
  mutate(log_SVL = log(SVL), log_TL = log(TL))
fit_log <- lm(log_TL ~ log_SVL, data = df_clean)
summary(fit_log)
#plot(fit_log)

##############################
### PROPAGATE MH log(SVL) POSTERIOR -> log(TL) ###
##############################

# Regression coefficients and their standard errors
coef_log <- coef(fit_log)
intercept_est <- coef_log[1]
slope_est <- coef_log[2]
se <- summary(fit_log)$coefficients[,2]
intercept_se <- se[1]
slope_se <- se[2]

# Draw regression coefficients with uncertainty
n_samples <- nrow(mh_samples_svl)
sampled_intercepts <- rnorm(n_samples, mean = intercept_est, sd = intercept_se)
sampled_slopes     <- rnorm(n_samples, mean = slope_est, sd = slope_se)

# Use log(SVL) samples from MH posterior
sampled_log_SVL <- mh_samples_svl$y_star

# Residual SD from model
residual_sd <- summary(fit_log)$sigma

# Generate log(TL) samples with residual variation
sampled_log_TL <- sampled_intercepts + sampled_slopes * sampled_log_SVL +
  rnorm(n_samples, mean = 0, sd = residual_sd)

# Back-transform to TL
sampled_TL <- exp(sampled_log_TL)

# Summary statistics
# Compute MAP estimate from posterior density
posterior_kde_tl <- density(sampled_TL)
map_TL <- posterior_kde_tl$x[which.max(posterior_kde_tl$y)]
map_log_TL <- log(map_TL)
TL_CI     <- quantile(sampled_TL, probs = c(0.025, 0.975))
TL_HPDI   <- hdi(sampled_TL, credMass = 0.95)

cat("\n--- TL endpoint (from MH-based SVL) ---\n")
cat("MAP TL:", round(map_TL, 2), "cm\n")
cat("95% Equal-tailed PI:", round(TL_CI[1], 2), "-", round(TL_CI[2], 2), "cm\n")
cat("95% HPDI:", round(TL_HPDI[1], 2), "-", round(TL_HPDI[2], 2), "cm\n")

##############################
### PLOT POSTERIOR FOR TL ENDPOINT ###
##############################
posterior_kde_tl <- density(sampled_TL)

p_tl <- ggplot(data.frame(TL = sampled_TL), aes(x = TL)) +
  geom_histogram(aes(y = ..density..), bins = 40, fill = "lightblue", color = "black", alpha = 0.6) +
  geom_line(data = data.frame(x = posterior_kde_tl$x, y = posterior_kde_tl$y), 
            aes(x = x, y = y), color = "darkblue", linewidth = 1.2) +
  geom_vline(xintercept = map_TL, color = "purple", linetype = "dashed", linewidth = 0.8) +
  geom_vline(xintercept = TL_HPDI, color = "orange", linetype = "dotdash", linewidth = 0.8) +
  geom_vline(xintercept = max(df_males$TL), color = "gray", linetype = "dashed") +
  geom_vline(xintercept = 450, color = "black", linetype = "solid", linewidth = 1.0) +
  labs(x = "TL endpoint (cm)", y = "Posterior density") +
  theme_science_polished
p_tl
# Save
save_plot(p_tl, "fig_alligator_tl_posterior_MH_density")

p_tl
save_plot(p_tl, "fig_alligator_tl_posterior_MH")

##############################
### PLOT log(TL) ~ log(SVL) with sex-based coloring ###
##############################

# Generate regression prediction line with extension
log_SVL_seq <- seq(min(df_clean$log_SVL), max(df_clean$log_SVL) + 0.3, length.out = 100)
predicted_log_TL <- predict(fit_log, newdata = data.frame(log_SVL = log_SVL_seq))

line_data <- data.frame(log_SVL = log_SVL_seq, log_TL = predicted_log_TL)
line_data$linetype <- ifelse(log_SVL_seq <= max(df_clean$log_SVL), "solid", "dashed")
p6 <- ggplot(df_clean, aes(x = log_SVL, y = log_TL, color = Sex)) +
  geom_point(alpha = 0.5) +
  geom_line(data = subset(line_data, linetype == "solid"), 
            aes(x = log_SVL, y = log_TL), inherit.aes = FALSE, color = "blue") +
  geom_line(data = subset(line_data, linetype == "dashed"), 
            aes(x = log_SVL, y = log_TL), inherit.aes = FALSE, color = "blue", linetype = "dashed") +
  geom_errorbar(aes(x = map_bayes_svl, ymin = log(TL_HPDI[1]), ymax = log(TL_HPDI[2])), 
                width = 0.02, color = "darkgreen") +
  geom_errorbarh(aes(y = map_log_TL, xmin = hpd_bayes_svl[1], xmax = hpd_bayes_svl[2]), 
                 height = 0.02, color = "darkgreen") +
  geom_point(aes(x = map_bayes_svl, y = map_log_TL), 
             inherit.aes = FALSE, color = "darkgreen", size = 2.5) +
  geom_point(aes(x = log(236), y = log(450)), 
             inherit.aes = FALSE, color = "gray", size = 2.5) +
  geom_vline(xintercept = threshold_opt, linetype = "dashed", color = "black") +
  #annotate("text", x = map_bayes_svl + 0.09, y = map_log_TL - 0.06, label = "EVT extrapolation", color = "darkgreen", size = 3) +
  #annotate("text", x = log(max(df_males$SVL)) - 0.1, y = log(max(df_males$TL)) + 0.03, label = "Max observed", color = "red", size = 3) +
  #annotate("text", x = log(236) - 0.05, y = log(450) + 0.04, label = "Stokes alligator", color = "black", size = 3) +
  annotate("text", x = threshold_opt + 0.01, y = mean(df_clean$log_TL), label = 'threshold', color = "black", size = 3) +
  scale_color_manual(values = c("F" = "deeppink", "M" = "dodgerblue")) +
  labs(color = "Sex", x = "log Snout-Vent Length", y = "log Total Length") +
  coord_fixed() +
  theme_science_polished
p6
save_plot(p6, "fig_alligator_regression_panel", width = 12, height = 7)


##############################
### EXCEEDANCE PROBABILITY: SVL > STOKES
##############################

# Define SVL of Stokes alligator in log scale
log_SVL_stokes <- log(236)

# Number of exceedances and total observations
n_exceed <- sum(log_SVL > threshold_opt)
n_total  <- length(log_SVL)

# Empirical exceedance rate
P1_empirical <- n_exceed / n_total

# Compute conditional exceedance probabilities using posterior samples
excess_stokes <- exp(log_SVL_stokes) - exp(threshold_opt)
y_star_post_exp <- exp(mh_samples_svl$y_star)
xi_post         <- mh_samples_svl$xi

P2_post <- ((y_star_post_exp - exp(log_SVL_stokes)) / 
              (y_star_post_exp - exp(threshold_opt)))^(-1 / xi_post)

# Full exceedance probability
P_exceed_samples <- P1_empirical * P2_post
P_exceed_stokes <- mean(P_exceed_samples, na.rm=TRUE)

# 95% HPDI for the exceedance probability
HPDI_exceed <- hdi(P_exceed_samples, credMass = 0.95)

cat("\n--- Exceedance Probability: SVL > 236 cm ---\n")
cat("Point estimate (in percent):", signif(P_exceed_stokes, 4)*100, "\n")
cat("95% HPDI (in percent): [", signif(HPDI_exceed[1], 4)*100, ",", signif(HPDI_exceed[2], 4)*100, "]\n")



