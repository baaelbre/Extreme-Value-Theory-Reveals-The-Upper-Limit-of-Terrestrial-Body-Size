# Load libraries
library(ggplot2)
library(extRemes)
library(readxl)
library(dplyr)

####################
### POT ANALYSIS ###
####################

data <- read_excel("manual_extracted_data_MALES.xlsx")
x_log <- log(data$`Total Length`)

u_0 <- 5.5
thresholds <- seq(u_0, max(x_log) - 0.1, length.out = 30)
scale_params <- shape_params <- scale_params_lower <- scale_params_upper <- shape_params_lower <- shape_params_upper <- numeric(length(thresholds))

# Fit GPD across thresholds
for (i in seq_along(thresholds)) {
  fit <- fevd(x_log, threshold = thresholds[i], type = "GP")
  scale_params[i] <- fit$results$par["scale"]
  shape_params[i] <- fit$results$par["shape"]
  tryCatch({
    ci_vals <- ci(fit, type = "parameter")
    scale_params_lower[i] <- ci_vals[1,1]
    scale_params_upper[i] <- ci_vals[1,3]
    shape_params_lower[i] <- ci_vals[2,1]
    shape_params_upper[i] <- ci_vals[2,3]
  }, error = function(e) {
    scale_params_lower[i] <- scale_params_upper[i] <- scale_params[i]
    shape_params_lower[i] <- shape_params_upper[i] <- shape_params[i]
  })
}

# Plot threshold diagnostics
ggplot(data.frame(thresholds, scale_params, scale_params_lower, scale_params_upper),
       aes(x = thresholds, y = scale_params)) +
  geom_point() +
  geom_errorbar(aes(ymin = scale_params_lower, ymax = scale_params_upper), width = 0.1, color = "blue") +
  labs(title = "Scale Parameter vs. Threshold", x = "Threshold", y = "Scale") +
  theme_minimal()

ggplot(data.frame(thresholds, shape_params, shape_params_lower, shape_params_upper),
       aes(x = thresholds, y = shape_params)) +
  geom_point() +
  geom_errorbar(aes(ymin = shape_params_lower, ymax = shape_params_upper), width = 0.1, color = "blue") +
  labs(title = "Shape Parameter vs. Threshold", x = "Threshold", y = "Shape") +
  theme_minimal()

# Choose threshold
threshold_opt <- 5.75

################################
### PROFILE LIKELIHOOD METHOD ###
################################

profile_likelihood <- function(y_star, y, threshold) {
  neg_log_lik <- function(xi) {
    if (xi >= 0 || y_star <= max(y)) return(Inf)
    exceedances <- y[y > threshold]
    T <- length(exceedances)
    term1 <- (T / xi) * log(y_star - threshold)
    term2 <- -T * log(-xi)
    term3 <- -(1 / xi + 1) * sum(log(y_star - exceedances))
    return(-(term1 + term2 + term3))
  }
  opt <- optim(par = -0.5, fn = neg_log_lik, method = "Brent", lower = -1.5, upper = -1e-6)
  return(-opt$value)
}

# Grid and evaluation
y_star_grid <- seq(max(x_log) + 0.01, max(x_log) + 1, length.out = 100)
pll_vals <- sapply(y_star_grid, function(y_star) profile_likelihood(y_star, x_log, threshold_opt))
pll_vals_scaled <- 2 * (pll_vals - max(pll_vals))

# Extract CI and MLE
ci_cutoff <- -qchisq(0.90, df = 1)
ci_idx <- which(pll_vals_scaled >= ci_cutoff)
ci_profile <- range(y_star_grid[ci_idx])
mle_y_star <- y_star_grid[which.max(pll_vals)]

# Plot: profile likelihood (log scale)
ggplot(data.frame(y_star = y_star_grid, pll = pll_vals_scaled), aes(x = y_star, y = pll)) +
  geom_line(color = "darkred", linewidth = 1.2) +
  geom_hline(yintercept = ci_cutoff, linetype = "dashed", color = "blue") +
  geom_vline(xintercept = ci_profile, linetype = "dotted", color = "black") +
  geom_vline(xintercept = mle_y_star, linetype = "solid", color = "darkgreen") +
  labs(title = "Profile likelihood for total length endpoint (log scale)", x = "log(TL)", y = "Deviance") +
  theme_minimal()

cat("Profile 90% CI for log total length:", ci_profile, "\n")
cat("Profile 90% CI for total length (cm):", exp(ci_profile), "\n")

#################################
### JOINT LIKELIHOOD SAMPLING ###
#################################
# --- Define joint log-likelihood ---
joint_loglik <- function(y_star, xi, threshold, data) {
  if (xi >= 0 || y_star <= max(data)) return(-Inf)
  y <- data[data > threshold]
  T <- length(y)
  term1 <- (T / xi) * log(y_star - threshold)
  term2 <- -T * log(-xi)
  term3 <- -(1 / xi + 1) * sum(log(y_star - y))
  return(term1 + term2 + term3)
}

# --- Evaluate log-likelihood grid ---
y_star_vals <- seq(max(x_log) + 0.01, max(x_log) + 1, length.out = 100)
xi_vals <- seq(-0.6, -0.1, length.out = 100)
joint_grid <- outer(y_star_vals, xi_vals, Vectorize(function(y, xi) joint_loglik(y, xi, threshold_opt, x_log)))

# --- Convert to flat data frame ---
grid_df <- expand.grid(y_star = y_star_vals, xi = xi_vals)
grid_df$loglik <- as.vector(joint_grid)

# --- Find MLE and log-likelihood threshold ---
mle_index <- which.max(grid_df$loglik)
mle_y_star <- grid_df$y_star[mle_index]
mle_xi <- grid_df$xi[mle_index]
ll_max <- grid_df$loglik[mle_index]

# Compute likelihood ratio region
chi2_threshold <- qchisq(0.90, df = 2)
grid_df$inside_CI <- 2 * (ll_max - grid_df$loglik) <= chi2_threshold

# --- Compute profile likelihood ridge (argmax xi at each y_star) ---
profile_ridge <- sapply(y_star_vals, function(y_star) {
  xi_grid <- seq(-0.6, -0.05, length.out = 100)
  ll_vals <- sapply(xi_grid, function(xi) joint_loglik(y_star, xi, threshold_opt, x_log))
  xi_grid[which.max(ll_vals)]
})
profile_df <- data.frame(y_star = y_star_vals, xi = profile_ridge)

# --- Display MLE ---
cat("MLE y* (log scale):", mle_y_star, "\n")
cat("MLE y* (cm):", round(exp(mle_y_star), 2), "\n")
cat("MLE xi:", mle_xi, "\n")

# --- Plot joint log-likelihood with MLE, CI, and profile ridge ---
ggplot(grid_df, aes(x = y_star, y = xi)) +
  geom_tile(aes(fill = loglik)) +  # Only color the tiles based on loglik
  geom_contour(aes(z = as.numeric(inside_CI)), breaks = 0.5, color = "blue", linewidth = 1.2) +
  geom_line(data = profile_df, aes(x = y_star, y = xi), color = "green", linetype = "dashed", linewidth = 1) +
  geom_point(aes(x = mle_y_star, y = mle_xi), shape = 4, size = 4, color = "green", stroke = 1.5) +  # cross symbol
  annotate("text", x = mle_y_star, y = mle_xi, label = "MLE", vjust = -1, color = "green") +
  labs(title = "Joint LLH Surface with 90% CI (blue) & profile llh ridge (green)",
       x = "log(total length)", y = expression(xi), fill = "Log-Lik") +
  theme_minimal()



# Sampling from joint surface
loglik_flat <- grid_df$loglik
loglik_flat <- loglik_flat - max(loglik_flat)
weights <- exp(loglik_flat)
weights <- weights / sum(weights)

set.seed(123)
indices <- sample(seq_along(weights), size = 10000, replace = TRUE, prob = weights)
joint_samples <- grid_df[indices, ]
joint_samples$length_cm <- exp(joint_samples$y_star)

# Compute CI from joint samples
ci_joint <- quantile(joint_samples$y_star, c(0.05, 0.95))
ci_joint_cm <- exp(ci_joint)
cat("Joint 90% CI for log total length:", round(ci_joint, 3), "\n")
cat("Joint 90% CI for total length (cm):", round(ci_joint_cm, 2), "\n")

##############################
### PLOTS ON ORIGINAL SCALE ###
##############################

# Profile plot in cm
profile_df_cm <- data.frame(length_cm = exp(y_star_grid), pll = pll_vals_scaled)
ci_profile_cm <- exp(ci_profile)
mle_cm <- exp(mle_y_star)
alabama_gator_cm <- 434

ggplot(profile_df_cm, aes(x = length_cm, y = pll)) +
  geom_line(color = "darkred", linewidth = 1.2) +
  geom_hline(yintercept = ci_cutoff, linetype = "dashed", color = "blue") +
  geom_vline(xintercept = ci_profile_cm, linetype = "dotted", color = "black") +
  geom_vline(xintercept = mle_cm, color = "darkgreen") +
  geom_vline(xintercept = alabama_gator_cm, color = "purple", linetype = "dashed") +
  annotate("text", x = alabama_gator_cm + 10, y = max(profile_df_cm$pll), label = "Alabama Gator", color = "purple", hjust = 0) +
  labs(title = "Profile Likelihood (Original Scale)", x = "Total Length (cm)", y = "2 × Log-Likelihood Difference") +
  theme_minimal()

# Joint sample histogram in cm
ggplot(joint_samples, aes(x = length_cm)) +
  geom_histogram(aes(y = ..density..), bins = 50, fill = 'lightblue', color = 'black') +
  geom_vline(xintercept = ci_joint_cm, linetype = "dashed", color = "blue") +
  geom_vline(xintercept = alabama_gator_cm, linetype = "dashed", color = "purple") +
  annotate("text", x = alabama_gator_cm + 10, y = 0.02, label = "Alabama Gator", color = "purple", hjust = 0) +
  labs(title = "Joint Likelihood Sampling (Original Scale)", x = "Total Length (cm)", y = "Density") +
  theme_minimal()
