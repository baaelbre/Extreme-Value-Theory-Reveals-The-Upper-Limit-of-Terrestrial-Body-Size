library(ggplot2)
library(extRemes)
library(ismev)
library(numDeriv)
library(dplyr)
library(MASS)
library(readxl)
library(pracma)  # For trapz
library(coda)    # For HPDI
library(HDInterval)

##############################
### DATA PREPARATION ###
##############################

data_farlow <-  c(
  440, 610, 620, 690, 720, 766, 805, 817, 840, 865, 865, 910, 912, 965, 970, 1000, 
  1040, 1120, 1130, 1295, 1420, 1422, 1445, 1560, 1613, 1689, 1765, 1790, 1791, 1854,
  1854, 1913, 1940, 1956, 1981, 2032, 2032, 2032, 2045, 2051, 2057, 2057, 2083, 2134,
  2159, 2159, 2184, 2215, 2223, 2261, 2237, 2340, 2349, 2426, 2435, 2438, 2464, 2470,
  2515, 2565, 2565, 2572, 2584, 2591, 2516, 2657, 2667, 2685, 2686, 2699, 2700, 2705,
  2743, 2743, 2769, 2769, 2769, 2769, 2769, 2794, 2794, 2807, 2819, 2845, 2845, 2845,
  2857, 2857, 2857, 2870, 2921, 2946, 2962, 3073, 3302, 3518, 3581, 3607, 3607, 3632,
  3670, 3700, 3785, 3810, 3810, 3823, 3873, 4039
)

males_farlow <- c(
  610, 805, 840, 865, 865, 965, 1000, 1120, 1295, 1422, 1613, 1765, 1790, 1791, 1913,
  1956, 2032, 2032, 2032, 2083, 2134, 2584, 2591, 3073, 3302, 3518, 3581, 3607, 3607,
  3632, 3670, 3700, 3785, 3810, 3810, 3823, 3873, 4039
)/10

data <- read_excel("manual_extracted_data_MALES.xlsx")
x_log <- sort(c(log(data$`Total Length`), log(males_farlow)))
x_log <- sort(c(log(data$`Total Length`)))

##############################
### THRESHOLD SELECTION ###
##############################

u_0 <- 5.67
thresholds <- seq(u_0, 5.95, length.out = 50)

scale_params <- shape_params <- shape_params_lower <- shape_params_upper <- numeric(length(thresholds))
scale_params_lower <- scale_params_upper <- adjusted_scale_params <- adjusted_scale_params_lower <- adjusted_scale_params_upper <- numeric(length(thresholds))

for (i in seq_along(thresholds)) {
  fit <- fevd(x_log, threshold = thresholds[i], type = "GP")
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

number <- 22
threshold_opt <- thresholds[number]

# Stability plots
ggplot(data.frame(thresholds, adjusted_scale_params), aes(x = thresholds, y = adjusted_scale_params)) +
  geom_point() +
  geom_errorbar(aes(ymin = adjusted_scale_params_lower, ymax = adjusted_scale_params_upper), width = 0.1, color = "blue") +
  geom_vline(xintercept = threshold_opt, color = "red", linetype = "dashed") +
  labs(title = "Adjusted Scale Parameter vs. Threshold", x = "Threshold", y = "Adjusted Scale Parameter") +
  theme_minimal()

ggplot(data.frame(thresholds, shape_params, shape_params_lower, shape_params_upper), aes(x = thresholds, y = shape_params)) +
  geom_point() +
  geom_errorbar(aes(ymin = shape_params_lower, ymax = shape_params_upper), width = 0.1, color = "blue") +
  geom_vline(xintercept = threshold_opt, color = "red", linetype = "dashed") +
  labs(title = "Shape Parameter vs. Threshold", x = "Threshold", y = "Shape Parameter") +
  theme_minimal()

fit_final <- fevd(x_log, threshold = threshold_opt, type = "GP")
plot(fit_final, type = "qq")

##############################
### PROFILE LIKELIHOOD ###
##############################

profile_likelihood <- function(z_star, y, threshold) {
  excesses <- y - threshold
  neg_log_likelihood <- function(sigma) {
    xi <- sigma / (threshold - z_star)
    if (xi >= 0) return(Inf)
    -sum(devd(excesses, loc = 0, scale = sigma, shape = xi, log = TRUE))
  }
  opt <- optim(par = 1, fn = neg_log_likelihood, method = "Brent", lower = 1e-5, upper = 1)
  return(-opt$value)
}

y_data <- x_log[x_log > threshold_opt]
z_grid <- seq(max(y_data) + 1e-4, 9.9, length.out = 100)
loglik_profile <- sapply(z_grid, function(z) profile_likelihood(z, y_data, threshold_opt))
loglik_profile <- (loglik_profile - max(loglik_profile)) * 2

ci_threshold <- -qchisq(0.9, df = 1)
z_star_CI <- z_grid[loglik_profile >= ci_threshold]
z_star_lower <- min(z_star_CI)
z_star_upper <- max(z_star_CI)
z_star_hat <- z_grid[which.max(loglik_profile)]

ggplot(data.frame(z_grid, loglik_profile), aes(x = z_grid, y = loglik_profile)) +
  geom_line(size = 1.2) +
  geom_hline(yintercept = ci_threshold, col = "red", linetype = "dashed") +
  geom_vline(xintercept = c(z_star_lower, z_star_upper), col = "blue", linetype = "dashed") +
  theme_minimal()

cat("Profile Likelihood Estimate:\n")
cat("z* =", z_star_hat, "log cm\n")
cat("90% CI:", z_star_lower, "-", z_star_upper, "log cm\n")
cat("Exponentiated:", exp(z_star_hat), "cm (", exp(z_star_lower), "-", exp(z_star_upper), ")\n")

##############################
### FULL BAYESIAN ANALYSIS ###
##############################

loglik <- function(y_star, xi) {
  if (xi >= 0 || any(y_star <= y_data)) return(-Inf)
  term1 <- length(y_data) * log(y_star - threshold_opt) / xi
  term2 <- -(1 / xi + 1) * sum(log(y_star - y_data))
  term3 <- -length(y_data) * log(-xi)
  return(term1 + term2 + term3)
}

y_grid <- seq(max(y_data) + 1e-4, 9.9, length.out = 1000)
xi_grid <- seq(-0.6, -1e-3, length.out = 1000)

marginal_log_posterior <- sapply(y_grid, function(y_star) {
  loglik_vals <- sapply(xi_grid, function(xi) loglik(y_star, xi))
  log_integrand <- loglik_vals
  integrand <- exp(log_integrand - max(log_integrand))
  log_integral <- log(trapz(xi_grid, integrand)) + max(log_integrand)
  return(log_integral)
})

prior <- rep(1, length(y_grid))
prior[y_grid <= max(y_data)] <- 0
unnormalized_posterior <- exp(marginal_log_posterior - max(marginal_log_posterior)) * prior
posterior_density <- unnormalized_posterior / trapz(y_grid, unnormalized_posterior)

n_bayes <- 100000
samples_bayes <- sample(y_grid, size = n_bayes, replace = TRUE, prob = posterior_density)

ci_bayes <- quantile(samples_bayes, probs = c(0.05, 0.95))
hpd_bayes <- hdi(samples_bayes, credMass = 0.90)
map_bayes <- y_grid[which.max(posterior_density)]

cat("\n--- Full Bayesian Inference ---\n")
cat("MAP estimate for y*:", round(map_bayes, 4), "\n")
cat("90% Equal-Tailed Credible Interval for y*: (", round(ci_bayes[1], 4), ",", round(ci_bayes[2], 4), ")\n")
cat("90% HPDI for y*: (", round(hpd_bayes[1], 4), ",", round(hpd_bayes[2], 4), ")\n")
cat("Exponentiated MAP:", round(exp(map_bayes), 2), "cm\n")
cat("Exponentiated Equal-Tailed CI: (", round(exp(ci_bayes[1]), 2), ",", round(exp(ci_bayes[2]), 2), ") cm\n")
cat("Exponentiated HPDI: (", round(exp(hpd_bayes[1]), 2), ",", round(exp(hpd_bayes[2]), 2), ") cm\n")

##############################
### PLOT BAYESIAN POSTERIOR ###
##############################

df_bayes <- data.frame(y_star = y_grid, posterior = posterior_density)

ggplot(df_bayes, aes(x = y_star, y = posterior)) +
  geom_line(color = "darkgreen", size = 1.2) +
  geom_histogram(data = data.frame(samples = samples_bayes), aes(x = samples, y = ..density..),
                 fill = "lightgreen", color = "black", bins = 40, alpha = 0.4) +
  geom_vline(xintercept = map_bayes, color = "purple", linetype = "dashed") +
  geom_vline(xintercept = ci_bayes, color = "darkgreen", linetype = "dotted") +
  geom_vline(xintercept = hpd_bayes, color = "blue", linetype = "dashed") +
  geom_vline(xintercept = max(y_data), color = "gray", linetype = "dotted") +
  labs(title = "Bayesian Posterior for y*", 
       x = "y* (log max length)", 
       y = "Posterior Density") +
  theme_minimal()
