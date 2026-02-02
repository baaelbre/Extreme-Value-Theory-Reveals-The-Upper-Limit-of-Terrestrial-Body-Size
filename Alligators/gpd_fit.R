library(ggplot2)
library(extRemes)
####################
### POT ANALYSIS ###
####################

### On logs ###
x_log <- readRDS("Data/logCFH")

# Set range of thresholds
u_0 <- 5.9*log10(exp(1))
thresholds <- seq(u_0, 7.4*log10(exp(1)), length.out = 30)

# Initialize vectors to store the scale and shape parameters
scale_params <- numeric(length(thresholds))
shape_params <- numeric(length(thresholds))
shape_params_lower <- numeric(length(thresholds))
shape_params_upper <- numeric(length(thresholds))
scale_params_lower <- numeric(length(thresholds))
scale_params_upper <- numeric(length(thresholds))
adjusted_scale_params <- numeric(length(thresholds))
adjusted_scale_params_lower <- numeric(length(thresholds))
adjusted_scale_params_upper <- numeric(length(thresholds))

# Loop through thresholds and fit models to extract scale and shape
for (i in seq_along(thresholds)) {
  fit <- fevd(x_log, threshold = thresholds[i], type = "GP")
  scale_params[i] <- fit$results$par["scale"]
  shape_params[i] <- fit$results$par["shape"]
  adjusted_scale_params[i] <- scale_params[i] - shape_params[i] * (thresholds[i] - u_0)
  
  shape_params_lower[i] <- shape_params[i]
  shape_params_upper[i] <- shape_params[i]
  scale_params_lower[i] <- scale_params[i]
  scale_params_upper[i] <- scale_params[i]
  adjusted_scale_params_lower[i] <- adjusted_scale_params[i]
  adjusted_scale_params_upper[i] <- adjusted_scale_params[i]
  
  # Try to calculate the confidence intervals, and handle any errors
  tryCatch({
    ci <- ci(fit, type='parameter')
    shape_params_lower[i] <- ci[2,1]
    shape_params_upper[i] <- ci[2,3]
    scale_params_lower[i] <- ci[1,1]
    scale_params_upper[i] <- ci[1,3]
    
  }, error = function(e) {
    # Do nothing on error
  })
}

library(ggplot2)

# Plot adjusted scale parameter vs. threshold
ggplot(data.frame(thresholds, adjusted_scale_params), aes(x = thresholds, y = adjusted_scale_params)) +
  geom_point() +
  labs(title = "Adjusted Scale Parameter vs. Threshold", x = "Threshold", y = "Adjusted Scale Parameter") +
  theme_minimal()

# Plot shape parameter vs. threshold with error bars


ggplot(data.frame(thresholds, shape_params, shape_params_lower, shape_params_upper), aes(x = thresholds, y = shape_params)) +
  geom_point() +
  geom_errorbar(aes(ymin = shape_params_lower, ymax = shape_params_upper), width = 0.1, color = "blue") +
  geom_vline(xintercept = thresholds[18], color = "red", linetype = "dashed") +
  labs(title = "Shape Parameter vs. Threshold", x = "Threshold", y = "Shape Parameter") +
  theme_minimal() +
  geom_hline(yintercept = shape_params[18], color = "red", linetype = "dashed") +
  annotate("text", x = thresholds[18]-.25, y = shape_params[18], label = paste("Optimal Shape:", round(shape_params[18], 2)), color = "red", vjust = -1)

# Plot scale parameter vs. threshold with error bars
ggplot(data.frame(thresholds, scale_params, scale_params_lower, scale_params_upper), aes(x = thresholds, y = scale_params)) +
  geom_point() +
  geom_errorbar(aes(ymin = scale_params_lower, ymax = scale_params_upper), width = 0.1, color = "blue") +
  geom_vline(xintercept = thresholds[18], color = "red", linetype = "dashed") +
  labs(title = "Scale Parameter vs. Threshold", x = "Threshold", y = "Scale Parameter") +
  theme_minimal() +
  geom_segment(aes(x = thresholds[18], y = scale_params[18], xend = thresholds[30], yend = scale_params[18] + (scale_params[22] - scale_params[18]) / (thresholds[22] - thresholds[18]) * (thresholds[30] - thresholds[18])), color = "red", linetype = "solid", size = 1.5) +
  annotate("text", x = (thresholds[18] + thresholds[26]) / 2-.25, y = (scale_params[18] + (scale_params[22] - scale_params[18]) / (thresholds[22] - thresholds[18]) * (thresholds[26] - thresholds[18])) / 2, label = "Optimal Scale Line", color = "red", vjust = -1)

threshold_opt <- thresholds[18]
fitGP <- fevd(x_log, threshold = threshold_opt, type = 'GP')
# already
fitGP
par(mfrow=c(1,1))
plot(fitGP)
ci(fitGP, type="parameter")

modsum <- summary(fitGP)
# endpoint
scale <- as.double(fitGP$results[[1]][1])
shape <- as.double(fitGP$results[[1]][2])
scale_se <- as.double(modsum$se.theta[['scale']])
shape_se <- as.double(modsum$se.theta[['shape']])

# endpoint is 100% quantile , of dus sf^-1(0)

endpoint <- threshold_opt-scale/shape # point estimate

# Set the number of bootstrap samples
n_bootstrap <- 100000

# Initialize vectors to store bootstrap estimates
bootstrap_endpoints <- numeric(n_bootstrap)

# Perform bootstrapping
set.seed(123) # for reproducibility
for (i in 1:n_bootstrap) {
  # Generate bootstrap samples for scale and shape
  scale_bootstrap <- rnorm(1, mean = scale, sd = scale_se)
  shape_bootstrap <- rnorm(1, mean = shape, sd = shape_se)
  
  # Calculate the bootstrap endpoint
  bootstrap_endpoints[i] <- threshold_opt - scale_bootstrap / shape_bootstrap
}

# Calculate the 95% confidence interval
ci_lower <- quantile(bootstrap_endpoints, 0.05)
ci_upper <- quantile(bootstrap_endpoints, 0.95)

# Print the confidence interval
cat("90% Confidence Interval for the log endpoint: [", ci_lower, ", ", ci_upper, "]\n")
cat("90% Confidence Interval for the endpoint: [", 10^ci_lower, ", ", 10^ci_upper, "]\n")
fitGP$threshold

qevd(0.999, loc=fitGP$threshold, scale=fitGP$results$par["scale"], shape=fitGP$results$par["shape"])
prob_exceed_argentinosaurus <- 1 - pevd(log10(2016.222), loc=fitGP$threshold, scale=fitGP$results$par["scale"], shape=fitGP$results$par["shape"])
prob_exceed_argentinosaurus

prob_exceed_brachiosaurus <- 1 - pevd(log10(1870), loc=fitGP$threshold, scale=fitGP$results$par["scale"], shape=fitGP$results$par["shape"])
prob_exceed_brachiosaurus

# two dinosaurs of 1870 and 2016 are lucky shots and
# should not be expected to be found in the future

ggplot() +
  geom_histogram(aes(x = bootstrap_endpoints, y = after_stat(density)), bins = 30, fill = "lightblue", color = "black") +
  geom_vline(xintercept = endpoint, color = "red", linetype = "dashed") +
  geom_vline(xintercept = ci_lower, color = "blue", linetype = "dashed") +
  geom_vline(xintercept = ci_upper, color = "blue", linetype = "dashed") +
  labs(title = "Bootstrap Distribution of the Endpoint", x = "Endpoint of ln C", y = "Density") +
  annotate("text", x = endpoint, y = 0.1, label = paste("Endpoint:", round(endpoint, 2)), color = "red", vjust = -1) +
  annotate("text", x = ci_lower, y = 0.1, label = paste("CI Lower:", round(ci_lower, 2)), color = "blue", vjust = -1) +
  annotate("text", x = ci_upper, y = 0.1, label = paste("CI Upper:", round(ci_upper, 2)), color = "blue", vjust = -1)
# Now plot the distribution of exp(endpoint) instead, to get the maximal bones circumference
ggplot() +
  geom_histogram(aes(x = 10^(bootstrap_endpoints), y = after_stat(density)), bins = 30, fill = "lightblue", color = "black") +
  geom_vline(xintercept = 10^(endpoint), color = "red", linetype = "dashed") +
  geom_vline(xintercept = 10^(ci_lower), color = "blue", linetype = "dashed") +
  geom_vline(xintercept = 10^(ci_upper), color = "blue", linetype = "dashed") +
  labs(title = "Bootstrap Distribution of the maximal femur + humerus Circumference", x = "$C_{F+H}$ (mm)", y = "Density")+
  annotate("text", x = 10^(endpoint), y = 0.0015, label = paste("Endpoint:", round(10^(endpoint), 2)), color = "red", vjust = 0.2) +
  annotate("text", x = 10^(ci_lower), y = 0.0015, label = paste("CI Lower:", round(10^(ci_lower), 2)), color = "blue", vjust = 1.8) +
  annotate("text", x = 10^(ci_upper), y = 0.0015, label = paste("CI Upper:", round(10^(ci_upper), 2)), color = "blue", vjust = 1.8)

# Model diagnostics
# QQ plot for model diagnostics
plot(fitGP, type='qq')
plot(fitGP, type='probprob')
