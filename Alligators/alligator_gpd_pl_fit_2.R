library(ggplot2)
library(extRemes)
library(ismev)
library(numDeriv)
library(dplyr)
library(MASS)
library(readxl)

##############################
### POT ANALYSIS ON LAG SVL ###
##############################

# Load the RData file
load("gator_growth.rdata")

# Extract g.in$lagSVL (snout-vent length)
x_lagSVL <- g.in$lagSVL

# oke probleem, slechts 2 mannetjes groter dan 3,6 meter. Slechte dataset (South Carolina)
# Florida has the big mammas.

# Original scale plot
ggplot(data.frame(x), aes(x = x)) +
  geom_histogram(bins = 30, fill = "lightblue", color = "black") +
  labs(title = "Histogram of lagSVL (Original Scale)",
       x = "lagSVL",
       y = "Frequency") +
  theme_minimal()

# Log scale plot
ggplot(data.frame(log_lagSVL = log(x)), aes(x = log_lagSVL)) +
  geom_histogram(bins = 30, fill = "lightgreen", color = "black") +
  labs(title = "Histogram of log(lagSVL)",
       x = "log(lagSVL)",
       y = "Frequency") +
  theme_minimal()

### Find optimal threshold by tracking parameter stability ###
# Set range of thresholds
u_0 <- exp(5.13)
thresholds <- seq(u_0, exp(5.27), length.out = 30) # Adjust range

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
  fit <- fevd(x_lagSVL, threshold = thresholds[i], type = "GP")
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
    ci <- ci(fit, type='parameter', alpha=0.1)
    cov <- summary(fit)$cov.theta
    shape_params_lower[i] <- ci[2,1]
    shape_params_upper[i] <- ci[2,3]
    scale_params_lower[i] <- ci[1,1]
    scale_params_upper[i] <- ci[1,3]
    var_adjusted_scale_param <- cov[1,1]+cov[2,2]*(thresholds[i] - u_0)^2
    adjusted_scale_params_lower[i] <- adjusted_scale_params[i]-qnorm(0.95)*var_adjusted_scale_param^0.5
    adjusted_scale_params_upper[i] <- adjusted_scale_params[i]+qnorm(0.95)*var_adjusted_scale_param^0.5
  }, error = function(e) {
    # Do nothing on error
  })
}

threshold_opt <- thresholds[which.min(abs(adjusted_scale_params - median(adjusted_scale_params)))] # Find threshold where adjusted scale is close to median

# Plot adjusted scale parameter vs. threshold
ggplot(data.frame(thresholds, adjusted_scale_params), aes(x = thresholds, y = adjusted_scale_params)) +
  geom_point() +
  geom_errorbar(aes(ymin = adjusted_scale_params_lower, ymax = adjusted_scale_params_upper), width=0.1, color="blue")+
  geom_vline(xintercept = threshold_opt, color="red", linetype="dashed")+
  labs(title = "Adjusted Scale Parameter vs. Threshold", x = "Threshold", y = "Adjusted Scale Parameter") +
  theme_minimal()+
  geom_hline(yintercept = adjusted_scale_params[which(thresholds == threshold_opt)], color="red", linetype="dashed")

# Plot shape parameter vs. threshold
ggplot(data.frame(thresholds, shape_params, shape_params_lower, shape_params_upper), aes(x = thresholds, y = shape_params)) +
  geom_point() +
  geom_errorbar(aes(ymin = shape_params_lower, ymax = shape_params_upper), width = 0.1, color = "blue") +
  geom_vline(xintercept = threshold_opt, color = "red", linetype = "dashed") +
  labs(title = "Shape Parameter vs. Threshold", x = "Threshold", y = "Shape Parameter") +
  theme_minimal() +
  geom_hline(yintercept = shape_params[which(thresholds == threshold_opt)], color = "red", linetype = "dashed")

# QQ plot
qq_plot <- fevd(x_lagSVL, threshold = threshold_opt, type = "GP")
plot(qq_plot, type = "qq")

#############################
### PROFILE LIKELIHOOD CI ###
#############################

# Define function to compute the profile likelihood for a given endpoint z*
profile_likelihood <- function(z_star, y, threshold) {
  # Fit a GPD model to the threshold exceedances
  excesses <- y-threshold
  
  # Define negative log-likelihood function
  neg_log_likelihood <- function(sigma) {
    xi <- sigma/(threshold-z_star)  # Reparameterization
    if (xi >= 0) return(Inf)  # Ensure finite endpoint constraint
    -sum(devd(excesses, loc = 0, scale = sigma, shape = xi, log = TRUE))
  }
  # Optimize sigma for a fixed z* (i.e. minimize neg_log_likelihood wrt sigma)
  opt <- optim(par = 1, fn = neg_log_likelihood, method = "Brent", 
               lower = 1e-5, upper = max(y)-threshold)  # Adjust upper bound
  
  return(-opt$value)
}

# Set range for endpoint estimates
x_lagSVL_exceed <- x_lagSVL[x_lagSVL > threshold_opt]
ez_star_values <- seq(max(x_lagSVL_exceed)+0.01, max(x_lagSVL_exceed) + (max(x_lagSVL_exceed)-threshold_opt)*1.5, length.out = 100) #Adjust max change

likelihood_values <- sapply(ez_star_values, function(z) profile_likelihood(z, x_lagSVL_exceed, threshold_opt))

# Normalize likelihoods
likelihood_values <- (likelihood_values - max(likelihood_values))*2

# Compute likelihood ratio statistic for confidence interval
ci_threshold <- -qchisq(0.90, df = 1) 
z_star_CI <- ez_star_values[likelihood_values >= ci_threshold]

# Extract lower and upper bounds of confidence interval
z_star_lower <- min(z_star_CI)
z_star_upper <- max(z_star_CI)

# Plot profile likelihood
plot(ez_star_values, likelihood_values, type = "l", lwd = 2, xlab = "z* (Maximum lagSVL)", ylab = "Profile Log-Likelihood")
abline(h = ci_threshold, col = "red", lty = 2)
abline(v = c(z_star_lower, z_star_upper), col = "blue", lty = 2)

# Print results
cat("Estimated endpoint (z*):", ez_star_values[which.max(likelihood_values)], "\n")
cat("90% Confidence Interval: (", z_star_lower, ",", z_star_upper, ")\n")

# Plot profile likelihood for original scale using ggplot2
df <- data.frame(ez_star_values = ez_star_values, likelihood_values = likelihood_values)

ggplot(df, aes(x = ez_star_values, y = likelihood_values)) +
  geom_line(size = 1.2) +
  geom_hline(yintercept = ci_threshold, col = "red", linetype = "dashed") +
  geom_vline(xintercept = c(z_star_lower, z_star_upper), col = "blue", linetype = "dashed") +
  labs(title = "Profile Log-Likelihood for Original Scale", x = "y* (Maximum lagSVL)", y = "Profile Log-Likelihood") +
  theme_minimal() +
  annotate("text", x = ez_star_values[which.max(likelihood_values)], y = max(likelihood_values), label = paste("Estimated endpoint:", round(ez_star_values[which.max(likelihood_values)], 2)), color = "black", vjust = -1) +
  annotate("text", x = z_star_lower, y = ci_threshold, label = paste("Lower CI:", round(z_star_lower, 2)), color = "blue", vjust = -1) +
  annotate("text", x = z_star_upper, y = ci_threshold, label = paste("Upper CI:", round(z_star_upper, 2)), color = "blue", vjust = -1)

# Print results
cat("Estimated endpoint (y*):", ez_star_values[which.max(likelihood_values)], "\n")
cat("90% Confidence Interval: (", z_star_lower, ",", z_star_upper, ")\n")

################################################
# Profile likelihood for xi
################################################

profile_likelihood_xi <- function(xi, y, threshold) {
  # Negative log-likelihood function for given xi
  neg_log_likelihood <- function(z_star) {
    sigma <- (threshold-z_star) * xi  # Reparameterization
    if (xi >= 0) return(Inf)  # Ensure valid parameters
    -sum(devd(y, loc = threshold, scale = sigma, shape = xi, log = TRUE))
  }
  
  # Optimize z_star
  opt <- optim(par = max(y) + 0.1, fn = neg_log_likelihood, 
               method = "Brent", lower = max(y) - (max(y)-threshold), upper = max(y) + (max(y)-threshold)*1.5) #Adjust bounds
  
  return(-opt$value)  # Return max log-likelihood for this xi
}

# Compute profile likelihood values for xi
xi_values <- seq(-0.6, -0.05, length.out = 100)
likelihood_values_xi <- sapply(xi_values, function(xi) profile_likelihood_xi(xi, x_lagSVL_exceed, threshold_opt))

# Normalize likelihoods
likelihood_values_xi <- (likelihood_values_xi - max(likelihood_values_xi)) * 2

# Compute confidence intervals
ci_threshold <- -qchisq(0.90, df = 1)

xi_CI <- xi_values[likelihood_values_xi >= ci_threshold]

xi_lower <- min(xi_CI)
xi_upper <- max(xi_CI)

# Plot profile likelihood for xi
df_xi <- data.frame(xi_values, likelihood_values_xi)

ggplot(df_xi, aes(x = xi_values, y = likelihood_values_xi)) +
  geom_line(size = 1.2) +
  geom_hline(yintercept = ci_threshold, col = "red", linetype = "dashed") +
  geom_vline(xintercept = c(xi_lower, xi_upper), col = "blue", linetype = "dashed") +
  labs(title = "Profile Log-Likelihood for xi", x = "xi (Shape parameter)", y = "Profile Log-Likelihood") +
  theme_minimal()
cat("Estimated xi:", xi_values[which.max(likelihood_values_xi)], "\n")
cat("90% Confidence Interval for xi: (", xi_lower, ",", xi_upper, ")\n")

################
# DEBUGGING
################

# Debug a single z_star value
profile_likelihood(ez_star_values[which.max(likelihood_values)], x_lagSVL_exceed, threshold_opt)

# Debug a single xi value
profile_likelihood_xi(xi_values[which.max(likelihood_values_xi)], x_lagSVL_exceed, threshold_opt)

plot(ez_star_values, likelihood_values, type = "l")

################################################
# Obtain joint confidence regions for z^* and xi
################################################

likelihood_2D <- function(y_star, xi, y, threshold) {
  # Compute scale parameter based on the reparameterization
  sigma <- -(y_star - threshold) * xi
  if (sigma <= 0 || xi >= 0) return(-Inf)  # Ensure valid GPD parameters
  
  # Compute the log-likelihood
  return(sum(devd(y, loc = threshold, scale = sigma, shape = xi, log = TRUE)))
}

# Define parameter ranges
y_star_values <- seq(min(ez_star_values), max(ez_star_values), length.out = 100)  # Finer grid for smoothness
xi_values <- seq(min(xi_values), max(xi_values), length.out = 100)  

# Compute profile likelihoods over the grid
likelihood_matrix <- outer(y_star_values, xi_values, Vectorize(function(y, xi) {
  likelihood_2D(y, xi, x_lagSVL_exceed, threshold_opt)
}))

# Normalize likelihoods
likelihood_matrix <- (likelihood_matrix - max(likelihood_matrix)) * 2

# Compute confidence region threshold
ci_threshold_2D <- -qchisq(0.90, df = 2)  # 90% confidence region threshold

# Convert data to dataframe for ggplot
df_likelihood <- expand.grid(y_star = y_star_values, xi = xi_values)
df_likelihood$likelihood <- as.vector(likelihood_matrix)

# Plot: Fill only the 90% confidence contour region
ggplot(df_likelihood, aes(x = y_star, y = xi, z = likelihood)) +
  geom_contour_filled(breaks = c(ci_threshold_2D, 0), alpha=0.5) +  # Fill only inside CI
  scale_fill_manual(values="hotpink", name = "90% Confidence Region") +  # Pink fill
  geom_contour(color = "black", breaks = c(ci_threshold_2D)) +  # Black border
  labs(title = "Joint 90% Confidence Region for (y*, ξ)",
       x = "y* (max lagSVL)", y = "ξ (shape parameter)") +
  theme_minimal()

################################################
# Exceedance probability using joint confidence region
################################################

# Function to compute exceedance probability
exceedance_prob <- function(z_star, xi, x, u) {
  ((z_star - x) / (z_star - u))^(-1/xi)
}

# Threshold and exceedance level
threshold <- threshold_opt
x_arg <- max(x_lagSVL_exceed) + (max(x_lagSVL_exceed)-threshold_opt)*0.5 #Example test value. Adjust as needed

# Filter for points within the confidence region
confidence_points <- df_likelihood[df_likelihood$likelihood >= -qchisq(0.99, df = 2), ]

# Number of samples
n_sim <- 10000

# Sample from the confidence region
sample_indices <- sample(1:nrow(confidence_points), n_sim, replace = TRUE)
sampled_y_star <- confidence_points$y_star[sample_indices]
sampled_xi <- confidence_points$xi[sample_indices]

# Compute exceedance probabilities for all samples
p_samples <- exceedance_prob(sampled_y_star, sampled_xi, x_arg, threshold)

# Compute exceedance probability estimate and confidence interval
p_hat <- mean(p_samples, na.rm = TRUE)
p_CI <- quantile(p_samples, c(0.05, 0.95), na.rm=TRUE) # 90% CI

# Print results
cat("Estimated exceedance probability P(X > ", x_arg, ") (from joint CI):", p_hat, "\n")
cat("90% Confidence Interval (from joint CI): (", p_CI[1], ",", p_CI[2], ")\n")

# Create a data frame for ggplot
df_plot <- data.frame(p_samples = p_samples)

# Create the ggplot histogram
ggplot(df_plot, aes(x = p_samples)) +
  geom_histogram(aes(y = after_stat(density)), # Use density for probability
                 bins = 20,
                 fill = "lightblue",
                 color = "black") +
  labs(title = paste("Probability of finding alligator with lagSVL >", round(x_arg, 2)),
       x = paste("P(X >", round(x_arg, 2), ")"),
       y = "Density") +
  theme_minimal() +
  
  # Add vertical lines for the confidence interval
  geom_vline(xintercept = p_CI, color = "red", linetype = "dashed", size = 1) +
  
  # Add vertical line for the mean exceedance probability
  geom_vline(xintercept = p_hat, color = "blue", linetype = "solid", size = 1) +
  
  # Add text labels
  annotate(
    "text",
    x = p_hat,
    y = max(hist(p_samples, plot = FALSE)$density) * 0.9, # Approximate y
    label = paste0(round(p_hat * 100, 2), "%"),
    color = "blue",
    hjust = 0, # Adjust text position
    vjust = 0  # Adjust text position
  ) +
  annotate(
    "text",
    x = p_CI[1],
    y = max(hist(p_samples, plot = FALSE)$density) * 0.9, # Approximate y
    label = paste0(round(p_CI[1] * 100, 2), "%"),
    color = "red",
    hjust = 0,
    vjust = 0
  ) +
  annotate(
    "text",
    x = p_CI[2],
    y = max(hist(p_samples, plot = FALSE)$density) * 0.9, # Approximate y
    label = paste0(round(p_CI[2] * 100, 2), "%"),
    color = "red",
    hjust = 0,
    vjust = 0
  )