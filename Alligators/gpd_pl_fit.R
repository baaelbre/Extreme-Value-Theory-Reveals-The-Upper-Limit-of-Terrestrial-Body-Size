library(ggplot2)
library(extRemes)
library(ismev)
library(numDeriv)
library(dplyr)
library(MASS)

##############################
### POT ANALYSIS ON LOG CFH ###
##############################

# Load the x_log object from the "Data/logCFH" file
x_log <- readRDS("Data/logCFH")


x_log_ <- sort(x_log)[-length(x_log)]  # Exclude Argentinosaurus

### Find optimal threshold by tracking parameter stability ###
# Set range of thresholds
u_0 <- 2.56
thresholds <- seq(u_0, 3.21, length.out = 30)

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

threshold_opt <- thresholds[18] #2.94

# Plot adjusted scale parameter vs. threshold
ggplot(data.frame(thresholds, adjusted_scale_params), aes(x = thresholds, y = adjusted_scale_params)) +
  geom_point() +
  geom_errorbar(aes(ymin = adjusted_scale_params_lower, ymax = adjusted_scale_params_upper), width=0.1, color="blue")+
  geom_vline(xintercept = thresholds[18], color="red", linetype="dashed")+
  labs(title = "Adjusted Scale Parameter vs. Threshold", x = "Threshold", y = "Adjusted Scale Parameter") +
  theme_minimal()+
  geom_hline(yintercept = adjusted_scale_params[18]-.05, color="red", linetype="dashed")

# Plot shape parameter vs. threshold
ggplot(data.frame(thresholds, shape_params, shape_params_lower, shape_params_upper), aes(x = thresholds, y = shape_params)) +
  geom_point() +
  geom_errorbar(aes(ymin = shape_params_lower, ymax = shape_params_upper), width = 0.1, color = "blue") +
  geom_vline(xintercept = thresholds[18], color = "red", linetype = "dashed") +
  labs(title = "Shape Parameter vs. Threshold", x = "Threshold", y = "Shape Parameter") +
  theme_minimal() +
  geom_hline(yintercept = shape_params[18]+0.04, color = "red", linetype = "dashed")

# Plot scale parameter vs. threshold with error bars
ggplot(data.frame(thresholds, scale_params, scale_params_lower, scale_params_upper), aes(x = thresholds, y = scale_params)) +
  geom_point() +
  geom_errorbar(aes(ymin = scale_params_lower, ymax = scale_params_upper), width = 0.1, color = "blue") +
  geom_vline(xintercept = thresholds[18], color = "red", linetype = "dashed") +
  labs(title = "Scale Parameter vs. Threshold", x = "Threshold", y = "Scale Parameter") +
  theme_minimal() +
  geom_segment(aes(x = thresholds[18], y = scale_params[18], xend = thresholds[30], yend = scale_params[18] + (scale_params[22] - scale_params[18]) / (thresholds[22] - thresholds[18]) * (thresholds[30] - thresholds[18])), color = "red", linetype = "solid", size = 1.5) +
  annotate("text", x = (thresholds[18] + thresholds[26]) / 2-.25, y = (scale_params[18] + (scale_params[22] - scale_params[18]) / (thresholds[22] - thresholds[18]) * (thresholds[26] - thresholds[18])) / 2, label = "Optimal Scale Line", color = "red", vjust = -1)

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
  opt <- optim(par = 1, fn = neg_log_likelihood, method = "Brent", lower = 1e-5, upper = 1)
  return(-opt$value)
}

# Set range for endpoint estimates
ez_star_values <- seq(3.305, 3.5, length.out = 100)
likelihood_values <- sapply(ez_star_values, function(z) profile_likelihood(z, x_log[x_log>threshold_opt], threshold_opt))

# Normalize likelihoods
likelihood_values <- (likelihood_values - max(likelihood_values))*2

# Compute likelihood ratio statistic for confidence interval
ci_threshold <- -qchisq(0.90, df = 1) 
z_star_CI <- ez_star_values[likelihood_values >= ci_threshold]

# Extract lower and upper bounds of confidence interval
z_star_lower <- min(z_star_CI)
z_star_upper <- max(z_star_CI)

# Plot profile likelihood
plot(ez_star_values, likelihood_values, type = "l", lwd = 2, xlab = "z* (Maximum log femur+humerus Circumference)", ylab = "Profile Log-Likelihood")
abline(h = ci_threshold, col = "red", lty = 2)
abline(v = c(z_star_lower, z_star_upper), col = "blue", lty = 2)

# Print results
cat("Estimated endpoint (z*):", ez_star_values[which.max(likelihood_values)], "\n")
cat("90% Confidence Interval: (", z_star_lower, ",", z_star_upper, ")\n")

# Plot profile likelihood for original scale using ggplot2
df <- data.frame(ez_star_values = 10^ez_star_values, likelihood_values = likelihood_values)

ggplot(df, aes(x = ez_star_values, y = likelihood_values)) +
  geom_line(size = 1.2) +
  geom_hline(yintercept = ci_threshold, col = "red", linetype = "dashed") +
  geom_vline(xintercept = c(10^z_star_lower, 10^z_star_upper), col = "blue", linetype = "dashed") +
  labs(title = "Profile Log-Likelihood for Original Scale", x = "10^y* (Maximum femur+humerus Circumference)", y = "Profile Log-Likelihood") +
  theme_minimal() +
  annotate("text", x = 10^ez_star_values[which.max(likelihood_values)], y = max(likelihood_values), label = paste("Estimated endpoint:", round(10^ez_star_values[which.max(likelihood_values)], 2)), color = "black", vjust = -1) +
  annotate("text", x = 10^z_star_lower, y = ci_threshold, label = paste("Lower CI:", round(10^z_star_lower, 2)), color = "blue", vjust = -1) +
  annotate("text", x = 10^z_star_upper, y = ci_threshold, label = paste("Upper CI:", round(10^z_star_upper, 2)), color = "blue", vjust = -1)

# Print results
cat("Estimated endpoint (10^z*):", 10^ez_star_values[which.max(likelihood_values)], "\n")
cat("90% Confidence Interval: (", 10^z_star_lower, ",", 10^z_star_upper, ")\n")

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
  opt <- optim(par = max(y) + 0.1, fn = neg_log_likelihood, method = "Brent", lower = max(y), upper = max(y) + 2)
  
  return(-opt$value)  # Return max log-likelihood for this xi
}

# Compute profile likelihood values for xi
xi_values <- seq(-0.6, -0.1, length.out = 100)
likelihood_values_xi <- sapply(xi_values, function(xi) profile_likelihood_xi(xi, x_log[x_log>threshold_opt], threshold_opt))

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
y_star_values <- seq(3.305, 3.55, length.out = 100)  # Finer grid for smoothness
xi_values <- seq(-0.6, -0.1, length.out = 100)  

# Compute profile likelihoods over the grid
likelihood_matrix <- outer(y_star_values, xi_values, Vectorize(function(y, xi) {
  likelihood_2D(y, xi, x_log[x_log>threshold_opt], threshold_opt)
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
       x = "y* (log max circumference)", y = "ξ (shape parameter)") +
  theme_minimal()

################################################
# Exceedance probability using joint confidence region
################################################

# Function to compute exceedance probability
exceedance_prob <- function(z_star, xi, x, u) {
  ((z_star - x) / (z_star - u))^(-1/xi)
}

# Threshold and exceedance level
threshold <- 2.94
x_arg <- 3.30

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
p_hat <- mean(p_samples)
p_CI <- quantile(p_samples, c(0.05, 0.95)) # 90% CI

# Print results
cat("Estimated exceedance probability P(X > 3.30) (from joint CI):", p_hat, "\n")
cat("90% Confidence Interval (from joint CI): (", p_CI[1], ",", p_CI[2], ")\n")

# Create a data frame for ggplot
df_plot <- data.frame(p_samples = p_samples)

# Create the ggplot histogram
ggplot(df_plot, aes(x = p_samples)) +
  geom_histogram(aes(y = after_stat(density)), # Use density for probability
                 bins = 20,
                 fill = "lightblue",
                 color = "black") +
  labs(title = "Probability of finding dinosaurs bigger than Argentinosaurus",
       x = "P(X > 3.30)",
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

###########################
# Mass estimation
###########################
library(readxl)
# Estimate regression coefficients
file_path <- "Data/campione_evans_2012.xlsx"
df <- read_excel(file_path, sheet = "Table 1.csv")

femur_col <- "Femur Circumference"
humerus_col <- "Humerus Circumference"
mass_col <- "Body Mass (g)"

# Remove rows with NA values
df <- na.omit(df)

# Compute the sum of femur and humerus circumference
df$sum_circumference <- df[[femur_col]] + df[[humerus_col]]

# Log-transform the sum of circumferences and mass
df$log_sum_circ <- log10(df$sum_circumference)
df$log_mass <- log10(df[[mass_col]])

# Fit linear and quadratic models
lin_model <- lm(log_mass ~ log_sum_circ, data = df)
df$log_sum_circ_sq <- df$log_sum_circ^2
quad_model <- lm(log_mass ~ log_sum_circ + log_sum_circ_sq, data = df)

sampled_lin_alpha <- rnorm(n_samples, mean = coef(lin_model)[1], sd = summary(lin_model)$coefficients[1, 2])
sampled_lin_beta <- rnorm(n_samples, mean = coef(lin_model)[2], sd = summary(lin_model)$coefficients[2, 2])

sampled_quad_alpha <- rnorm(n_samples, mean = coef(quad_model)[1], sd = summary(quad_model)$coefficients[1, 2])
sampled_quad_beta <- rnorm(n_samples, mean = coef(quad_model)[2], sd = summary(quad_model)$coefficients[2, 2])
sampled_quad_gamma <- rnorm(n_samples, mean = coef(quad_model)[3], sd = summary(quad_model)$coefficients[3, 2])

### Linear relationship ###

# Regression coefficients & standard errors
alpha_lin <- -1.104
beta <- 2.749
sigma_sq <- 0.134

alpha_lin_se <- 0.0339
beta_se <- 0.0197
x <- log10(2016)

# Theoretical variance using propagation rule
# for fixed measurements (i.e the Argentinosaurus)
log_mass_var_theory <- alpha_lin_se^2 + (x^2 * beta_se^2)+sigma_sq
log_mass <- alpha_lin+beta*x
log_mass_lwr <- log_mass-qnorm(0.95)*log_mass_var_theory^0.5
log_mass_upr <- log_mass+qnorm(0.95)*log_mass_var_theory^0.5
mass <- 10^log_mass/10^6
mass_lwr <- 10^log_mass_lwr/10^6
mass_upr <- 10^log_mass_upr/10^6

# Number of samples
n_samples <- 10000  

# Convert log-likelihood to sampling probabilities
df_likelihood$prob_weights <- df_likelihood$prob_weights / sum(df_likelihood$prob_weights)

# Sample endpoints proportional to likelihood
sampled_endpoints <- sample(df_likelihood$y_star, 
                            size = n_samples, replace = TRUE, 
                            prob = df_likelihood$prob_weights)

# Sample regression coefficients from normal distributions
sampled_alpha_lin <- rnorm(n_samples, mean = alpha_lin, sd = alpha_lin_se)
sampled_beta <- rnorm(n_samples, mean = beta, sd = beta_se)

# Compute log body mass
sampled_log_mass <- sampled_alpha_lin + sampled_beta * sampled_endpoints

# Convert to mass in tons (M in grams, so divide by 10^6)
sampled_mass_tons <- 10^(sampled_log_mass) / 1e6

# Compute 90% confidence intervals
log_mass_CI <- quantile(sampled_log_mass, probs = c(0.05, 0.95))
mass_CI <- quantile(sampled_mass_tons, probs = c(0.05, 0.95))

# Compute median estimates
log_mass_median <- median(sampled_log_mass)
mass_median <- median(sampled_mass_tons)

# Plot histogram of log mass
p1 <- ggplot(data.frame(log_mass = sampled_log_mass), aes(x = log_mass)) +
  geom_histogram(binwidth = 0.1, fill = "blue", alpha = 0.5, color = "black") +
  geom_vline(xintercept = log_mass_median, color = "red", linetype = "dashed", linewidth = 1.2) +
  geom_vline(xintercept = log_mass_CI, color = "black", linetype = "dotted", linewidth = 1.2) +
  annotate("text", x = log_mass_median, y = max(table(cut(sampled_log_mass, breaks = 30))) * 0.8, 
           label = paste0("Median: ", round(log_mass_median, 2)), color = "red", hjust = -0.1) +
  annotate("text", x = log_mass_CI[1], y = max(table(cut(sampled_log_mass, breaks = 30))) * 0.6, 
           label = paste0("90% CI: [", round(log_mass_CI[1], 2), ", ", round(log_mass_CI[2], 2), "]"), 
           color = "black", hjust = 0) +
  labs(title = "Distribution of Log Body Mass", x = "Log Mass (log10(grams))", y = "Frequency") +
  theme_minimal()

# Plot histogram of mass in tons
p2 <- ggplot(data.frame(mass_tons = sampled_mass_tons), aes(x = mass_tons)) +
  geom_histogram(binwidth = 5, fill = "purple", alpha = 0.5, color = "black") +
  geom_vline(xintercept = mass_median, color = "red", linetype = "dashed", linewidth = 1.2) +
  geom_vline(xintercept = mass_CI, color = "black", linetype = "dotted", linewidth = 1.2) +
  annotate("text", x = mass_median, y = max(table(cut(sampled_mass_tons, breaks = 30))) * 0.8, 
           label = paste0("Median: ", round(mass_median, 2), " tons"), color = "red", hjust = -0.1) +
  annotate("text", x = mass_CI[1], y = max(table(cut(sampled_mass_tons, breaks = 30))) * 0.6, 
           label = paste0("90% CI: [", round(mass_CI[1], 2), ", ", round(mass_CI[2], 2), "] tons"), 
           color = "black", hjust = 0) +
  labs(title = "Distribution of Maximum Body Mass (tons)", x = "Mass (tons)", y = "Frequency") +
  theme_minimal()

# Print plots
print(p1)
print(p2)

# Print results
cat("Estimated log mass (median):", log_mass_median, "\n")
cat("90% CI for log mass:", log_mass_CI[1], "to", log_mass_CI[2], "\n")
cat("Estimated mass (median):", mass_median, "tons\n")
cat("90% CI for mass:", mass_CI[1], "to", mass_CI[2], "tons\n")





### Campione (2017) ###
### Quadratic relationship ###

# Regression coefficients & standard errors
alpha_quad <- -1.25
beta <- 2.923
gamma <- -0.049
sigma_sq <- 0.134

alpha_quad_se <- 0.138
beta_se <- 0.161
gamma_se <- 0.0447

# Theoretical variance using propagation rule
log_mass_var_theory <- alpha_quad_se^2 + (x^2 * beta_se^2) + (x^4 * gamma_se^2)
log_mass <- alpha_quad+beta*x+gamma*x^2
log_mass_lwr <- log_mass-qnorm(0.95)*log_mass_var_theory^0.5
log_mass_upr <- log_mass+qnorm(0.95)*log_mass_var_theory^0.5
mass <- 10^log_mass/10^6
mass_lwr <- 10^log_mass_lwr/10^6
mass_upr <- 10^log_mass_upr/10^6

# Number of samples
n_samples <- 100000

# Convert log-likelihood to sampling probabilities
df_likelihood$prob_weights <- exp(df_likelihood$likelihood - max(df_likelihood$likelihood))

# Sample endpoints proportional to likelihood
sampled_endpoints <- sample(df_likelihood$y_star, 
                            size = n_samples, replace = TRUE, 
                            prob = df_likelihood$prob_weights)

# Sample regression coefficients from normal distributions
sampled_alpha_quad <- rnorm(n_samples, mean = alpha_quad, sd = alpha_quad_se)
sampled_beta <- rnorm(n_samples, mean = beta, sd = beta_se)
sampled_gamma <- rnorm(n_samples, mean = gamma, sd = gamma_se)

# Compute log body mass
sampled_log_mass <- sampled_alpha_quad + sampled_beta * sampled_endpoints + sampled_gamma * sampled_endpoints^2

# Convert to mass in tons (M in grams, so divide by 10^6)
sampled_mass_tons <- 10^(sampled_log_mass) / 1e6

# Compute 90% confidence intervals
log_mass_CI <- quantile(sampled_log_mass, probs = c(0.05, 0.95))
mass_CI <- quantile(sampled_mass_tons, probs = c(0.05, 0.95))

# Compute median estimates
log_mass_median <- median(sampled_log_mass)
mass_median <- median(sampled_mass_tons)

# Plot histogram of log mass
p1 <- ggplot(data.frame(log_mass = sampled_log_mass), aes(x = log_mass)) +
  geom_histogram(binwidth = 0.1, fill = "blue", alpha = 0.5, color = "black") +
  geom_vline(xintercept = log_mass_median, color = "red", linetype = "dashed", linewidth = 1.2) +
  geom_vline(xintercept = log_mass_CI, color = "black", linetype = "dotted", linewidth = 1.2) +
  annotate("text", x = log_mass_median, y = max(table(cut(sampled_log_mass, breaks = 30))) * 0.8, 
           label = paste0("Median: ", round(log_mass_median, 2)), color = "red", hjust = -0.1) +
  annotate("text", x = log_mass_CI[1], y = max(table(cut(sampled_log_mass, breaks = 30))) * 0.6, 
           label = paste0("90% CI: [", round(log_mass_CI[1], 2), ", ", round(log_mass_CI[2], 2), "]"), 
           color = "black", hjust = 0) +
  labs(title = "Distribution of Log Body Mass", x = "Log Mass (log10(grams))", y = "Frequency") +
  theme_minimal()

# Plot histogram of mass in tons
p2 <- ggplot(data.frame(mass_tons = sampled_mass_tons), aes(x = mass_tons)) +
  geom_histogram(binwidth = 5, fill = "purple", alpha = 0.5, color = "black") +
  geom_vline(xintercept = mass_median, color = "red", linetype = "dashed", linewidth = 1.2) +
  geom_vline(xintercept = mass_CI, color = "black", linetype = "dotted", linewidth = 1.2) +
  annotate("text", x = mass_median, y = 6000, 
           label = paste0("Median: ", round(mass_median, 2), " tons"), color = "red", hjust = -0.1) +
  annotate("text", x = mass_CI[1], y = 5000, 
           label = paste0("90% CI: [", round(mass_CI[1], 2), ", ", round(mass_CI[2], 2), "] tons"), 
           color = "black", hjust = 0) +
  scale_x_continuous(limits = c(0, 500))+
  scale_y_continuous(limits=c(0,7500))+
  labs(title = "Distribution of Maximum Body Mass, quadratic correction (tons)", x = "Mass (tons)", y = "Frequency") +
  theme_minimal()

# Print plots
print(p1)
print(p2)

# Print results
cat("Estimated log mass (median):", log_mass_median, "\n")
cat("90% CI for log mass:", log_mass_CI[1], "to", log_mass_CI[2], "\n")
cat("Estimated mass (median):", mass_median, "tons\n")
cat("90% CI for mass:", mass_CI[1], "to", mass_CI[2], "tons\n")


#####################################################################
# Now repeat everything for the dataset with argentinosaurus excluded
#####################################################################

#############################
### PROFILE LIKELIHOOD CI ###
#############################

# Define function to compute the profile likelihood for a given endpoint z*
profile_likelihood <- function(z_star, y, threshold) {
  
  # Fit a GPD model to the threshold exceedances
  excesses <- y-threshold
  
  # Define negative log-likelihood function
  neg_log_likelihood <- function(sigma) {
    xi <- sigma/(threshold-z_star)   # Reparameterization
    if (xi >= 0) return(Inf)   # Ensure finite endpoint constraint
    -sum(devd(excesses, loc = 0, scale = sigma, shape = xi, log = TRUE))
  }
  # Optimize sigma for a fixed z* (i.e. minimize neg_log_likelihood wrt sigma)
  opt <- optim(par = 1, fn = neg_log_likelihood, method = "Brent", lower = 1e-5, upper = 1)
  return(-opt$value)
}

# Set range for endpoint estimates
ez_star_values <- seq(3.27, 3.5, length.out = 100)
likelihood_values <- sapply(ez_star_values, function(z) profile_likelihood(z, x_log_, threshold_opt)) # Use x_log_

# Normalize likelihoods
likelihood_values <- (likelihood_values - max(likelihood_values))*2

# Compute likelihood ratio statistic for confidence interval
ci_threshold <- -qchisq(0.90, df = 1)
z_star_CI <- ez_star_values[likelihood_values >= ci_threshold]

# Extract lower and upper bounds of confidence interval
z_star_lower <- min(z_star_CI)
z_star_upper <- max(z_star_CI)

# Plot profile likelihood
plot(ez_star_values, likelihood_values, type = "l", lwd = 2, xlab = "z* (Maximum log femur+humerus Circumference)", ylab = "Profile Log-Likelihood")
abline(h = ci_threshold, col = "red", lty = 2)
abline(v = c(z_star_lower, z_star_upper), col = "blue", lty = 2)

# Print results
cat("Estimated endpoint (z*):", ez_star_values[which.max(likelihood_values)], "\n")
cat("90% Confidence Interval: (", z_star_lower, ",", z_star_upper, ")\n")

# Plot profile likelihood for original scale using ggplot2
df <- data.frame(ez_star_values = 10^ez_star_values, likelihood_values = likelihood_values)

ggplot(df, aes(x = ez_star_values, y = likelihood_values)) +
  geom_line(size = 1.2) +
  geom_hline(yintercept = ci_threshold, col = "red", linetype = "dashed") +
  geom_vline(xintercept = c(10^z_star_lower, 10^z_star_upper), col = "blue", linetype = "dashed") +
  labs(title = "Profile Log-Likelihood for Original Scale", x = "10^y* (Maximum femur+humerus Circumference)", y = "Profile Log-Likelihood") +
  theme_minimal() +
  annotate("text", x = 10^ez_star_values[which.max(likelihood_values)], y = max(likelihood_values), label = paste("Estimated endpoint:", round(10^ez_star_values[which.max(likelihood_values)], 2)), color = "black", vjust = -1) +
  annotate("text", x = 10^z_star_lower, y = ci_threshold, label = paste("Lower CI:", round(10^z_star_lower, 2)), color = "blue", vjust = -1) +
  annotate("text", x = 10^z_star_upper, y = ci_threshold, label = paste("Upper CI:", round(10^z_star_upper, 2)), color = "blue", vjust = -1)

# Print results
cat("Estimated endpoint (10^z*):", 10^ez_star_values[which.max(likelihood_values)], "\n")
cat("90% Confidence Interval: (", 10^z_star_lower, ",", 10^z_star_upper, ")\n")

################################################
# Profile likelihood for xi
################################################

profile_likelihood_xi <- function(xi, y, threshold) {
  # Negative log-likelihood function for given xi
  neg_log_likelihood <- function(z_star) {
    sigma <- (threshold-z_star) * xi  # Reparameterization
    if (xi >= 0) return(Inf)   # Ensure valid parameters
    -sum(devd(y, loc = threshold, scale = sigma, shape = xi, log = TRUE))
  }
  
  # Optimize z_star
  opt <- optim(par = max(y) + 0.1, fn = neg_log_likelihood, method = "Brent", lower = max(y), upper = max(y) + 2)
  
  return(-opt$value)  # Return max log-likelihood for this xi
}

# Compute profile likelihood values for xi
xi_values <- seq(-0.6, -0.1, length.out = 100)
likelihood_values_xi <- sapply(xi_values, function(xi) profile_likelihood_xi(xi, x_log_, threshold_opt)) # Use x_log_

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

################################################
# Obtain joint confidence regions for z^* and xi
################################################

likelihood_2D <- function(y_star, xi, y, threshold) {
  # Compute scale parameter based on the reparameterization
  sigma <- -(y_star - threshold) * xi
  if (sigma <= 0 || xi >= 0) return(-Inf)   # Ensure valid GPD parameters
  
  # Compute the log-likelihood
  return(sum(devd(y, loc = threshold, scale = sigma, shape = xi, log = TRUE)))
}

# Define parameter ranges
y_star_values <- seq(3.305, 3.5, length.out = 100)   # Finer grid for smoothness
xi_values <- seq(-0.6, -0.01, length.out = 100)

# Compute profile likelihoods over the grid
likelihood_matrix <- outer(y_star_values, xi_values, Vectorize(function(y, xi) {
  likelihood_2D(y, xi, x_log_, threshold_opt) # Use x_log_
}))

# Normalize likelihoods
likelihood_matrix <- (likelihood_matrix - max(likelihood_matrix)) * 2

# Compute confidence region threshold
ci_threshold_2D <- -qchisq(0.90, df = 2)   # 90% confidence region threshold

# Convert data to dataframe for ggplot
df_likelihood <- expand.grid(y_star = y_star_values, xi = xi_values)
df_likelihood$likelihood <- as.vector(likelihood_matrix)

# Plot: Fill only the 90% confidence contour region
ggplot(df_likelihood, aes(x = y_star, y = xi, z = likelihood)) +
  geom_contour_filled(breaks = c(ci_threshold_2D, 0), alpha=0.5) +   # Fill only inside CI
  scale_fill_manual(values="hotpink", name = "90% Confidence Region") +   # Pink fill
  geom_contour(color = "black", breaks = c(ci_threshold_2D)) +   # Black border
  labs(title = "Joint 90% Confidence Region for (y*, ξ)",
       x = "y* (log max circumference)", y = "ξ (shape parameter)") +
  theme_minimal()

################################################
# Exceedance probability using joint confidence region
################################################

# Function to compute exceedance probability
exceedance_prob <- function(z_star, xi, x, u) {
  ((z_star - x) / (z_star - u))^(-1/xi)
}

# Threshold and exceedance level
threshold <- 2.94
x_arg <- 3.30

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
p_hat <- mean(p_samples)
p_CI <- quantile(p_samples, c(0.05, 0.95)) # 90% CI

# Print results
cat("Estimated exceedance probability P(X > 3.30) (from joint CI):", p_hat, "\n")
cat("90% Confidence Interval (from joint CI): (", p_CI[1], ",", p_CI[2], ")\n")

# Create a data frame for ggplot
df_plot <- data.frame(p_samples = p_samples)

# Create the ggplot histogram
ggplot(df_plot, aes(x = p_samples)) +
  geom_histogram(aes(y = after_stat(density)),   # Use density for probability
                 bins = 20,
                 fill = "lightblue",
                 color = "black") +
  labs(title = "Probability of finding dinosaurs bigger than Argentinosaurus (Excluding Argentinosaurus)",
       x = "P(X > 3.30)",
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
    y = max(hist(p_samples, plot = FALSE)$density) * 0.9,   # Approximate y
    label = paste0(round(p_hat * 100, 2), "%"),
    color = "blue",
    hjust = 0,   # Adjust text position
    vjust = 0    # Adjust text position
  ) +
  annotate(
    "text",
    x = p_CI[1],
    y = max(hist(p_samples, plot = FALSE)$density) * 0.9,   # Approximate y
    label = paste0(round(p_CI[1] * 100, 2), "%"),
    color = "red",
    hjust = 0,
    vjust = 0
  ) +
  annotate(
    "text",
    x = p_CI[2],
    y = max(hist(p_samples, plot = FALSE)$density) * 0.9,   # Approximate y
    label = paste0(round(p_CI[2] * 100, 2), "%"),
    color = "red",
    hjust = 0,
    vjust = 0
  )
