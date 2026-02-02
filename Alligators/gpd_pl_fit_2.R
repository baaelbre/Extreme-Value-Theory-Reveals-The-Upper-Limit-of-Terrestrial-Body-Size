library(ggplot2)
library(extRemes)
library(ismev)

#######################
### PROFILE LIKELIHOOD FOR (z*, xi) ###
#######################

# Define function to compute the profile likelihood for given (z*, xi)
profile_likelihood <- function(z_star, xi, y, threshold) {
  excesses <- y[y > threshold] - threshold
  if (xi >= 0) return(-Inf)  # Ensure finite endpoint constraint
  
  # Compute scale parameter from reparametrized form
  sigma <- -xi * (threshold - z_star)
  if (sigma <= 0) return(-Inf)
  
  # Compute negative log-likelihood
  nll <- -sum(devd(excesses, loc = 0, scale = sigma, shape = xi, log = TRUE))
  return(nll)
}

# Set up a grid of (z*, xi) values
z_star_values <- seq(3.305, 3.5, length.out = 50)  # Adjust range if needed
xi_values <- seq(-0.5, -0.01, length.out = 50)  # Ensure negative xi
likelihood_matrix <- matrix(NA, nrow = length(z_star_values), ncol = length(xi_values))

# Compute likelihood for each pair (z*, xi)
for (i in seq_along(z_star_values)) {
  for (j in seq_along(xi_values)) {
    likelihood_matrix[i, j] <- profile_likelihood(z_star_values[i], xi_values[j], x_log, threshold_opt)
  }
}

# Normalize likelihood
max_likelihood <- max(likelihood_matrix, na.rm = TRUE)
likelihood_matrix <- likelihood_matrix - max_likelihood

# Compute confidence region threshold (-qchisq(0.95, df = 2) / 2)
ci_threshold <- -qchisq(0.95, df = 2) / 2

# Extract confidence region
conf_region <- which(likelihood_matrix >= ci_threshold, arr.ind = TRUE)

# Convert indices to parameter values
conf_z_star <- z_star_values[conf_region[, 1]]
conf_xi <- xi_values[conf_region[, 2]]

# Find MLE estimates
mle_index <- which(likelihood_matrix == 0, arr.ind = TRUE)
z_star_mle <- z_star_values[mle_index[1]]
xi_mle <- xi_values[mle_index[2]]

# Plot profile likelihood contours
likelihood_df <- expand.grid(z_star = z_star_values, xi = xi_values)
likelihood_df$likelihood <- as.vector(likelihood_matrix)

p <- ggplot(likelihood_df, aes(x = z_star, y = xi, z = likelihood)) +
  geom_contour_filled(breaks = seq(min(likelihood_matrix, na.rm = TRUE), 0, length.out = 20)) +
  geom_contour(breaks = ci_threshold, color = "red", size = 1.2) +
  geom_point(aes(x = z_star_mle, y = xi_mle), color = "black", size = 3) +
  labs(title = "Profile Likelihood Contours for (z*, xi)",
       x = "z* (Maximum log femur+humerus Circumference)",
       y = "Shape Parameter (xi)") +
  theme_minimal()

print(p)

# Print results
cat("Estimated (z*, xi):", z_star_mle, ",", xi_mle, "\n")
cat("95% Confidence Region for (z*, xi) shown in red contour.\n")
