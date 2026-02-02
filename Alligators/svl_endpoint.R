### LIBRARIES ###
library(ggplot2)
library(extRemes)
library(numDeriv)
library(dplyr)
library(MASS)
library(readxl)
library(pracma)
library(HDInterval)
library(gridExtra)

### SETTINGS ###
ci_level <- 0.95
n_samples <- 10000
svl_stokes <- 236
tl_stokes <- 450

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

### LOAD DATA ###
df <- read_excel("Data/experimental_alligator_harvest_woodward.xlsx")
df <- df %>% filter(!is.na(SVL))
df_clean <- df %>% filter(Sex == "M", Deform == 0, !is.na(SVL), !is.na(TL))
df_clean$log_SVL <- log(df_clean$SVL)
df_clean$log_TL <- log(df_clean$TL)
df$log_SVL <- log(df$SVL)
df$log_TL <- log(df$TL)

# Threshold sequence
u_0 <- log(165)
thresholds_male <- seq(u_0 - 0.3, u_0 + 0.15, length.out = 50)

# Preallocate
scale_params <- shape_params <- shape_params_lower <- shape_params_upper <- numeric(length(thresholds_male))
scale_params_lower <- scale_params_upper <- adjusted_scale_params <- adjusted_scale_params_lower <- adjusted_scale_params_upper <- numeric(length(thresholds_male))

for (i in seq_along(thresholds_male)) {
  tryCatch({
    fit <- fevd(df_clean$log_SVL, threshold = thresholds_male[i], type = "GP")
    scale_params[i] <- fit$results$par["scale"]
    shape_params[i] <- fit$results$par["shape"]
    adjusted_scale_params[i] <- scale_params[i] - shape_params[i] * (thresholds_male[i] - u_0)
    
    ci <- ci(fit, type = 'parameter', alpha = 0.1)
    cov <- summary(fit)$cov.theta
    shape_params_lower[i] <- ci[2,1]
    shape_params_upper[i] <- ci[2,3]
    scale_params_lower[i] <- ci[1,1]
    scale_params_upper[i] <- ci[1,3]
    
    var_adj <- cov[1,1] + cov[2,2] * (thresholds_male[i] - u_0)^2
    adjusted_scale_params_lower[i] <- adjusted_scale_params[i] - qnorm(0.95) * sqrt(var_adj)
    adjusted_scale_params_upper[i] <- adjusted_scale_params[i] + qnorm(0.95) * sqrt(var_adj)
  }, error = function(e) {})
}

p_adj_scale_male <- ggplot(data.frame(thresholds_male, adjusted_scale_params, 
                                      adjusted_scale_params_lower, adjusted_scale_params_upper),
                           aes(x = thresholds_male, y = adjusted_scale_params)) +
  geom_point() +
  geom_errorbar(aes(ymin = adjusted_scale_params_lower, ymax = adjusted_scale_params_upper), 
                width = 0.05, color = "blue") +
  geom_vline(xintercept = u_0, color = "red", linetype = "dashed") +
  geom_hline(yintercept = adjusted_scale_params[which.min(abs(thresholds_male - u_0))], 
             color = "red", linetype = "dashed") +
  labs(x = "Threshold (log SVL)", y = "Adjusted scale parameter") +
  theme_science_polished

p_adj_scale_male
ggsave("Figures_alligators/fig_param_adj_scale_male.tiff", plot = p_adj_scale_male, dpi = 600,
       width = 7, height = 5, units = "in", compression = "lzw")
ggsave("Figures_alligators/fig_param_adj_scale_male.png", plot = p_adj_scale_male, dpi = 600,
       width = 7, height = 5, units = "in")

p_shape_male <- ggplot(data.frame(thresholds_male, shape_params, 
                                  shape_params_lower, shape_params_upper),
                       aes(x = thresholds_male, y = shape_params)) +
  geom_point() +
  geom_errorbar(aes(ymin = shape_params_lower, ymax = shape_params_upper), 
                width = 0.05, color = "blue") +
  geom_vline(xintercept = u_0, color = "red", linetype = "dashed") +
  geom_hline(yintercept = shape_params[which.min(abs(thresholds_male - u_0))], 
             color = "red", linetype = "dashed") +
  labs(x = "Threshold (log SVL)", y = "Shape parameter") +
  theme_science_polished

p_shape_male
ggsave("Figures_alligators/fig_param_shape_male.tiff", plot = p_shape_male, dpi = 600,
       width = 7, height = 5, units = "in", compression = "lzw")
ggsave("Figures_alligators/fig_param_shape_male.png", plot = p_shape_male, dpi = 600,
       width = 7, height = 5, units = "in")



### EVT SETUP ###
threshold <- log(165)

male_svl <- df_clean$log_SVL
male_ecdf <- ecdf(male_svl)
quantile_male <- male_ecdf(threshold)

# Fit GPD
gpd_fit <- fevd(df_clean$log_SVL, threshold = threshold, type = "GP")
scale_hat <- gpd_fit$results$par["scale"]
shape_hat <- gpd_fit$results$par["shape"]
shape_se <- summary(gpd_fit)$se["shape"]

# Wald z-test statistic for H0: ξ >= 0 vs H1: ξ < 0
z_value <- (shape_hat - 0) / shape_se
p_value <- pnorm(z_value)  # lower-tail test (ξ < 0)

cat("\n--- One-sided test for ξ < 0 ---\n")
cat("Estimated ξ:", round(shape_hat, 4), "\n")
cat("Standard Error:", round(shape_se, 4), "\n")
cat("Z statistic:", round(z_value, 4), "\n")
cat("One-sided p-value:", signif(p_value, 4), "\n")

# Exceedances
y_excess <- df_clean$log_SVL[df_clean$log_SVL > threshold]
excesses <- y_excess - threshold
n <- length(excesses)

# Sorted empirical values and probabilities
sorted_excesses <- sort(excesses)
probs <- ppoints(n)

# Theoretical quantiles for Q-Q plot
if (abs(shape_hat) > 1e-6) {
  theo_q <- threshold + scale_hat / shape_hat * ((probs^(-shape_hat)) - 1)
} else {
  theo_q <- threshold - scale_hat * log(probs)
}

qq_df <- data.frame(
  Theoretical = rev(theo_q),
  Empirical = sort(y_excess)
)

p_qq_gpd <- ggplot(qq_df, aes(x = Theoretical, y = Empirical)) +
  geom_point(color = "steelblue") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  labs(title = "Q-Q Plot for GPD Fit (Males)",
       x = "Theoretical Quantiles",
       y = "Empirical Quantiles") +
  theme_science_polished
p_qq_gpd

# Theoretical and empirical probabilities for P-P plot
F_empirical <- (1:n) / n
if (abs(shape_hat) > 1e-6) {
  F_theoretical <- 1 - (1 + shape_hat * excesses / scale_hat)^(-1 / shape_hat)
} else {
  F_theoretical <- 1 - exp(-excesses / scale_hat)
}

pp_df <- data.frame(
  Theoretical = sort(F_theoretical),
  Empirical = F_empirical
)

p_pp_gpd <- ggplot(pp_df, aes(x = Theoretical, y = Empirical)) +
  geom_point(color = "darkgreen") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  labs(
       x = "Theoretical CDF",
       y = "Empirical CDF") +
  theme_science_polished
p_pp_gpd

########################
# Kolmogorov–Smirnov Test
########################

# Step 1: Compute empirical CDF values under fitted GPD
F_hat <- pgpd(excesses, loc = 0, scale = scale_hat, shape = shape_hat)
ks_obs <- ks.test(F_hat, "punif")$statistic

# Step 2: Visual check – histogram with reference uniform line
p_ks <- ggplot(data.frame(F_hat = F_hat), aes(x = F_hat)) +
  geom_histogram(aes(y = ..density..), bins = 20, fill = "skyblue", color = "black", alpha = 0.7) +
  geom_hline(yintercept = 1, color = "red", linetype = "dashed") +
  geom_density(color = "darkblue", size = 1.2, adjust = 1.5) +
  labs(title = "Uniformity Check of Transformed Exceedances",
       x = expression(hat(F)(y)),
       y = "Density") +
  theme_science_polished
p_ks

# Step 3: Parametric bootstrap of KS statistic
B <- 1000
ks_boot <- numeric(B)

set.seed(42)
for (b in 1:B) {
  sim_data <- rgpd(n, loc = 0, scale = scale_hat, shape = shape_hat)
  fit_b <- tryCatch(gpd.fit(sim_data + threshold, threshold, show = FALSE), error = function(e) NULL)
  
  if (!is.null(fit_b)) {
    sigma_b <- fit_b$mle[1]
    xi_b <- fit_b$mle[2]
    F_b <- pgpd(sim_data, loc = 0, scale = sigma_b, shape = xi_b)
    ks_boot[b] <- ks.test(F_b, "punif")$statistic
  } else {
    ks_boot[b] <- NA
  }
}

# Step 4: Clean and compute bootstrap p-value
ks_boot <- na.omit(ks_boot)
p_value <- mean(ks_boot >= ks_obs)

# Output
cat("Observed KS statistic (male alligators):", ks_obs, "\n")
cat("Bootstrap p-value:", p_value, "\n")


############################
# Obtain posterior of y*
###########################
y_excess <- df$log_SVL[df$log_SVL > threshold]
n_excess <- length(y_excess)
Y_max <- max(y_excess)
eps <- max(1e-6, 1e-4 * Y_max)
y_grid <- seq(Y_max + eps, Y_max + 5, length.out = 1000)
xi_grid <- seq(-1, -5e-2, length.out = 1000)

loglik <- function(y_star, xi) {
  if (xi >= 0 || any(y_star <= y_excess)) return(-Inf)
  term1 <- n_excess * log(y_star - threshold) / xi
  term2 <- -(1 / xi + 1) * sum(log(y_star - y_excess))
  term3 <- -n_excess * log(-xi)
  return(term1 + term2 + term3)
}

marginal_log_post <- sapply(y_grid, function(y_star) {
  loglik_vals <- sapply(xi_grid, function(xi) loglik(y_star, xi))
  integrand <- exp(loglik_vals - max(loglik_vals))
  log(trapz(xi_grid, integrand)) + max(loglik_vals)
})

# Posterior on y*
prior <- ifelse(y_grid > Y_max, 1, 0)
posterior_unnorm <- exp(marginal_log_post - max(marginal_log_post)) * prior
posterior <- posterior_unnorm / trapz(y_grid, posterior_unnorm)
samples <- sample(y_grid, size = n_samples, replace = TRUE, prob = posterior)
map <- y_grid[which.max(posterior)]
hpdi <- hdi(samples, credMass = ci_level)

# Convert to cm
samples_cm <- exp(samples)
map_cm <- exp(map)
hpdi_cm <- exp(hpdi)

cat("\n--- SVL endpoint (cm) ---\n")
cat("MAP:", round(map_cm, 2), "\n")
cat("HPDI:", round(hpdi_cm[1], 2), "-", round(hpdi_cm[2], 2), "\n")

samples_df <- data.frame(SVL = samples_cm)

p_svl <- ggplot(samples_df, aes(x = SVL)) +
  geom_histogram(aes(y = ..density..), bins = 50, fill = "skyblue", color = "black", alpha = 0.6) +
  geom_density(color = "darkblue", size = 1.2) +
  geom_vline(xintercept = map_cm, linetype = "dashed", color = "purple", size = 1.2) +
  geom_vline(xintercept = hpdi_cm, linetype = "dotdash", color = "orange", size = 1.2) +
  geom_vline(xintercept = svl_stokes, color = "black", size = 1.0) +
  labs(x = "SVL endpoint (cm)", y = "Density") +
  scale_x_continuous(limits = c(200, 350)) +
  theme_science_polished

p_svl

### REGRESSION log_TL ~ log_SVL ###
fit <- lm(log_TL ~ log_SVL, data = df)
summary(fit)
a <- coef(fit)[1]; b <- coef(fit)[2]
covmat <- vcov(fit)
sigma_resid <- summary(fit)$sigma

### CONVOLVE TO TL ###
log_tl_grid <- seq(min(df$log_TL), max(df$log_TL) + 1, length.out = 500)
tl_grid <- exp(log_tl_grid)
posterior_tl <- numeric(length(tl_grid))

for (i in seq_along(tl_grid)) {
  log_tl <- log(tl_grid[i])
  integrand <- numeric(length(y_grid))
  
  for (j in seq_along(y_grid)) {
    mu <- a + b * y_grid[j]
    var <- sigma_resid^2 + covmat[1,1] + y_grid[j]^2 * covmat[2,2]
    integrand[j] <- dnorm(log_tl, mean = mu, sd = sqrt(var)) * posterior[j]
  }
  posterior_tl[i] <- trapz(y_grid, integrand)
}
posterior_tl <- posterior_tl / trapz(tl_grid, posterior_tl)

# Posterior summaries
tl_samples <- sample(tl_grid, size = n_samples, replace = TRUE, prob = posterior_tl)
map_tl <- tl_grid[which.max(posterior_tl)]
hpdi_tl <- hdi(tl_samples, credMass = ci_level)

### PLOT TL POSTERIOR ###
tl_df <- data.frame(TL = tl_grid, density = posterior_tl)
tl_samp_df <- data.frame(TL = tl_samples)

p_tl <- ggplot() +
  geom_histogram(data = tl_samp_df, aes(x = TL, y = ..density..),
                 bins = 40, fill = "lightgreen", color='black',alpha = 0.6) +
  geom_line(data = tl_df, aes(x = TL, y = density),
            color = "darkgreen", size = 1.2) +
  geom_vline(xintercept = map_tl, color = "purple", linetype = "dashed", size=1.2) +
  geom_vline(xintercept = hpdi_tl, color = "orange", linetype = "dotdash", size=1.2) +
  geom_vline(xintercept = tl_stokes, color = "black", size = 1.0) +
  labs(x = "TL endpoint (cm)", y = "Density") +
  xlim(350, 700) +
  theme_science_polished
p_tl


cat("\n--- TL endpoint (cm) ---\n")
cat("MAP:", round(map_tl, 2), "\n")
cat("HPDI:", round(hpdi_tl[1], 2), "-", round(hpdi_tl[2], 2), "\n")

# Add female data
df_females <- df %>% filter(Sex == "F", !is.na(SVL), !is.na(TL), Deform==0)
df_females$log_SVL <- log(df_females$SVL)
df_females$log_TL <- log(df_females$TL)
df_females$Sex <- "Female"
df_clean$Sex <- "Male"

df_both <- bind_rows(df_clean, df_females)

intercept_est <- coef(fit)[1]
slope_est <- coef(fit)[2]

# Generate a sequence of log(SVL) values for plotting the regression line
log_SVL_seq <- seq(min(df_clean$log_SVL), max(df_clean$log_SVL) + 0.3, length.out = 200)

# Use the model fit to get predicted log TL values
predicted_log_TL <- intercept_est + slope_est * log_SVL_seq

# Update regression line data for legend
line_data <- data.frame(
  log_SVL = log_SVL_seq,
  log_TL = predicted_log_TL,
  linetype = ifelse(log_SVL_seq <= max(df_clean$log_SVL), "solid", "dashed")
)

map_bayes <- log(map_cm)      
hpd_bayes <- log(hpdi_cm) 

log_TL_estimated <- log(map_tl)     # log of MAP TL endpoint
TL_HPDI <- hpdi_tl                  # TL HPDI bounds (already in cm)


p_loglog <- ggplot(df_both, aes(x = log_SVL, y = log_TL, color = Sex)) +
  geom_point(alpha = 0.5) +
  geom_line(data = subset(line_data, linetype == "solid"),
            aes(x = log_SVL, y = log_TL), inherit.aes = FALSE, color = "blue") +
  geom_line(data = subset(line_data, linetype == "dashed"),
            aes(x = log_SVL, y = log_TL), inherit.aes = FALSE,
            color = "blue", linetype = "dashed") +
  geom_errorbar(aes(x = map_bayes, ymin = log(TL_HPDI[1]), ymax = log(TL_HPDI[2])),
                inherit.aes = FALSE, width = 0.02, color = "darkgreen") +
  geom_errorbarh(aes(y = log_TL_estimated, xmin = hpd_bayes[1], xmax = hpd_bayes[2]),
                 inherit.aes = FALSE, height = 0.02, color = "darkgreen") +
  geom_point(aes(x = map_bayes, y = log_TL_estimated), inherit.aes = FALSE,
             color = "darkgreen", size = 2) +
  geom_point(aes(x = log(svl_stokes), y = log(tl_stokes)),
             inherit.aes = FALSE, color = "black", size = 2) +
  geom_vline(xintercept = threshold, linetype = "dashed", color = "black") +
  labs(
    x = "Log(Snout–Vent Length [cm])",
    y = "Log(Total Length [cm])",
    color = "Sex"
  ) +
  scale_color_manual(values = c("Male" = "blue", "Female" = "deeppink")) +
  coord_fixed() +
  theme_science_polished +
  theme(legend.position = "right")

p_loglog

### INPUTS ###
y <- log(svl_stokes)
u <- threshold
n_total <- nrow(df)
n_exceed <- length(y_excess)
P_exceed_u <- n_exceed / n_total

### INPUT ###
y <- log(svl_stokes)         # log(236)
u <- threshold               # log(187)
n_total <- nrow(df)
n_exceed <- length(y_excess)
P_exceed_u <- n_exceed / n_total  # empirical rate

### RESAMPLE FROM JOINT POSTERIOR ###
set.seed(42)
y_star_grid <- y_grid

posterior_matrix <- posterior / sum(posterior)  # normalize to sum to 1

# Draw 20,000 (y*, xi) pairs from the joint posterior
sample_indices <- sample(length(posterior_matrix), size = 20000, replace = TRUE, prob = as.vector(posterior_matrix))
sample_y_idx <- ((sample_indices - 1) %% length(y_star_grid)) + 1
sample_xi_idx <- ((sample_indices - 1) %/% length(y_star_grid)) + 1

sample_y_star <- y_star_grid[sample_y_idx]
sample_xi <- xi_grid[sample_xi_idx]

# Compute conditional exceedance probability
P_cond <- ((sample_y_star - y) / (sample_y_star - u))^(-1 / sample_xi)
P_cond[!is.finite(P_cond) | P_cond <= 0] <- NA

# Final posterior samples of exceedance probability
P_samples <- na.omit(P_exceed_u * P_cond)

### SUMMARY STATISTICS ###
# HPDI
hpdi_prob <- hdi(P_samples, credMass = 0.95)

# MAP using histogram mode (to avoid kernel smoothing artifacts)
dens <- density(P_samples, adjust = 1.1)
map_prob <- dens$x[which.max(dens$y)]

# Optional: Mean/median for reference
mean_prob <- mean(P_samples)
median_prob <- median(P_samples)

cat("Posterior exceedance probability P(Y > 236):\n")
cat("  MAP:", signif(map_prob, 5), "\n")
cat("  95% HPDI:", signif(hpdi_prob[1], 5), "-", signif(hpdi_prob[2], 5), "\n")
cat("  Mean:", signif(mean_prob, 5), "\n")
cat("  Median:", signif(median_prob, 5), "\n")

### PLOT POSTERIOR ###
p_exceed_df <- data.frame(P_exceed = P_samples)

p_exceed_plot <- ggplot(p_exceed_df, aes(x = P_exceed)) +
  geom_histogram(aes(y = ..density..), bins = 50, fill = "skyblue", color = "black", alpha = 0.6) +
  geom_density(color = "darkblue", size = 1.2, adjust = 1.1) +
  geom_vline(xintercept = map_prob, linetype = "dashed", color = "purple", size = 1.2) +
  geom_vline(xintercept = hpdi_prob, linetype = "dotdash", color = "orange", size = 1.2) +
  scale_x_log10(limits = c(1e-3, 0.05)) +
  labs(x = expression(P(Y > 236)), y = "Posterior density") +
  theme_science_polished
p_exceed_plot

output_dir <- "Figures_alligators"
if (!dir.exists(output_dir)) dir.create(output_dir)

plot_list <- list(
  p_adj_scale_male,                # threshold–stability: adjusted scale
  p_shape_male,                    # threshold–stability: shape
  p_qq_gpd,                        # Q-Q plot
  p_pp_gpd,                        # P-P plot
  p_ks,                            # uniformity check of F̂(y)
  p_svl,                           # posterior of SVL endpoint
  p_tl,                            # posterior of TL endpoint
  p_loglog,                        # log–log regression with endpoints
  p_exceed_plot                    # posterior P(Y > 236 cm)
)

file_names <- c(
  "param_adj_scale_male",
  "param_shape_male",
  "diagnostic_qq_alligator",
  "diagnostic_pp_alligator",
  "diagnostic_uniform_alligator",
  "posterior_svl",
  "posterior_tl",
  "loglog_regression",
  "posterior_exceedance_probability"
)

for (i in seq_along(plot_list)) {
  ggsave(
    filename = file.path(output_dir, paste0(file_names[i], ".png")),
    plot     = plot_list[[i]],
    dpi      = 600,
    width    = 7, height = 5, units = "in"
  )
  ggsave(
    filename = file.path(output_dir, paste0(file_names[i], ".pdf")),
    plot     = plot_list[[i]],
    device   = cairo_pdf,        # robust font embedding
    width    = 7, height = 5, units = "in"
  )
}


