library(ggplot2)
library(readxl)
n_samples <-1000 
# Polished scientific plotting theme
theme_science_polished <- theme_minimal(base_family = "Arial", base_size = 12) +
  theme(
    axis.title = element_text(size = 16, face = "bold"),
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
# Data from campione and evans 2012 paper
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
df$log_sum_circ <- log(df$sum_circumference)
df$log_mass <- log(df[[mass_col]])
df$log_fem_circ <- log(df$`Femur Circumference`)

df$log_hum_circ <- log(df$`Humerus Circumference`)
# Fit linear and quadratic models
lin_model <- lm(log_mass ~ log_sum_circ, data = df)

# Compute squared term for the quadratic model
df$log_sum_circ_sq <- df$log_sum_circ^2

# Fit the quadratic model
quad_model <- lm(log_mass ~ log_sum_circ + log_sum_circ_sq, data = df)

# Predict fitted values across the range of the predictor
x_vals <- seq(min(df$log_sum_circ), max(df$log_sum_circ), length.out = 300)
pred_df <- data.frame(
  log_sum_circ = x_vals,
  log_sum_circ_sq = x_vals^2
)
pred_df$fit <- predict(quad_model, newdata = pred_df)

# Styled plot in same template as alligator regression
p_quad <- ggplot(df, aes(x = log_sum_circ, y = log_mass)) +
  geom_point(alpha = 0.5, color = "hotpink", size = 2) +
  geom_line(data = pred_df, aes(x = log_sum_circ, y = fit), 
            inherit.aes = FALSE, color = "blue", linewidth = 1.2) +
  labs(
    x = expression(log(C[F+H])~"[log mm]"),
    y = expression(log(Mass)~"[log g]"),
  ) +
  theme_science_polished  # Custom theme should already be defined

# Display and save
print(p_quad)
ggsave("Figures_sauropods/fig_mass_vs_circ_quadratic.png", plot = p_quad, width = 7, height = 5, dpi = 600)
ggsave("Figures_sauropods/fig_mass_vs_circ_quadratic.tiff", plot = p_quad, width = 7, height = 5, dpi = 600, compression = "lzw")

quad_model <- lm(log_mass ~ log_sum_circ + log_sum_circ_sq, data = df)

# Fit linear model: femur ~ humerus (log-log)
fem_hum_lm <- lm(log_fem_circ ~ log_hum_circ, data = df)

# Prepare newdata for prediction
hum_Arg <- 906
fem_Arg <- 1110
log_hum_arg <- log(hum_Arg)
log_fem_arg <- log(fem_Arg)
newdata <- data.frame(log_hum_circ = log_hum_arg)
pred <- predict(fem_hum_lm, newdata, interval = "prediction", level = 0.95)

# Prediction band over range
x_vals <- seq(min(df$log_hum_circ), log(906) + 0.1, length.out = 300)
pred_band <- predict(fem_hum_lm, newdata = data.frame(log_hum_circ = x_vals), interval = "prediction", level = 0.95)
pred_band_df <- data.frame(log_hum_circ = x_vals, pred_band)

# Prepare Argentinosaurus prediction and actual values
pred_df <- data.frame(log_hum_circ = log_hum_arg, pred)
arg_data <- data.frame(log_hum_circ = log_hum_arg, log_fem_circ = log_fem_arg)

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
# Final plot
p1 <- ggplot() +
  geom_point(data = df, aes(x = log_hum_circ, y = log_fem_circ), color = "gray30") +
  geom_line(data = pred_band_df, aes(x = log_hum_circ, y = fit), color = "blue", linewidth = 1) +
  geom_ribbon(data = pred_band_df, aes(x = log_hum_circ, ymin = lwr, ymax = upr), fill = "lightblue", alpha = 0.4) +
  geom_point(data = arg_data, aes(x = log_hum_circ, y = log_fem_circ), color = "darkred", size = 3) +
  labs(
    x = expression(log("Humerus circumference [mm]")),
    y = expression(log("Femur circumference [mm]")),
    title = "Femur vs Humerus Circumference (log scale)",
    subtitle = "Argentinosaurus: red = predicted; dark red = observed"
  ) +
  theme_science_polished
p1
ggsave("fig_femur_humerus_sauropods.png", plot=p1, width = 7, height = 5, dpi = 600)
ggsave("fig_femur_humerus_sauropods.tiff", plot=p1, width = 7, height = 5, dpi = 600, compression = "lzw")

## ================================================
## Top-10 CFH → 95% prediction intervals via lm()
##  - Uses lin_model and quad_model (or quad_model_centered if available)
##  - Exponentiate predictions at the end → tons
## ================================================

library(dplyr)
library(tibble)

# Your top 10 CFH (mm)
sum_circ_mm <- c(1639, 1672, 1687, 1692, 1694, 1695, 1704, 1827, 1870, 2016)
new_lin <- tibble(CFH_mm = sum_circ_mm,
                  log_sum_circ = log(as.double(CFH_mm)))

# ---- Linear model PI (requires: lin_model fit as lm(log_mass ~ log_sum_circ, data=df)) ----
pred_lin <- predict(lin_model, newdata = new_lin,
                    interval = "prediction", level = 0.95)

# ---- Quadratic model PI ----
quad_has_centered <- exists("quad_model_centered")

if (quad_has_centered) {
  # centered quadratic was fitted as: lm(log_mass ~ log_sum_circ_centered + log_sum_circ_centered_sq, data=df)
  c0 <- mean(df$log_sum_circ, na.rm = TRUE)
  new_cq <- tibble(
    log_sum_circ_centered    = new_lin$log_sum_circ - c0,
    log_sum_circ_centered_sq = (new_lin$log_sum_circ - c0)^2
  )
  pred_cq <- predict(quad_model_centered, newdata = new_cq,
                     interval = "prediction", level = 0.95)
} else {
  # non-centered quadratic was fitted as: lm(log_mass ~ log_sum_circ + log_sum_circ_sq, data=df)
  new_quad <- mutate(new_lin, log_sum_circ_sq = log_sum_circ^2)
  pred_cq  <- predict(quad_model, newdata = new_quad,
                      interval = "prediction", level = 0.95)
}

# ---- Exponentiate to tons and assemble tidy table ----
out <- tibble(
  CFH_mm               = new_lin$CFH_mm,
  Linear_pred_tons     = exp(pred_lin[,"fit"]) / 1e6,
  Linear_PI95_lo_t     = exp(pred_lin[,"lwr"]) / 1e6,
  Linear_PI95_hi_t     = exp(pred_lin[,"upr"]) / 1e6,
  Quadratic_pred_tons  = exp(pred_cq[,"fit"]) / 1e6,
  Quadratic_PI95_lo_t  = exp(pred_cq[,"lwr"]) / 1e6,
  Quadratic_PI95_hi_t  = exp(pred_cq[,"upr"]) / 1e6
) %>%
  mutate(across(-CFH_mm, ~round(., 1)))

print(out, n = nrow(out))

# Optional: write to CSV
# write.csv(out, "top10_CFH_mass_PI_lm_prediction.csv", row.names = FALSE)



###############
# Centered Quadratic Model Fitting and Saving
###############

# Center the predictor
log_sum_circ_mean <- mean(df$log_sum_circ)
df$log_sum_circ_centered <- df$log_sum_circ - log_sum_circ_mean
df$log_sum_circ_centered_sq <- df$log_sum_circ_centered^2

# Fit the centered quadratic model
quad_model_centered <- lm(log_mass ~ log_sum_circ_centered + log_sum_circ_centered_sq, data = df)

# Extract coefficients and standard errors
coef_vals <- coef(quad_model_centered)
coef_summary <- summary(quad_model_centered)$coefficients
vcov_mat <- vcov(quad_model_centered)
###########################
# Save Linear Model Coefficients
###########################

# Extract linear model parameters
lin_coef_vals <- coef(lin_model)
lin_coef_summary <- summary(lin_model)$coefficients
lin_vcov_mat <- vcov(lin_model)

lin_coefficients <- list(alpha=lin_alpha,beta=lin_beta,alpha_se=lin_alpha_se, beta_se=lin_beta_se, resid_sd=lin_sigma)

saveRDS(lin_coefficients, file = "linear_model_coefficients.rds")


# Prepare list of all relevant regression parameters
centered_quad_coefficients <- list(
  alpha    = coef_vals[1],
  beta     = coef_vals[2],
  gamma    = coef_vals[3],
  alpha_se = coef_summary[1, 2],
  beta_se  = coef_summary[2, 2],
  gamma_se = coef_summary[3, 2],
  cov_ab   = vcov_mat[1, 2],
  cov_ag   = vcov_mat[1, 3],
  cov_bg   = vcov_mat[2, 3],
  resid_sd = summary(quad_model_centered)$sigma,
  mean_log_sum_circ = log_sum_circ_mean
)

# Save to RDS
saveRDS(centered_quad_coefficients, file = "centered_quadratic_coefficients.rds")

# Center the predictor value for Argentinosaurus
x_centered <- x - log_sum_circ_mean

# Predict log mass
log_mass_hat_centered <- centered_quad_coefficients$alpha + 
  centered_quad_coefficients$beta * x_centered + 
  centered_quad_coefficients$gamma * x_centered^2

# Variance propagation for 95% prediction interval (add residual variance!)
var_log_mass_centered <- centered_quad_coefficients$alpha_se^2 + 
  (x_centered^2) * centered_quad_coefficients$beta_se^2 + 
  (x_centered^4) * centered_quad_coefficients$gamma_se^2 +
  centered_quad_coefficients$resid_sd^2

# Compute interval and mass on original scale
z <- qnorm(0.975)
ci_lower_centered <- log_mass_hat_centered - z * sqrt(var_log_mass_centered)
ci_upper_centered <- log_mass_hat_centered + z * sqrt(var_log_mass_centered)

mass_hat_centered <- exp(log_mass_hat_centered) / 1e6
mass_ci_lower_centered <- exp(ci_lower_centered) / 1e6
mass_ci_upper_centered <- exp(ci_upper_centered) / 1e6

cat("---- Centered quadratic (prediction interval) ----\n")
cat("Log mass estimate:", log_mass_hat_centered, "\n")
cat("95% PI for log mass: [", ci_lower_centered, ",", ci_upper_centered, "]\n")
cat("Mass estimate (tons):", mass_hat_centered, "\n")
cat("95% PI for mass (tons): [", mass_ci_lower_centered, ",", mass_ci_upper_centered, "]\n")

# Plot log-mass densities for Linear and Centered Quadratic
df_logmass_centered <- data.frame(
  model = c("Linear", "Centered quadratic"),
  mean = c(log_mass_hat_lin, log_mass_hat_centered),
  sd = c(sqrt(var_log_mass_lin), sqrt(var_log_mass_centered))
)

# Build data frame to plot both normal densities manually
x_vals <- seq(10, 20, length.out = 1000)

df_density <- data.frame(
  log_mass = rep(x_vals, 2),
  model = rep(c("Linear", "Centered quadratic"), each = length(x_vals)),
  density = c(
    dnorm(x_vals, mean = log_mass_hat_lin, sd = sqrt(var_log_mass_lin)),
    dnorm(x_vals, mean = log_mass_hat_centered, sd = sqrt(var_log_mass_centered))
  )
)

# Plot the two predictive log-mass densities correctly
p3 <- ggplot(df_density, aes(x = log_mass, y = density, color = model)) +
  geom_line(size = 1.2) +
  geom_vline(xintercept = c(ci_lower_lin, ci_upper_lin), linetype = "dashed", color = "blue") +
  geom_vline(xintercept = c(ci_lower_centered, ci_upper_centered), linetype = "dashed", color = "darkgreen") +
  scale_color_manual(values = c("Linear" = "blue", "Centered quadratic" = "darkgreen")) +
  labs(
    x = expression(log("Mass [g]")),
    y = "Normal density"
  ) +
  theme_science_polished
p3
ggsave("fig_linear_vs_centered_quadratic.png", plot=p3, width = 7, height = 5, dpi = 600)
ggsave("fig_linear_vs_centered_quadratic.tiff", plot=p3, width = 7, height = 5, dpi = 600, compression = "lzw")


#####################################
# Compare with mPPE
#####################################

circ_mm <- c(200,400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2016)
log_circ <- log(circ_mm)
x_centered_range <- log_circ - log_sum_circ_mean

# Predict log-mass and variance (with residual variance)
log_mass_hat_range <- centered_quad_coefficients$alpha + 
  centered_quad_coefficients$beta * x_centered_range + 
  centered_quad_coefficients$gamma * x_centered_range^2

var_log_mass_range <- centered_quad_coefficients$alpha_se^2 + 
  (x_centered_range^2) * centered_quad_coefficients$beta_se^2 + 
  (x_centered_range^4) * centered_quad_coefficients$gamma_se^2 + 
  centered_quad_coefficients$resid_sd^2  # include prediction variance

# z-values
z_95 <- qnorm(0.975)
z_40 <- qnorm(0.70)

# 95% PI
ci_lower_log_95 <- log_mass_hat_range - z_95 * sqrt(var_log_mass_range)
ci_upper_log_95 <- log_mass_hat_range + z_95 * sqrt(var_log_mass_range)

# 40% PI (for comparison with mPPE)
ci_lower_log_40 <- log_mass_hat_range - z_40 * sqrt(var_log_mass_range)
ci_upper_log_40 <- log_mass_hat_range + z_40 * sqrt(var_log_mass_range)

# Convert to mass (tons)
mass_est <- exp(log_mass_hat_range) / 1e6
ci_lower_mass_95 <- exp(ci_lower_log_95) / 1e6
ci_upper_mass_95 <- exp(ci_upper_log_95) / 1e6
ci_lower_mass_40 <- exp(ci_lower_log_40) / 1e6
ci_upper_mass_40 <- exp(ci_upper_log_40) / 1e6

# mPPE heuristic
mass_fixed_lower <- mass_est * 0.75
mass_fixed_upper <- mass_est * 1.25

# Combine into data frame
compare_df <- data.frame(
  Circumference = circ_mm,
  Mass_Estimate = mass_est,
  CI_Lower_40 = ci_lower_mass_40,
  CI_Upper_40 = ci_upper_mass_40,
  CI_Lower_95 = ci_lower_mass_95,
  CI_Upper_95 = ci_upper_mass_95,
  CI_Lower_25pct = mass_fixed_lower,
  CI_Upper_25pct = mass_fixed_upper
)

# Plot
p4 <- ggplot(compare_df, aes(x = Circumference)) +
  geom_line(aes(y = CI_Lower_95, color = "Model PI (95%)"), linewidth = 1.2) +
  geom_line(aes(y = CI_Upper_95, color = "Model PI (95%)"), linewidth = 1.2) +
  geom_line(aes(y = CI_Lower_40, color = "Model PI (40%)"), linewidth = 1.2, linetype = "twodash") +
  geom_line(aes(y = CI_Upper_40, color = "Model PI (40%)"), linewidth = 1.2, linetype = "twodash") +
  geom_line(aes(y = CI_Lower_25pct, color = "mPPE"), linetype = "dashed", linewidth = 1.2) +
  geom_line(aes(y = CI_Upper_25pct, color = "mPPE"), linetype = "dashed", linewidth = 1.2) +
  geom_point(aes(y = Mass_Estimate), color = "black", shape = 21, fill = "white", size = 3) +
  scale_color_manual(values = c(
    "Model PI (95%)" = "blue",
    "Model PI (40%)" = "darkgreen",
    "mPPE" = "red"
  )) +
  labs(
    x = "Femur + Humerus Circumference (mm)",
    y = "Mass Estimate [tons]",
    color = "Interval Type"
  ) +
  theme_science_polished
p4
ggsave("fig_PI_vs_mPPE.png", plot=p4, width = 7, height = 5, dpi = 600)
ggsave("fig_PI_vs_mPPE.tiff", plot=p4, width = 7, height = 5, dpi = 600, compression = "lzw")


# -------------------------------
# Plot corresponding confidence level for ±25% heuristic
# -------------------------------
# Compute z and confidence level required to match ±25% rule
z_target <- log(1.25) / sqrt(var_log_mass_range)
alpha_range <- 2 * (1 - pnorm(z_target))
conf_level_range <- 1 - alpha_range

# Create data frame
conf_vs_circ_df <- data.frame(
  Circumference = circ_mm,
  Confidence_Level = conf_level_range
)

# required confidence level to match the rule
# Plot
p5 <- ggplot(conf_vs_circ_df, aes(x = Circumference, y = Confidence_Level)) +
  geom_line(color = "steelblue", linewidth = 1.3) +
  geom_point(color = "black", size = 2) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0.3, 0.6)) +
  labs(
    x = "Femur + Humerus Circumference (mm)",
    y = "Confidence Level"
  ) +
  theme_science_polished
p5
ggsave("Figures_sauropods/fig_required_confidence_level.png", plot=p5, width = 7, height = 5, dpi = 600)
ggsave("FIgures_sauropods/fig_required_confidence_level.tiff", plot=p5, width = 7, height = 5, dpi = 600, compression = "lzw")









