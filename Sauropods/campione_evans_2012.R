# Load libraries
library(readxl)

# Load the dataset
df <- read_excel("Data/campione_evans_2012.xlsx", sheet = "Table 1.csv")
df <- na.omit(df)

# Compute sum of femur + humerus circumferences
df$sum_circumference <- df$`Femur Circumference` + df$`Humerus Circumference`
df$log_sum_circ <- log10(df$sum_circumference)
df$log_mass <- log10(df$`Body Mass (g)`)

#############################
# Linear Model
#############################

lin_model <- lm(log_mass ~ log_sum_circ, data = df)

lin_summary <- summary(lin_model)
lin_coeffs <- lin_summary$coefficients
lin_cov <- vcov(lin_model)
lin_sigma2 <- lin_summary$sigma^2  # residual variance

lin_coef_list <- list(
  alpha = lin_coeffs[1, 1],
  beta = lin_coeffs[2, 1],
  alpha_se = lin_coeffs[1, 2],
  beta_se = lin_coeffs[2, 2],
  cov = lin_cov,
  resid_var = lin_sigma2
)

saveRDS(lin_model, "linear_model.rds")
saveRDS(lin_coef_list, "linear_coefficients.rds")

#############################
# Centered Quadratic Model
#############################

# Center the predictor
mean_log_circ <- mean(df$log_sum_circ)
df$log_sum_circ_centered <- df$log_sum_circ - mean_log_circ
df$log_sum_circ_centered_sq <- df$log_sum_circ_centered^2

quad_model <- lm(log_mass ~ log_sum_circ_centered + log_sum_circ_centered_sq, data = df)

quad_summary <- summary(quad_model)
quad_coeffs <- quad_summary$coefficients
quad_cov <- vcov(quad_model)
quad_sigma2 <- quad_summary$sigma^2  # residual variance

quad_coef_list <- list(
  alpha = quad_coeffs[1, 1],
  beta = quad_coeffs[2, 1],
  gamma = quad_coeffs[3, 1],
  alpha_se = quad_coeffs[1, 2],
  beta_se = quad_coeffs[2, 2],
  gamma_se = quad_coeffs[3, 2],
  cov = quad_cov,
  mean_log_sum_circ = mean_log_circ,
  resid_var = quad_sigma2
)

saveRDS(quad_model, "centered_quadratic_model.rds")
saveRDS(quad_coef_list, "centered_quadratic_coefficients.rds")
