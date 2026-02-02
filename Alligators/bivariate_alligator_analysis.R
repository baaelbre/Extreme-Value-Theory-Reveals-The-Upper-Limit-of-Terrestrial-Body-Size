library(ggplot2)
library(extRemes)
library(ismev)
library(numDeriv)
library(dplyr)
library(MASS)
library(readxl)
library(pracma)
library(HDInterval)
library(evd)
library(copula)

##############################
### LOAD DATA ###
##############################

df <- read_excel("Data/experimental_alligator_harvest_woodward.xlsx")

# Keep males only, remove NAs
df_males <- df %>% filter(Sex == "M", Deform==0,!is.na(SVL), !is.na(TL))
df_males$log_SVL <- log(df_males$SVL)
df_males$log_TL <- log(df_males$TL)
df_males$log_M <- log(df_males$WTkg)

log_SVL <- log(df_males$SVL)
log_TL <- log(df_males$TL)
log_M <- log(df_males$WTkg)

threshold_stability <- function(log_measurement, u_0, thresholds) {
  n <- length(thresholds)
  df <- data.frame(
    threshold = thresholds,
    scale = NA, shape = NA,
    scale_lower = NA, scale_upper = NA,
    shape_lower = NA, shape_upper = NA,
    adjusted_scale = NA,
    adjusted_scale_lower = NA, adjusted_scale_upper = NA
  )
  
  for (i in seq_along(thresholds)) {
    fit <- fevd(log_measurement, threshold = thresholds[i], type = "GP")
    df$scale[i] <- fit$results$par["scale"]
    df$shape[i] <- fit$results$par["shape"]
    df$adjusted_scale[i] <- df$scale[i] - df$shape[i] * (thresholds[i] - u_0)
    
    tryCatch({
      ci <- ci(fit, type = "parameter", alpha = 0.1)
      cov <- summary(fit)$cov.theta
      
      df$scale_lower[i] <- ci[1,1]
      df$scale_upper[i] <- ci[1,3]
      df$shape_lower[i] <- ci[2,1]
      df$shape_upper[i] <- ci[2,3]
      
      var_adj <- cov[1,1] + cov[2,2] * (thresholds[i] - u_0)^2
      df$adjusted_scale_lower[i] <- df$adjusted_scale[i] - qnorm(0.95) * sqrt(var_adj)
      df$adjusted_scale_upper[i] <- df$adjusted_scale[i] + qnorm(0.95) * sqrt(var_adj)
    }, error = function(e) {})
  }
  
  return(df)
}

# Apply to both SVL and TL
u_0_svl <- 4.9
thresholds_svl <- seq(u_0_svl, 5.3, length.out = 50)
df_svl <- threshold_stability(log_SVL, u_0_svl, thresholds_svl)

u_0_tl <- 5.7
thresholds_tl <- seq(u_0_tl, 5.9, length.out = 50)
df_tl <- threshold_stability(log_TL, u_0_tl, thresholds_tl)

# Choose threshold index for highlighting (e.g., from parameter stability)
highlight_svl <- 39
highlight_tl <- 42

# SVL: Adjusted Scale vs Threshold
ggplot(df_svl, aes(x = threshold, y = adjusted_scale)) +
  geom_point() +
  geom_errorbar(aes(ymin = adjusted_scale_lower, ymax = adjusted_scale_upper),
                width = 0.1, color = "blue") +
  geom_vline(xintercept = df_svl$threshold[highlight_svl], color = "red", linetype = "dashed") +
  geom_hline(yintercept = df_svl$adjusted_scale[highlight_svl], color = "red", linetype = "dashed") +
  labs(x = "Threshold (log SVL)", y = "Adjusted scale parameter") +
  theme_minimal()

# SVL: Shape vs Threshold
ggplot(df_svl, aes(x = threshold, y = shape)) +
  geom_point() +
  geom_errorbar(aes(ymin = shape_lower, ymax = shape_upper),
                width = 0.1, color = "blue") +
  geom_vline(xintercept = df_svl$threshold[highlight_svl], color = "red", linetype = "dashed") +
  geom_hline(yintercept = df_svl$shape[highlight_svl], color = "red", linetype = "dashed") +
  labs(x = "Threshold (log SVL)", y = "Shape parameter") +
  theme_minimal()

# TL: Adjusted Scale vs Threshold
ggplot(df_tl, aes(x = threshold, y = adjusted_scale)) +
  geom_point() +
  geom_errorbar(aes(ymin = adjusted_scale_lower, ymax = adjusted_scale_upper),
                width = 0.1, color = "blue") +
  geom_vline(xintercept = df_tl$threshold[highlight_tl], color = "red", linetype = "dashed") +
  geom_hline(yintercept = df_tl$adjusted_scale[highlight_tl], color = "red", linetype = "dashed") +
  labs(x = "Threshold (log TL)", y = "Adjusted scale parameter") +
  theme_minimal()

# TL: Shape vs Threshold
ggplot(df_tl, aes(x = threshold, y = shape)) +
  geom_point() +
  geom_errorbar(aes(ymin = shape_lower, ymax = shape_upper),
                width = 0.1, color = "blue") +
  geom_vline(xintercept = df_tl$threshold[highlight_tl], color = "red", linetype = "dashed") +
  geom_hline(yintercept = df_tl$shape[highlight_tl], color = "red", linetype = "dashed") +
  labs(x = "Threshold (log TL)", y = "Shape parameter") +
  theme_minimal()

# Define thresholds
threshold_log_SVL <- df_svl$threshold[highlight_svl]
threshold_log_TL  <- df_tl$threshold[highlight_tl]

# Classify points based on threshold exceedances
df_males <- df_males %>%
  mutate(
    exceedance = case_when(
      log_SVL >= threshold_log_SVL & log_TL >= threshold_log_TL ~ "Exceeds both",
      log_SVL >= threshold_log_SVL ~ "Exceeds SVL only",
      log_TL >= threshold_log_TL ~ "Exceeds TL only",
      TRUE ~ "Below both"
    )
  )

# Plot with color-coded regions
ggplot(df_males, aes(x = log_SVL, y = log_TL, color = exceedance)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_vline(xintercept = threshold_log_SVL, linetype = "dashed", color = "red") +
  geom_hline(yintercept = threshold_log_TL, linetype = "dashed", color = "red") +
  scale_color_manual(
    values = c(
      "Exceeds both" = "darkgreen",
      "Exceeds SVL only" = "blue",
      "Exceeds TL only" = "purple",
      "Below both" = "gray"
    )
  ) +
  labs(
    x = "log(SVL)",
    y = "log(TL)",
    color = "Exceedance status"
  ) +
  coord_fixed() +
  theme_minimal()



# Assume df_males is already filtered for deform == 0 and has log_SVL and log_TL columns

# Set thresholds (from stability analysis)
u_svl <- 5.21  # threshold for log(SVL)
u_tl  <- 5.86  # threshold for log(TL)

# Exceedances
exceedances <- df_males %>%
  filter(log_SVL > u_svl | log_TL > u_tl) %>%
  mutate(
    x = log_SVL - u_svl,
    y = log_TL - u_tl
  )

# Fit margins (GPD) for log SVL and log TL separately
fit_svl <- fpot(exceedances$x, threshold = 0, model = "gpd")
fit_tl  <- fpot(exceedances$y, threshold = 0, model = "gpd")

# Standardize marginals using fitted GPD CDFs (transform to unit Fréchet or uniform scale)
svl_pp <- pgpd(exceedances$x, loc = 0, scale = fit_svl$estimate[1], shape = fit_svl$estimate[2])
tl_pp  <- pgpd(exceedances$y, loc = 0, scale = fit_tl$estimate[1], shape = fit_tl$estimate[2])

# Convert to unit Fréchet margins
u_frechet <- -1 / log(svl_pp)
v_frechet <- -1 / log(tl_pp)

# Bivariate extreme value copula fitting (e.g., logistic model)
data_uv <- cbind(u_frechet, v_frechet)
fit_ev <- fgev(data_uv)  # Not appropriate directly; use logistic dependence structure:

# Estimate logistic model via censored likelihood
loglik_logistic <- function(alpha) {
  if (alpha <= 0 || alpha > 1) return(Inf)
  z <- (u_frechet^(-1/alpha) + v_frechet^(-1/alpha))^alpha
  log_lik <- -sum(log((1/alpha) * z^(1/alpha - 2) * (u_frechet^(-1/alpha - 1)) * (v_frechet^(-1/alpha - 1))))
  return(log_lik)
}

opt <- optimize(loglik_logistic, interval = c(0.01, 1))
alpha_hat <- opt$minimum

cat("Estimated logistic dependence parameter α:", round(alpha_hat, 3), "\n")

# Define bivariate logistic extreme value copula survival function
bivariate_survival <- function(u, v, alpha) {
  z <- (u^(-1/alpha) + v^(-1/alpha))^alpha
  return(exp(-z))
}

# Estimate joint survival function on a grid
u_seq <- seq(1.1, 10, length.out = 100)
v_seq <- seq(1.1, 10, length.out = 100)
grid <- expand.grid(u = u_seq, v = v_seq)

grid$surv <- bivariate_survival(grid$u, grid$v, alpha_hat)

# Plot joint survival surface (unit Fréchet scale)
library(plotly)

plot_ly(
  x = ~grid$u, y = ~grid$v, z = ~grid$surv,
  type = "scatter3d", mode = "markers",
  marker = list(size = 2, color = ~grid$surv, colorscale = "Viridis")
) %>%
  layout(
    scene = list(
      xaxis = list(title = "Unit Fréchet SVL"),
      yaxis = list(title = "Unit Fréchet TL"),
      zaxis = list(title = "Joint Survival")
    ),
    title = "Joint Survival Function under Bivariate GPD (Logistic)"
  )


##############################
### log(TL) ~ log(SVL) REGRESSION ###
##############################
# Filter only individuals with no deformities and valid SVL + TL
df_clean <- df %>%
  filter(Deform == 0, !is.na(SVL), !is.na(TL)) %>%
  mutate(log_SVL = log(SVL), log_TL = log(TL), log_M=log(WTkg))
fit_log <- lm(log_M ~ log_TL+log_SVL, data = df_clean)
summary(fit_log)
plot(fit_log)

