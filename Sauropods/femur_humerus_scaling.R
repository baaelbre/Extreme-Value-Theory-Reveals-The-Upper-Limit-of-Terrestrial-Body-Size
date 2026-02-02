library(readxl)
library(dplyr)
library(ggplot2)

# Load sauropod data
saur_df <- read_excel("Data/DEMic23_updated_Supplemental_Data_withPubYear.xlsx", 
                      sheet = "humCirc~femCirc")

# Keep only rows with both femur and humerus circumference
saur_df <- saur_df %>%
  filter(!is.na(`femur circ (mm)`), !is.na(`humerus circ (mm)`)) %>%
  mutate(
    log_fem_circ = log10(`femur circ (mm)`),
    log_hum_circ = log10(`humerus circ (mm)`),
    taxon = `genus and species`
  )

# Fit linear model on sauropods only
saur_lm <- lm(log_fem_circ ~ log_hum_circ, data = saur_df)

# Add Argentinosaurus manually
arg_df <- data.frame(
  taxon = "Argentinosaurus",
  humerus = 906,
  femur = 1110,
  log_hum_circ = log10(906),
  log_fem_circ = log10(1110)
)

# Generate prediction band from sauropod model
x_vals <- seq(min(saur_df$log_hum_circ), max(c(saur_df$log_hum_circ, arg_df$log_hum_circ)) + 0.1, length.out = 300)
pred_band <- predict(saur_lm, 
                     newdata = data.frame(log_hum_circ = x_vals), 
                     interval = "prediction", level = 0.95)
pred_band_df <- data.frame(log_hum_circ = x_vals, pred_band)

# Combine and plot
ggplot() +
  geom_point(data = saur_df, aes(x = log_hum_circ, y = log_fem_circ), color = "darkgreen", size = 2.5) +
  geom_point(data = arg_df, aes(x = log_hum_circ, y = log_fem_circ), 
             color = "darkred", fill = "red", size = 4, shape = 21, stroke = 1.4) +
  geom_line(data = pred_band_df, aes(x = log_hum_circ, y = fit), color = "blue", linewidth = 1) +
  geom_ribbon(data = pred_band_df, aes(x = log_hum_circ, ymin = lwr, ymax = upr), 
              fill = "lightblue", alpha = 0.4) +
  labs(
    x = expression(log[10]("Humerus circumference [mm]")),
    y = expression(log[10]("Femur circumference [mm]")),
    title = "Sauropod Femur vs Humerus Circumference (log scale)",
    subtitle = "Red = Argentinosaurus (manually added)"
  ) +
  theme_minimal(base_size = 14)

# Predict at Argentinosaurus log humerus
arg_pred <- predict(saur_lm, 
                    newdata = data.frame(log_hum_circ = arg_df$log_hum_circ),
                    interval = "prediction", level = 0.95)

cat("Argentinosaurus observed femur (log):", arg_df$log_fem_circ, "\n")
cat("Predicted 95% PI:", arg_pred[1, "lwr"], "-", arg_pred[1, "upr"], "\n")

if (arg_df$log_fem_circ >= arg_pred[1, "lwr"] && arg_df$log_fem_circ <= arg_pred[1, "upr"]) {
  cat("✅ Argentinosaurus lies within the 95% prediction interval.\n")
} else {
  cat("❌ Argentinosaurus lies outside the 95% prediction interval.\n")
}
