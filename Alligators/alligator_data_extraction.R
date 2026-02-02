library(imager)
library(digitize)
library(ggplot2)

# Load the PNG image
image_path <- "males_2.png"  # Change this to your image path
img <- load.image(image_path)

# Open a new window to ensure correct resolution (Optional)
dev.new(width = 8, height = 8)  # Adjust size if needed

# Display the image with correct aspect ratio
par(mar=c(0,0,0,0))  # Remove margins
plot(img, main = "Original Image", asp = 1)

# Step 1: Calibrate the axes (Select known reference points)
data_points <- digitize(image_path)

# Print extracted data points
print(data_points)

# Overlay the extracted points on the image
points(data_points, col = "red", pch = 19, cex = 1.5)  # Plot extracted points

# Save manually extracted data
write.csv(data_points, "manual_extracted_data.csv", row.names = FALSE)

# Close the plotting window
dev.off()
# Assuming `data_points` is already a dataframe with columns "x" and "y"
ggplot(data_points, aes(x = x, y = y)) + 
  geom_point(color = "blue", size = 1) +  # Blue points for better visibility
  labs(
    x = "Head Length", 
    y = "Total Length", 
    title = "Scatter Plot for Male Alligators Measured in Florida (1987-1992)"
  ) +
  theme_minimal()  # Clean theme
# Regular Histogram of Total Length
ggplot(data_points, aes(x = y)) + 
  geom_histogram(binwidth = 10, fill = "blue", color = "black", alpha = 0.7) +
  labs(
    x = "Total Length", 
    y = "Count", 
    title = "Histogram of Total Length (Male Alligators, FL 1987-1992)"
  ) +
  theme_minimal()

# Log-transformed Histogram of Total Length
ggplot(data_points, aes(x = log(y))) + 
  geom_histogram(binwidth = 0.1, fill = "red", color = "black", alpha = 0.7) +
  labs(
    x = "Log(Total Length)", 
    y = "Count", 
    title = "Histogram of Log-Transformed Total Length (Male Alligators, FL 1987-1992)"
  ) +
  theme_minimal() 
