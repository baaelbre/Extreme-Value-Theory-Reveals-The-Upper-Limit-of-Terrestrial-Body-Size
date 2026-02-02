import cv2
import numpy as np
import matplotlib.pyplot as plt

# Load the image
image_path = "males.png"  # Change to "males.png" if needed
img = cv2.imread(image_path, cv2.IMREAD_GRAYSCALE)

# Apply adaptive thresholding to enhance contrast
thresh = cv2.adaptiveThreshold(img, 255, cv2.ADAPTIVE_THRESH_GAUSSIAN_C, 
                               cv2.THRESH_BINARY_INV, 15, 4)

# Morphological operations to remove small noise
kernel = np.ones((3,3), np.uint8)
morph = cv2.morphologyEx(thresh, cv2.MORPH_OPEN, kernel, iterations=1)

# Detect contours
contours, _ = cv2.findContours(morph, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)

# Set up plot and image dimensions
h, w = img.shape
filtered_contours = []

# Region filtering to exclude text and axes (assumes text is near the top)
for cnt in contours:
    x, y, w_cnt, h_cnt = cv2.boundingRect(cnt)
    
    # Heuristic to remove text, large objects, and elements near axes
    if 10 < w_cnt < 50 and 10 < h_cnt < 50:  # Filtering based on reasonable data point size
        if y > 50 and x > 50:  # Avoid labels near top-left
            filtered_contours.append(cnt)

# Create an output image with detected points
output = cv2.cvtColor(img, cv2.COLOR_GRAY2BGR)
cv2.drawContours(output, filtered_contours, -1, (0, 0, 255), 2)  # Draw detected points in red

# Show results
plt.figure(figsize=(6,6))
plt.imshow(output, cmap='gray')
plt.title("Detected Data Points")
plt.axis("off")
plt.show()
