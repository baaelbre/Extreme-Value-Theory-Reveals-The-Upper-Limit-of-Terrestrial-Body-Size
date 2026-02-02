## ------------------------------------------------------------
## 0. Packages & data import
## ------------------------------------------------------------
library(readr)
library(dplyr)
library(ggplot2)

# Read data from Data/ folder
gators <- read_csv("Data/alligator_three_lakes_harvest.csv")

str(gators)
# columns: year, lake, total_harvest, avg_length_ft, avg_length_in

## ------------------------------------------------------------
## 1. Time series of harvests per lake
## ------------------------------------------------------------

p_by_lake <- ggplot(gators,
                    aes(x = year, y = total_harvest, colour = lake)) +
  geom_line(linewidth = 1) +
  geom_point(size = 1.8) +
  labs(x = "Year",
       y = "Total harvest",
       colour = "Lake",
       title = "Annual alligator harvest per lake") +
  scale_x_continuous(breaks = pretty(gators$year)) +
  theme_minimal(base_size = 12)

print(p_by_lake)
# ggsave("Figures/harvest_by_lake.png", p_by_lake, width = 7, height = 4.5, dpi = 300)

## ------------------------------------------------------------
## 2. Time series of summed harvests (3 lakes combined)
## ------------------------------------------------------------

gators_total <- gators %>%
  group_by(year) %>%
  summarise(total_harvest = sum(total_harvest), .groups = "drop")

p_total <- ggplot(gators_total,
                  aes(x = year, y = total_harvest)) +
  geom_line(linewidth = 1, colour = "steelblue") +
  geom_point(size = 1.8, colour = "steelblue") +
  labs(x = "Year",
       y = "Total harvest (3 lakes combined)",
       title = "Annual alligator harvest in Lochloosa, Newnans & Orange Lakes") +
  scale_x_continuous(breaks = pretty(gators_total$year)) +
  theme_minimal(base_size = 12)

print(p_total)
# ggsave("Figures/harvest_total_3lakes.png", p_total, width = 7, height = 4.5, dpi = 300)

## ------------------------------------------------------------
## 3. Mean harvest from 2007 onwards (for return-period scaling)
## ------------------------------------------------------------

gators_total_2007 <- gators_total %>%
  filter(year >= 2007)

mean_harvest_2007plus <- mean(gators_total_2007$total_harvest)

mean_harvest_2007plus
# This is the N_year you can plug into your return-period calculation.

# Optional: show this mean as a horizontal line on the plot
p_total_2007 <- ggplot(gators_total_2007,
                       aes(x = year, y = total_harvest)) +
  geom_line(linewidth = 1, colour = "darkgreen") +
  geom_point(size = 1.8, colour = "darkgreen") +
  geom_hline(yintercept = mean_harvest_2007plus,
             linetype = "dashed", colour = "red") +
  annotate("text",
           x = min(gators_total_2007$year),
           y = mean_harvest_2007plus,
           vjust = -0.5, hjust = 0,
           label = paste0("Mean (≥2007): ",
                          round(mean_harvest_2007plus, 1)),
           size = 3.4) +
  labs(x = "Year",
       y = "Total harvest (3 lakes combined)",
       title = "Total harvest since 2007 with mean level") +
  scale_x_continuous(breaks = pretty(gators_total_2007$year)) +
  theme_minimal(base_size = 12)

print(p_total_2007)
# ggsave("Figures/harvest_total_3lakes_since2007.png", p_total_2007,
#        width = 7, height = 4.5, dpi = 300)
