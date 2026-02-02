##############################################################################
# Author: Bastiaan A. Van Velthoven
# Date: 18/01/2026
# Title: "Extreme value theory reveals the upper limit of terrestrial body size"
# Description: Code to generate time seris of sauropod discoveries and histograms
# using ggplot2
##############################################################################

library(readxl) # excel
library(dplyr) # dataframe manipulations
library(ggplot2) # plotting
library(scales) # for coloring the top 10 sauropods and making nice x-ticks
library(grid) # plot multiple plots onto a figure.

set.seed(42)

# ---------------------------
# Config
# ---------------------------
# point to data file
data_xlsx <- "Data/DEmic23_updated_Supplemental_Data_withPubYear update Nov25_2025.xlsx"
# save figures to ...
fig_dir   <- "Figure1" 
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

top_n                 <- 10L          # number of top taxa you want to use
use_all_masses        <- TRUE         # histogram uses all masses (even though the year may be missing)
out <- "Figure1" # name of the figure -> + _***.pdf

# ---------------------------
# General stylistic stuff
# ---------------------------
theme_science <- theme_minimal(base_family = "Arial", base_size = 12) +
  theme(
    axis.title       = element_text(size = 14, face = "bold"),
    axis.text        = element_text(size = 12),
    legend.title     = element_text(size = 10, face = "bold"),
    legend.text      = element_text(size = 10),
    panel.grid.major = element_line(color = "gray80", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.6),
    axis.ticks.length= unit(0.20, "cm"),
    axis.ticks       = element_line(color = "black", linewidth = 0.4),
    plot.margin      = margin(5, 5, 5, 5),
    legend.position  = "right"
  )

# ---------------------------
# Load data
# ---------------------------
df <- read_excel(data_xlsx)

df <- df %>%
  transmute(
    species = as.character(`genus and species`),
    year    = as.integer(`publication year`),
    m_kg    = as.numeric(`Campione est mass quadratic (kg)`),
    lo_kg   = as.numeric(`quadratic mass est low (kg)`),
    hi_kg   = as.numeric(`quadratic mass est high (kg)`)
  ) %>%
  mutate(
    m_t  = m_kg / 1000, # to tons 
    lo_t = lo_kg / 1000,
    hi_t = hi_kg / 1000
  )

# Time plot needs publication years
df_plot <- df %>% filter(is.finite(year))

# Histogram data: either all masses or only those for which year is available
df_hist <- if (use_all_masses) {
  df 
} else {
  df_plot
}

# y limits on linear and log scale
mass_range_lin <- c(0, 100)
mass_range_log <- c(0.1, 100)

# ---------------------------
# Top N taxa selection
# ---------------------------
top_species <- df_plot$species[order(df_plot$m_t, decreasing = TRUE)][seq_len(top_n)]

### coloring + legend ###
# if species in top 10, give it a color according to the scales::hue_pal, otherwise make it (light)grey
# relevel the factor
df_plot <- df_plot %>%
  mutate(
    highlight = factor(
      if_else(species %in% top_species, species, "Other"),
      levels = c(top_species, "Other")
    )
  )

# hue_pal is a color palette generator -> make a palet (list)
top_colors <- hue_pal()(top_n)
names(top_colors) <- top_species
col_map <- c(Other = "grey75", top_colors)


# ---------------------------
# Scatter plot (linear y)
# ---------------------------
p_time <- ggplot(df_plot, aes(x = year, y = m_t)) +
  geom_linerange(
    data = df_plot,
    aes(ymin = lo_t, ymax = hi_t, colour = highlight),
    linewidth = 0.55, alpha = 0.90,
    show.legend = FALSE
  ) +
  geom_point(
    aes(colour = highlight),
    size = 2.6, alpha = 0.95
  ) +
  scale_colour_manual(
    values = col_map,
    breaks = top_species,
    name   = paste("Top", length(top_species))
  ) +
  # pretty breaks from the scales package make sure it prints 1900,1920 etc instead of 1903, 1923
  
  scale_x_continuous(breaks = pretty_breaks(n = 10)) +
  scale_y_continuous(limits = mass_range_lin) +
  labs(x = "Publication year", y = "Mass estimate (tons)") +
  theme_science # +theme(legend.position = "none") # if you want the legend removed

print(p_time)
ggsave(filename = paste(fig_dir,"mass_time_linear.pdf", sep="/"),  plot = p_time, device = cairo_pdf,
       width = 7.2, height = 4.8, units = "in")

# ---------------------------
# Scatter plot (log y)
# ---------------------------
p_time_log <- ggplot(df_plot, aes(x = year, y = m_t)) +
  geom_linerange(
    data = df_plot,
    aes(ymin = lo_t, ymax = hi_t, colour = highlight), # with the color map we made
    linewidth = 0.55, alpha = 0.90,
    show.legend = FALSE
  ) +
  geom_point(
    aes(colour = highlight),
    size = 2.6, alpha = 0.95
  ) +
  scale_colour_manual(
    values = col_map, # here is the colour map
    breaks = top_species,
    name   = paste("Top", length(top_species))
  ) +
  scale_x_continuous(breaks = pretty_breaks(n = 10)) + 
  
  scale_y_log10(limits = mass_range_log) +
  labs(x = "Publication year", y = "Mass estimate (tons, log scale)") +
  theme_science +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1)))
  # + theme(legend.position = "none") # if you want a legend

print(p_time_log)
ggsave(filename=paste(fig_dir,"mass_time_log.pdf", sep="/"), plot = p_time_log, device = cairo_pdf,
       width = 7.2, height = 4.8)

# ============================================================
# Histograms (linear and log)
# ============================================================

p_hist <- ggplot(df_hist, aes(x = m_t)) +
  geom_histogram(bins = 30, color = "black", fill = "grey80") +
  scale_x_continuous(limits = mass_range_lin) +
  coord_flip() + # coordinate flip because I want them right of the time series
  labs(y = "Count", x = "Mass estimate (tons)") +
  theme_science +
  theme(legend.position = "none")
print(p_hist)
ggsave(paste(fig_dir,"hist.pdf", sep="/"), plot = p_hist, device = cairo_pdf,
       width = 3.0, height = 4.8)

p_hist_log <- ggplot(df_hist, aes(x = m_t)) +
  geom_histogram(bins = 30, color = "black", fill = "grey80") +
  scale_x_log10(limits = mass_range_log) +
  coord_flip() +
  labs(y = "Count", x = "Mass estimate (tons, log scale)") +
  theme_science +
  theme(legend.position = "none")
print(p_hist_log)
ggsave(paste(fig_dir,"hist_log.pdf", sep="/"), plot = p_hist_log, device = cairo_pdf,
       width = 3.0, height = 4.8)

# ============================================================
# Combo plots (linear andlog)
# ============================================================

# remove y label (mass label already on scatter), and top 10 legend
p_hist <- p_hist + labs(x=NULL, y = NULL)  + theme(legend.position='None')
p_hist_log <- p_hist_log + labs(x=NULL, y = NULL)  + theme(legend.position='None')
p_time <- p_time  + theme(legend.position='None')
p_time_log <- p_time_log  + theme(legend.position='None')

### now save p_time and p_hist into one figure with a shared y-axis ###
grid.newpage() # create blank canvas

# layout: 1 row, 2 columns with relative "null" proportions (3:1)
# check ?grid.layout 
lay <- grid.layout(
  nrow = 1, ncol = 2,
  widths = grid::unit(c(3, 1), "null")
)

# start drawing on it
pushViewport(viewport(layout=lay))
print(p_time, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(p_hist, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
ggsave(paste(fig_dir,"time_hist.pdf", sep="/"), grid.grab(),  # save plot from current grid
       device = cairo_pdf, width = 7.2, height = 4.8)

# for the log scale figures
grid.newpage() 
lay <- grid.layout(
  nrow = 1, ncol = 2,
  widths = grid::unit(c(3, 1), "null")
)
pushViewport(viewport(layout=lay))
print(p_time_log, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(p_hist_log, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))

ggsave(paste(fig_dir,"time_hist_log.pdf", sep="/"),   plot = grid.grab(),
       device = cairo_pdf, width = 7.2, height = 4.8)

