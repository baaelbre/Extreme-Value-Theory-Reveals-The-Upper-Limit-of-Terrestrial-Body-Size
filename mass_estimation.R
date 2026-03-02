##############################################################################
# Author: Bastiaan A. Van Velthoven
# Date: 18/01/2026
# Title: "Extreme value theory reveals the upper limit of terrestrial body size"
##############################################################################

library(readxl)
library(dplyr)
library(ggplot2)
library(scales)
library(zoo)

set.seed(42)

data_xlsx <- "Data/DEmic23_updated_Supplemental_Data_withPubYear update Nov25_2025.xlsx"
fig_dir   <- "Fig1"
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

top_n <- 10L

theme_nature <- theme_minimal(base_family = "Arial", base_size = 12) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 13, face = "bold"),
    axis.text  = element_text(size = 11),
    axis.ticks = element_line(linewidth = 0.4),
    axis.ticks.length = unit(2.5, "mm"),
    panel.border = element_rect(fill = NA, linewidth = 0.6),
    plot.title  = element_blank(),
    plot.margin = margin(8, 8, 8, 8),
    legend.text = element_text(size=9),
    legend.position = c(0.03, 0.03),
    legend.justification = c(0, 0),
    legend.direction = "horizontal",
    legend.background = element_rect(fill = alpha("white", 0.5), colour = NA),
    legend.key = element_rect(fill = NA, colour = NA),
    legend.key.size = unit(3.0, "mm"),       
    legend.key.width = unit(4.0, "mm"),
    legend.key.height = unit(3.0, "mm"),
    legend.spacing.x = unit(1.0, "mm"),
    legend.spacing.y = unit(0.5, "mm"),

  )

df <- read_excel(data_xlsx) %>%
  transmute(
    species = as.character(`genus and species`),
    year    = as.integer(`publication year`),
    m_kg    = as.numeric(`Campione est mass quadratic (kg)`),
    lo_kg   = as.numeric(`quadratic mass est low (kg)`),
    hi_kg   = as.numeric(`quadratic mass est high (kg)`)
  ) %>%
  mutate(
    m_t  = m_kg  / 1000,
    lo_t = lo_kg / 1000,
    hi_t = hi_kg / 1000
  ) %>%
  filter(is.finite(year), is.finite(m_t), is.finite(lo_t), is.finite(hi_t))

# ---------------------------
# Top-N taxa + colors
# ---------------------------
top_species <- df %>%
  arrange(desc(m_t)) %>%
  distinct(species) %>%
  slice_head(n = top_n) %>%
  pull(species)

df <- df %>%
  mutate(
    highlight = factor(
      if_else(species %in% top_species, species, "Other"),
      levels = c(top_species, "Other")
    )
  )

top_colors <- hue_pal()(length(top_species))
names(top_colors) <- top_species
col_map <- c(Other = "grey75", top_colors)

breaks <- c(top_species[1], top_species[6],
                    top_species[2], top_species[7],
                    top_species[3], top_species[8],
                    top_species[4], top_species[9],
                    top_species[5], top_species[10])

labels <- parse(text = paste0("italic('", breaks, "')"))

# ---------------------------
# Plot
# ---------------------------
p_time <- ggplot(df, aes(x = year, y = m_t)) +
  geom_linerange(
    aes(ymin = lo_t, ymax = hi_t, colour = highlight),
    linewidth = 0.55, alpha = 0.90,
    show.legend = FALSE
  ) +
  geom_hline(yintercept=100, col='darkblue')+ # remove this
  geom_hline(yintercept = 35, col = "darkblue") +  # remove this
  geom_point(
    aes(colour = highlight),
    size = 2.6, alpha = 0.95
  ) +
  scale_colour_manual(
    values = col_map,
    breaks = breaks,
    labels = labels,
    name   = ""
  ) +
  scale_x_continuous(breaks = pretty_breaks(n = 10)) +
  scale_y_log10(limits = c(0.1, 600)) +  # scale_y_continuous(limits = c(0, 100)) if on original scale
  labs(x = "Publication year", y = "Mass estimate (tons)") +
  theme_nature +
  guides(colour = guide_legend(
    ncol = 2,
    byrow = TRUE,   
  ))

print(p_time)

ggsave(
  filename = file.path(fig_dir, "mass_time.pdf"),
  plot = p_time,
  device = cairo_pdf)

# ---------------------------
# Discoveries per year
# ---------------------------
disc_per_year <- df %>%
  filter(is.finite(year)) %>%
  count(year, name = "n_discoveries") %>%
  arrange(year)

mean(tail(disc_per_year, 25)$n_discoveries, na.rm = TRUE)
