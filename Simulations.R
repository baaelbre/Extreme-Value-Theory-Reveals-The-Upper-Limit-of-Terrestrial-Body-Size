# ============================================================
# EVT animal size demo figures (complete script)
# - 4 figures + 2x2 panel (A–D) with tags UNDER each plot
# - Plotly 3D axis titles (Trait 1 / Trait 2 / Density) made BIGGER
# - Exports: PDF (vector container), SVG (if svglite), PNG (300 dpi)
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(MASS)       # kde2d
  library(plotly)     # interactive 3D
  library(evd)        # GPD density
  library(purrr)
  library(tibble)
  library(truncnorm)
  
  # panel/export helpers
  library(patchwork)
  library(grid)
  library(png)
})

set.seed(1)

# ------------------------------------------------------------
# Global options and theoretical endpoint
# ------------------------------------------------------------
show_true_endpoint <- TRUE

b_star      <- 7
Trait1_star <- b_star
Trait2_star <- 1.2 * b_star + 0.3   # = 8.7

endpoint_df <- tibble(
  Trait1 = Trait1_star,
  Trait2 = Trait2_star
)

# ------------------------------------------------------------
# 0) Example data (replace with your real traits)
# ------------------------------------------------------------
n  <- 1000
T1 <- rtruncnorm(n, a = -Inf, b = b_star, mean = 5, sd = 0.5)
T2 <- 1.2 * T1 + rnorm(n, 0, 0.3) + 0.3

dat <- tibble(Trait1 = T1, Trait2 = T2)

# Thresholds (e.g. 0.95-quantiles)
u1 <- as.numeric(quantile(dat$Trait1, 0.95))
u2 <- as.numeric(quantile(dat$Trait2, 0.95))

# Domain ranges for rectangles / 3D surfaces
x_min      <- min(dat$Trait1)
x_max_data <- max(dat$Trait1)
x_max      <- max(x_max_data, Trait1_star)

y_min      <- min(dat$Trait2)
y_max_data <- max(dat$Trait2)
y_max      <- max(y_max_data, Trait2_star)

# Regime classification
dat <- dat |>
  mutate(
    E1 = Trait1 >= u1,
    E2 = Trait2 >= u2,
    regime = case_when(
      E1 & E2  ~ "(E1,E2)",
      E1 & !E2 ~ "(E1,¬E2)",
      !E1 & E2 ~ "(¬E1,E2)",
      TRUE     ~ "none"
    ),
    regime = factor(regime, levels = c("(¬E1,E2)", "(E1,¬E2)", "(E1,E2)", "none"))
  )

regime_cols <- c(
  "(¬E1,E2)" = "#1b9e77",
  "(E1,¬E2)" = "#377eb8",
  "(E1,E2)"  = "#e41a1c",
  "none"     = "grey75"
)

# ------------------------------------------------------------
# 1) 2D scatter with three tail regimes + theoretical endpoint
# ------------------------------------------------------------
p_scatter <- ggplot(dat, aes(x = Trait1, y = Trait2)) +
  annotate(
    "rect",
    xmin = x_min, xmax = u1,
    ymin = u2,    ymax = y_max,
    fill = regime_cols["(¬E1,E2)"], alpha = 0.04
  ) +
  annotate(
    "rect",
    xmin = u1,    xmax = x_max,
    ymin = y_min, ymax = u2,
    fill = regime_cols["(E1,¬E2)"], alpha = 0.04
  ) +
  annotate(
    "rect",
    xmin = u1,    xmax = x_max,
    ymin = u2,    ymax = y_max,
    fill = regime_cols["(E1,E2)"], alpha = 0.07
  ) +
  geom_point(aes(color = regime), size = 1.8, alpha = 0.8) +
  scale_color_manual(values = regime_cols, guide = "none") +
  geom_vline(xintercept = u1, linetype = "dashed", linewidth = 0.7, colour = "red3") +
  geom_hline(yintercept = u2, linetype = "dashed", linewidth = 0.7, colour = "red3") +
  { if (show_true_endpoint)
    geom_point(
      data = endpoint_df,
      aes(x = Trait1, y = Trait2),
      inherit.aes = FALSE,
      shape = 4, size = 4, stroke = 1.2,
      colour = "red3"
    ) else NULL } +
  labs(x = "Trait 1", y = "Trait 2") +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank()
  )

ggsave("scatter_thresholds_regimes.png", p_scatter, width = 6, height = 5, dpi = 300)

# ------------------------------------------------------------
# 2) 3D trait–trait–density surface + EVT tail extension
# ------------------------------------------------------------
kde <- MASS::kde2d(dat$Trait1, dat$Trait2, n = 60)

# Demo tail model: independent GPD × GPD, extended to endpoints
xi1_tail <- -0.3
xi2_tail <- -0.3

sigma1_tail <- (u1 - Trait1_star) * xi1_tail
sigma2_tail <- (u2 - Trait2_star) * xi2_tail

x_tail_grid <- seq(u1, Trait1_star, length.out = 60)
y_tail_grid <- seq(u2, Trait2_star, length.out = 60)

f1_tail <- evd::dgpd(x_tail_grid, loc = u1, scale = sigma1_tail, shape = xi1_tail)
f2_tail <- evd::dgpd(y_tail_grid, loc = u2, scale = sigma2_tail, shape = xi2_tail)

p_joint     <- mean(dat$Trait1 >= u1 & dat$Trait2 >= u2)
z_tail_grid <- p_joint * outer(f1_tail, f2_tail, "*")

# Base-plane z=0 rectangle template
z0 <- matrix(0, nrow = 2, ncol = 2)

p_3d <- plotly::plot_ly() |>
  add_surface(
    x = ~kde$x, y = ~kde$y, z = ~kde$z,
    opacity    = 0.5,
    colorscale = "Greys",
    showscale  = FALSE,
    showlegend = FALSE
  ) |>
  add_surface(
    x = ~x_tail_grid, y = ~y_tail_grid, z = ~z_tail_grid,
    opacity    = 0.8,
    colorscale = "Reds",
    showscale  = FALSE,
    showlegend = FALSE
  ) |>
  # base-plane colored regimes
  add_surface(
    x = c(x_min, u1), y = c(u2, y_max), z = z0,
    opacity    = 0.25,
    colorscale = list(
      c(0, unname(regime_cols["(¬E1,E2)"])),
      c(1, unname(regime_cols["(¬E1,E2)"]))
    ),
    showscale  = FALSE,
    showlegend = FALSE
  ) |>
  add_surface(
    x = c(u1, x_max), y = c(y_min, u2), z = z0,
    opacity    = 0.25,
    colorscale = list(
      c(0, unname(regime_cols["(E1,¬E2)"])),
      c(1, unname(regime_cols["(E1,¬E2)"]))
    ),
    showscale  = FALSE,
    showlegend = FALSE
  ) |>
  add_surface(
    x = c(u1, x_max), y = c(u2, y_max), z = z0,
    opacity    = 0.45,
    colorscale = list(
      c(0, unname(regime_cols["(E1,E2)"])),
      c(1, unname(regime_cols["(E1,E2)"]))
    ),
    showscale  = FALSE,
    showlegend = FALSE
  ) |>
  # base-plane points
  add_markers(
    data = filter(dat, regime == "none"),
    x = ~Trait1, y = ~Trait2, z = 0,
    marker = list(size = 3, color = regime_cols["none"]),
    showlegend = FALSE
  ) |>
  add_markers(
    data = filter(dat, regime == "(¬E1,E2)"),
    x = ~Trait1, y = ~Trait2, z = 0,
    marker = list(size = 4, color = regime_cols["(¬E1,E2)"]),
    showlegend = FALSE
  ) |>
  add_markers(
    data = filter(dat, regime == "(E1,¬E2)"),
    x = ~Trait1, y = ~Trait2, z = 0,
    marker = list(size = 4, color = regime_cols["(E1,¬E2)"]),
    showlegend = FALSE
  ) |>
  add_markers(
    data = filter(dat, regime == "(E1,E2)"),
    x = ~Trait1, y = ~Trait2, z = 0,
    marker = list(size = 4.5, color = regime_cols["(E1,E2)"]),
    showlegend = FALSE
  ) |>
  # threshold lines on z=0
  add_lines(
    x = c(u1, u1), y = c(y_min, y_max), z = 0,
    line = list(width = 5, dash = "dash", color = "red"),
    showlegend = FALSE
  ) |>
  add_lines(
    x = c(x_min, x_max), y = c(u2, u2), z = 0,
    line = list(width = 5, dash = "dash", color = "red"),
    showlegend = FALSE
  )

if (show_true_endpoint) {
  p_3d <- p_3d |>
    add_markers(
      x = Trait1_star, y = Trait2_star, z = 0,
      marker = list(symbol = "x", size = 5, color = "red"),
      showlegend = FALSE
    )
}

# BIGGER axis titles + ticks (Trait 1 / Trait 2 / Density)
axis_title_size <- 28  # <-- increase this further if needed
axis_tick_size  <- 16

p_3d <- p_3d |>
  layout(
    scene = list(
      xaxis  = list(
        title     = "Trait 1",
        titlefont = list(size = axis_title_size),
        tickfont  = list(size = axis_tick_size)
      ),
      yaxis  = list(
        title     = "Trait 2",
        titlefont = list(size = axis_title_size),
        tickfont  = list(size = axis_tick_size)
      ),
      zaxis  = list(
        title     = "Density",
        titlefont = list(size = axis_title_size),
        tickfont  = list(size = axis_tick_size)
      ),
      camera = list(eye = list(x = 1.7, y = 1.4, z = 1.2))
    ),
    showlegend = FALSE,
    title = ""
  )

# ------------------------------------------------------------
# 3) GPD shape parameter ξ: tail and endpoint illustration
# ------------------------------------------------------------
u     <- 0
sigma <- 1
xs    <- seq(0, 6, length.out = 400)
xis   <- c(-0.3, 0, 0.3)

dens_df <- purrr::map_dfr(xis, function(xi) {
  support <- if (xi < 0) xs[xs < -sigma/xi] else xs
  tibble(
    x    = support,
    dens = evd::dgpd(support + u, loc = u, scale = sigma, shape = xi),
    xi   = factor(xi)
  )
})

xi_cols <- c(`-0.3` = "#D55E00", `0` = "#000000", `0.3` = "#0072B2")

end_df <- tibble(
  xi    = factor(-0.3),
  x_end = -sigma / (-0.3)
)

p_xi <- ggplot(dens_df, aes(x = x, y = dens, colour = xi)) +
  geom_rect(
    data = end_df,
    aes(xmin = x_end, xmax = max(xs), ymin = -Inf, ymax = Inf),
    inherit.aes = FALSE,
    fill = "grey90", alpha = 0.6
  ) +
  geom_line(linewidth = 1.1) +
  geom_vline(
    data = end_df,
    aes(xintercept = x_end),
    linetype = "dashed",
    linewidth = 0.8,
    colour = xi_cols["-0.3"]
  ) +
  annotate(
    "text",
    x = end_df$x_end,
    y = 0.1,
    label = "y*",
    hjust = -0.1,
    size  = 4,
    colour = xi_cols["-0.3"],
    parse = FALSE
  ) +
  scale_colour_manual(
    values = xi_cols,
    name   = "Shape ξ",
    labels = c("ξ = -0.3", "ξ = 0", "ξ = 0.3")
  ) +
  labs(
    x = expression("Exceedance " * (y - u)),
    y = "GPD density"
  ) +
  theme_minimal(base_size = 13) +
  theme(panel.grid.minor = element_blank())

# ------------------------------------------------------------
# 4) Log–log tail survival for both traits (slope comparison)
# ------------------------------------------------------------
make_tail_curve <- function(ex, n_grid = 30, name = "Trait") {
  ex <- ex[ex > 0] |> sort()
  probs_grid <- seq(0.1, 0.99, length.out = n_grid)
  x_grid     <- as.numeric(quantile(ex, probs_grid))
  S_grid <- vapply(x_grid, function(x) mean(ex >= x), numeric(1L))
  tibble(exceed = x_grid, S = S_grid, Trait = name)
}

ex1 <- dat$Trait1[dat$Trait1 >= u1] - u1
ex2 <- dat$Trait2[dat$Trait2 >= u2] - u2

curves <- bind_rows(
  make_tail_curve(ex1, 30, "Trait 1"),
  make_tail_curve(ex2, 30, "Trait 2")
)

trait_cols <- c("Trait 1" = "#377eb8", "Trait 2" = "#1b9e77")

p_tail_loglog <- ggplot(curves) +
  geom_point(aes(x = exceed, y = S, colour = Trait), size = 2, alpha = 1) +
  scale_colour_manual(values = trait_cols, name = "Trait") +
  scale_x_log10() +
  scale_y_log10() +
  labs(
    x = expression("Exceedance " * (y - u)),
    y = "SF(x) = P(Y > u + x | Y > u)"
  ) +
  theme_minimal(base_size = 13) +
  theme(panel.grid.minor = element_blank())

ggsave("tail_loglog_traits.png", p_tail_loglog, width = 6, height = 5, dpi = 300)

# ------------------------------------------------------------
# 5) Unconditional survival vs distance to empirical maximum
# ------------------------------------------------------------
y1_tail <- dat$Trait1[dat$Trait1 >= u1]
y2_tail <- dat$Trait2[dat$Trait2 >= u2]

y1_max <- max(y1_tail)
y2_max <- max(y2_tail)

t1 <- (y1_max - y1_tail)
t2 <- (y2_max - y2_tail)

t1 <- t1[t1 > 0]
t2 <- t2[t2 > 0]

p_tail1 <- mean(dat$Trait1 > u1)
p_tail2 <- mean(dat$Trait2 > u2)

make_endpoint_curve <- function(t, y_tail, y_max, p_tail, n_grid = 30, name = "Trait") {
  t <- sort(t)
  probs_grid <- seq(0.01, 0.99, length.out = n_grid)
  t_grid     <- as.numeric(quantile(t, probs_grid))
  S_cond <- vapply(t_grid, function(tt) mean(y_tail > (y_max - tt)), numeric(1L))
  tibble(t = t_grid, S = p_tail * S_cond, Trait = name)
}

curves_ep <- bind_rows(
  make_endpoint_curve(t1, y1_tail, y1_max, p_tail1, 30, "Trait 1"),
  make_endpoint_curve(t2, y2_tail, y2_max, p_tail2, 30, "Trait 2")
)

p_endpoint_survival <- ggplot(curves_ep, aes(x = t, y = S, colour = Trait)) +
  geom_point(size = 2, alpha = 0.9) +
  scale_colour_manual(values = trait_cols, name = "Trait") +
  scale_x_log10() +
  scale_y_log10() +
  labs(
    x = expression(hat(t) == y[max] - y),
    y = expression(hat(S)(t) == P(Y > y[max] - t))
  ) +
  theme_minimal(base_size = 13) +
  theme(panel.grid.minor = element_blank())

ggsave("endpoint_distance_survival_uncond.png", p_endpoint_survival, width = 6, height = 5, dpi = 300)

# ------------------------------------------------------------
# 6) Make a static PNG of the plotly 3D plot (for panel)
# ------------------------------------------------------------
p3d_png  <- "p_3d_static.png"
p3d_html <- "p_3d_static.html"
ok <- FALSE

# Preferred: plotly::save_image (Kaleido)
ok <- tryCatch({
  plotly::save_image(p_3d, file = p3d_png, width = 2000, height = 1600, scale = 2)
  file.exists(p3d_png)
}, error = function(e) FALSE)

# Fallback: htmlwidget + webshot2
if (!ok) {
  if (!requireNamespace("htmlwidgets", quietly = TRUE)) install.packages("htmlwidgets")
  if (!requireNamespace("webshot2", quietly = TRUE))    install.packages("webshot2")
  
  ok <- tryCatch({
    htmlwidgets::saveWidget(p_3d, file = p3d_html, selfcontained = TRUE)
    webshot2::webshot(url = p3d_html, file = p3d_png, vwidth = 2000, vheight = 1600, zoom = 2)
    file.exists(p3d_png)
  }, error = function(e) FALSE)
}

if (!ok) stop("Could not export p_3d to PNG. Try installing/updating 'plotly' or 'webshot2'.")

img <- png::readPNG(p3d_png)
p_3d_static <- patchwork::wrap_elements(full = grid::rasterGrob(img, interpolate = TRUE))

# ------------------------------------------------------------
# 7) 2x2 panel with tags UNDER each plot + export
# ------------------------------------------------------------
out_base <- "EVT_demo_panel_2x2"
panel_w  <- 10   # inches
panel_h  <- 8    # inches

tag_strip <- function(tag, size = 14) {
  ggplot() +
    annotate(
      "text", x = 0.5, y = 0.5,
      label = paste0("(", tag, ")"),
      fontface = "bold",
      size = size / 3
    ) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
    theme_void() +
    theme(plot.margin = margin(0, 0, 0, 0))
}

tight_bottom <- theme(plot.margin = margin(5.5, 5.5, 0, 5.5))

A <- p_3d_static / tag_strip("A") + plot_layout(heights = c(1, 0.08))
B <- (p_scatter + tight_bottom) / tag_strip("B") + plot_layout(heights = c(1, 0.08))
C <- (p_xi + tight_bottom)      / tag_strip("C") + plot_layout(heights = c(1, 0.08))
D <- (p_endpoint_survival + tight_bottom) / tag_strip("D") + plot_layout(heights = c(1, 0.08))

panel_2x2 <- wrap_plots(A, B, C, D, ncol = 2, byrow = TRUE) +
  plot_layout(widths = c(1, 1), heights = c(1, 1))

# Vector exports (container is vector; panel includes a raster for A)
ggsave(paste0(out_base, ".pdf"), panel_2x2, width = panel_w, height = panel_h, device = cairo_pdf)
if (requireNamespace("svglite", quietly = TRUE)) {
  ggsave(paste0(out_base, ".svg"), panel_2x2, width = panel_w, height = panel_h, device = svglite::svglite)
}
ggsave(paste0(out_base, ".png"), panel_2x2, width = panel_w, height = panel_h, dpi = 300)
