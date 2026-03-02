# ============================================================
# Title : Extreme Value Theory Reveals the Upper Limit of
#         Terrestrial Body Size: Alligators
# Author: Bastiaan A. Van Velthoven
#
# Update:
# - Replace PP plot with a QQ plot at the end.
# - QQ plot shows ACTUAL quantiles on the original scale (cm),
#   i.e. theoretical tail quantiles Y = u + qgpd(...) vs empirical Y (for y>u).
# ============================================================

suppressPackageStartupMessages({
  library(readxl)
  library(POT)
  library(dplyr)
  library(ggplot2)
  library(tibble)
})

set.seed(42)

# ---------------------------
# USER OPTIONS
# ---------------------------
ci_level <- 0.95
q_u <- 0.95
QQ_LINE_QUANTILES <- c(0.1, 0.9)  # anchor line through these QQ-quantiles (of the QQ cloud)

# ---------------------------
# Paths
# ---------------------------
DATA_XLSX <- "Data/alligators_woodward.xlsx"
FIG_DIR   <- "Figures/Alligators"
if (!dir.exists(FIG_DIR)) dir.create(FIG_DIR, recursive = TRUE)

# ---------------------------
# Read data
# ---------------------------
df_raw <- read_excel(DATA_XLSX)

TL <- suppressWarnings(as.numeric(df_raw$TL))
TL <- TL[is.finite(TL) & TL > 0]

y <- TL
xlab_size <- "Total length (cm)"

n <- length(y)
cat("\nN =", n, "\n")
cat("Summary:\n"); print(summary(y))

# ---------------------------
# Threshold diagnostics (POT)
# ---------------------------
u_range  <- as.numeric(quantile(y, probs = c(0.80, 0.98), na.rm = TRUE))

mrlplot(y, u.range = u_range, nt = 30, conf = ci_level,
        xlab = "u", ylab = "mean excess")
tcplot(y, u.range = u_range, nt = 30, conf = ci_level, which = 1, ask = FALSE)
tcplot(y, u.range = u_range, nt = 30, conf = ci_level, which = 2, ask = FALSE)

# ---------------------------
# Choose threshold u and fit tail model (POT)
# ---------------------------
u <- as.numeric(quantile(y, probs = q_u, na.rm = TRUE))
cat("\nChosen threshold: q_u =", q_u, " -> u =", signif(u, 6), "\n")

fit <- fitgpd(y, threshold = u, est = "mle")
print(fit)
summary(fit)

sigma_hat <- unname(fit$param["scale"])
xi_hat    <- unname(fit$param["shape"])

p_u <- mean(y > u)
n_u <- sum(y > u)

cat("\nTail fit (MLE):\n")
cat("  sigma_hat =", signif(sigma_hat, 6), "\n")
cat("  xi_hat    =", signif(xi_hat, 6), "\n")
cat("  exceedances n_u =", n_u, " (p_u =", signif(p_u, 6), ")\n")

# ---------------------------
# Endpoint estimate + delta-method
# ---------------------------
ystar_hat <- u - sigma_hat / xi_hat

V <- fit$var.cov
if (!is.null(rownames(V)) && all(c("scale", "shape") %in% rownames(V))) {
  V <- V[c("scale", "shape"), c("scale", "shape")]
}

grad <- c(-1 / xi_hat, sigma_hat / (xi_hat^2))
var_ystar <- as.numeric(t(grad) %*% V %*% grad)
se_ystar  <- sqrt(var_ystar)
ystar_ucb <- ystar_hat + qnorm(ci_level) * se_ystar

cat("\nFinite endpoint implied:\n")
cat("  y* =", signif(ystar_hat, 7), "\n")
cat("  se(y*) =", signif(se_ystar, 6), "\n")
cat("  one-sided", 100 * ci_level, "% UCB =", signif(ystar_ucb, 7), "\n")

# ---------------------------
# One-sided Wald test for xi < 0 
# ---------------------------
se_xi <- sqrt(V[2,2])
z_xi  <- xi_hat / se_xi                
p_xi  <- pnorm(z_xi)  
p_xi

# ---------------------------
# Exceedance probability for Stokes
# ---------------------------
y_stokes <- 450
N_year   <- 7600

p_stokes <- pgpd(y_stokes,
                 loc = u,
                 scale = sigma_hat,
                 shape = xi_hat,
                 lower.tail = FALSE,
                 lambda = 1 - p_u)

cat("\n================= Exceedance probability =================\n")
cat("Target y =", y_stokes, "cm\n")
cat("lambda  = P(Y<=u) =", signif(1 - p_u, 6), " (so p_u =", signif(p_u, 6), ")\n")
cat("p_hat   = P(Y>y)  =", signif(p_stokes, 8), "\n")
cat("Expected waiting time (years) ~", signif((1 / p_stokes) / N_year, 6), "\n")
cat("=========================================================\n")

# ============================================================
# ggplot styling helpers
# ============================================================
theme_nature <- theme_minimal(base_family = "Arial", base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.25),
    axis.title = element_text(size = 13, face = "bold"),
    axis.text  = element_text(size = 11),
    plot.title = element_blank(),
    plot.margin = margin(8, 8, 8, 8),
    legend.position = "none"
  )

col_tail   <- "#2C7FB8"
col_max    <- "#444444"
col_stokes <- "#D95F0E"
col_endpt  <- "#7B3294"

# ============================================================
# FIG 2 — Distribution of alligator TL
# ============================================================
df_y <- tibble(y = y)

p2 <- ggplot(df_y, aes(x = y)) +
  annotate("rect", xmin = u, xmax = Inf, ymin = -Inf, ymax = Inf,
           fill = col_tail, alpha = 0.08) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 40, fill = "grey82", color = "white") +
  geom_vline(xintercept = u, linewidth = 1.1, color = col_tail) +
  coord_cartesian(xlim = c(100, 600)) +
  labs(x = xlab_size, y = "Density") +
  theme_nature
p2

ggsave(file.path(FIG_DIR, "Fig2.png"),
       p2, dpi = 600, width = 6.6, height = 4.6, units = "in")
ggsave(
  file.path(FIG_DIR, "Fig2.pdf"),,
  plot = p2,
  device = cairo_pdf,
  width = 7.2, height = 4.8, units = "in"
)


# ============================================================
# FIG 3a — Empirical and fitted unconditional survival
# ============================================================
y_sorted <- sort(y)
k_all <- seq_along(y_sorted)
S_uncond_emp <- (n - k_all + 0.5) / (n + 1)

df_emp_all <- tibble(y = y_sorted, S = S_uncond_emp) %>%
  filter(y >= u)

y_grid_max <- max(y, na.rm = TRUE) + 30
if (is.finite(ystar_hat)) y_grid_max <- min(ystar_hat - 1e-6, y_grid_max)
y_grid_max <- max(y_grid_max, u + 1e-3)

y_grid <- seq(u, y_grid_max, length.out = 600)

S_uncond_fit <- pgpd(y_grid,
                     loc = u,
                     scale = sigma_hat,
                     shape = xi_hat,
                     lower.tail = FALSE,
                     lambda = 1 - p_u)

df_fit <- tibble(y = y_grid, S = S_uncond_fit)

p3a <- ggplot() +
  geom_point(data = df_emp_all, aes(x = y, y = S),
             size = 2.1, alpha = 0.55, color = "grey25") +
  geom_line(data = df_fit, aes(x = y, y = S),
            linewidth = 1.2, color = col_tail) +
  geom_vline(xintercept = u,        linewidth = 1.0, color = col_tail) +
  geom_vline(xintercept = y_stokes, linewidth = 1.0, color = col_stokes) +
  { if (is.finite(ystar_hat)) geom_vline(xintercept = ystar_hat, linewidth = 1.0, color = col_endpt, linetype = "dotdash") } +
  scale_y_log10() +
  labs(x = xlab_size, y = "Survival probability P(Y>y)") +
  theme_nature
p3a

ggsave(file.path(FIG_DIR, "Fig3a.png"),
       p3a, dpi = 600, width = 6.6, height = 4.6, units = "in")

# ============================================================
# FIG 3b — QQ plot (ACTUAL quantiles on original scale, cm)
#   x-axis: theoretical GPD quantiles for Y (u + exceedance quantile)
#   y-axis: empirical tail observations Y (for y>u)
# ============================================================
y_u <- sort(y[y > u])
n_exc <- length(y_u)

p_plot <- (seq_len(n_exc) - 0.5) / n_exc

q_exc <- qgpd(p_plot,
              loc = 0,
              scale = sigma_hat,
              shape = xi_hat,
              lambda = 0)   # conditional exceedance quantiles

q_Y <- u + q_exc            # convert exceedance quantiles to Y-scale

df_qq <- tibble(
  q_theo = q_Y,
  y_emp  = y_u
) %>% filter(is.finite(q_theo), is.finite(y_emp))

# Reference line anchored at chosen quantiles of the QQ cloud
a  <- as.numeric(quantile(df_qq$q_theo, QQ_LINE_QUANTILES[1], na.rm = TRUE))
b  <- as.numeric(quantile(df_qq$q_theo, QQ_LINE_QUANTILES[2], na.rm = TRUE))
ya <- as.numeric(quantile(df_qq$y_emp,  QQ_LINE_QUANTILES[1], na.rm = TRUE))
yb <- as.numeric(quantile(df_qq$y_emp,  QQ_LINE_QUANTILES[2], na.rm = TRUE))

slope_line <- (yb - ya) / (b - a)
inter_line <- ya - slope_line * a

p3b <- ggplot(df_qq, aes(x = q_theo, y = y_emp)) +
  geom_abline(intercept = inter_line, slope = slope_line,
              linewidth = 1.1, color = col_tail) +
  geom_point(size = 2.1, alpha = 0.60, color = "grey25") +
  labs(
    x = "Theoretical quantiles",
    y = "Empirical quantiles"
  ) +
  theme_nature
p3b

ggsave(file.path(FIG_DIR, "Fig3b.png"),
       p3b, dpi = 600, width = 6.6, height = 4.6, units = "in")

cat("\nSaved figures:\n",
    " - ", file.path(FIG_DIR, "Fig2.png"), "\n",
    " - ", file.path(FIG_DIR, "Fig3a.png"), "\n",
    " - ", file.path(FIG_DIR, "Fig3b.png"), "\n", sep = "")
