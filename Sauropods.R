# ============================================================
# Title : Extreme Value Theory Reveals the Upper Limit of
#         Terrestrial Body Size: Sauropods (univariate POT on mass)
# Author: Bastiaan A. Van Velthoven
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
REMOVE_LARGEST <- FALSE     # TRUE = exclude largest from training
u <- 34000                  # default sauropod threshold in kg
ci_level <- 0.95                   # one-sided UCB level for endpoint
QQ_LINE_QUANTILES <- c(0.05, 0.95) # anchor line through these QQ-quantiles (of the QQ cloud)

# ---------------------------
# Paths
# ---------------------------
DATA_XLSX <- "Data/sauropods_demic.xlsx"
FIG_DIR   <- "Figures/Sauropods"
if (!dir.exists(FIG_DIR)) dir.create(FIG_DIR, recursive = TRUE)

# ---------------------------
# Read data
# ---------------------------
col_mass <- "Campione est mass quadratic (kg)"
col_spec <- "genus and species"   # optional

df_raw <- read_excel(DATA_XLSX)

df_mass <- df_raw |>
  transmute(
    specimen = if (col_spec %in% names(df_raw)) as.character(.data[[col_spec]]) else NA_character_,
    mass_kg  = suppressWarnings(as.numeric(.data[[col_mass]]))
  ) |>
  filter(is.finite(mass_kg), mass_kg > 0)

# Analysis scale (kg)
y_all <- df_mass$mass_kg

# Remove the single largest training point by value (optional)
if (REMOVE_LARGEST) {
  idx_max <- which.max(y_all)
  cat("\nRemoving largest specimen from training: mass_kg =",
      signif(y_all[idx_max], 7), "kg\n")
  df_train <- df_mass[-idx_max, , drop = FALSE]
} else {
  df_train <- df_mass
}

y <- df_train$mass_kg
n <- length(y)

cat("\nN (training) =", n, "\n")
cat("Summary (kg):\n"); print(summary(y))

# ---------------------------
# Threshold diagnostics (POT) — still on kg
# ---------------------------
u_range <- as.numeric(quantile(y, probs = c(0.80, 0.95), na.rm = TRUE))
mrlplot(y, u.range = u_range, nt = 30, conf = ci_level,
        xlab = "u (kg)", ylab = "mean excess")
tcplot(y, u.range = u_range, nt = 30, conf = ci_level, which = 1, ask = FALSE)
tcplot(y, u.range = u_range, nt = 30, conf = ci_level, which = 2, ask = FALSE)

# ---------------------------
# Choose threshold u and fit tail model (POT) — in kg
# ---------------------------
cat("\nChosen threshold (fixed): u =", signif(u, 7), "kg\n")


fit <- fitgpd(y, threshold = u, est = "mle")
print(fit)
summary(fit)

sigma_hat <- unname(fit$param["scale"])
xi_hat    <- unname(fit$param["shape"])



p_u <- mean(y > u)
n_u <- sum(y > u)

cat("\nTail fit (MLE):\n")
cat("  sigma_hat =", signif(sigma_hat, 7), "\n")
cat("  xi_hat    =", signif(xi_hat, 7), "\n")
cat("  exceedances n_u =", n_u, " (p_u =", signif(p_u, 6), ")\n")

# ---------------------------
# Endpoint estimate + delta-method SE
# ---------------------------

ystar_hat <- u - sigma_hat / xi_hat
V <- fit$var.cov

grad <- c(-1 / xi_hat, sigma_hat / (xi_hat^2))
var_ystar <- as.numeric(t(grad) %*% V %*% grad)

se_ystar  <- sqrt(var_ystar)
ystar_ucb <- ystar_hat + qnorm(ci_level) * se_ystar

cat("\nFinite endpoint implied:\n")
cat("  y* =", signif(ystar_hat, 8), "kg\n")
cat("  se(y*) =", signif(se_ystar, 7), "\n")
cat("  one-sided", 100 * ci_level, "% UCB =", signif(ystar_ucb, 8), "kg\n")

# ---------------------------
# One-sided Wald test for xi < 0 
# ---------------------------
se_xi <- sqrt(V[2,2])
z_xi  <- xi_hat / se_xi                
p_xi  <- pnorm(z_xi)  
p_xi

# ---------------------------
# Exceedance probability for Argentinosaurus
# ---------------------------
y_arg <- 62000

p_arg <- pgpd(y_arg,
              loc = u,
              scale = sigma_hat,
              shape = xi_hat,
              lower.tail = FALSE,
              lambda = 1 - p_u)

cat("\n================= Exceedance probability =================\n")
cat("Target y =", y_arg, "kg (Argentinosaurus)\n")
cat("lambda = P(Y<=u) =", signif(1 - p_u, 6), " (so p_u =", signif(p_u, 6), ")\n")
cat("p_hat  = P(Y>y)  =", signif(p_arg, 8), "\n")
cat("=========================================================\n")

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
col_arg    <- "#D95F0E"
col_endpt  <- "#7B3294"

kg_to_t <- function(x) x / 1000

u_t       <- u/1000
y_t       <- u/1000
y_arg_t   <- u/1000
ystar_t   <- if (is.finite(ystar_hat)) kg_to_t(ystar_hat) else NA_real_

# ============================================================
# FIG 4 — Distribution 
# ============================================================
df_y <- tibble(y_t = y_t)

p4 <- ggplot(df_y, aes(x = y_t)) +
  annotate("rect", xmin = u_t, xmax = Inf, ymin = -Inf, ymax = Inf,
           fill = col_tail, alpha = 0.08) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 35, fill = "grey82", color = "white") +
  geom_vline(xintercept = u_t, linewidth = 1.1, color = col_tail) +
  coord_cartesian(xlim = c(0, max(200, max(y_t, na.rm = TRUE) * 1.15))) +
  labs(x = "Body mass (tons)", y = "Density") +
  theme_nature
p4

ggsave(file.path(FIG_DIR, "Fig4_sauropod_mass_distribution_tons.png"),
       p4, dpi = 600, width = 6.6, height = 4.6, units = "in")

# ============================================================
# FIG 3c — Tail survival P(Y>y) 
# ============================================================
y_sorted <- sort(y)                      # kg
k_all <- seq_along(y_sorted)
S_uncond_emp <- (n - k_all + 0.5) / (n + 1)

df_emp_all <- tibble(
  y_t = kg_to_t(y_sorted),
  S   = S_uncond_emp
) %>% filter(y_t >= u_t)

# grid in kg for probabilities, then convert for x-axis
y_grid_max <- max(y, na.rm = TRUE) + 20000
if (is.finite(ystar_hat)) y_grid_max <- min(ystar_hat - 1e-6, y_grid_max)
y_grid_max <- max(y_grid_max, u + 1e-3)

y_grid <- seq(u, y_grid_max, length.out = 600)      # kg
S_uncond_fit <- pgpd(y_grid,
                     loc = u,
                     scale = sigma_hat,
                     shape = xi_hat,
                     lower.tail = FALSE,
                     lambda = 1 - p_u)

df_fit <- tibble(y_t = kg_to_t(y_grid), S = S_uncond_fit)

p3c <- ggplot() +
  geom_point(data = df_emp_all, aes(x = y_t, y = S),
             size = 2.1, alpha = 0.55, color = "grey25") +
  geom_line(data = df_fit, aes(x = y_t, y = S),
            linewidth = 1.2, color = col_tail) +
  geom_vline(xintercept = u_t,     linewidth = 1.0, color = col_tail) +
  geom_vline(xintercept = ystar_t, linewidth = 1.0, color = col_endpt, linetype = "dotdash") +
  scale_y_log10() +
  labs(x = "Body mass (tons)", y =  "Survival probability P(Y>y)") +
  theme_nature
p3c

ggsave(file.path(FIG_DIR, "Fig3c.png"),
       p3c, dpi = 600, width = 6.6, height = 4.6, units = "in")

# ============================================================
# FIG 3d — QQ plot for points y>u
#   x-axis: theoretical quantiles of Y (NOT exceedances), i.e. u + qgpd(...)
#   y-axis: empirical tail observations Y (NOT exceedances)
#   Reference line anchored at quantiles of the QQ cloud
# ============================================================
y_u <- sort(y[y > u])        # kg, actual tail observations
n_exc <- length(y_u)

p_plot <- (seq_len(n_exc) - 0.5) / n_exc

# theoretical quantiles for exceedances, then add threshold to get Y-quantiles
q_exc <- qgpd(p_plot, loc = 0, scale = sigma_hat, shape = xi_hat, lambda = 0)  # exceedances (kg)
q_Y   <- u + q_exc                                                            # Y-quantiles (kg)

df_qq <- tibble(
  q_theo_t = kg_to_t(q_Y),
  y_emp_t  = kg_to_t(y_u)
) %>% filter(is.finite(q_theo_t), is.finite(y_emp_t))

# line anchored at chosen quantiles of the QQ cloud
a  <- as.numeric(quantile(df_qq$q_theo_t, QQ_LINE_QUANTILES[1], na.rm = TRUE))
b  <- as.numeric(quantile(df_qq$q_theo_t, QQ_LINE_QUANTILES[2], na.rm = TRUE))
ya <- as.numeric(quantile(df_qq$y_emp_t,  QQ_LINE_QUANTILES[1], na.rm = TRUE))
yb <- as.numeric(quantile(df_qq$y_emp_t,  QQ_LINE_QUANTILES[2], na.rm = TRUE))

slope_line <- (yb - ya) / (b - a)
inter_line <- ya - slope_line * a

p3d <- ggplot(df_qq, aes(x = q_theo_t, y = y_emp_t)) +
  geom_abline(intercept = inter_line, slope = slope_line,
              linewidth = 1.1, color = col_tail) +
  geom_point(size = 2.1, alpha = 0.60, color = "grey25") +
  labs(
    x = "Theoretical quantiles",
    y = "Empirical quantiles"
  ) +
  theme_nature
p3d

ggsave(file.path(FIG_DIR, "Fig3d.png"),
       p3d, dpi = 600, width = 6.6, height = 4.6, units = "in")

cat("\nSaved sauropod figures:\n",
    " - ", file.path(FIG_DIR, "Fig4_sauropod_mass_distribution_tons.png"), "\n",
    " - ", file.path(FIG_DIR, "Fig5a_sauropod_tail_survival_tons.png"), "\n",
    " - ", file.path(FIG_DIR, "Fig5b_sauropod_QQplot_tons.png"), "\n", sep = "")

