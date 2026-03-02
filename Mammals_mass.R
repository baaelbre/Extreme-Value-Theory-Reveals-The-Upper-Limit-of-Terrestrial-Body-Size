# ============================================================
# Title : Extreme Value Theory Reveals the Upper Limit of
#         Terrestrial Body Size: Mammals (PHYLACINE 1.2)
# Author: Bastiaan A. Van Velthoven
# File  : phylacine_mass_pot_endpoint_POT_delta_only.R
#
# Mirrors the "alligator POT version" but DELTA ONLY:
# - Threshold selection: mrlplot + tcplot (shape + modified scale) + lmomplot
# - Fit: POT::fitgpd (MLE)
# - Endpoint: y* = u - sigma/xi (finite iff xi<0)
# - Uncertainty: Delta-method one-sided UCB for y* (no profile likelihood)
# - Diagnostics: PP, QQ, dens (POT)
#
# Data: PHYLACINE Trait_data.csv (species-level adult mass)
# Columns used:
#   Binomial.1.2, Order.1.2, Terrestrial, Marine, Freshwater, Aerial,
#   Mass.g, IUCN.Status.1.2, Diet.Plant
# ============================================================

suppressPackageStartupMessages({
  library(POT)
  library(numDeriv)  # Hessian for delta method
  library(dplyr)
})

set.seed(42)

# (Optional) If plots don't show because a file device is open, close them:
# while (!is.null(dev.list())) dev.off()

# ---------------------------
# Read data
# ---------------------------
DATA_CSV <- "Data/Trait_data.csv"   # <-- set to your path
df_raw   <- read.csv(DATA_CSV, stringsAsFactors = FALSE, check.names = FALSE)

# ---------------------------
# Filters (edit as desired)
# ---------------------------
terrestrial_only   <- TRUE   # keep strictly terrestrial taxa
herbivore_only     <- TRUE   # diet filter using Diet.Plant (%)
diet_plant_min     <- 80     # e.g. 80 or 90 for "mostly herbivorous"
mega_only          <- FALSE  # optional: focus on big herbivores only
mega_mass_min_kg   <- 500    # only used if mega_only=TRUE

include_extinct    <- TRUE   # keep EX/EW/EP in PHYLACINE
# (PHYLACINE uses IUCN.Status.1.2 with codes incl EX, EW, EP)

# ---------------------------
# Choose analysis scale
# ---------------------------
use_log_mass <- FALSE  # TRUE -> analyze log(mass_kg)

# ---------------------------
# Columns (PHYLACINE Trait_data.csv)
# ---------------------------
col_species <- "Binomial.1.2"
col_order   <- "Order.1.2"
col_mass_g  <- "Mass.g"
col_terr    <- "Terrestrial"
col_marine  <- "Marine"
col_fresh   <- "Freshwater"
col_aerial  <- "Aerial"
col_iucn    <- "IUCN.Status.1.2"
col_dietP   <- "Diet.Plant"

stopifnot(all(c(col_species, col_order, col_mass_g, col_terr, col_marine, col_fresh, col_aerial, col_iucn, col_dietP) %in% names(df_raw)))

# ---------------------------
# Extract + clean
# ---------------------------
df_mass <- df_raw |>
  transmute(
    species   = as.character(.data[[col_species]]),
    order     = as.character(.data[[col_order]]),
    iucn      = as.character(.data[[col_iucn]]),
    terr      = suppressWarnings(as.numeric(.data[[col_terr]])),
    marine    = suppressWarnings(as.numeric(.data[[col_marine]])),
    fresh     = suppressWarnings(as.numeric(.data[[col_fresh]])),
    aerial    = suppressWarnings(as.numeric(.data[[col_aerial]])),
    dietPlant = suppressWarnings(as.numeric(.data[[col_dietP]])),
    mass_g    = suppressWarnings(as.numeric(.data[[col_mass_g]]))
  ) |>
  filter(is.finite(mass_g), mass_g > 0) |>
  mutate(
    mass_kg  = mass_g / 1000,
    log_mass = log(mass_kg)
  )

# ---------------------------
# Apply filters
# ---------------------------
if (isTRUE(terrestrial_only)) {
  # strict terrestrial: terrestrial==1 and not marine/freshwater/aerial
  df_mass <- df_mass |>
    filter(terr == 1, marine == 0, fresh == 0, aerial == 0)
}

if (!isTRUE(include_extinct)) {
  df_mass <- df_mass |>
    filter(!iucn %in% c("EX", "EW", "EP"))
}

if (isTRUE(herbivore_only)) {
  df_mass <- df_mass |>
    filter(is.finite(dietPlant), dietPlant >= diet_plant_min)
}

if (isTRUE(mega_only)) {
  df_mass <- df_mass |>
    filter(mass_kg >= mega_mass_min_kg)
}

# ---------------------------
# Analysis vector
# ---------------------------
y    <- if (use_log_mass) df_mass$log_mass else df_mass$mass_kg
ylab <- if (use_log_mass) "log(mass_kg)" else "mass_kg"

cat("\n================ PHYLACINE POT on", ylab, "================\n")
cat("N =", length(y), "\n")
cat("Summary (analysis scale):\n"); print(summary(y))

cat("\nTop 15 by mass (kg):\n")
print(df_mass |>
        arrange(desc(mass_kg)) |>
        select(species, order, iucn, mass_kg, dietPlant) |>
        head(15))

if (use_log_mass) {
  cat("\nSummary (kg scale):\n"); print(summary(df_mass$mass_kg))
}
cat("=========================================================\n\n")

stopifnot(length(y) >= 30)

# ---------------------------
# Threshold selection (POT)
# ---------------------------
ci_level <- 0.90  # one-sided UCB level (also used for plots)

# For mammals, you often want higher thresholds than 0.70–0.95; edit if needed
u_range <- as.numeric(quantile(y, probs = c(0.80, 0.99), na.rm = TRUE))

par(mfrow = c(2, 2))

POT::mrlplot(y, u.range = u_range, nt = 30, conf = ci_level,
             xlab = "u", ylab = "mean excess")
title(main = paste("MRL:", ylab))

POT::tcplot(y, u.range = u_range, nt = 30, conf = ci_level, which = 1, ask = FALSE)
title(main = paste("tcplot (shape):", ylab))

POT::tcplot(y, u.range = u_range, nt = 30, conf = ci_level, which = 2, ask = FALSE)
title(main = paste("tcplot (modified scale):", ylab))

POT::lmomplot(y, u.range = u_range, nt = 30, identify = FALSE)
title(main = paste("lmomplot:", ylab))

par(mfrow = c(1, 1))

# ---------------------------
# Choose threshold u and fit GPD (POT)
# ---------------------------
# IMPORTANT: if use_log_mass=TRUE, u must be on the log scale.
u <- NA_real_  # <-- set this after inspecting plots, OR leave NA to use a default

if (!is.finite(u)) {
  u <- as.numeric(quantile(y, probs = 0.97, na.rm = TRUE))
  cat(glue("No u supplied -> using default u = quantile(y, 0.97) = {signif(u,7)} ({ylab})\n"))
}

cat("\nThreshold u =", signif(u, 7), " (", ylab, ")\n", sep = "")
cat("Exceedances:", sum(y > u), "\n")

fit <- POT::fitgpd(y, threshold = u, est = "mle")
print(summary(fit))

sigma_hat <- unname(fit$param["scale"])
xi_hat    <- unname(fit$param["shape"])

cat("\nGPD MLE:\n")
cat("  sigma_hat =", signif(sigma_hat, 7), "\n")
cat("  xi_hat    =", signif(xi_hat, 7), "\n")

# ---------------------------
# Endpoint estimate + DELTA UCB (finite iff xi < 0)
# y* = u - sigma/xi
# ---------------------------
if (!is.finite(xi_hat) || xi_hat >= 0) {
  cat("\nNo finite endpoint implied (xi_hat >= 0). No endpoint CI.\n")
} else {
  
  ystar_hat <- u - sigma_hat / xi_hat
  cat("\nFinite endpoint implied (on", ylab, "scale): y* =", signif(ystar_hat, 8), "\n")
  if (use_log_mass) cat("Back-transformed endpoint (kg):", signif(exp(ystar_hat), 8), "\n")
  
  # exceedances
  x <- y[y > u] - u
  n <- length(x)
  
  # negative log-likelihood for exceedances (GPD with loc=0)
  nll_gpd <- function(par) {
    sig <- par[1]; xi <- par[2]
    if (!is.finite(sig) || !is.finite(xi) || sig <= 0) return(Inf)
    
    # exponential limit
    if (abs(xi) < 1e-8) {
      return(n * log(sig) + sum(x) / sig)
    }
    
    t <- 1 + xi * x / sig
    if (any(t <= 0) || any(!is.finite(t))) return(Inf)
    
    n * log(sig) + (1/xi + 1) * sum(log(t))
  }
  
  # observed information (numerical Hessian of nll)
  H <- try(numDeriv::hessian(nll_gpd, x = c(sigma_hat, xi_hat)), silent = TRUE)
  
  if (inherits(H, "try-error")) {
    warning("Hessian failed; delta-method UCB not available.")
  } else {
    
    V <- try(solve(H), silent = TRUE)
    
    if (inherits(V, "try-error")) {
      warning("Hessian not invertible; delta-method UCB not available.")
    } else {
      
      # gradient of y* = u - sigma/xi
      g <- c(d_sigma = -1 / xi_hat,
             d_xi    =  sigma_hat / (xi_hat^2))
      
      var_ystar <- as.numeric(t(g) %*% V %*% g)
      
      if (!is.finite(var_ystar) || var_ystar <= 0) {
        warning("Delta-method variance not positive; delta-method UCB not available.")
      } else {
        
        se_ystar <- sqrt(var_ystar)
        z <- qnorm(ci_level)  # one-sided upper
        ystar_ucb_delta <- ystar_hat + z * se_ystar
        
        cat("\nDelta method (one-sided upper):\n")
        cat("  se(y*) =", signif(se_ystar, 7), "\n")
        cat("  UCB_", 100*ci_level, "% = ", signif(ystar_ucb_delta, 8),
            " (", ylab, ")\n", sep = "")
        if (use_log_mass) cat("  UCB mass (kg) =", signif(exp(ystar_ucb_delta), 8), "\n")
        
        # --- Delta-method normal approx plot for y* ---
        xx <- seq(ystar_hat - 4*se_ystar, ystar_hat + 4*se_ystar, length.out = 400)
        yy <- dnorm(xx, mean = ystar_hat, sd = se_ystar)
        
        plot(xx, yy, type = "l",
             main = paste0("Delta-method approx: y* (", ylab, ")"),
             xlab = paste0("y* (", ylab, ")"), ylab = "density")
        abline(v = ystar_hat, lty = 2, lwd = 2)
        abline(v = ystar_ucb_delta, lty = 3, lwd = 2)
        legend("topright",
               legend = c("Normal approx", "MLE y*", paste0("UCB ", 100*ci_level, "%")),
               lty = c(1, 2, 3), bty = "n")
      }
    }
  }
}

# ---------------------------
# Quick goodness-of-fit diagnostics (POT)
# ---------------------------
par(mfrow = c(1, 3))
POT::pp(fit, ci = TRUE)
POT::qq(fit, ci = TRUE)
POT::dens(fit)
par(mfrow = c(1, 1))