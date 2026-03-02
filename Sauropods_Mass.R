# ============================================================
# Title : Extreme Value Theory Reveals the Upper Limit of
#         Terrestrial Body Size: Sauropods (univariate POT on mass)
# Author: Bastiaan A. Van Velthoven
# ============================================================
library(readxl)
library(POT)
library(numDeriv)
library(dplyr)

set.seed(42)


# ---------------------------
# Read data
# ---------------------------
DATA_XLSX <- "Data/DEmic23_updated_Supplemental_Data_withPubYear update Nov25_2025.xlsx"
df_raw    <- read_excel(DATA_XLSX)

col_mass <- "Campione est mass quadratic (kg)"
col_spec <- "genus and species"   # optional

df_mass <- df_raw |>
  transmute(
    specimen = if (col_spec %in% names(df_raw)) as.character(.data[[col_spec]]) else NA_character_,
    mass_kg  = suppressWarnings(as.numeric(.data[[col_mass]]))
  ) |>
  filter(is.finite(mass_kg), mass_kg > 0)

y    <- df_mass$mass_kg
ylab <- "Mass (kg)"

cat("\nN =", length(y), "\n")
cat("Summary (analysis scale):\n"); print(summary(y))

# ---------------------------
# Threshold selection
# ---------------------------
ci_level <- 0.95

u_range <- as.numeric(quantile(y, probs = c(0.60, 0.95), na.rm = TRUE))

mrlplot(y, u.range = u_range, nt = 30, conf = ci_level,
             xlab = "u", ylab = "mean excess")
tcplot(y, u.range = u_range, nt = 30, conf = ci_level, which = 1, ask = FALSE)
tcplot(y, u.range = u_range, nt = 30, conf = ci_level, which = 2, ask = FALSE)

u <- 34000 

cat("\nThreshold u =", u, " (", ylab, ")\n", sep = "")
cat("Exceedances:", sum(y > u), "\n")

fit <- POT::fitgpd(y, threshold = u, est = "mle")
print(summary(fit))

sigma_hat <- unname(fit$param["scale"])
xi_hat    <- unname(fit$param["shape"])

cat("\nGPD MLE:\n")
cat("  sigma_hat =", signif(sigma_hat, 6), "\n")
cat("  xi_hat    =", signif(xi_hat, 6), "\n")

# ---------------------------
# Endpoint estimate (and delta appr.)
# y* = u - sigma/xi
# ---------------------------

ystar_hat <- u - sigma_hat / xi_hat
cat("\nFinite endpoint implied (on", ylab, "scale): y* =", signif(ystar_hat, 7), "\n")

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
      cat("  se(y*) =", signif(se_ystar, 6), "\n")
      cat("  UCB_", 100*ci_level, "% = ", signif(ystar_ucb_delta, 7),
          " (", ylab, ")\n", sep = "")
    }
  }
}

# ---------------------------
# goodness-of-fit
# ---------------------------
POT::pp(fit)

