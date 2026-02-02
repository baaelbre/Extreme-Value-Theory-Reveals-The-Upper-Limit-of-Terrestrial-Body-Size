// =============================================================
// hier_gp_weibull_sauropod_flexible.stan
// Collapsed EVT model with optional truncation for y* prior
// =============================================================

functions {
  // ---- Derived mass for reporting (optional) ----
  real mass_of_ystar(real ystar, real alpha, real beta, real gamma, real c0) {
    real z = ystar - c0;
    real m_log = alpha + beta * z + gamma * z * z; // log grams
    return exp(m_log) / 1e6;                       // tons
  }
}

data {
  // ---- Exceedances ----
  int<lower=1> N;
  vector[N] y;
  real u;
  real ymax;                 // sample max among exceedances

  // ---- Prior support for endpoint ----
  real y_min;                // lower bound (user or data-based)
  real y_max;                // legacy/user cap (finite or Inf)
  real y_cap;                // hard cap for numerical safety

  // ---- Direct priors ----
  real mu_y0;
  real<lower=0> sd_mu_y;

  real mu_xi0;
  real<lower=0> sd_xi0;

  // ---- Allometry for reporting ----
  real alpha;
  real beta;
  real gamma;
  real c0;

  // ---- Toggle for truncation (1 = truncated, 0 = untruncated) ----
  int<lower=0, upper=1> use_truncation;
}

transformed data {
  // Lower bound for finite-endpoint support
  real ystar_lb = fmax(u, ymax) + 1e-6;
  real y_lower  = fmax(ystar_lb, y_min);
  // Finite effective upper bound (needed for parameter declaration)
  real y_upper_eff = fmin(y_cap, is_inf(y_max) ? 1e6 : y_max);

  if (y_upper_eff <= y_lower)
    reject("Invalid bounds: y_upper must exceed y_lower. Check y_cap/y_max vs u,ymax,y_min.");
}

parameters {
  // ---- Endpoint and shape parameters ----
  real<lower=y_lower, upper=y_upper_eff> ystar;
  real<lower=-0.999, upper=-1e-6> xi;   // finite-endpoint (Weibull) case
}

model {
  // ---- Prior for y* ----
  if (use_truncation == 1) {
    // Truncated Normal on [y_lower, y_upper_eff]
    target += normal_lpdf(ystar | mu_y0, sd_mu_y)
              - log_diff_exp(normal_lcdf(y_upper_eff | mu_y0, sd_mu_y),
                             normal_lcdf(y_lower     | mu_y0, sd_mu_y));
  } else {
    // Untruncated Normal
    target += normal_lpdf(ystar | mu_y0, sd_mu_y);
  }

  // ---- Prior for xi (always truncated) ----
  target += normal_lpdf(xi | mu_xi0, sd_xi0)
            - log_diff_exp(normal_lcdf(0   | mu_xi0, sd_xi0),
                           normal_lcdf(-1  | mu_xi0, sd_xi0));

  // ---- Likelihood: finite-endpoint GPD (Weibull case) ----
  {
    real sigma = (u - ystar) * xi;            // > 0 since xi<0 and ystar>u
    vector[N] z = 1 + (y - u) / (u - ystar);  // = (ystar - y)/(ystar - u) in (0,1)
    target += - N * log(sigma) + (-(1/xi + 1)) * sum(log(z));
  }
}

generated quantities {
  real CFHstar_mm = exp(ystar);
  real Mstar_out  = mass_of_ystar(ystar, alpha, beta, gamma, c0);
}
