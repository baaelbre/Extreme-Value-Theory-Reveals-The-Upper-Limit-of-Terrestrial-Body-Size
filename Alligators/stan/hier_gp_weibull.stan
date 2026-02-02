data {
  // --- Two groups (Sergio = 1, Allan = 2), log-scale data already prepared ---
    int<lower=1> N1;
    vector[N1] y1;          // exceedances (log TL) with y1 > u1
    real u1;                // threshold (log TL) for group 1
    real ymax1;             // max(y1)
    
    int<lower=1> N2;
    vector[N2] y2;          // exceedances (log TL) with y2 > u2
    real u2;                // threshold (log TL) for group 2
    real ymax2;             // max(y2)
    
    // --- Hyperprior hyperparameters (same as your R config) ---
      real mu_y0;
    real<lower=0> sd_y0;
    real mu_xi0;
    real<lower=0> sd_xi0;
    real<lower=0> hc_scale_y;
    real<lower=0> hc_scale_xi;
}
transformed data {
  // Lower bounds for endpoints: must exceed both threshold and sample max
  real ystar1_lb = fmax(u1, ymax1) + 1e-6;
  real ystar2_lb = fmax(u2, ymax2) + 1e-6;
}
parameters {
  // Group-specific endpoints y*_j and shapes xi_j (Weibull domain)
  real<lower=ystar1_lb> ystar1;
  real<lower=ystar2_lb> ystar2;
  
  real<lower=-0.999, upper=-1e-6> xi1;
  real<lower=-0.999, upper=-1e-6> xi2;
  
  // Population (hierarchical) parameters
  real mu_y;
  real<lower=0> tau_y;
  
  real mu_xi;
  real<lower=0> tau_xi;
}
model {
  // ----- Hyperpriors -----
    mu_y  ~ normal(mu_y0,  sd_y0);
  tau_y ~ cauchy(0, hc_scale_y);           // half-Cauchy via <lower=0>
    
    mu_xi  ~ normal(mu_xi0, sd_xi0);
  tau_xi ~ cauchy(0, hc_scale_xi);         // half-Cauchy via <lower=0>
    
    // ----- Hierarchical priors -----
    ystar1 ~ normal(mu_y,  tau_y);
  ystar2 ~ normal(mu_y,  tau_y);
  
  // Truncated-Normal prior on xi in (-1, 0): add normalization term explicitly
  target += normal_lpdf(xi1 | mu_xi, tau_xi)
  - log_diff_exp(normal_lcdf(0 | mu_xi, tau_xi),
                 normal_lcdf(-1 | mu_xi, tau_xi));
  target += normal_lpdf(xi2 | mu_xi, tau_xi)
  - log_diff_exp(normal_lcdf(0 | mu_xi, tau_xi),
                 normal_lcdf(-1 | mu_xi, tau_xi));
  
  // ----- Likelihood: finite-endpoint GP for exceedances (Weibull domain) -----
    {
      // sigma_j = (u_j - y*_j) * xi_j  (positive since both factors < 0)
      real sigma1 = (u1 - ystar1) * xi1;
      real sigma2 = (u2 - ystar2) * xi2;
      
      // Vectorized log-likelihoods
      // z = 1 + xi*(y-u)/sigma = 1 + (y - u)/(u - y*)  (independent of xi)
      vector[N1] z1 = 1 + (y1 - u1) / (u1 - ystar1);
      vector[N2] z2 = 1 + (y2 - u2) / (u2 - ystar2);
      
      // Support is guaranteed by parameter bounds, but guard numerically
      target += - N1 * log(sigma1) + (-(1/xi1 + 1)) * sum(log(z1));
      target += - N2 * log(sigma2) + (-(1/xi2 + 1)) * sum(log(z2));
    }
}
generated quantities {
  // Convenience transforms for plotting/summaries (TL scale in cm)
  real TLstar1_cm = exp(ystar1);
  real TLstar2_cm = exp(ystar2);
  
  // "Population TL*" as exp(mu_y) — matches your original definition
  real TLstar_pop_cm = exp(mu_y);
}
