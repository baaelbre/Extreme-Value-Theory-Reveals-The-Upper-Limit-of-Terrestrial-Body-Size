functions {
  // Log prior contribution for tau > 0 under different choices.
  real tau_prior_lp(
      real tau,
      int prior_code,
      // half-Cauchy
      real hc_scale,
      // half-Normal
      real hn_scale,
      // half-t_nu
      real ht_nu, real ht_scale,
      // PC prior (exponential on SD)
      real pc_lambda,
      // Log-Normal on SD
      real ln_mu, real ln_sd,
      // Gamma on precision (shape, rate)
      real ga_shape, real ga_rate,
      // Inv-Gamma on variance (shape, scale)
      real ig_shape, real ig_scale
  ) {
    real lp = 0;
    if (prior_code == 1) {
      // Half-Cauchy(0, hc_scale): symmetric Cauchy truncated to tau>0
      // Add log(2) so it's the proper half distribution (constant wrt params otherwise).
      lp = cauchy_lpdf(tau | 0, hc_scale) + log(2);
    } else if (prior_code == 2) {
      // Half-Normal(0, hn_scale)
      lp = normal_lpdf(tau | 0, hn_scale) + log(2);
    } else if (prior_code == 3) {
      // Half-t_nu(0, ht_scale)
      lp = student_t_lpdf(tau | ht_nu, 0, ht_scale) + log(2);
    } else if (prior_code == 4) {
      // PC prior toward tau=0: Exponential(lambda) on SD
      lp = exponential_lpdf(tau | pc_lambda);
    } else if (prior_code == 5) {
      // Log-Normal on SD
      lp = lognormal_lpdf(tau | ln_mu, ln_sd);
    } else if (prior_code == 6) {
      // Gamma on precision phi = 1/tau^2 with Jacobian
      // phi = tau^{-2}, dphi/dtau = -2 * tau^{-3} => |Jac| = 2 / tau^3
      real phi = inv_square(tau);
      lp = gamma_lpdf(phi | ga_shape, ga_rate) + log(2) - 3 * log(tau);
    } else if (prior_code == 7) {
      // Inv-Gamma on variance v = tau^2 with Jacobian
      // v = tau^2, dv/dtau = 2 * tau => log|Jac| = log(2) + log(tau)
      real v = square(tau);
      lp = inv_gamma_lpdf(v | ig_shape, ig_scale) + log(2) + log(tau);
    }
    return lp;
  }
}

data {
  // --- One group: sauropod exceedances on log CFH (y > u) ---
  int<lower=1> N;
  vector[N] y;
  real u;                   // threshold on log CFH
  real ymax;                // max(y) among exceedances

  // --- Hyperprior hyperparameters for means ---
  real mu_y0;
  real<lower=0> sd_y0;
  real mu_xi0;
  real<lower=0> sd_xi0;

  // --- Choice of hyperprior for tau_y and tau_xi ---
  // 1=half-Cauchy, 2=half-Normal, 3=half-t_nu, 4=PC (Exp),
  // 5=LogNormal, 6=Gamma on precision, 7=Inv-Gamma on variance
  int<lower=1, upper=7> tauy_prior;
  int<lower=1, upper=7> tauxi_prior;

  // -------- Hyperparameters for each choice (fill only those you use) --------
  // half-Cauchy scales
  real<lower=0> hc_scale_y;
  real<lower=0> hc_scale_xi;

  // half-Normal scales
  real<lower=0> hn_scale_y;
  real<lower=0> hn_scale_xi;

  // half-t: dof and scales
  real<lower=0> ht_nu_y;
  real<lower=0> ht_scale_y;
  real<lower=0> ht_nu_xi;
  real<lower=0> ht_scale_xi;

  // PC prior (Exponential) calibration via lambda
  // If you prefer (U, alpha), compute lambda=-log(alpha)/U outside and pass here.
  real<lower=0> pc_lambda_y;
  real<lower=0> pc_lambda_xi;

  // Log-Normal parameters for ln(tau)
  real ln_mu_y;
  real<lower=0> ln_sd_y;
  real ln_mu_xi;
  real<lower=0> ln_sd_xi;

  // Gamma on precision phi=1/tau^2  (shape, rate)
  real<lower=0> ga_shape_y;
  real<lower=0> ga_rate_y;
  real<lower=0> ga_shape_xi;
  real<lower=0> ga_rate_xi;

  // Inv-Gamma on variance v = tau^2  (shape, scale)
  real<lower=0> ig_shape_y;
  real<lower=0> ig_scale_y;
  real<lower=0> ig_shape_xi;
  real<lower=0> ig_scale_xi;
}

transformed data {
  // Endpoint must exceed both threshold and sample max
  real ystar_lb = fmax(u, ymax) + 1e-6;
}

parameters {
  // Group-specific endpoint and shape (Weibull domain)
  real<lower=ystar_lb> ystar;
  real<lower=-0.999, upper=-1e-6> xi;

  // Population (hierarchical) parameters
  real mu_y;
  real<lower=0> tau_y;

  real mu_xi;
  real<lower=0> tau_xi;
}

model {
  // ----- Hyperpriors on means -----
  mu_y  ~ normal(mu_y0,  sd_y0);
  mu_xi ~ normal(mu_xi0, sd_xi0);

  // ----- Hyperpriors on spreads (choose per code) -----
  target += tau_prior_lp(
              tau_y, tauy_prior,
              hc_scale_y,
              hn_scale_y,
              ht_nu_y, ht_scale_y,
              pc_lambda_y,
              ln_mu_y, ln_sd_y,
              ga_shape_y, ga_rate_y,
              ig_shape_y, ig_scale_y);

  target += tau_prior_lp(
              tau_xi, tauxi_prior,
              hc_scale_xi,
              hn_scale_xi,
              ht_nu_xi, ht_scale_xi,
              pc_lambda_xi,
              ln_mu_xi, ln_sd_xi,
              ga_shape_xi, ga_rate_xi,
              ig_shape_xi, ig_scale_xi);

  // ----- Hierarchical priors -----
  ystar ~ normal(mu_y,  tau_y);

  // Truncated-Normal prior on xi in (-1, 0): add normalization explicitly
  target += normal_lpdf(xi | mu_xi, tau_xi)
            - log_diff_exp(normal_lcdf(0   | mu_xi, tau_xi),
                           normal_lcdf(-1  | mu_xi, tau_xi));

  // ----- Likelihood: finite-endpoint GP for exceedances (Weibull domain) -----
  {
    // sigma = (u - y*) * xi; (u - y*) < 0 and xi < 0 ⇒ sigma > 0
    real sigma = (u - ystar) * xi;

    // z = 1 + xi*(y-u)/sigma = 1 + (y - u)/(u - y*)  (independent of xi)
    vector[N] z = 1 + (y - u) / (u - ystar);

    // Support ensured by bounds; still guard numerically via construction
    target += - N * log(sigma) + (-(1/xi + 1)) * sum(log(z));
  }
}

generated quantities {
  // Convenience transforms for summaries/plots
  real CFHstar_mm     = exp(ystar);
  real CFHstar_pop_mm = exp(mu_y);  // "population" endpoint scale
}
