functions {
  // ---- Derived mass for reporting (optional) ----
  real mass_of_ystar(real ystar, real alpha, real beta, real gamma, real c0) {
    real z = ystar - c0;
    real m_log = alpha + beta * z + gamma * z * z; // log grams
    return exp(m_log) / 1e6;                       // tons
  }

  // ---- SD prior switch for tau > 0 ----
  real tau_prior_lp(
      real tau,
      int prior_code,
      // 1 half-Cauchy
      real hc_scale,
      // 2 half-Normal
      real hn_scale,
      // 3 half-t (nu, scale)
      real ht_nu, real ht_scale,
      // 4 PC prior (Exponential on SD)
      real pc_lambda,
      // 5 Log-Normal on SD
      real ln_mu, real ln_sd,
      // 6 Gamma on precision phi = 1/tau^2  (shape, rate)
      real ga_shape, real ga_rate,
      // 7 Inv-Gamma on variance v = tau^2   (shape, scale)
      real ig_shape, real ig_scale
  ) {
    real lp = 0;
    if (prior_code == 1) {
      lp = cauchy_lpdf(tau | 0, hc_scale) + log(2);          // half-Cauchy
    } else if (prior_code == 2) {
      lp = normal_lpdf(tau | 0, hn_scale) + log(2);          // half-Normal
    } else if (prior_code == 3) {
      lp = student_t_lpdf(tau | ht_nu, 0, ht_scale) + log(2);// half-t
    } else if (prior_code == 4) {
      lp = exponential_lpdf(tau | pc_lambda);                // PC(Exp on SD)
    } else if (prior_code == 5) {
      lp = lognormal_lpdf(tau | ln_mu, ln_sd);               // LogNormal on SD
    } else if (prior_code == 6) {
      real phi = inv_square(tau);                            // precision
      lp = gamma_lpdf(phi | ga_shape, ga_rate) + log(2) - 3*log(tau);
    } else if (prior_code == 7) {
      real v = square(tau);                                  // variance
      lp = inv_gamma_lpdf(v | ig_shape, ig_scale) + log(2) + log(tau);
    }
    return lp;
  }
}

data {
  // ---- Exceedances ----
  int<lower=1> N;
  vector[N] y;
  real u;
  real ymax;                 // sample max among exceedances

  // ---- Prior support for endpoint ----
  real y_min;                // user lower cap for y*
  real y_max;                // legacy/user cap (kept for safety)
  real y_cap;                // hard cap for y* prior (set to 8.62 for 1000 t)

  // ---- Hyperprior on mu_y*: Normal(mu_y0, sd_mu_y) ----
  real mu_y0;                // set to 8.15
  real<lower=0> sd_mu_y;     // choose e.g. 0.3–0.6

  // ---- Hyperpriors for xi location ----
  real mu_xi0;
  real<lower=0> sd_xi0;

  // ---- Choice of SD priors ----
  int<lower=1, upper=7> tauy_prior;
  int<lower=1, upper=7> tauxi_prior;

  // ---- Hyperparams for SD priors ----
  real<lower=0> hc_scale_y;
  real<lower=0> hc_scale_xi;

  real<lower=0> hn_scale_y;
  real<lower=0> hn_scale_xi;

  real<lower=0> ht_nu_y;
  real<lower=0> ht_scale_y;
  real<lower=0> ht_nu_xi;
  real<lower=0> ht_scale_xi;

  real<lower=0> pc_lambda_y;
  real<lower=0> pc_lambda_xi;

  real ln_mu_y;
  real<lower=0> ln_sd_y;
  real ln_mu_xi;
  real<lower=0> ln_sd_xi;

  real<lower=0> ga_shape_y;
  real<lower=0> ga_rate_y;
  real<lower=0> ga_shape_xi;
  real<lower=0> ga_rate_xi;

  real<lower=0> ig_shape_y;
  real<lower=0> ig_scale_y;
  real<lower=0> ig_shape_xi;
  real<lower=0> ig_scale_xi;

  // ---- Allometry for reporting (optional) ----
  real alpha;
  real beta;
  real gamma;
  real c0;
}

transformed data {
  // EVT support: y* must exceed both threshold and sample max
  real ystar_lb = fmax(u, ymax) + 1e-6;
  // Effective lower & upper bounds for y*
  real y_lower  = fmax(ystar_lb, y_min);
  real y_upper  = fmin(y_cap, y_max);  // typically y_cap = 8.62

  // quick safety checks on bounds
  if (y_upper <= y_lower)
    reject("Invalid bounds: y_upper must exceed y_lower. Check y_cap/y_max vs u,ymax,y_min.");
}

parameters {
  // ---- Hyperparameters ----
  real mu_y;                 // now ~ Normal(mu_y0, sd_mu_y)
  real<lower=0> tau_y;

  real mu_xi;                // xi location
  real<lower=0> tau_xi;      // xi spread

  // ---- Dataset-level ----
  real<lower=y_lower, upper=y_upper> ystar;
  real<lower=-0.999, upper=-1e-6> xi;
}

model {
  // ---- Hyperpriors ----
  mu_y  ~ normal(mu_y0, sd_mu_y); // centered at 8.15
  target += tau_prior_lp(tau_y,  tauy_prior,
                         hc_scale_y, hn_scale_y, ht_nu_y, ht_scale_y,
                         pc_lambda_y, ln_mu_y, ln_sd_y,
                         ga_shape_y, ga_rate_y, ig_shape_y, ig_scale_y);

  mu_xi ~ normal(mu_xi0, sd_xi0);
  target += tau_prior_lp(tau_xi, tauxi_prior,
                         hc_scale_xi, hn_scale_xi, ht_nu_xi, ht_scale_xi,
                         pc_lambda_xi, ln_mu_xi, ln_sd_xi,
                         ga_shape_xi, ga_rate_xi, ig_shape_xi, ig_scale_xi);

  // ---- Priors ----
  // y* ~ Normal(mu_y, tau_y) truncated to [y_lower, y_upper]
  target += normal_lpdf(ystar | mu_y, tau_y)
            - log_diff_exp(normal_lcdf(y_upper | mu_y, tau_y),
                           normal_lcdf(y_lower | mu_y, tau_y));

  // xi ~ Normal(mu_xi, tau_xi) truncated to (-1, 0)
  target += normal_lpdf(xi | mu_xi, tau_xi)
            - log_diff_exp(normal_lcdf(0   | mu_xi, tau_xi),
                           normal_lcdf(-1  | mu_xi, tau_xi));

  // ---- Likelihood: finite-endpoint GPD (Weibull) ----
  {
    real sigma = (u - ystar) * xi;            // > 0 by construction (xi<0, ystar>u)
    vector[N] z = 1 + (y - u) / (u - ystar);  // = (ystar - y)/(ystar - u) in (0,1)
    target += - N * log(sigma) + (-(1/xi + 1)) * sum(log(z));
  }
}

generated quantities {
  real CFHstar_mm = exp(ystar);
  real Mstar_out  = mass_of_ystar(ystar, alpha, beta, gamma, c0);
}
