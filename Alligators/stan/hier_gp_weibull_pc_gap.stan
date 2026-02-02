functions {
  real gp_weibull_logpdf(vector y, real u, real ystar, real xi) {
    int n = num_elements(y);
    real denom = ystar - u;            // > 0
    vector[n] t;
    for (i in 1:n) {
      if (y[i] >= ystar) return negative_infinity();
      t[i] = (ystar - y[i]) / denom;   // in (0,1]
    }
    real sigma = (u - ystar) * xi;     // > 0 because xi<0 and u<ystar
    return - n * log(sigma) - (1/xi + 1) * sum(log(t));
  }
}
data {
  // Dataset 1
  int<lower=1> N1;
  vector[N1] y1;
  real u1;
  real ymax1;

  // Dataset 2
  int<lower=1> N2;
  vector[N2] y2;
  real u2;
  real ymax2;

  // PC hyper: lambda ~ Gamma(a_lambda, b_lambda)
  real<lower=0> a_lambda;
  real<lower=0> b_lambda;

  // Shape hierarchy (unchanged)
  real mu_xi0;
  real<lower=0> sd_xi0;
  real<lower=0> hc_scale_xi;

  // For population-predictive endpoint (anchor)
  real y_max_pool;
}
parameters {
  // Endpoint gaps (PC prior)
  real<lower=0> delta1;
  real<lower=0> delta2;

  // Shapes (Weibull domain)
  real<lower=-1, upper=0> xi1;
  real<lower=-1, upper=0> xi2;
  real mu_xi;
  real<lower=0> tau_xi;

  // Shared PC rate
  real<lower=0> lambda;
}
transformed parameters {
  real ystar1 = ymax1 + delta1;
  real ystar2 = ymax2 + delta2;
}
model {
  // Hyperprior on rate (Gamma shape-rate)
  lambda ~ gamma(a_lambda, b_lambda);

  // PC priors on gaps
  delta1 ~ exponential(lambda);
  delta2 ~ exponential(lambda);

  // Shape hierarchy
  mu_xi  ~ normal(mu_xi0, sd_xi0);
  tau_xi ~ cauchy(0, hc_scale_xi);               // half-Cauchy via <lower=0>
  xi1    ~ normal(mu_xi, tau_xi) T[-1, 0];
  xi2    ~ normal(mu_xi, tau_xi) T[-1, 0];

  // Likelihoods
  target += gp_weibull_logpdf(y1, u1, ystar1, xi1);
  target += gp_weibull_logpdf(y2, u2, ystar2, xi2);
}
generated quantities {
  // Dataset endpoints (cm)
  real TLstar1_cm = exp(ystar1);
  real TLstar2_cm = exp(ystar2);

  // Posterior-predictive population endpoint gap (same lambda)
  real delta_pred = exponential_rng(lambda);
  real ystar_pop  = y_max_pool + delta_pred;
  real TLstar_pop_cm = exp(ystar_pop);
}
