// gpd_endpoint.stan
data {
  int<lower=1> N;
  real<lower=0> y[N];         // Exceedances: y = log(SVL) - threshold
  real<lower=0> thresh;       // Threshold (log scale)
  real<lower=0> y_max;        // max(log_SVL)
}

parameters {
  real<lower=0> y_star;       // Endpoint of log SVL
  real psi;                   // Transformed xi: xi = -exp(psi)
}

transformed parameters {
  real xi;
  xi = -exp(psi);
}

model {
  // Impose constraint: y_star must be greater than max(log_SVL)
  if (y_star <= y_max) {
    reject("y_star must be greater than the largest data point");
  }
  
  // GPD likelihood (in log-SVL space)
  for (n in 1:N) {
    real z = y_star - y[n];
    target += log(y_star - thresh) / xi;
    target += -(1 / xi + 1) * log(z);
    target += -log(-xi);
  }
}
