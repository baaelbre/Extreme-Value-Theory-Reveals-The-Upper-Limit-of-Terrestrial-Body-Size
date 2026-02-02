
data {
  // Dimensions
  int<lower=1> J;                    // number of traits (e.g., 4: CF, CH, LF, LH)
  int<lower=1> N_spec;               // number of specimens
  int<lower=1> N_obs;                // number of observed bone measurements

  // Observation indices in long format
  array[N_obs] int<lower=1, upper=N_spec> obs_specimen;
  array[N_obs] int<lower=1, upper=J>      obs_trait;
  vector[N_obs]                            y;      // ln X_ij (natural log)

  // Mass information
  array[N_spec] int<lower=0, upper=1> has_mass;    // 1 if lnM observed for specimen i
  vector[N_spec]                         lnM_data; // observed lnM (ignored where has_mass=0)

  // Optional volumetric priors per specimen (on lnM)
  array[N_spec] int<lower=0, upper=1> has_vol_prior;
  vector[N_spec]                         vol_mu;   // mean lnM from volumetric
  vector<lower=0>[N_spec]                vol_sd;   // sd lnM from volumetric

  // Priors for allometry
  vector[J] alpha_hat;       // prior centers for alpha_j
  vector[J] beta_hat;        // prior centers for beta_j
  real<lower=0> sd_alpha;    // e.g., 0.05
  real<lower=0> sd_beta;     // e.g., 0.02

  // Inverse-Gamma prior on sigma_j^2
  real<lower=0> a0;          // e.g., 2.1
  real<lower=0> b0;          // e.g., 0.1

  // Prior for latent lnM when not observed and no volumetric prior
  real mu_M;                 // e.g., 10
  real<lower=0> tau_M;       // e.g., 5

  // Optional prediction targets (for missing bones, etc.)
  int<lower=0> N_pred;                               // 0 allowed
  array[N_pred] int<lower=1, upper=N_spec> pred_specimen;
  array[N_pred] int<lower=1, upper=J>      pred_trait;
}

parameters {
  vector[J] alpha;                         // intercepts
  vector[J] beta;                          // slopes
  vector<lower=0>[J] sigma2;               // residual variances per trait
  vector[N_spec] lnM;                      // ln mass per specimen
}

transformed parameters {
  vector[J] sigma;
  for (j in 1:J) sigma[j] = sqrt(sigma2[j]);
}

model {
  // Priors on trait parameters
  alpha ~ normal(alpha_hat, sd_alpha);
  beta  ~ normal(beta_hat,  sd_beta);
  for (j in 1:J) sigma2[j] ~ inv_gamma(a0, b0);

  // Priors / constraints on lnM per specimen
  for (i in 1:N_spec) {
    if (has_mass[i] == 1) {
      // Treat observed lnM as known (tight normal)
      lnM[i] ~ normal(lnM_data[i], 1e-6);
    } else if (has_vol_prior[i] == 1) {
      // Specimen-specific volumetric prior
      lnM[i] ~ normal(vol_mu[i], vol_sd[i]);
    } else {
      // Default weakly-informative prior
      lnM[i] ~ normal(mu_M, tau_M);
    }
  }

  // Likelihood
  for (n in 1:N_obs) {
    int i = obs_specimen[n];
    int j = obs_trait[n];
    y[n] ~ normal(alpha[j] + beta[j] * lnM[i], sigma[j]);
  }
}

generated quantities {
  // Posterior predictive for observed rows
  vector[N_obs] y_rep;
  for (n in 1:N_obs) {
    int i = obs_specimen[n];
    int j = obs_trait[n];
    y_rep[n] = normal_rng(alpha[j] + beta[j] * lnM[i], sigma[j]);
  }

  // Posterior predictive for requested predictions
  vector[N_pred] y_pred;
  for (k in 1:N_pred) {
    int i = pred_specimen[k];
    int j = pred_trait[k];
    y_pred[k] = normal_rng(alpha[j] + beta[j] * lnM[i], sigma[j]);
  }
}


