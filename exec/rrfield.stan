data {
  int<lower=1> nKnots;
  int<lower=1> nLocs;
  int<lower=1> nT;
  int<lower=1> N;
  int<lower=1> stationID[N];
  int<lower=1> yearID[N];
  real y[N]; // y for normal and gamma obs. model
  int y_int[N]; // y for NB2 obs. model
  real prior_gp_scale[3];
  real prior_gp_sigma[3];
  real prior_sigma[3];
  real prior_rw_sigma[3];
  real prior_intercept[3];
  real prior_beta[3];
  real prior_ar[3];
  matrix[nKnots,nKnots] distKnots;
  matrix[nLocs,nKnots] distKnots21;
  int<lower=0> nCov;
  matrix[N,nCov] X;
  int<lower=0,upper=1> cov_func; # 0 exp, 1 = sq_exp, 2 = matern
  int<lower=0,upper=1> est_df;
  int<lower=0,upper=1> est_ar;
  int<lower=0,upper=1> norm_params;
  int<lower=0,upper=1> gamma_params;
  int<lower=0,upper=1> nb2_params;
  int<lower=0,upper=6> obs_model;
  real<lower=1> fixed_df_value;
  real fixed_ar_value;
  int<lower=0,upper=1> est_temporalRE;
  int<lower=0> n_year_effects;
  int<lower=0> lower_truncation;
  int<lower=0,upper=1> fixed_intercept;
  real matern_kappa;
}
parameters {
  real<lower=0> gp_scale;
  real<lower=0> gp_sigma;
  real<lower=2> df[est_df];
  real<lower=0> sigma[norm_params];
  real<lower=0> CV[gamma_params];
  real<lower=0> nb2_phi[nb2_params];
  real yearEffects[n_year_effects];
  real<lower=0> year_sigma[est_temporalRE];
  vector[nKnots] spatialEffectsKnots[nT];
  vector[nCov] B;
  real<lower=-1, upper=1> ar[est_ar];
  real<lower=0> W[nT];
}
transformed parameters {
  vector[nKnots] muZeros;
  vector[nLocs] spatialEffects[nT];
  matrix[nKnots, nKnots] SigmaKnots;
  matrix[nKnots, nKnots] transformed_dist;
  matrix[nLocs, nKnots] transformed_dist21;
  matrix[nLocs, nKnots] SigmaOffDiag;
  matrix[nLocs, nKnots] invSigmaKnots;
  vector[N] y_hat;
  real<lower=0> gammaA[gamma_params];
  real<lower=0> gp_sigma_sq;
  gp_sigma_sq = pow(gp_sigma, 2.0);

  // allow user to switch between squared-exponential and exponential covariance
  if(cov_func == 0) {
    // cov matrix between knots
    SigmaKnots = gp_sigma_sq * exp(-distKnots / gp_scale);
    // cov matrix between knots and projected locs
    SigmaOffDiag = gp_sigma_sq * exp(-distKnots21 / gp_scale);
  }
  if (cov_func == 1) {
    // cov matrix between knots:
    SigmaKnots = gp_sigma_sq *
      exp(-inv(2.0 * pow(gp_scale, 2.0)) * distKnots); // dist^2 as data
    // cov matrix between knots and projected locs:
    SigmaOffDiag = gp_sigma_sq *
      exp(-inv(2.0 * pow(gp_scale, 2.0)) * distKnots21); // dist^2 as data
  }
  if (cov_func == 2) {
    if(matern_kappa == 1.5) {
      // cov matrix between knots
      transformed_dist = sqrt(3) * distKnots / gp_scale;
      SigmaKnots = gp_sigma_sq * (1 + transformed_dist) * exp (-transformed_dist);
      // cov matrix between knots and projected locs
      transformed_dist21 = sqrt(3) * distKnots21 / gp_scale;
      SigmaOffDiag = gp_sigma_sq * (1 + transformed_dist21) * exp (-transformed_dist21);
    }
    if(matern_kappa == 2.5) {
      // cov matrix between knots
      transformed_dist = sqrt(5) * distKnots / gp_scale;
      SigmaKnots = gp_sigma_sq * (1 + transformed_dist + (transformed_dist .* transformed_dist)/3) * exp (-transformed_dist);
      // cov matrix between knots and projected locs
      transformed_dist21 = sqrt(5) * distKnots21 / gp_scale;
      SigmaOffDiag = gp_sigma_sq * (1 + transformed_dist21 + (transformed_dist21 .* transformed_dist21)/3) * exp (-transformed_dist21);
    }
  }

  for(k in 1:nKnots) {
    muZeros[k] = 0;
  }
  // multiply and invert once, used below:
  SigmaOffDiag = SigmaOffDiag * inverse_spd(SigmaKnots);
  for(t in 1:nT) {
    spatialEffects[t] = SigmaOffDiag * spatialEffectsKnots[t];
  }

  // calculate predicted value of each observation
  for(i in 1:N) {
    if(est_temporalRE == 0) {
      if(fixed_intercept == 0) {
        y_hat[i] = X[i] * B + spatialEffects[yearID[i], stationID[i]];
      } else {
        y_hat[i] = spatialEffects[yearID[i], stationID[i]];
      }
    } else {
      y_hat[i] = spatialEffects[yearID[i], stationID[i]] + yearEffects[yearID[i]];
    }
  }

  if(obs_model==0) {
    gammaA[1] = inv(pow(CV[1], 2.0));
  }

}
model {
  // priors:
  gp_scale ~ student_t(prior_gp_scale[1], prior_gp_scale[2], prior_gp_scale[3]);
  gp_sigma ~ student_t(prior_gp_sigma[1], prior_gp_sigma[2], prior_gp_sigma[3]);

  if(est_ar == 1) {
    ar ~ student_t(prior_ar[1], prior_ar[2], prior_ar[3]);
  }

  if(nCov >= 1) {
    // global intercept, absorbed into year re [1] if those estimated
    B[1] ~ student_t(prior_intercept[1], prior_intercept[2], prior_intercept[3]);
  }
  if (nCov >= 2) {
    for (i in 2:nCov) {
      // coefficients associated with non-intercept covariates
      B[i] ~ student_t(prior_beta[1], prior_beta[2], prior_beta[3]);
    }
  }

  // temporal random effects, if estimated global intercept = effect in first year
  if(est_temporalRE == 1) {
    year_sigma ~ student_t(prior_rw_sigma[1], prior_rw_sigma[2], prior_rw_sigma[3]);
    // random walk in year terms
    yearEffects[1] ~ student_t(prior_intercept[1], prior_intercept[2], prior_intercept[3]);
    for(t in 2:nT) {
      yearEffects[t] ~ normal(yearEffects[t-1], year_sigma);
    }
  }

  // if est_df == 1 estimate MVT degrees of freedom, otherwise use fixed df
  if (est_df == 1) {
    W ~ scaled_inv_chi_square(df[1], 1);
    df ~ gamma(2, 0.1);
  } else {
    W ~ scaled_inv_chi_square(fixed_df_value, 1);
  }

  // spatial deviates in first time slice
  spatialEffectsKnots[1] ~ multi_normal(muZeros, W[1]*SigmaKnots);

  // spatial deviates in remaining time slices
  for(t in 2:nT) {
    if(est_ar == 1) {
      spatialEffectsKnots[t] ~ multi_normal(ar[1] * spatialEffectsKnots[t-1],
          W[t]*SigmaKnots);
    } else {
      spatialEffectsKnots[t] ~ multi_normal(fixed_ar_value * spatialEffectsKnots[t-1],
          W[t]*SigmaKnots);
    }
  }

  // switch between observation error models: normal (1), gamma (0), NB2 (2), binomial (4), poisson (5)
  // lognormal (6)
  // (tweedie (3) is in a branch)
  if(obs_model == 5) {
    y_int ~ poisson_log(y_hat);
  }
  if(obs_model == 4) {
    y_int ~ bernoulli_logit(y_hat);
  }
  if(obs_model == 2) {
    nb2_phi[1] ~ student_t(prior_sigma[1], prior_sigma[2], prior_sigma[3]);
    if (lower_truncation == 0) {
      y_int ~ neg_binomial_2_log(y_hat, nb2_phi[1]);
    } else {
      for (i in 1:N) {
        y_int[i] ~ neg_binomial_2(exp(y_hat[i]), nb2_phi[1]) T[lower_truncation, ];
      }
    }
  }
  if(obs_model == 1) {
    sigma[1] ~ student_t(prior_sigma[1], prior_sigma[2], prior_sigma[3]);
    y ~ normal(y_hat, sigma[1]);
  }
  if(obs_model == 0) {
    // prior on CV of gamma obs error, gamma shape 'a' is derived parameter
    CV[1] ~ student_t(prior_sigma[1], prior_sigma[2], prior_sigma[3]);
    y ~ gamma(gammaA[1], gammaA[1] ./ exp(y_hat));
  }
  if(obs_model == 6) {
    sigma[1] ~ student_t(prior_sigma[1], prior_sigma[2], prior_sigma[3]);
    y ~ lognormal(y_hat, sigma[1]);
  }
}
generated quantities {
  vector[N] log_lik;

  for (i in 1:N) {
    if(obs_model == 0) {
      log_lik[i] = gamma_lpdf(y[i] | gammaA[1], gammaA[1] ./ exp(y_hat[i]));
    }
    if(obs_model == 1) {
      log_lik[i] = normal_lpdf(y[i] | y_hat[i], sigma[1]);
    }
    if(obs_model == 2) {
      if (lower_truncation == 0) {
        log_lik[i] = neg_binomial_2_log_lpmf(y_int[i] | y_hat[i], nb2_phi[1]);
      } else {
        // TODO!!! Note that I had to remove T[lower_truncation, ] from the following line
        // and I think that will make this calculation incorrect
        // in the case of truncated negative binomial
        // the package will issue a warning
        log_lik[i] = neg_binomial_2_lpmf(y_int[i] | exp(y_hat[i]), nb2_phi[1]);
      }
    }
    if(obs_model == 4) {
      log_lik[i] = bernoulli_logit_lpmf(y_int[i] | y_hat[i]);
    }
    if(obs_model == 5) {
      log_lik[i] = poisson_log_lpmf(y_int[i] | y_hat[i]);
    }
    if(obs_model == 6) {
      log_lik[i] = lognormal_lpdf(y[i] | y_hat, sigma[1]);
    }
  }
}
