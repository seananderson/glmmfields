data {
  int<lower=1> nKnots;
  int<lower=1> nLocs;
  int<lower=1> nT;
  int<lower=1> N;
  int<lower=1> stationID[N];
  int<lower=1> yearID[N];
  real y[N];
  real prior_gp_scale[3];
  real prior_gp_sigma[3];
  real prior_sigma[3];
//  real prior_ar[3];
  matrix[nKnots,nKnots] distKnotsSq;
  matrix[nLocs,nKnots] distKnots21Sq;
  int<lower=1> nCov;
  matrix[N,nCov] X;
  int<lower=0,upper=1> gauss_cor;
  int<lower=0,upper=1> est_df;
  int<lower=0,upper=1> obs_model;
  real<lower=2> fixed_df_value;
}
parameters {
  real<lower=0> gp_scale;
  real<lower=0> gp_sigma;
  real<lower=2> df[est_df];
  real<lower=0> sigma;
  real<lower=0> CV;
  vector[nKnots] spatialEffectsKnots[nT];
  vector[nCov] B;
}
transformed parameters {
	vector[nKnots] muZeros;
	vector[nLocs] spatialEffects[nT];
  matrix[nKnots, nKnots] SigmaKnots;
  matrix[nLocs, nKnots] SigmaOffDiag;
  matrix[nLocs, nKnots] invSigmaKnots;
  vector[N] y_hat;
  real<lower=0> gp_sigmaSq;
  real<lower=0> gammaA;

  // allow user to switch between gaussian and exponential covariance
  if (gauss_cor == 1) {
    gp_sigmaSq = gp_sigma^2;
  } else {
    gp_sigmaSq = gp_sigma;
  }

  // cov matrix between knots:
  SigmaKnots = gp_sigmaSq * exp(-gp_scale * distKnotsSq);
  // cov matrix between knots and projected locs:
  SigmaOffDiag = gp_sigmaSq * exp(-gp_scale * distKnots21Sq);
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
	  y_hat[i] = X[i] * B + spatialEffects[yearID[i], stationID[i]];
	}

	if(obs_model==0) {
	  gammaA = 0;
	} else {
	  gammaA = 1/(CV*CV);
	}
}
model {
  // priors:
  gp_scale ~ student_t(prior_gp_scale[1], prior_gp_scale[2], prior_gp_scale[3]);
  gp_sigma ~ student_t(prior_gp_sigma[1], prior_gp_sigma[2], prior_gp_sigma[3]);

  B ~ normal(0, 1);

  // if est_df = 1 estimate degrees of freedom for MVT, otherwise fit MVT with fixed value
  if (est_df == 1) {
    df ~ gamma(2, 0.1);
    for(t in 2:nT) {
      spatialEffectsKnots[t] ~ multi_student_t(df[1], muZeros, SigmaKnots);
    }
  } else {
    for(t in 2:nT) {
      spatialEffectsKnots[t] ~ multi_student_t(fixed_df_value, muZeros, SigmaKnots);
    }
  }

  // switch between normal observation error (1) and gamma (2)
  if(obs_model == 1) {
    sigma ~ student_t(prior_sigma[1], prior_sigma[2], prior_sigma[3]);
    y ~ normal(y_hat, sigma);
  } else {
    // prior on CV of gamma obs error, gamma shape 'a' is derived parameter
    CV ~ student_t(prior_sigma[1], prior_sigma[2], prior_sigma[3]);
    for(i in 1:N) {
       // STAN needs this to be done in loop, argument to gamma can't be vector
       y ~ gamma(gammaA, gammaA/exp(y_hat[i]));
    }
  }

}
