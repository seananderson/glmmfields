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
  real prior_ar[3];
  matrix[nKnots,nKnots] distKnotsSq;
  matrix[nLocs,nKnots] distKnots21Sq;
  int<lower=1> nCov;
  matrix[N,nCov] X;
}
parameters {
  real<lower=0> gp_scale;
  real<lower=0> gp_sigma;
  real<lower=2> df;
  real<lower=0> sigma;
  real yearEffects[nT];
  real<lower=0> year_sigma;
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
  // flag for gaussian
  gp_sigmaSq = gp_sigma^2;

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
	for(i in 1:N) {
	  y_hat[i] = X[i] * B + spatialEffects[yearID[i], stationID[i]];
	}
}
model {
  // priors:
  gp_scale ~ student_t(prior_gp_scale[1], prior_gp_scale[2], prior_gp_scale[3]);
  gp_sigma ~ student_t(prior_gp_sigma[1], prior_gp_sigma[2], prior_gp_sigma[3]);
  sigma ~ student_t(prior_sigma[1], prior_sigma[2], prior_sigma[3]);
  df ~ gamma(2, 0.1);

  // regression coefficients:
  B ~ normal(0, 1);

  for(t in 2:nT) {
    spatialEffectsKnots[t] ~ multi_student_t(df, muZeros, SigmaKnots);
  }

  y ~ normal(y_hat, sigma);
}
