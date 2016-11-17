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
}
parameters {
  real<lower=0> gp_scale;
  real<lower=0> gp_sigma;
  real<lower=2> df;
  real<lower=0> sigma;
  real<lower=-1,upper=1> ar;
  real yearEffects[nT];
  real<lower=0> year_sigma;
  vector[nKnots] spatialEffectsKnots[nT];
}
transformed parameters {
	vector[nKnots] muZeros;
	vector[nLocs] spatialEffects[nT];
  matrix[nKnots, nKnots] SigmaKnots;
  matrix[nLocs, nKnots] SigmaOffDiag;
  matrix[nLocs, nKnots] invSigmaKnots;
  vector[N] y_hat;
  real<lower=0> gp_sigmaSq;

  gp_sigmaSq = gp_sigma^2;

  // cov matrix between knots:
  SigmaKnots = gp_sigmaSq * exp(-gp_scale * distKnotsSq);
  // cov matrix between knots and projected locs:
  SigmaOffDiag = gp_sigmaSq * exp(-gp_scale * distKnots21Sq);
	for(i in 1:nKnots) {
		muZeros[i] = 0;
	}
	// multiply and invert once, used below:
	SigmaOffDiag = SigmaOffDiag * inverse_spd(SigmaKnots);
	for(i in 1:nT) {
    spatialEffects[i] = SigmaOffDiag * spatialEffectsKnots[i];
	}
	for(i in 1:N) {
	  y_hat[i] = yearEffects[yearID[i]] + spatialEffects[yearID[i],stationID[i]];
	}
}
model {
  // priors:
  gp_scale ~ student_t(prior_gp_scale[1], prior_gp_scale[2], prior_gp_scale[3]);
  gp_sigma ~ student_t(prior_gp_sigma[1], prior_gp_sigma[2], prior_gp_sigma[3]);
  sigma ~ student_t(prior_sigma[1], prior_sigma[2], prior_sigma[3]);
  ar ~ student_t(prior_ar[1], prior_ar[2], prior_ar[3]);
  df ~ gamma(2, 0.1);
  year_sigma ~ student_t(3, 0, 2);

  // random walk in year terms:
  yearEffects[1] ~ student_t(3, 0, 2);
  for(t in 2:nT) {
    yearEffects[t] ~ normal(yearEffects[t-1], year_sigma);
  }
  spatialEffectsKnots[1] ~ multi_student_t(df, muZeros, SigmaKnots);
  for(t in 2:nT) {
    spatialEffectsKnots[t] ~ multi_student_t(df,
      ar*spatialEffectsKnots[t-1], SigmaKnots);
  }

  y ~ normal(y_hat, sigma);
}
