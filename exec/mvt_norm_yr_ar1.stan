data {
  int<lower=1> nKnots;
  int<lower=1> nLocs;
  int<lower=1> nT;
  int<lower=1> N;
  int<lower=1> stationID[N];
  int<lower=1> yearID[N];
  real y[N];
  real x[N];
  matrix[nKnots,nKnots] distKnotsSq;
  matrix[nLocs,nKnots] distKnots21Sq;
}
parameters {
  real<lower=0> gp_scale;
  real<lower=0> gp_sigmaSq;
  real<lower=2> scaledf;
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
  matrix[nLocs,nKnots] SigmaOffDiag;
  matrix[nLocs,nKnots] invSigmaKnots;
  vector[N] y_hat;

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
  # priors:
  gp_scale ~ student_t(3, 0, 10);
  gp_sigmaSq ~ student_t(3, 0, 2);
  sigma ~ student_t(3, 0, 2);
  ar ~ normal(0, 1);

  // https://github.com/stan-dev/stan/wiki/Prior-Choice-Recommendations:
  scaledf ~ gamma(2, 0.1);
  year_sigma ~ student_t(3, 0, 2);

  // random walk in year terms:
  yearEffects[1] ~ student_t(3, 0, 2);
  for(t in 2:nT) {
    yearEffects[t] ~ normal(yearEffects[t-1], year_sigma);
  }
  spatialEffectsKnots[1] ~ multi_student_t(scaledf, muZeros, SigmaKnots);
  for(t in 2:nT) {
    spatialEffectsKnots[t] ~ multi_student_t(scaledf,
      ar*spatialEffectsKnots[t-1], SigmaKnots);
  }

  y ~ normal(y_hat, sigma);
}
