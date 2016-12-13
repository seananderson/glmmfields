# model for presence-absence data
data {
  int<lower=1> nKnots;
  int<lower=1> nLocs;
  int<lower=1> nT;
  int<lower=1> N;
  int<lower=1> stationID[N];
  int<lower=1> yearID[N];
  real y[N];
  real x[N];
  #matrix[nT,nLocs] y;
  matrix[nKnots,nKnots] distKnotsSq;
  matrix[nLocs,nKnots] distKnots21Sq;
}
parameters {
  real<lower=0> gp_scale;
  real<lower=0> gp_sigmaSq;
  real<lower=1> scaledf;
  real<lower=0> CV;
  real B;
  real ar;
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
  real gammaA;
  SigmaKnots = gp_sigmaSq * exp(-gp_scale * distKnotsSq);# cov matrix between knots
  SigmaOffDiag = gp_sigmaSq * exp(-gp_scale * distKnots21Sq);# cov matrix between knots and projected locs
	for(i in 1:nKnots) {
		muZeros[i] = 0;
	}
	SigmaOffDiag = SigmaOffDiag * inverse_spd(SigmaKnots); # multiply and invert once, used below
	for(i in 1:nT) {
  spatialEffects[i] = SigmaOffDiag * spatialEffectsKnots[i];
	}
	gammaA = 1/(CV*CV);
}
model {
  # priors on parameters for covariances, etc
  gp_scale ~ student_t(3, 0, 20);
  gp_sigmaSq ~ student_t(3, 0, 2);
  CV ~ normal(0,1);#student_t(3, 0, 2);
  #2 * (ar - 0.5) ~ beta(2, 2);
  ar ~ normal(0, 1);
  scaledf ~ gamma(2, 0.1); # https://github.com/stan-dev/stan/wiki/Prior-Choice-Recommendations
  year_sigma ~ student_t(3, 0, 2);

  # random walk / random effects in year terms
  yearEffects[1] ~ student_t(3, 0, 2);
  for(t in 2:nT) {
    yearEffects[t] ~ normal(yearEffects[t-1],year_sigma);
  }

  spatialEffectsKnots[1] ~ multi_student_t(scaledf, muZeros, SigmaKnots);
  for(t in 2:nT) {
  spatialEffectsKnots[t] ~ multi_student_t(scaledf, ar*spatialEffectsKnots[t-1], SigmaKnots);
  }

  # not sure this can be efficiently vectorized. See: http://stackoverflow.com/questions/35271398/simple-gamma-glm-in-stan
  for(i in 1:N) {
    y[i] ~ gamma(gammaA, gammaA/exp(B*x[i] + yearEffects[yearID[i]] + spatialEffects[yearID[i],stationID[i]]));
  }
}

