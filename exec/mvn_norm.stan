# model for presence-absence data
data {
  int<lower=1> nKnots;
  int<lower=1> nLocs;
  int<lower=1> nT;
  matrix[nT,nLocs] y;
  matrix[nKnots,nKnots] distKnotsSq;
  matrix[nLocs,nKnots] distKnots21Sq;
}
parameters {
  real<lower=0> gp_scale;
  real<lower=0> gp_sigmaSq;
  vector[nKnots] spatialEffectsKnots[nT];
}
transformed parameters {
	vector[nKnots] muZeros;
	vector[nLocs] spatialEffects[nT];
  matrix[nKnots, nKnots] SigmaKnots;
  matrix[nLocs,nKnots] SigmaOffDiag;
	for(i in 1:nKnots) {
		muZeros[i] = 0;
	}

  SigmaKnots = gp_sigmaSq * exp(-gp_scale * distKnotsSq);# cov matrix between knots
  SigmaOffDiag = gp_sigmaSq * exp(-gp_scale * distKnots21Sq);# cov matrix between knots and projected
	SigmaOffDiag = SigmaOffDiag * inverse_spd(SigmaKnots); # multiply and invert once, used below

	for(t in 1:nT) {
  spatialEffects[t] = SigmaOffDiag * spatialEffectsKnots[t];
	}
}
model {
  # priors on parameters for covariances, etc
  gp_scale ~ cauchy(0,5);
  gp_sigmaSq ~ cauchy(0,5);

  for(t in 1:nT) {
    spatialEffectsKnots[t] ~ multi_normal(muZeros, SigmaKnots);
  }

  for(t in 1:nT) {
    y[t] ~ normal(spatialEffects[t], 0.1);
  }

}
