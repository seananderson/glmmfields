
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
  real<lower=0> sigma;
  vector[nKnots] spatialEffectsKnots[nT];

  real<lower=2> scaledf;

  real<lower=0> CV;

}

transformed parameters {
  vector[nKnots] muZeros;
  vector[nLocs] spatialEffects[nT];
  cov_matrix[nKnots] SigmaKnots;
  matrix[nLocs,nKnots] SigmaOffDiag;
  matrix[nLocs,nKnots] invSigmaKnots;

real<lower=0> gammaA;

  // cov matrix between knots:
  SigmaKnots = gp_sigmaSq * exp(-gp_scale * distKnotsSq);

  // cov matrix between knots and projected locs:
  SigmaOffDiag = gp_sigmaSq * exp(-gp_scale * distKnots21Sq);

  for(i in 1:nKnots) {
    muZeros[i] = 0;
  }

  // multiply and invert   once, used below
  SigmaOffDiag = SigmaOffDiag * inverse_spd(SigmaKnots);

  for(i in 1:nT) {
    spatialEffects[i] = SigmaOffDiag * spatialEffectsKnots[i];
  }

}

model {
  gp_scale ~ student_t(3, 0, 10);
  gp_sigmaSq ~ student_t(3, 0, 2);
  sigma ~ student_t(3, 0, 2);
  scaledf ~ gamma(2,0.1);

  CV ~ normal(0, 1);

  for(t in 1:nT) {
    spatialEffectsKnots[t] ~ multi_student_t(scaledf, muZeros, SigmaKnots);
  }

  for(t in 1:nT) {
    for(l in 1:nLocs) {
      y[t,l] ~ gamma(gammaA, gammaA/exp(spatialEffects[t,l]));
    }
  }
}

}

