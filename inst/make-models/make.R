
make_stan <- function(obs = c("normal", "gamma")) {
s <- list()

s[[1]] <-
"
data {
  int<lower=1> nKnots;
  int<lower=1> nLocs;
  int<lower=1> nT;
  matrix[nT,nLocs] y;
  matrix[nKnots,nKnots] distKnotsSq;
  matrix[nLocs,nKnots] distKnots21Sq;
"

s[[2]] <- ""

s[[3]] <-
"
}
"

s[[4]] <-
"
parameters {
  real<lower=0> gp_scale;
  real<lower=0> gp_sigmaSq;
  real<lower=0> sigma;
  vector[nKnots] spatialEffectsKnots[nT];
"

s[[5]] <- "
  real<lower=2> scaledf;
"

s[[6]] <- switch(obs[[1]],
  normal = "",
  gamma =
"
  real<lower=0> CV;
")

s[[7]] <-
"
}
"

s[[8]] <-
"
transformed parameters {
  vector[nKnots] muZeros;
  vector[nLocs] spatialEffects[nT];
  cov_matrix[nKnots] SigmaKnots;
  matrix[nLocs,nKnots] SigmaOffDiag;
  matrix[nLocs,nKnots] invSigmaKnots;
"

s[[9]] <- switch(obs[[1]],
  normal = "",
  gamma =
"
real<lower=0> gammaA;
")

s[[10]] <-
"
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
"

s[[11]] <- ""

s[[12]] <-
"
}
"

s[[13]] <-
"
model {
  gp_scale ~ student_t(3, 0, 10);
  gp_sigmaSq ~ student_t(3, 0, 2);
  sigma ~ student_t(3, 0, 2);
  scaledf ~ gamma(2,0.1);
"

s[[14]] <- switch(obs[[1]],
  normal = "",
  gamma =
"
  CV ~ normal(0, 1);
")

s[[15]] <- "
  for(t in 1:nT) {
    spatialEffectsKnots[t] ~ multi_student_t(scaledf, muZeros, SigmaKnots);
  }
"

s[[16]] <- switch(obs[[1]],
  normal =
"
  for(t in 1:nT) {
    y[t] ~ normal(spatialEffects[t], sigma);
  }
",
  gamma =
"
  for(t in 1:nT) {
    for(l in 1:nLocs) {
      y[t,l] ~ gamma(gammaA, gammaA/exp(spatialEffects[t,l]));
    }
  }
}
"
)

s[[17]] <-
"
}
"
s
}

x <- make_stan()
xx <- paste(x, collapse = "")
writeLines(xx, "inst/make-models/test.stan")

x <- make_stan(obs = "gamma")
xx <- paste(x, collapse = "")
writeLines(xx, "inst/make-models/test-gamma.stan")
