// from https://github.com/LuZhangstat/NNGP_STAN/blob/master/nngp_latent.stan
functions{
  real nngp_w_lpdf(vector w, real sigmasq, real gp_theta, matrix NN_dist,
                   matrix NN_distM, int[,] NN_ind, int N, int M,
                   int cov_func){

    vector[N] V;
    vector[N] I_Aw = w;
    int dim;
    int h;

    for (i in 2:N) {

      matrix[ i < (M + 1)? (i - 1) : M, i < (M + 1)? (i - 1): M]
      iNNdistM;
      matrix[ i < (M + 1)? (i - 1) : M, i < (M + 1)? (i - 1): M]
      iNNCholL;
      vector[ i < (M + 1)? (i - 1) : M] iNNcorr;
      vector[ i < (M + 1)? (i - 1) : M] v;
      row_vector[i < (M + 1)? (i - 1) : M] v2;

      dim = (i < (M + 1))? (i - 1) : M;

      if(dim == 1){iNNdistM[1, 1] = 1;}
      else{
        h = 0;
        for (j in 1:(dim - 1)){
          for (k in (j + 1):dim){
            h = h + 1;
            if(cov_func==0) {
              iNNdistM[j, k] = exp(- NN_distM[(i - 1), h] / gp_theta);
            }
            if(cov_func==1) {
              iNNdistM[j, k] = exp(-inv(2.0 * pow(gp_theta, 2.0)) *  NN_distM[(i - 1), h]);
            }
            iNNdistM[k, j] = iNNdistM[j, k];
          }
        }
        for(j in 1:dim){
          iNNdistM[j, j] = 1;
        }
      }

      iNNCholL = cholesky_decompose(iNNdistM);
      if(cov_func==0) {
        iNNcorr = to_vector(exp(- NN_dist[(i - 1), 1:dim] / gp_theta));
      }
      if(cov_func==1) {
        iNNcorr = to_vector(exp(-inv(2.0 * pow(gp_theta, 2.0)) * NN_dist[(i - 1), 1:dim]));
      }
      v = mdivide_left_tri_low(iNNCholL, iNNcorr);

      V[i] = 1 - dot_self(v);

      v2 = mdivide_right_tri_low(v', iNNCholL);

                                 I_Aw[i] = I_Aw[i] - v2 * w[NN_ind[(i - 1), 1:dim]];

    }
                                 V[1] = 1;
                                 return - 0.5 * ( 1 / sigmasq * dot_product(I_Aw, (I_Aw ./ V)) +
                                 sum(log(V)) + N * log(sigmasq));
  }
}


data {
  int<lower=1> N;
  int<lower=0> nCov;
  real y[N]; // y for normal and gamma obs. model
  int y_int[N]; // y for NB2 or poisson or binomial obs. model
  matrix[N, nCov] X;
  int<lower=0, upper=1> norm_params;
  int<lower=0, upper=1> gamma_params;
  int<lower=0, upper=1> nb2_params;
  int<lower=0, upper=6> obs_model;
  int<lower=0> lower_truncation;
  real prior_sigma[3];
  real prior_intercept[3];
  real prior_beta[3];
  real prior_gp_theta[3];
  real prior_gp_sigma[3];
  int<lower=0, upper=2> cov_func; // 0 = exp, 1 = sq_exp, 2 = matern

  int<lower=1> M;
  int NN_ind[N - 1, M];
  matrix[N - 1, M] NN_dist;
  matrix[N - 1, (M * (M - 1) / 2)] NN_distM;
  real ss;
  real st;
  real ap;
  real bp;
}
parameters{
  vector[nCov] B;
  vector[N] w;
  real<lower=0> sigma[norm_params];
  real<lower=0> CV[gamma_params];
  real<lower=0> nb2_phi[nb2_params];
  real<lower = 0> gp_theta;
  real<lower = 0> gp_sigma;

  real<lower = 0> sigmasq;
  real<lower = 0> tau;

}
transformed parameters {
  //real sigmasq = square(sigma);
  vector[N] y_hat; // vector to hold predicted values
  real tausq = square(tau);
  real<lower=0> gammaA[gamma_params];
    if (obs_model==0) {
    gammaA[1] = inv(pow(CV[1], 2.0));
  }
  // predicted
  y_hat = X * B + w;
}
model{
  B ~ normal(0,1);//multi_normal_cholesky(uB, L_VB);
  tau ~ normal(0, st);
  w ~ nngp_w(gp_sigma^2, gp_theta, NN_dist, NN_distM, NN_ind, N, M, cov_func);

  // priors:
  gp_theta ~ student_t(prior_gp_theta[1], prior_gp_theta[2], prior_gp_theta[3]);
  gp_sigma ~ student_t(prior_gp_sigma[1], prior_gp_sigma[2], prior_gp_sigma[3]);

  if (nCov >= 1) {
    // global intercept, absorbed into year re [1] if those estimated
    B[1] ~ student_t(prior_intercept[1], prior_intercept[2], prior_intercept[3]);
  }
  if (nCov >= 2) {
    for (i in 2:nCov) {
      // coefficients associated with non-intercept covariates
      B[i] ~ student_t(prior_beta[1], prior_beta[2], prior_beta[3]);
    }
  }

  // switch between observation error models:
  // gamma (0), normal (1), NB2 (2), binomial (4), poisson (5), lognormal (6)
  // where is 3? tweedie (3) is in a branch and is too slow to be practical
  if (obs_model == 0) {
    // prior on CV of gamma obs error, gamma shape 'a' is derived parameter
    CV[1] ~ student_t(prior_sigma[1], prior_sigma[2], prior_sigma[3]);
    y ~ gamma(gammaA[1], gammaA[1] ./ exp(y_hat));
  }
  if (obs_model == 1) {
    sigma[1] ~ student_t(prior_sigma[1], prior_sigma[2], prior_sigma[3]);
    y ~ normal(y_hat, sigma[1]);
  }
  if (obs_model == 2) {
    nb2_phi[1] ~ student_t(prior_sigma[1], prior_sigma[2], prior_sigma[3]);
    if (lower_truncation == 0) {
      y_int ~ neg_binomial_2_log(y_hat, nb2_phi[1]);
    } else {
      for (i in 1:N) {
        y_int[i] ~ neg_binomial_2(exp(y_hat[i]), nb2_phi[1]) T[lower_truncation, ];
      }
    }
  }
  if (obs_model == 4) {
    y_int ~ bernoulli_logit(y_hat);
  }
  if (obs_model == 5) {
    y_int ~ poisson_log(y_hat);
  }
  if (obs_model == 6) {
    sigma[1] ~ student_t(prior_sigma[1], prior_sigma[2], prior_sigma[3]);
    y ~ lognormal(y_hat, sigma[1]);
  }

}
generated quantities {
  // log_lik is for use with the loo package
  vector[N] log_lik;

  for (i in 1:N) {
    if (obs_model == 0) {
      log_lik[i] = gamma_lpdf(y[i] | gammaA[1], gammaA[1] ./ exp(y_hat[i]));
    }
    if (obs_model == 1) {
      log_lik[i] = normal_lpdf(y[i] | y_hat[i], sigma[1]);
    }
    if (obs_model == 2) {
      if (lower_truncation == 0) {
        log_lik[i] = neg_binomial_2_log_lpmf(y_int[i] | y_hat[i], nb2_phi[1]);
      } else {
        // Note that I had to remove T[lower_truncation, ] from the following line
        // and I think that will make this calculation incorrect
        // in the case of truncated negative binomial
        // the package will issue a warning
        log_lik[i] = neg_binomial_2_lpmf(y_int[i] | exp(y_hat[i]), nb2_phi[1]);
      }
    }
    if (obs_model == 4) {
      log_lik[i] = bernoulli_logit_lpmf(y_int[i] | y_hat[i]);
    }
    if (obs_model == 5) {
      log_lik[i] = poisson_log_lpmf(y_int[i] | y_hat[i]);
    }
    if (obs_model == 6) {
      log_lik[i] = lognormal_lpdf(y[i] | y_hat, sigma[1]);
    }
  }
}
