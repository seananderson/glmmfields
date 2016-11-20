#' @export

predict_rrfield <- function(fitted_model, new_data, mcmc_draws, time="time") {

  # newdata is df with time, y, lon, lat

  knots = fitted_model$knots
  n_knots = nrow(knots)

  distKnots = as.matrix(dist(knots))
  distKnotsSq = distKnots^2 # squared distances

  # Calculate distance from knots to grid
  dist_all = as.matrix(stats::dist(rbind(new_data[, c(fitted_model$lon, fitted_model$lat)], knots)))
  n_locs = nrow(new_data)

  # this is the transpose of the lower left corner
  dist_knots21 <- t(
    dist_all[-c(seq_len(n_locs)), -c((n_locs + 1):ncol(dist_all))])


  # extract mcmc pars
  pars = rstan::extract(fitted_model$model, permuted = TRUE)

  mcmc.i = sample(1:length(pars$lp__), size=mcmc_draws, replace=F)
  pred_values = matrix(NA, n_locs, mcmc_draws)
  for(i in 1:mcmc_draws) {
    # create cov matrix @ knots
    if(fitted_model$correlation=="exponential") {
      covmat = pars$gp_sigma[mcmc.i[i]] * exp(-distKnots/pars$gp_scale[mcmc.i[i]])
      covmat21 = pars$gp_sigma[mcmc.i[i]] * exp(-dist_knots21/pars$gp_scale[mcmc.i[i]])
    } else {
      covmat = pars$gp_sigma[mcmc.i[i]] * exp(-2*(distKnots^2)/(pars$gp_scale[mcmc.i[i]]^2))
      covmat21 = pars$gp_sigma[mcmc.i[i]] * exp(-2*(distKnots21^2)/(pars$gp_scale[mcmc.i[i]]^2))
    }

    # these are projected spatial effects, dim = new data points x time
    spat_effects = covmat21 %*% solve(covmat) %*% t(pars$spatialEffectsKnots[mcmc.i[i],,])

    rows = seq_len(n_locs)
    cols = new_data$time
    pred_values[,i] = pars$B[mcmc.i[i],1] + spat_effects[cbind(rows,cols)] # check this for > 1 year. B will also have to be modified
  }

  pred_values
}
