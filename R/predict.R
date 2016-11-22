#' Predict from a rrfield model
#'
#' @param object An object returned by \code{\link{rrfield}}.
#' @param newdata A data frame to predict on
#' @param mcmc_draws The number of MCMC samples to draw from the posterior
#' @param ... Ignored currently
#'
#' @export
predict.rrfield <- function(object, newdata, mcmc_draws, ...) {

  # newdata is df with time, y, lon, lat

  time = object$time
  knots = object$knots
  n_knots = nrow(knots)

  distKnots = as.matrix(dist(knots))
  distKnotsSq = distKnots^2 # squared distances

  # Calculate distance from knots to grid
  dist_all = as.matrix(stats::dist(rbind(newdata[, c(object$lon, object$lat)], knots)))
  n_locs = nrow(newdata)

  # this is the transpose of the lower left corner
  dist_knots21 <- t(
    dist_all[-c(seq_len(n_locs)), -c((n_locs + 1):ncol(dist_all))])


  # extract mcmc pars
  pars = rstan::extract(object$model, permuted = TRUE)

  mcmc.i = sample(1:length(pars$lp__), size=mcmc_draws, replace=F)
  pred_values = matrix(NA, n_locs, mcmc_draws)
  for(i in 1:mcmc_draws) {
    # create cov matrix @ knots
    if(object$correlation=="exponential") {
      covmat = pars$gp_sigma[mcmc.i[i]] * exp(-distKnots/pars$gp_scale[mcmc.i[i]])
      covmat21 = pars$gp_sigma[mcmc.i[i]] * exp(-dist_knots21/pars$gp_scale[mcmc.i[i]])
    } else {
      covmat = pars$gp_sigma[mcmc.i[i]] * exp(-2*(distKnots^2)/(pars$gp_scale[mcmc.i[i]]^2))
      covmat21 = pars$gp_sigma[mcmc.i[i]] * exp(-2*(dist_knots21^2)/(pars$gp_scale[mcmc.i[i]]^2))
    }

    # these are projected spatial effects, dim = new data points x time
    spat_effects = covmat21 %*% solve(covmat) %*% t(pars$spatialEffectsKnots[mcmc.i[i],,])

    rows = seq_len(n_locs)
    cols = newdata[,time]
    pred_values[,i] = pars$B[mcmc.i[i],1] + spat_effects[cbind(rows,cols)] # check this for > 1 year. B will also have to be modified
  }

  pred_values
}
