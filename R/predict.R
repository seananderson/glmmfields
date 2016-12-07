#' Predict from a rrfield model
#'
#' @param object An object returned by \code{\link{rrfield}}.
#' @param newdata Optionally, a data frame to predict on
#' @param mcmc_draws The number of MCMC samples to draw from the posterior
#' @param obs_model The observation error (0 = gamma, 1 = normal, 2 = neg bin)
#' @param ... Ignored currently
#'
#' @importFrom stats median quantile
#'
#' @export
predict.rrfield <- function(object, newdata = NULL, mcmc_draws, obs_model = 0, ...) {

  # newdata is df with time, y, lon, lat
  # if null, defaults to data used to fit model
  if(is.null(newdata)) newdata = object$data
  # create model.matrix() as in fitting function, only with newdata
  X = model.matrix(object$formula, model.frame(object$formula, newdata))

  time = object$time
  knots = object$knots
  n_knots = nrow(knots)

  distKnots = as.matrix(dist(knots))

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
    if(object$covariance=="exponential") {
      covmat = pars$gp_sigma[mcmc.i[i]] *
        exp(-distKnots/pars$gp_scale[mcmc.i[i]])
      covmat21 = pars$gp_sigma[mcmc.i[i]] *
        exp(-dist_knots21/pars$gp_scale[mcmc.i[i]])
    } else {
      covmat = pars$gp_sigma[mcmc.i[i]] *
        exp(-2*(distKnots^2)/(pars$gp_scale[mcmc.i[i]]^2))
      covmat21 = pars$gp_sigma[mcmc.i[i]] *
        exp(-2*(dist_knots21^2)/(pars$gp_scale[mcmc.i[i]]^2))
    }

    # these are projected spatial effects, dim = new data points x time
    spat_effects = covmat21 %*% solve(covmat) %*% t(pars$spatialEffectsKnots[mcmc.i[i],,])

    rows = seq_len(n_locs)
    cols = as.numeric(as.factor(newdata[,time]))
    # check this for > 1 year. B will also have to be modified
    if(object$year_re == FALSE) {
      pred_values[,i] = X %*% matrix(pars$B[mcmc.i[i],], nrow=1) +
        spat_effects[cbind(rows,cols)]
    } else {
      pred_values[,i] = spat_effects[cbind(rows,cols)] + pars$yearEffects[mcmc.i[i],][cols]
    }
  }

  if(obs_model %in% c(0,2)) {
    pred_values = exp(pred_values)
  }

  if(obs_model == 0) {
    # gamma, CV is returned. gammaA = 1/(CV*CV)
    pred_values_obs = pred_values
    for(i in 1:mcmc_draws) {
    pred_values_obs[,i] = rgamma(nrow(pred_values), shape = 1/(pars$CV[mcmc.i[i]]^2),
      rate = 1/(pars$CV[mcmc.i[i]]^2)/pred_values[,i])
    }
  }
  if(obs_model == 1) {
    # normal
  }
  if(obs_model == 2) {
    # negative binomial
  }

  summary_mat = data.frame("mean" = apply(pred_values, 1, mean),
    "median" = apply(pred_values, 1, median),
    "predictive_lower2.5" = apply(pred_values, 1, quantile, 0.025),
    "predictive_upper97.5" = apply(pred_values, 1, quantile, 0.975),
    "confint_lower2.5" = apply(pred_values_obs, 1, quantile, 0.025),
    "confint_upper97.5" = apply(pred_values_obs, 1, quantile, 0.975))

  return(list(predictions = pred_values, summary = summary_mat))

}
