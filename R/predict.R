#' Predict from an rrfield model
#'
#' @param object An object returned by \code{\link{rrfield}}.
#' @param newdata Optionally, a data frame to predict on
#' @param interval Type of interval calculation. Same as \code{\link[stats]{predict.lm}}
#' @param estimate_method Method for computing point estimate ("mean" or median")
#' @param conf_level Probability level for CI
#' @param type Whether the predictions are returned on "link" scale or "response" scale (Same as \code{\link[stats]{predict.glm}})
#' @param ... Ignored currently
#'
#' @importFrom stats median quantile predict rgamma rnbinom
#' @importFrom assertthat assert_that
#'
#' @export
predict.rrfield <- function(object, newdata = NULL,
  estimate_method = c("median", "mean"), conf_level = 0.95,
  interval = c("confidence", "prediction"), type = c("link", "response"), ...) {

  assert_that(is.character(estimate_method[[1]]))
  assert_that(is.character(interval[[1]]))
  assert_that(is.character(type[[1]]))
  assert_that(is.numeric(conf_level))
  assert_that(identical(length(conf_level), 1L))
  assert_that(conf_level > 0 & conf_level < 1)
  assert_that(identical(class(object), "rrfield"))
  assert_that(type[[1]] %in% c("link", "response"))
  assert_that(estimate_method[[1]] %in% c("median", "mean"))
  assert_that(interval[[1]] %in% c("confidence", "prediction"))

  obs_model <- object$obs_model

  # newdata is df with time, y, lon, lat
  # if null, defaults to data used to fit model
  if(is.null(newdata)) newdata <- object$data
  # create model.matrix() as in fitting function, only with newdata
  X <- model.matrix(object$formula, model.frame(object$formula, newdata))

  time <- object$time
  knots <- object$knots
  n_knots <- nrow(knots)

  # pred_values <- t(rstan::extract(object$model, pars = "y_hat")[[1]])

  distKnots <- as.matrix(dist(knots))

  # Calculate distance from knots to grid
  dist_all <- as.matrix(stats::dist(rbind(newdata[, c(object$lon, object$lat)], knots)))
  n_locs <- nrow(newdata)

  # this is the transpose of the lower left corner
  dist_knots21 <- t(
    dist_all[-c(seq_len(n_locs)), -c((n_locs + 1):ncol(dist_all))])

  # extract mcmc pars
  pars <- rstan::extract(object$model, permuted = TRUE)
  ##
  mcmc.i <- seq_len(length(pars$lp__))
  mcmc_draws <- max(mcmc.i)
  pred_values <- matrix(NA, n_locs, mcmc_draws)
  for(i in 1:mcmc_draws) {
    # create cov matrix @ knots
    if(object$covariance == "exponential") {
      covmat <- pars$gp_sigma[mcmc.i[i]] *
        exp(-distKnots/pars$gp_scale[mcmc.i[i]])
      covmat21 <- pars$gp_sigma[mcmc.i[i]] *
        exp(-dist_knots21/pars$gp_scale[mcmc.i[i]])
    } else {
      covmat <- pars$gp_sigma[mcmc.i[i]] *
        exp(-(distKnots^2)/(2 * pars$gp_scale[mcmc.i[i]]^2))
      covmat21 <- pars$gp_sigma[mcmc.i[i]] *
        exp(-(dist_knots21^2)/(2 * pars$gp_scale[mcmc.i[i]]^2))
    }

    # these are projected spatial effects, dim = new data points x time
    spat_effects <- covmat21 %*% solve(covmat) %*% t(pars$spatialEffectsKnots[mcmc.i[i],,])

    rows <- seq_len(n_locs)
    cols <- as.numeric(as.factor(newdata[, time][[1]]))
    # check this for > 1 year. B will also have to be modified
    if(object$year_re == FALSE) {
      if (!object$fixed_intercept) {
        pred_values[,i] <- X %*% matrix(pars$B[mcmc.i[i],], ncol = 1) +
          spat_effects[cbind(rows,cols)]
      } else {
        pred_values[,i] <- spat_effects[cbind(rows,cols)]
      }
    } else {
      pred_values[,i] <- spat_effects[cbind(rows,cols)] + pars$yearEffects[mcmc.i[i],][cols]
    }
  }

  mcmc_draws <- ncol(pred_values)

  # if type == link, don't include observation/data model.

  # If predictions other than on link scale, use observation model and link to
  # generate (1) confidence intervals on mean or (2) prediction intervals including obs error

  if(type[[1]] == "response") {
    if(obs_model %in% c(0, 2, 5, 6)) pred_values <- exp(pred_values) # gamma or NB2 or poisson

    if(obs_model == 1) {
      # normal, sigma is returned
      pp <- t(apply(pred_values, 1, function(x) stats::rnorm(mcmc_draws, mean = x, sd = pars$sigma[,1])))
    }

    if(obs_model == 4) pred_values <- stats::plogis(pred_values) # binomial (plogis = inverse logit)

    if(obs_model == 0) {
      # gamma, CV is returned; gammaA = 1/(CV*CV)
      pp <- t(apply(pred_values, 1, function(x) stats::rgamma(mcmc_draws, shape = 1/(pars$CV[,1]^2),
        rate = 1/(pars$CV[,1]^2)/x)))
    }
    if(obs_model == 2) {
      # negative binomial, phi returned
      pp <- t(apply(pred_values, 1, function(x) stats::rnbinom(mcmc_draws, mu = x, size = pars$nb2_phi[,1])))
    }
    if(obs_model == 4) {
      # binomial
      pp <- t(apply(pred_values, 1, function(x) stats::rbinom(mcmc_draws, size = 1, prob = x)))
    }
    if(obs_model == 5) {
      # poisson
      pp <- t(apply(pred_values, 1, function(x) stats::rpois(mcmc_draws, lambda = x)))
    }
    if(obs_model == 6) {
      # lognormal, sigma is returned
      pp <- t(apply(pred_values, 1, function(x) stats::rlnorm(mcmc_draws, meanlog = x, sdlog = pars$sigma[,1])))
    }
  }

  est_method <- switch(estimate_method[[1]], median = median, mean = mean)
  out <- data.frame(estimate = apply(pred_values, 1, est_method))

  if (interval[[1]] == "confidence") {
    out$conf_low <- apply(pred_values, 1, quantile, probs = (1 - conf_level) / 2)
    out$conf_high <- apply(pred_values, 1, quantile, probs = 1 - (1 - conf_level) / 2)
  }
  if (interval[[1]] == "prediction" & (type[[1]] == "response" | (type[[1]] != "response" & obs_model==1))) {
    out$conf_low <- apply(pp, 1, quantile, probs = (1 - conf_level) / 2)
    out$conf_high <- apply(pp, 1, quantile, probs = 1 - (1 - conf_level) / 2)
  }

  dplyr::as.tbl(out)
}
