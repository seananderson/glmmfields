#' Predict from a glmmfields model
#'
#' These functions extract posterior draws or credible intervals. The helper
#' functions are named to match those in the \pkg{rstanarm} package and call the
#' function `predict()` with appropriate argument values.
#'
#' @param object An object returned by [glmmfields()].
#' @param newdata Optionally, a data frame to predict on
#' @param interval Type of interval calculation. Same as for
#'   [stats::predict.lm()].
#' @param estimate_method Method for computing point estimate ("mean" or
#'   "median")
#' @param conf_level Probability level for the credible intervals.
#' @param type Whether the predictions are returned on "link" scale or
#'   "response" scale (Same as for [stats::predict.glm()]).
#' @param return_mcmc Logical. Should the full MCMC draws be returned for the
#'   predictions?
#' @param iter Number of MCMC iterations to draw. Defaults to all.
#' @param ... Ignored currently
#'
#' @importFrom stats median quantile rgamma rnbinom
#' @importFrom assertthat assert_that
#'
#' @examples
#' \donttest{
#' library(ggplot2)
#'
#' # simulate:
#' set.seed(1)
#' s <- sim_glmmfields(
#'   n_draws = 12, n_knots = 12, gp_theta = 2.5,
#'   gp_sigma = 0.2, sd_obs = 0.1
#' )
#'
#' # fit:
#' # options(mc.cores = parallel::detectCores()) # for parallel processing
#' m <- glmmfields(y ~ 0,
#'   data = s$dat, time = "time",
#'   lat = "lat", lon = "lon",
#'   nknots = 12, iter = 800, chains = 1
#' )
#'
#' # Predictions:
#' # Link scale credible intervals:
#' p <- predict(m, type = "link", interval = "confidence")
#' head(p)
#'
#' # Prediction intervals on new observations (include observation error):
#' p <- predictive_interval(m)
#' head(p)
#'
#' # Posterior prediction draws:
#' p <- posterior_predict(m, iter = 100)
#' dim(p) # rows are iterations and columns are data elements
#'
#' # Draws from the linear predictor (not in link space):
#' p <- posterior_linpred(m, iter = 100)
#' dim(p) # rows are iterations and columns are data elements
#'
#' # Use the `tidy` method to extract parameter estimates as a data frame:
#' head(tidy(m, conf.int = TRUE, conf.method = "HPDinterval"))
#'
#' # Make predictions on a fine-scale spatial grid:
#' pred_grid <- expand.grid(
#'   lat = seq(min(s$dat$lat), max(s$dat$lat), length.out = 25),
#'   lon = seq(min(s$dat$lon), max(s$dat$lon), length.out = 25),
#'   time = unique(s$dat$time)
#' )
#' pred_grid$prediction <- predict(m,
#'   newdata = pred_grid, type = "response", iter = 100,
#'   estimate_method = "median"
#' )$estimate
#'
#' ggplot(pred_grid, aes(lon, lat, fill = prediction)) +
#'   facet_wrap(~time) +
#'   geom_raster() +
#'   scale_fill_gradient2()
#' }
#' @name predict
NULL

#' @name predictive_interval
#' @rdname predict
#' @export
#' @importFrom rstantools predictive_interval
NULL

#' @name posterior_linpred
#' @rdname predict
#' @export
#' @importFrom rstantools posterior_linpred
NULL

#' @name posterior_predict
#' @rdname predict
#' @export
#' @importFrom rstantools posterior_predict
NULL

#' @rdname predict
#' @export
predictive_interval.glmmfields <- function(object, ...) {
  predict.glmmfields(object, interval = "prediction", type = "response",
    return_mcmc = FALSE, ...)
}

#' @rdname predict
#' @export
posterior_linpred.glmmfields <- function(object, ...) {
  predict.glmmfields(object, interval = "confidence", type = "response",
    return_mcmc = TRUE, ...)
}

#' @rdname predict
#' @export
posterior_predict.glmmfields <- function(object, ...) {
  predict.glmmfields(object, interval = "prediction", type = "response",
    return_mcmc = TRUE, ...)
}

#' @importFrom stats predict
#' @rdname predict
#' @export
predict.glmmfields <- function(object, newdata = NULL,
                               estimate_method = c("median", "mean"),
                               conf_level = 0.95,
                               interval = c("confidence", "prediction"),
                               type = c("link", "response"),
                               return_mcmc = FALSE,
                               iter = "all", ...) {
  estimate_method <- match.arg(estimate_method)
  interval <- match.arg(interval)
  type <- match.arg(type)

  assert_that(is.numeric(conf_level))
  assert_that(identical(length(conf_level), 1L))
  assert_that(conf_level > 0 && conf_level < 1)
  assert_that(identical(class(object), "glmmfields"))

  if (interval == "prediction" && type != "response") {
    stop("type must be 'response' if interval is 'prediction")
  }

  obs_model <- object$obs_model

  # newdata is df with time, y, lon, lat
  # if null, defaults to data used to fit model
  if (is.null(newdata)) newdata <- object$data

  response <- all.vars(nlme::getResponseFormula(object$formula))
  newdata[[response]] <- 1.0

  newdata <- tibble::as_tibble(newdata)
  # create model.matrix() as in fitting function, only with newdata
  X <- model.matrix(object$formula,
    model.frame(object$formula, newdata, na.action = na.omit))

  .missing <- colnames(object$X)[!colnames(object$X) %in% colnames(X)]
  if (length(.missing) > 0L)
    stop(paste(.missing, collapse = ", "), " are missing in `newdata`.")

  if (nrow(X) < nrow(newdata))
    stop("Some predictors in newdata had NA values.")

  time <- object$time
  knots <- object$knots

  dist_knots <- as.matrix(dist(knots))

  # Calculate distance from knots to grid
  dist_all <- as.matrix(stats::dist(rbind(
    newdata[, c(object$lon, object$lat)],
    knots
  )))
  n_locs <- nrow(newdata)

  # this is the transpose of the lower left corner
  dist_knots21 <- t(
    dist_all[-seq_len(n_locs), -seq(n_locs + 1, ncol(dist_all))]
  )

  # extract mcmc pars
  pars <- rstan::extract(object$model, permuted = TRUE)
  ##
  if (iter == "all") {
    mcmc.i <- seq_len(length(pars$lp__))
  } else {
    mcmc.i <- base::sample(seq_len(length(pars$lp__)), size = iter)
  }
  mcmc_draws <- length(mcmc.i)
  pred_values <- matrix(NA, n_locs, mcmc_draws)
  for (i in seq_len(mcmc_draws)) {
    # create cov matrix at knots
    if (object$covariance == "exponential") {
      covmat <- pars$gp_sigma[mcmc.i[i]] *
        exp(-dist_knots / pars$gp_theta[mcmc.i[i]])
      covmat21 <- pars$gp_sigma[mcmc.i[i]] *
        exp(-dist_knots21 / pars$gp_theta[mcmc.i[i]])
    }
    if (object$covariance == "squared-exponential") {
      covmat <- pars$gp_sigma[mcmc.i[i]] *
        exp(-(dist_knots^2) / (2 * pars$gp_theta[mcmc.i[i]]^2))
      covmat21 <- pars$gp_sigma[mcmc.i[i]] *
        exp(-(dist_knots21^2) / (2 * pars$gp_theta[mcmc.i[i]]^2))
    }
    if (object$covariance == "matern") {
      if (object$matern_kappa == 1.5) {
        transformed_dist <- sqrt(3) * dist_knots / pars$gp_theta[mcmc.i[i]]
        covmat <-
          pars$gp_sigma[mcmc.i[i]] *
            (1 + transformed_dist) * exp(-transformed_dist)

        transformed_dist <- sqrt(3) * dist_knots21 / pars$gp_theta[mcmc.i[i]]
        covmat21 <-
          pars$gp_sigma[mcmc.i[i]] *
            (1 + transformed_dist) * exp(-transformed_dist)
      }
      if (object$matern_kappa == 2.5) {
        transformed_dist <- sqrt(5) * dist_knots / pars$gp_theta[mcmc.i[i]]
        covmat <-
          pars$gp_sigma[mcmc.i[i]] *
            (1 + transformed_dist + (transformed_dist^2) / 3) * exp(-transformed_dist)

        transformed_dist <- sqrt(5) * dist_knots21 / pars$gp_theta[mcmc.i[i]]
        covmat21 <-
          pars$gp_sigma[mcmc.i[i]] *
            (1 + transformed_dist + (transformed_dist^2) / 3) * exp(-transformed_dist)
      }
    }
    # these are projected spatial effects, dim = new data points x time
    spat_eff_knots_i <- pars$spatialEffectsKnots[mcmc.i[i], , ]
    if (is.matrix(spat_eff_knots_i)) {
      spat_eff_knots_i <- t(spat_eff_knots_i)
    } else { # these are for one time slice and are a vector
      spat_eff_knots_i <- t(t(spat_eff_knots_i))
    }
    spat_effects <- covmat21 %*% solve(covmat) %*% spat_eff_knots_i

    rows <- seq_len(n_locs)

    if (identical(time, "null_time_")) {
      newdata$null_time_ <- 1
      time <- "null_time_"
    }
    cols <- as.numeric(as.factor(newdata[, time][[1]]))
    # check this for > 1 year. B will also have to be modified
    if (!object$year_re) {
      if (!object$fixed_intercept) {
        pred_values[, i] <- X %*% matrix(pars$B[mcmc.i[i], ], ncol = 1) +
          spat_effects[cbind(rows, cols)]
      } else {
        pred_values[, i] <- spat_effects[cbind(rows, cols)]
      }
    } else {
      pred_values[, i] <- spat_effects[cbind(rows, cols)] +
        pars$yearEffects[mcmc.i[i], ][cols]
    }
  }

  mcmc_draws <- ncol(pred_values)

  # if type == link, don't include observation/data model.

  # If predictions other than on link scale, use observation model and link to
  # generate (1) confidence intervals on mean or (2) prediction intervals
  # including obs error

  if (type == "response") {
    # gamma or NB2 or poisson:
    if (obs_model %in% c(0, 2, 5, 6)) pred_values <- exp(pred_values)

    if (obs_model == 1) {
      # normal, sigma is returned
      pp <- t(apply(pred_values, 1, function(x)
        stats::rnorm(mcmc_draws, mean = x, sd = pars$sigma[, 1])))
    }

    # binomial (plogis = inverse logit):
    if (obs_model == 4) pred_values <- stats::plogis(pred_values)

    if (obs_model == 0) {
      # gamma, CV is returned; gammaA = 1/(CV*CV)
      pp <- t(apply(pred_values, 1, function(x)
        stats::rgamma(mcmc_draws,
          shape = 1 / (pars$CV[, 1]^2),
          rate = 1 / (pars$CV[, 1]^2) / x
        )))
    }
    if (obs_model == 2) {
      # negative binomial, phi returned
      pp <- t(apply(pred_values, 1, function(x)
        stats::rnbinom(mcmc_draws, mu = x, size = pars$nb2_phi[, 1])))
    }
    if (obs_model == 4) {
      # binomial
      pp <- t(apply(pred_values, 1, function(x)
        stats::rbinom(mcmc_draws, size = 1, prob = x)))
    }
    if (obs_model == 5) {
      # poisson
      pp <- t(apply(pred_values, 1, function(x)
        stats::rpois(mcmc_draws, lambda = x)))
    }
    if (obs_model == 6) {
      # lognormal, sigma is returned
      pp <- t(apply(pred_values, 1, function(x)
        stats::rlnorm(mcmc_draws,
          meanlog = log(x),
          sdlog = pars$sigma[, 1]
        )))
    }
  }

  est_method <- switch(estimate_method[[1]], median = median, mean = mean)
  out <- data.frame(estimate = apply(pred_values, 1, est_method))

  if (interval == "confidence") {
    out$conf_low <- apply(pred_values, 1, quantile,
      probs = (1 - conf_level) / 2
    )
    out$conf_high <- apply(pred_values, 1, quantile,
      probs = 1 - (1 - conf_level) / 2
    )
  }
  if (interval == "prediction" && type == "response") {
    out$conf_low <- apply(pp, 1, quantile, probs = (1 - conf_level) / 2)
    out$conf_high <- apply(pp, 1, quantile, probs = 1 - (1 - conf_level) / 2)
  }

  out <- tibble::as_tibble(out)
  if (return_mcmc) out <- t(pred_values) # to match rstanarm generic methods
  out
}
