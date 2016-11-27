#' Return a vector of parameters
#'
#' @param obs_error The observation error distribution
#' @param estimate_df Logical indicating whether the degrees of freedom
#'   parameter should be estimated
#' @param est_temporalRE Logical: estimate a random walk for the time variable?
stan_pars <- function(obs_error, estimate_df = TRUE, est_temporalRE = FALSE) {
  p <- c("gp_sigma",
    "gp_scale",
    "B",
    switch(obs_error[[1]], normal = "sigma", gamma = "CV", nb2 = "nb2_phi"),
    "spatialEffectsKnots")
  if (estimate_df) p <- c("df", p)
  if (est_temporalRE) {
    p <- c("year_sigma", "yearEffects", p)
    p <- p[!p=="B"] # no main effects if random walk for now
  }
  p
}

#' Fit a robust spatiotemporal random fields model
#'
#' @param formula The model formula
#' @param data A data frame
#' @param time A character object giving the name of the time column
#' @param lon A character object giving the name of the longitude column
#' @param lat A character object giving the name of the latitude column
#' @param station A numeric vector giving the integer ID of the station
#' @param nknots The number of knots to use in the predictive process model
#' @param prior_gp_scale The prior on the Gaussian Process scale parameter. Must
#'   be declared with \code{\link{half_t}}.
#' @param prior_gp_sigma The prior on the Gaussian Process sigma parameter. Must
#'   be declared with \code{\link{half_t}}.
#' @param prior_sigma The prior on the observation process scale parameter. Must
#'   be declared with \code{\link{half_t}}. This acts as a
#'   substitute for the scale parameter in whatever observation distribution is
#'   being used. I.e. the CV for the Gamma or the dispersion parameter for
#'   the negative binomial.
#' @param prior_intercept The prior on the intercept parameter. Must be declared
#'   with \code{\link{student_t}}.
#' @param prior_beta The prior on the slope parameters (if any). Must be
#'   declared with \code{\link{student_t}}.
#' @param fixed_df_value The fixed value for the student-t degrees of freedom
#'   parameter if the degrees of freedom parameter is fixed. If the degrees of
#'   freedom parameter is estimated then this argument is ignored.
#' @param estimate_df Logical: should the degrees of freedom perimeter be
#'   estimated?
#' @param obs_error Character object indicating the observation process
#'   distribution.
#' @param covariance Character object describing the covariance
#'   function of the Gaussian Process.
#' @param algorithm Character object describing whether the model should be fit
#'   with full NUTS MCMC or via the variational inference mean-field approach.
#'   See \code{\link[rstan]{vb}}. Note that the variational inference approach
#'   should not be trusted for final inference and is much more likely to give
#'   incorrect inference than MCMC.
#' @param year_re Logical: estimate a random walk for the time variable? If
#'   \code{TRUE}, then no fixed effects (B coefficients) will be estimated.
#' @param ... Any other arguments to pass to \code{\link[rstan]{sampling}}.
#'
#' @export
#' @importFrom rstan sampling vb
#' @import Rcpp
#' @importFrom stats dist model.frame model.matrix model.response rnorm runif

rrfield <- function(formula, data, time, lon, lat, station = "", nknots = 25L,
  prior_gp_scale = student_t(3, 0, 5),
  prior_gp_sigma = student_t(3, 0, 5),
  prior_sigma = student_t(3, 0, 5),
  prior_intercept = student_t(3, 0, 10),
  prior_beta = student_t(3, 0, 2),
  fixed_df_value = 5,
  estimate_df = TRUE,
  obs_error = c("normal", "gamma", "nb2"),
  covariance = c("squared-exponential", "exponential"),
  algorithm = c("sampling", "meanfield"),
  year_re = FALSE,
  ...) {

  mf <- model.frame(formula, data)
  X <- model.matrix(formula, mf)
  y <- model.response(mf, "numeric")

  # user inputs raw data. this function formats it for STAN
  data_list <- format_data(data = data, y = y, X = X, time = time,
    lon = lon, lat = lat, station=station, nknots = nknots, covariance = covariance)
  stan_data = data_list$spatglm_data
  data_knots = data_list$knots

  obs_model <- switch(obs_error[[1]], normal = 1L, gamma = 0L, nb2 = 2L, 1L)

  est_temporalRE <- ifelse(year_re, 1L, 0L)

  stan_data <- c(stan_data,
    list(prior_gp_scale = parse_t_prior(prior_gp_scale),
      prior_gp_sigma = parse_t_prior(prior_gp_sigma),
      prior_sigma = parse_t_prior(prior_sigma),
      prior_intercept = parse_t_prior(prior_intercept),
      prior_beta = parse_t_prior(prior_beta),
      sqexp_cov = switch(covariance[[1]], `squared-exponential` = 1L,
        exponential = 0L, 1L),
      obs_model = obs_model,
      est_df = as.integer(estimate_df),
      gamma_params = ifelse(obs_error[[1]] == "gamma", 1L, 0L),
      norm_params = ifelse(obs_error[[1]] == "normal", 1L, 0L),
      nb2_params = ifelse(obs_error[[1]] == "nb2", 1L, 0L),
      fixed_df_value = fixed_df_value,
      est_temporalRE = est_temporalRE,
      n_year_effects = ifelse(year_re, stan_data$nT, 0L)))

  if (obs_model == 2) { # NB2 obs model
    stan_data <- c(stan_data, list(y_int = stan_data$y))
  } else {
    stan_data <- c(stan_data, list(y_int = rep(0L, stan_data$N)))
  }

  sampling_args <- list(
    object = stanmodels$rrfield,
    data = stan_data,
    pars = stan_pars(obs_error = obs_error, estimate_df = estimate_df, est_temporalRE = est_temporalRE),
    ...)

  if (algorithm[[1]] == "meanfield") {
    sampling_args$chains <- NULL
    m <- do.call(vb, sampling_args)
  } else {
    m <- do.call(sampling, sampling_args)
  }

  out <- list(model = m, knots = data_knots, y = y, X = X, data = data, formula = formula,
    covariance = covariance[[1]], lon = lon, lat = lat, time = time, year_re = year_re,
    station=data_list$stationID)
  out <- structure(out, class = "rrfield")
}
