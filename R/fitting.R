#' Return a vector of parameters
#'
#' @param obs_error The observation error distribution
#' @param estimate_df Logical indicating whether the degrees of freedom
#'   parameter should be estimated
#' @param est_temporalRE Logical: estimate a random walk for the time variable?
#' @param estimate_ar Logical indicating whether the ar
#'   parameter should be estimated
#' @param fixed_intercept Should the intercept be fixed?
#' @param save_log_lik Logical: should the log likelihood for each data point be
#'   saved so that information criteria such as LOOIC or WAIC can be calculated?
#'   Defaults to \code{FALSE} so that the size of model objects is smaller.
stan_pars <- function(obs_error, estimate_df = TRUE, est_temporalRE = FALSE,
  estimate_ar = FALSE, fixed_intercept = FALSE, save_log_lik = FALSE) {
  p <- c("gp_sigma",
    "gp_scale",
    "B",
    switch(obs_error[[1]], normal = "sigma", gamma = "CV", nb2 = "nb2_phi"),
    "spatialEffectsKnots")
  if (estimate_df) p <- c("df", p)
  if (estimate_ar) p <- c("ar", p)
  if (est_temporalRE) {
    p <- c("year_sigma", "yearEffects", p)
    p <- p[!p=="B"] # no main effects if random walk for now
  }
  if (fixed_intercept) p <- p[p != "B"]
  if (save_log_lik) p <- c(p, "log_lik")
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
#' @param prior_rw_sigma The prior on the standard deviation parameter of the
#'   random walk process (if specified). Must be declared with
#'   \code{\link{half_t}}.
#' @param prior_intercept The prior on the intercept parameter. Must be declared
#'   with \code{\link{student_t}}.
#' @param prior_beta The prior on the slope parameters (if any). Must be
#'   declared with \code{\link{student_t}}.
#' @param fixed_df_value The fixed value for the student-t degrees of freedom
#'   parameter if the degrees of freedom parameter is fixed. If the degrees of
#'   freedom parameter is estimated then this argument is ignored. Must be 1 or
#'   greater.
#' @param estimate_df Logical: should the degrees of freedom perameter be
#'   estimated?
#' @param estimate_ar Logical: should the ar parameter be
#'   estimated?
#' @param fixed_ar_value The fixed value for temporal autoregressive parameter,
#'   between random fields at time(t) and time(t-1). If the ar parameter
#'   is estimated then this argument is ignored.
#' @param obs_error Character object indicating the observation process
#'   distribution (i.e. the GLM "family"). Links are hardcoded. Gamma, NB2,
#'   and Poisson have a log link. Binomial has a logit link.
#' @param covariance Character object describing the covariance
#'   function of the Gaussian Process.
#' @param algorithm Character object describing whether the model should be fit
#'   with full NUTS MCMC or via the variational inference mean-field approach.
#'   See \code{\link[rstan]{vb}}. Note that the variational inference approach
#'   should not be trusted for final inference and is much more likely to give
#'   incorrect inference than MCMC.
#' @param year_re Logical: estimate a random walk for the time variable? If
#'   \code{TRUE}, then no fixed effects (B coefficients) will be estimated.
#' @param nb_lower_truncation For NB2: lower truncation value. E.g. 0 for no
#'   truncation, 1 for 1 and all values above
#' @param control List to pass to \code{\link[rstan]{sampling}}
#' @param save_log_lik Logical: should the log likelihood for each data point be
#'   saved so that information criteria such as LOOIC or WAIC can be calculated?
#'   Defaults to \code{FALSE} so that the size of model objects is smaller.
#' @param ... Any other arguments to pass to \code{\link[rstan]{sampling}}.
#'
#' @export
#' @importFrom rstan sampling vb
#' @import Rcpp
#' @importFrom stats dist model.frame model.matrix model.response rnorm runif
#' @importFrom assertthat assert_that is.count is.number

rrfield <- function(formula, data, time, lon, lat, station = NULL, nknots = 25L,
  prior_gp_scale = half_t(3, 0, 5),
  prior_gp_sigma = half_t(3, 0, 5),
  prior_sigma = half_t(3, 0, 5),
  prior_rw_sigma = half_t(3, 0, 5),
  prior_intercept = student_t(3, 0, 10),
  prior_beta = student_t(3, 0, 2),
  fixed_df_value = 5,
  fixed_ar_value = 0,
  estimate_df = TRUE,
  estimate_ar = FALSE,
  obs_error = c("normal", "gamma", "poisson", "nb2", "binomial", "lognormal"),
  covariance = c("squared-exponential", "exponential"),
  algorithm = c("sampling", "meanfield"),
  year_re = FALSE,
  nb_lower_truncation = 0,
  control = list(adapt_delta = 0.9),
  save_log_lik = FALSE,
  ...) {

  # argument checks:
  is.count(nb_lower_truncation)
  assert_that(nb_lower_truncation >= 0)
  assert_that(fixed_df_value >= 1)
  is.number(fixed_ar_value)
  assert_that(all(unlist(lapply(list(obs_error, covariance, algorithm), is.character))))
  assert_that(
    obs_error[[1]] %in% c("normal", "gamma", "poisson", "nb2", "binomial", "lognormal"))
  assert_that(covariance[[1]] %in% c("squared-exponential", "exponential"))
  assert_that(algorithm[[1]] %in% c("sampling", "meanfield"))
  assert_that(is.logical(save_log_lik))
  assert_that(is.logical(estimate_df))
  assert_that(is.logical(estimate_ar))
  assert_that(is.logical(year_re))
  assert_that(is.list(control))

  if (nb_lower_truncation > 0)
    warning(paste0(c("Lower truncation with negative binomial has not been ",
      "extensively tested and calculation of log likelihood for information ",
      "criteria purposes has not been checked.")))

  mf <- model.frame(formula, data)
  X <- model.matrix(formula, mf)
  y <- model.response(mf, "numeric")
  fixed_intercept <- ifelse(ncol(X) == 0, TRUE, FALSE)

  # user inputs raw data. this function formats it for STAN
  data_list <- format_data(data = data, y = y, X = X, time = time,
    lon = lon, lat = lat, station=station, nknots = nknots, covariance = covariance,
    fixed_intercept = fixed_intercept)
  stan_data = data_list$spatglm_data
  data_knots = data_list$knots

  obs_model <- switch(obs_error[[1]], normal = 1L, gamma = 0L, nb2 = 2L, binomial = 4L,
    poisson = 5L, lognormal = 6L, stop(paste("observation model", obs_error[[1]], "is not defined.")))

  est_temporalRE <- ifelse(year_re, 1L, 0L)

  stan_data <- c(stan_data,
    list(prior_gp_scale = parse_t_prior(prior_gp_scale),
      prior_gp_sigma = parse_t_prior(prior_gp_sigma),
      prior_sigma = parse_t_prior(prior_sigma),
      prior_intercept = parse_t_prior(prior_intercept),
      prior_rw_sigma = parse_t_prior(prior_rw_sigma),
      prior_beta = parse_t_prior(prior_beta),
      sqexp_cov = switch(covariance[[1]], `squared-exponential` = 1L, exponential = 0L,
        stop(paste("covariance function", covariance[[1]], "is not defined."))),
      obs_model = obs_model,
      est_df = as.integer(estimate_df),
      est_ar = as.integer(estimate_ar),
      gamma_params = ifelse(obs_error[[1]] == "gamma", 1L, 0L),
      norm_params = ifelse(obs_error[[1]] %in% c("normal", "lognormal"), 1L, 0L),
      nb2_params = ifelse(obs_error[[1]] == "nb2", 1L, 0L),
      fixed_df_value = fixed_df_value,
      fixed_ar_value = fixed_ar_value,
      est_temporalRE = est_temporalRE,
      n_year_effects = ifelse(year_re, stan_data$nT, 0L),
      lower_truncation = nb_lower_truncation,
      fixed_intercept = as.integer(fixed_intercept)))

  if (obs_model %in% c(2L, 4L, 5L)) { # NB2 or binomial or poisson obs model
    stan_data <- c(stan_data, list(y_int = stan_data$y))
  } else {
    stan_data <- c(stan_data, list(y_int = rep(0L, stan_data$N)))
  }

  sampling_args <- list(
    object = stanmodels$rrfield,
    data = stan_data,
    pars = stan_pars(obs_error = obs_error, estimate_df = estimate_df,
      est_temporalRE = est_temporalRE, estimate_ar = estimate_ar,
      fixed_intercept = fixed_intercept, save_log_lik = save_log_lik),
    control = control, ...)

  if (algorithm[[1]] == "meanfield") {
    sampling_args$chains <- NULL
    sampling_args$control <- NULL
    m <- do.call(vb, sampling_args)
  } else {
    m <- do.call(sampling, sampling_args)
  }

  out <- list(model = m, knots = dplyr::as.tbl(as.data.frame(data_knots)), y = y, X = X,
    data = dplyr::as.tbl(data), formula = formula,
    covariance = covariance[[1]], lon = lon, lat = lat, time = time, year_re = year_re,
    station = data_list$stationID, obs_model = obs_model, fixed_intercept = fixed_intercept)
  out <- structure(out, class = "rrfield")
}
