#' Fit a spatiotemporal random fields GLMM
#'
#' Fit a spatiotemporal random fields model that optionally uses the MVT
#' distribution instead of a MVN distribution to allow for spatial extremes
#' through time. It is also possible to fit a spatial random fields model
#' without a time component.
#'
#' @param formula The model formula
#' @param data A data frame
#' @param time A character object giving the name of the time column
#' @param lon A character object giving the name of the longitude column
#' @param lat A character object giving the name of the latitude column
#' @param nknots The number of knots to use in the predictive process model.
#'   Smaller values will be faster but may not adequately represent the shape
#'   of the spatial pattern.
#' @param prior_gp_scale The prior on the Gaussian Process scale parameter. Must
#'   be declared with \code{\link{half_t}}. Here, and throughout, priors that
#'   are normal or half-normal can be implemented by setting the first
#'   parameter in the half-t or student-t distribution to a large value.
#' @param prior_gp_sigma The prior on the Gaussian Process sigma parameter. Must
#'   be declared with \code{\link{half_t}}.
#' @param prior_sigma The prior on the observation process scale parameter. Must
#'   be declared with \code{\link{half_t}}. This acts as a substitute for the
#'   scale parameter in whatever observation distribution is being used. I.e.
#'   the CV for the Gamma or the dispersion parameter for the negative
#'   binomial.
#' @param prior_rw_sigma The prior on the standard deviation parameter of the
#'   random walk process (if specified). Must be declared with
#'   \code{\link{half_t}}.
#' @param prior_intercept The prior on the intercept parameter. Must be declared
#'   with \code{\link{student_t}}.
#' @param prior_beta The prior on the slope parameters (if any). Must be
#'   declared with \code{\link{student_t}}.
#' @param prior_ar The prior on the AR parameter. Must be
#'   declared with \code{\link{student_t}}.
#' @param fixed_df_value The fixed value for the student-t degrees of freedom
#'   parameter if the degrees of freedom parameter is fixed in the MVT. If the
#'   degrees of freedom parameter is estimated then this argument is ignored.
#'   Must be 1 or greater. Very large values (e.g. the default value)
#'   approximate the normal distribution. If the value is >=1000 then a true
#'   MVN distribution will be fit.
#' @param estimate_df Logical: should the degrees of freedom perameter be
#'   estimated?
#' @param estimate_ar Logical: should the AR (autoregressive) parameter be
#'   estimated? Here, this refers to a autoregressive process in the evolution
#'   of the spatial field through time.
#' @param fixed_ar_value The fixed value for temporal autoregressive parameter,
#'   between random fields at time(t) and time(t-1). If the AR parameter
#'   is estimated then this argument is ignored.
#' @param family Family object describing the observation model. Note that only
#'   one link is implemented for each distribution. Gamma, negative binomial
#'   (specified as \code{nbinom2(link = "log")}, and Poisson must have a log
#'   link. Binomial must have a logit link. Also implemented is
#'   \code{lognormal(link = "log")}. Besides the negative binomial and
#'   lognormal, other families are specified as shown in
#'   \code{\link[stats]{family}}.
#' @param covariance The covariance function of the Gaussian Process.
#'   One of "squared-exponential", "exponential", or "matern".
#' @param matern_kappa Optional parameter for the Matern covariance function.
#'   Optional values are 1.5 or 2.5. Values of 0.5 are equivalent to exponential.
#' @param algorithm Character object describing whether the model should be fit
#'   with full NUTS MCMC or via the variational inference mean-field approach.
#'   See \code{\link[rstan]{vb}}. Note that the variational inference approach
#'   should not be trusted for final inference and is much more likely to give
#'   incorrect inference than MCMC.
#' @param year_re Logical: estimate a random walk for the time variable? If
#'   \code{TRUE}, then no fixed effects (B coefficients) will be estimated.
#'   \code{TRUE}, then \code{prior_intercept} will be used as the prior for
#'   the initial value in time.
#' @param nb_lower_truncation For NB2 only: lower truncation value. E.g. 0 for
#'   no truncation, 1 for 1 and all values above. Note that estimation is
#'   likely to be considerably slower with lower truncation because the
#'   sampling is not vectorized. Also note that the log likelihood values
#'   returned for estimating quantities like LOOIC will not be correct if
#'   lower truncation is implemented.
#' @param control List to pass to \code{\link[rstan]{sampling}}
#' @param save_log_lik Logical: should the log likelihood for each data point be
#'   saved so that information criteria such as LOOIC or WAIC can be calculated?
#'   Defaults to \code{FALSE} so that the size of model objects is smaller.
#' @param gp_sigma_scaling_factor An optional factor to scale the spatial
#'   variance parameter gp_sigma by to help sampling. Sometimes, especially in
#'   the case of an auto regressive model, gp_sigma could be quite small and
#'   Stan can struggle with tiny parameters. If so, try setting the scaling
#'   factor to a value less than 1.
#' @param ... Any other arguments to pass to \code{\link[rstan]{sampling}}.
#'
#' @details
#' Note that there is no guarantee the priors will remain the same in future
#' versions. Therefore it is important that you specify any priors that are
#' used in your model, even if they replicate the defaults in the package.
#'
#' @export
#' @importFrom rstan sampling vb
#' @import Rcpp
#' @importFrom stats dist model.frame model.matrix model.response rnorm runif
#' @importFrom assertthat assert_that is.count is.number
#' @importFrom stats gaussian

glmmfields <- function(formula, data, lon, lat,
  time = NULL,
  nknots = 15L,
  prior_gp_scale = half_t(3, 0, 5),
  prior_gp_sigma = half_t(3, 0, 5),
  prior_sigma = half_t(3, 0, 5),
  prior_rw_sigma = half_t(3, 0, 5),
  prior_intercept = student_t(3, 0, 10),
  prior_beta = student_t(3, 0, 3),
  prior_ar = student_t(1000, 0, 0.5),
  fixed_df_value = 1000,
  fixed_ar_value = 0,
  estimate_df = FALSE,
  estimate_ar = FALSE,
  family = gaussian(link = "identity"),
  covariance = c("squared-exponential", "exponential", "matern"),
  matern_kappa = 0.5,
  algorithm = c("sampling", "meanfield"),
  year_re = FALSE,
  nb_lower_truncation = 0,
  control = list(adapt_delta = 0.9),
  save_log_lik = FALSE,
  gp_sigma_scaling_factor = 1,
  ...) {

  # argument checks:
  is.count(nb_lower_truncation)
  assert_that(nb_lower_truncation >= 0)
  assert_that(fixed_df_value >= 1)
  is.number(fixed_ar_value)
  assert_that(all(unlist(lapply(list(covariance, algorithm), is.character))))
  assert_that(covariance[[1]] %in% c("squared-exponential", "exponential", "matern"))
  assert_that(algorithm[[1]] %in% c("sampling", "meanfield"))
  assert_that(matern_kappa %in% c(0.5, 1.5, 2.5))
  assert_that(is.logical(save_log_lik))
  assert_that(is.logical(estimate_df))
  assert_that(is.logical(estimate_ar))
  assert_that(is.logical(year_re))
  assert_that(is.list(control))

  family <- check_family(family)
  obs_error <- tolower(family$family)
  assert_that(obs_error[[1]] %in% c("gaussian", "gamma", "poisson", "nbinom2",
    "binomial", "lognormal"))

  if(covariance[[1]] == "matern" & matern_kappa %in% c(1.5, 2.5) == FALSE) {
    warning(paste0(c("Matern covariance specified, but Matern kappa not 1.5 or 2.5",
      ": defaulting to 0.5, or exponential")))
    covariance[[1]] = "exponential"
  }

  if (nb_lower_truncation > 0)
    warning(paste0(c("Lower truncation with negative binomial has not been ",
      "extensively tested and calculation of log likelihood for information ",
      "criteria purposes is likely to be incorrect.")))

  mf <- model.frame(formula, data)
  X <- model.matrix(formula, mf)
  y <- model.response(mf, "numeric")
  fixed_intercept <- ifelse(ncol(X) == 0, TRUE, FALSE)

  if (is.null(time)) {
    data$time <- 1
    time <- "time"
  }

  if ("station" %in% names(list(...)))
    stop(paste("The 'station' argument is no longer needed when calling glmmfields().",
      "Please remove it."))
  data$station_ <- paste(data[[lon]], data[[lat]])

  # user inputs raw data. this function formats it for STAN
  data_list <- format_data(data = data, y = y, X = X, time = time,
    lon = lon, lat = lat, station = "station_", nknots = nknots,
    covariance = covariance,
    fixed_intercept = fixed_intercept)
  stan_data = data_list$spatglm_data
  data_knots = data_list$knots

  obs_model <- switch(obs_error[[1]], gaussian = 1L, gamma = 0L, nbinom2 = 2L,
    binomial = 4L, poisson = 5L, lognormal = 6L,
    stop(paste("observation model", obs_error[[1]], "is not defined.")))

  est_temporalRE <- ifelse(year_re, 1L, 0L)

  stan_data <- c(stan_data,
    list(prior_gp_scale = parse_t_prior(prior_gp_scale),
      prior_gp_sigma = parse_t_prior(prior_gp_sigma),
      prior_sigma = parse_t_prior(prior_sigma),
      prior_intercept = parse_t_prior(prior_intercept),
      prior_rw_sigma = parse_t_prior(prior_rw_sigma),
      prior_beta = parse_t_prior(prior_beta),
      prior_ar = parse_t_prior(prior_ar),
      cov_func = switch(covariance[[1]], exponential = 0L, `squared-exponential` = 1L,
        matern = 2L,
        stop(paste("covariance function", covariance[[1]], "is not defined."))),
      obs_model = obs_model,
      est_df = as.integer(estimate_df),
      est_ar = as.integer(estimate_ar),
      gamma_params = ifelse(obs_error[[1]] == "gamma", 1L, 0L),
      norm_params = ifelse(obs_error[[1]] %in% c("gaussian", "lognormal"), 1L, 0L),
      nb2_params = ifelse(obs_error[[1]] == "nbinom2", 1L, 0L),
      fixed_df_value = fixed_df_value,
      fixed_ar_value = fixed_ar_value,
      est_temporalRE = est_temporalRE,
      n_year_effects = ifelse(year_re, stan_data$nT, 0L),
      lower_truncation = nb_lower_truncation,
      fixed_intercept = as.integer(fixed_intercept),
      matern_kappa = matern_kappa,
      gp_sigma_scaling_factor = gp_sigma_scaling_factor,
      nW = ifelse(fixed_df_value[[1]] > 999 & !estimate_df, 0L, stan_data$nT)))

  if (obs_model %in% c(2L, 4L, 5L)) { # integers: NB2 or binomial or poisson obs model
    stan_data <- c(stan_data, list(y_int = stan_data$y))
  } else {
    stan_data <- c(stan_data, list(y_int = rep(0L, stan_data$N)))
  }

  sampling_args <- list(
    object = stanmodels$glmmfields,
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

  out <- list(model = m,
    knots = dplyr::as.tbl(as.data.frame(data_knots)),
    y = y, X = X,
    data = dplyr::as.tbl(data), formula = formula,
    covariance = covariance[[1]], matern_kappa = matern_kappa,
    lon = lon, lat = lat,
    time = time, year_re = year_re,
    station = data_list$stationID, obs_model = obs_model,
    fixed_intercept = fixed_intercept, family = family)
  out <- structure(out, class = "glmmfields")
}
