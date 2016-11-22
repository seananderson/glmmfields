stan_pars <- function(obs_error, estimate_df = TRUE) {
  p <- c("gp_sigma",
    "gp_scale",
    "B",
    "spatialEffectsKnots",
    switch(obs_error[[1]], normal = "sigma", gamma = "CV", nb2 = "nb2_phi")
  )
  if (estimate_df) p <- c("df", p)
  p
}

parse_t_prior <- function(x) {
  as.vector(unlist(x)[-1], mode = "numeric")
}

#' Fit a robust spatiotemporal random fields model
#'
#' @export
#' @importFrom rstanarm student_t normal
#' @importFrom rstan sampling vb
#' @import Rcpp
#' @importFrom stats dist model.frame model.matrix model.response rnorm runif
rrfield <- function(formula, data, time, lon, lat, nknots = 25L,
  prior_gp_scale = student_t(3, 0, 5),
  prior_gp_sigma = student_t(3, 0, 5),
  prior_sigma = student_t(3, 0, 5),
  fixed_df_value = 5,
  estimate_df = TRUE,
  obs_error = c("normal", "gamma", "nb2"),
  correlation = c("gaussian", "exponential"),
  algorithm = c("sampling", "meanfield"),
  ...) {

  mf <- model.frame(formula, data)
  X <- model.matrix(formula, mf)
  y <- model.response(mf, "numeric")

  # user inputs raw data. this function formats it for STAN
  data_list <- format_data(data = data, y = y, X = X, time = time,
    lon = lon, lat = lat, nknots = nknots, correlation = correlation)
  stan_data = data_list$spatglm_data
  data_knots = data_list$knots

  obs_model <- switch(obs_error[[1]], normal = 1L, gamma = 0L, nb2 = 2L, 1L)

  stan_data <- c(stan_data,
    list(prior_gp_scale = parse_t_prior(prior_gp_scale),
      prior_gp_sigma = parse_t_prior(prior_gp_sigma),
      prior_sigma = parse_t_prior(prior_sigma),
      gauss_cor = switch(correlation[[1]], gaussian = 1L, exponential = 0L, 1L),
      obs_model = obs_model,
      est_df = as.integer(estimate_df),
      gamma_params = ifelse(obs_error[[1]] == "gamma", 1L, 0L),
      norm_params = ifelse(obs_error[[1]] == "normal", 1L, 0L),
      nb2_params = ifelse(obs_error[[1]] == "nb2", 1L, 0L),
      fixed_df_value = fixed_df_value))

  if (obs_model == 2) { # NB2 obs model
    stan_data <- c(stan_data, list(y_int = stan_data$y))
  } else {
    stan_data <- c(stan_data, list(y_int = rep(0L, stan_data$N)))
  }

  sampling_args <- list(
    object = stanmodels$rrfield,
    data = stan_data,
    pars = stan_pars(obs_error = obs_error, estimate_df = estimate_df),
    ...)

  if (algorithm[[1]] == "meanfield") {
    sampling_args$chains <- NULL
    m <- do.call(vb, sampling_args)
  } else {
    m <- do.call(sampling, sampling_args)
  }

  out <- list(model = m, knots = data_knots, y = y, X = X,
    correlation = correlation[[1]], lon = lon, lat = lat, time = time)
  out <- structure(out, class = "rrfield")
}
