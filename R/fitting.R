#' @export

format_data <- function(data, y, X, time, lon = "lon", lat = "lat", nknots = 25L) {

  knots = cluster::pam(data[, c(lon, lat)], nknots)$medoids

  distKnots = as.matrix(dist(knots))
  distKnotsSq = distKnots^2 # squared distances

  # Calculate distance from knots to grid
  distAll = as.matrix(stats::dist(rbind(data[, c(lon, lat)], knots)))^2
  nLocs = nrow(data)
  # this is the transpose of the lower left corner
  distKnots21Sq = t(distAll[-c(1:nLocs), -c((nLocs + 1):ncol(distAll))])

  yearID = as.numeric(as.factor(data[,time]))
  stationID = seq(1, nrow(data))

  # create list for STAN
  spatglm_data = list(
    nKnots = nknots,
    nLocs = nLocs,
    nT = length(unique(yearID)),
    N = length(y),
    stationID = stationID,
    yearID = yearID,
    y = y,
    distKnotsSq = distKnotsSq,
    distKnots21Sq = distKnots21Sq,
    X = X,
    nCov = ncol(X))
  return(list(spatglm_data = spatglm_data, knots = knots))
}

stan_pars <- function(obs_error) {
  c(
    "df",
#    "yearEffects",
    "gp_sigma",
    "gp_scale",
#    "year_sigma",
#    "ar",
    "B",
    "spatialEffectsKnots",
    ifelse(obs_error==1,"sigma","CV")
  )
}

parse_t_prior <- function(x) {
  as.vector(unlist(x)[-1], mode = "numeric")
}

#' Fit a robust spatiotemporal random fields model
#'
#' @export
#' @importFrom rstanarm student_t normal
#' @importFrom rstan sampling
#' @import Rcpp
rrfield <- function(formula, data, time, lon, lat, nknots = 25L,
  prior_gp_scale = student_t(3, 0, 10),
  prior_gp_sigma = student_t(3, 0, 2),
  prior_sigma = student_t(3, 0, 2),
  fixed_df_value = 5,
  estimate_df = TRUE,
  obs_error = c("normal", "gamma"),
  correlation = c("gaussian", "exponential"),
  ...) {

  mf <- model.frame(formula, data)
  X <- model.matrix(formula, mf)
  y <- model.response(mf, "numeric")

  # user inputs raw data. this function formats it for STAN
  data_list <- format_data(data = data, y = y, X = X, time = time,
    lon = lon, lat = lat, nknots = nknots)
  stan_data = data_list$spatglm_data
  data_knots = data_list$knots

  gauss_cor <- switch(correlation[[1]], gaussian = 1L, exponential = 0L, 1L)

  obs_model <- switch(obs_error[[1]], normal = 1L, gamma = 0L, 1L)

  stan_data <- c(stan_data,
    list(prior_gp_scale = parse_t_prior(prior_gp_scale),
      prior_gp_sigma = parse_t_prior(prior_gp_sigma),
      prior_sigma = parse_t_prior(prior_sigma),
      gauss_cor = gauss_cor,
      obs_model = obs_model,
      est_df = as.integer(estimate_df),
      fixed_df_value = fixed_df_value))

  sampling_args <- list(
    object = stanmodels$rrfield,
    data = stan_data,
    pars = stan_pars(obs_model),
    ...)

  m <- do.call(sampling, sampling_args)
  return(list(model = m, knots = data_knots))
}
