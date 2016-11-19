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
  spatglm_data
}

stan_pars <- function() {
  c(
    "df",
    "yearEffects",
    "sigma",
    "gp_sigma",
    "gp_scale",
    "year_sigma",
#    "ar",
    "spatialEffectsKnots"
  )
}

parse_t_prior <- function(x) {
  as.vector(unlist(x)[-1], mode = "numeric")
}

#' Fit a robust spatiotemporal random fields model
#'
#' @export
#' @importFrom rstanarm student_t normal
rrfield <- function(formula, data, time, lon, lat, nknots = 25L,
  prior_gp_scale = rstanarm::student_t(3, 0, 10),
  prior_gp_sigma = rstanarm::student_t(3, 0, 2),
  prior_sigma = rstanarm::student_t(3, 0, 2),
  prior_ar = rstanarm::student_t(100, 0, 1),
  estimate_df = TRUE,
  obs_error = c("normal","gamma"),
  correlation = c("gaussian", "exponential"),
  ...) {

  mf <- model.frame(formula, data)
  X <- model.matrix(formula, mf)

  # user inputs raw data. this function formats it for STAN
  stan_data <- format_data(data = data, y = mf[,1], X = X, time = time,
    lon = lon, lat = lat, nknots = nknots)

  gauss_cor <- switch(correlation[[1]], gaussian = 1, exponential = 0, 1)

  stan_data <- c(stan_data,
    list(prior_gp_scale = parse_t_prior(prior_gp_scale),
      prior_gp_sigma = parse_t_prior(prior_gp_sigma),
      prior_sigma = parse_t_prior(prior_sigma),
      prior_ar = parse_t_prior(prior_ar),
      gauss_cor = gauss_cor))

  m <- rstan::sampling(stanmodels$mvt_norm_yr_ar1, data = stan_data,
    pars = stan_pars(), ...)
  m
}
