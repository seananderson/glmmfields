

#' @export

format_data <- function(data, y, time, lon = "lon", lat = "lat", nKnots = 25L) {

  knots = cluster::pam(data[, c(lon, lat)], nKnots)$medoids

  distKnots = as.matrix(dist(knots))
  distKnotsSq = distKnots^2 # squared distances

  # Calculate distance from knots to grid
  distAll = as.matrix(stats::dist(rbind(data[, c(lon, lat)], knots)))^2
  nLocs = nrow(data)
  # this is the transpose of the lower left corner
  distKnots21Sq = t(distAll[-c(1:nLocs), -c((nLocs + 1):ncol(distAll))])

  Y = as.numeric(data[,y])
  yearID = as.numeric(as.factor(data[,time]))
  stationID = seq(1, nrow(data))

  # create list for STAN
  spatglm_data = list(
    nKnots = nKnots,
    nLocs = nLocs,
    nT = length(unique(yearID)),
    N = length(Y),
    stationID = stationID,
    yearID = yearID,
    y = Y,
    distKnotsSq = distKnotsSq,
    distKnots21Sq = distKnots21Sq)
  spatglm_data
}

#' @export
stan_pars <- function() {
  spatglm_pars = c(
    "df",
    "yearEffects",
    "sigma",
    "gp_sigma",
    "gp_scale",
    "year_sigma",
    "ar",
    "spatialEffectsKnots"
  )
  spatglm_pars
}

parse_t_prior <- function(x) {
  as.vector(unlist(x)[-1], mode = "integer")
}

#' @export
#' @importFrom rstanarm student_t normal
rrfield <- function(data, pars = stan_pars(),
  prior_gp_scale = student_t(3, 0, 10),
  prior_gp_sigma = student_t(3, 0, 2),
  prior_sigma = student_t(3, 0, 2),
  prior_ar = student_t(100, 0, 1),
  ...) {

  data <- c(data,
    list(prior_gp_scale = parse_t_prior(prior_gp_scale)),
    list(prior_gp_scale = parse_t_prior(prior_gp_sigma)),
    list(prior_gp_scale = parse_t_prior(prior_sigma)),
    list(prior_gp_scale = parse_t_prior(prior_ar))
  )
  m <- rstan::sampling(stanmodels$mvt_norm_yr_ar1, data = data, pars = pars, ...)
  m
}
