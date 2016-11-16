

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
    distKnots21Sq = distKnots21Sq,
    x = rep(0, length(Y))
  )
  spatglm_data
}

#' @export
stan_pars <- function() {
  spatglm_pars = c(
    "scaledf",
    "yearEffects",
    "sigma",
    "gp_sigmaSq",
    "gp_scale",
    "year_sigma",
    "ar",
    "spatialEffectsKnots"
  )
  spatglm_pars
}

#' @export
rrfield <- function(data, pars = stan_pars(), ...) {
  m <- rstan::sampling(stanmodels$mvt_norm_yr_ar1, data = data, pars = pars, ...)
  m
}
