#' @export

format_data <- function(data, y, X, time, lon = "lon", lat = "lat", nKnots = 25L) {

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
    X = X,
    nCov = ncol(X))
  spatglm_data
}

stan_pars <- function() {
  spatglm_pars = c(
  "muZeros",
  "spatialEffects",
  "SigmaKnots",
  "SigmaOffDiag",
  "invSigmaKnots",
  "y_hat",
  "gp_sigmaSq")
}

parse_t_prior <- function(x) {
  as.vector(unlist(x)[-1], mode = "numeric")
}

#' @export
#' @importFrom rstanarm student_t normal
rrfield <- function(formula,
  time,
  lon,
  lat,
  data,
  y,
  X,
  pars = stan_pars(),
  nKnots = 25L,
  prior_gp_scale = rstanarm::student_t(3, 0, 10),
  prior_gp_sigma = rstanarm::student_t(3, 0, 2),
  prior_sigma = rstanarm::student_t(3, 0, 2),
  prior_ar = rstanarm::student_t(100, 0, 1),
  estimate_df = TRUE,
  year_effects = c("fixed","zero"),
  obs_error = c("normal","gamma"),
  correlation = c("gaussian", "exponential"),
  ...) {

  # user inputs raw data. this function formats it for STAN
  data = format_data(data=data, y=y, X=X, time=time, lon=lon, lat=lat, nKnots=nKnots)

  data <- c(data,
    list(prior_gp_scale = parse_t_prior(prior_gp_scale)),
    list(prior_gp_scale = parse_t_prior(prior_gp_sigma)),
    list(prior_gp_scale = parse_t_prior(prior_sigma)),
    list(prior_gp_scale = parse_t_prior(prior_ar))
  )

  m <- rstan::sampling(stanmodels$mvt_norm_yr_ar1, data = data, pars = pars, include=FALSE, ...)
  m
}
