#' Format data for fitting an rrfield model
#'
#' @param data A data frame to be formatted
#' @param y A numeric vector of the response
#' @param X A matrix of the predictors
#' @param time A character object giving the name of the time column
#' @param lon A character object giving the name of the longitude column
#' @param lat A character object giving the name of the latitude column
#' @param nknots The number of knots
#' @param covariance The type of covariance function
#'
#' @export
format_data <- function(data, y, X, time, lon = "lon", lat = "lat", nknots = 25L,
  covariance = "squared-exponential") {

  knots = cluster::pam(data[, c(lon, lat)], nknots)$medoids

  distKnots = as.matrix(dist(knots))

  # Calculate distance from knots to grid
  distAll = as.matrix(stats::dist(rbind(data[, c(lon, lat)], knots)))
  nLocs = nrow(data)

  if (covariance[[1]] == "squared-exponential") {
    distKnots = distKnots^2 # squared distances
    distAll = distAll^2 # squared distances
  }

  # this is the transpose of the lower left corner
  distKnots21 = t(distAll[-c(1:nLocs), -c((nLocs + 1):ncol(distAll))])

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
    distKnots = distKnots,
    distKnots21 = distKnots21,
    X = X,
    nCov = ncol(X))
  list(spatglm_data = spatglm_data, knots = knots)
}
