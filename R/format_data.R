#' @export
format_data <- function(data, y, X, time, lon = "lon", lat = "lat", nknots = 25L,
  correlation = "gaussian") {

  knots = cluster::pam(data[, c(lon, lat)], nknots)$medoids

  distKnots = as.matrix(dist(knots))

  # Calculate distance from knots to grid
  distAll = as.matrix(stats::dist(rbind(data[, c(lon, lat)], knots)))
  nLocs = nrow(data)

  if (correlation[[1]] == "gaussian") {
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
