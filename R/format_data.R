#' Format data for fitting a glmmfields model
#'
#' @param data A data frame to be formatted
#' @param y A numeric vector of the response
#' @param X A matrix of the predictors
#' @param time A character object giving the name of the time column
#' @param lon A character object giving the name of the longitude column
#' @param lat A character object giving the name of the latitude column
#' @param station A numeric vector giving the integer ID of the station
#' @param nknots The number of knots
#' @param covariance The type of covariance function
#' @param fixed_intercept Should the intercept be fixed?
#' @param cluster The type of clustering algorithm used to determine the not locations.
#'   \code{"pam"} = \code{\link[cluster]{pam}}. \code{kmeans} is faster for large datasets.
format_data <- function(data, y, X, time,
                        lon = "lon", lat = "lat",
                        station = NULL, nknots = 25L,
                        covariance = c("squared-exponential",
                          "exponential", "matern"),
                        fixed_intercept = FALSE,
                        cluster = c("pam", "kmeans")) {
  data <- as.data.frame(data)
  cluster <- match.arg(cluster)
  covariance <- match.arg(covariance)

  if (is.null(time)) {
    data$time <- 1
    time <- "time"
  }
  yearID <- as.numeric(data[, time, drop = TRUE])
  yearID <- yearID - min(yearID) + 1 # convert to 1, ..., nT
  if (is.null(station)) {
    stationID <- seq_len(nrow(data))
  } else {
    stationID <- as.numeric(forcats::as_factor(data[, station, drop = TRUE]))
  }
  data$stationID <- stationID

  # if stationID is duplicated, perform clustering on the subset of data
  if (length(unique(stationID)) < length(stationID)) {
    first_instance <- which(!duplicated(stationID))

    if (cluster == "pam") {
      knots <- cluster::pam(data[first_instance, c(lon, lat), drop = FALSE], nknots)$medoids
    } else {
      if (cluster == "kmeans") {
        knots <- stats::kmeans(data[first_instance, c(lon, lat), drop = FALSE], nknots)$centers
      } else {
        # each point is unique, predictive process not used
        knots <- data[first_instance, c(lon, lat)]
      }
    }

    distKnots <- as.matrix(dist(knots))
    ix <- sort(data[first_instance, "stationID"], index.return = TRUE)$ix

    # Calculate distance from knots to grid
    distAll <- as.matrix(stats::dist(rbind(data[first_instance, c(lon, lat)][ix, ], knots)))
    nLocs <- length(first_instance)
  } else {
    if (cluster == "pam") {
      knots <- cluster::pam(data[, c(lon, lat), drop = FALSE], nknots)$medoids
    } else {
      if (cluster == "kmeans") {
        knots <- stats::kmeans(data[, c(lon, lat), drop = FALSE], nknots)$centers
      } else {
        # each point is unique, predictive process not used
        knots <- data[, c(lon, lat), drop = FALSE]
      }
    }
    distKnots <- as.matrix(dist(knots))

    # Calculate distance from knots to grid
    distAll <- as.matrix(stats::dist(rbind(data[, c(lon, lat), drop = FALSE], knots)))
    nLocs <- nrow(data)
  }

  if (covariance == "squared-exponential") {
    distKnots <- distKnots^2 # squared distances
    distAll <- distAll^2 # squared distances
  }

  # this is the transpose of the lower left corner
  distKnots21 <- t(distAll[-seq_len(nLocs), -seq(nLocs + 1, ncol(distAll))])

  # create list for Stan
  spatglm_data <- list(
    nKnots = nknots,
    nLocs = nLocs,
    nT = max(yearID),
    N = length(y),
    stationID = stationID,
    yearID = yearID,
    y = y,
    distKnots = distKnots,
    distKnots21 = distKnots21,
    X = X,
    nCov = if(fixed_intercept) 0 else ncol(X)
  )
  list(spatglm_data = spatglm_data, knots = knots)
}
