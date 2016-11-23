#' Simulate a random field with a MVT distribution
#'
#' @param n_knots The number of knots
#' @param n_draws The number of draws (for example, the number of years)
#' @param gp_scale The Gaussian Process scale parameter
#' @param gp_sigma The Gaussian Process sigma parameter
#' @param mvt Logical: MVT? (vs. MVN)
#' @param df The degrees of freedom parameter for the MVT distribution
#' @param seed The random seed value
#' @param nDataPoints The number of data points
#' @param sd_obs The observation process scale parameter
#' @param correlation The correlation or covariance function
#' @param obs_error The observation error distribution
#' @param b0 The intercept parameter
#'
#' @export
#' @importFrom ggplot2 ggplot aes_string facet_wrap geom_point scale_color_gradient2
sim_rrfield <- function(n_knots = 15, n_draws = 10, gp_scale = 0.5,
  gp_sigma = 0.2, mvt = TRUE, df = 4, seed = NULL, nDataPoints = 100,
  sd_obs = 0.1, correlation = "gaussian",
  obs_error = c("normal", "gamma", "nb2"), B = c(0),
  X = rep(1, n_draws * nDataPoints)) {

  g <- data.frame(lon = runif(nDataPoints, 0, 10),
    lat = runif(nDataPoints, 0, 10))
  n_pts <- nrow(g)

  if (!is.null(seed)) {
    set.seed(seed)
  }

  # cluster analysis to determine knot locations
  knots <- jitter(cluster::pam(g, n_knots)$medoids)
  distKnots <- as.matrix(dist(knots))


  if (!correlation[[1]] %in% c("gaussian", "exponential")) {
    stop(paste(correlation[[1]], "not implemented"))
  }
  if (correlation[[1]] == "gaussian") {
    dist_knots_sq <- distKnots^2 # squared distances
    cor_knots <- exp(-dist_knots_sq / (2 * gp_scale^2))
  }
  if (correlation[[1]] == "exponential") {
    dist_knots_sq <- distKnots # NOT squared distances despite name
    cor_knots <- exp(-dist_knots_sq / (gp_scale))
  }

  sigma_knots <- gp_sigma^2 * cor_knots
  invsigma_knots <- base::solve(sigma_knots)

  # this is the transpose of the lower left corner
  if (correlation[[1]] == "gaussian") {
    # calculate distance from knots to grid
    dist_all <- as.matrix(dist(rbind(g, knots)))^2
    dist_knots21_sq <- t(
      dist_all[-c(seq_len(n_pts)), -c((n_pts + 1):ncol(dist_all))])
    sigma21 <- gp_sigma^2 * exp(-dist_knots21_sq / (2 * gp_scale^2))
  }
  if (correlation[[1]] == "exponential") {
    # calculate distance from knots to grid
    dist_all <- as.matrix(dist(rbind(g, knots)))
    dist_knots21_sq <- t( # NOT squared distances despite name
      dist_all[-c(seq_len(n_pts)), -c((n_pts + 1):ncol(dist_all))])
    sigma21 <- gp_sigma^2 * exp(-dist_knots21_sq / (gp_scale))
  }

  # generate vector of random effects
  # each 'draw' here is hypothetical draw from posterior
  if (mvt)
    re_knots <- mvtnorm::rmvt(n_draws, sigma = sigma_knots, df = df)
  if (!mvt)
    re_knots <- mvtnorm::rmvnorm(n_draws, sigma = sigma_knots)

  # project random effects to locations of the data
  proj <- t((sigma21 %*% invsigma_knots) %*% t(re_knots))

  # multiply coefficients by the design matrix
  eta <- as.vector(B %*% t(X), mode = "double")
  eta_mat <- matrix(eta, nrow = nrow(proj), byrow = TRUE)

  # add the observation process:
  N <- ncol(proj) * nrow(proj)
  if (obs_error[[1]] == "normal") {
    y <- proj + eta_mat + matrix(data = stats::rnorm(N, 0, sd_obs),
      ncol = ncol(proj), nrow = nrow(proj))
  }
  if (obs_error[[1]] == "nb2") {
    proj <- proj + eta_mat
    y <- matrix(data = stats::rnbinom(N, mu = exp(proj), size = sd_obs),
      ncol = ncol(proj), nrow = nrow(proj))
  }
  if (obs_error[[1]] == "gamma") {
    gamma_a = 1/(sd_obs^2) # sd_obs means CV here
    gamma_b = gamma_a/exp(proj + eta_mat)
    y <- matrix(data = stats::rgamma(N, shape = gamma_a, rate = gamma_b),
      ncol = ncol(proj), nrow = nrow(proj))
  }

  # Reshape for output
  out <- reshape2::melt(y)
  names(out) <- c("time", "pt", "y")
  out <- dplyr::arrange_(out, "time", "pt")
  out$lon <- rep(g$lon, n_draws)
  out$lat <- rep(g$lat, n_draws)

  plot <- ggplot(out, aes_string(x = "lon", y = "lat", colour = "y")) +
    facet_wrap(~ time) +
    geom_point(size = 2) +
    scale_color_gradient2()

  list(
    knots = knots,
    re_knots = re_knots,
    proj = proj,
    dist_knots_sq = dist_knots_sq,
    dist_knots21_sq = dist_knots21_sq,
    sigma_knots = sigma_knots,
    g = g,
    plot = plot,
    dat = out)
}
