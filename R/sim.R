#' Simulate a random field with a MVT distribution
#'
#' @param n_knots The number of knots
#' @param n_draws The number of draws (for example, the number of years)
#' @param gp_theta The Gaussian Process scale parameter
#' @param gp_sigma The Gaussian Process variance parameter
#' @param mvt Logical: MVT? (vs. MVN)
#' @param df The degrees of freedom parameter for the MVT distribution
#' @param seed The random seed value
#' @param n_data_points The number of data points per draw
#' @param sd_obs The observation process scale parameter
#' @param covariance The covariance function of the Gaussian process
#'   ("squared-exponential", "exponential", "matern")
#' @param matern_kappa The optional matern parameter. Can be 1.5 or 2.5. Values
#'   of 0.5 equivalent to exponential model.
#' @param obs_error The observation error distribution
#' @param B A vector of parameters. The first element is the intercept
#' @param phi The auto regressive parameter on the mean of the random field knots
#' @param X The model matrix
#' @param g Grid of points
#'
#' @export
#' @importFrom ggplot2 ggplot aes_string facet_wrap geom_point scale_color_gradient2
sim_glmmfields <- function(n_knots = 15, n_draws = 10, gp_theta = 0.5,
                           gp_sigma = 0.2, mvt = TRUE, df = 1e6,
                           seed = NULL, n_data_points = 100,
                           sd_obs = 0.1,
  covariance = c("squared-exponential", "exponential", "matern"),
                           matern_kappa = 0.5,
                           obs_error = c("normal", "gamma", "poisson", "nb2", "binomial", "lognormal"),
                           B = c(0), phi = 0, X = rep(1, n_draws * n_data_points),
                           g = data.frame(
                             lon = runif(n_data_points, 0, 10),
                             lat = runif(n_data_points, 0, 10)
                           )) {

  obs_error <- match.arg(obs_error)
  covariance <- match.arg(covariance)

  station_id <- seq_len(n_data_points)
  n_pts <- nrow(g)

  if (!is.null(seed)) {
    set.seed(seed)
  }

  # cluster analysis to determine knot locations
  knots <- jitter(cluster::pam(g, n_knots)$medoids)
  distKnots <- as.matrix(dist(knots))

  if (covariance == "matern") {
    if (matern_kappa %in% c(1.5, 2.5) == FALSE) {
      matern_kappa <- 0.5
      covariance[[1]] <- "exponential"
    }
    else {
      if (matern_kappa == 1.5) {
        dist_knots_sq <- distKnots # NOT squared distances despite name
        transformed_dist <- sqrt(3) * dist_knots_sq / gp_theta
        cor_knots <- (1 + transformed_dist) * exp(-transformed_dist)
      }
      if (matern_kappa == 2.5) {
        dist_knots_sq <- distKnots # NOT squared distances despite name
        transformed_dist <- sqrt(5) * dist_knots_sq / gp_theta
        cor_knots <- (1 + transformed_dist + (transformed_dist^2) / 3) *
          exp(-transformed_dist)
      }
    }
  }
  if (covariance == "squared-exponential") {
    dist_knots_sq <- distKnots^2 # squared distances
    cor_knots <- exp(-dist_knots_sq / (2 * gp_theta^2))
  }
  if (covariance == "exponential") {
    dist_knots_sq <- distKnots # NOT squared distances despite name
    cor_knots <- exp(-dist_knots_sq / (gp_theta))
  }

  sigma_knots <- gp_sigma^2 * cor_knots
  invsigma_knots <- base::solve(sigma_knots)

  # this is the transpose of the lower left corner
  if (covariance == "squared-exponential") {
    # calculate distance from knots to grid
    dist_all <- as.matrix(dist(rbind(g, knots)))^2
    dist_knots21_sq <- t(
      dist_all[-c(seq_len(n_pts)), -c((n_pts + 1):ncol(dist_all))]
    )
    sigma21 <- gp_sigma^2 * exp(-dist_knots21_sq / (2 * gp_theta^2))
  }
  if (covariance == "exponential") {
    # calculate distance from knots to grid
    dist_all <- as.matrix(dist(rbind(g, knots)))
    dist_knots21_sq <- t( # NOT squared distances despite name
      dist_all[-c(seq_len(n_pts)), -c((n_pts + 1):ncol(dist_all))]
    )
    sigma21 <- gp_sigma^2 * exp(-dist_knots21_sq / (gp_theta))
  }
  if (covariance[[1]] == "matern") {
    # calculate distance from knots to grid
    dist_all <- as.matrix(dist(rbind(g, knots)))
    dist_knots21_sq <- t( # NOT squared distances despite name
      dist_all[-c(seq_len(n_pts)), -c((n_pts + 1):ncol(dist_all))]
    )
    if (matern_kappa == 1.5) {
      transformed_dist <- sqrt(3) * dist_knots21_sq / gp_theta
      sigma21 <- gp_sigma^2 * (1 + transformed_dist) * exp(-transformed_dist)
    }
    if (matern_kappa == 2.5) {
      transformed_dist <- sqrt(5) * dist_knots21_sq / gp_theta
      sigma21 <- gp_sigma^2 * (1 + transformed_dist + (transformed_dist^2) / 3) *
        exp(-transformed_dist)
    }
  }

  # generate vector of random effects
  # each 'draw' here is hypothetical draw from posterior
  # initialize:
  re_knots <- matrix(ncol = n_knots, nrow = n_draws)
  if (mvt) {
    re_knots[1, ] <- mvtnorm::rmvt(1, sigma = sigma_knots, df = df)
  }
  if (!mvt) {
    re_knots[1, ] <- mvtnorm::rmvnorm(1, sigma = sigma_knots)
  }
  # potentially with AR process:
  if (n_draws > 1) {
    for (i in seq(2, n_draws)) {
      if (mvt) {
        re_knots[i, ] <- mvtnorm::rmvt(1,
          # delta = ar * (re_knots[i - 1, ] - mean(re_knots[i - 1, ])),
          delta = phi * (re_knots[i - 1, ]),
          sigma = sigma_knots, df = df
        )
      }
      if (!mvt) {
        re_knots[i, ] <- mvtnorm::rmvnorm(1,
          # mean = ar * (re_knots[i - 1, ] - mean(re_knots[i - 1, ])),
          mean = phi * (re_knots[i - 1, ]),
          sigma = sigma_knots
        )
      }
    }
  }

  # project random effects to locations of the data
  proj <- t((sigma21 %*% invsigma_knots) %*% t(re_knots))

  # multiply coefficients by the design matrix
  eta <- as.vector(B %*% t(X), mode = "double")
  eta_mat <- matrix(eta, nrow = nrow(proj), byrow = TRUE)

  # add the observation process:
  N <- ncol(proj) * nrow(proj)
  if (obs_error == "normal") {
    y <- proj + eta_mat + matrix(
      data = stats::rnorm(N, 0, sd_obs),
      ncol = ncol(proj), nrow = nrow(proj)
    )
  }
  if (obs_error == "nb2") {
    y <- matrix(
      data = stats::rnbinom(N, mu = exp(proj + eta_mat), size = sd_obs),
      ncol = ncol(proj), nrow = nrow(proj)
    )
  }
  if (obs_error == "gamma") {
    gamma_a <- 1 / (sd_obs^2) # sd_obs means CV here
    gamma_b <- gamma_a / exp(proj + eta_mat)
    y <- matrix(
      data = stats::rgamma(N, shape = gamma_a, rate = gamma_b),
      ncol = ncol(proj), nrow = nrow(proj)
    )
  }
  if (obs_error == "binomial") { # plogis = inverse_logit
    y <- matrix(data = stats::rbinom(N,
      size = 1,
      prob = stats::plogis(proj + eta_mat)
    ), ncol = ncol(proj), nrow = nrow(proj))
  }
  if (obs_error == "poisson") {
    y <- matrix(
      data = stats::rpois(N, lambda = exp(proj + eta_mat)),
      ncol = ncol(proj), nrow = nrow(proj)
    )
  }
  if (obs_error == "lognormal") {
    y <- matrix(
      data = stats::rlnorm(N, meanlog = proj + eta_mat, sdlog = sd_obs),
      ncol = ncol(proj), nrow = nrow(proj)
    )
  }
  # Reshape for output
  out <- reshape2::melt(y)
  names(out) <- c("time", "pt", "y")
  out <- dplyr::arrange_(out, "time", "pt")
  out$lon <- rep(g$lon, n_draws)
  out$lat <- rep(g$lat, n_draws)
  out$station_id <- rep(station_id, n_draws)

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
    dat = out
  )
}
