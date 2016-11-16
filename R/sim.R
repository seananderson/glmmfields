#' @export

sim_mvt_rf <- function(n_knots = 15, n_draws = 10, gp_scale = 0.5,
  sigma_t = 0.2, mvt = TRUE, df = 3, seed = NULL, nDataPoints = 100) {

  g <- data.frame(lon = runif(nDataPoints, 0, 10),
    lat = runif(nDataPoints, 0, 10))
  n_pts <- nrow(g)

  if (!is.null(seed)) {
    set.seed(seed)
  }
  # cluster analysis to determine knot locations
  knots <- jitter(cluster::pam(g, n_knots)$medoids)
  distKnots <- as.matrix(dist(knots))
  dist_knots_sq <- distKnots^2 # squared distances

  cor_knots <- exp(-gp_scale * dist_knots_sq)
  sigma_knots <- sigma_t * sigma_t * cor_knots
  invsigma_knots <- base::solve(sigma_knots)

  # calculate distance from knots to grid
  dist_all <- as.matrix(dist(rbind(g, knots)))^2

  # this is the transpose of the lower left corner
  dist_knots21_sq <- t(
    dist_all[-c(seq_len(n_pts)), -c((n_pts + 1):ncol(dist_all))])
  sigma21 <- exp(-gp_scale * dist_knots21_sq) * sigma_t * sigma_t

  # generate vector of random effects
  # each 'draw' here is hypothetical draw from posterior
  if (mvt)
    re_knots <- mvtnorm::rmvt(n_draws, sigma = sigma_knots, df = df)
  if (!mvt)
    re_knots <- mvtnorm::rmvnorm(n_draws, sigma = sigma_knots)

  # project random effects to locations of the data
  proj <- t((sigma21 %*% invsigma_knots) %*% t(re_knots))

  return(
    list(
      knots = knots,
      re_knots = re_knots,
      proj = proj,
      dist_knots_sq = dist_knots_sq,
      dist_knots21_sq = dist_knots21_sq,
      sigma_knots = sigma_knots,
      g = g
    )
  )
}
