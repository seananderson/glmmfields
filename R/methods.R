#' @export
#' @import methods
print.glmmfields <- function(x, pars = c("spatialEffectsKnots", "log_lik"),
                             include = FALSE, ...) {
  print(x$model, pars = pars, include = include, ...)
}

#' Tidy model output
#'
#' @param x Output from [glmmfields()]
#' @param ... Other arguments
#' @export
#' @rdname tidy
tidy <- function(x, ...) {
  UseMethod("tidy")
}

#' @importFrom broom tidy
#' @export
#' @rdname tidy
tidy.glmmfields <- function(x, ...) {
  broom::tidyMCMC(x$model, ...)
}

#' Return LOO information criteria
#'
#' Extract the LOOIC (leave-one-out information criterion) using
#' [loo::loo()].
#'
#' @param object Output from [glmmfields()].
#'   Must be fit with `save_log_lik = TRUE`, which is *not* the default.
#' @param cores Number of cores to use for parallelization.
#'
#' @export
#' @examples
#' \donttest{
#' set.seed(1)
#' s <- sim_glmmfields(n_draws = 12, n_knots = 12, gp_theta = 1.5,
#' gp_sigma = 0.2, sd_obs = 0.2)
#' options(mc.cores = parallel::detectCores()) # for parallel processing
#'
#' # save_log_lik defaults to FALSE to save space but is needed for loo():
#' m <- glmmfields(y ~ 0, time = "time",
#'  lat = "lat", lon = "lon", data = s$dat,
#'  nknots = 12, iter = 1000, chains = 4, seed = 1,
#'  save_log_lik = TRUE)
#' loo(m)
#' }
#' @rdname loo
loo.glmmfields <- function(object, cores = getOption("mc.cores", 1L)) {
  log_lik <- loo::extract_log_lik(object$model, merge_chains = FALSE)
  rel_eff <- loo::relative_eff(exp(log_lik), cores = cores)
  loo::loo.array(log_lik,
    r_eff = rel_eff,
    cores = cores,
    save_psis = FALSE)
}

#' @name loo
#' @rdname loo
#' @export
#' @importFrom loo loo
NULL
