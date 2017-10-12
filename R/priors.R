#' Student-t and half-t priors
#'
#' Student-t and half-t priors. Note that this can be used to represent an
#' effectively normal distribution prior by setting the first argument (the
#' degrees of freedom parameter) to a large value (roughly 50 or above).
#'
#' @param df Degrees of freedom parameter
#' @param location Location parameter
#' @param scale Scale parameter
#' @export
#' @rdname priors
#' @examples
#' student_t(3, 0, 1)
student_t <- function(df = 3, location = 0, scale = 1) {
  stopifnot(is.numeric(df), is.numeric(location), is.numeric(scale),
    df >= 1, scale > 0)
  list(dist = "t", df = df, location = location, scale = scale)
}

#' @export
#' @rdname priors
#' @examples
#' half_t(3, 0, 1)
half_t <- function(df = 3, location = 0, scale = 1) {
  if(location != 0) warning("half-t location != 0")
  ht <- student_t(df, location, scale)
  ht[[1]] <- "half-t"
  ht
}

parse_t_prior <- function(x) {
  as.vector(unlist(x)[-1], mode = "numeric")
}
