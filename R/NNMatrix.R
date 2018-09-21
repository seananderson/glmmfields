# Code developed by Lu Zhang licensed under BSD-3 (3-clause), https://opensource.org/licenses/BSD-3-Clause

#' distance matrix for location i and its neighbors
#' @param i The index of point to use
#' @param neighbor_index The indices of neighbors
#' @param s the locations of neighbors
i_dist <- function(i, neighbor_index, s){
  dist(s[c(i, neighbor_index[[i - 1]]), ])
}

#' Function
#'
#'  @param ind
#'  @param ind_distM_d
get_NN_distM <- function (ind, ind_distM_d) {
  if (ind < M ){l = ind } else {l = M};
  M_i <- rep(0, M * (M - 1) / 2);
  if (l == 1) {}
  else{
    M_i[1: (l * (l - 1) / 2)] <-
      c(ind_distM_d[[ind]])[(l + 1): (l * (l + 1) / 2)]
  }
  return(M_i)
}

#' Function to get nearest neighbor distances
#'
#'  @param ind
#'  @param ind_dist_M_d
get_NN_dist <- function (ind, ind_distM_d) {
  if (ind < M ){l = ind } else {l = M};
  D_i <- rep(0, M);
  D_i[1:l] <- c(ind_distM_d[[ind]])[1:l]
  return(D_i)
}

#' Function to get nearest neighbor indices
#'
#'  @param ind
#'  @param ind_dist_M_i
get_NN_ind <- function (ind, ind_distM_i) {
  if (ind < M ){l = ind } else {l = M};
  D_i <- rep(0, M);
  D_i[1:l] <- c(ind_distM_i[[ind]])[1:l]
  return(D_i)
}


#' NNMatrix: A wrapper of spConjNNGP to build Nearest Neighbor matrics ####
#'
#' @param coords: An n x 2 matrix of the observation coordinates in R^2
#' @param n.neighbors: Number of neighbors used in the NNGP.
#' @param n.omp.threads: A positive integer indicating the number of threads to use
#' for SMP parallel processing.
#' @param search.type:   a quoted keyword that specifies type of nearest neighbor
#' search algorithm. Supported method key words are: "tree" and "brute"both will
#' yield the same solution but "tree" should be much faster.
#' @importFrom spNNGP spConjNNGP
NNMatrix <- function(coords, n.neighbors, n.omp.threads = 2,
  search.type = "tree"){
  N <- nrow(coords)
  m.c <- spNNGP::spConjNNGP(rep(0, N) ~ 1, coords = coords,
    n.neighbors = n.neighbors, k.fold = 1,
    theta.alpha = c("phi" = 5, "alpha" = 0.5),
    n.omp.threads = n.omp.threads,
    search.type = search.type,
    return.neighbors = T, sigma.sq.IG = c(2, 1),
    cov.model = "exponential", verbose = F)

  NN_ind <- t(sapply(1: (N - 1), get_NN_ind, m.c$n.indx[-1]))
  neighbor_dist <- sapply(2:N, i_dist, m.c$n.indx[-1], m.c$coords.ord)
  NN_distM <- t(sapply(1: (N - 1), get_NN_distM, neighbor_dist))
  NN_dist <- t(sapply(1: (N - 1), get_NN_dist, neighbor_dist))

  return(list(ord = m.c$ord, coords.ord = m.c$coords.ord,
    NN_ind = NN_ind, NN_distM = NN_distM, NN_dist = NN_dist))
}
