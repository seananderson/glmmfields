library(rstan)
library(tweedie)
rstan_options(auto_write = TRUE)

set.seed(1234)
N <- 400
Y <- tweedie::rtweedie(N, power = 1.5, mu = 2.8, phi = 1.2)
# M = number of series expansions, more = more accurate
# see Dunn and Smyth (2005) Table 1, depends on phi and power (theta)
M <- 20

plot(Y)
data <- list(N = N, M = M, Y = Y)
fit <- stan("inst/tests/tweedie.stan", data = data, iter = 300, chains = 3, cores = 3)
fit
bayesplot::mcmc_areas(as.matrix(fit), pars = c("mu", "phi", "theta"))
