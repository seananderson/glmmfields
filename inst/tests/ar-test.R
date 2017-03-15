library(rrfields)
set.seed(1)

gp_sigma <- 0.2
sigma <- 0.01
df <- 6
gp_scale <- 1.8
n_draws <- 20
nknots <- 10
ar <- 1
B <- vector(mode = "double", length = n_draws)
B[1] <- 6
year_sigma <- 0.4
for (i in 2:length(B)) {
  B[i] <- B[i-1] + rnorm(1, 0, year_sigma) # random walk
}

s <- sim_rrfield(df = df, n_draws = n_draws, gp_scale = gp_scale,
  gp_sigma = gp_sigma, sd_obs = sigma, n_knots = nknots, ar = ar,
  B = B, X = model.matrix(~ a - 1, data.frame(a = gl(n_draws, 100))))
# print(s$plot)

# year effects:
plot(B, pch = 19)

# Empirical mean values at each year:
library(dplyr)
means <- group_by(s$dat, time) %>%
  summarise(m = mean(y))
points(means$m, col = "red")

# mean of each random field draw:
re_knots <- apply(s$re_knots, 1, mean)
print(re_knots)

# this gets closer... but maybe I'm still missing something?
points(B + re_knots, col = "blue")
