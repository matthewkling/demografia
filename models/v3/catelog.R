
library(cmdstanr)


### version directly from stan docs #####

dat <- list(
  K = 3,
  N = 10,
  D = 2,
  y = sample(1:3, 10, replace = T),
  x = matrix(runif(10*2), 10, 2)
)

model <- cmdstan_model("models/v3/catelog.stan")

# run MCMC
fit <- model$sample(
  data = dat,
  chains = 2
)

fit



### implementation with multiple outcomes per data observation ####

dat <- list(
  K = 3,
  N = 10,
  D = 2,
  y = matrix(sample(1:3, 10*3, replace = T), 10),
  x = matrix(runif(10*2), 10, 2)
)

model <- cmdstan_model("models/v3/catelog2.stan")

# run MCMC
fit <- model$sample(
  data = dat,
  chains = 2
)

fit




### implementation with full transition matrix ####

dat <- list(
  K = 3,
  N = 10,
  D = 2,
  y = array(sample(1:3, 3*10*3, replace = T), dim = c(3, 10, 3)),
  x = matrix(runif(10*2), 10, 2)
)

model <- cmdstan_model("models/v3/catelog3.stan")

# run MCMC
fit <- model$sample(
  data = dat,
  chains = 2
)

fit
