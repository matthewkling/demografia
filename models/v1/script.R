
# this version includes a very basic set of markov transitions



library(tidyverse)
library(stranger)
library(cmdstanr)


## simulate dataset ########

stages <- c("s", "b", "j", "a") # seed-baby-juvenile-adult
k <- length(stages)

# true projection matrix
p <- matrix(0, k, k, dimnames = list(stages, stages))
p["b", "s"] <- .003
p["b", "b"] <- .5
p["j", "b"] <- .01
p["j", "j"] <- .98
p["a", "a"] <- .995
p <- rbind(p, x = apply(p, 2, function(x) 1 - sum(x)))
stages1 <- rownames(p)

fecundity <- 400



# simulate
n_i <- 10000 # number of individuals

d0 <- sample(stages, n_i, replace = T)
d1 <- apply(p[,d0], 2, function(x) sample(stages1, 1, prob = x))

d <- tibble(s0 = sample(stages, n_i, replace = T),
            s1 = apply(p[,s0], 2, function(x) sample(stages1, 1, prob = x))) %>%
  count(s0, s1) %>%
  full_join(expand_grid(s0 = colnames(p),
                        s1 = rownames(p))) %>%
  mutate(n = ifelse(is.na(n), 0, n)) %>%
  spread(s0, n) %>%
  column_to_rownames("s1") %>%
  as.matrix()
d <- d[rownames(p), colnames(p)]

sum(d[,"a"]) * fecundity



# hyperparameters for transition probability priors
alpha <- c(s = c(0, 1, 0, 0, 10),
           b = c(0, 2, 1, 0, 2),
           j = c(0, 0, 10, 1, 1),
           a = c(0, 0, 0, 10, 1)) %>%
  matrix(ncol = 4, byrow = F)
alpha[alpha == 0] <- .1


## fit stan model ############

# data inputs
dat <- list(
  K = ncol(d),
  y = t(d),
  alpha = t(alpha)
)

# compile
model <- cmdstan_model("models/v1/model.stan")

# run MCMC
fit <- model$sample(
  data = dat,
  chains = 2
)



## assess fit #########

f <- fit$draws() %>% 
  apply(3, median) %>%
  as.data.frame() %>%
  rownames_to_column("param") %>%
  rename(value = ".") %>%
  filter(param != "lp__",
         str_detect(param, "theta")) %>%
  mutate(param = str_remove_all(param, "theta|\\[|\\]")) %>%
  separate(param, c("s0", "s1"), sep = ",") %>%
  mutate(s0 = colnames(p)[as.integer(s0)],
         s1 = rownames(p)[as.integer(s1)]) %>%
  spread(s0, value) %>%
  column_to_rownames("s1") %>%
  as.matrix()
f <- f[rownames(p), colnames(p)]

plot(p, f)
