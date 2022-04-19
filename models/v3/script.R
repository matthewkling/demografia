
# this version adds environmental (also effectively density dependence) effects. 
# substantial reengineering of model required.


library(tidyverse)
library(cmdstanr)


## simulate dataset ########

stages <- c("s", "b", "j", "a") # seed-baby-juvenile-adult
k <- length(stages)

# true projection matrix
p <- matrix(0, k, k, dimnames = list(stages, stages))
p["b", "s"] <- .003
p["b", "b"] <- .5
p["j", "b"] <- .01
p["j", "j"] <- .9
p["a", "j"] <- .05
p["a", "a"] <- .995
p <- rbind(p, x = apply(p, 2, function(x) 1 - sum(x)))
p <- p[2:nrow(p),]
stages1 <- rownames(p)

# fecundity <- 400



# simualation specs
n_i <- 100 # number of individuals per site
n_sites <- 100 # number of sites

# environmental predictors
x <- matrix(runif(n_sites*2), n_sites, 2)
colnames(x) <- c("v1", "v2")

# specify predictor effects
betas <- array(0, 
               dim = c(ncol(p), ncol(x), nrow(p)),
               dimnames = list(colnames(p), colnames(x), rownames(p)))
betas["b", "v1", "b"] <- 3
betas["b", "v2", "b"] <- -2
betas["a", "v1", "x"] <- -1
betas["a", "v2", "a"] <- -1
betas["j", "v1", "j"] <- -.5
betas["j", "v2", "a"] <- 2

# binary indicator of which betas to infer (passed to model)
beta_switch <- betas
beta_switch[] <- as.integer(betas[] != 0)

# transforms
logit <- function(x) log(x / (1-x))
softmax <- function(x) exp(x) / sum(exp(x))
# softmax(log(c(.1, .2, .3, .4))) # demo that log is the inverse of softmax


# outcome counts
y <- array(dim = c(k, k, n_sites),
           dimnames = list(colnames(p),
                           rownames(p),
                           paste0("s", 1:n_sites)))
for(s in 1:n_sites){
  
  xs <- x[s,] # local environment
  ps <- p # to, from
  
  for(i in 1:k){ # from
    for(j in 1:k){ # to
      for(v in 1:ncol(x)){ # var
        b <- betas[i,v,j]
        if(b == 0) next()
        # if(b == -1) stop()
        lp <- log(ps[,i]) # convert to logit scale
        lp[j] = lp[j] + betas[i, v, j] * xs[[v]] # apply effect
        ps[,i] <- softmax(lp) # convert back to probability
      }
    }
  }
  
  d <- tibble(s0 = sample(stages, n_i, replace = T),
              s1 = apply(ps[,s0], 2, function(x) sample(stages1, 1, prob = x))) %>%
    count(s0, s1) %>%
    full_join(expand_grid(s0 = colnames(p),
                          s1 = rownames(p))) %>%
    mutate(n = ifelse(is.na(n), 0, n)) %>%
    spread(s0, n) %>%
    column_to_rownames("s1") %>%
    as.matrix()
  d <- d[rownames(p), colnames(p)]
  y[,,s] <- t(d)
}

# # check that covarites had the desired effect
# plot(x[,1], y["a","x",])
# plot(x[,2], y["j","a",])

alpha_switch <- ceiling(p) %>% t()



## fit stan model ############

dat <- list(
  K = k,
  N = n_sites,
  D = ncol(x),
  y = y, 
  x = x, # site, var
  alpha_switch = alpha_switch, # from, to
  beta_switch = beta_switch # from, var, to
)

# compile
model <- cmdstan_model("models/v3/model.stan")

# run MCMC
fit <- model$sample(
  data = dat,
  # init = function() list(alpha = t(apply(log(p)+1e-10, 2, softmax))),
  chains = 1,
  iter_warmup = 200, iter_sampling = 100, refresh = 30
)

fit




## assess fit #########


## alphas ##

pd <- p %>%
  as.data.frame() %>% as_tibble() %>%
  rownames_to_column("s1") %>%
  mutate(s1 = rownames(p)[as.integer(s1)]) %>%
  gather(s0, true, -s1)

f <- fit$summary() %>%
  filter(str_detect(variable, "alpha_switched")) %>%
  mutate(variable = str_remove_all(variable, "alpha_switched|\\[|\\]")) %>%
  separate(variable, c("s0", "s1"), sep = ",") %>%
  mutate(s0 = colnames(p)[as.integer(s0)],
         s1 = rownames(p)[as.integer(s1)]) %>%
  left_join(pd) %>%
  group_by(s0) %>%
  mutate(pred = softmax(median))


l <- 6
ggplot(f, aes(log(true), log(pred), 
              ymin = logit(q5), ymax = logit(q95),
              color = s0, fill = s1)) + 
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_linerange(size = .25) +
  geom_point(shape = 21, size = 3) +
  coord_cartesian(xlim = c(-l, 0),
                  ylim = c(-l, 0)) +
  theme_minimal()


## betas ##

bd <- betas %>%
  as.data.frame() %>% 
  rownames_to_column("s0") %>% as_tibble() %>%
  gather(par, true, -s0) %>%
  separate(par, c("v", "s1"), sep = "\\.")

f <- fit$summary() %>%
  filter(str_detect(variable, "beta_switched")) %>%
  mutate(variable = str_remove_all(variable, "beta_switched|\\[|\\]")) %>%
  separate(variable, c("s0", "v", "s1"), sep = ",") %>%
  mutate(s0 = colnames(p)[as.integer(s0)],
         s1 = rownames(p)[as.integer(s1)],
         v = colnames(x)[as.integer(v)]) %>%
  left_join(bd)

f %>%
  ggplot(aes(true, median, ymin = q5, ymax = q95)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_pointrange()
