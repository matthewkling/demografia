
# this version adds correction for multiyear survey intervals
# via the project_simplex functions on the R and Stan sides




library(tidyverse)
library(cmdstanr)


## simulate dataset ########

stages <- c("p", "s", "j", "a", "x") # propagule-seedling-juvenile-adult-dead
k <- length(stages)

# true projection matrix -- origin in COLS, outcome in ROWS
p <- matrix(0, k, k, dimnames = list(stages, stages))
p["s", "p"] <- .003
p["s", "s"] <- .5
p["j", "s"] <- .01
p["j", "j"] <- .9
p["a", "j"] <- .05
p["a", "a"] <- .995
p["x",] <- apply(p, 2, function(x) 1 - sum(x))


# fecundity <- 400

# function that converts an annual simplex to a multiyear simplex
project_simplex <- function(p, # projection matrix 
                            t, # number of time steps
                            i){ # life stage to project
  n <- matrix(0, 1, ncol(p))
  n[i] <- 1
  for(j in 1:t) n <- n %*% t(p)
  return(as.vector(n))
}

# simulation specs
n_i <- 1000 # number of individuals per site
n_sites <- 1000 # number of sites

# environmental predictors
x <- matrix(runif(n_sites*2), n_sites, 2)
colnames(x) <- c("v1", "v2")

# specify predictor effects
betas <- array(0, 
               dim = c(ncol(p), ncol(x), nrow(p)),
               dimnames = list(colnames(p), colnames(x), rownames(p)))
betas["s", "v1", "s"] <- 3
betas["s", "v2", "s"] <- -2
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

# timesteps between years
t <- sample(5:10, n_sites, replace = T)
# t <- rep(1, n_sites)

for(s in 1:n_sites){
  
  xs <- x[s,] # local environment
  ps <- p # to, from
  
  # apply covariate effects
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
  
  probs <- sapply(1:k, function(i) project_simplex(ps, t[s], i))
  colnames(probs) <- colnames(ps)
  
  d <- tibble(s0 = sample(stages, n_i, replace = T),
              s1 =  apply(probs[,s0], 2, function(x) sample(stages, 1, prob = x))) %>%
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


# specify which transitions to infer
alpha_switch <- ceiling(p) %>% t()



## fit stan model ############

dat <- list(
  K = k,
  N = n_sites,
  D = ncol(x),
  y = y, 
  x = x, # site, var
  t = t,
  alpha_switch = alpha_switch, # from, to
  # beta_switch = beta_switch # from, var, to
  beta_switch = aperm(beta_switch, c(2, 1, 3)) # var, from, to
)

# compile
model <- cmdstan_model("models/v4/model.stan")

# fit <- model$optimize(
#   data = dat,
#   init = function() list(alpha = t(log(apply(log(p+1e-5), 2, softmax))))
# )
 
fit <- model$variational(
  data = dat,
  iter = 1e4
)

# # run MCMC
# fit <- model$sample(
#   data = dat,
#   # init = function() list(alpha = t(apply(log(p)+1e-10, 2, softmax))),
#   chains = 1,
#   iter_warmup = 200, iter_sampling = 100, refresh = 30
# )
# 
# fit




## assess fit #########

fs <- fit$summary %>% mutate(mean = estimate)

## alphas ##

pd <- p %>%
  as.data.frame() %>% 
  rownames_to_column("s1") %>% as_tibble() %>%
  gather(s0, true, -s1)

f <- fit$summary() %>%
  filter(str_detect(variable, "alpha_switched")) %>%
  mutate(variable = str_remove_all(variable, "alpha_switched|\\[|\\]")) %>%
  separate(variable, c("s0", "s1"), sep = ",") %>%
  mutate(s0 = colnames(p)[as.integer(s0)],
         s1 = rownames(p)[as.integer(s1)]) %>%
  left_join(pd) %>%
  group_by(s0) %>%
  mutate(pred = softmax(mean),
         q5 = softmax(q5),
         q95 = softmax(q95))


l <- 8
ggplot(f, aes(logit(true), logit(pred), 
              ymin = logit(q5), ymax = logit(q95))) + 
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  # geom_linerange(size = .25) +
  geom_point(aes(color = s1, shape = "t + 1"), size = 3) +
  geom_point(aes(color = s0, shape = "t"), size = 5) +
  scale_shape_manual(values = c(10, 16)) +
  coord_cartesian(xlim = c(-l, l),
                  ylim = c(-l, l)) +
  theme_minimal() +
  labs(color = "stage",
       shape = NULL)


l <- 8
ggplot(f %>% filter(true > 0, paste0(s0, s1) != "xx"), 
       aes(logit(true), logit(pred), 
           ymin = logit(q5), ymax = logit(q95))) + 
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_linerange(size = .25) +
  geom_point() +
  ggrepel::geom_text_repel(aes(label = paste0(s0, s1)), direction = "y") +
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
  separate(variable, c("v", "s0", "s1"), sep = ",") %>%
  mutate(s0 = colnames(p)[as.integer(s0)],
         s1 = rownames(p)[as.integer(s1)],
         v = colnames(x)[as.integer(v)]) %>%
  left_join(bd)

f %>%
  ggplot(aes(true, mean,
             # ymin = q5, ymax = q95
             )) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_point() +
  ggrepel::geom_text_repel(aes(label = paste0(v, "-", s0, s1)), direction = "y") +
  theme_minimal()
