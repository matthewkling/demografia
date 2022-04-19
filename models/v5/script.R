
# this version incorporates fecundity and unknown seedling turnover



library(tidyverse)
library(cmdstanr)
library(patchwork)


## utils #########

logit <- function(x) log(x / (1-x))
softmax <- function(x) exp(x) / sum(exp(x))
# softmax(log(c(.1, .2, .3, .4))) # demo that log is the inverse of softmax


## simulate dataset ########

stages <- c("p", "s", "j", "a", "x") # propagule-seedling-juvenile-adult-dead
k <- length(stages)

# true projection matrix -- origin in COLS, outcome in ROWS
p <- matrix(0, k, k, dimnames = list(stages, stages))
p["s", "p"] <- .05
p["s", "s"] <- .3
p["j", "s"] <- .01
p["j", "j"] <- .95
p["a", "j"] <- .01
p["a", "a"] <- .995
p["x",] <- apply(p, 2, function(x) 1 - sum(x))
p <- apply(p, 2, function(x) softmax(log(x)))

fecundity <- 400

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
n_i <- 1000 # number of individuals per block
n_blocks <- 100 # number of blocks (where a block is e.g. one resurvey for a site)

# environmental predictors
x <- matrix(runif(n_blocks*2), n_blocks, 2)
colnames(x) <- c("v1", "v2")

# specify predictor effects
betas <- array(0, 
               dim = c(ncol(p), ncol(x), nrow(p)),
               dimnames = list(colnames(p), colnames(x), rownames(p)))
# betas["s", "v1", "s"] <- 3
# betas["s", "v2", "s"] <- -2
betas["a", "v1", "x"] <- -1
betas["a", "v2", "a"] <- -1
betas["j", "v1", "j"] <- -.5
betas["j", "v2", "a"] <- 2

# binary indicator of which betas to infer (passed to model)
zero_beta <- betas
zero_beta[] <- as.integer(betas[] != 0)

# individual-level outcome counts
y <- array(dim = c(k, k, n_blocks),
           dimnames = list(colnames(p),
                           rownames(p),
                           paste0("b", 1:n_blocks)))

# population frequency counts
m <- array(dim = c(2, k, n_blocks),
           dimnames = list(c("t0", "tt"),
                           stages,
                           dimnames(y)[[3]]))

# timesteps between years
t <- sample(5:10, n_blocks, replace = T)
# t <- rep(1, n_blocks)

# simulate
for(s in 1:n_blocks){
  
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
  rownames(probs) <- colnames(probs) <- colnames(ps)
  
  d <- tibble(s0 = sample(stages, n_i, replace = T),
              s1 = apply(probs[, s0], 2, function(x) sample(stages, 1, prob = x))) %>%
    count(s0, s1) %>%
    full_join(expand_grid(s0 = colnames(p),
                          s1 = rownames(p))) %>%
    mutate(n = ifelse(is.na(n), 0, n)) %>%
    spread(s0, n) %>%
    column_to_rownames("s1") %>%
    as.matrix()
  d <- d[rownames(p), colnames(p)]
  y[ , , s] <- t(d)
  
  # full matrix projection
  A <- probs
  A["p", "a"] <- fecundity
  mi <- matrix(apply(d, 2, sum), 1,
               dimnames = list(NA, stages))
  mi[, "p"] <- mi[, "a"] * fecundity
  m["t0", , s] <- mi
  for(i in 1:t[s]) mi <- mi %*% t(A)
  m["tt", , s] <- mi
  
}

# # check that covarites had the desired effect
# plot(x[,1], y["a","x",])
# plot(x[,2], y["j","a",])


# specify which transitions to infer
zero_alpha <- ceiling(p) %>% t()

# for FIA, have to set ONE of these two params to zero (or another arbitrary value)
# because we don't have the individual-level data to infer both
zero_alpha["s", "x"] <- 0
zero_alpha["s", "s"] <- 0




## process to mimic FIA ########

## stage distribution changes ##

# cull to match FIA, with no propagule or dead counts
m[, c("p", "x"), ] <- 0

# add propagules, based on adults
# assumes adults at t = -1 == adults at t = 0, which we'll need to do for FIA
m[1, "p", ] <- m[1, "a", ] * fecundity

# specify which stages to evaluate in model likelihood calculation
m_eval <- 2:4

# convert to integers
# m <- round(m)


## individual-level transitions ##
y0_eval <- match(c("j", "a"), stages)
y1_eval <- match(c("j", "a", "x"), stages)
y <- y[y0_eval, y1_eval, ]




## fit stan model ############

dat <- list(
  K = k,
  N = n_blocks,
  D = ncol(x),
  yk0 = min(y0_eval),
  yK0 = length(y0_eval),
  yK1 = length(y1_eval),
  y = y, # from, to, block
  m = m,
  mk0 = min(m_eval),
  mk1 = max(m_eval),
  x = x, # block, var
  t = t,
  "F" = fecundity,
  Fk = match("a", stages),
  zero_alpha = zero_alpha, # from, to
  zero_beta = aperm(zero_beta, c(2, 1, 3)) # var, from, to
)

# compile
model <- cmdstan_model("models/v5/model.stan")

# fit <- model$optimize(
#   data = dat,
#   init = function() list(alpha = t(log(apply(log(p+1e-5), 2, softmax))))
# )

# fit <- model$variational(
#   data = dat
# )

# # run MCMC
fit <- model$sample(
  data = dat,
  # init = function() list(alpha = t(apply(log(p)+1e-10, 2, softmax))),
  chains = 1,
  iter_warmup = 200, iter_sampling = 100, refresh = 30
)



# fit$profiles()



## assess fit #########

# fs <- fit$summary() %>% mutate(mean = estimate)

## alphas ##

# f <- fit$draws()[,1,] %>%
#   as.data.frame() %>% as_tibble() %>%
#   mutate(i = 1:nrow(.)) %>%
#   gather(variable, value, -i) %>%
#   filter(str_detect(variable, "zeroed_alpha")) %>%
#   mutate(variable = str_remove_all(variable, "zeroed_alpha|\\[|\\]|1\\.")) %>%
#   separate(variable, c("s0", "s1"), sep = ",") %>%
#   mutate(s0 = colnames(p)[as.integer(s0)],
#          s1 = rownames(p)[as.integer(s1)]) %>%
#   group_by(i, s0) %>%
#   mutate(value = softmax(value)) %>%
#   ungroup() %>%
#   unite(trans, s0, s1) %>%
#   spread(trans, value)
# 
# fx <- f %>% select(starts_with("p_")) %>%
#   apply(1, softmax) %>%
#   t() %>%
#   as.data.frame() %>% as_tibble()
# 
# ggplot(f, aes(p_s, s_x, color = p_x)) +
#   geom_point() +
#   geom_smooth(method = lm, se = F, color = "black") +
#   # scale_y_log10() +
#   scale_color_viridis_c(trans = "log10")



pd <- p %>%
  as.data.frame() %>% 
  rownames_to_column("s1") %>% as_tibble() %>%
  gather(s0, true, -s1)

fa <- fit$summary() %>%
  filter(str_detect(variable, "zeroed_alpha")) %>%
  mutate(variable = str_remove_all(variable, "zeroed_alpha|\\[|\\]")) %>%
  separate(variable, c("s0", "s1"), sep = ",") %>%
  mutate(s0 = colnames(p)[as.integer(s0)],
         s1 = rownames(p)[as.integer(s1)]) %>%
  left_join(pd) %>%
  group_by(s0) %>%
  mutate(pred = softmax(mean),
         q5 = softmax(q5),
         q95 = softmax(q95))

# l <- 8
# ggplot(f, aes(logit(true), logit(pred), 
#               ymin = logit(q5), ymax = logit(q95))) + 
#   geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
#   # geom_linerange(size = .25) +
#   geom_point(aes(color = s1, shape = "t + 1"), size = 3) +
#   geom_point(aes(color = s0, shape = "t"), size = 5) +
#   scale_shape_manual(values = c(10, 16)) +
#   coord_cartesian(xlim = c(-l, l),
#                   ylim = c(-l, l)) +
#   theme_minimal() +
#   labs(color = "stage",
#        shape = NULL)

p_alpha <- ggplot(fa %>% filter(true > 0, paste0(s0, s1) != "xx"), 
                  aes(logit(true), logit(pred), 
                      ymin = logit(q5), ymax = logit(q95))) + 
  geom_abline(slope = 1, intercept = 0, alpha = .3) +
  geom_linerange(size = .25) +
  geom_point() +
  ggrepel::geom_text_repel(aes(label = paste0(s0, s1)), direction = "x") +
  theme_minimal() +
  labs(title = expression(alpha))



## betas ##

bd <- betas %>%
  as.data.frame() %>% 
  rownames_to_column("s0") %>% as_tibble() %>%
  gather(par, true, -s0) %>%
  separate(par, c("v", "s1"), sep = "\\.")

f <- fit$summary() %>%
  filter(str_detect(variable, "zeroed_beta")) %>%
  mutate(variable = str_remove_all(variable, "zeroed_beta|\\[|\\]")) %>%
  separate(variable, c("v", "s0", "s1"), sep = ",") %>%
  mutate(s0 = colnames(p)[as.integer(s0)],
         s1 = rownames(p)[as.integer(s1)],
         v = colnames(x)[as.integer(v)]) %>%
  left_join(bd)

p_beta <- f %>%
  ggplot(aes(true, mean,
             ymin = q5, ymax = q95
             )) +
  geom_abline(slope = 1, intercept = 0, alpha = .3) +
  geom_linerange(size = .25) +
  geom_point() +
  ggrepel::geom_text_repel(aes(label = paste0(v, "-", s0, s1)), direction = "x") +
  theme_minimal() +
  labs(title = expression(beta),
       y = "pred")


p <- p_alpha + p_beta +
  plot_layout(nrow = 1) &
  theme(plot.title = element_text(hjust = .5, size= 20))
ggsave("models/v5/alpha_beta.png", p, width = 6, height = 3.5, units = "in")




## simulate true vs predicted population trajectories #############

f <- fit$draws()[,1,] %>%
  as.data.frame() %>% as_tibble() %>%
  mutate(i = 1:nrow(.)) %>%
  gather(variable, value, -i) %>%
  filter(str_detect(variable, "zeroed_alpha")) %>%
  mutate(variable = str_remove_all(variable, "zeroed_alpha|\\[|\\]|1\\.")) %>%
  separate(variable, c("s0", "s1"), sep = ",") %>%
  mutate(s0 = colnames(p)[as.integer(s0)],
         s1 = rownames(p)[as.integer(s1)]) %>%
  group_by(i, s0) %>%
  mutate(value = softmax(value)) %>%
  ungroup()

sim <- function(x, p, fecundity, iter = 50){
  # x <- matrix(c(0, 0, 0, 10, 0), 1)
  # p[1,4] <- fecundity
  x <- as.vector(x)
  n <- matrix(NA, 1 + iter, ncol(p),
              dimnames = list(NULL, colnames(p)))
  n[1,] <- x
  for(i in 1:iter){
    # x <- x %*% t(p)
    xi <- sapply(1:length(x), function(i){
      rmultinom(1, x[i], p[,i])
    }) %>%
      apply(1, sum)
    xi[1] <- x[4] * fecundity
    # apply(p, 2, function(z) rmultinom(as.vector(x), 1, z))
    n[i+1,] <- xi
    x <- xi
  }
  
  n %>% as.data.frame() %>% as_tibble() %>%
    mutate(t = 0:iter)
}


nt <- 25
# n0 <- matrix(c(10*fecundity, 300, 30, 10, 0), 1)
n0 <- matrix(c(10*fecundity, 1000, 1000, 1000, 0), 1)
s <- map_df(unique(f$i), function(i) sim(n0, p, fecundity, nt) %>%
              mutate(i = i)) %>%
  mutate(model = "true")
# s <- sim(n0, p, fecundity, nt) %>%
#   mutate(i = 0)

x <- f %>%
  split(.$i) %>%
  map_df(function(z){
    y <- z %>%
      select(-i) %>%
      spread(s0, value) %>%
      column_to_rownames("s1")
    y <- y[match(rownames(p), rownames(y)), 
           match(colnames(p), colnames(y))]
    sim(n0, y, fecundity, nt) %>%
      mutate(i = pull(z, i)[1])
  }) %>%
  mutate(model = "pred") %>%
  bind_rows(s) %>%
  gather(stage, n, -t, -i, -model) %>%
  mutate(stage = factor(stage, levels = colnames(p)))

# x %>%
#   filter(stage != "x") %>%
#   spread(model, n) %>%
#   group_by(t, stage) %>%
#   summarize(p = mean(sign(pred - true))) %>%
#   ggplot(aes(t, p, color = stage)) +
#   geom_line() +
#   ylim(c(-1, 1)) +
#   theme_minimal() +
#   labs(y = "mean(sign(pred - true))")
# 
# x %>%
#   filter(stage != "x") %>%
#   ggplot(aes(t, n, color = model, group = paste(model, i))) +
#   facet_wrap(~stage, scales = "free") +
#   geom_line(alpha = .25) +
#   theme_bw()

x %>%
  filter(stage != "x") %>%
  group_by(stage, model, t) %>%
  mutate(mid = median(n),
         low = quantile(n, .05), 
         high = quantile(n, .95)) %>%
  ggplot(aes(t, mid, ymin = low, ymax = high, 
             color = model, fill = model)) +
  facet_wrap(~stage, scales = "free", nrow = 1) +
  geom_ribbon(alpha = .25, color = NA) +
  geom_line() +
  theme_bw() +
  labs(y = "n")
