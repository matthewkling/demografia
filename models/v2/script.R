
# this version adds adds sparse temporal sampling.


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
stages1 <- rownames(p)

# fecundity <- 400



# simulate
n_i <- 1000 # number of individuals
n_steps <- 10

d0 <- sample(stages, n_i, replace = T)
di <- d0
d <- d0
for(i in 1:n_steps){
  alive <- di != "x"
  di[alive] <- apply(p[,di[alive]], 2, function(x) sample(stages1, 1, prob = x))
  d <- cbind(d, di)
}



d <- tibble(s0 = d[,1],
            s1 = d[,ncol(d)]) %>%
  count(s0, s1) %>%
  full_join(expand_grid(s0 = colnames(p),
                        s1 = rownames(p))) %>%
  mutate(n = ifelse(is.na(n), 0, n)) %>%
  spread(s0, n) %>%
  column_to_rownames("s1") %>%
  as.matrix()
d <- d[rownames(p), colnames(p)]

# sum(d[,"a"]) * fecundity




# hyperparameters for transition probability priors
alpha <- c(s = c(0, 1, 0, 0, 10),
           b = c(0, 2, 1, 0, 2),
           j = c(0, 0, 10, 1, 1),
           a = c(0, 0, 0, 10, 1)) %>%
  matrix(ncol = 4, byrow = F)
alpha[alpha == 0] <- .01
rownames(alpha) <- rownames(p)
colnames(alpha) <- colnames(p)

# alpha[] <- 1 # just use a flat prior


## fit stan model ############

# data inputs
dat <- list(
  K = ncol(d),
  t = n_steps,
  y = t(d),
  alpha = t(alpha)
)

# compile
model <- cmdstan_model("models/v2/model.stan")

# run MCMC
fit <- model$sample(
  data = dat,
  chains = 2
)



## assess fit #########

pd <- p %>%
  as.data.frame() %>% as_tibble() %>%
  rownames_to_column("s1") %>%
  mutate(s1 = rownames(p)[as.integer(s1)]) %>%
  gather(s0, true, -s1)

ad <- alpha %>%
  as.data.frame() %>% as_tibble() %>%
  rownames_to_column("s1") %>%
  mutate(s1 = rownames(p)[as.integer(s1)]) %>%
  gather(s0, alpha, -s1)

f <- fit$summary() %>%
  filter(str_detect(variable, "theta")) %>%
  mutate(variable = str_remove_all(variable, "theta|\\[|\\]")) %>%
  separate(variable, c("s0", "s1"), sep = ",") %>%
  mutate(s0 = colnames(p)[as.integer(s0)],
         s1 = rownames(p)[as.integer(s1)]) %>%
  left_join(pd) %>%
  left_join(ad)


logit <- function(p) log(p / (1-p))
l <- 7
ggplot(f, aes(logit(true), logit(median), 
              ymin = logit(q5), ymax = logit(q95),
              color = s0)) + 
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_pointrange(size = .25) +
  coord_cartesian(xlim = c(-l, l),
                  ylim = c(-l, l)) +
  theme_minimal()

filter(f, alpha != min(alpha), true == 0)


p
round(f, 3)
alpha
