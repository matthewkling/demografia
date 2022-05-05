
# this version tackles density dependence
# it also implements proper elliptic paraboloids for betas

# environmental effects: mortality
# density effects: transition
# reference class: survival




## setup #################

library(tidyverse)
library(cmdstanr)
library(patchwork)
library(popbio)

here <- "models/v7"

## utils ####################

logit <- function(x) log(x / (1-x))

softmax <- function(x) exp(x) / sum(exp(x))
# softmax(log(c(.1, .2, .3, .4))) # demo that log is the inverse of softmax

project_multinom <- function(x, p){
   sapply(1:length(x), function(j){
      xj <- x[j]
      if(xj <= .Machine$integer.max) return(rmultinom(1, xj, p[,j]))
      if(xj > .Machine$integer.max) return(round(xj * p[,j]))
   }) %>%
      apply(1, sum)
}

sim <- function(x, p, fecundity, iter = 50){
   # x <- matrix(c(0, 0, 0, 10, 0), 1)
   # p[1,4] <- fecundity
   x <- as.vector(x)
   k <- length(x)
   n <- matrix(NA, 1 + iter, ncol(p),
               dimnames = list(NULL, colnames(p)))
   n[1,] <- x
   for(i in 1:iter){
      # x <- x %*% t(p)
      xi <- project_multinom(x, p)
      xi[1] <- x[4] * fecundity
      n[i+1,] <- xi[1:k]
      x <- xi[1:k]
   }
   
   n %>% as.data.frame() %>% as_tibble() %>%
      mutate(t = 0:iter)
}

project_simplex <- function(p, # projection matrix 
                            t, # number of time steps
                            i){ # life stage to project
   n <- matrix(0, 1, nrow(p))
   n[i] <- 1
   for(j in 1:t) n <- n %*% t(p)
   return(as.vector(n))
}








## priors #############

# see alpha_priors.r for testing and vis

make_priors <- function(){
   
   stages <- c("p", "s", "j", "a") # propagule-seedling-juvenile-adult-dead
   labels <- c("propagules", "seedlings", "saplings", "adults")
   names(labels) <- stages
   k <- length(stages)
   
   
   ## alphas: transition probabilities under mean environmental conditions ##
   # origin in COLS, outcome in ROWS
   
   p <- matrix(0, k+1, k, dimnames = list(c(stages, "x"), stages))
   p["s", "p"] <- .005
   p["s", "s"] <- .3
   p["j", "s"] <- .05
   p["j", "j"] <- .9
   p["a", "j"] <- .05
   p["a", "a"] <- .995
   p["x",] <- apply(p, 2, function(x) 1 - sum(x))
   p <- apply(p, 2, function(x) softmax(log(x)))
   
   # boolean indicating nonzero transitions
   alpha_i <- ceiling(t(p)) 
   
   # reference classes for each source stage (will have alphas fixed, not random)
   alpha_ref <- alpha_i
   alpha_ref[] <- 0
   alpha_ref["p", "s"] <- 1
   alpha_ref["s", "s"] <- 1
   alpha_ref["j", "j"] <- 1
   alpha_ref["a", "a"] <- 1
   alpha_rand <- alpha_i - alpha_ref
   
   mu <- log(p)
   sigma <- mu
   sigma[] <- rep(c(2, 1, 1, .75), each = nrow(p))
   
   alpha_mu <- mu[t(alpha_rand) == 1]
   alpha_sigma <- sigma[t(alpha_rand) == 1]
   
   
   ## betas: effect of environmental variables on log prob ##
   
   n_pred <- 2 # number of predictors
   pred_names <- paste0("v", 1:n_pred)
   
   beta <- array(0, 
                 dim = c(n_pred, ncol(p), nrow(p)),
                 dimnames = list(pred_names, colnames(p), rownames(p)))
   
   # flag betas to be nonzero
   # beta[, "a", "a"] <- 1
   # beta[, "j", "j"] <- 1
   # beta[, "j", "a"] <- 1
   # beta[, "s", "j"] <- 1
   # beta[, "s", "x"] <- 1
   beta[, "j", "x"] <- 1
   beta[, "a", "x"] <- 1
   
   # prior standard deviation of beta location terms
   beta_loc_sd <- 2
   
   # prior exponential scale parameter for beta curvature terms
   beta_shape <- 3
   beta_rate <- 10
   
   
   ## gammas: density dependence ##
   
   gamma <- array(0, 
                  dim = c(ncol(p), ncol(p), nrow(p)),
                  dimnames = list(colnames(p), colnames(p), rownames(p)))
   
   # flag gammas to be nonzero
   gamma["s", "s", "j"] <- 1
   gamma["j", "s", "j"] <- 1
   gamma["j", "j", "a"] <- 1
   gamma["a", "j", "a"] <- 1
   # gamma["a", "a", "a"] <- 1
   
   # exponential scale parameter for gammas
   gamma_scale <- 1000
   
   
   
   ## bundle ##
   list(p = p,
        mu = mu,
        sigma = sigma,
        alpha_i = alpha_i,
        alpha_rand = alpha_rand,
        alpha_ref = alpha_ref,
        alpha_mu = alpha_mu,
        alpha_sigma = alpha_sigma,
        beta = beta,
        beta_loc_sd = beta_loc_sd,
        beta_shape = beta_shape,
        beta_rate = beta_rate,
        gamma = gamma,
        gamma_scale = gamma_scale,
        stages = labels)
}




## species ###########

# create true species parameters by sampling from prior
make_species <- function(prior){
   
   fecundity <- 500
   
   # alpha
   p <- prior$p
   for(i in 1:ncol(p)){
      nz <- which(p[,i] > 0)
      p[nz, i] <- softmax(rnorm(length(nz), prior$mu[nz, i], prior$sigma[nz, i]))
   }
   
   # betas location and curvature
   b_crv <- b_loc <- prior$beta
   b_loc[b_loc == 1] <- rnorm(sum(b_loc), 0, prior$beta_loc_sd)
   # b_shp[b_shp == 1] <- rexp(sum(b_shp), prior$beta_shp_scale)
   b_crv[b_crv == 1] <- rgamma(sum(b_crv), prior$beta_shape, prior$beta_rate)
   
   # gamma
   g <- prior$gamma
   g[g == 1] <- 0 - rexp(sum(g), prior$gamma_scale)
   
   # bundle
   list(p = p,
        fecundity = fecundity,
        beta_crv = b_crv,
        beta_loc = b_loc,
        gamma = g)
}

# visualize environmental effects
plot_betas <- function(sp, normalize = T){
   
   p <- log(sp$p)
   bl <- sp$beta_loc
   bc <- sp$beta_crv
   
   p <- p %>%
      as.data.frame() %>%
      rownames_to_column("s1") %>%
      gather(s0, lp, -s1) %>%
      filter(is.finite(lp)) %>%
      expand_grid(v1 = seq(-3, 3, .2),
                  v2 = seq(-3, 3, .2))
   
   bl <- bl %>%
      as.data.frame() %>%
      rownames_to_column("v") %>%
      gather(trans, b, -v) %>%
      separate(trans, c("s0", "s1")) %>%
      mutate(v = paste0(v, "_bl")) %>%
      spread(v, b)
   
   bc <- bc %>%
      as.data.frame() %>%
      rownames_to_column("v") %>%
      gather(trans, b, -v) %>%
      separate(trans, c("s0", "s1")) %>%
      mutate(v = paste0(v, "_bc")) %>%
      spread(v, b)
   
   d <- p %>%
      left_join(bl) %>%
      left_join(bc) %>%
      mutate(s0 = factor(s0, levels = colnames(sp$p)),
             s1 = factor(s1, levels = rownames(sp$p))) %>%
      mutate(lp = lp + (v1 - v1_bl)^2 * v1_bc^2 + (v2 - v2_bl)^2 * v2_bc^2)
   
   if(normalize){
      d <- d %>%
         group_by(v1, v2, s0) %>%
         mutate(p = softmax(lp))
      
      ggplot(d, aes(v1, v2, fill = logit(p), z = logit(p))) +
         facet_grid(s1 ~ s0) +
         geom_raster() +
         geom_contour(color = "white", alpha = .5,
                      breaks = seq(min(logit(d$p)), max(logit(d$p)), length.out = 50)) +
         scale_fill_viridis_c() +
         theme_bw() +
         theme(panel.grid = element_blank()) +
         scale_x_continuous(expand = c(0, 0)) +
         scale_y_continuous(expand = c(0, 0)) 
      
   }else{
      ggplot(d, 
             aes(v1, v2, fill = lp, z = lp)) +
         facet_grid(s1 ~ s0) +
         geom_raster() +
         geom_contour(color = "white", alpha = .5,
                      breaks = seq(min(d$lp), max(d$lp), length.out = 50)) +
         scale_fill_viridis_c() +
         theme_bw() +
         theme(panel.grid = element_blank()) +
         scale_x_continuous(expand = c(0, 0)) +
         scale_y_continuous(expand = c(0, 0)) 
   }
   
}


# visualize density effects as surfaces
plot_gammas <- function(sp){
   
   p <- log(sp$p) %>%
      as.data.frame() %>%
      rownames_to_column("s1") %>%
      gather(s0, lp, -s1) %>%
      filter(is.finite(lp))
   
   g <- sp$gamma %>%
      as.data.frame() %>%
      rownames_to_column("se") %>%
      gather(trans, g, -se) %>%
      separate(trans, c("s0", "s1")) %>%
      filter(g != 0) %>%
      left_join(p) %>%
      expand_grid(n = 10 ^ seq(log10(1), log10(10000), length.out = 100)) %>%
      mutate(p = exp(lp + g * n),
             trans = paste0(s0, s1),
             effect = paste0(se, "-", s0, s1))
   
   g %>%
      ggplot(aes(n, p, color = se)) +
      facet_wrap(~trans, scales = "free") +
      geom_line() +
      ylim(0, NA) +
      theme_bw() +
      labs(x = "N",
           y = "transition probability",
           color = "density\nclass")
   
   # tibble(n = 0:10000, p = 1/n) %>% ggplot(aes(n, p)) + geom_line() + scale_y_log10()
}





# visualize density effects as time series
plot_simulated_ts <- function(sp, n0 = matrix(c(1000, 100, 10, 1), 1), t = 1000, type = "lines"){
   
   p <- log(sp$p)
   g <- sp$gamma
   bl <- sp$beta_loc
   bc <- sp$beta_crv
   
   d <- expand_grid(v1 = -3:3, v2 = -3:3) %>%
      pmap_df(function(v1, v2){
         env <- c(v1, v2)
         for(d in 1:length(env)) p <- p + t((env[d] - bl[d,,])^2 * bc[d,,]^2)
         
         n <- matrix(NA, t+1, length(n0))
         nt <- n0
         n[1,] <- nt
         for(i in 1:t){
            pt <- p
            for(k in 1:dim(g)[1]) pt <- pt + nt[k] * t(g[k,,])
            for(k in 1:ncol(pt)) pt[,k] <- softmax(pt[,k])
            pt[1,4] <- sp$fecundity
            nt <- nt[1,1:4] %*% t(pt)
            n[i+1,] <- nt[1,1:4]
         }
         
         n %>%
            as.data.frame() %>%
            setNames(c("seed", "seedling", "sapling", "adult")) %>%
            as_tibble() %>%
            mutate(t = 1:nrow(.))  %>%
            gather(stage, N, -t) %>%
            mutate(stage = factor(stage, levels = c("seed", "seedling", "sapling", "adult")),
                   v1 = v1,
                   v2 = v2)
      })
   
   
   if(type == "lines"){
      plt <- d %>%
         ggplot(aes(t, N, color = v1, alpha = v2, group = paste(v1, v2))) +
         facet_wrap(~stage, scales = "free", nrow = 1) +
         geom_line() +
         scale_color_gradientn(colors = c("black", "blue", "red", "orange")) +
         scale_alpha_continuous(range = c(.3, 1)) +
         theme_bw() +
         ylim(0, NA)
   }
   
   if(type == "heatmap"){
      d <- d %>% 
         ungroup() %>%
         filter(t == max(t)) %>%
         group_by(stage) %>%
         mutate(NN = N / max(N))
      l <- d %>%
         group_by(stage) %>%
         filter(N == max(N))
      plt <- d %>%
         ggplot(aes(v1, v2, fill = NN)) +
         facet_wrap(~stage, scales = "free", nrow = 1) +
         geom_raster() +
         geom_point(data = l) +
         ggrepel::geom_text_repel(data = l, aes(label = round(N))) +
         scale_fill_viridis_c() +
         theme_bw() +
         scale_x_continuous(expand = c(0, 0)) +
         scale_y_continuous(expand = c(0, 0)) +
         labs(fill = paste0("N / max(N)\n(t = ", t, ")"))
      
   }
   
   return(plt)
}

# plot_simulated_ts(sp, n0 = matrix(rep(100, 4), 1), type = "lines")
# plot_simulated_ts(sp, n0 = matrix(rep(100, 4), 1), type = "heatmap")





## simulate dataset ############################

simulate_data <- function(sp, n_blocks = 100, t = 5, n0 = 10000){
   
   io_stages <- 3:4 # stages for which individual outcomes are tracked
   
   p <- sp$p
   f <- sp$fecundity
   bl <- sp$beta_loc
   bc <- sp$beta_crv
   g <- sp$gamma
   
   # predictor data
   n_pred <- dim(bl)[1]
   x <- matrix(rnorm(n_blocks * n_pred), 
               n_blocks, n_pred, 
               dimnames = list(paste0("b", 1:n_blocks),
                               dimnames(bl)[[1]]))
   
   # containers
   n <- array(dim = c(n_blocks, ncol(p), 2),
              dimnames = list(rownames(x),
                              colnames(p),
                              c("t0", "tt")))
   io <- array(dim = c(n_blocks, length(io_stages), nrow(p)),
               dimnames = list(rownames(x),
                               colnames(p)[io_stages],
                               rownames(p)))
   
   # simulate for each block
   prog = txtProgressBar(min = 0, max = n_blocks, initial = 0) 
   for(i in 1:n_blocks){
      setTxtProgressBar(prog,i)
      
      # assemble projection matrix from alpha, beta, and x
      pi <- log(p)
      for(s0 in 1:ncol(p)){
         for(s1 in 1:nrow(p)){
            for(v in 1:ncol(x)){
               # pi[s1, s0] <- pi[s1, s0] + b[v, s0, s1] * x[i, v] + b2[v, s0, s1] * x[i, v]^2
               pi[s1, s0] <- pi[s1, s0] + (x[i, v] - bl[v, s0, s1])^2 * bc[v, s0, s1]^2
            }
         }
         pi[, s0] = softmax(pi[, s0])
      }
      
      # starting population values from stable stage distribution + noise
      A <- cbind(pi, x = c(rep(0, 4), 1))
      A[1, 4] <- f
      eig <- eigen.analysis(A[1:4, 1:4])
      ssd <- eig$stable.stage[1:4]
      ssd <- ssd / sum(ssd)
      n_t0 <- ceiling(softmax(log(ssd) + rnorm(length(ssd), 0, 1)) * n0) + 1
      n_t0 <- c(1000, 100, 10, 10) * 10
      n_t0 <- ceiling(n_t0 * 10 ^ runif(4, -1, 1))
      n[i, , 1] <- n_t0
      
      
      nt <- matrix(n_t0, 1)
      M <- matrix(0, 2, 5)
      for(k in 1:length(io_stages)) M[k,io_stages[k]] <- 1
      
      for(ti in 1:t){
         
         # apply density effects
         pit <- log(pi)
         for(s in 1:ncol(p)) pit <- pit + t(g[s,,]) * nt[s]
         for(s0 in 1:ncol(p)) pit[, s0] <- softmax(pit[, s0])
         
         # project population
         propagules <- f * nt[4]
         nt <- matrix(project_multinom(nt[1,], pit), 1)
         nt <- matrix(nt[1, 1:4], 1)
         nt[1, 1] <- propagules
         
         # individual outcomes
         pit <- cbind(pit, c(rep(0, 4), 1))
         pit[1, 4] <- 0
         for(k in 1:length(io_stages)) M[k,] <- M[k,] %*% t(pit)
         
      }
      
      # collective counts
      n[i, , 2] <- nt
      
      # individual outcomes
      io[i, , ] <- io_stages %>%
         sapply(function(j) rmultinom(1, n_t0[j], M[match(j, io_stages),])) %>%
         t()
   }
   close(prog)
   
   # bundle
   list(n = n,
        io = io,
        x = x,
        io_stages = io_stages,
        t = t,
        K = ncol(sp$p),
        fecundity = sp$fecundity)
}


plot_individual_outcomes <- function(d){
   
   x <- d$x %>%
      as.data.frame() %>%
      rownames_to_column("block") %>%
      as_tibble()
   
   f <- d$io %>%
      as.data.frame() %>%
      rownames_to_column("block") %>%
      as_tibble() %>%
      gather(trans, value, -block) %>%
      separate(trans, c("s0", "s1")) %>%
      group_by(block, s0) %>%
      mutate(p = value / sum(value)) %>%
      left_join(x) %>%
      mutate(s0 = factor(s0, levels = c("p", "s", "j", "a")),
             s1 = factor(s1, levels = c("p", "s", "j", "a", "x")))
   
   b <- f %>%
      mutate(v1 = round(v1, 0),
             v2 = round(v2, 0)) %>%
      group_by(s0, s1, v1, v2) %>%
      summarize(value = sum(value)) %>%
      group_by(v1, v2, s0) %>%
      mutate(p = value / sum(value))
   
   ggplot(b, aes(v1, v2, fill = logit(p))) +
      facet_grid(s1 ~ s0) +
      geom_raster() +
      scale_fill_viridis_c(na.value = "black") +
      theme_bw() +
      theme(panel.grid = element_blank()) +
      scale_x_continuous(expand = c(0, 0)) +
      scale_y_continuous(expand = c(0, 0))
}


# 
# 
# prior <- make_priors()
# sp <- make_species(prior)
# sp$p[,3] <- c(0, 0, .9, .1, .1)
# sp$p[,4] <- c(0, 0, 0, .99, .01)
# sp$gamma["a", "j", "a"] <- -.01
# sp$gamma["a", "a", "a"] <- -.01 
# d <- simulate_data(sp, n_blocks = 1000, n0 = 1e5, t = 1)
# 
# 
# tibble(n0 = d$n[,4,1],
#            nt = d$n[,4,2],
#            dn = nt - n0) %>%
#    ggplot(aes(n0, dn)) +
#    geom_jitter()





## fit stan model ###################

fit_model <- function(d, prior, algorithm = "sample"){
   
   # prepare input data
   dat <- list(
      K = d$K,
      N = nrow(d$x),
      D = ncol(d$x),
      n0 = d$n[ , , 1],
      nt = d$n[ , , 2],
      nk = 2,
      io = d$io,
      iok = d$io_stages,
      T = d$t,
      "F" = d$fecundity,
      Fk = match("a", dimnames(d$n)[[2]]),
      x = d$x,
      
      # alpha_i = prior$alpha_i,
      alpha_rand = prior$alpha_rand,
      alpha_ref = prior$alpha_ref,
      n_alpha = sum(prior$alpha_rand),
      alpha_mu = prior$alpha_mu,
      alpha_sigma = prior$alpha_sigma,
      
      n_beta = sum(prior$beta),
      beta_i = prior$beta,
      beta_loc_sd = prior$beta_loc_sd,
      beta_shape = prior$beta_shape,
      beta_rate = prior$beta_rate,
      
      n_gamma = sum(prior$gamma),
      gamma_i = prior$gamma,
      gamma_scale = prior$gamma_scale
   )
   
   model <- cmdstan_model(paste0(here, "/model.stan"))
   
   if(algorithm == "sample"){
      fit <- model$sample(
         data = dat,
         chains = 1,
         # adapt_delta = .99
         iter_warmup = 200, iter_sampling = 100, refresh = 30
      )
   }
   
   if(algorithm == "variational"){
      fit <- model$variational(data = dat)
   }
   
   if(algorithm == "optimize"){
      fit <- model$optimize(data = dat)
   }
   
   return(fit)
}





## assess fit ##############################

## true vs fitted parameters ##

plot_param_fits <- function(fit, sp, 
                           save = F, tag = NULL){
   
   ## alphas ##
   
   pd <- sp$p %>% 
      as.data.frame() %>% 
      rownames_to_column("s1") %>% as_tibble() %>%
      gather(s0, true, -s1)
   
   pd <- fit$draws() %>%
      as.data.frame() %>% as_tibble() %>%
      mutate(i = 1:nrow(.)) %>%
      gather(variable, value, -i) %>%
      filter(str_detect(variable, "A")) %>%
      mutate(variable = str_remove_all(variable, "A|\\[|\\]")) %>%
      separate(variable, c("s0", "s1"), sep = ",") %>%
      mutate(s0 = colnames(sp$p)[as.integer(s0)],
             s1 = rownames(sp$p)[as.integer(s1)]) %>%
      group_by(s0, i) %>%
      mutate(prob = softmax(value)) %>%
      group_by(s0, s1) %>%
      summarize(pred = mean(prob),
                q5 = quantile(prob, .05),
                q95 = quantile(prob, .95)) %>%
      ungroup() %>%
      mutate(trans = paste0(s0, s1),
             type = ifelse(str_detect(trans, "j|a"), "individual", "collective")) %>% 
      left_join(pd) %>%
      filter(true > 0, 
             trans != "xx")
   
   alpha <- pd %>%
      ggplot(aes(logit(true), logit(pred), ymin = logit(q5), ymax = logit(q95),
                 color = type)) + 
      geom_abline(slope = 1, intercept = 0, alpha = .3) +
      ggrepel::geom_text_repel(aes(label = trans), 
                               direction = "x", segment.size = .25) +
      geom_linerange(size = .5) +
      geom_point() +
      scale_color_manual(values = c("darkred", "darkblue")) +
      theme_minimal() +
      theme(legend.position = "bottom") +
      labs(title = expression(alpha),
           x = "true logit transition probability",
           y = "posterior logit transition probability\n(median & 95% CI)",
           color = "data type")
   
   ## betas ##
   
   # location term #
   
   bl <- sp$beta_loc %>%
      as.data.frame() %>% 
      rownames_to_column("v") %>% as_tibble() %>%
      gather(par, true, -v) %>%
      separate(par, c("s0", "s1"), sep = "\\.")
   
   bl <- fit$summary() %>%
      filter(str_detect(variable, "BL")) %>%
      mutate(variable = str_remove_all(variable, "BL|\\[|\\]")) %>%
      separate(variable, c("v", "s0", "s1"), sep = ",") %>%
      mutate(s0 = colnames(sp$p)[as.integer(s0)],
             s1 = rownames(sp$p)[as.integer(s1)],
             v = dimnames(sp$beta_loc)[[1]][as.integer(v)]) %>%
      left_join(bl) %>%
      filter(true != 0) %>%
      mutate(label = paste0(v, "-", s0, s1))
   
   beta_loc <- bl %>%
      ggplot(aes(true, median, ymin = q5, ymax = q95)) +
      geom_abline(slope = 1, intercept = 0, alpha = .3) +
      ggrepel::geom_text_repel(aes(label = label), direction = "x", segment.size = .25) +
      geom_linerange(size = .5, color = "darkgreen") +
      geom_point(color = "darkgreen") +
      theme_minimal() +
      labs(title = expression(beta[optimum]),
           x = "true effect",
           y = "posterior effect\n(median & 95% CI)")
   
   
   # curvature term #
   
   bc <- sp$beta_crv %>%
      as.data.frame() %>% 
      rownames_to_column("v") %>% as_tibble() %>%
      gather(par, true, -v) %>%
      separate(par, c("s0", "s1"), sep = "\\.")
   
   bc <- fit$summary() %>%
      filter(str_detect(variable, "BC")) %>%
      mutate(variable = str_remove_all(variable, "BC|\\[|\\]")) %>%
      separate(variable, c("v", "s0", "s1"), sep = ",") %>%
      mutate(s0 = colnames(sp$p)[as.integer(s0)],
             s1 = rownames(sp$p)[as.integer(s1)],
             v = dimnames(sp$beta_crv)[[1]][as.integer(v)]) %>%
      left_join(bc) %>%
      filter(true != 0) %>%
      mutate(label = paste0(v, "-", s0, s1))
   
   beta_crv <- bc %>%
      ggplot(aes(true, median, ymin = q5, ymax = q95)) +
      geom_vline(xintercept = 0, alpha = .3) +
      geom_hline(yintercept = 0, alpha = .3) +
      geom_abline(slope = 1, intercept = 0, alpha = .3) +
      ggrepel::geom_text_repel(aes(label = label), direction = "x", segment.size = .25) +
      geom_linerange(size = .5, color = "darkgreen") +
      geom_point(color = "darkgreen") +
      theme_minimal() +
      labs(title = expression(beta[slope]),
           x = "true effect",
           y = "posterior effect\n(median & 95% CI)")
   
   
   
   ## gammas ##
   
   gd <- sp$gamma %>%
      as.data.frame() %>%
      rownames_to_column("sx") %>% as_tibble() %>%
      gather(par, true, -sx) %>%
      separate(par, c("s0", "s1"), sep = "\\.")
   
   gd <- fit$summary() %>%
      filter(str_detect(variable, "G")) %>%
      mutate(variable = str_remove_all(variable, "G|\\[|\\]")) %>%
      separate(variable, c("sx", "s0", "s1"), sep = ",") %>%
      mutate(s0 = colnames(sp$p)[as.integer(s0)],
             s1 = rownames(sp$p)[as.integer(s1)],
             sx = colnames(sp$p)[as.integer(sx)]) %>%
      left_join(gd) %>%
      filter(true != 0)
   
   gamma <- gd %>%
      ggplot(aes(true, median, ymin = q5, ymax = q95)) +
      geom_vline(xintercept = 0, alpha = .3) +
      geom_hline(yintercept = 0, alpha = .3) +
      geom_abline(slope = 1, intercept = 0, alpha = .3) +
      ggrepel::geom_text_repel(aes(label = paste0(sx, "-", s0, s1)),
                               direction = "x", segment.size = .25) +
      geom_linerange(size = .5, color = "darkorchid") +
      geom_point(color = "darkorchid") +
      xlim(NA, 0) +
      ylim(NA, 0) +
      theme_minimal() +
      labs(title = expression(gamma),
           x = "true effect",
           y = "posterior effect\n(median & 95% CI)")
   
   
   ## assemble ##
   
   p_abg <- alpha + gamma + beta_crv + beta_loc + 
      plot_layout(nrow = 2) &
      theme(plot.title = element_text(hjust = .5, size = 30))
   
   if(save) ggsave(paste0(here, "/ABG", tag, ".pdf"), p_abg, width = 8, height = 4, units = "in")
   
   p_abg
}




## simulated projections under true vs fitted params ##

plot_trajectories <- function(fit, sp, prior, d, 
                            t = 100, reps = 100, 
                            seed = NULL, n0 = NULL, 
                            save = F, tag = NULL){
   
   # starting population for simulation
   set.seed(seed)
   n0sim <- d$n[sample(1:(dim(d$n)[1]), 1), , 1]
   # n0sim[1] <- n0sim[1] * d$fecundity / sp$fecundity
   set.seed(NULL)
   if(!is.null(n0)) n0sim <- n0
   
   # format posterior values
   f <- fit$draws() %>% #[,1,] %>%
      as.data.frame() %>% as_tibble() %>%
      mutate(i = 1:nrow(.)) %>%
      gather(variable, value, -i) %>%
      filter(str_detect(variable, "A")) %>%
      mutate(variable = str_remove_all(variable, "A|\\[|\\]|1\\.")) %>%
      separate(variable, c("s0", "s1"), sep = ",") %>%
      mutate(s0 = colnames(sp$p)[as.integer(s0)],
             s1 = rownames(sp$p)[as.integer(s1)]) %>%
      ungroup() %>%
      mutate(value = ifelse(s0 == "a" & s1 == "p", 0, value)) %>%
      filter(i %in% sample(max(i), min(max(i), reps)))
   
   # posterior sim
   post <- f %>%
      split(.$i) %>%
      map_df(function(z){
         y <- z %>%
            select(-i) %>%
            spread(s0, value) %>%
            column_to_rownames("s1")
         y <- y[match(rownames(sp$p), rownames(y)), 
                match(colnames(sp$p), colnames(y))]
         yy <- try(sim(n0sim, y, d$fecundity, t) %>% # use d$fecundity rather than sp$fecundity in case they differ
                      mutate(i = pull(z, i)[1]))
         if(class(yy) == "try-error") return(NULL)
         yy
      }) %>%
      mutate(model = "posterior") 
   
   # prior sim
   pryr <- map_df(unique(f$i), 
                  function(i){
                     yy <- try(sim(n0sim, make_species(prior)$p, sp$fecundity, t) %>%
                                  mutate(i = i))
                     if(class(yy) == "try-error") return(NULL)
                     yy
                  }) %>%
      mutate(model = "prior")
   
   # true sim
   true <- map_df(unique(f$i), 
                  function(i) sim(n0sim, sp$p, sp$fecundity, t) %>%
                     mutate(i = i)) %>%
      mutate(model = "true")
   
   
   
   # combine and summarize
   x <- bind_rows(true, post, pryr) %>%
      gather(stage, n, -t, -i, -model) %>%
      mutate(stage = factor(stage, levels = names(prior$stages), labels = prior$stages),
             model = factor(model, levels = c("prior", "posterior", "true"))) %>%
      filter(stage != "x") %>%
      group_by(stage, model, t) %>%
      mutate(mid = median(n),
             low = quantile(n, .025), 
             high = quantile(n, .975))
   
   # plot
   p_traj <- x %>%
      ggplot(aes(t, mid, ymin = low, ymax = high, 
                 color = model, fill = model)) +
      facet_wrap(~stage, scales = "free_x", nrow = 1) +
      # facet_grid(stage ~ ., scales = "free") +
      geom_ribbon(alpha = .25, color = NA) +
      geom_line() +
      scale_y_log10(breaks = c(1, 10, 100, 1000, 10000, 100000),
                    minor_breaks = c(1:9,
                                     seq(10, 90, 10),
                                     seq(100, 900, 100),
                                     seq(1000, 9000, 1000),
                                     seq(10000, 90000, 10000),
                                     seq(100000, 1000000, 100000))) +
      scale_x_continuous(expand = c(0, 0)) +
      coord_cartesian(ylim = c(NA, max(x$high[x$model != "prior"]))) +
      scale_color_manual(values = c("gray40", "darkred", "dodgerblue")) +
      scale_fill_manual(values = c("gray40", "darkred", "dodgerblue")) +
      theme_bw() +
      theme(strip.background = element_rect(fill = "black"),
            strip.text = element_text(color = "white", size = 16),
            legend.position = "top") +
      labs(y = paste0("N(t) across ", reps, " simulations\n(median & 95% CI)"),
           x = "t",
           color = NULL, fill = NULL)
   
   if(save) ggsave(paste0(here, "/trajectory", tag, ".pdf"), p_traj, width = 8, height = 4, units = "in")
   
   p_traj
}


as_species <- function(fit){
   
   stages <- c("p", "s", "j", "a", "x")
   
   ii <- sample(nrow(fit$draws()), 1)
   
   fd <- fit$draws() %>%
      as.data.frame() %>% as_tibble() %>%
      mutate(i = 1:nrow(.)) %>%
      gather(variable, value, -i) %>%
      filter(i == ii) %>%
      select(-i)
   
   # alpha
   p <- fd %>%
      filter(str_detect(variable, "A")) %>%
      mutate(variable = str_remove_all(variable, "A|\\[|\\]")) %>%
      separate(variable, c("s0", "s1"), sep = ",") %>%
      mutate(s0 = colnames(sp$p)[as.integer(s0)],
             s1 = rownames(sp$p)[as.integer(s1)]) %>%
      group_by(s0) %>%
      mutate(value = softmax(value)) %>%
      ungroup() %>%
      spread(s0, value) %>%
      column_to_rownames("s1")
   p <- p[match(stages, rownames(p)), match(stages[1:4], colnames(p))]
   
   # beta curvature
   bc <- fd %>%
      filter(str_detect(variable, "BC")) %>%
      mutate(variable = str_remove_all(variable, "BC|\\[|\\]")) %>%
      separate(variable, c("v", "s0", "s1"), sep = ",") %>%
      mutate(s0 = colnames(p)[as.integer(s0)],
             s1 = rownames(p)[as.integer(s1)],
             v = dimnames(sp$beta_crv)[[1]][as.integer(v)])
   pred_names <- unique(bc$v)
   beta_crv <- array(NA,
               dim = c(n_pred, ncol(p), nrow(p)),
               dimnames = list(pred_names, colnames(p), rownames(p)))
   for(i in dimnames(beta_crv)[[1]]){
      for(j in dimnames(beta_crv)[[2]]){
         for(k in dimnames(beta_crv)[[3]]){
            beta_crv[i,j,k] <- bc %>% filter(v == i, s0 == j, s1 == k) %>% pull(value)
         }
      }
   }
   
   # beta location
   bc <- fd %>%
      filter(str_detect(variable, "BL")) %>%
      mutate(variable = str_remove_all(variable, "BL|\\[|\\]")) %>%
      separate(variable, c("v", "s0", "s1"), sep = ",") %>%
      mutate(s0 = colnames(p)[as.integer(s0)],
             s1 = rownames(p)[as.integer(s1)],
             v = dimnames(sp$beta_crv)[[1]][as.integer(v)])
   pred_names <- unique(bc$v)
   beta_loc <- array(NA,
                     dim = c(n_pred, ncol(p), nrow(p)),
                     dimnames = list(pred_names, colnames(p), rownames(p)))
   for(i in dimnames(beta_loc)[[1]]){
      for(j in dimnames(beta_loc)[[2]]){
         for(k in dimnames(beta_loc)[[3]]){
            beta_loc[i,j,k] <- bc %>% filter(v == i, s0 == j, s1 == k) %>% pull(value)
         }
      }
   }
   
   # gamma
   bc <- fd %>%
      filter(str_detect(variable, "G")) %>%
      mutate(variable = str_remove_all(variable, "G|\\[|\\]")) %>%
      separate(variable, c("sx", "s0", "s1"), sep = ",") %>%
      mutate(s0 = colnames(p)[as.integer(s0)],
             s1 = rownames(p)[as.integer(s1)],
             sx = colnames(p)[as.integer(sx)])
   gamma <- array(0, 
                  dim = c(ncol(p), ncol(p), nrow(p)),
                  dimnames = list(colnames(p), colnames(p), rownames(p)))
   for(i in dimnames(gamma)[[1]]){
      for(j in dimnames(gamma)[[2]]){
         for(k in dimnames(gamma)[[3]]){
            gamma[i,j,k] <- bc %>% filter(sx == i, s0 == j, s1 == k) %>% pull(value)
         }
      }
   }
   
   
   list(p = p,
        fecundity = 500,
        beta_crv = beta_crv,
        beta_loc = beta_loc,
        gamma = gamma)
}



## execute ######################

# setup
prior <- make_priors()
sp <- make_species(prior)
d <- simulate_data(sp, n_blocks = 100, n0 = 1e5, t = 5)

# visualize inputs
plot_betas(sp) | plot_individual_outcomes(d)
plot_gammas(sp)
plot_simulated_ts(sp, n0 = matrix(c(1000, 100, 10, 10), 1), type = "lines")
plot_simulated_ts(sp, n0 = matrix(c(1000, 100, 10, 10), 1), type = "heatmap")

# fit model
fit <- fit_model(d, prior, algorithm = "variational")

# visualize outputs
plot_param_fits(fit, sp)
plot_trajectories(fit, sp, prior, d, reps = 100)

plot_betas(sp) | plot_betas(as_species(fit))
plot_gammas(sp) / plot_gammas(as_species(fit))






# what if we got fecundity wrong?
d2 <- d
d2$fecundity <- sp$fecundity / 10 # fecundity used in fitting (d) versus in simulating (sp)
fit <- fit_model(d2, prior, algorithm = "variational")
plot_param_fits(fit, sp)
plot_trajectories(fit, sp, prior, d2, reps = 100)
