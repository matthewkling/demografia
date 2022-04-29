
# this version reintroduces environmental dependence
# and implements other streamlines and bugfixes


## setup #################

library(tidyverse)
library(cmdstanr)
library(patchwork)
library(popbio)

here <- "models/v6"

## utils ####################

logit <- function(x) log(x / (1-x))

softmax <- function(x) exp(x) / sum(exp(x))
# softmax(log(c(.1, .2, .3, .4))) # demo that log is the inverse of softmax

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
      xi <- sapply(1:k, function(j){
         xj <- x[j]
         if(xj <= .Machine$integer.max) return(rmultinom(1, xj, p[,j]))
         if(xj > .Machine$integer.max) return(round(xj * p[,j]))
      }) %>%
         apply(1, sum)
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
   p["s", "p"] <- .002
   p["s", "s"] <- .3
   p["j", "s"] <- .02
   p["j", "j"] <- .95
   p["a", "j"] <- .01
   p["a", "a"] <- .995
   p["x",] <- apply(p, 2, function(x) 1 - sum(x))
   p <- apply(p, 2, function(x) softmax(log(x)))
   
   alpha_i <- ceiling(t(p)) # boolean indicating nonzero transitions
   
   mu <- log(p)
   sigma <- mu
   sigma[] <- rep(c(2, 1, 2, .75), each = nrow(p))
   
   alpha_mu <- mu[is.finite(mu)]
   alpha_sigma <- sigma[is.finite(mu)]
   
   
   ## betas: effect of environmental variables on log prob ##
   
   n_pred <- 2 # number of predictors
   pred_names <- paste0("v", 1:n_pred)
   
   beta <- array(0, 
                 dim = c(n_pred, ncol(p), nrow(p)),
                 dimnames = list(pred_names, colnames(p), rownames(p)))
   
   # flag betas to be nonzero
   beta[, "a", "a"] <- 1
   beta[, "j", "j"] <- 1
   beta[, "j", "a"] <- 1
   
   # standard deviation of betas (mean is zero)
   beta_sd <- 1
   
   # bundle
   list(p = p,
        mu = mu,
        sigma = sigma,
        alpha_i = alpha_i,
        alpha_mu = alpha_mu,
        alpha_sigma = alpha_sigma,
        beta = beta,
        beta_sd = beta_sd,
        stages = labels)
}





## true species parameters, sampled from prior ###########

make_species <- function(prior){
   
   fecundity <- 500
   
   # alpha
   p <- prior$p
   for(i in 1:ncol(p)){
      nz <- which(p[,i] > 0)
      p[nz, i] <- softmax(rnorm(length(nz), prior$mu[nz, i], prior$sigma[nz, i]))
   }
   
   # beta
   b <- prior$beta
   b[b == 1] <- rnorm(sum(b), 0, prior$beta_sd)
   
   # bundle
   list(p = p,
        fecundity = fecundity,
        beta = b)
}




## simulate dataset ############################

simulate_data <- function(sp, n_blocks = 100, t = 5, n0 = 10000){
   
   io_stages <- 3:4 # stages for which individual outcomes are tracked
   
   p <- sp$p
   f <- sp$fecundity
   b <- sp$beta
   
   # predictor data
   n_pred <- dim(b)[1]
   x <- matrix(rnorm(n_blocks * n_pred), 
               n_blocks, 2, 
               dimnames = list(paste0("b", 1:n_blocks),
                               dimnames(b)[[1]]))
   
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
               pi[s1, s0] <- pi[s1, s0] + b[v, s0, s1] * x[i, v]
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
      # n_t0 <- c(1000, 100, 10, 10)
      n[i, , 1] <- n_t0
      
      # collective counts
      n_ts <- sim(n_t0, pi, f, t)
      n_tt <- n_ts[t+1, ] %>% select(-t) %>% as.matrix() %>% as.vector()
      n[i, , 2] <- n_tt
      
      # individual outcomes
      A[1, 4] <- 0
      io[i, , ] <- io_stages %>%
         sapply(function(j) rmultinom(1, n_t0[j], project_simplex(A, t, j))) %>%
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
      t = d$t,
      "F" = d$fecundity,
      Fk = match("a", dimnames(d$n)[[2]]),
      x = d$x,
      alpha_i = prior$alpha_i,
      n_alpha = length(prior$alpha_mu),
      alpha_mu = prior$alpha_mu,
      alpha_sigma = prior$alpha_sigma,
      n_beta = sum(prior$beta),
      beta_i = prior$beta,
      beta_sigma = prior$beta_sd
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

alphabet_plot <- function(fit, sp, 
                          save = F, tag = NULL){
   
   ## alphas ##
   
   pd <- sp$p %>% 
      as.data.frame() %>% 
      rownames_to_column("s1") %>% as_tibble() %>%
      gather(s0, true, -s1)
   
   pd <- fit$summary() %>%
      filter(str_detect(variable, "A")) %>%
      mutate(variable = str_remove_all(variable, "A|\\[|\\]")) %>%
      separate(variable, c("s0", "s1"), sep = ",") %>%
      mutate(s0 = colnames(sp$p)[as.integer(s0)],
             s1 = rownames(sp$p)[as.integer(s1)]) %>%
      left_join(pd) %>%
      group_by(s0) %>%
      mutate(pred = mean,
             trans = paste0(s0, s1),
             type = ifelse(str_detect(trans, "j|a"), "individual", "collective")) %>% 
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
   
   bd <- sp$beta %>%
      as.data.frame() %>% 
      rownames_to_column("v") %>% as_tibble() %>%
      gather(par, true, -v) %>%
      separate(par, c("s0", "s1"), sep = "\\.")
   
   bd <- fit$summary() %>%
      filter(str_detect(variable, "B")) %>%
      mutate(variable = str_remove_all(variable, "B|\\[|\\]")) %>%
      separate(variable, c("v", "s0", "s1"), sep = ",") %>%
      mutate(s0 = colnames(sp$p)[as.integer(s0)],
             s1 = rownames(sp$p)[as.integer(s1)],
             v = dimnames(sp$beta)[[1]][as.integer(v)]) %>%
      left_join(bd) %>%
      filter(true != 0)
   
   beta <- bd %>%
      ggplot(aes(true, median, ymin = q5, ymax = q95)) +
      geom_abline(slope = 1, intercept = 0, alpha = .3) +
      ggrepel::geom_text_repel(aes(label = paste0(v, "-", s0, s1)), 
                               direction = "x", segment.size = .25) +
      geom_linerange(size = .5, color = "darkgreen") +
      geom_point(color = "darkgreen") +
      theme_minimal() +
      labs(title = expression(beta),
           x = "true effect",
           y = "posterior effect\n(median & 95% CI)")
   
   p_ab <- alpha + beta +
      plot_layout(nrow = 1) &
      theme(plot.title = element_text(hjust = .5, size = 30))
   
   if(save) ggsave(paste0(here, "/alpha_beta", tag, ".pdf"), p_ab, width = 8, height = 4, units = "in")
   
   p_ab
}




## simulated projections under true vs fitted params ##

trajectory_plot <- function(fit, sp, prior, d, 
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
         yy <- try(sim(n0sim, y, d$fecundity, t) %>% # using d$fecundity rather than sp$fecundity in case they differ
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


## run ######################

prior <- make_priors()
sp <- make_species(prior)
d <- simulate_data(sp, n_blocks = 100, n0 = 1e3, t = 5)
fit <- fit_model(d, prior, algorithm = "variational")
alphabet_plot(fit, sp)
trajectory_plot(fit, sp, prior, d, reps = 100)



# what if we got fecundity wrong?
d2 <- d
d2$fecundity <- sp$fecundity / 10 # fecundity used in fitting (d) versus in simulating (sp)
fit <- fit_model(d2, prior, algorithm = "variational")
alphabet_plot(fit, sp)
trajectory_plot(fit, sp, prior, d2, reps = 100)
