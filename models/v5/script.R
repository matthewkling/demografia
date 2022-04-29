
# this version incorporates fecundity and unknown seedling turnover
# but drops environmental effects for the time being

## setup #################
library(tidyverse)
library(cmdstanr)
library(patchwork)

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
                  rmultinom(1, x[j], p[,j])
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





## projection matrix priors #############

# see alpha_priors.r for testing and vis

make_priors <- function(stages){
      
      stages <- c("p", "s", "j", "a") # propagule-seedling-juvenile-adult-dead
      labels <- c("propagules", "seedlings", "saplings", "adults")
      names(labels) <- stages
      k <- length(stages)
      
      # prior expectations for projection matrix -- origin in COLS, outcome in ROWS
      pp <- matrix(0, k+1, k, dimnames = list(c(stages, "x"), stages))
      pp["s", "p"] <- .002
      pp["s", "s"] <- .3
      pp["j", "s"] <- .02
      pp["j", "j"] <- .95
      pp["a", "j"] <- .01
      pp["a", "a"] <- .995
      pp["x",] <- apply(pp, 2, function(x) 1 - sum(x))
      pp <- apply(pp, 2, function(x) softmax(log(x)))
      
      alpha_i <- ceiling(t(pp)) # boolean indicating nonzero transitions
      
      mu <- log(pp)
      sigma <- mu
      sigma[] <- rep(c(2, 1, 2, 2), each = nrow(pp))
      
      alpha_mu <- mu[is.finite(mu)]
      alpha_sigma <- sigma[is.finite(mu)]
      
      list(p = pp,
           mu = mu,
           sigma = sigma,
           alpha_i = alpha_i,
           alpha_mu = alpha_mu,
           alpha_sigma = alpha_sigma,
           stages = labels)
}

prior <- make_priors(stages)




## true matrix, sampled from prior ###########

make_species <- function(prior){
      p <- prior$p
      for(i in 1:ncol(p)){
            nz <- which(p[,i] > 0)
            p[nz, i] <- softmax(rnorm(length(nz), prior$mu[nz, i], prior$sigma[nz, i]))
      }
      
      fecundity <- 500
      
      list(p = p,
           fecundity = fecundity)
}

sp <- make_species(prior)



## simulate dataset ############################

simulate_data <- function(sp, t = 5){
      
      # collective counts
      n0 <- c(300, 10, 10, 3) * 10000
      nts <- sim(n0, sp$p, sp$fecundity, t)
      nt <- nts[t+1,] %>% select(-t) %>% as.matrix() %>% as.vector()
      
      # individual outcomes
      io_stages <- 3:4
      p_sym <- cbind(sp$p, x = c(0, 0, 0, 0, 1))
      io <- io_stages %>%
            sapply(function(i) rmultinom(1, n0[i], project_simplex(p_sym, t, i))) %>%
            t()
      
      list(n0 = n0,
           nt = nt,
           io = io,
           io_stages = io_stages,
           t = t)
}

d <- simulate_data(sp)



## fit stan model ###################

fit_model <- function(sp, d, algorithm = "sample"){
      
      # prepare input data
      dat <- list(
            K = ncol(sp$p),
            n0 = d$n0,
            nt = d$nt,
            nk = 2,
            io = d$io,
            iok = d$io_stages,
            t = d$t,
            "F" = sp$fecundity,
            Fk = match("a", rownames(sp$p)),
            alpha_i = prior$alpha_i,
            n_alpha = length(prior$alpha_mu),
            alpha_mu = prior$alpha_mu,
            alpha_sigma = prior$alpha_sigma
      )
      
      model <- cmdstan_model("models/v5/model.stan")
      
      if(algorithm == "sample"){
            fit <- model$sample(
                  data = dat,
                  chains = 1
                  # adapt_delta = .99
                  # iter_warmup = 200, iter_sampling = 100, refresh = 30
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

fit <- fit_model(sp, d)




## assess fit ##############################

## true vs fitted parameters ##

alphabet_plot <- function(fit, sp){
      
      pd <- sp$p %>% 
            as.data.frame() %>% 
            rownames_to_column("s1") %>% as_tibble() %>%
            gather(s0, true, -s1)
      
      fa <- fit$summary() %>%
            filter(str_detect(variable, "A")) %>%
            mutate(variable = str_remove_all(variable, "A|\\[|\\]")) %>%
            separate(variable, c("s0", "s1"), sep = ",") %>%
            mutate(s0 = colnames(sp$p)[as.integer(s0)],
                   s1 = rownames(sp$p)[as.integer(s1)]) %>%
            left_join(pd) %>%
            group_by(s0) %>%
            mutate(pred = mean) 
      
      ggplot(fa %>% filter(true > 0, paste0(s0, s1) != "xx"), 
             aes(logit(true), logit(pred), ymin = logit(q5), ymax = logit(q95))) + 
            geom_abline(slope = 1, intercept = 0, alpha = .3) +
            ggrepel::geom_text_repel(aes(label = paste0(s0, s1)), 
                                     direction = "x", segment.size = .25) +
            geom_linerange(size = .5, color = "darkred") +
            geom_point(color = "darkred") +
            theme_minimal() +
            labs(title = expression(alpha))
}

alphabet_plot(fit, sp)


## simulated projections under true vs fitted params ##

trajectory_plot <- function(fit, sp, prior, t = 100, reps = 100){
      
      # starting population for simulation
      n0sim <- d$n0 / 1000
      
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
                  yy <- try(sim(n0sim, y, sp$fecundity, t) %>%
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
                  strip.text = element_text(color = "white", size = 16)) +
            labs(y = paste0("N(t): 95% CI across ", reps, " simulations"),
                 x = "t",
                 color = NULL, fill = NULL)
      p_traj
}

trajectory_plot(fit, sp, prior, reps = 100)




