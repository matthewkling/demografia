




get_FIA_data <- function(sp = "Acer saccharum"){
  
  require(tidyverse)
  require(raster)
  select <- dplyr::select
  
  
  ## generic FIA data prep ##
  
  t <- readRDS("data/derived/fia_trees.rds") %>%
    mutate(tree_id = paste(subplot_id, designcd, tree)) %>%
    filter(manual >= 2.0)
  
  s <- readRDS("data/derived/fia_seedlings.rds") %>%
    filter(manual >= 2.0) # seedlings only counted with M2+
  
  
  ## individual outcomes ##
  
  # note: still need to figure out whether there are implicit deaths,
  # in addition to the explicit deaths counted here
  
  io <- t %>% 
    filter(!is.na(dia),
           paste(genus, species) == sp) %>%
    mutate(status = recode(statuscd, 
                           "1" = "live", "2" = "dead", "3" = "harvested",
                           .default = NA_character_),
           class = ifelse(dia >= 5, "tree", "sapling")) %>%
    group_by(tree_id) %>%
    arrange(invyr) %>%
    mutate(dia_next = lead(dia),
           invyr_next = lead(invyr),
           class_next = lead(class),
           status_next = lead(status)) %>%
    ungroup() %>%
    filter(is.finite(invyr),
           is.finite(invyr_next),
           !is.na(status),
           status == "live",
           status_next != "harvested") %>%
    mutate(trans = case_when(status == status_next & class == class_next ~ "survived",
                             status == "live" & status_next == "dead" ~ "died",
                             class == "sapling" & class_next == "tree" ~ "graduated",
                             class == "tree" & class_next == "sapling" ~ "reverted",
                             TRUE ~ "other"),
           years = invyr_next - invyr) %>%
    count(plot_id, invyr, years, class, trans) %>%
    mutate(s0 = case_when(class == "tree" ~ "a",
                          class == "sapling" ~ "j"),
           s1 = case_when(trans == "died" ~ "x",
                          trans == "graduated" & s0 == "j" ~ "a",
                          trans %in% c("survived", "reverted") ~ s0)) %>%
    select(plot_id, invyr, t = years, s0, s1, n) %>%
    mutate(inv = paste(plot_id, invyr))
  
  
  ## collective counts ##
  
  ss <- s %>%
    filter(paste(genus, species) == sp) %>%
    select(plot_id, invyr, treecount) %>%
    mutate(class = "s") %>%
    group_by(plot_id, invyr, class) %>%
    summarize(n = sum(treecount))
    
  cc <- t %>%
    filter(!is.na(dia),
           paste(genus, species) == sp) %>%
    mutate(status = recode(statuscd, 
                           "1" = "live", "2" = "dead", "3" = "harvested",
                           .default = NA_character_),
           class = ifelse(dia >= 5, "a", "j")) %>%
    group_by(plot_id) %>%
    filter(!any(status == "harvested")) %>% # exclude all plots with any harvesting
    ungroup() %>%
    filter(status == "live") %>%
    group_by(plot_id, invyr, class) %>%
    summarize(n = n()) %>%
    ungroup() %>%
    bind_rows(ss)
  
  cc <- expand_grid(plot_id = unique(cc$plot_id),
              invyr = unique(cc$invyr),
              class = unique(cc$class)) %>%
    left_join(cc) %>%
    mutate(n = ifelse(is.na(n), 0, n)) %>%
    group_by(plot_id, invyr) %>%
    filter(sum(n) > 0) %>%
    group_by(plot_id) %>%
    filter(length(unique(invyr)) > 1) %>%
    ungroup()
  
  cc <- cc %>%
    arrange(plot_id, class, invyr) %>%
    group_by(plot_id, class) %>%
    mutate(n_next = lead(n),
           invyr_next = lead(invyr)) %>%
    filter(is.finite(n_next)) %>%
    ungroup() %>%
    mutate(inv = paste(plot_id, invyr),
           t = invyr_next - invyr) %>%
    select(-invyr_next)
  
  
  
  ## common inventories ##
  
  master_inv <- intersect(unique(io$inv), unique(cc$inv))
  io <- filter(io, inv %in% master_inv)
  cc <- filter(cc, inv %in% master_inv)
  
  
  
  ## climate data ##
  
  ll <- t %>% 
    filter(plot_id %in% unique(io$plot_id)) %>%
    dplyr::select(plot_id, lon, lat) %>% 
    distinct()
  
  climate <- list.files("f:/chelsa/v2/climatology/raw/", pattern = "_bio1",
                        full.names = T) %>%
    stack() %>%
    extract(dplyr::select(ll, lon, lat) %>% as.matrix()) %>%
    as.data.frame() %>%
    as_tibble() %>%
    setNames(c("temp", "prec")) %>%
    mutate(plot_id = ll$plot_id)
  
  io <- left_join(io, climate)
  cc <- left_join(cc, climate)
  
  
  
  
  ## summarize to modeling blocks ##
  
  bcc <- cc %>%
    rename(t0 = n,
           tt = n_next) %>%
    gather(time, n, tt, t0) %>%
    unite(var, class, time) %>%
    distinct() %>%
    select(-plot_id, -invyr) %>%
    spread(var, n)
  
  bio <- io %>%
    select(-plot_id, -invyr) %>%
    unite(var, s0, s1) %>%
    group_by(inv, var, t, temp, prec) %>%
    summarize(n = sum(n)) %>%
    ungroup() %>%
    spread(var, n, fill = 0)
  
  
  log_bin <- function(x, res = .5) round(10 ^ plyr::round_any(log10(x), res))
  
  clim_bins <- 5
  
  b <- inner_join(bcc, bio) %>%
    mutate(temp = plyr::round_any(temp, (max(temp) - min(temp)) / clim_bins),
           prec = plyr::round_any(prec, (max(prec) - min(prec)) / clim_bins)) %>%
    mutate(n0_a = log_bin(a_t0),
           n0_j = log_bin(j_t0),
           n0_s = log_bin(s_t0)) %>%
    group_by(temp, prec, t, n0_a, n0_j, n0_s) %>%
    select(-inv) %>%
    mutate(n_inv = 1) %>%
    summarize_all(sum)
  
  
  ## add hectares surveyed ##
  
  # area surveyed per PLOT, in hectares 
  # (these values /4 give areas surveyed per subplot/microplot)
  area <- tibble(class = c("s_area", "j_area", "a_area"),
                 ha = base::pi * (c(6.8, 6.8, 24) / 3.28084) ^ 2 / 10000 * 4) %>%
    spread(class, ha)
  
  b <- b %>%
    bind_cols(area) %>%
    mutate(a_area = a_area * n_inv,
           j_area = j_area * n_inv,
           s_area = s_area * n_inv) %>%
    ungroup()
  
  
  ## format for model ##
  
  # environment
  x <- b %>% 
    select(v1 = temp, v2 = prec) %>%
    as.matrix()
  x_scale <- list(v1_mean = mean(x[,1]),
                  v2_mean = mean(x[,2]),
                  v1_sd = sd(x[,1]),
                  v2_sd = sd(x[,2]))
  x <- scale(x)
  
  # individual outcomes
  iod <- b %>%
    select(a_a:j_x) %>%
    mutate(a_j = 0,
           block = 1:nrow(.)) %>%
    gather(var, n, -block) %>%
    separate(var, c("s0", "s1")) %>%
    arrange(s0, s1, block)
  io <- array(0, 
              c(max(iod$block), 2, 5),
              dimname = list(paste0("b", 1:max(iod$block)),
                             c("j", "a"),
                             c("p", "s", "j", "a", "x")))
  for(is0 in unique(iod$s0)){
    for(is1 in unique(iod$s1)){
      io[ ,is0, is1] <- iod %>% filter(s0 == is0, s1 == is1) %>% arrange(block) %>% pull(n)
    }
  }
  
  # collective counts
  ccd <- b %>%
    select(a_t0:s_tt) %>%
    mutate(block = 1:nrow(.)) %>%
    gather(var, n, -block) %>%
    separate(var, c("class", "t")) %>%
    arrange(class, t, block)
  cc <- array(0, 
              c(max(ccd$block), 3, 2),
              dimname = list(paste0("b", 1:max(ccd$block)),
                             c("s", "j", "a"),
                             c("t0", "tt")))
  for(cl in unique(ccd$class)){
    for(tt in unique(ccd$t)){
      cc[ ,cl, tt] <- ccd %>% filter(class == cl, t == tt) %>% arrange(block) %>% pull(n)
    }
  }
  
  # time intervals
  t <- b$t
  
  # land area
  area <- b %>%
    select(s_area, j_area, a_area) %>%
    as.matrix()
    
  
  
  ## bundle for export ##
  
  list(K = 4,
       x = x,
       x_scale = x_scale,
       n = cc,
       io = io,
       io_stages = 3:4,
       t = t,
       fecundity = 500,
       area = area)
  
}

d <- get_FIA_data("Acer saccharum") %>% 
  saveRDS("data/derived/species/Acer saccharum.rds")
