
t <- readRDS("data/derived/fia_trees.rds") %>%
  mutate(tree_id = paste(subplot_id, designcd, tree))

s <- readRDS("data/derived/fia_seedlings.rds")


sp <- 318 # sugar maple





# identify (sub?)plots surveyed in 2+ years under same designcd
p <- bind_rows(t, s) %>%
  select(plot_id, lon, lat, invyr, spcd, designcd) %>%
  distinct() %>%
  mutate(series = paste(plot_id, designcd)) %>%
  group_by(series, lon, lat) %>%
  filter(sp %in% spcd) %>%
  filter(length(unique(invyr)) > 1) %>%
  ungroup() %>%
  select(-spcd) %>%
  distinct()

mpd <- map_data("state")

plt <- p %>%
  count(series, lon, lat) %>%
  # summarize(n = length(unique(invyr)),
  #           lon = mean(lon),
  #           lat = mean(lat)) %>%
  ggplot(aes(lon, lat, color = n)) +
  geom_polygon(data = mpd, aes(long, lat, group = group),
               fill = "gray80", color = NA) +
  geom_point() +
  scale_color_viridis_c() +
  coord_cartesian(xlim = range(p$lon), ylim = range(p$lat)) +
  theme_minimal() +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank()) +
  labs(color = "n surveys",
       title = "FIA plots that contain sugar maple and were surveyed 2+ times under a single designcd")
ggsave("figures/acer_saccharum/map_nyears.png", plt,
       width = 8, height = 6, units = "in")

x <- p %>%
  group_by(series) %>%
  filter(length(unique(invyr)) > 1) %>%
  arrange(invyr) %>%
  summarize(years = paste(invyr, collapse = "_")) %>%
  ungroup() %>%
  select(years) %>%
  distinct() %>%
  separate(years, paste0("y", 1:6), remove = F) %>%
  gather(i, year, y1:y6) %>% #filter(years == years[1])
  na.omit() %>%
  mutate(year = as.integer(year)) %>%
  filter(year < 2023) %>%
  group_by(years) %>%
  mutate(n = length(year))
plt <- ggplot(x, aes(year, years, group = years)) +
  facet_wrap(~n, scales = "free", nrow = 1, labeller = label_both) +
  geom_line(size = .25) + 
  geom_point(size = .75) +
  theme_bw() +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank()) +
  labs(title = "unique visit schedules for FIA plots containing sugar maple")
ggsave("figures/acer_saccharum/invyr.png", plt,
       width = 8, height = 5, units = "in")


# tree-seedling co-observation patterns (tangent)
bind_rows(t %>% mutate(stage = "tree"), 
          s %>% mutate(stage = "seedling")) %>%
  mutate(series = paste(plot_id, designcd)) %>%
  filter(series %in% unique(p$series)) %>%
  group_by(series, stage, invyr) %>%
  summarize(pres = sp %in% spcd) %>%
  ungroup() %>%
  unite(occ, stage, pres) %>%
  group_by(series, invyr) %>%
  summarize(occ = paste(occ, collapse = ", ")) %>%
  ungroup() %>%
  count(occ)



x <- t %>% 
  filter(!is.na(dia),
         spcd == sp) %>%
  filter(tree_id %in% sample(tree_id, 1000)) %>%
  count(tree_id)
pd <- t %>% 
  filter(tree_id %in% sample(x$tree_id[x$n >= 2], 100)) %>%
  mutate(status = recode(statuscd, 
                         "1" = "live", "2" = "dead", "3" = "harvested"))
plt <- pd %>%
  ggplot(aes(invyr, dia, group = tree_id, color = status)) +
  geom_hline(yintercept = 5, color = "dodgerblue") +
  geom_line(size = .25) + 
  geom_point(size = .75) +
  scale_y_log10() +
  scale_color_manual(values = c("red", "green", "black"), 
                     na.value = "yellow") +
  theme_minimal() +
  labs(title = "DBH tracjectories for random sugar maples")
ggsave("figures/acer_saccharum/dbh_trajectories.png", plt,
       width = 7, height = 5, units = "in")



# note: still need to figure out whether there are implicit deaths,
# in addition to the explicit deaths counted here


x <- t %>% 
  filter(!is.na(dia),
         spcd == sp) %>%
  mutate(status = recode(statuscd, "1" = "live", "2" = "dead", "3" = "harvested"),
         class = ifelse(dia >= 5, "tree", "sapling")) %>%
  group_by(tree_id) %>%
  arrange(dia) %>%
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
         years = invyr_next - invyr) 

x %>%
  group_by(status, class, trans) %>%
  summarize(n_events = n(),
            n_years = sum(years)) %>%
  group_by(status, class) %>%
  mutate(p_event = n_events / sum(n_events),
         p_year = n_events / sum(n_years),
         p_year = ifelse(trans == "survived", 
                         1-sum(p_year[trans != "survived"]), 
                         p_year))


## model fits ##############

library(nnet)
library(raster)

# climate data
ll <- x %>% dplyr::select(plot_id, lon, lat) %>% distinct()
climate <- list.files("f:/chelsa/v2/climatology/raw/", pattern = "_bio1",
                      full.names = T) %>%
  stack()
climate <- extract(climate, dplyr::select(ll, lon, lat) %>% as.matrix()) %>%
  as.data.frame() %>%
  as.tibble() %>%
  setNames(c("temp", "prec")) %>%
  mutate(plot_id = ll$plot_id)
md <- x %>% 
  filter(class == "sapling") %>%
  left_join(climate) %>%
  mutate(rtrans = relevel(factor(trans), ref = "survived"))

# density data
density <- t %>% 
  filter(!is.na(dia)) %>%
  group_by(subplot_id, invyr) %>%
  summarize(n_trees = length(dia[dia >= 5 & spcd == sp]) %>% sqrt(),
            n_saplings = length(dia[dia < 5 & spcd == sp]) %>% sqrt())
md <- left_join(density, md)

# fit multinomial model
fit <- multinom(trans ~ dia + temp + prec + n_saplings + n_trees, 
                data = md)

pred <- md %>% 
  ungroup() %>%
  select(dia, temp, prec, n_saplings, n_trees) %>%
  na.omit() %>%
  bind_cols(predict(fit, ., type = "probs") %>%
              as.data.frame()) %>%
  gather(trans, value, died:survived) %>%
  mutate(n_saplings = n_saplings ^ 2,
         n_trees = n_trees ^ 2)

ggplot(pred %>% filter(round(dia, 0) == 4), 
       aes(n_trees, value, color = n_saplings)) +
  facet_wrap(~trans, scales = "free") +
  geom_point() +
  scale_color_viridis_c() +
  scale_x_sqrt() +
  theme_bw()




pred <- md %>% 
  ungroup() %>%
  filter(is.finite(temp), is.finite(dia)) %>%
  mutate(temp = median(temp),
         prec = median(prec),
         dia = 4,
         n_trees = round(n_trees, 1),
         n_saplings = round(n_saplings, 1)) %>%
  select(temp, prec, dia, n_trees, n_saplings) %>%
  distinct() %>%
  bind_cols(predict(fit, ., type = "probs") %>%
              as.data.frame()) %>%
  mutate(n_saplings = n_saplings ^ 2,
         n_trees = n_trees ^ 2) %>%
  gather(trans, value, died:survived)


plt <- ggplot(pred, 
       aes(n_trees, value, alpha = n_saplings,
           color = trans, group = paste(trans, n_saplings))) +
  facet_wrap(~ trans, scales = "free", nrow = 1) +
  geom_line(size = 1) +
  scale_alpha_continuous(trans = "sqrt", range = c(.3, 1)) +
  theme_bw() +
  labs(color = "transition",
       y = "probability",
       title = "multiniomial logistic fits of sugar maple sapling outcomes")
ggsave("figures/acer_saccharum/sapling_multinom_density.png", plt,
       width = 8, height = 5, units = "in")


########

pred <- md %>% mutate(temp = round(temp),
                      prec = plyr::round_any(prec, 500),
                      dia = round(dia, 0)) %>%
  select(temp, prec, dia) %>%
  distinct() %>%
  bind_cols(predict(fit, ., type = "probs") %>%
              as.data.frame()) %>%
  gather(trans, value, died:survived) %>%
  rename(dbh = dia)

plt <- ggplot(pred, aes(temp, value, alpha = prec, color = trans,
                        group = paste(trans, prec))) +
  facet_grid(. ~ dbh, scales = "free", labeller = label_both) +
  geom_line(size = 1) +
  scale_alpha_continuous(range = c(.25, 1)) +
  scale_color_manual(values = c("darkred", "darkgreen", "darkblue")) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.box = "vertical") +
  labs(color = "transition",
       alpha = "precipitation",
       x = "temperature",
       y = "probability",
       title = "multiniomial logistic fits of sugar maple sapling outcomes")
ggsave("figures/acer_saccharum/sapling_multinom.png", plt,
       width = 7, height = 5, units = "in")




pred <- climate %>%
  left_join(select(md, plot_id, lon, lat)) %>%
  mutate(dia = 4) %>%
  bind_cols(predict(fit, ., type = "probs") %>%
              as.data.frame()) %>%
  gather(trans, value, died:survived)

plt <- pred %>%
  filter(trans != "survived") %>%
  ggplot(aes(lon, lat, color = value)) +
  facet_wrap(~ trans) +
  geom_polygon(data = mpd, aes(long, lat, group = group),
               fill = "gray80", color = NA) +
  geom_point() +
  scale_color_viridis_c() +
  coord_cartesian(xlim = range(p$lon), ylim = range(p$lat)) +
  theme_minimal() +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank()) +
  labs(color = "prob",
       title = "sugar maple saplings: predicted transition probabilities based on climate")
ggsave("figures/acer_saccharum/sapling_multinom_map.png", plt,
       width = 10, height = 4, units = "in")




####################

# if we assume that seedlings are in equilibrium:
# - we can use their counts, plus appearances of new saplings, to calculate mean graduation rates
# - we can use their mean age to calculate turnover rates (but virtually never reported)


##  find instances of seedling-sapling transition ##

# identify all seedling survey events that were followed by sapling survey events
# count seedlings in those surveys (including on plots that had no sapling appearance), 
# and saplings that newly appeared in the following survey


# subplot time series with seedlings for focal species
seedlings <- s %>%
  filter(spcd == sp) %>%
  mutate(series_id = paste(subplot_id, designcd)) %>%
  select(series_id, invyr, treecount) %>%
  distinct()

# subplot time series in which saplings of any species were counted
series <- t %>%
  filter(dia < 5, !is.na(dia)) %>%
  mutate(series_id = paste(subplot_id, designcd)) %>%
  select(series_id, invyr) %>%
  distinct() %>%
  group_by(series_id) %>%
  summarize(series_start = min(invyr))

# new appearances of saplings of focal species
saplings <- t %>% 
  filter(spcd == sp,
         dia < 5, !is.na(dia)) %>%
  mutate(series_id = paste(subplot_id, designcd)) %>%
  select(series_id, tree_id, invyr, dia) %>%
  group_by(series_id, tree_id) %>%
  filter(invyr ==  min(invyr)) %>%
  rename(sapling_start = invyr) %>%
  left_join(series) %>%
  filter(sapling_start > series_start)


# seedling-sapling transitions
ss <- left_join(seedlings, saplings) %>%
  filter(is.finite(treecount),
         sapling_start > invyr | is.na(sapling_start)) %>%
  group_by(series_id) %>%
  mutate(n = n()) %>%
  ungroup() %>%
  group_by(series_id) %>%
  mutate(next_year = unique(c(invyr, sapling_start))[
    match(invyr, unique(c(invyr, sapling_start))) + 1]) %>%
  group_by(series_id, invyr) %>%
  summarize(seedlings = treecount[1],
            new_saplings = sum(sapling_start == next_year &
                                 !is.na(sapling_start))) %>%
  ungroup()

sum(ss$new_saplings) / sum(ss$seedlings)
