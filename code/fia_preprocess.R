library(data.table)
library(tidyverse)
library(janitor)




## load FIA tables ############

# PLOT table
# primary key: CN
# unique key: STATECD, INVYR, UNITCD, COUNTYCD, PLOT
plot <- fread("data/raw/fia/ENTIRE/PLOT.csv", stringsAsFactors = F, 
              colClasses = c(CN = "character")) %>% 
  select(PLT_CN = CN, MANUAL, DESIGNCD, PLOT_STATUS_CD,
         STATECD, INVYR, UNITCD, COUNTYCD, PLOT, 
         LAT, LON, ELEV) %>% as_tibble()

# SUBPLOT table
# primary key: CN
# unique key: PLT_CN, SUBP
# When PLOT.MANUAL <1.0, the field crew measured slope at a condition level
# When PLOT.MANUAL >=1.0, slope is collected at a subplot level, and then the aspect 
# from the subplot representing the greatest proportion of the condition is assigned as a surrogate
subplot <- fread("data/raw/fia/ENTIRE/SUBPLOT.csv", stringsAsFactors = F, 
                 colClasses = c(CN = "character", PLT_CN = "character")) %>% 
  select(PLT_CN, SUBP, SUBP_CN = CN, PREV_SBP_CN, POINT_NONSAMPLE_REASN_CD,
         SLOPE, ASPECT) %>% as_tibble()

# CONDITION table (used to identify and exclude planted plots)
# primary key: CN
# unique key: PLT_CN, CONDID
cond <- fread("data/raw/fia/ENTIRE/COND.csv", stringsAsFactors = F, 
              colClasses = c(CN = "character", PLT_CN = "character")) %>% 
  select(PLT_CN, CONDID, PROP_BASIS, STDORGCD, PHYSCLCD, COND_STATUS_CD, CONDPROP_UNADJ,
         SLOPE_COND = SLOPE, ASPECT_COND = ASPECT) %>% as_tibble()
planted <- cond %>%
  filter(STDORGCD == 1) %>%
  select(PLT_CN) %>%
  distinct() %>%
  mutate(planted = TRUE)

# TREE table
# primary key: CN
# unique key: PLT_CN, SUBP, TREE
# Tree numbers (TREE) can be used to track trees when PLOT.DESIGNCD is the same between inventories.
tree <- fread("data/raw/fia/ENTIRE/TREE.csv", 
              stringsAsFactors = F, 
              select = c("PLT_CN", "PLOT", "SUBP", "INVYR", "TREE", "SPCD", 
                         "DIA", "DIAHTCD", "STATUSCD", "DIACHECK"),
              colClasses = c(PLT_CN = "character")) %>% as_tibble()

# SEEDLING table
# primary key: CN
# unique key: PLT_CN, SUBP, CONDID, SPCD
seedling <- fread("data/raw/fia/ENTIRE/SEEDLING.csv", 
              stringsAsFactors = F, 
              select = c("PLT_CN", "PLOT", "SUBP", "INVYR", "SPCD", 
                         "TREECOUNT", "TOTAGE", "STOCKING"),
              colClasses = c(PLT_CN = "character")) %>% as_tibble()


# TREE_WOODLAND_STEMS table -- placeholder

# SPECIES table
species <- fread("data/raw/fia/REF_SPECIES.csv", stringsAsFactors = F) %>%
  select(SPCD, COMMON_NAME, GENUS, SPECIES) %>% as_tibble()




## join tables ########

d <- plot %>%
  left_join(planted) %>%
  left_join(subplot) %>%
  left_join(tree) %>%
  # left_join(seedling) %>%
  left_join(species) %>% 
  mutate(PLOT_ID = paste(STATECD, UNITCD, COUNTYCD, PLOT),
         SUBPLOT_ID = paste(PLOT_ID, SUBP)) %>%
  filter(is.na(planted),
         DESIGNCD %in% c(1, 111:118, 230:242, 311:323, 328, 501:506), # 4-subplot design
         SUBP <= 4) %>%  # remove extraneous subplots
  # select(PLT_CN, STATECD, UNITCD, COUNTYCD, PLOT, SUBP, INVYR, DESIGNCD,
  #        PLOT_ID, SUBPLOT_ID, LON, LAT, ELEV, PLOT_STATUS_CD, POINT_NONSAMPLE_REASN_CD,
  #        SLOPE, ASPECT, GENUS, SPECIES) %>%
  clean_names() %>%
  mutate(lon = ifelse(lon>0, lon-360, lon)) %>%
  filter(lat > 24)

saveRDS(d, "data/derived/fia_trees.rds")


s <- plot %>%
  left_join(planted) %>%
  left_join(subplot) %>%
  # left_join(tree) %>%
  left_join(seedling) %>%
  left_join(species) %>% 
  mutate(PLOT_ID = paste(STATECD, UNITCD, COUNTYCD, PLOT),
         SUBPLOT_ID = paste(PLOT_ID, SUBP)) %>%
  filter(is.na(planted),
         DESIGNCD %in% c(1, 111:118, 230:242, 311:323, 328, 501:506), # 4-subplot design
         SUBP <= 4) %>%  # remove extraneous subplots
  clean_names() %>%
  mutate(lon = ifelse(lon>0, lon-360, lon)) %>%
  filter(lat > 24)

saveRDS(s, "data/derived/fia_seedlings.rds")

