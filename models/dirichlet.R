

library(MCMCpack)
library(ggtern)


rdirichlet(1000, c(10, 1, .1)) %>%
  as.data.frame() %>%
  as_tibble() %>%
  ggtern(aes(V1, V2, V3)) +
  geom_point(alpha = .25) +
  # theme_minimal() +
  theme_showarrows()
