# -------------------------------------- #
#### Script for calculating d' ####
# Author: XXXXX
# -------------------------------------- #

librarian::shelf(bipartite, tidyverse, glmmTMB, ggeffects)

# load webs
webs.matrix <- readRDS(file = file.path("data", "webs.matrix.RData"))


### calculating species level d' ####
net_d <- lapply(webs.matrix, specieslevel, index = "d", level = "both")
net_d.df <- lapply(net_d, unlist)
net_d.df <- purrr::imap_dfr(net_d.df, ~{
  as.data.frame(as.list(unlist(.x))) %>%
    mutate(unique_ID = .y)
  })

net_d.df <- net_d.df %>%
  dplyr::select(unique_ID, everything()) %>%
  pivot_longer(cols = 2:91, names_to = "species", values_to = "d", values_drop_na = TRUE) %>%
  mutate(level = case_when(
    str_detect(species, "higher.level") ~ "higher.level",
    str_detect(species, "lower.level")  ~ "lower.level",
    .default = "ERROR"
  ))

#### higher level d' ####
higher_level_d <- net_d.df %>%
  filter(level == "higher.level") %>%
  dplyr::select(-species) %>%
  group_by(unique_ID) %>%
  summarise(mean_d_higher = mean(d))

#### higher level d' ####
lower_level_d <- net_d.df %>%
  filter(level == "lower.level") %>%
  dplyr::select(-species) %>%
  group_by(unique_ID) %>%
  summarise(mean_d_lower = mean(d))

# all together
net_d_mean <- higher_level_d %>%
  left_join(lower_level_d, by = "unique_ID") %>%
  separate(unique_ID, into = c("block", "patch"), sep = "\\.", remove = F)

# write.csv(net_d_mean, file = file.path("data", "L2_boot_metrics", "net_d_mean.csv"), row.names = FALSE)
