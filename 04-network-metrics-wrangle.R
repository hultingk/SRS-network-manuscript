# -------------------------------------- #
#### Script for wrangling output of network calculations ####
# Author: XXXXX
# -------------------------------------- #

librarian::shelf(tidyverse, plyr, data.table)

# reading in files
lst_h2 <- readRDS(file = file.path("data", "L2_boot_metrics", "h2.RData"))
lst_links <- readRDS(file = file.path("data", "L2_boot_metrics", "links.RData"))
lst_niche <- readRDS(file = file.path("data", "L2_boot_metrics", "niche.RData"))
lst_nodf <- readRDS(file = file.path("data", "L2_boot_metrics", "nodf.RData"))
lst_robustness <- readRDS(file = file.path("data", "L2_boot_metrics", "robustness.RData"))
lst_shannon <- readRDS(file = file.path("data", "L2_boot_metrics", "shannon.RData"))
net_d_mean <- read.csv(file = file.path("data", "L2_boot_metrics", "net_d_mean.csv"))

##### H2 ####
h2_boot <- lst_h2[["H2"]][["stats_df"]]
h2_boot <- h2_boot %>%
  filter(spl_size == 158) %>% # subsetting at 158 interactions to compare networks
  separate(web, into = c("block", "patch")) %>%
  dplyr::rename(h2 = mean) %>%
  dplyr::select(block, patch, h2)

##### links ####
links_boot <- lst_links[["links per species"]][["stats_df"]]
links_boot <- links_boot %>%
  filter(spl_size == 158) %>% # subsetting at 158 interactions to compare networks
  separate(web, into = c("block", "patch")) %>%
  dplyr::rename(links_per_sp = mean) %>%
  dplyr::select(block, patch, links_per_sp)

##### niche ####
HL_niche_boot <- lst_niche[["niche.overlap.HL"]][["stats_df"]]
LL_niche_boot <- lst_niche[["niche.overlap.LL"]][["stats_df"]]

HL_niche_boot <- HL_niche_boot %>%
  filter(spl_size == 158) %>% # subsetting at 158 interactions to compare networks
  separate(web, into = c("block", "patch")) %>%
  dplyr::rename(HL_niche = mean) %>%
  dplyr::select(block, patch, HL_niche)

LL_niche_boot <- LL_niche_boot %>%
  filter(spl_size == 158) %>% # subsetting at 158 interactions to compare networks
  separate(web, into = c("block", "patch")) %>%
  dplyr::rename(LL_niche = mean) %>%
  dplyr::select(block, patch, LL_niche)


#### nodf ####
nodf_boot <- lst_nodf[["NODF"]][["stats_df"]]
nodf_boot <- nodf_boot %>%
  filter(spl_size == 158) %>% # subsetting at 158 interactions to compare networks
  separate(web, into = c("block", "patch")) %>%
  dplyr::rename(nodf = mean) %>%
  dplyr::select(block, patch, nodf)

#### robustness ####
HL_robustness_boot <- lst_robustness[["robustness.HL"]][["stats_df"]]
LL_robustness_boot <- lst_robustness[["robustness.LL"]][["stats_df"]]

HL_robustness_boot <- HL_robustness_boot %>%
  filter(spl_size == 158) %>% # subsetting at 158 interactions to compare networks
  separate(web, into = c("block", "patch")) %>%
  dplyr::rename(HL_robustness = mean) %>%
  dplyr::select(block, patch, HL_robustness)

LL_robustness_boot <- LL_robustness_boot %>%
  filter(spl_size == 158) %>% # subsetting at 158 interactions to compare networks
  separate(web, into = c("block", "patch")) %>%
  dplyr::rename(LL_robustness = mean) %>%
  dplyr::select(block, patch, LL_robustness)


##### shannon diversity of interactions ####
shannon_boot <- lst_shannon[["Shannon diversity"]][["stats_df"]]
shannon_boot <- shannon_boot %>%
  filter(spl_size == 158) %>% # subsetting at 158 interactions to compare networks
  separate(web, into = c("block", "patch")) %>%
  dplyr::rename(shannon = mean) %>%
  dplyr::select(block, patch, shannon)





#### All together ####
network_metrics_boot <- nodf_boot %>%
  left_join(h2_boot, by = c("block", "patch")) %>%
  left_join(links_boot, by = c("block", "patch")) %>%
  left_join(HL_niche_boot, by = c("block", "patch")) %>%
  left_join(LL_niche_boot, by = c("block", "patch")) %>%
  left_join(HL_robustness_boot, by = c("block", "patch")) %>%
  left_join(LL_robustness_boot, by = c("block", "patch")) %>%
  left_join(shannon_boot, by = c("block", "patch")) %>%
  left_join(net_d_mean, by = c("block", "patch"))


# write.csv(network_metrics_boot, file = file.path("data", "L2_boot_metrics", "network_metrics_boot.csv"), row.names = FALSE)

