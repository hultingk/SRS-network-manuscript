# -------------------------------------- #
#### Analysis script for calculating species turnover and interaction rewiring ####
# Author: XXXXX
# -------------------------------------- #

# loading libraries
librarian::shelf(tidyverse, plyr, vegan, bipartite, data.table, betapart, glmmTMB, ggeffects)

source(here::here(file.path("scripts", "00_functions.R")))

# loading data 
pollinator <- read.csv(file = file.path("data", "L1_data", "SRS-plant-pollinator.csv"))

# creating unique ID
pollinator <- pollinator %>%
  mutate(unique_ID = paste(block, patch, sep = ".")) %>%
  dplyr::select(c("unique_ID", "sampling_round", "order", "family", "pollinator_species", "flower_species")) 


# subsetting interatively adding a sampling round
round1 <- pollinator %>%
  filter(sampling_round == 1) %>%
  dplyr::count(unique_ID, pollinator_species, flower_species) 

round1_2 <- pollinator %>%
  filter(sampling_round %in% c(1, 2)) %>%
  dplyr::count(unique_ID, pollinator_species, flower_species) 

round1_3 <- pollinator %>%
  filter(sampling_round %in% c(1, 2, 3)) %>%
  dplyr::count(unique_ID, pollinator_species, flower_species) 

round1_4 <- pollinator %>%
  filter(sampling_round %in% c(1, 2, 3, 4)) %>%
  dplyr::count(unique_ID, pollinator_species, flower_species) 

round1_5 <- pollinator %>%
  filter(sampling_round %in% c(1, 2, 3, 4, 5)) %>%
  dplyr::count(unique_ID, pollinator_species, flower_species) 

round1_6 <- pollinator %>%
  filter(sampling_round %in% c(1, 2, 3, 4, 5, 6)) %>%
  dplyr::count(unique_ID, pollinator_species, flower_species) 

round1_7 <- pollinator %>%
  filter(sampling_round %in% c(1, 2, 3, 4, 5, 6, 7)) %>%
  dplyr::count(unique_ID, pollinator_species, flower_species) 


# splitting by patch
round1_split <- round1 %>%
  dplyr::group_by(unique_ID) %>%
  group_split() 
round1_2_split <- round1_2 %>%
  dplyr::group_by(unique_ID) %>%
  group_split() 
round1_3_split <- round1_3 %>%
  dplyr::group_by(unique_ID) %>%
  group_split() 
round1_4_split <- round1_4 %>%
  dplyr::group_by(unique_ID) %>%
  group_split() 
round1_5_split <- round1_5 %>%
  dplyr::group_by(unique_ID) %>%
  group_split()
round1_6_split <- round1_6 %>%
  dplyr::group_by(unique_ID) %>%
  group_split()
round1_7_split <- round1_7 %>%
  dplyr::group_by(unique_ID) %>%
  group_split()

# adding patch names to each network
webs.names <- c("10.B","10.W","52.B","52.W", "53N.B", "53N.W",
                "53S.B", "53S.W", "54S.B", "54S.W", "57.B", "57.W", "8.B", "8.W")

# getting each network into correct format
round1_webs <- round1_split %>%
  lapply(prepare_matrix)
names(round1_webs) <- webs.names

round1_2_webs <- round1_2_split %>%
  lapply(prepare_matrix)
names(round1_2_webs) <- webs.names

round1_3_webs <- round1_3_split %>%
  lapply(prepare_matrix)
names(round1_3_webs) <- webs.names

round1_4_webs <- round1_4_split %>%
  lapply(prepare_matrix)
names(round1_4_webs) <- webs.names

round1_5_webs <- round1_5_split %>%
  lapply(prepare_matrix)
names(round1_5_webs) <- webs.names

round1_6_webs <- round1_6_split %>%
  lapply(prepare_matrix)
names(round1_6_webs) <- webs.names

round1_7_webs <- round1_7_split %>%
  lapply(prepare_matrix)
names(round1_7_webs) <- webs.names


## turnover and rewiring
round1_network_dissimilarity <- betalinkr_multi(webs2array(list("10.B" = as.matrix(round1_webs[[1]]), "10.W" = as.matrix(round1_webs[[2]]),
                                                         "52.B" = as.matrix(round1_webs[[3]]), "52.W" = as.matrix(round1_webs[[4]]),
                                                         "53N.B" = as.matrix(round1_webs[[5]]), "53N.W" = as.matrix(round1_webs[[6]]),
                                                         "53S.B" = as.matrix(round1_webs[[7]]), "53S.W" = as.matrix(round1_webs[[8]]),
                                                         "54S.B" = as.matrix(round1_webs[[9]]), "54S.W" = as.matrix(round1_webs[[10]]),
                                                         "57.B" = as.matrix(round1_webs[[11]]), "57.W" = as.matrix(round1_webs[[12]]),
                                                         "8.B" = as.matrix(round1_webs[[13]]), "8.W" = as.matrix(round1_webs[[14]]))), 
                                         index = "sorensen", partitioning="commondenom")
round1_network_dissimilarity <- round1_network_dissimilarity %>%
  mutate(block_pair = paste(i, j, sep = "-")) %>%
  mutate(real_pairs = case_when(
    block_pair %in% c("10.B-10.W", "52.B-52.W", "53N.B-53N.W", "53S.B-53S.W",
                      "54S.B-54S.W", "57.B-57.W", "8.B-8.W") ~ "real", 
    .default = "random"
  )) %>%
  filter(real_pairs == "real") %>%
  pivot_longer(cols = S:ST, names_to = "type", values_to = "dissimilarity") %>%
  mutate(max_sampling_round = 1)


round1_2_network_dissimilarity <- betalinkr_multi(webs2array(list("10.B" = as.matrix(round1_2_webs[[1]]), "10.W" = as.matrix(round1_2_webs[[2]]),
                                                                "52.B" = as.matrix(round1_2_webs[[3]]), "52.W" = as.matrix(round1_2_webs[[4]]),
                                                                "53N.B" = as.matrix(round1_2_webs[[5]]), "53N.W" = as.matrix(round1_2_webs[[6]]),
                                                                "53S.B" = as.matrix(round1_2_webs[[7]]), "53S.W" = as.matrix(round1_2_webs[[8]]),
                                                                "54S.B" = as.matrix(round1_2_webs[[9]]), "54S.W" = as.matrix(round1_2_webs[[10]]),
                                                                "57.B" = as.matrix(round1_2_webs[[11]]), "57.W" = as.matrix(round1_2_webs[[12]]),
                                                                "8.B" = as.matrix(round1_2_webs[[13]]), "8.W" = as.matrix(round1_2_webs[[14]]))), 
                                                index = "sorensen", partitioning="commondenom")
round1_2_network_dissimilarity <- round1_2_network_dissimilarity %>%
  mutate(block_pair = paste(i, j, sep = "-")) %>%
  mutate(real_pairs = case_when(
    block_pair %in% c("10.B-10.W", "52.B-52.W", "53N.B-53N.W", "53S.B-53S.W",
                      "54S.B-54S.W", "57.B-57.W", "8.B-8.W") ~ "real", 
    .default = "random"
  )) %>%
  filter(real_pairs == "real") %>%
  pivot_longer(cols = S:ST, names_to = "type", values_to = "dissimilarity") %>%
  mutate(max_sampling_round = 2)

round1_3_network_dissimilarity <- betalinkr_multi(webs2array(list("10.B" = as.matrix(round1_3_webs[[1]]), "10.W" = as.matrix(round1_3_webs[[2]]),
                                                                  "52.B" = as.matrix(round1_3_webs[[3]]), "52.W" = as.matrix(round1_3_webs[[4]]),
                                                                  "53N.B" = as.matrix(round1_3_webs[[5]]), "53N.W" = as.matrix(round1_3_webs[[6]]),
                                                                  "53S.B" = as.matrix(round1_3_webs[[7]]), "53S.W" = as.matrix(round1_3_webs[[8]]),
                                                                  "54S.B" = as.matrix(round1_3_webs[[9]]), "54S.W" = as.matrix(round1_3_webs[[10]]),
                                                                  "57.B" = as.matrix(round1_3_webs[[11]]), "57.W" = as.matrix(round1_3_webs[[12]]),
                                                                  "8.B" = as.matrix(round1_3_webs[[13]]), "8.W" = as.matrix(round1_3_webs[[14]]))), 
                                                  index = "sorensen", partitioning="commondenom")
round1_3_network_dissimilarity <- round1_3_network_dissimilarity %>%
  mutate(block_pair = paste(i, j, sep = "-")) %>%
  mutate(real_pairs = case_when(
    block_pair %in% c("10.B-10.W", "52.B-52.W", "53N.B-53N.W", "53S.B-53S.W",
                      "54S.B-54S.W", "57.B-57.W", "8.B-8.W") ~ "real", 
    .default = "random"
  )) %>%
  filter(real_pairs == "real") %>%
  pivot_longer(cols = S:ST, names_to = "type", values_to = "dissimilarity") %>%
  mutate(max_sampling_round = 3)



round1_4_network_dissimilarity <- betalinkr_multi(webs2array(list("10.B" = as.matrix(round1_4_webs[[1]]), "10.W" = as.matrix(round1_4_webs[[2]]),
                                                                  "52.B" = as.matrix(round1_4_webs[[3]]), "52.W" = as.matrix(round1_4_webs[[4]]),
                                                                  "53N.B" = as.matrix(round1_4_webs[[5]]), "53N.W" = as.matrix(round1_4_webs[[6]]),
                                                                  "53S.B" = as.matrix(round1_4_webs[[7]]), "53S.W" = as.matrix(round1_4_webs[[8]]),
                                                                  "54S.B" = as.matrix(round1_4_webs[[9]]), "54S.W" = as.matrix(round1_4_webs[[10]]),
                                                                  "57.B" = as.matrix(round1_4_webs[[11]]), "57.W" = as.matrix(round1_4_webs[[12]]),
                                                                  "8.B" = as.matrix(round1_4_webs[[13]]), "8.W" = as.matrix(round1_4_webs[[14]]))), 
                                                  index = "sorensen", partitioning="commondenom")
round1_4_network_dissimilarity <- round1_4_network_dissimilarity %>%
  mutate(block_pair = paste(i, j, sep = "-")) %>%
  mutate(real_pairs = case_when(
    block_pair %in% c("10.B-10.W", "52.B-52.W", "53N.B-53N.W", "53S.B-53S.W",
                      "54S.B-54S.W", "57.B-57.W", "8.B-8.W") ~ "real", 
    .default = "random"
  )) %>%
  filter(real_pairs == "real") %>%
  pivot_longer(cols = S:ST, names_to = "type", values_to = "dissimilarity") %>%
  mutate(max_sampling_round = 4)



round1_5_network_dissimilarity <- betalinkr_multi(webs2array(list("10.B" = as.matrix(round1_5_webs[[1]]), "10.W" = as.matrix(round1_5_webs[[2]]),
                                                                  "52.B" = as.matrix(round1_5_webs[[3]]), "52.W" = as.matrix(round1_5_webs[[4]]),
                                                                  "53N.B" = as.matrix(round1_5_webs[[5]]), "53N.W" = as.matrix(round1_5_webs[[6]]),
                                                                  "53S.B" = as.matrix(round1_5_webs[[7]]), "53S.W" = as.matrix(round1_5_webs[[8]]),
                                                                  "54S.B" = as.matrix(round1_5_webs[[9]]), "54S.W" = as.matrix(round1_5_webs[[10]]),
                                                                  "57.B" = as.matrix(round1_5_webs[[11]]), "57.W" = as.matrix(round1_5_webs[[12]]),
                                                                  "8.B" = as.matrix(round1_5_webs[[13]]), "8.W" = as.matrix(round1_5_webs[[14]]))), 
                                                  index = "sorensen", partitioning="commondenom")
round1_5_network_dissimilarity <- round1_5_network_dissimilarity %>%
  mutate(block_pair = paste(i, j, sep = "-")) %>%
  mutate(real_pairs = case_when(
    block_pair %in% c("10.B-10.W", "52.B-52.W", "53N.B-53N.W", "53S.B-53S.W",
                      "54S.B-54S.W", "57.B-57.W", "8.B-8.W") ~ "real", 
    .default = "random"
  )) %>%
  filter(real_pairs == "real") %>%
  pivot_longer(cols = S:ST, names_to = "type", values_to = "dissimilarity") %>%
  mutate(max_sampling_round = 5)

round1_6_network_dissimilarity <- betalinkr_multi(webs2array(list("10.B" = as.matrix(round1_6_webs[[1]]), "10.W" = as.matrix(round1_6_webs[[2]]),
                                                                  "52.B" = as.matrix(round1_6_webs[[3]]), "52.W" = as.matrix(round1_6_webs[[4]]),
                                                                  "53N.B" = as.matrix(round1_6_webs[[5]]), "53N.W" = as.matrix(round1_6_webs[[6]]),
                                                                  "53S.B" = as.matrix(round1_6_webs[[7]]), "53S.W" = as.matrix(round1_6_webs[[8]]),
                                                                  "54S.B" = as.matrix(round1_6_webs[[9]]), "54S.W" = as.matrix(round1_6_webs[[10]]),
                                                                  "57.B" = as.matrix(round1_6_webs[[11]]), "57.W" = as.matrix(round1_6_webs[[12]]),
                                                                  "8.B" = as.matrix(round1_6_webs[[13]]), "8.W" = as.matrix(round1_6_webs[[14]]))), 
                                                  index = "sorensen", partitioning="commondenom")
round1_6_network_dissimilarity <- round1_6_network_dissimilarity %>%
  mutate(block_pair = paste(i, j, sep = "-")) %>%
  mutate(real_pairs = case_when(
    block_pair %in% c("10.B-10.W", "52.B-52.W", "53N.B-53N.W", "53S.B-53S.W",
                      "54S.B-54S.W", "57.B-57.W", "8.B-8.W") ~ "real", 
    .default = "random"
  )) %>%
  filter(real_pairs == "real") %>%
  pivot_longer(cols = S:ST, names_to = "type", values_to = "dissimilarity") %>%
  mutate(max_sampling_round = 6)

round1_7_network_dissimilarity <- betalinkr_multi(webs2array(list("10.B" = as.matrix(round1_7_webs[[1]]), "10.W" = as.matrix(round1_7_webs[[2]]),
                                                                  "52.B" = as.matrix(round1_7_webs[[3]]), "52.W" = as.matrix(round1_7_webs[[4]]),
                                                                  "53N.B" = as.matrix(round1_7_webs[[5]]), "53N.W" = as.matrix(round1_7_webs[[6]]),
                                                                  "53S.B" = as.matrix(round1_7_webs[[7]]), "53S.W" = as.matrix(round1_7_webs[[8]]),
                                                                  "54S.B" = as.matrix(round1_7_webs[[9]]), "54S.W" = as.matrix(round1_7_webs[[10]]),
                                                                  "57.B" = as.matrix(round1_7_webs[[11]]), "57.W" = as.matrix(round1_7_webs[[12]]),
                                                                  "8.B" = as.matrix(round1_7_webs[[13]]), "8.W" = as.matrix(round1_7_webs[[14]]))), 
                                                  index = "sorensen", partitioning="commondenom")
round1_7_network_dissimilarity <- round1_7_network_dissimilarity %>%
  mutate(block_pair = paste(i, j, sep = "-")) %>%
  mutate(real_pairs = case_when(
    block_pair %in% c("10.B-10.W", "52.B-52.W", "53N.B-53N.W", "53S.B-53S.W",
                      "54S.B-54S.W", "57.B-57.W", "8.B-8.W") ~ "real", 
    .default = "random"
  )) %>%
  filter(real_pairs == "real") %>%
  pivot_longer(cols = S:ST, names_to = "type", values_to = "dissimilarity") %>%
  mutate(max_sampling_round = 7)


## all together
network_dissimilarity <- rbind(round1_network_dissimilarity, round1_2_network_dissimilarity,
                               round1_3_network_dissimilarity, round1_4_network_dissimilarity,
                               round1_5_network_dissimilarity, round1_6_network_dissimilarity,
                               round1_7_network_dissimilarity)
network_dissimilarity$type <- factor(network_dissimilarity$type, levels=c("S", "WN", "ST", "OS"))

# plotting over sampling round additions
dissimilarity_rounds <- network_dissimilarity %>%
  filter(type != "S") %>%
  ggplot() +
  geom_jitter(aes(max_sampling_round, dissimilarity, color = type), size = 5, alpha = 0.55, width = 0.10, height = 0) +
  geom_smooth(aes(max_sampling_round, dissimilarity, fill = type, color = type)) +
  theme_classic(base_size = 22) +
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        axis.text = element_text(size = 22),
        axis.ticks = element_line(color = "black", linewidth = 0.7),
        strip.text.x = element_text(hjust = -0.05),
        panel.background = element_rect(fill = "transparent", color = NA), # Inside axes
        plot.background = element_rect(fill = "transparent", color = NA))  +
  scale_color_manual(values=c("#C2697FFF","#E68E54FF", "#DCB254"), labels = c(expression(beta[WN]), expression(beta[ST]), expression(beta[OS])), name = "Component") +
  scale_fill_manual(values=c("#C2697FFF","#E68E54FF", "#DCB254"), labels = c(expression(beta[WN]), expression(beta[ST]), expression(beta[OS])), name = "Component") +
  xlab("Number of sampling rounds included") +
  theme(legend.position = "right") + 
  ylab(expression(atop("Dissimilarity between", paste("patch pairs"))))
dissimilarity_rounds
# exporting
# pdf(file = file.path("plots", "dissimilarity_rounds.pdf"), width = 12.5, height = 5.5)
# dissimilarity_rounds
# dev.off()


#### all rounds together ####
# filtering for just interaction rewiring and species turnover, filtering for pairs within a block
round1_7_network_dissimilarity <- round1_7_network_dissimilarity %>%
  filter(!type %in% c("S"))

# filtering just species turnover and iteraction rewiring
round1_7_network_dissimilarity_filtered <- round1_7_network_dissimilarity %>%
  filter(!type %in% c("WN")) %>%
  pivot_wider(names_from = type, values_from = dissimilarity)
# t test
m.round1_7_network_dissimilarity <- t.test(round1_7_network_dissimilarity_filtered$OS, round1_7_network_dissimilarity_filtered$ST, paired = TRUE, data = round1_7_network_dissimilarity_filtered)
m.round1_7_network_dissimilarity

# model predictions for plotting
m.round1_7_network_dissimilarity_plot <- glmmTMB(dissimilarity ~ type,
                                            data = round1_7_network_dissimilarity)
round1_7_network_dissimilarity.predict <- ggpredict(m.round1_7_network_dissimilarity_plot, terms = c("type"), back_transform = TRUE)

round1_7_network_dissimilarity$type <- factor(round1_7_network_dissimilarity$type, levels=c("WN", "ST", "OS"))
round1_7_network_dissimilarity.predict$x <- factor(round1_7_network_dissimilarity.predict$x, levels=c("WN", "ST", "OS"))

# plotting
dissimilarity_plot <- round1_7_network_dissimilarity %>%
  ggplot() +
  geom_jitter(aes(x = type, y = dissimilarity, color = type), data = round1_7_network_dissimilarity, size = 5, alpha = 0.55,
              width = 0.08, height = 0) +
  geom_errorbar(aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high, fill = x), color = "black",
                data = round1_7_network_dissimilarity.predict, width = 0, linewidth = 2.5) +
  #geom_line(aes(x = x, y = predicted, group = group), data = network_dissimilarity.predict, linewidth = 2, linetype = 1) +
  geom_point(aes(x = x, y = predicted, fill = x), data = round1_7_network_dissimilarity.predict, 
             size = 6, colour="black", pch=21, stroke = 2) +
  scale_color_manual(values=c("#C2697FFF","#E68E54FF", "#DCB254")) +
  scale_fill_manual(values=c("#C2697FFF","#E68E54FF", "#DCB254")) +
  scale_x_discrete(labels = c(expression(beta[WN]), expression(beta[ST]), expression(beta[OS]))) +
  theme_classic(base_size = 22) +
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        axis.text = element_text(size = 22),
        axis.ticks = element_line(color = "black", linewidth = 0.7),
        strip.text.x = element_text(hjust = -0.05),
        panel.background = element_rect(fill = "transparent", color = NA), # Inside axes
        plot.background = element_rect(fill = "transparent", color = NA))  +
  xlab(NULL) +
  theme(legend.position = "none") + 
  ylab(expression(atop("Dissimilarity between", paste("patch pairs"))))
dissimilarity_plot

# # exporting
# pdf(file = file.path("plots", "dissimilarity_plot.pdf"), width = 6, height = 5)
# dissimilarity_plot
# dev.off()


# exporting all together
si_turnover_plot <- cowplot::plot_grid(dissimilarity_plot, dissimilarity_rounds,
                                       labels = c("(A)", "(B)"), rel_widths = c(1, 1.5),
                                       label_size = 26, label_x = c(0.80, 0.57), label_y = 0.95)
si_turnover_plot

# pdf(file = file.path("plots", "si_turnover_plot.pdf"), width = 13, height = 4.6)
# si_turnover_plot
# dev.off()







