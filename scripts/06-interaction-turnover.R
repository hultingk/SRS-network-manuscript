# -------------------------------------- #
#### Analysis script for calculating species turnover and interaction rewiring ####
# Author: XXXXX
# -------------------------------------- #

# loading libraries
librarian::shelf(tidyverse, plyr, vegan, bipartite, data.table, betapart, glmmTMB, ggeffects)

source(here::here(file.path("scripts", "00-functions.R")))

# loading data 
pollinator <- read.csv(file = file.path("data", "L1_data", "SRS-plant-pollinator.csv"))

# creating unique ID
pollinator <- pollinator %>%
  mutate(unique_ID = paste(block, patch, sep = ".")) %>%
  dplyr::select(c("unique_ID", "sampling_round", "order", "family", "pollinator_species", "flower_species")) 


round1_7 <- pollinator %>%
  filter(sampling_round %in% c(1, 2, 3, 4, 5, 6, 7)) %>%
  dplyr::count(unique_ID, pollinator_species, flower_species) 

# splitting by patch
round1_7_split <- round1_7 %>%
  dplyr::group_by(unique_ID) %>%
  group_split()

# adding patch names to each network
webs.names <- c("10.B","10.W","52.B","52.W", "53N.B", "53N.W",
                "53S.B", "53S.W", "54S.B", "54S.W", "57.B", "57.W", "8.B", "8.W")

# getting each network into correct format
round1_7_webs <- round1_7_split %>%
  lapply(prepare_matrix)
names(round1_7_webs) <- webs.names


## turnover and rewiring
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








