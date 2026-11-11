# -------------------------------------- #
#### Diversity metrics excluding Apis ####
# Author: XXXX
# -------------------------------------- #


# loading libraries
librarian::shelf(tidyverse, vegan, glmmTMB, iNEXT, DHARMa)

# loading data 
pollinator <- read.csv(file = file.path("data", "L1_data", "SRS-plant-pollinator.csv")) # plant-pollinator interaction data


# creating unique ID
pollinator <- pollinator %>%
  mutate(unique_ID = paste(block, patch, sep = ".")) %>%
  dplyr::select(c("unique_ID", "order", "family", "pollinator_species", "flower_species")) 

# No Apis - floral Hill numbers #
floral_noApis_inext <- pollinator %>%
  filter(!pollinator_species %in% c("Apis mellifera", "Poecilognathus sulphureus")) %>%
  dplyr::count(unique_ID, flower_species) %>% # counting # of each floral species per patch
  pivot_wider(names_from = unique_ID, values_from = n, values_fill = 0) %>% # pivoting into wider format for flower species by patch matrix
  column_to_rownames("flower_species")

floral_noApis_hill <- iNEXT(floral_noApis_inext, q = c(0,1,2), datatype = "abundance") # hill numbers
floral_noApis_est <- floral_noApis_hill$iNextEst[["coverage_based"]] # getting coveraged based estimates
floral_noApis_target_cov <- floral_noApis_est %>% # figuring out what is the lowest max sampling coverage - will use estimates from that coverage for all patches
  filter(Method == "Observed") %>%
  group_by(Assemblage) %>%
  summarise(max_cov = max(SC)) %>%
  summarise(min(max_cov)) %>%
  pull()

floral_noApis_est <- floral_noApis_est %>% # getting estimates as close to the target sampling coverage as possible
  group_by(Assemblage, Order.q) %>%
  slice_min(abs(SC - floral_noApis_target_cov), n = 1) %>%
  ungroup()

floral_noApis_est <- floral_noApis_est %>% # wrangling into format
  separate(Assemblage, into = c("block", "patch"), remove = FALSE) %>%
  mutate(Order.q = as.factor(Order.q)) %>%
  select(Assemblage, block, patch, Order.q, qD) %>%
  pivot_wider(names_from = Order.q, values_from = qD) %>%
  dplyr::rename(unique_ID = Assemblage, floral_noApis_0 = `0`, floral_noApis_1 = `1`, floral_noApis_2 = `2`)

# floral no Apis models
m.floral_0_noApis <- glmmTMB(floral_noApis_0 ~ patch + (1|block), 
                             data = floral_noApis_est)
summary(m.floral_0_noApis)
plot(simulateResiduals(m.floral_0_noApis))

m.floral_1_noApis <- glmmTMB(floral_noApis_1 ~ patch + (1|block), 
                             data = floral_noApis_est)
summary(m.floral_1_noApis)
plot(simulateResiduals(m.floral_1_noApis))


# No Apis - pollinator Hill numbers #
pollinator_noApis_inext <- pollinator %>%
  filter(!pollinator_species %in% c("Apis mellifera", "Poecilognathus sulphureus")) %>%
  dplyr::count(unique_ID, pollinator_species) %>% # counting # of each pollinator species per patch
  pivot_wider(names_from = unique_ID, values_from = n, values_fill = 0) %>% # pivoting into wider format for flower species by patch matrix
  column_to_rownames("pollinator_species")

pollinator_noApis_hill <- iNEXT(pollinator_noApis_inext, q = c(0,1,2), datatype = "abundance") # hill numbers
pollinator_noApis_est <- pollinator_noApis_hill$iNextEst[["coverage_based"]] # getting coveraged based estimates
pollinator_noApis_target_cov <- pollinator_noApis_est %>% # figuring out what is the lowest max sampling coverage - will use estimates from that coverage for all patches
  filter(Method == "Observed") %>%
  group_by(Assemblage) %>%
  summarise(max_cov = max(SC)) %>%
  summarise(min(max_cov)) %>%
  pull()

pollinator_noApis_est <- pollinator_noApis_est %>% # getting estimates as close to the target sampling coverage as possible
  group_by(Assemblage, Order.q) %>%
  slice_min(abs(SC - pollinator_noApis_target_cov), n = 1) %>%
  ungroup()

pollinator_noApis_est <- pollinator_noApis_est %>% # wrangling into format
  separate(Assemblage, into = c("block", "patch"), remove = FALSE) %>%
  mutate(Order.q = as.factor(Order.q)) %>%
  select(Assemblage, block, patch, Order.q, qD) %>%
  pivot_wider(names_from = Order.q, values_from = qD) %>%
  dplyr::rename(unique_ID = Assemblage, pollinator_noApis_0 = `0`, pollinator_noApis_1 = `1`, pollinator_noApis_2 = `2`)


# pollinator no Apis models
m.pollinator_0_noApis <- glmmTMB(pollinator_noApis_0 ~ patch + (1|block), 
                             data = pollinator_noApis_est)
summary(m.pollinator_0_noApis)
plot(simulateResiduals(m.pollinator_0_noApis))

m.pollinator_1_noApis <- glmmTMB(pollinator_noApis_1 ~ patch + (1|block), 
                             data = pollinator_noApis_est)
summary(m.pollinator_1_noApis)
plot(simulateResiduals(m.pollinator_1_noApis))




#### plotting diversity ####
# floral plot
# model predictions for plotting
m.floral_1_noApis.df <- ggpredict(m.floral_1_noApis, terms = c("patch"), back_transform = TRUE)
# plotting
floral_1_noApis.pred <- m.floral_1_noApis.df %>%
  ggplot() +
  geom_jitter(aes(x = patch, y = floral_noApis_1, color = patch), data = floral_noApis_est, size = 6, alpha = 0.55,
              width = 0.08, height = 0) +
  geom_errorbar(aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high, fill = x), color = "black",
                data = m.floral_1_noApis.df, width = 0, linewidth = 2.5) +
  geom_line(aes(x = x, y = predicted, group = group), linewidth = 2, linetype = 2) +
  geom_point(aes(x = x, y = predicted, fill = x), size = 8, colour="black", pch=21, stroke = 2) +
  scale_x_discrete(labels = c('Connected', 'Unconnected')) +
  theme_classic(base_size = 26) +
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        # panel.grid.major = element_line(linetype = 2, linewidth = 0.7, color = "grey85"), 
        panel.grid.minor = element_blank(),
        axis.ticks = element_line(color = "black", linewidth = 0.7),
        strip.text.x = element_text(hjust = -0.05),
        panel.background = element_rect(fill = "transparent", color = NA), # Inside axes
        plot.background = element_rect(fill = "transparent", color = NA)) +
  scale_color_manual(values=c("#5B859EFF","#DCB254")) +
  scale_fill_manual(values=c("#5B859EFF","#DCB254")) +
  xlab("Patch type") +
  ylab(expression(atop("Plant diversity excluding two", paste("dominant pollinators (q = 1)")))) +
  theme(legend.position = "none") 
floral_1_noApis.pred


# pollinator plot
# model predictions for plotting
m.pollinator_1_noApis.df <- ggpredict(m.pollinator_1_noApis, terms = c("patch"), back_transform = TRUE)
# plotting
pollinator_1_noApis.pred <- m.pollinator_1_noApis.df %>%
  ggplot() +
  geom_jitter(aes(x = patch, y = pollinator_noApis_1, color = patch), data = pollinator_noApis_est, size = 6, alpha = 0.55,
              width = 0.08, height = 0) +
  geom_errorbar(aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high, fill = x), color = "black",
                data = m.pollinator_1_noApis.df, width = 0, linewidth = 2.5) +
  geom_line(aes(x = x, y = predicted, group = group), linewidth = 2, linetype = 2) +
  geom_point(aes(x = x, y = predicted, fill = x), size = 8, colour="black", pch=21, stroke = 2) +
  scale_x_discrete(labels = c('Connected', 'Unconnected')) +
  theme_classic(base_size = 26) +
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        # panel.grid.major = element_line(linetype = 2, linewidth = 0.7, color = "grey85"), 
        panel.grid.minor = element_blank(),
        axis.ticks = element_line(color = "black", linewidth = 0.7),
        strip.text.x = element_text(hjust = -0.05),
        panel.background = element_rect(fill = "transparent", color = NA), # Inside axes
        plot.background = element_rect(fill = "transparent", color = NA)) +
  scale_color_manual(values=c("#5B859EFF","#DCB254")) +
  scale_fill_manual(values=c("#5B859EFF","#DCB254")) +
  xlab("Patch type") +
  ylab(expression(atop("Pollinator diversity excluding two", paste("dominant pollinators (q = 1)")))) +
  theme(legend.position = "none") 
pollinator_1_noApis.pred



si_noApis_diversity_plot <- cowplot::plot_grid(floral_1_noApis.pred, pollinator_1_noApis.pred, 
                                     labels = c("A)", "B)"), label_size = 26, label_x = 0.25, label_y = 0.95)
si_noApis_diversity_plot

# pdf(file = file.path("plots", "diversity-plots", "si_noApis_diversity_plot.pdf"), width = 13, height = 6)
# si_noApis_diversity_plot
# dev.off()



