# -------------------------------------- #
#### Analysis script for network metrics ####
# Author: XXXXX
# -------------------------------------- #

# loading libraries
librarian::shelf(glmmTMB, DHARMa, tidyverse, ggeffects, car, kableExtra)

# loading data
network_metrics_boot <- read.csv(file = file.path("data", "L2_boot_metrics", "network_metrics_boot.csv"))


#### models ####
# nestedness 
m_nodf <- glmmTMB(nodf ~ patch + (1|block),
                  data = network_metrics_boot)
summary(m_nodf)
plot(simulateResiduals(m_nodf))
anova.nodf <- Anova(m_nodf, type = "III")


# links per species 
m_links <- glmmTMB(links_per_sp ~ patch + (1|block),
                   data = network_metrics_boot)
summary(m_links)
plot(simulateResiduals(m_links))
anova.links <- Anova(m_links, type = "III")
-0.10589/1.28683 * 100 # 8.228748% decrease
confint(m_links)
-0.2243528/1.2030633 * 100 # 18.64846% decrease
0.01256507/1.37058953 * 100 # 0.9167639% decrease


# shannon diversity 
m_shannon <- glmmTMB(shannon ~ patch + (1|block),
                     data = network_metrics_boot)
summary(m_shannon)
plot(simulateResiduals(m_shannon))
anova.shannon <- Anova(m_shannon, type = "III")
-0.8589/4.0814 * 100 # 21.04425% decrease
confint(m_shannon)
-1.402504/3.689878 * 100 # 38.0095% decrease
-0.3153479/4.4729143 * 100 # 7.050166% decrease


# d' higher level 
m.d_higher <- glmmTMB(mean_d_higher ~ patch + (1|block),
                      data = network_metrics_boot)
summary(m.d_higher)
anova.d_higher <- Anova(m.d_higher, type = "III")


# d' lower level 
m.d_lower <- glmmTMB(mean_d_lower ~ patch + (1|block),
                     data = network_metrics_boot)
summary(m.d_lower)
anova.d_lower <- Anova(m.d_lower, type = "III")
0.1106/0.5298 * 100 # 20.8758% decrease
confint(m.d_lower)
0.04577496/0.48390382 * 100 # 9.459516% decrease
0.1755210/0.5756481 * 100 # 30.49102% decrease


##### tables #####
anova.shannon.df <- as.data.frame(anova.shannon)
anova.shannon.df <- anova.shannon.df %>%
  rownames_to_column("model_term") %>%
  mutate(`Response variable` = "Shannon diversity of interactions")

anova.links.df <- as.data.frame(anova.links)
anova.links.df <- anova.links.df %>%
  rownames_to_column("model_term") %>%
  mutate(`Response variable` = "Links per species")

anova.nodf.df <- as.data.frame(anova.nodf)
anova.nodf.df <- anova.nodf.df %>%
  rownames_to_column("model_term") %>%
  mutate(`Response variable` = "Nestedness (NODF)")

anova.d_lower.df <- as.data.frame(anova.d_lower)
anova.d_lower.df <- anova.d_lower.df %>%
  rownames_to_column("model_term") %>%
  mutate(`Response variable` = "Plant d'")

anova.d_higher.df <- as.data.frame(anova.d_higher)
anova.d_higher.df <- anova.d_higher.df %>%
  rownames_to_column("model_term") %>%
  mutate(`Response variable` = "Pollinator d'")

m.net_anova_all <- rbind(
  anova.shannon.df, anova.links.df, anova.nodf.df, anova.d_lower.df, anova.d_higher.df
)


rename_variable_anova <- tibble(model_term = c("(Intercept)", "patch"),
                                Predictor = c("Intercept", "Patch Type"))

m.net_anova_all <- m.net_anova_all %>%
  left_join(rename_variable_anova, by = "model_term") %>%
  dplyr::select(`Response variable`, Predictor, Chisq, Df, `Pr(>Chisq)`) %>%
  dplyr::rename(p.value = `Pr(>Chisq)`, df = Df)

table_net_anova <- m.net_anova_all %>%
  kbl(digits = 3) %>%
  kable_classic(full_width = T) %>%
  kable_styling(html_font = "Times New Roman",
                font_size = 16) %>%
  collapse_rows(columns = c(1), target = 1) %>%
  row_spec(0, extra_css = "border-bottom: 5px double;") %>%
  row_spec(1:nrow(m.net_anova_all), extra_css = "border-bottom: 1px solid;") %>%
  row_spec(0:nrow(m.net_anova_all), extra_css = "padding-bottom: 5px;")
table_net_anova

# exporting
#save_kable(table_net_anova, file = file.path("tables", "table_net_anova.html"))





##### plotting #####
# nestedness plot
# plotting nodf per species
m.nodf.df <- ggpredict(m_nodf, terms = c("patch"), back_transform = TRUE)
# plotting
nodf.pred <- m.nodf.df %>%
  ggplot() +
  geom_jitter(aes(x = patch, y = nodf, color = patch), data = network_metrics_boot, size = 6, alpha = 0.55,
              width = 0.08, height = 0) +
  geom_errorbar(aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high, fill = x), color = "black",
                data = m.nodf.df, width = 0, linewidth = 2.5) +
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
  ylab(expression("Nestedness (NODF)")) +
  theme(legend.position = "none") 
nodf.pred

# pdf(file = file.path("plots", "network-plots", "nodf.pred.pdf"), width = 7, height = 6)
# nodf.pred
# dev.off()



# plotting links per species
m.links.df <- ggpredict(m_links, terms = c("patch"), back_transform = TRUE)
# plotting
links.pred <- m.links.df %>%
  ggplot() +
  geom_jitter(aes(x = patch, y = links_per_sp, color = patch), data = network_metrics_boot, size = 6, alpha = 0.55,
              width = 0.08, height = 0) +
  geom_errorbar(aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high, fill = x), color = "black",
                data = m.links.df, width = 0, linewidth = 2.5) +
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
  ylab(expression("Links per species")) +
  theme(legend.position = "none") 
links.pred

# pdf(file = file.path("plots", "network-plots", "links.pred.pdf"), width = 7, height = 6)
# links.pred
# dev.off()


# plotting shannon diversity of interactions
m.shannon.df <- ggpredict(m_shannon, terms = c("patch"), back_transform = TRUE)
# plotting
shannon.pred <- m.shannon.df %>%
  ggplot() +
  geom_jitter(aes(x = patch, y = shannon, color = patch), data = network_metrics_boot, size = 6, alpha = 0.55,
              width = 0.08, height = 0) +
  geom_errorbar(aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high, fill = x), color = "black",
                data = m.shannon.df, width = 0, linewidth = 2.5) +
  geom_line(aes(x = x, y = predicted, group = group), linewidth = 2, linetype = 1) +
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
  ylab(expression(atop("Shannon diversity", paste("of interactions")))) +
  theme(legend.position = "none") 
shannon.pred

# pdf(file = file.path("plots", "network-plots", "shannon.pred.pdf"), width = 7, height = 6)
# shannon.pred
# dev.off()

# plotting
m.d_higher.df <- ggpredict(m.d_higher, terms = c("patch"), back_transform = TRUE)
# plotting
d_higher.pred <- m.d_higher.df %>%
  ggplot() +
  geom_jitter(aes(x = patch, y = mean_d_higher, color = patch), data = network_metrics_boot, size = 6, alpha = 0.55,
              width = 0.08, height = 0) +
  geom_errorbar(aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high, fill = x), color = "black",
                data = m.d_higher.df, width = 0, linewidth = 2.5) +
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
  ylab(expression("Pollinator d'")) +
  theme(legend.position = "none") 
d_higher.pred

# pdf(file = file.path("plots", "network-plots", "d_higher.pred.pdf"), width = 7, height = 6)
# d_higher.pred
# dev.off()



# plotting
m.d_lower.df <- ggpredict(m.d_lower, terms = c("patch"), back_transform = TRUE)
# plotting
d_lower.pred <- m.d_lower.df %>%
  ggplot() +
  geom_jitter(aes(x = patch, y = mean_d_lower, color = patch), data = network_metrics_boot, size = 6, alpha = 0.55,
              width = 0.08, height = 0) +
  geom_errorbar(aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high, fill = x), color = "black",
                data = m.d_lower.df, width = 0, linewidth = 2.5) +
  geom_line(aes(x = x, y = predicted, group = group), linewidth = 2, linetype = 1) +
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
  #scale_y_continuous(limits = c(0.45, 0.75), labels = label_number(accuracy = 0.1)) +
  xlab("Patch type") +
  #ylab(expression("Plant d'")) +
  ylab(expression(atop(" ", paste("Plant d'")))) +
  theme(legend.position = "none") 
d_lower.pred

# pdf(file = file.path("plots", "network-plots", "d_lower.pred.pdf"), width = 7, height = 6)
# d_lower.pred
# dev.off()



# all together
network_plot <- cowplot::plot_grid(shannon.pred, links.pred, d_lower.pred, d_higher.pred,
                                   labels = c("A)", "B)", "C)", "D)"), rel_widths = c(1.07, 1),
                                   label_size = 26, label_x = c(0.25, 0.20, 0.23, 0.19), label_y = 0.95)
network_plot

# pdf(file = file.path("plots", "network-plots", "network_plot.pdf"), width = 13, height = 11)
# network_plot
# dev.off()


