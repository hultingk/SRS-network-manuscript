# -------------------------------------- #
#### Analysis script for diversity metrics ####
# Author: XXXXX
# -------------------------------------- #

# loading libraries 
librarian::shelf(tidyverse, glmmTMB, DHARMa, performance, ggeffects, ggpubr, car, kableExtra)


# load data
diversity_metrics <- read.csv(file = file.path("data", "L3_metrics", "diversity_metrics.csv"))


#### abundance ####
# # abundance all interactions
m.abundance <- glmmTMB(abundance ~ patch + (1|block), 
                       data = diversity_metrics,
                       family = "nbinom2")
summary(m.abundance)
plot(simulateResiduals(m.abundance))
anova.abundance <- Anova(m.abundance, type = "III")


##### floral diversity #####
m.floral_0 <- glmmTMB(floral_0 ~ patch + (1|block), 
                      data = diversity_metrics)
summary(m.floral_0)
plot(simulateResiduals(m.floral_0))
anova.floral_0 <- Anova(m.floral_0, type = "III")

# shannon all species
m.floral_1 <- glmmTMB(floral_1 ~ patch + (1|block), 
                      data = diversity_metrics)
summary(m.floral_1)
plot(simulateResiduals(m.floral_1))
anova.floral_1 <- Anova(m.floral_1, type = "III")
confint(m.floral_1)



##### pollinator diversity #####
m.pollinator_0 <- glmmTMB(pollinator_0 ~ patch + (1|block),
                          data = diversity_metrics)
summary(m.pollinator_0)
plot(simulateResiduals(m.pollinator_0))
anova.pollinator_0 <- Anova(m.pollinator_0, type = "III")
confint(m.pollinator_0)



# shannon all species
m.pollinator_1 <- glmmTMB(pollinator_1 ~ patch + (1|block),
                      data = diversity_metrics)
summary(m.pollinator_1)
plot(simulateResiduals(m.pollinator_1))
anova.pollinator_1 <- Anova(m.pollinator_1, type = "III")
confint(m.pollinator_1)



#### tables ####
anova.abundance.df <- as.data.frame(anova.abundance)
anova.abundance.df <- anova.abundance.df %>%
  rownames_to_column("model_term") %>%
  mutate(`Response variable` = "Total abundance of interactions")

anova.floral_0.df <- as.data.frame(anova.floral_0)
anova.floral_0.df <- anova.floral_0.df %>%
  rownames_to_column("model_term") %>%
  mutate(`Response variable` = "Floral richness (q=0)")

anova.floral_1.df <- as.data.frame(anova.floral_1)
anova.floral_1.df <- anova.floral_1.df %>%
  rownames_to_column("model_term") %>%
  mutate(`Response variable` = "Floral diversity (q=1)")

anova.pollinator_0.df <- as.data.frame(anova.pollinator_0)
anova.pollinator_0.df <- anova.pollinator_0.df %>%
  rownames_to_column("model_term") %>%
  mutate(`Response variable` = "Pollinator richness (q=0)")

anova.pollinator_1.df <- as.data.frame(anova.pollinator_1)
anova.pollinator_1.df <- anova.pollinator_1.df %>%
  rownames_to_column("model_term") %>%
  mutate(`Response variable` = "Pollinator diversity (q=1)")


m.diversity_anova_all <- rbind(
  anova.abundance.df, anova.floral_0.df, anova.floral_1.df, anova.pollinator_0.df, anova.pollinator_1.df
)


rename_variable_anova <- tibble(model_term = c("(Intercept)", "patch"),
                                Predictor = c("Intercept", "Patch Type"))


m.diversity_anova_all <- m.diversity_anova_all %>%
  left_join(rename_variable_anova, by = "model_term") %>%
  dplyr::select(`Response variable`, Predictor, Chisq, Df, `Pr(>Chisq)`) %>%
  dplyr::rename(p.value = `Pr(>Chisq)`, df = Df)

table_diversity_anova <- m.diversity_anova_all %>%
  kbl(digits = 3) %>%
  kable_classic(full_width = T) %>%
  kable_styling(html_font = "Times New Roman",
                font_size = 16) %>%
  collapse_rows(columns = c(1), target = 1) %>%
  row_spec(0, extra_css = "border-bottom: 5px double;") %>%
  row_spec(1:nrow(m.diversity_anova_all), extra_css = "border-bottom: 1px solid;") %>%
  row_spec(0:nrow(m.diversity_anova_all), extra_css = "padding-bottom: 5px;")
table_diversity_anova

# exporting
# save_kable(table_diversity_anova, file = file.path("tables", "table_diversity_anova.html"))




#### plotting ####
# abundance plot
# model predictions for plotting
m.abundance.df <- ggpredict(m.abundance, terms = c("patch"), back_transform = TRUE)
# plotting
abundance.pred <- m.abundance.df %>%
  ggplot() +
  geom_jitter(aes(x = patch, y = abundance, color = patch), data = diversity_metrics, size = 6, alpha = 0.55,
              width = 0.08, height = 0) +
  geom_errorbar(aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high, fill = x), color = "black",
                data = m.abundance.df, width = 0, linewidth = 2.5) +
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
  ylab(expression("Interaction abundance")) +
  theme(legend.position = "none") 
abundance.pred

# pdf(file = file.path("plots", "diversity-plots", "abundance.pred.pdf"), width = 7, height = 6)
# abundance.pred
# dev.off()


# floral 0
# model predictions for plotting
m.floral_0.df <- ggpredict(m.floral_0, terms = c("patch"), back_transform = TRUE)
# plotting
floral_0.pred <- m.floral_0.df %>%
  ggplot() +
  geom_jitter(aes(x = patch, y = floral_0, color = patch), data = diversity_metrics, size = 6, alpha = 0.55,
              width = 0.08, height = 0) +
  geom_errorbar(aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high, fill = x), color = "black",
                data = m.floral_0.df, width = 0, linewidth = 2.5) +
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
  ylab(expression("Plant richness (q = 0)")) +
  theme(legend.position = "none") 
floral_0.pred

# pdf(file = file.path("plots", "diversity-plots", "floral_0.pred.pdf"), width = 7, height = 6)
# floral_0.pred
# dev.off()



# Shannon all species
# model predictions for plotting
m.floral_1.df <- ggpredict(m.floral_1, terms = c("patch"), back_transform = TRUE)
# plotting
floral_1.pred <- m.floral_1.df %>%
  ggplot() +
  geom_jitter(aes(x = patch, y = floral_1, color = patch), data = diversity_metrics, size = 6, alpha = 0.55,
              width = 0.08, height = 0) +
  geom_errorbar(aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high, fill = x), color = "black",
                data = m.floral_1.df, width = 0, linewidth = 2.5) +
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
  ylab(expression("Plant diversity (q = 1)")) +
  theme(legend.position = "none") 
floral_1.pred
# 
# pdf(file = file.path("plots", "diversity-plots", "floral_1.pred.pdf"), width = 7, height = 6)
# floral_1.pred
# dev.off()



# pollinator 
# model predictions for plotting
m.pollinator_0.df <- ggpredict(m.pollinator_0, terms = c("patch"), back_transform = TRUE)
# plotting
pollinator_0.pred <- m.pollinator_0.df %>%
  ggplot() +
  geom_jitter(aes(x = patch, y = pollinator_0, color = patch), data = diversity_metrics, size = 6, alpha = 0.55,
              width = 0.08, height = 0) +
  geom_errorbar(aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high, fill = x), color = "black",
                data = m.pollinator_0.df, width = 0, linewidth = 2.5) +
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
  ylab(expression("Pollinator richness (q = 0)")) +
  theme(legend.position = "none") 
pollinator_0.pred

# pdf(file = file.path("plots", "diversity-plots", "pollinator_0.pred.pdf"), width = 7, height = 6)
# pollinator_0.pred
# dev.off()



# Shannon all species
# model predictions for plotting
m.pollinator_1.df <- ggpredict(m.pollinator_1, terms = c("patch"), back_transform = TRUE)
# plotting
pollinator_1.pred <- m.pollinator_1.df %>%
  ggplot() +
  geom_jitter(aes(x = patch, y = pollinator_1, color = patch), data = diversity_metrics, size = 6, alpha = 0.55,
              width = 0.08, height = 0) +
  geom_errorbar(aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high, fill = x), color = "black",
                data = m.pollinator_1.df, width = 0, linewidth = 2.5) +
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
  ylab(expression("Pollinator diversity (q = 1)")) +
  theme(legend.position = "none") 
pollinator_1.pred

# pdf(file = file.path("plots", "diversity-plots", "pollinator_1.pred.pdf"), width = 7, height = 6)
# pollinator_1.pred
# dev.off()


# all diversity plots together
diversity_plot <- cowplot::plot_grid(floral_0.pred, pollinator_0.pred,floral_1.pred, pollinator_1.pred, 
                                     labels = c("A)", "B)", "C)", "D)"), label_size = 26, label_x = 0.19, label_y = 0.95)
diversity_plot

# pdf(file = file.path("plots", "diversity-plots", "diversity_plot.pdf"), width = 13, height = 11)
# diversity_plot
# dev.off()
# 
# 
# # supplemental plots
# pdf(file = file.path("plots", "diversity-plots", "si_abundance_plot.pdf"), width = 7, height = 6)
# abundance.pred
# dev.off()



