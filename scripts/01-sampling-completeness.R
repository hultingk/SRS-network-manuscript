# -------------------------------------- #
#### Sampling completeness script ####
# Author: XXXXX
# -------------------------------------- #

# loading libraries
librarian::shelf(tidyverse, iNEXT, bipartite, kableExtra)

# loading data 
pollinator <- read.csv(file = file.path("data", "L1_data", "SRS-plant-pollinator.csv"))


# getting plant-pollinator data into correct format for sampling completeness estimation
sc_pollinator <- pollinator %>%
  mutate(interaction = paste(pollinator_species, flower_species, sep = "-")) %>% # combining pollinator and flower into one column - creating interaction
  mutate(unique_ID = paste(block, patch, sep = ".")) %>% # creating unique ID for each patch
  dplyr::count(unique_ID, interaction) %>% # counting # of each interaction per patch
  pivot_wider(names_from = unique_ID, values_from = n, values_fill = 0) # pivoting into wider format for interaction by patch matrix

# excluding patch from interaction matrix
sc_nosite_pollinator <- as.data.frame(sc_pollinator[, 2:15])

# estimating sampling completeness
out <- iNEXT(sc_nosite_pollinator, q=c(1), datatype="abundance", endpoint=NULL, se = TRUE) 

# mean sampling completeness across all patches
mean(out[["DataInfo"]][["SC"]])
sc.table <- out[["DataInfo"]]
sc.table <- sc.table %>%
  dplyr::select("Assemblage", "n", "SC") %>%
  separate(Assemblage, into = c("Block", "Patch")) %>%
  dplyr::rename(`Sampling Coverage` = SC, `Number of Observations` = n) %>%
  mutate(Patch = if_else(Patch == "B", "Connected", "Unconnected"))

# table of sampling completeness
sc.table <- sc.table %>%
  kbl(digits = 3) %>%
  kable_classic(full_width = T) %>%
  kable_styling(html_font = "Times New Roman",
                font_size = 16) %>%
  row_spec(0, extra_css = "border-bottom: 5px double;") %>%
  row_spec(1:nrow(sc.table), extra_css = "border-bottom: 1px solid;")
sc.table

# exporting
# save_kable(sc.table, file = file.path("tables", "sc.table.html"))



# visualizing sampling completeness
ggiNEXT(out, type=1, facet.var="Order.q") # Sample-size-based R/E curve
ggiNEXT(out, type=2, facet.var="Order.q") # Sample completeness curve
ggiNEXT(out, type=3, facet.var="Order.q") # Coverage-based R/E curve


sampling_plot <- ggiNEXT(out, type=2, facet.var="Order.q") # Sample completeness curve
sampling_plot <- sampling_plot + theme_minimal(base_size = 16)
sampling_plot

r_e_curve <- ggiNEXT(out, type=1, facet.var="Order.q") # Sample-size-based R/E curve
r_e_curve <- r_e_curve + theme_minimal(base_size = 16) + ylab("Interaction diversity") + theme(legend.position = "none")

total_coverage_plot <- cowplot::plot_grid(r_e_curve, sampling_plot, cols = 1)
total_coverage_plot

# 
# pdf(file = file.path("plots", "total_coverage_plot.pdf"), width = 7, height = 11)
# total_coverage_plot
# dev.off()

