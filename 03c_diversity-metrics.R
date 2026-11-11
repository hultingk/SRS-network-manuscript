# -------------------------------------- #
#### Diversity metrics calculations script ####
# Author: XXXX
# -------------------------------------- #
# end result - summarized dataframe ready for analysis

# loading libraries
librarian::shelf(tidyverse, vegan, glmmTMB, iNEXT)

# loading data 
pollinator <- read.csv(file = file.path("data", "L1_data", "SRS-plant-pollinator.csv")) # plant-pollinator interaction data


# creating unique ID
pollinator <- pollinator %>%
  mutate(unique_ID = paste(block, patch, sep = ".")) %>%
  dplyr::select(c("unique_ID", "order", "family", "pollinator_species", "flower_species")) 

#### Pollinator abundance ####
# # total pollinator visitation
abundance <- pollinator %>%
  dplyr::count(unique_ID) %>%
  mutate(abundance = n) %>%
  separate(unique_ID, into = c("block", "patch"), remove = F) %>%
  dplyr::select(!c("n"))

#### Flower diversity ####
# Hill numbers #
floral_inext <- pollinator %>%
  dplyr::count(unique_ID, flower_species) %>% # counting # of each flower species per patch
  pivot_wider(names_from = unique_ID, values_from = n, values_fill = 0) %>% # pivoting into wider format for flower species by patch matrix
  column_to_rownames("flower_species")

floral_hill <- iNEXT(floral_inext, q = c(0,1,2), datatype = "abundance") # hill numbers
floral_est <- floral_hill$iNextEst[["coverage_based"]] # getting coveraged based estimates
floral_target_cov <- floral_est %>% # figuring out what is the lowest max sampling coverage - will use estimates from that coverage for all patches
  filter(Method == "Observed") %>%
  group_by(Assemblage) %>%
  summarise(max_cov = max(SC)) %>%
  summarise(min(max_cov)) %>%
  pull()

floral_est <- floral_est %>% # getting estimates as close to the target sampling coverage as possible
  group_by(Assemblage, Order.q) %>%
  slice_min(abs(SC - floral_target_cov), n = 1) %>%
  ungroup()

floral_est <- floral_est %>% # wrangling into format
  separate(Assemblage, into = c("block", "patch"), remove = FALSE) %>%
  mutate(Order.q = as.factor(Order.q)) %>%
  select(Assemblage, block, patch, Order.q, qD) %>%
  pivot_wider(names_from = Order.q, values_from = qD) %>%
  dplyr::rename(unique_ID = Assemblage, floral_0 = `0`, floral_1 = `1`, floral_2 = `2`)


# pollinator diversity #
pollinator_inext <- pollinator %>%
  dplyr::count(unique_ID, pollinator_species) %>% # counting # of each pollinator species per patch
  pivot_wider(names_from = unique_ID, values_from = n, values_fill = 0) %>% # pivoting into wider format for flower species by patch matrix
  column_to_rownames("pollinator_species")

pollinator_hill <- iNEXT(pollinator_inext, q = c(0,1,2), datatype = "abundance") # hill numbers
pollinator_est <- pollinator_hill$iNextEst[["coverage_based"]] # getting coveraged based estimates
pollinator_target_cov <- pollinator_est %>% # figuring out what is the lowest max sampling coverage - will use estimates from that coverage for all patches
  filter(Method == "Observed") %>%
  group_by(Assemblage) %>%
  summarise(max_cov = max(SC)) %>%
  summarise(min(max_cov)) %>%
  pull()

pollinator_est <- pollinator_est %>% # getting estimates as close to the target sampling coverage as possible
  group_by(Assemblage, Order.q) %>%
  slice_min(abs(SC - pollinator_target_cov), n = 1) %>%
  ungroup()

pollinator_est <- pollinator_est %>% # wrangling into format
  separate(Assemblage, into = c("block", "patch"), remove = FALSE) %>%
  mutate(Order.q = as.factor(Order.q)) %>%
  select(Assemblage, block, patch, Order.q, qD) %>%
  pivot_wider(names_from = Order.q, values_from = qD) %>%
  dplyr::rename(unique_ID = Assemblage, pollinator_0 = `0`, pollinator_1 = `1`, pollinator_2 = `2`)


# excluding two dominant pollinators: No Poecilognathus sulphureus or Apis mellifera - Hill numbers #
pollinator_noPoe_Apis_inext <- pollinator %>%
  filter(!pollinator_species %in% c("Poecilognathus sulphureus", "Apis mellifera")) %>%
  dplyr::count(unique_ID, pollinator_species) %>% # counting # of each pollinator species per patch
  pivot_wider(names_from = unique_ID, values_from = n, values_fill = 0) %>% # pivoting into wider format for flower species by patch matrix
  column_to_rownames("pollinator_species")

pollinator_noPoe_Apis_hill <- iNEXT(pollinator_noPoe_Apis_inext, q = c(0,1,2), datatype = "abundance") # hill numbers
pollinator_noPoe_Apis_est <- pollinator_noPoe_Apis_hill$iNextEst[["coverage_based"]] # getting coveraged based estimates
pollinator_noPoe_Apis_target_cov <- pollinator_noPoe_Apis_est %>% # figuring out what is the lowest max sampling coverage - will use estimates from that coverage for all patches
  filter(Method == "Observed") %>%
  group_by(Assemblage) %>%
  summarise(max_cov = max(SC)) %>%
  summarise(min(max_cov)) %>%
  pull()

pollinator_noPoe_Apis_est <- pollinator_noPoe_Apis_est %>% # getting estimates as close to the target sampling coverage as possible
  group_by(Assemblage, Order.q) %>%
  slice_min(abs(SC - pollinator_noPoe_Apis_target_cov), n = 1) %>%
  ungroup()

pollinator_noPoe_Apis_est <- pollinator_noPoe_Apis_est %>% # wrangling into format
  separate(Assemblage, into = c("block", "patch"), remove = FALSE) %>%
  mutate(Order.q = as.factor(Order.q)) %>%
  select(Assemblage, block, patch, Order.q, qD) %>%
  pivot_wider(names_from = Order.q, values_from = qD) %>%
  dplyr::rename(unique_ID = Assemblage, pollinator_noPoe_Apis_0 = `0`, pollinator_noPoe_Apis_1 = `1`, pollinator_noPoe_Apis_2 = `2`)



# putting all together
diversity_metrics <- abundance %>%
  left_join(floral_est, by = c("unique_ID", "block", "patch")) %>%
  left_join(pollinator_est, by = c("unique_ID", "block", "patch")) %>%
  left_join(pollinator_noPoe_Apis_est, by = c("unique_ID", "block", "patch"))
 
write.csv(diversity_metrics, file = file.path("data", "L3_metrics", "diversity_metrics.csv"))

