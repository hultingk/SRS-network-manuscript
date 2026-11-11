# -------------------------------------- #
#### Functions script ####
# Author: XXXXX
# -------------------------------------- #
library(parallel, tidyverse)

# function to get data into correct format for network analysis
prepare_matrix <- function(df) {
  df_wide <- df %>% 
    pivot_wider(names_from = pollinator_species, values_from = n, values_fill = 0) %>% # wide format
    dplyr::select(!c("unique_ID")) %>% # remove unique ID column
    column_to_rownames("flower_species") #convert years to rownames
}

####### NULL DISTRIBUTION FUNCTIONS #######
# null distributions
net.null.networklevel = function(nulls, fun, metric, level, web.names, cores = (detectCores() - 2)){
  net.null.metric <- mclapply(
    nulls,
    function(webs) {
      do.call(
        rbind,
        lapply(webs, fun, index = metric, level = level)
      )
    },
    mc.cores = cores)
  
  names(net.null.metric) <- webs.names
  return(net.null.metric)
}





# null distribution for extinction 
net.null.extinct <- function(nulls, participant, method){
  net.null.robust.metric  <- list()
  
  for (i in 1:length(nulls)) {
    net.null.extinct.metric <- lapply(nulls[[i]], second.extinct, participant = participant, method = method)
    net.null.robust.metric[[i]] = do.call('rbind', 
                                          lapply(net.null.extinct.metric, robustness))
    }
  names(net.null.robust.metric) <- webs.names
  return(net.null.robust.metric)
    
}




####### Z SCORE FUNCTIONS #########
# Z-score function for comparing different networks
net.zscore = function(obsval, nullval) {
  (obsval - mean(nullval))/sd(nullval)  
} 

zscore_metric = function(obsval, nullval, metric){
  zscore_values <- list() 
  for(i in 1:length(obsval)){
    zscore_values[[i]] = net.zscore(obsval[[i]][metric], 
                                         nullval[[i]][ ,metric])
  }
  names(zscore_values) <- webs.names
  return(zscore_values)
}
