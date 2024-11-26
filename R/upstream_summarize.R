###################################
### upstream_summarize function ###
###################################
# Description #### 
# This function filters a sfnetwork from a given node to all of its 
# upstream reaches and summarizes information of other specified nodes present 
# in these sections.

# Dependent packages:
# The function requires the dplyr and sfnetworks packages

# Arguments:
# net - the network from which the function calculates upstream reaches. Must be sfnetwork object.

# start - row number of node, from which upstream reache should be calculated.

# cols - columns from which the data should be summarized

# IDs - IDs of nodes, which should be counted. Not the row number of nodes!

# overwrite - should output be written in input columns? 

# Code ####
upstream_summarize <- function(net, start, cols, IDs, overwrite) {
  
  if (class(net)[1] != "sfnetwork") {
    stop("net must be an sfnetwork object!", call. = F)  # check for valid network input
  }
  
  if(typeof(overwrite) != "logical") {
    stop("overwrite must be TRUE or FALSE!")
  }
  
  # Perform the analysis and filtering steps
  tbl <- suppressWarnings(igraph::shortest_paths(net, from = start, to = igraph::V(net), mode = "in")) %>%
    unlist(.$vpath) %>%
    unique() %>%
    st_as_sf(net, "nodes")[.,]

  # Calculate the sum of each column in 'cols'
  if(overwrite == F){
  tbl_sum <- tbl %>% as.data.frame() %>% 
    summarise(across(.cols = all_of(cols), ~sum(.x, na.rm = T), .names = "sum_{.col}"),
              across(.cols = all_of(IDs), ~sum(!is.na(.x)), .names = "num_{.col}"))
  } else if(overwrite == T) {
    tbl_sum <- tbl %>% as.data.frame() %>% 
      summarise(across(.cols = all_of(cols), ~sum(.x, na.rm = T)),
                across(.cols = all_of(IDs), ~sum(!is.na(.x)), .names = "num_{.col}"))
                }
  
  # Return the summarized table with sums
  return(as.data.frame(tbl_sum))
}