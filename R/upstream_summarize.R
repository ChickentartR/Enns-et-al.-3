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

# col - column from which the data should be summarized

# stat - method of summarization. When stat = 'sum' then col must be numeric 

# Code ####
upstream_summarize <- function(net, start, col, stat, count) {
  col <- ensym(col)  # Ensure the column is properly quoted for non-standard evaluation

  if (class(net)[1] != "sfnetwork") {
    stop("net must be an sfnetwork object!", call. = F)  # check for valid network input
  }
  
  if (!(stat %in% c("count", "sum"))) {
    stop("stat must either be 'count' or 'sum'!", call. = F)  # check for valid stat input
  }
  
  # Perform the analysis and filtering steps
  tbl <- suppressWarnings(igraph::shortest_paths(net, from = start, to = igraph::V(net), mode = "in")) %>%
    unlist(.$vpath) %>%
    unique() %>%
    st_as_sf(net, "nodes")[.,] %>%
    filter(!is.na(!!col))  # Filter out NA values
  
  # Handle the case when 'stat' is either 'count' or 'sum'
  if (stat == "count") {
    result <- nrow(tbl)  # Count the rows
  } 
  else {
    result <- sum(tbl[[as.character(col)]], na.rm = TRUE)  # Sum the values in the specified column
  }
  
  return(result)
}