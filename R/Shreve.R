#######################
### Shreve function ###
#######################
# Description #### 
# This function calculates the SHreve stream magnitude, that is
# the number of sources upstream of a given point. The function
# is a simpler derivate of the upstream_summarize function.

# Dependent packages:
# The function requires the dplyr and sfnetworks packages

# Arguments:
# net - the network from which the function calculates upstream reaches. Must be sfnetwork object.

# start - row number of node, from which upstream reache should be calculated.

# Code ####

# disable traceback
options(error = NULL)

# Check for valid class of net input
if (!is.sfnetwork(net)) {
  stop("net must be an sfnetwork object!", call. = F)  # check for valid network input
}

# Check for presence of start node in net
if (start %in% rownames(st_as_sf(net, "nodes"))) {
  stop("start node not present in net!", call. = F)
}

Shreve <- function(net, start) {
  
  # Perform filtering steps
  nodes_us <- suppressWarnings(igraph::shortest_paths(net, from = start, to = igraph::V(net), mode = "in")) %>%
    unlist(.$vpath) %>%
    unique()
  
  sub_nodes <- st_as_sf(net, "nodes") %>% mutate(igraph_ID = seq.int(nrow(.))) %>% .[nodes_us,]
  
  # Calculate Shreve
  tab_shreve <- tibble(
    Shreve = st_as_sf(net, "edges") %>% filter(from %in% sub_nodes$igraph_ID & !(from %in% to)) %>% nrow()
  )
  
  # Return the summarized table with sums
  return(tab_shreve)
}