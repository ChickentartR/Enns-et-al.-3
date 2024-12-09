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

# node_cols - columns from which the node data should be summarized

# IDs - IDs of nodes, which should be counted. Not the row number of nodes!

# area - sf object with polygon geometry for which upstream reach coverage should be calculated 

# area_cols - columns from which the area data should be summarized

# threshold - shrinkage value for the area polygon to not intersect with adjacent polygons
#             not covered by upstream reach. Pay attention to CRS units!

# Code ####
upstream_summarize <- function(net, start, node_cols, IDs, area = NULL, area_cols = NULL, threshold) {
  
  # disable traceback
  options(error = NULL)
  
  # Check for valid class of net input
  if (class(net)[1] != "sfnetwork") {
    stop("net must be an sfnetwork object!", call. = F)  # check for valid network input
  }
  
  # Check for valid area input
  if (!is.null(area)) {
    if (st_crs(net) != st_crs(area)) {
      stop("net and area are not in the same CRS!", call. = F)}# Check for same net and area CRS
    if (class(area)[1] != "sf" | any(st_geometry_type(area) != "POLYGON")) {
      stop("area must be an sf object which contains only polygon geometries!", call. = F)}  # check for valid network input
  }
  
  # Check presence of columns in data
  if (any(!(node_cols %in% colnames(st_as_sf(net,"nodes"))))) {
    stop("Column name not found in net!")
  }
  
  if (!is.null(area) & is.null(area_cols)) {
    stop("must provide area_cols if area is given!")
  }
  
  if (!is.null(area_cols)) {
  if (any(!(area_cols %in% colnames(area)))) {
    stop("Column name not found in area!")}
  }
  
  # Perform the analysis and filtering steps
  nodes_us <- suppressWarnings(igraph::shortest_paths(net, from = start, to = igraph::V(net), mode = "in")) %>%
    unlist(.$vpath) %>%
    unique()
    
  tab_nodes <-  st_as_sf(net, "nodes")[nodes_us,]

  # Calculate the sum of each column in 'cols'
  tab_sum <- tab_nodes %>% as.data.frame() %>% 
    summarise(across(.cols = all_of(node_cols), ~sum(.x, na.rm = T)),
              across(.cols = all_of(IDs), ~sum(!is.na(.x)), .names = "num_{.col}"))
  
  # select watersheds
  if (!is.null(area)) {
    area_poly <- st_filter(area, st_shift(tab_nodes, st_centroid(st_union(tab_nodes)), 0.001), .predicate = st_intersects) 
    
    if(any(!st_is(area_poly, "POLYGON") | length(st_is(area_poly, "POLYGON") == 0))) {stop("there were non-polygon geometries created!")}
    else {area_poly <- area_poly %>% st_union() %>% st_remove_holes() %>% st_buffer(dist = -threshold)}
    
    # sample area by filled area polygon and summarize
    area_sum <- st_filter(area, area_poly, .predicate = st_intersects) %>% 
      as.data.frame() %>% summarise(across(.cols = all_of(area_cols), ~sum(.x, na.rm = T)))
    
    tab_sum <- bind_cols(tab_sum, area_sum)
  }
  
  # Return the summarized table with sums
  return(tab_sum)
}