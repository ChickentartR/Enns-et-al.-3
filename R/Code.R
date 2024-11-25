# Enns et al. 3
# library and working directory ####

library(readxl)
library(tidyverse)
library(stars)
library(sf)
library(sfnetworks)
library(tidygraph)
library(SSN2)
library(SSNbler)
library(geodata)
library(whitebox)
library(osmextract)
library(httr2)
library(tmap)

#install_whitebox()
wbt_init()
setwd("C:/Users/Daniel Enns/Documents/Promotion/MZB/Enns-et-al.-3/Data")


## 1 Stream network from HLNUG ####
## SSNbler Workflow

# load data
stream_net <- st_read("./stream_net/stream_net_rev.gpkg") %>% 
  select(geom, GESAMT, HP3, GEWKZ, ABS, NAME) %>% st_cast("LINESTRING") %>% 
  mutate(ABS = as.numeric(ABS)) 

# build initial LSN
lines_to_lsn(
  stream_net,
  "./SSNbler",
  check_topology = T,
  snap_tolerance = 3, 
  topo_tolerance = 30,
  overwrite = T
)

# indentify topological errors & restrictions

node_errors <- st_read("./node_errors.gpkg") %>% modify_if(is.character, as.factor)

summary(node_errors)

## Conclusions:
# It would take a lot of time to fix all network issues. In theory, the network issues would be fixed by 
# adjusting the parameters in the LSN building process and postprocessing using GIS. To allocate network-
# and simple features information to the biological sampling sites I use the stream- and segment number 
# system (see section X) instead of network routing algorithms. 

## 2. Bio sampling data ####

## select most recent sampling for macrobenthos ...
mzb <- read_delim("./EQC_MZB.csv", delim = ";") %>% 
  mutate(EQC = factor(EQC)) %>% group_by(ID_SITE) %>% arrange(desc(YEAR)) %>% slice_head() %>% 
  st_as_sf(coords = c("UTM_EAST", "UTM_NORTH"), crs = st_crs(25832)) 

# save as shapefile for watershed delineation
st_write(mzb, "./Bio_data/MZB.shp")

## ... for fish

## ... for diatoms

## ... for macrophytes

## check overlap with stream network

mzb <- cbind(mzb, stream_net[st_nearest_feature(mzb, stream_net),]) %>% select(-geometry.1)

## 3. WWTP & Storm overflow ####

# load data
wwtp <- read_delim("./WWTP.csv", delim = ";") %>% 
  st_as_sf(coords = c("UTM_EAST", "UTM_NORTH"), crs = st_crs(25832))

storm <- st_read()

# check overlap with stream network

wwtp <- cbind(wwtp, stream_net[st_nearest_feature(wwtp, stream_net),3:4]) %>% select(-geometry.1)

## 4. OSM Highways ####
### 1.5.1 query Highways from OSM ####

# buold SQL like query
hw_vectrotrans <- c(
  "-select", "osm_id,highway",
  "-where", "highway IN ('motorway')"
)

# read in data
highways <- oe_get("Hessen",
                   layer = "lines",
                   vectortranslate_options = mjr_vectrotrans,
                   boundary = st_as_sf(hess)
) %>% as_spatvector() 

# plot
plot(highways)

# save & load file
writeVector(stream_vec_osm,"highways.shp")
highways <- vect("highways.shp")

### 1.5.2 extract nodes from higway-stream net. crossings ####

## Y. Delineate catchments ####

## Burn DEM

## D8 flow accumulation

## D8 pointer grid

## snap pour points

## delineate watersheds

## X. Data allocation ####

end_points <- stream_net %>% filter(ABS == 1) %>% st_boundary() %>% st_cast("POINT") %>% group_by(GEWKZ) %>% slice_tail()

test <- st_nearest_feature(end_points, stream_net)

st_write(test,"./stream_net/test.shp")

# count number of WWTPs
mzb <- mzb %>% rowwise() %>%  mutate(NUM_WWTP =
        sum( ifelse( str_starts(wwtp$GEWKZ, GEWKZ), 1, 0))-
          sum(ifelse((wwtp$GEWKZ %in% GEWKZ) & (ABS >= wwtp$ABS), 1, 0))
)

## Allocate Data by routing ####

# round coordinates of stream network
str_net_rout <- st_set_precision(stream_net, 0.01)

# create sfnetwork object
net_rout <- as_sfnetwork(str_net_rout) %>% convert(to_spatial_simple) %>% convert(to_spatial_smooth)

# blend in Bio and stressor data as nodes
net_rout_blend <- st_network_blend(net_rout, mzb) %>% st_network_blend(.,wwtp)

#extract bio nodes for querying
nodes <- st_as_sf(net_rout_blend, "nodes")
data_allocated <- nodes %>% filter(!is.na(ID_SITE)) %>% 
  mutate(ID_NODE = row.names(nodes)[with(nodes, !is.na(ID_SITE))])

# write function to filter for upstream reaches and count/sum stressor variables
upstream_count <- function(net, start, col, stat) {
  col <- ensym(col)  # Ensure the column is properly quoted for non-standard evaluation
  
  if (class(net)[1] != "sfnetwork") {
    stop("net must be an sfnetwork object!") #check for valid network input
  }
  
  if (!(stat %in% c("count", "sum"))) {
    stop("stat must either be 'count' or 'sum'!") # check for valid stat input
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

upstream_count.2(net_rout_blend ,20348, CON_HOUSE, "count")

# summarize colmuns values
test <- st_as_sf(net_rout_blend, "nodes")[x,] %>% st_drop_geometry() %>% filter(!is.na(CON_HOUSE)) %>% summarize(sum(CON_HOUSE))

# count instances
test <- st_as_sf(net_rout_blend, "nodes")[x,] %>% st_drop_geometry() %>% filter(!is.na(CON_HOUSE)) %>% count()

## Z. watershed delineation ####

# load burned DEM
DEM <- read_stars("./DEM/SRTMGL1_30m_2px_burned.tif")

# load rasterized stream network
stream_net_dis <- st_read("./stream_net/stream_net_dis.gpkg")

# breach and fill pits in raster
wbt_breach_depressions_least_cost(
  dem = "./DEM/SRTMGL1_30m_2px_burned.tif", 
  output = "./DEM/srtm_breach.tif", 
  dist = 7,
  fill = T
)

# Create flow accumulation & pointer grids
wbt_d8_flow_accumulation(
  input = "./DEM/srtm_breach.tif",
  output = "./DEM/D8FA.tif"
)

wbt_d8_pointer(
  dem = "./DEM/srtm_breach.tif",
  output = "./DEM/D8pointer.tif"
)

# snap points to stream raster
wbt_jenson_snap_pour_points(
  pour_pts = "./Bio_data/MZB.shp",
  streams = "./DEM/HLNUG_raster.tif",
  output = "./Bio_data/MZB_snap.shp",
  snap_dist = 300
)

# !!! MASK D8POINTER BEFORE DELINEATION !!!

# delineate watersheds
wbt_watershed(
  d8_pntr = "./DEM/D8pointer.tif",
  pour_pts = "./Bio_data/MZB_snap.shp",
  output = "./DEM/watersheds.tif"
)

# load and plot watersheds
watershed <- read_stars("./DEM/watersheds.tif") %>% st_as_sf(merge = T)


# extract streams only for/ move to README
wbt_extract_streams(
  flow_accum = "./DEM/D8FA.tif",
  output = "./DEM/stream_raster.tif",
  threshold = 700
)
stream_raster <- read_stars("./DEM/stream_raster.tif")
