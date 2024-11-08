# Enns et al. 3
# library and working directory ####

library(readxl)
library(tidyverse)
library(stars)
library(sf)
library(sfnetworks)
library(SSN2)
library(SSNbler)
library(riverdist)
library(rivernet)
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
stream_net <- st_read("./stream_net/Gewässerstrukturgüte_Hessen_rev.shp") %>% 
  select(geometry, GESAMT, HP3, GEWKZ, ABS, NAME) %>% st_cast("LINESTRING") %>% 
  mutate(ABS = as.numeric(ABS)) 

# build initial LSN
lines_to_lsn(
  stream_net_dis,
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

# count number of WWTPs
mzb <- mzb %>% rowwise() %>%  mutate(NUM_WWTP =
        sum( ifelse( str_starts(wwtp$GEWKZ, GEWKZ), 1, 0))-
          sum(ifelse((wwtp$GEWKZ %in% GEWKZ) & (ABS >= wwtp$ABS), 1, 0))
)

## preparation of stream network

# group by stream ID (gwz) and unify, summing up shape_LEN
## Z. watershed delineation ####

# reduce stream network
stream_net_red7 <- stream_net %>% filter(str_length(GEWKZ) <= 7) %>% group_by(GEWKZ) %>% summarize(geometry = st_union(geometry)) %>% st_cast("LINESTRING")

# load in SRTM 30m raster
srtm_30 <- read_stars("./DEM/output_SRTMGL1.tif")

# union stream network
stream_net_dis <- stream_net %>% mutate(geometry = )

# rasterize stream network
test <- stream_net %>% st_transform(st_crs(srtm_30)) %>% st_rasterize(srtm_30, options = "ALL_TOUCHED=TRUE", align = T)