# Enns et al. 3
# library and working directory ####

# general
library(readxl)
library(tidyverse)

# Saptial analysis
library(geodata)
library(stars)
library(sf)
library(nngeo)
library(sfnetworks)
library(tidygraph)
library(SSN2)
library(SSNbler)
library(whitebox)
library(osmextract)
library(httr2)
library(tmap)

# modeling
library(caret)
library(mgcv)
library(xgboost)
library(pdp)
library(doParallel)

# install_whitebox()
wbt_init()
setwd("C:/Users/Daniel Enns/Documents/Promotion/MZB/Enns-et-al.-3/Data")

# source function for data allocation
source("C:/Users/Daniel Enns/Documents/Promotion/MZB/Enns-et-al.-3/R/upstream_summarize.R")
source("C:/Users/Daniel Enns/Documents/Promotion/MZB/Enns-et-al.-3/R/st_shift.R")

# I. Data generation ####
## 1 Stream network from HLNUG ####

# load stream network data
stream_net <- st_read("./stream_net/stream_net_rev.gpkg") %>% 
  select(geom, GESAMT, HP3, GEWKZ, ABS, NAME) %>% st_cast("LINESTRING") %>% 
  mutate(ABS = as.numeric(ABS)) 

# load outline of Hesse
hess <- gadm("Deu", level = 1, tempfile()) %>% subset(.$NAME_1 == "Hessen") %>% st_as_sf()

# SSNbler Workflow
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
# and simple features information to the biological sampling sites the precision of the stream network
# is lowered and a network is build for routing (section 7. Data allocation). 

## 2. Bio sampling data ####

# select most recent sampling for macrobenthos ...
mzb <- read_delim("./Bio_data/EQC_MZB.csv", delim = ";") 

# get average EQC & Taxa richness
mzb_avg <- mzb %>% group_by(ID_SITE) %>% summarize(EQC_mean = round(mean(EQC)), Taxa_mean = mean(NUM_TAXA))

# join with averages
mzb <- left_join(mzb, mzb_avg, "ID_SITE")

mzb <- mzb %>% mutate(EQC = factor(EQC)) %>% group_by(ID_SITE) %>% arrange(desc(YEAR)) %>% 
  slice_head() %>% st_as_sf(coords = c("UTM_EAST", "UTM_NORTH"), crs = st_crs(25832)) 

# get infos from stream network and remove sites from streams/rivers with uncomplete catchments
mzb <- cbind(mzb, stream_net[st_nearest_feature(mzb, stream_net),]) %>% select(-geom) %>% 
  filter(!GEWKZ  %in% c("2", "238", "24", "41", "4", "44", "239152", "23932", "2396", "23988", "24779892"), !ID_SAMPLE %in% c("1210718","1273528"))

# save as shapefile for watershed delineation
st_write(mzb, "./Bio_data/MZB.shp", append = F)

## ... for fish

## ... for diatoms

## ... for macrophytes

## 3. WWTP, Storm overflow & Dams ####

# load data
wwtp <- read_delim("./Antropo_data/WWTP.csv", delim = ";") %>% 
  st_as_sf(coords = c("UTM_EAST", "UTM_NORTH"), crs = st_crs(st_crs(stream_net)))

storm <- read_delim("./Antropo_data/SWOF.csv", delim = ";") %>% 
  st_as_sf(coords = c("Ostwert", "Nordwert"), crs = st_crs(st_crs(stream_net))) %>% 
  select(OBJECTID) %>% rename(SWOF_ID = OBJECTID)

barriers <- read_delim("./Antropo_data/AMBER/AMBER_all.csv", delim = ";") %>% 
  st_as_sf(coords = c("Longitude_WGS84", "Latitude_WGS84"), crs = st_crs(4326)) %>% 
  st_filter(hess) %>% select(GUID, type) %>% mutate(count = 1) %>% 
  pivot_wider(names_from = "type", values_from = count, values_fill = NA) %>% st_transform(st_crs(stream_net))

## 4. OSM Highways ####

# read osm data and filter for highways and railways
transport_net <- oe_read("./OSM/hessen-latest.osm.pbf", extra_tags = "railway") %>% 
select(osm_id, highway, railway, geometry) %>% filter(highway == "motorway" | railway == "rail") %>% 
  unite("type", highway, railway, na.rm = T) %>% st_transform(st_crs(stream_net))

# extract highway and railway crossings with stream net
crossings <- st_intersection(stream_net, transport_net) %>% select(osm_id, type, geom) %>% st_cast("POINT")
mtw <- crossings %>% filter(type == "motorway")
rail <- crossings %>% filter(type == "rail")

## 5. watershed delineation ####

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

# extract streams
wbt_extract_streams(
  flow_accum = "./DEM/D8FA.tif",
  output = "./DEM/stream_raster.tif",
  threshold = 700
)

# snap points to stream raster
wbt_jenson_snap_pour_points(
  pour_pts = "./Bio_data/MZB.shp",
  streams = "./DEM/stream_raster.tif",
  output = "./Bio_data/MZB_snap.shp",
  snap_dist = 300
)

# delineate watersheds
wbt_watershed(
  d8_pntr = "./DEM/D8pointer.tif",
  pour_pts = "./Bio_data/MZB_snap.shp",
  output = "./DEM/watersheds.tif"
)

# Convert watershed raster to polygon
ws <- read_stars("./DEM/watersheds.tif") %>% st_as_sf(as_points = F, merge = T)

## 6. sample CLC by watersheds ####

# load CLC
clc <- st_read("./CLC/CLC_2018.gpkg")

# intersect with watersheds and create dataframe
clc_ws <- st_intersection(clc, ws) %>% mutate(area = st_area(.)) %>% st_drop_geometry() %>% 
  as.data.frame() %>% select(-CLC_201) %>% filter(!is.na(type))

clc_tab <- clc_ws %>% group_by(watersheds.tif, type) %>% summarize(area = sum(area)) %>% 
  pivot_wider(names_from = type, values_from = area)

# combine watersheds with CLC area data
ws_clc <- left_join(ws, clc_tab, by = "watersheds.tif")

## 7. Data allocation ####

# round coordinates of stream network
str_net_rout <- st_set_precision(stream_net, 0.01)

# create sfnetwork object
net_rout <- as_sfnetwork(str_net_rout) %>% convert(to_spatial_simple) %>% convert(to_spatial_smooth)

# blend in Bio and stressor data as nodes
net_rout_blend <- st_network_blend(net_rout, mzb) %>% st_network_blend(.,wwtp) %>% 
  st_network_blend(.,storm) %>% st_network_blend(.,barriers) %>% st_network_blend(.,mtw) %>% 
  st_network_blend(.,rail)

#extract bio nodes for querying
nodes <- st_as_sf(net_rout_blend, "nodes")
bio_data <- nodes %>% filter(!is.na(ID_SITE)) %>% 
  mutate(ID_NODE = row.names(nodes)[with(nodes, !is.na(ID_SITE))])

data_all <- bio_data %>% rowwise() %>% 
  mutate(upstream_summarize(
    net = net_rout_blend,
    start = ID_NODE,
    node_cols = c("CON_HOUSE", "TOT_WASTE", "other", "weir", "culvert", "ford", "dam", "ramp"),
    IDs = c("ID_WWTP", "SWOF_ID", "GUID", "osm_id.x", "osm_id.y"),
    area = ws_clc,
    area_cols = c("Agriculture", "Urban", "semi-Natural", "Wetland"),
    threshold = 30)
    )

# II. Modeling ####

## 1. Data preparation ####

# data cleaning
data4model <- data_all %>% cbind(st_coordinates(.$geom)) %>% as.data.frame() %>%
  rename(cros_hw = num_osm_id.x, cros_rail = num_osm_id.y, semi_Natural = semi.Natural) %>% 
  select(EQC_mean, 
         GESAMT, 
         TOT_WASTE,
         other:ramp,
         num_ID_WWTP:Y,
         RIV_TYPE
  ) %>% 
  mutate(EQC_mean = replace(EQC_mean, EQC_mean == 1, 2),
         across(c(EQC_mean,GESAMT,RIV_TYPE), as.factor), 
         across(TOT_WASTE:Y, as.numeric),
         Urban_prop = Urban / rowSums(.[15:18]),
         Agriculture_prop = Agriculture / rowSums(.[15:18]),
         semi_Natural_prop = semi_Natural / rowSums(.[15:18]))

# create training and test datasets
set.seed(1234)
inTrain <- createDataPartition(
  y = data4model$EQC_mean,
  p = 0.8,
  list = F
)

train <- data4model[inTrain,] 
test <- data4model[-inTrain,] 

# plot Class counts
ggplot(train)+
  geom_histogram(aes(x = EQC_mean), stat = "count")+
  theme_bw()

## 2. GAM ####

# build initial model
gam_model <- gam(as.numeric(EQC_mean) ~ s(TOT_WASTE, by = interaction(RIV_TYPE,GESAMT)) + s(other, by = interaction(RIV_TYPE,GESAMT)) + s(weir, by = interaction(RIV_TYPE,GESAMT)) +
                   s(culvert, by = interaction(RIV_TYPE,GESAMT)) + s(ford, by = interaction(RIV_TYPE,GESAMT)) + s(dam, by = interaction(RIV_TYPE,GESAMT)) + s(ramp, by = interaction(RIV_TYPE,GESAMT)) + s(num_ID_WWTP, by = interaction(RIV_TYPE,GESAMT)) +
                   s(num_SWOF_ID, by = interaction(RIV_TYPE,GESAMT)) + s(Agriculture_prop, by = interaction(RIV_TYPE,GESAMT)) + s(Urban_prop, by = interaction(RIV_TYPE,GESAMT)) + s(semi_Natural_prop, by = interaction(RIV_TYPE,GESAMT)) +
                   s(cros_hw, by = interaction(RIV_TYPE,GESAMT)) + s(cros_rail, by = interaction(RIV_TYPE,GESAMT)) + s(X,Y, by = interaction(RIV_TYPE,GESAMT))+s(Agriculture_prop, Urban_prop, semi_Natural_prop, by = interaction(RIV_TYPE,GESAMT)) + RIV_TYPE * GESAMT,
                 select = T,
                 data = train,
                 family = ocat(R = 4)
                 )

# check model
gam.check(gam_model)

# summary and plot
summary(gam_model)
plot.gam(gam_model, trans = plogis, scheme = 2)

# check concurvity
concurvity(gam_model, full = T)



## 3. XGBoost ####

# extract target
target <- as.numeric(train$EQC_mean)-1

# one-hot encode rivtype
train <- train %>% mutate(rv_count = 1) %>% 
pivot_wider(names_from = RIV_TYPE, values_from = rv_count, values_fill = 0, names_prefix = "RT_")

# extract data & remove redundant info
train_ini <- train %>% 
  select(-EQC_mean,
         -Wetland) %>% 
  sapply(as.numeric)

# create xgb data matrix
dtrain <- xgb.DMatrix(data = train_ini, label = target)

# set model parameters
xgb_params <- list(
  "objective" = "multi:softprob",
  "eval_metric" = "mlogloss",
  "num_class" = 4
)

# build initial model
xgb_model <- xgb.cv(
  data = dtrain,
  params = xgb_params,
  nrounds = 250,
  nfold = 10,
  prediction = T
  )

# predictions & classification error
OOF_prediction <- data.frame(xgb_model$pred) %>%
  mutate(max_prob = max.col(., ties.method = "last"),
         label = target + 1)
head(OOF_prediction)

# confusion matrix
confusionMatrix(factor(OOF_prediction$max_prob),
                factor(OOF_prediction$label),
                mode = "everything")

# train model
train_model <- xgb.train(
  params = xgb_params,
  data = dtrain,
  nrounds = 250
)

# feature importance
xgb.importance(feature_names = colnames(train_ini), model = train_model) %>% 
xgb.ggplot.importance()

# Hyper-parameter tuning
xgbFit_4cl <- train %>% as.data.frame() %>% mutate(EQC_mean = factor(EQC_mean, levels = c("good", "moderate", "bad", "very_bad"))) %>% 
  caret::train(
    form = EQC_mean ~ .,
    data = .,
    method = "xgbTree",
    tuneLength = 10,
    trControl = trainControl(
      method = "cv", 
      number = 10,
      classProbs = T,
      summaryFunction = mnLogLoss
    ),
    objective = "multi:softprob",
    num_class = 4
  )

plot(varImp(xgbFit_4cl))
partial(xgbFit_4cl, pred.var = "X", train = train, prob = T, which.class = "good") %>% autoplot(smooth = T)

# feature importance
xgb.importance(feature_names = colnames(train_fin), model = final_model) %>% 
  xgb.ggplot.importance()

# Partial dependence plots
partial(final_model, pred.var = "semi_Natural_prop", train = train_fin, prob = T, which.class = 1) %>% autoplot(smooth = T)

# model performance with test data

## 4. Model comparisons ####

# wrangle test data
target_test <- as.numeric(test$EQC_mean)-1
dtest <- test %>% select(-EQC_mean) %>% as.numeric() %>% xgb.DMatrix(data = ., label = target_test)

# Confusion matrix
# GAM
gampred <- predict.gam(gam_model, test, type = "response") %>% as.data.frame() %>% 
  rowwise() %>% mutate(pred = colnames(.)[which.max(c_across(1:4))] %>% parse_number())
gampred$pred <- (gampred$pred + 1) %>% as_factor()
confusionMatrix(gampred$pred, test$EQC_mean)

# xgboost
confusionMatrix(
  predict(xgbFit_4cl, ),
)
