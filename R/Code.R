# Enns et al. 3
# I. library and working directory ####

# general
library(readxl)
library(tidyverse)
library(dtplyr)
library(data.table)

# Saptial analysis
library(geodata)
library(stars)
library(sf)
library(spdep)
library(nngeo)
library(sfnetworks)
library(tidygraph)
library(whitebox)
library(osmextract)
library(httr2)
library(tmap)
library(units)

# modeling
library(caret)
library(mlr)
library(xgboost)
library(pdp)
library(doParallel)
library(parallelMap)
library(GWmodel)

# install_whitebox()
wbt_init()
setwd("C:/Users/Daniel Enns/Documents/Promotion/MZB/Enns-et-al.-3/Data")

# source function for data allocation
source("C:/Users/Daniel Enns/Documents/Promotion/MZB/Enns-et-al.-3/R/upstream_summarize.R")
source("C:/Users/Daniel Enns/Documents/Promotion/MZB/Enns-et-al.-3/R/st_shift.R")
source("C:/Users/Daniel Enns/Documents/Promotion/MZB/Enns-et-al.-3/R/Shreve.R")

# II. Data generation ####
## 1 Stream network from HLNUG ####

# load stream network data
stream_net <- st_read("./stream_net/stream_net_rev.gpkg") %>% 
  select(geom, GESAMT, HP3, GEWKZ, ABS, NAME) %>% st_cast("LINESTRING") %>% 
  mutate(ABS = as.numeric(ABS)) 

# load outline of Hesse
hess <- gadm("Deu", level = 1, tempfile()) %>% subset(.$NAME_1 == "Hessen") %>% st_as_sf()

## 2. Bio sampling data ####

# read in macrobenthos sampling
mzb <- read_delim("./Bio_data/EQC_MZB.csv", delim = ";") 

# get average EQC & Taxa richness
mzb_avg <- mzb %>% group_by(ID_SITE) %>% summarize(EQC_mean = round(mean(EQC)), Taxa_mean = mean(NUM_TAXA), MMI_mean = mean(MMI))

# join with averages
mzb <- left_join(mzb, mzb_avg, "ID_SITE")

mzb <- mzb %>% mutate(EQC = factor(EQC)) %>% group_by(ID_SITE) %>% arrange(desc(YEAR)) %>% 
  slice_head() %>% st_as_sf(coords = c("UTM_EAST", "UTM_NORTH"), crs = st_crs(25832)) 

# get infos from stream network and remove sites from streams/rivers with uncomplete catchments
mzb <- cbind(mzb, stream_net[st_nearest_feature(mzb, stream_net),]) %>% select(-geom) %>% 
  filter(!GEWKZ  %in% c("2", "238", "24", "41", "4", "44", "239152", "23932", "2396", "23988", "24779892"), !ID_SAMPLE %in% c("1210718","1273528"))

# save as shapefile for watershed delineation
st_write(mzb, "./Bio_data/MZB.shp", append = F)

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
str_net_rout <- st_set_precision(stream_net, set_units(70, "m"))

# create sfnetwork object
net_rout <- as_sfnetwork(str_net_rout) %>% convert(to_spatial_simple) %>% convert(to_spatial_smooth)

# blend in Bio and stressor data as nodes
net_rout_blend <- st_network_blend(net_rout, mzb) %>% st_network_blend(.,wwtp) %>% 
  st_network_blend(.,storm) %>% st_network_blend(.,barriers) %>% st_network_blend(.,mtw) %>% 
  st_network_blend(.,rail)

# extract bio nodes for querying
nodes <- st_as_sf(net_rout_blend, "nodes")
bio_data <- nodes %>% filter(!is.na(ID_SITE)) %>% 
  mutate(ID_NODE = row.names(nodes)[with(nodes, !is.na(ID_SITE))])

# set up cluster
num_cores <- (detectCores()-1)
num_cores %>% makeCluster() %>% registerDoParallel()

# split data into chunks
bio_data_chunks <- bio_data %>% split(., ceiling(seq_along(row_number(.)) / (length(row_number(.)) %/% num_cores)))

# run function parallel over data chunks
data_all <- foreach(chunk = bio_data_chunks, .combine = rbind, .packages = c("dplyr", "dtplyr","sf","sfnetworks","tidygraph","nngeo")) %dopar% {
  chunk %>% rowwise() %>% 
  mutate(upstream_summarize(
    net = net_rout_blend,
    start = ID_NODE,
    node_cols = c("CON_HOUSE", "TOT_WASTE", "other", "weir", "culvert", "ford", "dam", "ramp"),
    IDs = c("ID_WWTP", "SWOF_ID", "GUID", "osm_id.x", "osm_id.y"),
    area = ws_clc,
    area_cols = c("Agriculture", "Urban", "semi-Natural", "Wetland"),
    dist = "all",
    threshold = 30
    )
    )
}

# calculate Shreve stream magnitude

# set up cluster
num_cores <- (detectCores()-1)
num_cores %>% makeCluster() %>% registerDoParallel()

# split data into chunks
bio_data_chunks <- data_all %>% split(., ceiling(seq_along(row_number(.)) / (length(row_number(.)) %/% num_cores)))

# run function parallel over data chunks
data_all <- foreach(chunk = bio_data_chunks, .combine = rbind, .packages = c("dplyr", "dtplyr","sf","sfnetworks","tidygraph")) %dopar% {
  chunk %>% rowwise() %>% 
    mutate(
      Shreve(net_rout_blend, ID_NODE)
    )
}    

# if R crashes trying to work with the output dataframe, load in from csv
data_all <- read.csv("data_all.csv", sep = ";")

  
# III. Modeling ####
## 1. Data cleaning ####

# renaming
data4model <- data_all %>% cbind(st_coordinates(.$geom)) %>% as.data.frame() %>%
  rename(cros_hw = num_osm_id.x, 
         cros_rail = num_osm_id.y,
         dist_cros_hw_min = dist_osm_id.x_min,
         dist_cros_rail_min = dist_osm_id.y_min,
         semi_Natural = semi.Natural)

# renaming when data_all is loaded from csv
data4model <- data_all %>% as.data.frame() %>%
  rename(cros_hw = num_osm_id.x, 
         cros_rail = num_osm_id.y,
         dist_cros_hw_min = dist_osm_id.x_min,
         dist_cros_rail_min = dist_osm_id.y_min,
         semi_Natural = semi.Natural)

# remove units from values
data4model <- data4model %>% 
  mutate(across(.cols = c("dist_ID_WWTP_min", "dist_SWOF_ID_min", "dist_GUID_min", "dist_cros_hw_min", "dist_cros_rail_min", "Agriculture":"Wetland"),
                .fns = ~ as.vector(.x)))

# change values in distance
data4model <- data4model %>% 
  mutate(across(.cols = c("dist_ID_WWTP_min", "dist_SWOF_ID_min", "dist_GUID_min", "dist_cros_hw_min", "dist_cros_rail_min"),
                .fns = ~case_when(
                  .x == 0 ~ 5,
                  .x == Inf | .x == -Inf | is.na(.x) ~ 0,
                  TRUE ~ .x
                )))

# transform distances by 1/dist
data4model <- data4model %>% mutate(
  across(.cols= c("dist_ID_WWTP_min", "dist_SWOF_ID_min", "dist_GUID_min", "dist_cros_hw_min", "dist_cros_rail_min"),
         .fns = ~case_when(
           .x != 0 ~ .x / max(.x),
           .x == 0 ~ .x
         )
  )
)

# select variables and change classes
data4model <- data4model %>% 
  select(EQC_mean, 
         GESAMT, 
         TOT_WASTE,
         other:ramp,
         num_ID_WWTP:cros_rail,
         RIV_TYPE,
         MMI_mean:Y,
         Agriculture:Wetland,
         matches("_min$")
         )

## 2. EDA ####

# summary
summarizeColumns(data4model)

# correlation among variables
varicor <- data4model %>% select(where(is.numeric)) %>% cor(.) %>% round(.,2) %>% as.data.frame()


## 3. Spatial autocorrelation ####

# k-nearest neighbours
suppressPackageStartupMessages(require(deldir))
data_nb_knn <- knn2nb(knearneigh(data_all,k = 1))

# distance neighbours
dsts <- unlist(nbdists(data_nb_knn, data_all))
summary(dsts)
max_1nn <- max(dsts)

data_nb_dst <- dnearneigh(data_all, d1 = 0, d2 = 0.75*max_1nn)

data_lw <- nb2listw(data_nb_dst, style = "W", zero.policy = T)

# Moran statistics
moran.test(data_all$MMI_mean, data_lw, alternative = "greater")
moran.plot(data_all$MMI_mean, data_lw)

## 4. Feature engeneering ####

data4model <- data4model %>% 
  mutate(EQC_mean = replace(EQC_mean, EQC_mean == 1, 2),
         across(c(EQC_mean,GESAMT,RIV_TYPE), as.factor), 
         Urban_prop = Urban / rowSums(select(., Agriculture:Wetland)),
         Agriculture_prop = Agriculture / rowSums(select(., Agriculture:Wetland)),
         semi_Natural_prop = semi_Natural / rowSums(select(., Agriculture:Wetland))
  )

# standardize variables
data4model_std <- normalizeFeatures(data4model, target = c("MMI_mean", "EQC_mean"))

# PCA  
PC_pnt <- data4model_std %>% select(TOT_WASTE:cros_rail) %>% princomp(scores = T)
PC_dst <- data4model_std %>% select(matches("_min")) %>% princomp(scores = T)
PC_lcc <- data4model_std %>% select(Agriculture:semi_Natural) %>% princomp(scores = T)

summary(PC_pnt)
summary(PC_dst)
summary(PC_lcc)

PC_loadings <- bind_rows(
  PC_pnt$loadings[,1:2] %>% as.data.frame() %>% mutate(var = rownames(.), PCA = "points") %>% pivot_longer(cols = Comp.1:Comp.2, names_to = "component", values_to = "scores"),
  PC_dst$loadings[,1:4] %>% as.data.frame() %>% mutate(var = rownames(.), PCA = "distances") %>% pivot_longer(cols = Comp.1:Comp.4, names_to = "component", values_to = "scores"),
  PC_lcc$loadings[,1:2] %>% as.data.frame() %>% select(-Comp.2) %>% mutate(var = rownames(.), PCA = "land_use") %>% pivot_longer(cols = Comp.1, names_to = "component", values_to = "scores")
  )

# plot loadings
ggplot(
  data = PC_loadings,
  aes(x = component, y = scores, fill = var)
  )+
  geom_bar(
    position = "fill",
    stat = "identity",
    colour="#666666"
  )+
  facet_grid(factor(PCA, levels = c("distances","points","land_use"))~., scales ="free", space = "free")+
  coord_flip()+
  theme_bw()+
  theme(
    axis.title = element_blank()
  )

# PC scores data
data4model_sc <- data4model %>% select(MMI_mean, EQC_mean, GESAMT, RIV_TYPE, X, Y) %>% 
  mutate(
    pnt_comp1 = PC_pnt$scores[,1],
    pnt_comp2 = PC_pnt$scores[,2],
    dst_comp1 = PC_dst$scores[,1],
    dst_comp2 = PC_dst$scores[,2],
    dst_comp3 = PC_dst$scores[,3],
    dst_comp4 = PC_dst$scores[,4],
    lcc_comp1 = PC_lcc$scores[,1]
  )

### 4.1 feature engeneering for 2nd version ####

# cahnge meters to kilometers 
data4model_2nd <- data4model %>% select(-Shreve) %>% 
  mutate(EQC_mean = replace(EQC_mean, EQC_mean == 1, 2),
         across(c(EQC_mean,GESAMT,RIV_TYPE), as.factor), 
         across(X:Y, ~ .x / 1000),
         across(Agriculture:Wetland, ~ .x / 10^6)
  )   

# standardize variables
data4model_std_2nd <- normalizeFeatures(data4model_2nd, target = c("MMI_mean", "EQC_mean"))

## 5. XGBoost ####
### 5.1 Training and Test sets ####
set.seed(1234)
inTrain <- createDataPartition(
  y = data4model$EQC_mean,
  p = 0.8,
  list = F
)

train_coords <-data4model[inTrain,] %>% select(X,Y) 
test_coords <-data4model[-inTrain,] %>% select(X,Y) 

train_reg <- data4model[inTrain,] %>% select(-EQC_mean, -X,-Y) %>% 
  createDummyFeatures(target = "MMI_mean", cols = c("GESAMT", "RIV_TYPE")) %>% makeRegrTask(data = ., target = "MMI_mean", coordinates = train_coords)
test_reg <- data4model[-inTrain,] %>% select(-EQC_mean, -X,-Y) %>% 
  createDummyFeatures(target = "MMI_mean", cols = c("GESAMT", "RIV_TYPE")) %>% makeRegrTask(data = ., target = "MMI_mean", coordinates = test_coords)

train_clas <- data4model[inTrain,] %>% select(-MMI_mean, -X,-Y) %>% 
  createDummyFeatures(target = "EQC_mean", cols = c("GESAMT", "RIV_TYPE")) %>% makeClassifTask(data = ., target = "EQC_mean", coordinates = train_coords)
test_clas <- data4model[-inTrain,] %>% select(-MMI_mean, -X,-Y) %>% 
  createDummyFeatures(target = "EQC_mean", cols = c("GESAMT", "RIV_TYPE")) %>% makeClassifTask(data = ., target = "EQC_mean", coordinates = test_coords)


# standardized data
train_coords <-data4model_std[inTrain,] %>% select(X,Y) 
test_coords <-data4model_std[-inTrain,] %>% select(X,Y) 

train_reg_std <- data4model_std[inTrain,] %>% select(-EQC_mean, -X,-Y) %>% 
  createDummyFeatures(target = "MMI_mean", cols = c("GESAMT", "RIV_TYPE")) %>% makeRegrTask(data = ., target = "MMI_mean", coordinates = train_coords)
test_reg_std <- data4model_std[-inTrain,] %>% select(-EQC_mean, -X,-Y) %>%  
  createDummyFeatures(target = "MMI_mean", cols = c("GESAMT", "RIV_TYPE")) %>% makeRegrTask(data = ., target = "MMI_mean", coordinates = test_coords)


train_clas_std <- data4model_std[inTrain,] %>% select(-MMI_mean, -X,-Y) %>% 
  createDummyFeatures(target = "EQC_mean", cols = c("GESAMT", "RIV_TYPE")) %>% makeClassifTask(data = ., target = "EQC_mean", coordinates = train_coords)
test_clas_std <- data4model_std[-inTrain,] %>% select(-MMI_mean, -X,-Y) %>% 
  createDummyFeatures(target = "EQC_mean", cols = c("GESAMT", "RIV_TYPE")) %>% makeClassifTask(data = ., target = "EQC_mean", coordinates = test_coords)

# PCA scores
train_coords <-data4model_sc[inTrain,] %>% select(X,Y) 
test_coords <-data4model_sc[-inTrain,] %>% select(X,Y) 

train_reg_sc <- data4model_sc[inTrain,] %>% select(-EQC_mean, -X,-Y) %>% 
  createDummyFeatures(target = "MMI_mean", cols = c("GESAMT", "RIV_TYPE")) %>% makeRegrTask(data = ., target = "MMI_mean", coordinates = train_coords)
test_reg_sc <- data4model_sc[-inTrain,] %>% select(-EQC_mean, -X,-Y) %>% 
  createDummyFeatures(target = "MMI_mean", cols = c("GESAMT", "RIV_TYPE")) %>% makeRegrTask(data = ., target = "MMI_mean", coordinates = test_coords)

train_clas_sc <- data4model_sc[inTrain,] %>% select(-MMI_mean, -X,-Y) %>% 
  createDummyFeatures(target = "EQC_mean", cols = c("GESAMT", "RIV_TYPE")) %>% makeClassifTask(data = ., target = "EQC_mean", coordinates = train_coords)
test_clas_sc <- data4model_sc[-inTrain,] %>% select(-MMI_mean, -X,-Y) %>% 
  createDummyFeatures(target = "EQC_mean", cols = c("GESAMT", "RIV_TYPE")) %>% makeClassifTask(data = ., target = "EQC_mean", coordinates = test_coords)

# create learners
xgb_reg_learner <- makeLearner(
  "regr.xgboost",
  predict.type = "response"
)

xgb_clas_learner <- makeLearner(
  "classif.xgboost",
  predict.type = "response"
)

### 5.2 Hyper-parameter tuning ####
parallelStartSocket(detectCores()-1)
tuned_reg <- tuneParams(
  learner = xgb_reg_learner,
  task = train_reg,
  par.set = makeParamSet(makeIntegerParam("nrounds", lower = 100, upper = 500), 
                         makeIntegerParam("max_depth", lower = 1, upper = 10), 
                         makeNumericParam("eta", lower = 0.1, upper = 0.5), 
                         makeNumericParam("lambda",lower = -1, upper = 0, trafo = function(x) 10^x)),
  resampling = makeResampleDesc("SpCV", iters = 5),
  control = makeTuneControlRandom(maxit = 200)
)
parallelStop()

parallelStartSocket(detectCores()-1)
tuned_clas <- tuneParams(
  learner = xgb_clas_learner,
  task = train_clas,
  par.set = makeParamSet(makeIntegerParam("nrounds", lower = 100, upper = 500), 
                         makeIntegerParam("max_depth", lower = 1, upper = 10), 
                         makeNumericParam("eta", lower = 0.1, upper = 0.5), 
                         makeNumericParam("lambda",lower = -1, upper = 0, trafo = function(x) 10^x)),
  resampling = makeResampleDesc("SpCV", iters = 5),
  control = makeTuneControlRandom(maxit = 200)
)
parallelStop()

# standardized data
parallelStartSocket(detectCores()-1)
tuned_reg_std <- tuneParams(
  learner = xgb_reg_learner,
  task = train_reg_std,
  par.set = makeParamSet(makeIntegerParam("nrounds", lower = 100, upper = 500), 
                         makeIntegerParam("max_depth", lower = 1, upper = 10), 
                         makeNumericParam("eta", lower = 0.1, upper = 0.5), 
                         makeNumericParam("lambda",lower = -1, upper = 0, trafo = function(x) 10^x)),
  resampling = makeResampleDesc("SpCV", iters = 5),
  control = makeTuneControlRandom(maxit = 200)
)
parallelStop()

parallelStartSocket(detectCores()-1)
tuned_clas_std <- tuneParams(
  learner = xgb_clas_learner,
  task = train_clas_std,
  par.set = makeParamSet(makeIntegerParam("nrounds", lower = 100, upper = 500), 
                         makeIntegerParam("max_depth", lower = 1, upper = 10), 
                         makeNumericParam("eta", lower = 0.1, upper = 0.5), 
                         makeNumericParam("lambda",lower = -1, upper = 0, trafo = function(x) 10^x)),
  resampling = makeResampleDesc("SpCV", iters = 5),
  control = makeTuneControlRandom(maxit = 200)
)
parallelStop()

# PCA scores
parallelStartSocket(detectCores()-1)
tuned_reg_sc <- tuneParams(
  learner = xgb_reg_learner,
  task = train_reg_sc,
  par.set = makeParamSet(makeIntegerParam("nrounds", lower = 100, upper = 500), 
                         makeIntegerParam("max_depth", lower = 1, upper = 10), 
                         makeNumericParam("eta", lower = 0.1, upper = 0.5), 
                         makeNumericParam("lambda",lower = -1, upper = 0, trafo = function(x) 10^x)),
  resampling = makeResampleDesc("SpCV", iters = 5),
  control = makeTuneControlRandom(maxit = 200)
)
parallelStop()

parallelStartSocket(detectCores()-1)
tuned_clas_sc <- tuneParams(
  learner = xgb_clas_learner,
  task = train_clas_sc,
  par.set = makeParamSet(makeIntegerParam("nrounds", lower = 100, upper = 500), 
                         makeIntegerParam("max_depth", lower = 1, upper = 10), 
                         makeNumericParam("eta", lower = 0.1, upper = 0.5), 
                         makeNumericParam("lambda",lower = -1, upper = 0, trafo = function(x) 10^x)),
  resampling = makeResampleDesc("SpCV", iters = 5),
  control = makeTuneControlRandom(maxit = 200)
)
parallelStop()

# Evaluate tuning
generateHyperParsEffectData(tuned_reg, partial.dep = T) %>% plotHyperParsEffect(x = "iteration", y = "mse.test.mean", plot.type = "line", partial.dep.learn = xgb_reg_learner)
generateHyperParsEffectData(tuned_clas, partial.dep = T) %>% plotHyperParsEffect(x = "iteration", y = "mmce.test.mean", plot.type = "line", partial.dep.learn = xgb_reg_learner)

generateHyperParsEffectData(tuned_reg_std, partial.dep = T) %>% plotHyperParsEffect(x = "iteration", y = "mse.test.mean", plot.type = "line", partial.dep.learn = xgb_reg_learner)
generateHyperParsEffectData(tuned_clas_std, partial.dep = T) %>% plotHyperParsEffect(x = "iteration", y = "mmce.test.mean", plot.type = "line", partial.dep.learn = xgb_reg_learner)

generateHyperParsEffectData(tuned_reg_sc, partial.dep = T) %>% plotHyperParsEffect(x = "iteration", y = "mse.test.mean", plot.type = "line", partial.dep.learn = xgb_reg_learner)
generateHyperParsEffectData(tuned_clas_sc, partial.dep = T) %>% plotHyperParsEffect(x = "iteration", y = "mmce.test.mean", plot.type = "line", partial.dep.learn = xgb_reg_learner)

### 5.3 Build models ####

# models
xgb_regmodel <- train(setHyperPars(learner = xgb_reg_learner, par.vals = tuned_reg$x), train_reg)
xgb_clasmodel <- train(setHyperPars(learner = xgb_clas_learner, par.vals = tuned_clas$x), train_clas)

xgb_regmodel_std <- train(setHyperPars(learner = xgb_reg_learner, par.vals = tuned_reg_std$x), train_reg_std)
xgb_clasmodel_std <- train(setHyperPars(learner = xgb_clas_learner, par.vals = tuned_clas_std$x), train_clas_std)

xgb_regmodel_sc <- train(setHyperPars(learner = xgb_reg_learner, par.vals = tuned_reg_sc$x), train_reg_sc)
xgb_clasmodel_sc <- train(setHyperPars(learner = xgb_clas_learner, par.vals = tuned_clas_sc$x), train_clas_sc)

# feature importance

featimp <- bind_rows(getFeatureImportance(xgb_regmodel)$res %>% mutate(handling = "raw", task = "regression"),
                     getFeatureImportance(xgb_regmodel_std)$res %>% mutate(handling = "standardized", task = "regression"),
                     getFeatureImportance(xgb_regmodel_sc)$res %>% mutate(handling = "PCA", task = "regression"),
                     getFeatureImportance(xgb_clasmodel)$res %>% mutate(handling = "raw", task = "classification"),
                     getFeatureImportance(xgb_clasmodel_std)$res %>% mutate(handling = "standardized", task = "classification"),
                     getFeatureImportance(xgb_clasmodel_sc)$res %>% mutate(handling = "PCA", task = "classification")
) %>% mutate(handling = factor(handling, levels = c("raw", "standardized", "PCA")))

featimp %>% filter(importance > 0.05) %>% 
  ggplot(aes(x = variable, y = importance, fill = handling))+
  geom_bar(stat = "identity")+
  facet_grid(handling ~ task, scales = "free", space = "free")+
  coord_flip()+
  theme_bw()+
  theme(
    axis.title.x = element_text(size = 20, face = "bold"),
    axis.title.y = element_blank(),
    axis.text = element_text(size = 16),
    stip.text = element_text(size = 20, face = "bold"),
    legend.position = "none"
  )


### 5.4 Predictions #### 
pred <- predict(xgb_regmodel, task = train_reg)
pred_std <- predict(xgb_regmodel_std, task = train_reg_std)
pred_sc <- predict(xgb_regmodel_sc, task = train_reg_sc)

# Performance
performance(pred)
performance(pred_std)
performance(pred_sc)

### 5.5 XGBoost 2nd version ####
# units converted to Km and Km^2
# standardized vs. unstandardized data

#### 5.5.1 Training and Test sets ####
set.seed(1234)
inTrain <- createDataPartition(
  y = data4model_2nd$EQC_mean,
  p = 0.7,
  list = F
)

train_coords_2nd <-data4model_2nd[inTrain,] %>% select(X,Y) 
test_coords_2nd <-data4model_2nd[-inTrain,] %>% select(X,Y) 

train_reg_2nd <- data4model_2nd[inTrain,] %>% select(-EQC_mean, -X,-Y) %>% 
  createDummyFeatures(target = "MMI_mean", cols = c("GESAMT", "RIV_TYPE")) %>% makeRegrTask(data = ., target = "MMI_mean", coordinates = train_coords_2nd)
test_reg_2nd <- data4model_2nd[-inTrain,] %>% select(-EQC_mean, -X,-Y) %>% 
  createDummyFeatures(target = "MMI_mean", cols = c("GESAMT", "RIV_TYPE")) %>% makeRegrTask(data = ., target = "MMI_mean", coordinates = test_coords_2nd)

train_clas_2nd <- data4model_2nd[inTrain,] %>% select(-MMI_mean, -X,-Y) %>% 
  createDummyFeatures(target = "EQC_mean", cols = c("GESAMT", "RIV_TYPE")) %>% makeClassifTask(data = ., target = "EQC_mean", coordinates = train_coords_2nd)
test_clas_2nd <- data4model_2nd[-inTrain,] %>% select(-MMI_mean, -X,-Y) %>% 
  createDummyFeatures(target = "EQC_mean", cols = c("GESAMT", "RIV_TYPE")) %>% makeClassifTask(data = ., target = "EQC_mean", coordinates = test_coords_2nd)


# standardized data
train_coords_std_2nd <-data4model_std_2nd[inTrain,] %>% select(X,Y) 
test_coords_std_2nd <-data4model_std_2nd[-inTrain,] %>% select(X,Y) 

train_reg_std_2nd <- data4model_std_2nd[inTrain,] %>% select(-EQC_mean, -X,-Y) %>% 
  createDummyFeatures(target = "MMI_mean", cols = c("GESAMT", "RIV_TYPE")) %>% makeRegrTask(data = ., target = "MMI_mean", coordinates = train_coords_std_2nd)
test_reg_std_2nd <- data4model_std_2nd[-inTrain,] %>% select(-EQC_mean, -X,-Y) %>%  
  createDummyFeatures(target = "MMI_mean", cols = c("GESAMT", "RIV_TYPE")) %>% makeRegrTask(data = ., target = "MMI_mean", coordinates = test_coords_std_2nd)

train_clas_std_2nd <- data4model_std_2nd[inTrain,] %>% select(-MMI_mean, -X,-Y) %>% 
  createDummyFeatures(target = "EQC_mean", cols = c("GESAMT", "RIV_TYPE")) %>% makeClassifTask(data = ., target = "EQC_mean", coordinates = train_coords_std_2nd)
test_clas_std_2nd <- data4model_std_2nd[-inTrain,] %>% select(-MMI_mean, -X,-Y) %>% 
  createDummyFeatures(target = "EQC_mean", cols = c("GESAMT", "RIV_TYPE")) %>% makeClassifTask(data = ., target = "EQC_mean", coordinates = test_coords_std_2nd)


# create learners
xgb_reg_learner <- makeLearner(
  "regr.xgboost",
  predict.type = "response"
)

xgb_clas_learner <- makeLearner(
  "classif.xgboost",
  predict.type = "response"
)

#### 5.5.2 Hyper-parameter tuning ####
parallelStartSocket(detectCores()-1)
tuned_reg_2nd <- tuneParams(
  learner = xgb_reg_learner,
  task = train_reg_2nd,
  par.set = makeParamSet(makeIntegerParam("nrounds", lower = 100, upper = 500), 
                         makeIntegerParam("max_depth", lower = 1, upper = 10), 
                         makeNumericParam("eta", lower = 0.1, upper = 0.5), 
                         makeNumericParam("lambda",lower = -1, upper = 0, trafo = function(x) 10^x)),
  resampling = makeResampleDesc("SpCV", iters = 5),
  control = makeTuneControlRandom(maxit = 200)
)
parallelStop()

parallelStartSocket(detectCores()-1)
tuned_clas_2nd <- tuneParams(
  learner = xgb_clas_learner,
  task = train_clas_2nd,
  par.set = makeParamSet(makeIntegerParam("nrounds", lower = 100, upper = 500), 
                         makeIntegerParam("max_depth", lower = 1, upper = 10), 
                         makeNumericParam("eta", lower = 0.1, upper = 0.5), 
                         makeNumericParam("lambda",lower = -1, upper = 0, trafo = function(x) 10^x)),
  resampling = makeResampleDesc("SpCV", iters = 5),
  control = makeTuneControlRandom(maxit = 200)
)
parallelStop()

# standardized data
parallelStartSocket(detectCores()-1)
tuned_reg_std_2nd <- tuneParams(
  learner = xgb_reg_learner,
  task = train_reg_std_2nd,
  par.set = makeParamSet(makeIntegerParam("nrounds", lower = 100, upper = 500), 
                         makeIntegerParam("max_depth", lower = 1, upper = 10), 
                         makeNumericParam("eta", lower = 0.1, upper = 0.5), 
                         makeNumericParam("lambda",lower = -1, upper = 0, trafo = function(x) 10^x)),
  resampling = makeResampleDesc("SpCV", iters = 5),
  control = makeTuneControlRandom(maxit = 200)
)
parallelStop()

parallelStartSocket(detectCores()-1)
tuned_clas_std_2nd <- tuneParams(
  learner = xgb_clas_learner,
  task = train_clas_std_2nd,
  par.set = makeParamSet(makeIntegerParam("nrounds", lower = 100, upper = 500), 
                         makeIntegerParam("max_depth", lower = 1, upper = 10), 
                         makeNumericParam("eta", lower = 0.1, upper = 0.5), 
                         makeNumericParam("lambda",lower = -1, upper = 0, trafo = function(x) 10^x)),
  resampling = makeResampleDesc("SpCV", iters = 5),
  control = makeTuneControlRandom(maxit = 200)
)
parallelStop()


# Evaluate tuning
generateHyperParsEffectData(tuned_reg_2nd, partial.dep = T) %>% plotHyperParsEffect(x = "iteration", y = "mse.test.mean", plot.type = "line", partial.dep.learn = xgb_reg_learner)
generateHyperParsEffectData(tuned_clas_2nd, partial.dep = T) %>% plotHyperParsEffect(x = "iteration", y = "mmce.test.mean", plot.type = "line", partial.dep.learn = xgb_reg_learner)

generateHyperParsEffectData(tuned_reg_std_2nd, partial.dep = T) %>% plotHyperParsEffect(x = "iteration", y = "mse.test.mean", plot.type = "line", partial.dep.learn = xgb_reg_learner)
generateHyperParsEffectData(tuned_clas_std_2nd, partial.dep = T) %>% plotHyperParsEffect(x = "iteration", y = "mmce.test.mean", plot.type = "line", partial.dep.learn = xgb_reg_learner)

#### 5.5.3 Build models ####

# models
xgb_regmodel_2nd <- train(setHyperPars(learner = xgb_reg_learner, par.vals = tuned_reg_2nd$x), train_reg_2nd)
xgb_clasmodel_2nd <- train(setHyperPars(learner = xgb_clas_learner, par.vals = tuned_clas_2nd$x), train_clas_2nd)

xgb_regmodel_std_2nd <- train(setHyperPars(learner = xgb_reg_learner, par.vals = tuned_reg_std_2nd$x), train_reg_std_2nd)
xgb_clasmodel_std_2nd <- train(setHyperPars(learner = xgb_clas_learner, par.vals = tuned_clas_std_2nd$x), train_clas_std_2nd)

# feature importance
featimp_2nd <- bind_rows(getFeatureImportance(xgb_regmodel_2nd)$res %>% mutate(handling = "raw", task = "regression"),
                     getFeatureImportance(xgb_regmodel_std_2nd)$res %>% mutate(handling = "standardized", task = "regression"),
                     getFeatureImportance(xgb_clasmodel_2nd)$res %>% mutate(handling = "raw", task = "classification"),
                     getFeatureImportance(xgb_clasmodel_std_2nd)$res %>% mutate(handling = "standardized", task = "classification"),
) %>% mutate(handling = factor(handling, levels = c("raw", "standardized")))

featimp_2nd %>% filter(importance > 0.05) %>% 
  ggplot(aes(x = variable, y = importance, fill = handling))+
  geom_bar(stat = "identity")+
  facet_grid(handling ~ task, scales = "free", space = "free")+
  coord_flip()+
  theme_bw()+
  theme(
    axis.title.x = element_text(size = 20, face = "bold"),
    axis.title.y = element_blank(),
    axis.text = element_text(size = 16),
    stip.text = element_text(size = 20, face = "bold"),
    legend.position = "none"
  )

#### 5.5.4 Predictions #### 
pred_2nd <- predict(xgb_regmodel_2nd, task = train_reg_2nd)
pred_std_2nd <- predict(xgb_regmodel_std_2nd, task = train_reg_std_2nd)

# Performance
performance(pred_2nd)
performance(pred_std_2nd)

## 6. GWR ####

# data from predictions
pred_data <- train_reg$env$data %>% bind_cols(pred$data, train_reg$coordinates) %>% 
  select(response, Urban_prop, semi_Natural_prop, Urban, TOT_WASTE, cros_hw, culvert, X, Y) %>% st_as_sf(coords = c("X","Y"))

# standardized data
pred_data_std <- train_reg_std$env$data %>% bind_cols(pred_std$data, train_reg_std$coordinates) %>% 
  select(response, Urban_prop, semi_Natural_prop, semi_Natural, TOT_WASTE, culvert, X, Y) %>% st_as_sf(coords = c("X","Y"))

# PCA scores
pred_data_sc <- train_reg_sc$env$data %>% bind_cols(pred_sc$data, train_reg_sc$coordinates) %>% 
  select(response, pnt_comp1, pnt_comp2, lcc_comp1, dst_comp2, GESAMT.4, GESAMT.6, GESAMT.7, X, Y) %>% st_as_sf(coords = c("X","Y"))


# bandwidths
bwG <- bw.gwr(response ~ Urban_prop + semi_Natural_prop + Urban + TOT_WASTE + cros_hw + culvert, 
         data = pred_data,
         dMat = gw.dist(as.matrix(train_reg$coordinates)),
         adaptive = T
         )
# standardized data
bwG_std <- bw.gwr(response ~ Urban_prop + semi_Natural_prop + semi_Natural + TOT_WASTE + culvert, 
              data = pred_data_std,
              dMat = gw.dist(as.matrix(train_reg_std$coordinates)),
              adaptive = T
              )
# PCA scores
bwG_sc <- bw.gwr(response ~ pnt_comp1 + pnt_comp2 + lcc_comp1 + dst_comp2 + GESAMT.4 + GESAMT.6 + GESAMT.7,
              data = pred_data_sc,
              dMat = gw.dist(as.matrix(train_reg_sc$coordinates)),
              adaptive = T
              )

# models
gwr_model <- gwr.basic(
  response ~ Urban_prop + semi_Natural_prop + Urban + TOT_WASTE + cros_hw + culvert,
  data = pred_data, 
  bw = bwG,
  adaptive = T,
  dMat = gw.dist(as.matrix(train_reg$coordinates))
)

# standardized data
gwr_model_std <- gwr.basic(
  response ~ Urban_prop + semi_Natural_prop + semi_Natural + TOT_WASTE + culvert,
  data = pred_data_std, 
  bw = bwG_std,
  adaptive = T,
  dMat = gw.dist(as.matrix(train_reg_std$coordinates))
)

# PCA scores
gwr_model_sc <- gwr.basic(
  response ~ pnt_comp1 + pnt_comp2 + lcc_comp1 + dst_comp2 + GESAMT.4 + GESAMT.6 + GESAMT.7,
  data = pred_data_sc, 
  bw = bwG_sc,
  adaptive = T,
  dMat = gw.dist(as.matrix(train_reg_sc$coordinates))
)

### 6.1. GWR 2nd version ####
# data from predictions
pred_data_2nd <- train_reg_2nd$env$data %>% bind_cols(pred_2nd$data, train_reg_2nd$coordinates) %>% 
  select(response, semi_Natural, Urban, TOT_WASTE, cros_hw, culvert, X, Y) %>% st_as_sf(coords = c("X","Y"))

# standardized data
pred_data_std_2nd <- train_reg_std_2nd$env$data %>% bind_cols(pred_std_2nd$data, train_reg_std_2nd$coordinates) %>% 
  select(response, Urban, semi_Natural, culvert, cros_hw, X, Y) %>% st_as_sf(coords = c("X","Y"))

# bandwidths
bwG_2nd <- bw.gwr(response ~ semi_Natural + Urban + TOT_WASTE + cros_hw + culvert, 
              data = pred_data_2nd,
              dMat = gw.dist(as.matrix(train_reg_2nd$coordinates)),
              adaptive = T
)
# standardized data
bwG_std_2nd <- bw.gwr(response ~ Urban + semi_Natural + cros_hw + culvert, 
                  data = pred_data_std_2nd,
                  dMat = gw.dist(as.matrix(train_reg_std_2nd$coordinates)),
                  adaptive = T
)

# models
gwr_model_2nd <- gwr.basic(
  response ~ semi_Natural + Urban + TOT_WASTE + cros_hw + culvert,
  data = pred_data_2nd, 
  bw = bwG_2nd,
  adaptive = T,
  dMat = gw.dist(as.matrix(train_reg_2nd$coordinates))
)

# standardized data
gwr_model_std_2nd <- gwr.basic(
  response ~ Urban + semi_Natural + cros_hw + culvert,
  data = pred_data_std_2nd, 
  bw = bwG_std_2nd,
  adaptive = T,
  dMat = gw.dist(as.matrix(train_reg_std_2nd$coordinates))
)


## 7. Model Performance ####
### 7.1 XGBoost ####
# test predictions
pred_test_xgb <- predict(xgb_regmodel, task = test_reg)
pred_test_std_xgb <- predict(xgb_regmodel_std, task = test_reg_std)
pred_test_sc_xgb <- predict(xgb_regmodel_sc, task = test_reg_sc)

pred_test_cl_xgb <- predict(xgb_clasmodel, task = test_clas)
pred_test_cl_std_xgb <- predict(xgb_clasmodel_std, task = test_clas_std)
pred_test_cl_sc_xgb <- predict(xgb_clasmodel_sc, task = test_clas_sc)

# performance
performance(pred_test_xgb)
performance(pred_test_std_xgb)
performance(pred_test_sc_xgb)

performance(pred_test_cl_xgb)
performance(pred_test_cl_std_xgb)
performance(pred_test_cl_sc_xgb)

# partial dependence
generatePartialDependenceData(xgb_clasmodel, test_clas) %>% plotPartialDependence()

### 7.2 GWR & XGBoost ####

# wrangle test data
test_gwr <- test_reg$env$data %>% bind_cols(test_reg$coordinates) %>% st_as_sf(coords = c("X","Y"))
test_gwr_std <- test_reg_std$env$data %>% bind_cols(test_reg_std$coordinates) %>% st_as_sf(coords = c("X","Y"))
test_gwr_sc <- test_reg_sc$env$data %>% bind_cols(test_reg_sc$coordinates) %>% st_as_sf(coords = c("X","Y"))

# Predictions
pred_test_gwr <- gwr.predict(
  response ~ Urban_prop + semi_Natural_prop + Urban + TOT_WASTE + cros_hw + culvert,
  data = pred_data, 
  predictdata = test_gwr,
  bw = bwG,
  adaptive = T
  )

# standardized data
pred_test_gwr_std <- gwr.predict(
  response ~ Urban_prop + semi_Natural_prop + semi_Natural + TOT_WASTE + culvert,
  data = pred_data_std, 
  predictdata = test_gwr,
  bw = bwG,
  adaptive = T
)

# PCA scores
pred_test_gwr_sc <- gwr.predict(
  response ~ pnt_comp1 + pnt_comp2 + lcc_comp1 + dst_comp2 + GESAMT.4 + GESAMT.6 + GESAMT.7,
  data = pred_data_sc, 
  predictdata = test_gwr_sc,
  bw = bwG_sc,
  adaptive = T
)

# MSE
sum((test_gwr$MMI_mean-pred_test_gwr$SDF$prediction)^2)/nrow(pred_test_gwr$SDF)
sum((test_gwr_std$MMI_mean-pred_test_gwr_std$SDF$prediction)^2)/nrow(pred_test_gwr_std$SDF)
sum((test_gwr_sc$MMI_mean-pred_test_gwr_sc$SDF$prediction)^2)/nrow(pred_test_gwr_sc$SDF)

### 7.3 classification from reg models ####

# xgb
xgb_cl_from_rg <- pred_test_std_xgb$data %>% mutate(
  across(
    .cols = c(truth, response),
    .fns = ~ case_when(
      .x <= 0.2 ~ 5,
      .x > 0.2 & .x <= 0.4 ~ 4,
      .x > 0.4 & .x <= 0.6 ~ 3,
      .x > 0.6 ~ 2
    ),
    .names = "{.col}_class"
  )
) %>% mutate(
  match = case_when(truth_class == response_class ~ 0, .default = 1)
)

# gwr_xgb
gw_xgb_cl_from_rg <- pred_test_gwr$SDF %>% 
  bind_cols(pred_test_std_xgb$data) %>%
  select(truth, prediction) %>% 
  mutate(
  across(
    .cols = c(truth, prediction),
    .fns = ~ case_when(
      .x <= 0.2 ~ 5,
      .x > 0.2 & .x <= 0.4 ~ 4,
      .x > 0.4 & .x <= 0.6 ~ 3,
      .x > 0.6 ~ 2
    ),
    .names = "{.col}_class"
  )
) %>% mutate(
  match = case_when(truth_class == prediction_class ~ 0, .default = 1)
)

# mmce
sum(xgb_cl_from_rg$match)/length(xgb_cl_from_rg$match)
sum(gw_xgb_cl_from_rg$match)/length(gw_xgb_cl_from_rg$match)

### 7.4 Model performance 2nd version ####
#### 7.4.1 XGBoost ####

pred_test_xgb_2nd <- predict(xgb_regmodel_2nd, task = test_reg_2nd)
pred_test_std_xgb_2nd <- predict(xgb_regmodel_std_2nd, task = test_reg_std_2nd)

pred_test_cl_xgb_2nd <- predict(xgb_clasmodel_2nd, task = test_clas_2nd)
pred_test_cl_std_xgb_2nd <- predict(xgb_clasmodel_std_2nd, task = test_clas_std_2nd)

# performance
performance(pred_test_xgb_2nd)
performance(pred_test_std_xgb_2nd)

performance(pred_test_cl_xgb_2nd)
performance(pred_test_cl_std_xgb_2nd)

#### 7.4.2 GWR & XGBoost ####

# wrangle test data
test_gwr_2nd <- test_reg_2nd$env$data %>% bind_cols(test_reg_2nd$coordinates) %>% st_as_sf(coords = c("X","Y"))
test_gwr_std_2nd <- test_reg_std_2nd$env$data %>% bind_cols(test_reg_std_2nd$coordinates) %>% st_as_sf(coords = c("X","Y"))

# Predictions
pred_test_gwr_2nd <- gwr.predict(
  response ~ Urban + semi_Natural + TOT_WASTE + cros_hw + culvert,
  data = pred_data_2nd, 
  predictdata = test_gwr_2nd,
  bw = bwG_2nd,
  adaptive = T
)

# standardized data
pred_test_gwr_std_2nd <- gwr.predict(
  response ~ Urban + semi_Natural + cros_hw + culvert,
  data = pred_data_std_2nd, 
  predictdata = test_gwr_2nd,
  bw = bwG_2nd,
  adaptive = T
)

# MSE
sum((test_gwr_2nd$MMI_mean-pred_test_gwr_2nd$SDF$prediction)^2)/nrow(pred_test_gwr_2nd$SDF)
sum((test_gwr_std_2nd$MMI_mean-pred_test_gwr_std_2nd$SDF$prediction)^2)/nrow(pred_test_gwr_std_2nd$SDF)

#### 7.4.3 classification from reg models ####

# xgb
xgb_cl_from_rg_2nd <- pred_test_std_xgb_2nd$data %>% mutate(
  across(
    .cols = c(truth, response),
    .fns = ~ case_when(
      .x <= 0.2 ~ 5,
      .x > 0.2 & .x <= 0.4 ~ 4,
      .x > 0.4 & .x <= 0.6 ~ 3,
      .x > 0.6 ~ 2
    ),
    .names = "{.col}_class"
  )
) %>% mutate(
  match = case_when(truth_class == response_class ~ 0, .default = 1)
)

# mmce
sum(xgb_cl_from_rg_2nd$match)/length(xgb_cl_from_rg_2nd$match)