# README
Daniel Enns
2024-10-17

# Study outline

Freshwater ecosystems exists in very heterogenous landscapes, including
not only natural, but also anthropogenic “spatial features”(Carlisle,
Falcone, and Meador 2008; Marzin, Verdonschot, and Pont 2012; Murphy and
Davy-Bowker 2005). Can these features, like barriers, wastewater
treatment plants, road & railroad crossings, etc. explain the ecological
status of the freshwater fauna? In this study, we try to answer this
question combining various spatial data sets and utilizing machine
learning algorithms (Dedman et al. 2015; Elith, Leathwick, and Hastie
2008; Heß et al. 2023; Yu, Cooper, and Infante 2020).

# Methods

## 1. Library

The following packages are needed for the analysis:

``` r
# general
library(readxl)
library(tidyverse)
library(doParallel)

# spatial analysis
library(stars)
library(sf)
library(sfnetworks)
library(tidygraph)
library(nngeo)
library(geodata)
library(whitebox)
library(osmextract)

# modeling
library(caret)
library(mgcv)
library(xgboost)
library(pdp)
```

## 2. Combining Spatial Data

### 2.1 Point features

For the spatial analysis all necessary shape and raster files are
projected into EPSG 25832.

Highway and railroad networks can be read in from OpenStreetMap pbf
files via the `osmextract` package. These can be intersected with the
stream network to create point features. In the osm highway and railroad
networks each lane / track is represented as an individual line feature,
so that multiple crossings are present in locations with high lane /
track density.

``` r
# Extract highway and railroad networks
transport_net <- oe_read("./hessen-latest.osm.pbf", extra_tags = "railway") %>%
  filter(highway == "motorway" | railway == "rail") %>% 
  st_transform(st_crs(stream_net))

# Intersect with stream network
crossings <- st_intersection(stream_net, transport_net) %>% 
  st_cast("POINT")
```

### 2.2 Watershed delineation

In order to delineate upstream watersheds for each bio sampling site the
SRTM GL1 30m (OpenTopography 2013) digital elevation model will be used
to extract flow accumulation and pointer grids, from which the
watersheds are calculated. However, the delineation algorithms struggles
in regions with low profiles (namely the Rhine basin), creating an
inaccurate stream raster and watersheds. To improve accuracy, the
existing stream network is burned into the DEM following these steps in
QGIS:

1.  Reproject raster network to stream CRS (EPSG: 25832)

2.  dissolve the stream network

3.  Run the Grass r.carve algorithm with a stream width of 60 m, depth
    of 3 m and with no flat areas in direction of stream flow allowed
    (Note: the algorithm does not work on latitude-longitude CRS, see
    [r.carve GRASS GIS
    manual](https://grass.osgeo.org/grass-stable/manuals/r.carve.html))

The algorithm will take some time to run (approximately 1 hour and 18
minutes with an Intel core i7, 16 GB RAM). Finally, watersheds can be
extracted using the `whitebox` package, following the instructions from
[Gannon,
2024](https://vt-hydroinformatics.github.io/Quarto_Book/15-Watershed-Delineation-and-Analysis.html).

``` r
# Breach and fill pits in raster
wbt_breach_depressions_least_cost(
  dem = "./SRTMGL1_30m_2px_burned.tif", 
  output = "./srtm_breach.tif", 
  dist = 7,
  fill = T
)

# Create flow accumulation & pointer grids
wbt_d8_flow_accumulation(
  input = "./srtm_breach.tif",
  output = "./D8FA.tif"
)

wbt_d8_pointer(
  dem = "./srtm_breach.tif",
  output = "./D8pointer.tif"
)

# Extract streams
wbt_extract_streams(
  flow_accum = "./D8FA.tif",
  output = "./stream_raster.tif",
  threshold = 700
)

# Snap points to stream raster
wbt_jenson_snap_pour_points(
  pour_pts = "./MZB.shp",
  streams = "./stream_raster.tif",
  output = "./MZB_snap.shp",
  snap_dist = 300
)

# Delineate watersheds
wbt_watershed(
  d8_pntr = "./D8pointer.tif",
  pour_pts = "./MZB_snap.shp",
  output = "./watersheds.tif"
)
```

Figure 1 demonstrates the structure of the watersheds, which were build
consecutively between multiple pour points.

![Figure 1: Example watershed. Each watershed is build consecutively
between specified
pourpoints.](README_files/figure-commonmark/watershed_example_figure-1.png)

The created watershed raster can be read in and converted to polygon
features, which subsequently be used to intersect the [CORINE Land
Cover](https://land.copernicus.eu/en/products/corine-land-cover)
polygons. The area of each intersection can be calculated using the
`st_area()` function.

``` r
# Read in Watersheds and convert to polygon 
ws_poly <- read_stars("./watersheds.tif") %>% 
  st_as_sf(as_points = F, merge = T)

# Create dataframe with land cover areas
landcover_ws <- st_intersection(landcover, ws_poly) %>% 
  mutate(area = st_area(.)) %>% st_drop_geometry() %>% 
  as.data.frame() %>% filter(!is.na(type))

# Summarize land cover types in each watershed
landcover_sum <- landcover_ws %>% group_by(watersheds.tif, type) %>% 
  summarize(area = sum(area)) %>% 
  pivot_wider(names_from = type, values_from = area)

# Combine watersheds with summarized land cover
watersheds_lc <- left_join(ws_poly, landcover_sum, by = "watersheds.tif")
```

### 2.3 Network building

The Hessian stream network is a collection of 100 m long line segments,
each containing a stream ID and a segment ID. The stream ID is designed
in such a way, that it consecutively builds upon the ID of a previous
stream segment, adding numbers to it, whenever the stream splits in
upstream direction. The segment number increases towards the source of a
stream (Fig. 2). Each 100 m segment already contains information on the
in-stream and surrounding habitats structural quality (1 - pristine, 7 -
completly altered).

![Figure 2: Example of stream- and segment ID system. The segment number
increases in the direction of the source. The stream ID of smaller
confluence resembels always that of the bigger stream, in which it flows
into, plus additional
numbers.](README_files/figure-commonmark/ID_example_figure-1.png)

The `sfnetwork` package provides many useful functions for network
analysis and routing operations, where the main function
`as_sfnetwork()` builds a network by connecting aligning nodes. An
initial build reveals that many stream segments visually seem connected,
but are separated from their nearest node or edge by a few meters. Thus,
to build a useful network for routing the precision of the stream
network `sf` object needs to be lowered (for more information see
Details on `st_set_precision()` and `st_as_binary()` help page). To what
precision to round depends, in general, on the maximum distance of
dangling nodes and if the network can be build by connecting nodes with
edges, or just by connecting nodes. The precision will be set to 70
meters (be aware of the object CRS units!).

In the resulting network, the right streams are connected (Fig. 3) and
it can be used for routing operations. Performance of routing can be
greatly increased by simplifying the network and removing pseudo-nodes
(more info
[here](https://luukvdmeer.github.io/sfnetworks/articles/sfn02_preprocess_clean.html)).

``` r
# Lower precision to 50 meters
str_net_50m <- st_set_precision(stream_net, set_units(70, "m"))

# Create sf_network object, simplify and smooth pseudo-nodes
network <- as_sfnetwork(str_net_50m) %>% convert(to_spatial_simple) %>% 
  convert(to_spatial_smooth)
```

![Figure 3: Example for how precision changes the connectedness of the
network. The numbers indicate from which node a stream
originates.](README_files/figure-commonmark/Precision_example_figure-1.png)

Sampling sites and point stressors can be blend in as nodes into the
created network.

``` r
# Blend in sampling sites, wastewater treatment plants and dams
network_blend <- st_network_blend(network, mzb) %>% st_network_blend(.,wwtp) %>% 
  st_network_blend(.,dams)
```

## 3. Spatial Data quantification

For the quantification of point- and polygon features and their
attributes I wrote two custom functions, `st_shift()` and
`upstream_summarize()`. The first function is required within the latter
and shifts a set of point geometries a proportional distance towards
another specified point. The `upstream_summarize()` function builds a
sub-network from a given point to all its vertices (stream sources) and
extracts the node data, from which it counts the number of specified
nodes and sums up specified attribute values. Further, it can calculate
the minimal, maximal, or average network distance between the start
point and a set of nodes specified by their IDs with the `IDs` argument.
If no such nodes are present in the sub-network, `inf` or `NaN` values
are returned. Additionally the Shreve stream magnitude, defined as the
sum of sources upstream from a specifiedy node, can be calculated. The
function can sum up attribute values of provided polygons the following
way: First, the set of nodes present in the sub-network are shifted
towards their centroid by 0.1% of their length, to avoid including
adjacent polygons. Then, it creates a filtering mask by selecting all
polygons which are touched by the shifted nodes and fills in ‘holes’ in
the set of polygons. By setting a threshold, the mask can be shrunken to
avoid selecting adjacent polygons. This mask is finally used to select
all polygons present in the sub-network, from which their specified
attributes are summarized. The function contains the following set of
arguments:

- `net` : The complete network with blended in points of interest
- `start` : Row name of node from which to rout upstream
- `node_cols` : Names of attribute columns to summarize, present in
  nodes
- `IDs` : ID columns for sets of nodes, which should be unique in each
  set
- `area` : `sf` object with only polygon geometries
- `area_cols` : names of attribute columns to summarize, present in
  `area`
- `dist` : Either `'min'`, `'max'`, `'mean'`, or `'all'`. Calculates
  distances between start and nodes, specified by `IDs`.
- `threshold` : Value (in polygon CRS units) by which the polygon mask
  should be shrunken.
- `Shreve` : logical, should Shreve stream magnitud be calculate?

Before this function can be used, the points from which the function
should rout must be extracted as nodes from the network and their row
name must be saved as a new column.

``` r
# Extract nodes of invertebrate sampling sites and add row names
mzb_nodes <- st_as_sf(network_blend, "nodes") %>% filter(!is.na(ID_SITE)) %>% 
  mutate(ID_NODE = row.names(nodes)[with(nodes, !is.na(ID_SITE))])
```

`upstream_summarize()` can be used in combination with rowwise() and
mutate() to perform the action over multiple points. Depending on the
number of nodes within the entire network, this task can take quite a
lot of time. Luckily this task can be parallelized to save a lot of
time.

``` r
# set up cluster
num_cores <- (detectCores()-1)
num_cores %>% makeCluster() %>% registerDoParallel()

# split data into chunks
data_chunks <- mzb_nodes %>% 
  split(., ceiling(seq_along(row_number(.)) / (length(row_number(.)) %/% num_cores)))

# Apply upstream_summarize parallel and row-wise over all sampling site nodes 
mzb_data_complete <- foreach(chunk = data_chunks, .combine = rbind, .packages = c("dplyr", "dtplyr","sf","sfnetworks","tidygraph","nngeo")) %dopar% {
  chunk %>% rowwise() %>% 
  mutate(upstream_summarize(
    net = network_blend,
    start = ID_NODE,
    node_cols = c("population_equivalents", "dam_discharge"),
    IDs = c("ID_WWTP", "ID_DAMS"),
    area = ws_clc,
    area_cols = c("Agriculture", "Urban", "semi-Natural"),
    dist = T,
    threshold = 30)
    )
}
```

## 4. Data analysis and handling

### 4.1 EDA

Before modeling it is good practice to get a sens of the target and
explanatory variables. Here, multiple models should perform two
different tasks, regression and multi-class classification. The
respective target variables are the multi metric index (MMI) and the
resulting ecological quality class (EQC). The MMI is a summarization of
multiple core metrics, compared to reference conditions (for more
details see (Hering et al. 2006)). It ranges from zero to one,
representing bad to good condition respectively.

![Figure 4: Histogram of the site averaged multi metric
index.](README_files/figure-commonmark/MMI_histogram-1.png)

![Figure 5: Histogram of the site averaged ecological quality
class.](README_files/figure-commonmark/EQC_histogram-1.png)

From the EQC histogram (Fig. 5) a strong class imbalance is observable.
To counter this and increase model performance, the data for the
minority class could be oversampled using \[…\]. Since the classes are
of ordered categorical, the minority could be merged with an adjacent
class, although that would mean a loss of information.

Some covariates follow the same method of quantification and are
therefore auto correlated by default. Sites at bigger reaches will have
more point stresors in the upstream watershed compared to sites further
upstream, therefore the different stressor covariates are automatically
correlated. The correlation can be displayed in a correlation plot (Fig.
6).

![Figure 6: Correlation plot of
covariates.](README_files/figure-commonmark/autocorrelation_plot-1.png)

Although tree ensemble learning algorithms are robust to auto-correlated
covariates, it is best practice to either drop some related covariates
or reduce the dimensionallity and use PCA scores as new variables.

### 4.2 Spatial autocorrelation

Toblers first law states that *“everything is related to everything
else, but near things are more related than distant things”*. This
degree can be quantified by the Moranes Index, which can be calculated
with distance-based weight matrix for points.

``` r
# k-nearest neighbours 
suppressPackageStartupMessages(require(deldir)) data_nb_knn <- knn2nb(knearneigh(data_all,k = 1))  
# distance neighbours 
dsts <- unlist(nbdists(data_nb_knn, data_all)) summary(dsts) max_1nn <- max(dsts)  data_nb_dst <- dnearneigh(data_all, d1 = 0, d2 = 0.75*max_1nn)  data_lw <- nb2listw(data_nb_dst, style = "W", zero.policy = T)  
# Moran statistics 
moran.test(data_all$MMI_mean, data_lw, alternative = "greater") moran.plot(data_all$MMI_mean, data_lw) 
```

       Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
        1.0   332.7   840.2  1197.1  1663.7  8718.0 


        Moran I test under randomisation

    data:  data_all$MMI_mean  
    weights: data_lw  
    n reduced by no-neighbour observations  

    Moran I statistic standard deviate = 39.571, p-value < 2.2e-16
    alternative hypothesis: greater
    sample estimates:
    Moran I statistic       Expectation          Variance 
         0.4357911968     -0.0005844535      0.0001216123 

The index shows a significant deviation from the expectation, signifying
strong spatial autocorrelation. This can also be visualized with a
Moranes plot.

![Figure 6: Moranes plot for the site averaged multi metric
index.](README_files/figure-commonmark/moranes_plot-1.png)

Consequently, a geographically weighted XGBoost model should be
calculated following \[Li (2019)\](Schimohr, Doebler, and Scheiner
2022), by using the geo-referenced predictions from the XGBoost model as
inputs for the geographically wheighted regression (GWR) model. Further,
spatial heterogenity should be kept in mind for model tuning and cross
validation (Schratz et al. 2019).

### 4.3 Data cleaning & feature engeneering

Tree ensemble learning algorithms can handle covariates of different
types and are robust to auto-correlated covariates. To see how the data
cleaning and feature engeneering influences the model performance, the
same data was prepared in three different ways. Modelling was performed
on the raw data, on data with standardized covariates and on scores
derived from three separate PCAs performed on point stressor, point
stressor distance and land use variables. Beforehand, distance variables
where transformed in two steps:

1.  A five meter tolerance was added to distances with zero meters and
    missing values were replaced by zero, representing absence of
    features
2.  For each stressor variable distances were

## 5. Modeling

### 5.1 Ensemble machine learning with XGBoost

We first separate the available data into a training and test set by a
ratio of 80% to 20% respectively.

``` r
set.seed(1234)

# create split
inTrain <- createDataPartition(
  y = data4model$EQC_mean,
  p = 0.8,
  list = F
)

# split data
```

### 5.2 GW-XGBoost

### 5.3 Model performance

The performance of regression models can be evaluated by comparison of
the mean squared errors. The classification models can be compared using
the missclassification rate. The predicted values from the regression
models can be translated to the EQC, which can also be evaluated via the
misclassification rate.

``` r
```

The results in Table XX reveal that model XX is the best performing
regression model while model XX is the best classifiaction model.
Further, the geographically wheighted

## 6. Data availability

Data on WFD invertebrate sampling and the Hessian stream network were
kindly provided by the Hessian state office for nature, environment and
geology (HLNUG). Data on wastewater treatment plants, stormwater
overflows are available for download from the
[WRRL-viewer](https://wrrl.hessen.de/mapapps/resources/apps/wrrl/index.html?lang=de).
Further, data on hydrological barriers are available from the [Amber
Barrier Atlas](https://amber.international/european-barrier-atlas/).
OpenStreetMap data can be downloaded as pbf files at
[geofabrik](https://download.geofabrik.de/). CORINE land cover
shapefiles and rasters are available at [Copernicus Land Monitoring
Service](https://land.copernicus.eu/en/products/corine-land-cover). The
SRTM GL1 30m digital elevation model can be downloaded from
[OpenTopography](https://portal.opentopography.org/raster?opentopoID=OTSRTM.082015.4326.1).

## 7. References

<div id="refs" class="references csl-bib-body hanging-indent"
entry-spacing="0">

<div id="ref-carlisle2008" class="csl-entry">

Carlisle, Daren M., James Falcone, and Michael R. Meador. 2008.
“Predicting the Biological Condition of Streams: Use of Geospatial
Indicators of Natural and Anthropogenic Characteristics of Watersheds.”
*Environmental Monitoring and Assessment* 151 (1-4): 143–60.
<https://doi.org/10.1007/s10661-008-0256-z>.

</div>

<div id="ref-dedman2015" class="csl-entry">

Dedman, Simon, Rick Officer, Deirdre Brophy, Maurice Clarke, and David
G. Reid. 2015. “Modelling Abundance Hotspots for Data-Poor Irish Sea
Rays.” *Ecological Modelling* 312 (September): 77–90.
<https://doi.org/10.1016/j.ecolmodel.2015.05.010>.

</div>

<div id="ref-elith2008" class="csl-entry">

Elith, J., J. R. Leathwick, and T. Hastie. 2008. “A Working Guide to
Boosted Regression Trees.” *Journal of Animal Ecology* 77 (4): 802–13.
<https://doi.org/10.1111/j.1365-2656.2008.01390.x>.

</div>

<div id="ref-hering2006" class="csl-entry">

Hering, Daniel, Christian K. Feld, Otto Moog, and Thomas Ofenböck. 2006.
“Cook Book for the Development of a Multimetric Index for Biological
Condition of Aquatic Ecosystems: Experiences from the European AQEM and
STAR Projects and Related Initiatives.” In, 311–24. Springer
Netherlands. <https://doi.org/10.1007/978-1-4020-5493-8_22>.

</div>

<div id="ref-heß2023" class="csl-entry">

Heß, Sebastian, Delia Hof, Matthias Oetken, and Andrea Sundermann. 2023.
“Effects of Multiple Stressors on Benthic Invertebrates Using Water
Framework Directive Monitoring Data.” *Science of The Total Environment*
878 (June): 162952. <https://doi.org/10.1016/j.scitotenv.2023.162952>.

</div>

<div id="ref-li2019" class="csl-entry">

Li, Lianfa. 2019. “Geographically Weighted Machine Learning and
Downscaling for High-Resolution Spatiotemporal Estimations of Wind
Speed.” *Remote Sensing* 11 (11): 1378.
<https://doi.org/10.3390/rs11111378>.

</div>

<div id="ref-marzin2012" class="csl-entry">

Marzin, Anahita, Piet F. M. Verdonschot, and Didier Pont. 2012. “The
Relative Influence of Catchment, Riparian Corridor, and Reach-Scale
Anthropogenic Pressures on Fish and Macroinvertebrate Assemblages in
French Rivers.” *Hydrobiologia* 704 (1): 375–88.
<https://doi.org/10.1007/s10750-012-1254-2>.

</div>

<div id="ref-murphy2005" class="csl-entry">

Murphy, John F., and John Davy-Bowker. 2005. “Spatial Structure in Lotic
Macroinvertebrate Communities in England and Wales: Relationship with
Physical, Chemical and Anthropogenic Stress Variables.” *Hydrobiologia*
534 (1-3): 151–64. <https://doi.org/10.1007/s10750-004-1451-8>.

</div>

<div id="ref-opentopography2013" class="csl-entry">

OpenTopography. 2013. “Shuttle Radar Topography Mission (SRTM) Global.”
<https://doi.org/10.5069/G9445JDF>.

</div>

<div id="ref-schimohr2022" class="csl-entry">

Schimohr, Katja, Philipp Doebler, and Joachim Scheiner. 2022.
“Prediction of Bike-Sharing Trip Counts: Comparing Parametric Spatial
Regression Models to a Geographically Weighted XGBoost Algorithm.”
*Geographical Analysis* 55 (4): 651–84.
<https://doi.org/10.1111/gean.12354>.

</div>

<div id="ref-schratz2019" class="csl-entry">

Schratz, Patrick, Jannes Muenchow, Eugenia Iturritxa, Jakob Richter, and
Alexander Brenning. 2019. “Hyperparameter Tuning and Performance
Assessment of Statistical and Machine-Learning Algorithms Using Spatial
Data.” *Ecological Modelling* 406 (August): 109–20.
<https://doi.org/10.1016/j.ecolmodel.2019.06.002>.

</div>

<div id="ref-yu2020" class="csl-entry">

Yu, Hao, Arthur R. Cooper, and Dana M. Infante. 2020. “Improving Species
Distribution Model Predictive Accuracy Using Species Abundance:
Application with Boosted Regression Trees.” *Ecological Modelling* 432
(September): 109202. <https://doi.org/10.1016/j.ecolmodel.2020.109202>.

</div>

</div>
