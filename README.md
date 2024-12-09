README
================
Daniel Enns
2024-10-17

# Study outline

Freshwater ecosystems exists in very heterogenous landscapes, including
not only natural, but also anthropogenic “spatial features”[^1]. Can
these features, like barriers, wastewater treatment plants, road &
railroad crossings, etc. explain the ecological status of the freshwater
fauna? In this study, we try to answer this question combining various
spatial data sets and utilizing machine learning algorithms[^2].

# Methods

## Combining Spatial Data

For the spatial analysis all necessary shape and raster files are
projected into EPSG . The HLNUG stream network is a collection of 100 m
long line segments, each containing a stream ID and a segment ID. The
stream ID is designed in such a way, that it represents the one from the
subsequent stream in which the current flows into, followed by
additional numbers. The segment number increases towards the source of a
stream (Fig. 1). Each 100 m segment already contains information on the
habitats structural quality (1 - pristine, 7 - completly altered).

<figure>
<img src="README_files/figure-gfm/ID_example_figure-1.png"
alt="Fig. 1: Example for stream- and segment ID. The segment number increases in the direction of the source. The stream ID of smaller confluence resembels always that of the bigger stream, in which it flows into, plus additional numbers." />
<figcaption aria-hidden="true">Fig. 1: Example for stream- and segment
ID. The segment number increases in the direction of the source. The
stream ID of smaller confluence resembels always that of the bigger
stream, in which it flows into, plus additional numbers.</figcaption>
</figure>

The stream network can be dissolved in QGIS using the dissolve tool and
subsequently an initial landscape network (LSN) can build using the
`lines_to_lsn()` function from the `SSNbler` package, where the snap
tolerance of nodes is set to 3 m and the topology tolerance to 30 m
(note that the digitized flow needs to lead from source to outlet
points). The summary of the node error report reveals that the network
is highly erroneous. (For more detail on the LSN parameters and error
types see [Peterson & Pearse
2024](https://github.com/pet221/SSNbler/blob/main/inst/tutorials/Topology_Editing/QGIS/TopologyEditing_QGIS.pdf))

    ##     pointid            nodecat                            error      
    ##  Min.   :    2   Confluence: 9226   Converging Node          :   54  
    ##  1st Qu.: 5206   Outlet    :10167   Dangling Node            :19172  
    ##  Median :10520                      Downstream Divergence    :   75  
    ##  Mean   :10518                      Intersection Without Node:   88  
    ##  3rd Qu.:15808                      Unsnapped Node           :    4  
    ##  Max.   :21137                                                       
    ##  NA's   :9149                                                        
    ##             geom      
    ##  POINT        :19393  
    ##  epsg:25832   :    0  
    ##  +proj=utm ...:    0  
    ##                       
    ##                       
    ##                       
    ## 

Fixing topological errors in the entire network, manually and using
GRASS GIS tools, would take too much time, therefore using this network
for routing operations is not feasible. To quantify the number of point
features upstream of a respective bio sampling site, the existing stream
and section ID system can be used . The assignment works as follows for
each stressor type:

1.  count all features on the entire stream and its confluences

2.  subtract the number of features with equal or smaller segment number
    of the bio sampling site

In order to delineate upstream watersheds for each bio sampling site a
digital elevation model will be used to extract flow accumulation and
pointer grids, from which the watershed is calculated. Using the SRTM
GL1 30m DEM[^3] for Hessen has the drawback of creating inaccurate flow
accumulation and pointer grids in areas with low profiles (namely the
Rhine basin). The existing stream network could be burned into the DEM
beforehand to improve accuracy, however the network needs to be fixed of
all topological errors and the DEM is too large.

<figure>
<img src="README_files/figure-gfm/stream_size_figure-1.png"
alt="Fig. 2: streams split by their size (stream ID digits)" />
<figcaption aria-hidden="true">Fig. 2: streams split by their size
(stream ID digits)</figcaption>
</figure>

<figure>
<img src="README_files/figure-gfm/plot_mzb_stream_count-1.png"
alt="Fig. 3: number of invertebrate sampling sites located on different stream sizes" />
<figcaption aria-hidden="true">Fig. 3: number of invertebrate sampling
sites located on different stream sizes</figcaption>
</figure>

Finally, watersheds can be extracted using the `whitebox` package,
following the instructions from [Gannon,
2024](https://vt-hydroinformatics.github.io/Quarto_Book/15-Watershed-Delineation-and-Analysis.html).

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

<div id="ref-heß2023" class="csl-entry">

Heß, Sebastian, Delia Hof, Matthias Oetken, and Andrea Sundermann. 2023.
“Effects of Multiple Stressors on Benthic Invertebrates Using Water
Framework Directive Monitoring Data.” *Science of The Total Environment*
878 (June): 162952. <https://doi.org/10.1016/j.scitotenv.2023.162952>.

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

<div id="ref-opentopography2013a" class="csl-entry">

OpenTopography. 2013. “Shuttle Radar Topography Mission (SRTM) Global.”
<https://doi.org/10.5069/G9445JDF>.

</div>

<div id="ref-yu2020" class="csl-entry">

Yu, Hao, Arthur R. Cooper, and Dana M. Infante. 2020. “Improving Species
Distribution Model Predictive Accuracy Using Species Abundance:
Application with Boosted Regression Trees.” *Ecological Modelling* 432
(September): 109202. <https://doi.org/10.1016/j.ecolmodel.2020.109202>.

</div>

</div>

[^1]: (Carlisle, Falcone, and Meador 2008); (Marzin, Verdonschot, and
    Pont 2012); (Murphy and Davy-Bowker 2005)

[^2]: (Dedman et al. 2015); (Elith, Leathwick, and Hastie 2008); (Heß et
    al. 2023); (Yu, Cooper, and Infante 2020)

[^3]: (OpenTopography 2013)
