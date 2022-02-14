suppressPackageStartupMessages({
    library(BiocManager)
    library(ggthemes)
    library(purrr)
    library(tidyr)
    library(data.table)
    library(dplyr)
    library(glue)
    library(tibble)
    library(Matrix)
    library(ggplot2)
    library(patchwork)
    library(ggrepel)
    library(presto)
    library(scales)
    library(ComplexHeatmap)
    library(viridis)
    
    ## Parallel 
    library(future)
    library(furrr)

    ## New spatial libraries
    library(sf)
    library(concaveman)
    library(smoothr)
    library(stars)
    library(geojson)
    library(classInt)
})

fig.size <- function(h, w) {
    options(repr.plot.width=w, repr.plot.height=h)
}