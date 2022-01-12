#################################################################################################
#### Code to prepare your R environment for all the subsequent analyses done in this project ####
#                                                                                               #
# Author: Pedro Henrique Pereira Braga                                                          #
# Last Update: "2021-09-21"                                                        #
#                                                                                               #
#################################################################################################

?set.seed()

# Begin with setting up non-CRAN packages

## try http:// if https:// URLs are not supported
## you may need to upgrade: biocLite("BiocUpgrade") 

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ggtree")
BiocManager::install("limma")

## Install wesanderson colour palettes
devtools::install_github("karthik/wesanderson")

# Load ggtree and wesanderson
library(ggtree); library(wesanderson)

## Or the development version from GitHub
# install.packages("devtools")
devtools::install_github("stamats/MKmisc")

# Check for needed packages, and install the missing ones
required.libraries <- c("ade4","subniche",
                        "RcppEigen", "snowfall", "rlist",
                        "letsR", 
                        "ape", "picante", "phytools",
                        "vegan", "geiger", "Rcpp", "adephylo", "phylobase",
                        "geiger", "hisse",
                        "lme4", "spgwr", "dglm", "coda",
                        "devtools", "rgeos", "dplyr", "gtable", "grid", "readxl",
                        "sf", "lwgeom", "sp","raster", "maptools", "splancs",
                        "ggplot2", "ggtree",  "viridis", "magrittr",
                        "rasterVis", "gridExtra", "PhyloMeasures",
                        "tidyverse", "wesanderson", "mosaic",
                        "broom", "tidyr", "scales",
                        "BiocManager", "ggpubr",
                        "kableExtra", "knitr",
                        "BAMMtools", "coda",
                        "rgdal", "tidyverse", "rbin", 
                        "quantreg", "WRTDStidal",
                        "ggrepel"
)

needed.libraries <- required.libraries[!(required.libraries %in% installed.packages()[,"Package"])]

if(length(needed.libraries)) install.packages(needed.libraries, Ncpus = 12)

# Load all required libraries at once

lapply(required.libraries, 
       require, 
       character.only = TRUE)

##### Import utility functions #############################################

source("scripts/SX_fun_CommWeightedMeans.R")
source("scripts/SX_fun_sf.ses.phylostr.R")
source("scripts/SX_fun_brokenStick.selection.R")

############################################################################
###### Utility function to create a hexagonal or square polygon grid #######
############################################################################

make_grid_sf <- function(x, type, cell_width, cell_area, clip = TRUE, propCover = 0.5) {
  if (!type %in% c("square", "hexagonal")) {
    stop("Type must be either 'square' or 'hexagonal'")
  }
  
  if (missing(cell_width)) {
    if (missing(cell_area)) {
      stop("Must provide cell_width or cell_area")
    } else {
      if (type == "square") {
        cell_width <- sqrt(cell_area)
      } else if (type == "hexagonal") {
        cell_width <- sqrt(2 * cell_area / sqrt(3))
      }
    }
  }
  # buffered extent of study area to define cells over
  ext <- as(extent(x) + cell_width, "SpatialPolygons")
  projection(ext) <- projection(x)
  # generate grid
  if (type == "square") {
    g <- raster(ext, resolution = cell_width)
    g <- as(g, "SpatialPolygons")
  } else if (type == "hexagonal") {
    
    # create a buffer to sample points
    ext.large <- gBuffer(x, byid = T, width = cell_width/2)
    
    # generate array of hexagon centers
    g <- spsample(ext.large, type = "hexagonal", cellsize = cell_width, offset = c(0, 0))
    
    # convert center points to hexagons
    g <- HexPoints2SpatialPolygons(g, dx = cell_width)
  }
  
  # clip to boundary of study area
  if (clip) {
    g.SF <- sf::st_as_sf(g) # convert x object to sf compatible
    x.SF <- sf::st_as_sf(x) # convert g to sf compatible
    
    g.Ins.SF <- sf::st_intersection(g.SF, x.SF) # intersect both polygons
    # keep hexagons that cover a certain percent of proportinal area cover
    if(propCover > 0) { 
      g.Area.SF <- st_area(g.SF) # 
      g.Ins.Area.SF <- st_area(g.Ins.SF)
      
      g.Prop.SF <- g.Ins.Area.SF/g.Area.SF
      
      g <- g.SF[as.numeric(g.Prop.SF) >= 0.3,]
    }
  } else {
    g <- g[x, ]
  }
  
  # clean up feature IDs
  row.names(g) <- as.character(1:nrow(g))
  return(g)
}


#############################################################################################
#### Function to match (and subset) species between phylogenetic and community data sets ####
#### Modified from the picante package                                                   ####   
#############################################################################################

match.phylo.comm. <- function(phy, comm) {
  if (!(is.data.frame(comm) | is.matrix(comm))) {
    stop("Community data should be a data.frame or matrix with samples in rows and taxa in columns")
  }
  res <- list()
  phytaxa <- phy$tip.label
  commtaxa <- colnames(comm)
  if (is.null(commtaxa)) {
    stop("Community data set lacks taxa (column) names, these are required to match phylogeny and community data")
  }
  if (!all(commtaxa %in% phytaxa)) {
    print("Dropping taxa from the community because they are not present in the phylogeny:")
    print(setdiff(commtaxa, phytaxa))
    comm <- comm[, intersect(commtaxa, phytaxa), drop = FALSE]
    commtaxa <- colnames(comm)
  }
  if (any(!(phytaxa %in% commtaxa))) {
    print("Dropping tips from the tree because they are not present in the community data:")
    print(setdiff(phytaxa, commtaxa))
    res$phy <- prune.sample(comm, phy)
  }
  else {
    res$phy <- phy
  }
  res$comm <- comm[, res$phy$tip.label, drop = FALSE]
  return(res)
}

####################################
####### Graphical parameters #######
####################################

# General theme for the map

theme_map <- function(...) {
  theme_minimal() +
    theme(
      text = element_text(family = "Ubuntu Regular", color = "#22211d"),
      axis.line = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      # panel.grid.minor = element_line(color = "#ebebe5", size = 0.2),
      panel.grid.major = element_line(color = "#ebebe5", size = 0.2),
      panel.grid.minor = element_blank(),
      plot.background = element_blank(), 
      panel.background = element_blank(), 
      legend.background = element_blank(),
      panel.border = element_blank(),
      ...
    )
}

