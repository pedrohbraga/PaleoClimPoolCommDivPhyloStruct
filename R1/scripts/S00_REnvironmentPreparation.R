#################################################################################################
#### Code to prepare your R environment for all the subsequent analyses done in this project ####
#                                                                                               #
# Author: Pedro Henrique Pereira Braga                                                          #
# Last Update: "2023-01-24"                                                                     #
#                                                                                               #
#################################################################################################

# Set seed #### 
# This seed was used for most parts of this project.

set.seed("213293")

# Install and load packages ####

# Begin with setting up non-CRAN packages

## try http:// if https:// URLs are not supported
## you may need to upgrade: biocLite("BiocUpgrade") 

if(!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ggtree")

BiocManager::install("limma")

# Load ggtree and wesanderson
library(ggtree); library(wesanderson)

## Or the development version from GitHub

# install.packages("devtools")

devtools::install_github("stamats/MKmisc")

devtools::install_github("cran/PhyloMeasures")

devtools::install_github("fawda123/WRTDStidal")

## Install wesanderson colour palettes
devtools::install_github("karthik/wesanderson")

# Check for needed packages, and install the missing ones
required.libraries <- c("ade4",
                        "RcppEigen", "snowfall", "rlist",
                        "ape", "picante", "phytools",
                        "vegan", "geiger", "Rcpp", "adephylo", "phylobase",
                        "lme4", "coda",
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
                        "ggrepel", "ggpubr", "wesanderson",
                        "latex2exp",
                        "naniar", "data.table",
                        "robust")

needed.libraries <- required.libraries[!(required.libraries %in% installed.packages()[,"Package"])]

# Install required libraries

if(length(needed.libraries)) install.packages(needed.libraries, Ncpus = 12)

# Load all required libraries

lapply(required.libraries, 
       require, 
       character.only = TRUE)

# Import utility functions ######

source("scripts/S00_fun_CommWeightedMeans.R")
source("scripts/S00_fun_sf.ses.phylostr.R")
source("scripts/S00_fun_brokenStick.selection.R")
source("scripts/S00_fun_make_grid_sf.R")
source("scripts/S00_fun_ggplot_theme_map.R")
source("scripts/S00_fun_match_phylo_comm.R")

