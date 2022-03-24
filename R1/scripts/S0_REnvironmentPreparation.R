#################################################################################################
#### Code to prepare your R environment for all the subsequent analyses done in this project ####
#                                                                                               #
# Author: Pedro Henrique Pereira Braga                                                          #
# Last Update: "2021-09-21"                                                                     #
#                                                                                               #
#################################################################################################

set.seed("213293")

# Begin with setting up non-CRAN packages

## try http:// if https:// URLs are not supported
## you may need to upgrade: biocLite("BiocUpgrade") 

if(!requireNamespace("BiocManager", quietly = TRUE))
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
                        "ggrepel",
                        "naniar", "data.table"
)

needed.libraries <- required.libraries[!(required.libraries %in% installed.packages()[,"Package"])]

if(length(needed.libraries)) install.packages(needed.libraries, Ncpus = 12)

# Load all required libraries at once

lapply(required.libraries, 
       require, 
       character.only = TRUE)

##### Import utility functions #############################################

source("scripts/SX_fun_CommWeightedMeans.R")
source("scripts/SX_fun_brokenStick.selection.R")
source("scripts/SX_fun_make_grid_sf.R")
source("scripts/SX_fun_ggplot_theme_map.R")
source("scripts/SX_fun_match_phylo_comm.R")
source("scripts/SX_fun_ses.opt.rarefaction.phylostr.R")
source("scripts/SX_fun_ses.phylostr.query.sf.R")

