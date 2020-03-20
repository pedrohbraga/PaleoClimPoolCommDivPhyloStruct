#################################################################################################
#### Code to prepare your R environment for all the subsequent analyses done in this project ####
#                                                                                               #
# Author: Pedro Henrique Pereira Braga                                                          #
# Last Update: "2019-01-30 15:29:30 EST"                                                        #
#                                                                                               #
#################################################################################################

# Check for needed packages, and install the missing ones
required.libraries <- c("ade4","subniche",
                        "RcppEigen", "snowfall", "rlist",
                        "letsR", 
                        "ape", "picante", "phytools",
                        "vegan", "geiger", "Rcpp", "adephylo", "phylobase",
                        "geiger", "hisse",
                        "lme4", "spgwr", "dglm", "coda",
                        "devtools", "rgeos", "dplyr", "gtable", "grid", "readxl",
                        "sp","raster", "rgdal", "maptools", "splancs",
                        "ggplot2", "ggtree",  "viridis", "magrittr",
                        "rasterVis", "gridExtra", "PhyloMeasures",
                        "tidyverse", "wesanderson", "mosaic", "velox",
                        "broom", "tidyr", "MKmisc",
                        "BiocManager", "ggpubr" # Have to be loaded at last
)
needed.libraries <- required.libraries[!(required.libraries %in% installed.packages()[,"Package"])]

if(length(needed.libraries)) install.packages(needed.libraries)

# Load all required libraries at once
lapply(required.libraries, 
       require, 
       character.only = TRUE)

### Install and load packages not available in CRAN

## try http:// if https:// URLs are not supported
## you may need to upgrade: biocLite("BiocUpgrade") 

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ggtree")
BiocManager::install("limma")

## Install wesanderson color palettes
devtools::install_github("karthik/wesanderson")

# Load ggtree and wesanderson

library(ggtree); library(wesanderson)

## Or the development version from GitHub
# install.packages("devtools")
devtools::install_github("stamats/MKmisc")

##### Import utility functions #############################################

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

###############################################################
### Modified functions to allow the computation in parallel ###
###############################################################

#### Mean pairwise distances (ses.MPD), heavily based on picante's original function

sf.ses.mpd <- function(samp, dis, abundance.weighted = FALSE, runs = 999, cores = 4){
  dis <- as.matrix(dis) # coerce's the distance matrix to matrix
  
  # calculate the observed mean phylogenetic distance in each observing unit
  mpd.obs <- mpd(samp, dis, abundance.weighted = abundance.weighted) 
  
  ## build null models in parallel
  sfInit(parallel = TRUE, 
         cpus = cores, 
         slaveOutfile="logAnalysis.txt")
  
  # load libraries in parallel
  sfLibrary(picante)
  sfLibrary(base)
  sfLibrary(mosaic)
  sfLibrary(snowfall)
  
  # export objects 
  sfExport('samp', 'dis', 'abundance.weighted', 'runs')
  
  # calculate expected mean phylogenetic distances
  mpd.rand <- t(sfSapply(1:runs, 
                         function(i){ 
                           sfCat(paste("Iteration ", i), sep="\n")
                           mpd(samp, taxaShuffle(dis),
                               abundance.weighted = abundance.weighted)
                         }))
  
  # 
  mpd.rand.mean <- sfApply(x = mpd.rand, margin = 2, fun = mean, 
                           na.rm = TRUE)
  
  mpd.rand.sd <- sfApply(x = mpd.rand, margin = 2, fun = sd, 
                         na.rm = TRUE)
  
  mpd.obs.rank <- sfApply(x = rbind(mpd.obs, mpd.rand), margin = 2, 
                          fun = rank)[1, ]
  
  sfStop()
  
  ###
  
  mpd.obs.z <- (mpd.obs - mpd.rand.mean)/mpd.rand.sd
  
  mpd.obs.rank <- ifelse(is.na(mpd.rand.mean), NA, mpd.obs.rank)
  
  data.frame(ntaxa = specnumber(samp), 
             mpd.obs, 
             mpd.rand.mean, 
             mpd.rand.sd, 
             mpd.obs.rank, 
             mpd.obs.z, 
             mpd.obs.p = mpd.obs.rank/(runs + 1), 
             runs = runs, 
             row.names = row.names(samp))
}

sf.ses.mntd <- function(samp, dis, abundance.weighted = FALSE, runs = 999, cores = 4){
  dis <- as.matrix(dis)
  mntd.obs <- mntd(samp, dis, abundance.weighted = abundance.weighted)
  
  ###
  sfInit(parallel = TRUE, cpus = cores, slaveOutfile="logAnalysis.txt")
  
  sfLibrary(picante)
  sfLibrary(base)
  sfLibrary(mosaic)
  sfLibrary(snowfall)
  
  sfExport('samp', 'dis', 'abundance.weighted', 'runs')
  
  mntd.rand <- t(sfSapply(1:runs, 
                          function(i){ 
                            sfCat(paste("Iteration ", i), sep="\n")
                            mntd(samp, taxaShuffle(dis),
                                 abundance.weighted = abundance.weighted)
                          }))
  
  mntd.rand.mean <- sfApply(x = mntd.rand, margin = 2, fun = mean, 
                            na.rm = TRUE)
  
  mntd.rand.sd <- sfApply(x = mntd.rand, margin = 2, fun = sd, 
                          na.rm = TRUE)
  
  mntd.obs.rank <- sfApply(x = rbind(mntd.obs, mntd.rand), margin = 2, 
                           fun = rank)[1, ]
  
  sfStop()
  
  ###
  
  mntd.obs.z <- (mntd.obs - mntd.rand.mean)/mntd.rand.sd
  
  mntd.obs.rank <- ifelse(is.na(mntd.rand.mean), NA, mntd.obs.rank)
  
  data.frame(ntaxa = specnumber(samp), 
             mntd.obs, 
             mntd.rand.mean, 
             mntd.rand.sd, 
             mntd.obs.rank, 
             mntd.obs.z, 
             mntd.obs.p = mntd.obs.rank/(runs + 1), 
             runs = runs, 
             row.names = row.names(samp))
}


#system.time(sf.ses.mpd(phylocom$sample, cophenetic(phylocom$phylo), runs = 1000))
#system.time(ses.mpd(phylocom$sample, cophenetic(phylocom$phylo), null.model= "taxa.labels", runs = 1000))

###############################################################################################
##############
########################
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


brokenStick.selection <- function(dudi.pca.object){
  ## First, choose how many axes are going to be kept
  ev <- dudi.pca.object$eig
  
  # Broken stick model, largely based on Pierre Legendre's class
  
  n <- length(dudi.pca.object$eig)
  bsm <- data.frame(j = seq(1:n), 
                    p = 0)
  bsm$p[1] <- 1/n
  
  for (i in 2:n){
    bsm$p[i] = bsm$p[i-1] + (1/(n + 1 - i))
  }
  
  bsm$p <- 100*bsm$p/n
  
  ## Plot eigenvalues and % of variance for each axis
  
  plot.new()
  
  par(mfrow = c(2,1))
  
  barplot(ev, 
          main="Eigenvalues", 
          col="bisque", las=2)
  
  abline(h=mean(ev), 
         col="red")		# average eigenvalue
  
  legend("topright", 
         "Average eigenvalue", 
         lwd=1, col=2, bty="n")
  
  
  brokenStick.matrix <- t(cbind(100*ev/sum(ev), bsm$p[n:1]))
  brokenStick.kept <- length(which(brokenStick.matrix[1,] > brokenStick.matrix[2,]))
  
  barplot(brokenStick.matrix, 
          beside= TRUE, 
          main = "% variance", col=c("bisque",2), las=2)
  
  legend("topright", 
         c("% eigenvalue", 
           "Broken stick model"), 
         pch= 15, 
         col= c("bisque", 2), 
         bty= "n")
  
  return(list(brokenStick.kept = brokenStick.kept,
              brokenStick.matrix = brokenStick.matrix))
}

##

# This function has been replaced ##

CWM_Std_TW <- function(Trait = tip.rates$lambda.avg, Distrib = Chiroptera.FaurSven.comm){
  
  Trait <- as.matrix(Trait)
  
  Distrib <- as.matrix(Distrib)
  Distrib <- Distrib[ , (colnames(Distrib) %in% row.names(Trait))]
  
  centering_mat <- function(X,w){ X - rep(1,length(w))%*%t(w)%*% X }
  
  standardize_w <- function(X,w){
    ones <- rep(1,length(w))
    Xc <- X - ones %*% t(w)%*% X
    Xc / ones%*%sqrt(t(ones)%*%(Xc*Xc*w)) 
  } 
  
  # check_L()
  rows <- seq_len(nrow(Distrib))
  cols <- seq_len(ncol(Distrib))
  rni <- which(rowSums(Distrib)==0)
  
  repeat {
    if (length(rni)) {Distrib <- Distrib[-rni,,drop = FALSE]; rows <-rows[-rni]}
    ksi <- which(colSums(Distrib)==0)
    if (length(ksi)) {Distrib <- Distrib[,-ksi, drop = FALSE]; cols <- cols[-ksi]}
    rni <- which(rowSums(Distrib)==0)
    if (length(rni)==0 & length(ksi)==0) { 
      break 
    }
  }
  
  # E <-E[rows,,drop = FALSE]
  
  Trait <- Trait[cols, , drop = FALSE]
  # end check_L()
  
  Distrib <- Distrib/sum(Distrib)
  
  # dimensions
  S <- ncol(Distrib) # number of species
  n <- nrow(Distrib) # number of Distribunities
  # p <- ncol(E) # number of environmental predictors
  q <- ncol(Trait) # number of traits
  
  # setting up matrices
  Wn <- rowSums(Distrib)
  Ws <- colSums(Distrib)
  # cor matrices are trait by environment
  CWM <- Distrib%*%Trait/Wn  # weighted means wrt to Trait
  # CWM.cor <- cor(CWM,E)
  
  # SNC <- t(Distrib)%*%E/Ws  # weighted means wrt to E
  # SNC.cor <- cor(Trait,SNC)
  
  CWMstd_w  <- standardize_w(CWM, Wn)
  # Estd_w <- standardize_w(E,Wn)
  # wCWM.cor <- t(t(Estd_w)%*%(CWMstd_w*Wn))
  
  # SNCstd_w <- standardize_w(SNC,Ws)
  Tstd_w <-  standardize_w(Trait, Ws)
  # wSNC.cor <- t(Tstd_w)%*%(SNCstd_w*Ws)
  
  # Fourth corner calculated as W_n weighted covariance between 
  # CWM and standardized Trait (trait)
  
  CWM_std_tw <- Distrib%*%Tstd_w/Wn #CWM wrt to standardized Trait (trait)
  
  return(
    CommWeightMeans <- data.frame(CWM = CWM,
                                  CWM_std_tw = CWM_std_tw,
                                  CWM_std_w = CWMstd_w))
}


####################################
####### Graphical parameters #######
####################################

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
      plot.background = element_rect(fill = "#f5f5f2", color = NA), 
      panel.background = element_rect(fill = "#f5f5f2", color = NA), 
      legend.background = element_rect(fill = "#f5f5f2", color = NA),
      panel.border = element_blank(),
      ...
    )
}

# library(broom)
# library(tidyr)
