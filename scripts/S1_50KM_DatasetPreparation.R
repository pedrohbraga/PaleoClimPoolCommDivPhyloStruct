######################################################################################
######## Code to prepare the dataset for all analyses performed in this project ######
#                                                                                    #
# Author: Pedro Henrique Pereira Braga                                               #
# Last Update: "2019-01-30 15:29:30 EST"                                             #
#                                                                                    #
# Alert: this script needs to be updated to integrate sf and exclude bioclim parts   #
######################################################################################

#############################
#### Dataset preparation ####
#############################

# Load world projection grid shapefile
world_grid_50km <- readOGR(dsn = "data/grids", 
                           layer = "world_grid_prev_50km")

# Load biogeographical realms map
biogeoRealms <- readOGR(dsn = "data/shapefiles", 
                        layer = "wwf_terr_realms")

# Load biogeographical realms map
worldPlates <- readOGR(dsn = "data/shapefiles", 
                       layer = "PB2002_plates")

# Load terrestrial ecoregions map
terrEcoregions <- readOGR(dsn = "data/shapefiles", 
                          layer = "wwf_terr_ecos") #%>% 
#  disaggregate %>% 
#  geometry

# Inspect polygon grid: plot(world_grid_50km)
# Inspect all other polygon grids with plot()

#### Define and match CRS for the spatial grid polygons ####
proj4string(world_grid_50km) <- c("+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 
                                  +datum=WGS84 +ellps=WGS84 +units=m +no_defs")

# Attribute an "id" variable
world_grid_50km@data$id = 1:length(world_grid_50km)

# Define and match CRS for the other shapefiles
terrEcoregions <- CRS("+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 
                      +datum=WGS84 +ellps=WGS84 +units=m +no_defs") %>% 
  spTransform(terrEcoregions, .)

worldPlates <- CRS("+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 
                   +datum=WGS84 +ellps=WGS84 +units=m +no_defs") %>% 
  spTransform(worldPlates, .)

biogeoRealms <- CRS("+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 
                    +datum=WGS84 +ellps=WGS84 +units=m +no_defs") %>% 
  spTransform(biogeoRealms, .)


### Extract coordinates of the centroids of each quadrat from the polygon grid ###

# Get the centroids and then convert them to a SpatialPointsDataFrame
world_grid_50km_centroids <- SpatialPointsDataFrame(gCentroid(world_grid_50km, 
                                                              byid = TRUE), 
                                                    world_grid_50km@data, 
                                                    match.ID = FALSE)

# Overlap coordinates and shapefile and extract required information
overlaping.regions.world_grid_50km <- data.frame(xx = over(world_grid_50km_centroids, 
                                                           terrEcoregions))

overlaping.plates.world_grid_50km <- data.frame(xx = over(world_grid_50km_centroids, 
                                                          worldPlates))

# Bind latitude and longitude and subset ecoregions, biomes and realms for cross-scale analyses
regions.LatLong <- cbind(world_grid_50km_centroids@data$id,
                         world_grid_50km_centroids@coords,
                         subset(overlaping.regions.world_grid_50km, 
                                select = c(xx.ECO_NAME, xx.BIOME, xx.REALM)),
                         subset(overlaping.plates.world_grid_50km, 
                                select = c(xx.PlateName, xx.Code)))

# Rename column names
colnames(regions.LatLong) <- c("ID", "Longitude", "Latitude",
                               "NAME_Ecoregion", 
                               "ID_Biome", 
                               "ID_Realm",
                               "ID_PlateName", 
                               "ID_PlateCode")

# Assign ID to rownames
rownames(regions.LatLong) <- regions.LatLong$ID

## Rename and reorder levels to make sense
levels(regions.LatLong$ID_Realm)[levels(regions.LatLong$ID_Realm)=="NA"] <- "Nearctic"
levels(regions.LatLong$ID_Realm)[levels(regions.LatLong$ID_Realm)=="PA"] <- "Palearctic"
levels(regions.LatLong$ID_Realm)[levels(regions.LatLong$ID_Realm)=="NT"] <- "Neotropical"
levels(regions.LatLong$ID_Realm)[levels(regions.LatLong$ID_Realm)=="AA"] <- "Australasian"
levels(regions.LatLong$ID_Realm)[levels(regions.LatLong$ID_Realm)=="AN"] <- "Antarctic"
levels(regions.LatLong$ID_Realm)[levels(regions.LatLong$ID_Realm)=="AT"] <- "Afrotropical"
levels(regions.LatLong$ID_Realm)[levels(regions.LatLong$ID_Realm)=="IM"] <- "Indomalay"
levels(regions.LatLong$ID_Realm)[levels(regions.LatLong$ID_Realm)=="OC"] <- "Oceanic"

## Do the same for the biomes, in a slightly different approach
Biomes.code <- c("Tropical_Subtropical_Moist_Broadleaf_Forests" = 1,
                 "Tropical_Subtropical_Dry_Broadleaf_Forests" = 2,
                 "Tropical_Subtropical_Coniferous_Forests" = 3,
                 "Temperate_Broadleaf_Mixed_Forests" = 4,
                 "Temperate_Conifer_Forests" = 5,
                 "Boreal_Forests_Taiga" = 6,
                 "Tropical_Subtropical_Grasslands_Savannas_Shrublands" = 7,
                 "Temperate_Grasslands_Savannas_Shrublands" = 8,
                 "Flooded_Grasslands_Savannas" = 9,
                 "Montane_Grasslands_Shrublands" = 10,
                 "Tundra" = 11,
                 "Mediterranean_Forests_Woodlands_Scrub" = 12,
                 "Deserts_Xeric_Shrublands" = 13,
                 "Mangroves" = 14)

regions.LatLong$ID_Biome <- names(Biomes.code)[match(regions.LatLong$ID_Biome, Biomes.code)]
regions.LatLong$ID_Biome <- as.factor(regions.LatLong$ID_Biome)

# Create a new variable that nests biomes within their respective realms
regions.LatLong$ID_Biome_Realm <- paste(regions.LatLong$ID_Biome, 
                                        "__", 
                                        regions.LatLong$ID_Realm, 
                                        sep = "")

#####################################################
### Prepare community and environmental datasets ####
##################################################### 

# Load environmental data

current.env.raw <- read.csv("data/matrices/world_Env_50KM.csv", row.names = 1, header=TRUE)
lgm.env.raw <- read.csv("data/matrices/world_cclgmbi_env_50KM.csv", row.names = 1, header=TRUE)
mid.holo.env.raw <- read.csv("data/matrices/world_ccmidbi_env_50KM.csv", row.names = 1, header=TRUE)

# Verify if both tables present the right number of rows. 
nrow(current.env.raw); nrow(lgm.env.raw); nrow(mid.holo.env.raw); nrow(regions.LatLong) # nrow(comm.raw); 

head(current.env.raw)

# Verify cor.plot of environmental variables
#plot(current.env.raw)

#### Applying transformations ####

# Select chosen variables
new.current.Env <- data.frame(bio_1 = current.env.raw$bio1,
                              bio_2 = current.env.raw$bio2,
                              bio_3 = current.env.raw$bio3,
                              bio_4 = current.env.raw$bio4,
                              bio_5 = current.env.raw$bio5,
                              bio_6 = current.env.raw$bio6,
                              bio_7 = current.env.raw$bio7,
                              bio_8 = current.env.raw$bio8,
                              bio_9 = current.env.raw$bio9,
                              bio_10 = current.env.raw$bio10,
                              bio_11 = current.env.raw$bio11,
                              bio_12 = base::log1p(current.env.raw$bio12),
                              bio_13 = base::log1p(current.env.raw$bio13),
                              bio_14 = base::log1p(current.env.raw$bio14),
                              bio_15 = current.env.raw$bio15,
                              bio_16 = base::log1p(current.env.raw$bio16),
                              bio_17 = base::log1p(current.env.raw$bio17),
                              bio_18 = base::log1p(current.env.raw$bio18),
                              bio_19 = base::log1p(current.env.raw$bio19),
                              row.names = row.names(regions.LatLong))

# Select chosen variables
new.lgm.Env <- data.frame(bio_1 = lgm.env.raw$cclgmbi1,
                          bio_2 = lgm.env.raw$cclgmbi2,
                          bio_3 = lgm.env.raw$cclgmbi3,
                          bio_4 = lgm.env.raw$cclgmbi4,
                          bio_5 = lgm.env.raw$cclgmbi5,
                          bio_6 = lgm.env.raw$cclgmbi6,
                          bio_7 = lgm.env.raw$cclgmbi7,
                          bio_8 = lgm.env.raw$cclgmbi8,
                          bio_9 = lgm.env.raw$cclgmbi9,
                          bio_10 = lgm.env.raw$cclgmbi10,
                          bio_11 = lgm.env.raw$cclgmbi11,
                          bio_12 = base::log1p(lgm.env.raw$cclgmbi12),
                          bio_13 = base::log1p(lgm.env.raw$cclgmbi13),
                          bio_14 = base::log1p(lgm.env.raw$cclgmbi14),
                          bio_15 = lgm.env.raw$cclgmbi15,
                          bio_16 = base::log1p(lgm.env.raw$cclgmbi16),
                          bio_17 = base::log1p(lgm.env.raw$cclgmbi17),
                          bio_18 = base::log1p(lgm.env.raw$cclgmbi18),
                          bio_19 = base::log1p(lgm.env.raw$cclgmbi19),
                          row.names = row.names(regions.LatLong))

# Select chosen variables
new.mid.holo.Env <- data.frame(bio_1 = mid.holo.env.raw$ccmidbi1,
                               bio_2 = mid.holo.env.raw$ccmidbi2,
                               bio_3 = mid.holo.env.raw$ccmidbi3,
                               bio_4 = mid.holo.env.raw$ccmidbi4,
                               bio_5 = mid.holo.env.raw$ccmidbi5,
                               bio_6 = mid.holo.env.raw$ccmidbi6,
                               bio_7 = mid.holo.env.raw$ccmidbi7,
                               bio_8 = mid.holo.env.raw$ccmidbi8,
                               bio_9 = mid.holo.env.raw$ccmidbi9,
                               bio_10 = mid.holo.env.raw$ccmidbi10,
                               bio_11 = mid.holo.env.raw$ccmidbi11,
                               bio_12 = base::log1p(mid.holo.env.raw$ccmidbi12),
                               bio_13 = base::log1p(mid.holo.env.raw$ccmidbi13),
                               bio_14 = base::log1p(mid.holo.env.raw$ccmidbi14),
                               bio_15 = mid.holo.env.raw$ccmidbi15,
                               bio_16 = base::log1p(mid.holo.env.raw$ccmidbi16),
                               bio_17 = base::log1p(mid.holo.env.raw$ccmidbi17),
                               bio_18 = base::log1p(mid.holo.env.raw$ccmidbi18),
                               bio_19 = base::log1p(mid.holo.env.raw$ccmidbi19),
                               row.names = row.names(regions.LatLong))

new.current.Env <- as.data.frame(zoo::na.spline(new.current.Env))

new.mid.holo.Env <- as.data.frame(zoo::na.spline(new.mid.holo.Env))

new.lgm.Env <- as.data.frame(zoo::na.spline(new.lgm.Env))

# Identify rows that have NA in all environmental datasets
(rows.NA <- c(row.names(which(is.na(new.current.Env), 
                              arr.ind=TRUE)), 
              row.names(which(is.na(new.mid.holo.Env), 
                              arr.ind=TRUE)),
              row.names(which(is.na(new.lgm.Env), 
                              arr.ind=TRUE))))

# (rows.NA <- c(rows.NA, row.names(which(is.na(regions.LatLong), arr.ind=TRUE))))

length(unique(rows.NA))

nrow(regions.LatLong) - length(rows.NA) # Difference between total original rows and number of rows with NA

# Remove rows from the current.Env matrix based on the rows.NA
current.Env <- new.current.Env[!(rownames(new.current.Env) %in% rows.NA), ] 
lgm.Env <- new.lgm.Env[!(rownames(new.lgm.Env) %in% rows.NA), ] 
mid.holo.Env <- new.mid.holo.Env[!(rownames(new.mid.holo.Env) %in% rows.NA), ] 

# Check if everything went well: unequal numbers flag problems!
nrow(current.Env) + length(unique(rows.NA)); nrow(new.current.Env)
nrow(lgm.Env) + length(unique(rows.NA)); nrow(new.lgm.Env)
nrow(mid.holo.Env) + length(unique(rows.NA)); nrow(new.mid.holo.Env)

#### Preparing the community data ####

# Load community data
comm.raw <- read.csv("data/matrices/world_chiroptera_PropCover_50KM.csv",header=TRUE)

# Remove the same rows from the community matrix 
new.Comm <- comm.raw[!(rownames(comm.raw) %in% rows.NA), ] 
nrow(new.Comm); nrow(regions.LatLong) # Verify lengths

# Remove species (i.e. columns) that have sum zero
Comm <- new.Comm[,colSums(new.Comm) >= 1] # Why do this??? No!
ncol(new.Comm); ncol(Comm) # Verify lengths

#### Rename community data ###

Chiroptera.Comm <- Comm

# Replace periods in species names by underscores
colnames(Chiroptera.Comm) <- gsub('\\.', "_", colnames(Chiroptera.Comm))

head(Chiroptera.Comm); nrow(Chiroptera.Comm)

### A few more levels of cleaning ####

# Remove rows.NA from the regions and coordinates data.frame

regions.LatLong <- regions.LatLong[!(rownames(regions.LatLong) %in% rows.NA), ] 
# regions.LatLong$NAME_Ecoregion <- paste("__", regions.LatLong$NAME_Ecoregion, sep = '')
# regions.LatLong$ID_Biome <- paste("__", regions.LatLong$ID_Biome, sep = '')

regions.LatLong$ID_Realm <- droplevels(regions.LatLong$ID_Realm)
regions.LatLong$ID_Biome <- droplevels(regions.LatLong$ID_Biome)

nrow(new.Comm); nrow(current.Env); nrow(regions.LatLong) # Everything should have the same number of rows

#########################################
### Prepare the phylogenetic dataset ####
#########################################

#### Load phylogenetic trees from Faurby and Svenning ####

Mammal.FaurSven.tree.1 <-  read.nexus("data/phylogenies/Fully_resolved_phylogeny_1.nex")
Mammal.FaurSven.tree.2 <-  read.nexus("data/phylogenies/Fully_resolved_phylogeny_2.nex")
Mammal.FaurSven.tree.3 <-  read.nexus("data/phylogenies/Fully_resolved_phylogeny_3.nex")

Mammal.FaurSven.trees <- c(Mammal.FaurSven.tree.1, 
                           Mammal.FaurSven.tree.2, 
                           Mammal.FaurSven.tree.3)

# Calculate a maximum credibility tree
Mammal.FaurSven.MCC.tree <- phangorn::maxCladeCred(Mammal.FaurSven.trees)

# Filter species from the phylogenetic tree

#ChiropteraPhylo <- drop.tip(Chiroptera.tree, 
#                               Chiroptera.tree$tip.label[-match(colnames(Chiroptera.Comm), 
#                                                                Chiroptera.tree$tip.label)]) #Prunning the tree
#length(ChiropteraPhylo$tip.label); ncol(Chiroptera.Comm)

Chiroptera.FaurSven.tree <- match.phylo.comm(Mammal.FaurSven.MCC.tree, Chiroptera.Comm)$phy
Chiroptera.FaurSven.comm <- match.phylo.comm(Mammal.FaurSven.MCC.tree, Chiroptera.Comm)$comm

# Save results

write.tree(Mammal.FaurbySven.tree, file = "data/phylogenies/Mammal.FaurbySven.tree.phy")
write.tree(Chiroptera.FaurSven.tree, file = "data/phylogenies/Chiroptera.FaurSven.tree.phy")

# For future analyses, load them from here
# Mammal.FaurSven.MCC.tree <- read.tree("data/phylogenies/Mammal.FaurbySven.tree.phy")
# Chiroptera.FaurSven.tree <- read.tree("data/phylogenies/Chiroptera.FaurSven.tree.phy")

# ChiropteraPhylo
is.ultrametric(Chiroptera.FaurSven.tree)

# Force tree to be ultrametric
# Chiroptera.FaurSven.Phylo.ultra <- force.ultrametric(Chiroptera.FaurSven.tree)

# Check if number of species matches between both trait and phylogenetic relationship datasets
length(Chiroptera.FaurSven.tree$tip.label); ncol(Chiroptera.FaurSven.comm)

###################

library(tidyverse)
library(ggtree)

# or add the entire scale to the x axis with theme_tree2()
# ggtree(Chiroptera.FaurSven.tree) + theme_tree2()

# Label the tips
# ggtree(Chiroptera.FaurSven.tree) + theme_tree2() + geom_tiplab()

#######################################################################################################

rm(Comm, 
   comm.raw, 
   current.env.raw, mid.holo.env.raw, lgm.env.raw, 
   new.Comm,
   new.current.Env, new.mid.holo.Env, new.lgm.Env,
   Mammal.FaurSven.tree.1,
   Mammal.FaurSven.tree.2,
   Mammal.FaurSven.tree.3,
   Mammal.FaurSven.trees,
   Mammal.FaurSven.MCC.tree)

#######################################################################################################
