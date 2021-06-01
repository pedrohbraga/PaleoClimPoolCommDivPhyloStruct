#### Dataset preparation ####

# Load Afrotropical region grid shapefile
world_grid_2 <- readOGR(dsn = "data/grids", layer = "world_grid_2")

# Load biogeographical realms map
biogeoRealms <- readOGR(dsn = "data/shapefiles", layer = "wwf_terr_realms")

# Load biogeographical realms map
terrBiomes <- readOGR(dsn = "data/shapefiles", layer = "terr_biomes")

# Load ecoregion map
terrEcoregions <- readOGR(dsn = "data/shapefiles", layer = "wwf_terr_ecos") #%>% 
#  disaggregate %>% 
#  geometry


# plot(world_grid_2)

# Change CRS

terrEcoregions <- CRS("+proj=laea +datum=WGS84 +no_defs") %>% 
  spTransform(terrEcoregions, .)

biogeoRealms <- CRS("+proj=laea +datum=WGS84 +no_defs") %>% 
  spTransform(biogeoRealms, .)

world_grid_2 <- CRS("+proj=laea +datum=WGS84 +no_defs") %>% 
  spTransform(world_grid_2, .)


### Extract coordinates
# Get the centroids and then convert them to a SpatialPointsDataFrame
world_grid_2_centroids <- SpatialPointsDataFrame(gCentroid(world_grid_2, byid=TRUE), 
                                                 world_grid_2@data, match.ID=FALSE)

# Overlap coordinates and shapefile and extract required information
overlaping.regions.world_grid_2 <- data.frame(xx=over(world_grid_2_centroids, terrEcoregions))


# Bind Latitude and Longitude and subset ecoregions, biomes and realms for cross-scale analyses
regions.LatLong <- cbind(world_grid_2_centroids@data$ID,
                         world_grid_2_centroids@coords,
                         subset(overlaping.regions.world_grid_2, 
                                select = c(xx.ECO_NAME, xx.BIOME, xx.REALM)))

colnames(regions.LatLong) <- c("ID", "Longitude", "Latitude",
                               "NAME_Ecoregion", "ID_Biome", "ID_Realm")

rownames(regions.LatLong) <- regions.LatLong$ID
levels(regions.LatLong$ID_Realm)[levels(regions.LatLong$ID_Realm)=="NA"] <- "Nearctic"
levels(regions.LatLong$ID_Realm)[levels(regions.LatLong$ID_Realm)=="PA"] <- "Palearctic"
levels(regions.LatLong$ID_Realm)[levels(regions.LatLong$ID_Realm)=="NT"] <- "Neotropical"
levels(regions.LatLong$ID_Realm)[levels(regions.LatLong$ID_Realm)=="AA"] <- "Australasian"
levels(regions.LatLong$ID_Realm)[levels(regions.LatLong$ID_Realm)=="AN"] <- "Antarctic"
levels(regions.LatLong$ID_Realm)[levels(regions.LatLong$ID_Realm)=="AT"] <- "Afrotropical"
levels(regions.LatLong$ID_Realm)[levels(regions.LatLong$ID_Realm)=="IM"] <- "Indomalay"
levels(regions.LatLong$ID_Realm)[levels(regions.LatLong$ID_Realm)=="OC"] <- "Oceanic"

## You could also have done that  by doing this!

Biomes.code <- c("Tropical Subtropical Moist Broadleaf Forests" = 1,
                 "Tropical Subtropical Dry Broadleaf Forests" = 2,
                 "Tropical Subtropical Coniferous Forests" = 3,
                 "Temperate Broadleaf Mixed Forests" = 4,
                 "Temperate Conifer Forests" = 5,
                 "Boreal Forests Taiga" = 6,
                 "Tropical Subtropical Grasslands Savannas Shrublands" = 7,
                 "Temperate Grasslands Savannas Shrublands" = 8,
                 "Flooded Grasslands Savannas" = 9,
                 "Montane Grasslands Shrublands" = 10,
                 "Tundra" = 11,
                 "Mediterranean Forests Woodlands Scrub" = 12,
                 "Deserts Xeric Shrublands" = 13,
                 "Mangroves" = 14)

regions.LatLong$ID_Biome <- names(Biomes.code)[match(regions.LatLong$ID_Biome, Biomes.code)]
regions.LatLong$ID_Biome <- as.factor(regions.LatLong$ID_Biome)

regions.LatLong$ID_Biome_Realm <- paste(regions.LatLong$ID_Biome, "__", regions.LatLong$ID_Realm, sep = "")

### Prepare community and environmental dataset

#############################################
############## Prepare datasets #############
############################################# 

paleo.veg.raw <- read.table("data/matrices/world_paleoveg_2.txt", row.names = 1, header=TRUE)
env.raw <- read.table("data/matrices/world_select_env_2.txt", row.names = 1, header=TRUE)
comm.raw <- read.table("data/matrices/world_bats_2_ALL_PAM.txt", row.names = 1, header=TRUE)

row.names(paleo.veg.raw) <- row.names(regions.LatLong)
row.names(env.raw) <- row.names(regions.LatLong)
row.names(comm.raw) <- row.names(regions.LatLong) 

# Verify if both tables present the right number of rows. 
nrow(env.raw); nrow(comm.raw); nrow(regions.LatLong); nrow(paleo.veg.raw)

head(env.raw)

# Verify cor.plot of environmental variables
plot(env.raw)

# Select chosen variables
new.Env <- data.frame(AnMnTemp_Bio1 = env.raw$bio_1,
                      TempSeasonal_Bio4 = env.raw$bio_4,
                      MnTempWetQuart_Bio8 = env.raw$bio_8,
                      TempSeasonal = env.raw$bio_12,
                      PrecSeasonal_Bio15 = env.raw$bio_15,
                      PrecWetQuart_Bio18 = env.raw$bio_18,
                      DEM = env.raw$dem,
                      Slope = env.raw$h_slope,
                      VegCoverLGM = paleo.veg.raw$vegcover_pct_CLIMAP_LGM_w_2,
                      VegCoverMOD = paleo.veg.raw$vegcover_pct_CLIMAP_MOD_w_2,
                      VegCoverDiff = paleo.veg.raw$diff_vegcover_pct_LGM_MOD,
                      row.names = rownames(env.raw))

# Identify rows that have NA
(rows.NA <- row.names(which(is.na(new.Env), arr.ind=TRUE)))
(rows.NA <- c(rows.NA, row.names(which(is.na(regions.LatLong), arr.ind=TRUE))))
length(unique(rows.NA))
nrow(new.Env) - length(rows.NA) # Difference between total original rows and number of rows with NA

# Remove rows from the Env matrix based on the rows.NA
Env <- new.Env[!(rownames(new.Env) %in% rows.NA), ] 

# Check if everything went well
nrow(Env) + length(unique(rows.NA)); nrow(new.Env)

# Remove the same rows from the community matrix 
new.Comm <- comm.raw[!(rownames(comm.raw) %in% rows.NA), ] 
nrow(new.Comm); nrow(Env) # Verify lengths
Comm <- new.Comm
# Remove species (i.e. columns) that have sum zero
## Comm <- new.Comm[,colSums(new.Comm) >= 2] # Why do this??? No!
ncol(new.Comm); ncol(Comm) # Verify lengths

# Remove rows.NA from the regions and coordinates data.frame

regions.LatLong <- regions.LatLong[!(rownames(regions.LatLong) %in% rows.NA), ] 
# regions.LatLong$NAME_Ecoregion <- paste("__", regions.LatLong$NAME_Ecoregion, sep = '')
# regions.LatLong$ID_Biome <- paste("__", regions.LatLong$ID_Biome, sep = '')

regions.LatLong$ID_Realm <- droplevels(regions.LatLong$ID_Realm)
regions.LatLong$ID_Biome <- droplevels(regions.LatLong$ID_Biome)
nrow(new.Comm); nrow(Env); nrow(regions.LatLong) # Everything should have the same number of rows

#### Rename community data ####

Chiroptera.Comm <- Comm

colnames(Chiroptera.Comm) <- gsub('\\.', "_", colnames(Chiroptera.Comm))
head(Chiroptera.Comm); nrow(Chiroptera.Comm)

#######################################################################################################

#### Load phylogenetic tree ####

Chiroptera.Jones.tree <- read.tree("data/phylogenies/tree.fran.phy")

# Filter species from the phylogenetic tree
#ChiropteraPhylo <- drop.tip(Chiroptera.tree, 
#                               Chiroptera.tree$tip.label[-match(colnames(Chiroptera.Comm), 
#                                                                Chiroptera.tree$tip.label)]) #Prunning the tree
#length(ChiropteraPhylo$tip.label); ncol(Chiroptera.Comm)

Chiroptera.Jones.tree <- match.phylo.comm(Chiroptera.Jones.tree, Chiroptera.Comm)$phy
Chiroptera.Jones.comm <- match.phylo.comm(Chiroptera.Jones.tree, Chiroptera.Comm)$comm
setdiff(colnames(Chiroptera.Comm), Chiroptera.Jones.tree$tip.label)
# ChiropteraPhylo
is.ultrametric(Chiroptera.Jones.tree)

# Force tree to be ultrametric
Chiroptera.Jones.Phylo.ultra <- force.ultrametric(Chiroptera.Jones.tree)

# Check if number of species matches between both trait and phylogenetic relationship datasets
length(Chiroptera.Jones.Phylo.ultra$tip.label); ncol(Chiroptera.Jones.comm)

#######################################################################################################

#### Load phylogenetic tree ####

Mammal.FaurSven.tree.1 <-  read.nexus("data/phylogenies/Fully_resolved_phylogeny_1.nex")
Mammal.FaurSven.tree.2 <-  read.nexus("data/phylogenies/Fully_resolved_phylogeny_2.nex")
Mammal.FaurSven.tree.3 <-  read.nexus("data/phylogenies/Fully_resolved_phylogeny_3.nex")

Mammal.FaurSven.trees <- c(Mammal.FaurSven.tree.1, Mammal.FaurSven.tree.2, Mammal.FaurSven.tree.3)

Mammal.FaurSven.MCC.tree <- phangorn::maxCladeCred(Mammal.FaurSven.trees)

#Chiroptera.tree

#Mammal.FaurSven.tree

# Filter species from the phylogenetic tree
#ChiropteraPhylo <- drop.tip(Chiroptera.tree, 
#                               Chiroptera.tree$tip.label[-match(colnames(Chiroptera.Comm), 
#                                                                Chiroptera.tree$tip.label)]) #Prunning the tree
#length(ChiropteraPhylo$tip.label); ncol(Chiroptera.Comm)

Chiroptera.FaurSven.tree <- match.phylo.comm(Mammal.FaurSven.MCC.tree, Chiroptera.Comm)$phy
Chiroptera.FaurSven.comm <- match.phylo.comm(Mammal.FaurSven.MCC.tree, Chiroptera.Comm)$comm

# ChiropteraPhylo
is.ultrametric(Chiroptera.FaurSven.tree)

# Force tree to be ultrametric
# Chiroptera.FaurSven.Phylo.ultra <- force.ultrametric(Chiroptera.FaurSven.tree)

# Check if number of species matches between both trait and phylogenetic relationship datasets
length(Chiroptera.FaurSven.tree$tip.label); ncol(Chiroptera.FaurSven.comm)

write.tree(Chiroptera.FaurSven.tree, file = "data/phylogenies/Chiroptera.FaurSven.tree.phy")

###################

library(tidyverse)
library(ggtree)

# or add the entire scale to the x axis with theme_tree2()
ggtree(Chiroptera.FaurSven.tree) + theme_tree2()

# Label the tips
ggtree(Chiroptera.FaurSven.tree) + theme_tree2() + geom_tiplab()


setdiff(colnames(Chiroptera.Comm), Mammal.FaurSven.MCC.tree$tip.label)


#######################################################################################################

rm(Comm, 
   comm.raw, 
   env.raw, 
   new.Comm,
   new.Env)
