################################################################################
### Code to prepare the dataset for all analyses performed in this project #####
#                                                                              #
# Author: Pedro Henrique Pereira Braga                                         #
# Last Update: "2022-01-12"                                                    #
#                                                                              #
################################################################################

### Preparing the spatial dataset ----------------------------------------------

# Load the world projection grid shapefile
world_grid_50km <- st_read("data/grids/world_grid_prev_50km.shp")
world_grid_50km$CRI <- NULL

equal_area_projection <- st_crs(world_grid_50km)

# Load the shapefile containing biogeographical realms
biogeoRealms <- st_read("data/shapefiles/wwf_terr_realms.shp")

# Load the shapefile defining tectonic plates
worldPlates <- st_read("data/shapefiles/PB2002_plates.shp")

# Load the terrestrial ecoregions map
terrEcoregions <- st_read("data/shapefiles/wwf_terr_ecos.shp")

# Inspect polygon grid: plot(world_grid_50km)
# Inspect all other polygon grids with plot()

# Inspect the CRS of each simple feature

st_crs(biogeoRealms); st_crs(worldPlates); st_crs(terrEcoregions)

# Transform the coordinate systems of all other shapefiles towards the Berhmann
# projection, which is equivalent to c("+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 
# +datum=WGS84 +ellps=WGS84 +units=m +no_defs")

biogeoRealms <- st_transform(biogeoRealms, crs = equal_area_projection)
worldPlates <- st_transform(worldPlates, crs = equal_area_projection)
terrEcoregions <- st_transform(terrEcoregions, crs = equal_area_projection)

# Inspect the CRS of each simple feature

st_crs(biogeoRealms); st_crs(worldPlates); st_crs(terrEcoregions)

### Manipulating the 50 km by 50 km grid ####

# Attribute an id variable to the world_grid_50km

world_grid_50km$id <- 1:nrow(world_grid_50km)

### Extract coordinates of the centroids for each quadrat from the polygon grid

# Obtain the centroids for the grid cells

world_grid_50km_centroids <- st_centroid(world_grid_50km)

# Create a new object for the world grid, which will contain categorizations for
# each cell

world_grid_50km_cat <- world_grid_50km

# To classify each cell according to the plates where it occurs, intersect
# plates with the grid cell centroids

intersects_plates_world_grid_50km_centroids <- st_intersects(worldPlates,
                                                             world_grid_50km_centroids, 
                                                             sparse = FALSE)

# The above result is boolean, and will require the names of each plate to be
# attributed to each cell. You may do this as follows.

world_grid_50km_cat$PlateName <- apply(intersects_plates_world_grid_50km_centroids,
                                       2,
                                       function(col) { 
                                         worldPlates[which(col), ]$PlateName
                                       })

# We can repeat the same steps above to classify each cell according to the
# realms, biomes, and terrestrial ecoregions where it occurs.

# Begin by intersecting plates with the grid cell centroids

intersects_terrEcoregions_world_grid_50km_centroids <- st_intersects(terrEcoregions,
                                                                     world_grid_50km_centroids, 
                                                                     sparse = FALSE)

# Then, convert the Boolean results from the intersects to region names

world_grid_50km_cat$ECO_NAME <- apply(intersects_terrEcoregions_world_grid_50km_centroids,
                                      2,
                                      function(col) { 
                                        terrEcoregions[which(col), ]$ECO_NAME
                                      })

world_grid_50km_cat$BIOME <- apply(intersects_terrEcoregions_world_grid_50km_centroids,
                                   2,
                                   function(col) { 
                                     terrEcoregions[which(col), ]$BIOME
                                   })

world_grid_50km_cat$REALM <- apply(intersects_terrEcoregions_world_grid_50km_centroids,
                                   2,
                                   function(col) { 
                                     terrEcoregions[which(col), ]$REALM
                                   })

world_grid_50km_centroids_lat_long <- world_grid_50km_centroids %>%
  st_coordinates() %>%
  as.data.frame() %>%
  rename(Longitude = X,
         Latitude = Y)

world_grid_50km_cat_sf <- world_grid_50km_cat %>%
  bind_cols(world_grid_50km_centroids_lat_long) %>%
  rename(ID = id,
         ID_Ecoregion = ECO_NAME,
         ID_Biome = BIOME,
         ID_Realm = REALM,
         ID_PlateName = PlateName) %>%
  unnest(cols = c(ID_Ecoregion, ID_Biome, ID_Realm),
         keep_empty = TRUE)

## Rename and reorder levels to make sense

world_grid_50km_cat_df <- as.data.frame(world_grid_50km_cat_sf) %>%
  mutate(ID_Realm = dplyr::recode_factor(ID_Realm, 
                                         "NA" = "Nearctic",
                                         "NT" = "Neotropical",
                                         "PA" = "Palearctic",
                                         "AT" = "Afrotropical",
                                         "IM" = "Indomalay",
                                         "OC" = "Oceanic",
                                         "AA" = "Australasian",
                                         "AN" = "Antarctic"))

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

world_grid_50km_cat_df$ID_Biome <- names(Biomes.code)[match(world_grid_50km_cat_df$ID_Biome, 
                                                            Biomes.code)]

world_grid_50km_cat_df$ID_Biome <- as.factor(world_grid_50km_cat_df$ID_Biome)

# Create a new variable that nests biomes within their respective realms
world_grid_50km_cat_df$ID_Biome_Realm <- paste(world_grid_50km_cat_df$ID_Biome, 
                                               "__",
                                               world_grid_50km_cat_df$ID_Realm, 
                                               sep = "")
# Define NA values

world_grid_50km_cat_df <- world_grid_50km_cat_df %>%
  naniar::replace_with_na(replace = list(ID_Biome_Realm = c("NA__NA")))


#### Preparing the bat community dataset ---------------------------------------

# Load community data ###

Chiroptera.Comm <- read.csv("data/matrices/world_chiroptera_PropCover_50KM.csv",
                            header = TRUE)

nrow(Chiroptera.Comm); nrow(world_grid_50km_cat_df) # Verify lengths

# Replace periods in species names by underscores
colnames(Chiroptera.Comm) <- gsub('\\.', 
                                  "_", 
                                  colnames(Chiroptera.Comm))

head(Chiroptera.Comm); nrow(Chiroptera.Comm)

# Remove species that do not occur in any community
# To identify them: colnames(Chiroptera.Comm[,colSums(Chiroptera.Comm) < 1])

# Remove species (i.e. columns) that have sum zero
Chiroptera.Comm <- Chiroptera.Comm[, colSums(Chiroptera.Comm) >= 1]

# Drop levels, but not working so far

world_grid_50km_cat_df$ID_Realm <- droplevels(world_grid_50km_cat_df$ID_Realm)
world_grid_50km_cat_df$ID_Biome <- droplevels(world_grid_50km_cat_df$ID_Biome)

world_grid_50km_cat_df <- world_grid_50km_cat_df %>%
  select(-QuadratID, -geometry)

### Preparing the phylogenetic dataset -----------------------------------------

#### Load phylogenetic trees from Faurby and Svenning, 2015 () #

Mammal.FaurSven.tree.1 <-  read.nexus("data/phylogenies/Fully_resolved_phylogeny_1.nex")
Mammal.FaurSven.tree.2 <-  read.nexus("data/phylogenies/Fully_resolved_phylogeny_2.nex")
Mammal.FaurSven.tree.3 <-  read.nexus("data/phylogenies/Fully_resolved_phylogeny_3.nex")

Mammal.FaurSven.trees <- c(Mammal.FaurSven.tree.1, 
                           Mammal.FaurSven.tree.2, 
                           Mammal.FaurSven.tree.3)

# Obtain the maximum credibility tree

Mammal.FaurSven.MCC.tree <- phangorn::maxCladeCred(Mammal.FaurSven.trees)

# Calculate the log clade credibility scores, which will be used to sample other
# trees to assess whether the results from the maximum credibility clade tree
# differs across different trees

Mammal.FaurSven.MCC.tree.scores <- phangorn::maxCladeCred(Mammal.FaurSven.trees, 
                                                          tree = FALSE)

# Represent log credibility scores in a histogram 

histogram <- hist(Mammal.FaurSven.MCC.tree.scores,
                  breaks = 25,
                  col = "gray",
                  main = "",
                  xlab = "log(Clade Credibility for Faurby and Svenning's (2015) tree)")

arrows(x0 = max(Mammal.FaurSven.MCC.tree.scores), 
       y0 = 0.5*max(histogram$counts),
       x1 = max(Mammal.FaurSven.MCC.tree.scores),
       y1 = 0,
       col = "red",
       lwd = 2,
       length = 0.15, 
       angle = 20)

text(x = max(Mammal.FaurSven.MCC.tree.scores),
     y = 0.5*max(histogram$counts)+ 1,
     "Clade credibility of the tree we used",
     adj = c(0, 0.3),
     srt = 90)

arrows(x0 = median(Mammal.FaurSven.MCC.tree.scores), 
       y0 = 0.5*median(histogram$counts),
       x1 = median(Mammal.FaurSven.MCC.tree.scores),
       y1 = 0,
       col = "blue",
       lwd = 2,
       length = 0.15,
       angle = 20)

Mammal.FaurSven.MCC.tree.scores

Mammal.FaurSven.MCC.tree.scores[Mammal.FaurSven.MCC.tree.scores < quantile(Mammal.FaurSven.MCC.tree.scores, 0.05)]

# Filter species from the phylogenetic tree

# ChiropteraPhylo <- drop.tip(Chiroptera.tree, 
#                               Chiroptera.tree$tip.label[-match(colnames(Chiroptera.Comm), 
#                                                                Chiroptera.tree$tip.label)]) #Prunning the tree
# length(ChiropteraPhylo$tip.label); ncol(Chiroptera.Comm)

Chiroptera.FaurSven.tree <- match.phylo.comm(Mammal.FaurSven.MCC.tree, 
                                             Chiroptera.Comm)$phy
Chiroptera.FaurSven.comm <- match.phylo.comm(Mammal.FaurSven.MCC.tree, 
                                             Chiroptera.Comm)$comm

# Save results

write.tree(Mammal.FaurbySven.tree, file = "data/phylogenies/Mammal.FaurbySven.tree.phy")
write.tree(Chiroptera.FaurSven.tree, file = "data/phylogenies/Chiroptera.FaurSven.tree.phy")

# For future analyses, load them from here
# Mammal.FaurSven.MCC.tree <- read.tree("data/phylogenies/Mammal.FaurbySven.tree.phy")
# Chiroptera.FaurSven.tree <- read.tree("data/phylogenies/Chiroptera.FaurSven.tree.phy")

# ChiropteraPhylo
is.ultrametric(Chiroptera.FaurSven.tree)

# Force tree to be ultrametric
Chiroptera.FaurSven.Phylo.ultra <- force.ultrametric(Chiroptera.FaurSven.tree)

Chiroptera.FaurSven.Phylo.ultra

# Check if number of species matches between both trait and phylogenetic relationship datasets
length(Chiroptera.FaurSven.tree$tip.label); ncol(Chiroptera.FaurSven.comm)

######################################################################################################

# library(tidyverse)
# library(ggtree)

# or add the entire scale to the x axis with theme_tree2()
# ggtree(Chiroptera.FaurSven.tree) + theme_tree2()

# Label the tips
# ggtree(Chiroptera.FaurSven.tree) + theme_tree2() + geom_tiplab()

######################################################################################################

rm(intersects_plates_world_grid_50km_centroids,
   intersects_terrEcoregions_world_grid_50km_centroids,
   Mammal.FaurSven.tree.1,
   Mammal.FaurSven.tree.2,
   Mammal.FaurSven.tree.3)

#######################################################################################################

