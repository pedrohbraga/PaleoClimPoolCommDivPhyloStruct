################################################################################
### Code to prepare the spatial dataset for  this project                  #####
#                                                                              #
# Author: Pedro Henrique Pereira Braga                                         #
# Last Update: "2022-01-12"                                                    #
#                                                                              #
################################################################################

### Preparing the spatial dataset ----------------------------------------------

# Load the world projection grid shapefile
world_grid_50km <- st_read("data/grids/world_grid_50km.shp")

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
# attributed to each cell. You may do this as follows. It may take a while.

world_grid_50km_cat$PlateName <- apply(intersects_plates_world_grid_50km_centroids,
                                       2,
                                       function(col) { 
                                         worldPlates[which(col), ]$PlateName
                                       })

# We can repeat the same steps above to classify each cell according to the
# realms, biomes, and terrestrial ecoregions where it occurs.

# Begin by intersecting plates with the grid cell centroids:

intersects_terrEcoregions_world_grid_50km_centroids <- st_intersects(terrEcoregions,
                                                                     world_grid_50km_centroids, 
                                                                     sparse = FALSE)

# Then, convert the Boolean results from the intersects to region names:

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

# Rename X to Longitude and Y to Latitude

world_grid_50km_centroids_lat_long <- world_grid_50km_centroids %>%
  st_coordinates() %>%
  as.data.frame() %>%
  rename(Longitude = X,
         Latitude = Y)

# Rename other columns and unnest data

world_grid_50km_cat_sf <- world_grid_50km_cat %>%
  bind_cols(world_grid_50km_centroids_lat_long) %>%
  rename(ID = id,
         ID_Ecoregion = ECO_NAME,
         ID_Biome = BIOME,
         ID_Realm = REALM,
         ID_PlateName = PlateName) %>%
  unnest(cols = c(ID_Ecoregion, ID_Biome, ID_Realm),
         keep_empty = TRUE)

# Compute the area of quadrats. Note they should be the same because we are using an
# equal area projection.

world_grid_50km_cat_sf$Quadrat_AREA <- world_grid_50km_cat_sf %>%
  st_area()

## Rename and reorder levels for realms

world_grid_50km_cat_df <- as.data.frame(world_grid_50km_cat_sf) %>%
  mutate(ID_Realm = dplyr::recode_factor(ID_Realm, 
                                         "NA" = "Nearctic",
                                         "NT" = "Neotropical",
                                         "PA" = "Palearctic",
                                         "AT" = "Afrotropical",
                                         "IM" = "Indomalay",
                                         "OC" = "Oceanic",
                                         "AA" = "Australasian",
                                         "AN" = "Antarctic")
  )

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
# Replace regions that were NA__NA (NA biomes and NA realms) with NA values

world_grid_50km_cat_df <- world_grid_50km_cat_df %>%
  naniar::replace_with_na(replace = list(ID_Biome_Realm = c("NA__NA"))) 

# The code below is not used. It just allowed me to quickly explore the area of each
# region used in this study.

# area_ID_Hemisphere <- world_grid_50km_cat_df %>%
#   mutate(ID_Region = ifelse(ID_Realm %in% c("Nearctic", "Neotropical"), 
#                             "New World",
#                             ifelse(ID_Realm %in% c("Palearctic",
#                                                    "Afrotropical",
#                                                    "Indomalay",
#                                                    "Oceanic",
#                                                    "Australasian"), "Old World",
#                                    NA))) %>%
#   group_by(ID_Region) %>%
#   dplyr::summarise(sum_Quadrat_AREA = 
#                      sum(Quadrat_AREA, na.rm = T)) %>%
#   mutate(ID_extent = rep("ID_Hemisphere"))
# 
# area_ID_Realm <- world_grid_50km_cat_df %>%
#   rename(ID_Region = ID_Realm) %>%
#   group_by(ID_Region) %>%
#   dplyr::summarise(sum_Quadrat_AREA = 
#                      sum(Quadrat_AREA, na.rm = T)) %>%
#   mutate(ID_extent = rep("ID_Realm"))
# 
# area_ID_PlateName <- world_grid_50km_cat_df %>%
#   rename(ID_Region = ID_PlateName) %>%
#   group_by(ID_Region) %>%
#   dplyr::summarise(sum_Quadrat_AREA = 
#                      sum(Quadrat_AREA, na.rm = T)) %>%
#   mutate(ID_extent = rep("ID_PlateName"))
# 
# area_ID_Biome_Realm <- world_grid_50km_cat_df %>%
#   rename(ID_Region = ID_Biome_Realm) %>%
#   group_by(ID_Region) %>%
#   dplyr::summarise(sum_Quadrat_AREA = 
#                      sum(Quadrat_AREA, na.rm = T)) %>%
#   
#   mutate(ID_extent = rep("ID_Biome_Realm"))
# 
# area_ID_Ecoregion <- world_grid_50km_cat_df %>%
#   rename(ID_Region = ID_Ecoregion) %>%
#   group_by(ID_Region) %>%
#   dplyr::summarise(sum_Quadrat_AREA = 
#                      sum(Quadrat_AREA, na.rm = T)) %>%
#   mutate(ID_extent = rep("ID_Ecoregion"))
# 
# area_ID_extent <- bind_rows(area_ID_Hemisphere,
#                             area_ID_Realm,
#                             area_ID_PlateName,
#                             area_ID_Biome_Realm,
#                             area_ID_Ecoregion)
# 
# area_ID_extent$sum_Quadrat_AREA <- units::set_units(area_ID_extent$sum_Quadrat_AREA, km^2)
# 
# area_ID_extent %>%
#   group_by(ID_extent) %>%
#   summarise(mean_area = mean(sum_Quadrat_AREA, na.rm =T)
#   )
# 
# ggplot(data = area_ID_extent,
#        aes(x = ID_extent,
#            y =  sum_Quadrat_AREA)) +
#   geom_point() 
# 
# bind_rows(area_ID_Hemisphere,
#           area_ID_Realm,
#           area_ID_PlateName,
#           area_ID_Biome_Realm,
#           area_ID_Ecoregion) 

# Checking levels

str(
  world_grid_50km_cat_df %>%
    droplevels()
)

######################################################################################

rm(
  intersects_plates_world_grid_50km_centroids,
  intersects_terrEcoregions_world_grid_50km_centroids,
)
