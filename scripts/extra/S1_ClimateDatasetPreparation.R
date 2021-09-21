# Decomissioned as velox does not exist anymore

list.files("data/rasters")
list.files("data/rasters")[1]

rasterFilenames <- list.files("data/rasters", 
           full.names = TRUE, pattern = ".zip")

for (i in 1:length(rasterFilenames)) {
  unzip(rasterFilenames[i],
        exdir = gsub('.zip', 
                     '', 
                     rasterFilenames[i]))
}

list.dirs("data/rasters")

##

require(raster)
require(velox)

## Load rasters as a RasterStack

for (j in 8:length(rasterFilenames)) {
  
# Generate a list of input rasters ("grids")
# Pattern = "*.asc$" - filters for main raster files only and skips any associated files (e.g. world files)
(bioclimRasters <- list.files(list.dirs("data/rasters")[j], 
                              pattern = "*.tif$"))

#create a raster stack from the input raster files 
s <- raster::stack(paste0(list.dirs("data/rasters")[j], "/", 
                          bioclimRasters))

# If you do not have the bioclimatic data yet, you can download it using this:
# s <- raster::getData('worldclim', 
#             var = 'bio', 
#             res = 2.5)

# In case the dataset you have has no CRS defined:
# proj4string(s) <- c("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

# Project your raster to Behrmann CRS, to match the polygon grid

s <- raster::projectRaster(s, 
                           crs = ('+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs'), # Behrmann CRS
                           bylayer = TRUE,  
                           method = 'bilinear')

# You may save the reprojected rasters, as this takes a while: 

# writeRaster(s,
#            filename = names(s), 
#            bylayer = TRUE, 
#            format = "EHdr") # .bil format

#### Fast raster value extraction to polygons ####

## Please note that all raster data will be saved to memory,
## and therefore, a lot of free memory space is required.

environmentalData <- matrix(NA, 
                            nrow(world_grid_50km), 
                            dim(s)[3])

colnames(environmentalData) <- names(s)

for(i in names(s)){
  
  ## Cast as VeloxRaster
  vx.i <- velox(s[[i]])
  
  ## Extract values and calculate mean
  environmentalData[, i] <- vx.i$extract(sp = world_grid_50km, 
                                         fun = function(x) mean(x, 
                                                                na.rm = TRUE),
                                         small = TRUE)
}

# Inspect
head(environmentalData)

#### Save extraction results ####
## If you use a large dataset, consider using ff::fwrite()

write.csv(environmentalData, paste0("data/matrices/world_", list.dirs("data/rasters", full.names = FALSE)[j], "_50km", ".csv"))
}


