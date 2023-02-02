############################################################
#### 
############################################################

library(fasterize)

#### Preparing data ####

## Manipulating the polygon grid

# Load polygon grid data
world_grid_50km <- st_read("/home/shared/pedro.braga/chapter-BiogHistPhyloRelatedSpatScales/data/shapefiles/world_grid_prev_50km.shp") %>%
  rename("Quadrat_ID" = "QuadratID")

# Determine the common CRS for our data 
equal_area_projection <- c("+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs")

# Set the polygon's grid CRS
st_crs(world_grid_50km) <- equal_area_projection

# Create a template base raster that will have the extent of the polygon grid
# and our wanted resolution

templateBase <- raster(world_grid_50km, 
                       res = 50000, 
                       ext = extent(world_grid_50km)) # 5000m resolution

st_area(world_grid_50km); area(templateBase)@data@values

# extent(templateBase) <- extent(world_grid_50km)

# You may rasterize the grid to the template base
templateBase.fasterized <- fasterize(world_grid_50km["Quadrat_ID"], 
                                     templateBase)


# Obtain file names of all rasters that are going to be processed

(rasterFilenames <- list.files(path = export_dir, 
#                               pattern='.tif$',
pattern='bio14',
                               all.files = TRUE, 
                               full.names = TRUE,
                               include.dirs = TRUE,
                               recursive = FALSE)
)

rasterFilenames <- c(sort(rasterFilenames[1:40], 
                          decreasing = TRUE),
                     rasterFilenames[41:45])




# You can subset one of the rasters as dummy to see how the following function is working
rasterFilename.i <- rasterFilenames[1]

# resampled.rasterStack.CHELSA_TraCE21k.cea.stack <- stack(templateBase.fasterized)
# 

# Setting a function that can work in parallel to process, project and extract the data
extractRasterSF <- function(rasterFilename.i){
  
  psnice(pid = Sys.getpid(), value = 19)
  
  rasterOptions(maxmemory = Inf,
                memfrac = 0.9,
                chunksize = 1e+20)
  
  # Progress follow-up 
  cat('Processing raster', basename(rasterFilename.i), 'of', length(rasterFilenames),'\n')
  
  # Note processing time to determine how long each iteration takes
  counter.proc.time <- proc.time()[3]
  
  # Load TIF as a rasterBrick
  rasterStack.CHELSA_TraCE21k.i <- brick(rasterFilename.i)
  
  # Progress follow-up
  cat('Projecting rasters to from', basename(rasterFilename.i), "to", equal_area_projection, '\n')
  
  
  #   # Project rasterBrick to the desired CRS and export it
  rasterStack.CHELSA_TraCE21k.cea.i <- projectRaster(rasterStack.CHELSA_TraCE21k.i,
                                                     to = templateBase,
                                                     crs = equal_area_projection,
                                                     filename = paste0(
                                                       export_dir,"/",
                                                       sub(pattern = "(.*)\\..*$",
                                                           replacement = "\\1",
                                                           basename(rasterFilename.i)
                                                       ),
                                                       "_cea30_50km", ".tif"),
                                                     format = "GTiff", 
                                                     overwrite = TRUE)
  
  # Perform the exact extraction of the raster data to the polygon grid
  extract.rasterStack.CHELSA_TraCE21k.i <- exactextractr::exact_extract(
    rasterStack.CHELSA_TraCE21k.cea.i, 
    world_grid_50km,
    include_xy = FALSE, # include coordinates
    include_cell = FALSE, # include cell IDs
    fun = "mean",
    force_df = TRUE,
    max_cells_in_memory = 3e+08# ,
    # include_cols = rasterStack.CHELSA_TraCE21k.cea.i@data@names
    # full_colnames = T
  )
  
  colnames(extract.rasterStack.CHELSA_TraCE21k.i) <- rasterStack.CHELSA_TraCE21k.cea.i@data@names
  
  
  # Manipulate output
  (CHELSA_TraCE21k_cea30_50km <- extract.rasterStack.CHELSA_TraCE21k.i %>% 
      as.data.frame() %>%
      mutate(Quadrat_ID = world_grid_50km$Quadrat_ID) %>%
      rename_at(.vars = vars(starts_with("CHELSA_TraCE21k_")),
                .funs = funs(sub("V1.0_", "", .))) %>%
      rename_all(.funs = funs(gsub(".", "-", ., fixed = TRUE))) %>%
      #      rename_all(.funs = funs(str_replace(., "_cea30_5km.", "_bio_"))) %>%
      dplyr::relocate(Quadrat_ID) %>%
      as.data.frame()
  )
  
  # CHELSA_TraCE21k_cea30_5km %>%
  # group_by(Quadrat_ID) %>% 
  # tally() %>% 
  # filter(n != 100) %>%
  # as.data.frame()
  
  data.table::fwrite(CHELSA_TraCE21k_cea30_50km,
                     file = paste0("~/PaleoClimPoolCommDivPhyloStruct/R1/data/matrices/CHELSA_TraCE21k_500yr_cea30_50km/", 
                                   gsub('-', '_', 
                                        gsub('\\.', '_', 
                                             sub(pattern = "(.*)\\..*$", 
                                                 replacement = "\\1", 
                                                 basename(rasterFilename.i)
                                             )
                                        )), 
                                   "_ext_50km", ".csv")
  )
  
  cat(paste0(gsub('\\.', '_', 
                  sub(pattern = "(.*)\\..*$", 
                      replacement = "\\1", 
                      basename(rasterFilename.i)
                  )
  ), 
  "_ext_50km", ".csv"), 'exported! This iteration took', proc.time()[3] - counter.proc.time, 'seconds', '\n')
}


## Init snowfall
library(snowfall)

sfInit(parallel = TRUE, 
       cpus = 40,
       slaveOutfile= "~/PaleoClimPoolCommDivPhyloStruct/R1/data/matrices/CHELSA_TraCE21k_cea30_50km/log.txt") ## we select 2 CPUs. If you have 8 CPUs, put 8

## Export packages
sfLibrary('raster', character.only=TRUE)
sfLibrary('terra', character.only=TRUE)
sfLibrary('rgdal', character.only = TRUE)
sfLibrary('exactextractr', character.only=TRUE)
sfLibrary('tidyverse', character.only=TRUE)
sfLibrary('data.table', character.only=TRUE)
sfLibrary('tools', character.only=TRUE)

## Export variables
sfExport('world_grid_50km')
sfExport('equal_area_projection')
sfExport('rasterFilenames')
sfExport('extractRasterSF')
sfExport('templateBase')
sfExport('export_dir')

# you may also use sfExportAll() to export all your workspace variables

## Do the run
mySFrasterExtraction <- sfLapply(rasterFilenames, 
                                 extractRasterSF)

## stop snowfall
sfStop(nostop = FALSE)

####

(cea.rasterFilenames <- list.files(
  path = paste0(export_dir, "/cea_50km"), 
  pattern='.tif$', 
  all.files = TRUE, 
  full.names = TRUE,
  include.dirs = TRUE,
  recursive = TRUE)
)

cea.rasterFilenames.bio01 <- cea.rasterFilenames[1:45]

c(sort(cea.rasterFilenames.bio01[1:40], decreasing = T), cea.rasterFilenames[41:45])

# Load TIF as a rasterBrick
rasterStack.cea.CHELSA_TraCE21k.bio_1 <- stack(c(sort(cea.rasterFilenames.bio01[1:40], decreasing = T), cea.rasterFilenames.bio01[41:45]))

tempTrend.CHELSA_TraCE21k.bio_1<- tempTrend(rasterStack.cea.CHELSA_TraCE21k.bio_1, th = 10)
spatGrad.CHELSA_TraCE21k.bio_1 <- spatGrad(rasterStack.cea.CHELSA_TraCE21k.bio_1, th = 0.0001, projected = TRUE)

#

plot(tempTrend.CHELSA_TraCE21k.bio_1)
plot(spatGrad.CHELSA_TraCE21k.bio_1)

#

gVoCC.CHELSA_TraCE21k.bio_1 <- gVoCC(tempTrend.CHELSA_TraCE21k.bio_1,
                                     spatGrad.CHELSA_TraCE21k.bio_1)
plot(gVoCC.CHELSA_TraCE21k.bio_1)

###


(cea.rasterFilenames <- list.files(
  path = paste0(export_dir, "/cea_50km"), 
  pattern='bio12', 
  all.files = TRUE, 
  full.names = TRUE,
  include.dirs = TRUE,
  recursive = TRUE)
)

cea.rasterFilenames.

cea.rasterFilenames.bio12 <- list.files(
  path = paste0(export_dir, "/cea_50km"), 
  pattern='bio12', 
  all.files = TRUE, 
  full.names = TRUE,
  include.dirs = TRUE,
  recursive = TRUE)

c(sort(cea.rasterFilenames.bio12[1:40], decreasing = T), cea.rasterFilenames.bio12[41:45])

# Load TIF as a rasterBrick
rasterStack.cea.CHELSA_TraCE21k.bio_12 <- stack(
  c(
    sort(cea.rasterFilenames.bio12[1:40], 
         decreasing = T), 
    cea.rasterFilenames.bio12[41:45]
  )
)

rasterStack.cea.CHELSA_TraCE21k.bio_12.log <- log1p(rasterStack.cea.CHELSA_TraCE21k.bio_12)

tempTrend.CHELSA_TraCE21k.bio_12.log <- tempTrend(rasterStack.cea.CHELSA_TraCE21k.bio_12.log, th = 10)
spatGrad.CHELSA_TraCE21k.bio_12.log <- spatGrad(rasterStack.cea.CHELSA_TraCE21k.bio_12.log, th = 0.0001, projected = TRUE)

plot(tempTrend.CHELSA_TraCE21k.bio_12.log)
plot(spatGrad.CHELSA_TraCE21k.bio_12.log)

gVoCC.CHELSA_TraCE21k.bio_12.log <- gVoCC(tempTrend.CHELSA_TraCE21k.bio_12.log,
                                          spatGrad.CHELSA_TraCE21k.bio_12.log)

plot(gVoCC.CHELSA_TraCE21k.bio_12.log)
#

tempTrend.CHELSA_TraCE21k.bio_12 <- tempTrend(rasterStack.cea.CHELSA_TraCE21k.bio_12, th = 10)
spatGrad.CHELSA_TraCE21k.bio_12 <- spatGrad(rasterStack.cea.CHELSA_TraCE21k.bio_12, th = 0.0001, projected = TRUE)


plot(tempTrend.CHELSA_TraCE21k.bio_12)
plot(spatGrad.CHELSA_TraCE21k.bio_12)

#

gVoCC.CHELSA_TraCE21k.bio_12 <- gVoCC(tempTrend.CHELSA_TraCE21k.bio_12,
                                      spatGrad.CHELSA_TraCE21k.bio_12)
plot(gVoCC.CHELSA_TraCE21k.bio_12)


stack(gVoCC.CHELSA_TraCE21k.bio_1,
      gVoCC.CHELSA_TraCE21k.bio_12,
      gVoCC.CHELSA_TraCE21k.bio_12.log)


# Perform the exact extraction of the raster data to the polygon grid
extract.gVoCC.CHELSA_TraCE21k.bio_1_bio_12 <- exactextractr::exact_extract(
  stack(gVoCC.CHELSA_TraCE21k.bio_1,
        gVoCC.CHELSA_TraCE21k.bio_12,
        gVoCC.CHELSA_TraCE21k.bio_12.log), 
  world_grid_50km,
  include_xy = FALSE, # include coordinates
  include_cell = FALSE, # include cell IDs
  fun = "mean",
  force_df = TRUE,
  max_cells_in_memory = 3e+08# ,
  # include_cols = rasterStack.CHELSA_TraCE21k.cea.i@data@names
  # full_colnames = T
)

colnames(extract.gVoCC.CHELSA_TraCE21k.bio_1_bio_12) <- c("voccMag_bio_1", "voccAng_bio_1",
                                                          "voccMag_bio_12", "voccAng_bio_12",
                                                          "voccMag_log_bio_12", "voccAng_log_bio_12")

  
  # Manipulate output
  (gVoCC.CHELSA_TraCE21k.bio_1_bio_12_50km <- extract.gVoCC.CHELSA_TraCE21k.bio_1_bio_12 %>% 
      as.data.frame() %>%
      mutate(Quadrat_ID = world_grid_50km$Quadrat_ID) %>%
      #      rename_all(.funs = funs(str_replace(., "_cea30_5km.", "_bio_"))) %>%
      dplyr::relocate(Quadrat_ID) %>%
      as.data.frame()
  )

# CHELSA_TraCE21k_cea30_5km %>%
# group_by(Quadrat_ID) %>% 
# tally() %>% 
# filter(n != 100) %>%
# as.data.frame()

data.table::fwrite(extract.gVoCC.CHELSA_TraCE21k.bio_1_bio_12,
                   file = paste0("~/PaleoClimPoolCommDivPhyloStruct/R1/data/matrices/", 
                                 "extract.gVoCC.CHELSA_TraCE21k.bio_1_bio_12_500yr",
                                 "_ext_50km", ".csv")
)

