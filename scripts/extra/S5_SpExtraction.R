# Routine to extract a community presence-absence matrix from IUCN species distributions
# using lets.presab.grid(). This was not used for this project.

require(sp)
require(rgdal)
require(snowfall)

mammalsExtentIUCN <- readOGR(dsn = "data/organisms", layer = "TERRESTRIAL_MAMMALS") 
world_grid_05 <- readOGR(dsn = "data/grids", layer = "world_05_grid")

world_grid_05$ID <- c(1:length(world_grid_05))

proj4string(world_grid_05) <- proj4string(mammalsExtentIUCN)

# projection(mammalsExtentIUCN) <- projection(world_grid_05)

uniqueSpecies <- as.vector(unique(mammalsExtentIUCN@data$binomial))

spExtraction <- function(sp.uniquenames){
  sp.i <- mammalsExtentIUCN[mammalsExtentIUCN$binomial == sp.uniquenames, ]
  
  presAbs.sp.i <- lets.presab.grid(sp.i, world_grid_05, "ID")$PAM
  
  write.table(sp.uniquenames, file=paste("data/progress/", sp.uniquenames, ".done", sep=""))
  return(as.data.frame(presAbs.sp.i))
}

# Check if it is working with:
# system.time(a <- lapply(uniqueSpecies, spExtraction))

library(snowfall)
sfInit(parallel = TRUE, cpus = 24) ## we select 2 CPUs. If you have 8 CPUs, put 8.

## Export packages
sfLibrary('letsR', character.only=TRUE)

## Export variables
sfExport('mammalsExtentIUCN')
sfExport('uniqueSpecies')
sfExport('world_grid_05')
# you may also use sfExportAll() to export all your workspace variables

## Do the run
mySFPresAbs <- sfLapply(uniqueSpecies, spExtraction)

## stop snowfall
sfStop(nostop=FALSE )

### Building presence-absence matrix
current_PresAbs <- as.data.frame(mySFPresAbs)

write.table(current_PresAbs, "data/matrices/world_mammals_PAM_05.txt")
