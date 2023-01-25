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

