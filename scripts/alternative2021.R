
#################################
########### Visualize ###########
#################################

# We will create a map using sp and maptools, because an issue prevented me from doing it with "sf"
# This implies in an slower procedure, but it will work.

# Load the polygon grid
polygonGrid_world_50km <- readOGR(dsn="data/grids/", 
                                  layer="world_grid_prev_50km")

# Load world projection grid shapefile
polygonGrid_world_50km <- st_read("data/grids/world_grid_prev_50km.shp")


# Define its projection
proj4string(polygonGrid_world_50km) <- "+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs"

# Attribute an "id" variable
polygonGrid_world_50km@data$id = 1:length(polygonGrid_world_50km)

# st_as_sf()
polygonGrid_world_50km.sf <- st_as_sf(polygonGrid_world_50km)
class(polygonGrid_world_50km.sf)



polygonGrid_world_50km_merged_sf<- polygonGrid_world_50km.sf %>%
  merge(filter(MPD.MNTD.LatLong.AllScales,
               is.na(ID_Realm) == FALSE,
               #               SamplingPool == "Global sampling"
  ) %>%
    mutate(QuadratID = ID -1),
  by = c("QuadratID" = "QuadratID")
  )

# Gather

long.polygonGrid_world_50km_merged_sf <- polygonGrid_world_50km_merged_sf %>%
  gather(key = "PhyloMetric",
         value = "PhyloStructure",
         c(NRI, NTI))

# Plot the information using geom_sf() and ggplot()

ggplot(data = long.polygonGrid_world_50km_merged_sf) +
  geom_sf(mapping = aes(fill = PhyloStructure),
          lwd = 0 # no quadrat boundaries  
  ) +
  scale_fill_gradientn(colours = extended.colPalette, 
                       name = "Phylogenetic Structure", 
                       guide = guide_colorbar(
                         direction = "horizontal",
                         barheight = unit(2, units = "mm"),
                         barwidth = unit(75, units = "mm"),
                         draw.ulim = F,
                         title.position = 'top',
                         # some shifting around
                         title.hjust = 0.5,
                         label.hjust = 0.5
                       ) # FALSE will suppress legend
  ) +
  facet_grid(SamplingPool ~ PhyloMetric) +
  # add theme
  theme_map() + 
  theme(legend.position="bottom")  

long.polygonGrid_world_50km_merged_sf

# Make it "ggplot2-friendly" by converting it to data.frame and using broom::tidy()
# polygonGrid_world_50km.points <- broom::tidy(polygonGrid_world_50km)

# Explore to understand the result of fortifying: 
# head(polygonGrid_world_50km.points)

CWM.Div.MPD.Chiroptera.Comm$QuadratID <-  (CWM.Div.MPD.Chiroptera.Comm$ID - 1)

# Join CRI data to the fortified data
polygonGrid_world_50km.df <- polygonGrid_world_50km.points %>% 
  mutate_each(list(as.numeric), id) %>%
  left_join(CWM.Div.MPD.Chiroptera.Comm %>%
              filter(SamplingPool == "Global sampling"), 
            by = c("id" = "QuadratID"))

head(polygonGrid_world_50km.df)
head(CWM.Div.MPD.Chiroptera.Comm)

# Inspect: head(polygonGrid_world_50km.df)

#### Straightfoward visualisation ####

#pal <- wes_palette("Zissou1", 100, type = "continuous")
#pal <- wes_palette("Moonrise3", 100, type = "continuous")

raw.colPalette <- brewer.pal(9, 'YlOrBr')
pie(rep(1,9), col = raw.colPalette)

newcol <- colorRampPalette(raw.colPalette)
ncols <- 100

extended.colPalette <- newcol(ncols) #apply the function to get 100 colours

pie(rep(1, ncols), 
    col = extended.colPalette, 
    border = NA, 
    labels = NA)

p <- ggplot() +
  # municipality polygons
  geom_polygon(data = polygonGrid_world_50km.df, aes(fill = netDiv_CWM_std_w, 
                                                     x = long, 
                                                     y = lat, 
                                                     group = group)) +
  # no outline
  #geom_path(data = map_data, aes(x = long, 
  #                               y = lat, 
  #                               group = group), 
  #          color = "white", size = 0.1) +
  
  coord_equal() +
  
  # add a basic theme
  theme_map() +
  
  theme(legend.position = "bottom") +
  
  scale_fill_gradientn(colours = extended.colPalette,
                       name = "Net Diversification Rates", 
                       guide = guide_colorbar(
                         direction = "horizontal",
                         barheight = unit(2, units = "mm"),
                         barwidth = unit(75, units = "mm"),
                         draw.ulim = F,
                         title.position = 'top',
                         # some shifting around
                         title.hjust = 0.5,
                         label.hjust = 0.5
                       )) + 
  # insert a North arrow
  ggsn::north(polygonGrid_world_50km.df, 
              symbol = 12, scale = 0.04,
              location = "bottomright",
              anchor = c(x = +14500000, y = -6850000)) +
  
  labs(x = NULL, 
       y = NULL, 
       title = "Net Diversification Rates", 
       subtitle = "Community Weighted Mean",
       caption = "Bats")
p
