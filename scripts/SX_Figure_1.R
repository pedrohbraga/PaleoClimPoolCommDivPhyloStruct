#################################################################################
###            ###
###                         ###
#                                                                               #
# Author: Pedro Henrique Pereira Braga                                          #
# Last Update: "2020-05-13"                                                     #
#                                                                               # 
#################################################################################

# Load world projection grid shapefile
polygonGrid_world_50km <- st_read("data/grids/world_grid_prev_50km.shp")

# Attribute an "id" variable
polygonGrid_world_50km$id = 1:length(polygonGrid_world_50km)

# Merge the MPD.MNTD data with the grid
polygonGrid_world_50km_merged_sf <- polygonGrid_world_50km %>%
  merge(filter(MPD.MNTD.LatLong.AllScales,
               is.na(ID_Realm) == FALSE,
               # SamplingPool == "Global sampling"
  ) %>%
    mutate(QuadratID = ID -1),
  by = c("QuadratID" = "QuadratID")
  )

# Gather variables towards columns

long.polygonGrid_world_50km_merged_sf <- polygonGrid_world_50km_merged_sf %>%
  gather(key = "PhyloMetric",
         value = "PhyloStructure",
         c(NRI, NTI)
  )

# head(long.polygonGrid_world_50km_merged_sf)

# Plot the information using geom_sf() and ggplot()

ggplot(data = long.polygonGrid_world_50km_merged_sf) +
  geom_sf(mapping = aes(fill = PhyloStructure),
          lwd = 0 # no quadrat boundaries 
  ) +
  scico::scale_colour_scico(palette = "acton", direction = -1) +
  # scale_fill_gradient2( # Divergent colour scale 
  #   name = "Phylogenetic Structure",
  #   low = muted("blue"),
  #   mid = "white",
  #   high = muted("red"),
  #   midpoint = 0,
  #   na.value = "lightgrey",
  #   breaks = breaks_extended(9),
  #   guide = guide_colorbar(
  #     direction = "horizontal",
  #     barheight = unit(2, units = "mm"),
  #     barwidth = unit(75, units = "mm"),
  #     draw.ulim = F, 
  #     title.position = 'top', 
  #     title.hjust = 0.5,  # shifting title height 
  #     label.hjust = 0.5   # shifting label height
  #   ) # FALSE will suppress legend
  # ) +
  facet_grid(SamplingPool ~ PhyloMetric,
             switch = "x") + # facet_grid top labels on bottom
  # add theme
  theme_map() + 
  theme(legend.position = "bottom",
        strip.text.y = element_text(size = 10),
        strip.text.x = element_text(size = 10)
        )

# Export figure

ggsave(filename = "Figure.1.NRI.NTI.Sampling.maps.png", 
       dpi = 300, 
       width = 8.5 , height = 11.5, 
       units = "in")

