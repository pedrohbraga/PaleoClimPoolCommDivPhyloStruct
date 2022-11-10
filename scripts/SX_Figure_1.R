#################################################################################
### Map representation of the globe-wise phylogenetic structure of bat        ###
### communities across species sampling pools                                 ###
##                                                                             ##
#                                                                               #
# Author: Pedro Henrique Pereira Braga                                          #
# Last Update: "2022-11-09"                                                     #
#                                                                               # 
#################################################################################

# Load world projection grid shapefile
polygonGrid_world_50km <- st_read("data/grids/world_grid_prev_50km.shp")

# Attribute an "id" variable
polygonGrid_world_50km$id = 1:length(polygonGrid_world_50km)

# Merge the MPD.MNTD data with the grid
polygonGrid_world_50km_merged_sf <- polygonGrid_world_50km %>%
  merge(filter(MPD.MNTD.LatLong.AllScales.raref.min.relative.worldClimate.diff.CWM.Div,
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
         c(nri, nti)
  )

# head(long.polygonGrid_world_50km_merged_sf)

# Choose colours, play with luminance and chroma values

show_col(c("red", "blue", "lightgray",
           muted("red", l = 20, c = 90), 
           muted("blue", l = 20, c = 90),
           "lightgray",
           muted("red"),
           muted("blue"),
           "lightgray"
)
)

# Plot the information using geom_sf() and ggplot()

nri.AllScales.map <- ggplot(data = long.polygonGrid_world_50km_merged_sf %>%
                              filter(PhyloMetric == "nri")) +
  geom_sf(mapping = aes(fill = PhyloStructure),
          colour = NA
          # lwd = 0, # no quadrat boundaries
  ) +
  # scico::scale_colour_scico(palette = "acton", direction = -1) +
  scale_fill_gradient2( # Divergent colour scale
    name = "NRI",
    low = muted("blue",
                l = 20,
                c = 90),
    mid = "white",
    high = muted("red",
                 l = 20,
                 c = 90),
    midpoint = 0,
    na.value = "lightgrey",
    breaks = breaks_extended(10),
    guide = guide_colorbar(
      direction = "horizontal",
      barheight = unit(2, units = "mm"),
      barwidth = unit(75, units = "mm"),
      draw.ulim = F,
      title.position = 'top',
      title.hjust = 0.5,  # shifting title height
      label.hjust = 0.5   # shifting label height
    ) # FALSE will suppress legend
  ) +
  # colorspace::scale_fill_continuous_diverging(palette = "Blue-Red 3", 
  #                                 l1 = 30, 
  #                                 l2 = 100, 
  #                                 p1 = .9, 
  #                                 p2 = 1.2) +
  facet_grid(SamplingPool ~ .,
             switch = "x") + # facet_grid top labels on bottom
  # add theme
  theme_minimal() +
  theme_map() + 
  theme(legend.position = "bottom",
        strip.text.y = element_text(size = 10),
        # strip.text.y = element_blank(),
        strip.text.x = element_text(size = 10),
        panel.grid.major = element_line(colour = "transparent")
  )


# Plot the information using geom_sf() and ggplot()

nti.AllScales.map <- ggplot(data = long.polygonGrid_world_50km_merged_sf %>%
                              filter(PhyloMetric == "nti")) +
  geom_sf(mapping = aes(fill = PhyloStructure),
          colour = NA
          # lwd = 0, # no quadrat boundaries
  ) +
  # scico::scale_colour_scico(palette = "acton", direction = -1) +
  scale_fill_gradient2( # Divergent colour scale
    name = "NTI",
    low = muted("blue",
                l = 20,
                c = 90),
    mid = "white",
    high = muted("red",
                 l = 20,
                 c = 90),
    midpoint = 0,
    na.value = "lightgrey",
    breaks = breaks_extended(9),
    guide = guide_colorbar(
      direction = "horizontal",
      barheight = unit(2, units = "mm"),
      barwidth = unit(75, units = "mm"),
      draw.ulim = F,
      title.position = 'top',
      title.hjust = 0.5,  # shifting title height
      label.hjust = 0.5   # shifting label height
    ) # FALSE will suppress legend
  ) +
  # colorspace::scale_fill_continuous_diverging(palette = "Blue-Red 3", 
  #                                 l1 = 30, 
  #                                 l2 = 100, 
  #                                 p1 = .9, 
  #                                 p2 = 1.2) +
  facet_grid(SamplingPool ~ .,
             switch = "x") + # facet_grid top labels on bottom
  # add theme
  theme_minimal() +
  theme_map() + 
  theme(legend.position = "bottom",
        strip.text.y = element_text(size = 10),
        strip.text.x = element_text(size = 10),
        panel.grid.major = element_line(colour = "transparent")
  )


ggarrange(nri.AllScales.map + 
            theme(strip.text.y =  element_text(size = 10,
                                               colour = "transparent")),
          nti.AllScales.map,
          ncol = 2,
          nrow = 1)

# Export Figure 1

ggplot2::ggsave(filename = "figures/Figure.1.NRI.NTI.Sampling.maps.png", 
       dpi = 600, 
       width = 9 , height = 13, 
       units = "in",
#       device = "png",
       type = "agg_png")
