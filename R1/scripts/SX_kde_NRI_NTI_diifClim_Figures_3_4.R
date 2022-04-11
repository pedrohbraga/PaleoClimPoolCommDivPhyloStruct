# Function

ks_kde_geom_contours_phyloStr_diffClim <- function(parent_data = MPD.MNTD.CWM.Div.diff.worldClimate.subset.data,
                                                   var_x = "diff.AnnTemp.LGM_cur",
                                                   var_y = "nri",
                                                   lab_x = c(expression(atop(
                                                     "Historical Change in Temperature (°C)",
                                                     scriptstyle("MAT"[Contemporary] - "MAT"[LGM])
                                                   ))),
                                                   lab_y = c(expression("NRI"["Global"]))
){
  
  phyloStr.diff.Clim.subset.data <- parent_data %>%
    select({{var_x}},
           {{var_y}}) %>%
    drop_na()
  
  
  subset.var_y.var_x.Hpi <- ks::Hpi(
    x =  phyloStr.diff.Clim.subset.data,
    pilot = "samse"
  )
  
  # subset.var_y.var_x.Hpi
  
  # Compute KDE
  
  kde.subset.var_y.var_x.Hpi <- ks::kde(
    x =  phyloStr.diff.Clim.subset.data,
    H = subset.var_y.var_x.Hpi,
    binned = TRUE,
    compute.cont = TRUE,
    approx.cont = TRUE
  )
  
  # Represent with ks::plot.kde()
  
  # "cont" specifies the density contours, which are upper percentages of highest
  # density regions. The default contours are at 25%, 50%, and 75%
  
  percentiles <- percent(c(
    seq(0.01, 0.99,
        length.out = 99
    )
  ),
  accuracy = 1)
  
  # percentiles <- percent(c(
  #   0.01,
  #   seq(0.01, 0.9,
  #       length.out = 98
  #   ),
  #   0.99
  # ),
  # accuracy = 1)
  # 
  
    # plot(kde.subset.var_y.var_x.Hpi,
  #      display = "filled.contour2",
  #      approx = FALSE,
  #      cont = as.numeric(sub("%", "", 
  #                            percentiles,
  #                            fixed = TRUE)))
  
  # plot(kde.subset.var_y.var_x.Hpi, display = "persp")
  
  ## Subset and prepare information to represent
  
  fhat.kde.subset.var_y.var_x.Hpi <- kde.subset.var_y.var_x.Hpi
  
  dimnames(fhat.kde.subset.var_y.var_x.Hpi[["estimate"]]) <- list(
    fhat.kde.subset.var_y.var_x.Hpi[["eval.points"]][[1]],
    fhat.kde.subset.var_y.var_x.Hpi[["eval.points"]][[2]]
  )
  
  molten.fhat.kde.subset.var_y.var_x.Hpi <- reshape2::melt(fhat.kde.subset.var_y.var_x.Hpi[["estimate"]],
                                                           value.name = "estimate") %>%
    rename(
      var_x = Var1,
      var_y = Var2
    )
  
  # Represent contours with ggplot() and geom_contour
  
  fig.kde.subset.var_y.var_x.Hpi <- ggplot(
    molten.fhat.kde.subset.var_y.var_x.Hpi,
    aes(
      x = var_x,
      y = var_y
    )
  ) +
    geom_contour_filled(aes(z = estimate),
                        breaks = 
                          fhat.kde.subset.var_y.var_x.Hpi[["cont"]][percentiles]
    ) +
    scico::scale_fill_scico_d(palette = "acton", 
                              direction = -1, 
                              drop = FALSE) +
    labs(
      x = lab_x,
      y = lab_y) +
    geom_hline(yintercept = 0, alpha = 0.25) +
    scale_y_continuous(
      breaks = pretty_breaks(6,
                             high.u.bias = 0.1,
                             eps.correct = 0,
                             min.n = 6
      ),
      limits = c(-0.5, 0.5) + range(molten.fhat.kde.subset.var_y.var_x.Hpi$var_y,
                                    na.rm = TRUE
      )
    ) +
    scale_x_continuous(
      limits = c(-0.1, 0.1) + range(molten.fhat.kde.subset.var_y.var_x.Hpi$var_x,
                                    na.rm = TRUE
      )) +
    theme_classic(base_size = 17 * 1.6) +
    theme(
      legend.position = "none",
      legend.direction = "horizontal",
      legend.key.size = unit(0.5, "cm"),
      legend.text = element_blank(),
      legend.title = element_blank(),
      legend.box = "horizontal",
      # axis.title.x = element_blank(),
      # axis.text.x = element_blank(),
      # axis.text.y = element_text(face = "bold", size = 14),
      axis.title.y = element_text(
        size = 17 * 1.7,
        face = "bold"
      ),
      axis.title.x = element_text(
        size = 17 * 1.6,
        face = "bold"
      )
    ) +
    guides(colour = guide_legend(nrow = 1))
  return(fig.kde.subset.var_y.var_x.Hpi)
}

# Global -------------------------------------------------

saveRDS(MPD.MNTD.Chiroptera.Comm.AllScales.CWM.Div,
        "data/matrices/MPD.MNTD.Chiroptera.Comm.AllScales.CWM.Div.rds")

saveRDS(MPD.MNTD.LatLong.AllScales.worldClimate.diff,
        "data/matrices/MPD.MNTD.LatLong.AllScales.worldClimate.diff.rds")

## NRI and Annual Temperature and Precipitation ----------

fig.kde.Global.NRI.diff.AnnTemp.LGM_cur.Hpi <- ks_kde_geom_contours_phyloStr_diffClim(parent_data = MPD.MNTD.Chiroptera.Comm.AllScales.CWM.Div %>%
                                                                                        filter(SamplingPool == "Global sampling") %>%
                                                                                        left_join(MPD.MNTD.LatLong.AllScales.worldClimate.diff %>%
                                                                                                    filter(SamplingPool == "Global sampling") %>%
                                                                                                    select(
                                                                                                      ID,
                                                                                                      diff.AnnTemp.LGM_cur:diff.log.AnnPrec.log.LGM_cur
                                                                                                    ),
                                                                                                  by = "ID"
                                                                                        ),
                                                                                      var_x = "diff.AnnTemp.LGM_cur",
                                                                                      var_y = "nri",
                                                                                      lab_x = c(expression(atop(
                                                                                        "Historical Change in Temperature (°C)",
                                                                                        scriptstyle("MAT"[Contemporary] - "MAT"[LGM])
                                                                                      ))),
                                                                                      lab_y = c(expression("NRI"["Global"]))
)


fig.kde.Global.NRI.diff.log.AnnPrec.log.LGM_cur.Hpi <- ks_kde_geom_contours_phyloStr_diffClim(parent_data = MPD.MNTD.Chiroptera.Comm.AllScales.CWM.Div %>%
                                                                                                filter(SamplingPool == "Global sampling") %>%
                                                                                                left_join(MPD.MNTD.LatLong.AllScales.worldClimate.diff %>%
                                                                                                            filter(SamplingPool == "Global sampling") %>%
                                                                                                            select(
                                                                                                              ID,
                                                                                                              diff.AnnTemp.LGM_cur:diff.log.AnnPrec.log.LGM_cur
                                                                                                            ),
                                                                                                          by = "ID"
                                                                                                ),
                                                                                              var_x = "diff.log.AnnPrec.log.LGM_cur",
                                                                                              var_y = "nri",
                                                                                              lab_x = c(expression(atop(
                                                                                                "Historical Change in Precipitation",
                                                                                                scriptstyle(log("MAP"[Contemporary]) - log("MAP"[LGM]))
                                                                                              ))),
                                                                                              lab_y = c(expression("NRI"["Global"]))
)


## NTI and Annual Temperature and Precipitation ----------


fig.kde.Global.NTI.diff.AnnTemp.LGM_cur.Hpi <- ks_kde_geom_contours_phyloStr_diffClim(parent_data = MPD.MNTD.Chiroptera.Comm.AllScales.CWM.Div %>%
                                                                                        filter(SamplingPool == "Global sampling") %>%
                                                                                        left_join(MPD.MNTD.LatLong.AllScales.worldClimate.diff %>%
                                                                                                    filter(SamplingPool == "Global sampling") %>%
                                                                                                    select(
                                                                                                      ID,
                                                                                                      diff.AnnTemp.LGM_cur:diff.log.AnnPrec.log.LGM_cur
                                                                                                    ),
                                                                                                  by = "ID"
                                                                                        ),
                                                                                      var_x = "diff.AnnTemp.LGM_cur",
                                                                                      var_y = "nti",
                                                                                      lab_x = c(expression(atop(
                                                                                        "Historical Change in Temperature (°C)",
                                                                                        scriptstyle("MAT"[Contemporary] - "MAT"[LGM])
                                                                                      ))),
                                                                                      lab_y = c(expression("NTI"["Global"]))
)


fig.kde.Global.NTI.diff.log.AnnPrec.log.LGM_cur.Hpi <- ks_kde_geom_contours_phyloStr_diffClim(parent_data = MPD.MNTD.Chiroptera.Comm.AllScales.CWM.Div %>%
                                                                                                filter(SamplingPool == "Global sampling") %>%
                                                                                                left_join(MPD.MNTD.LatLong.AllScales.worldClimate.diff %>%
                                                                                                            filter(SamplingPool == "Global sampling") %>%
                                                                                                            select(
                                                                                                              ID,
                                                                                                              diff.AnnTemp.LGM_cur:diff.log.AnnPrec.log.LGM_cur
                                                                                                            ),
                                                                                                          by = "ID"
                                                                                                ),
                                                                                              var_x = "diff.log.AnnPrec.log.LGM_cur",
                                                                                              var_y = "nti",
                                                                                              lab_x = c(expression(atop(
                                                                                                "Historical Change in Precipitation",
                                                                                                scriptstyle(log("MAP"[Contemporary]) - log("MAP"[LGM]))
                                                                                              ))),
                                                                                              lab_y = c(expression("NTI"["Global"]))
)

# Hemispheric ----------------------------------------------

## NRI and Annual Temperature and Precipitation ------------

fig.kde.Hemispheric.NRI.diff.AnnTemp.LGM_cur.Hpi <- ks_kde_geom_contours_phyloStr_diffClim(parent_data = MPD.MNTD.Chiroptera.Comm.AllScales.CWM.Div %>%
                                                                                             filter(SamplingPool == "Hemispheric sampling") %>%
                                                                                             left_join(MPD.MNTD.LatLong.AllScales.worldClimate.diff %>%
                                                                                                         filter(SamplingPool == "Hemispheric sampling") %>%
                                                                                                         select(
                                                                                                           ID,
                                                                                                           diff.AnnTemp.LGM_cur:diff.log.AnnPrec.log.LGM_cur
                                                                                                         ),
                                                                                                       by = "ID"
                                                                                             ),
                                                                                           var_x = "diff.AnnTemp.LGM_cur",
                                                                                           var_y = "nri",
                                                                                           lab_x = c(expression(atop(
                                                                                             "Historical Change in Temperature (°C)",
                                                                                             scriptstyle("MAT"[Contemporary] - "MAT"[LGM])
                                                                                           ))),
                                                                                           lab_y = c(expression("NRI"["Hemispheric"]))
)


fig.kde.Hemispheric.NRI.diff.log.AnnPrec.log.LGM_cur.Hpi <- ks_kde_geom_contours_phyloStr_diffClim(parent_data = MPD.MNTD.Chiroptera.Comm.AllScales.CWM.Div %>%
                                                                                                     filter(SamplingPool == "Hemispheric sampling") %>%
                                                                                                     left_join(MPD.MNTD.LatLong.AllScales.worldClimate.diff %>%
                                                                                                                 filter(SamplingPool == "Hemispheric sampling") %>%
                                                                                                                 select(
                                                                                                                   ID,
                                                                                                                   diff.AnnTemp.LGM_cur:diff.log.AnnPrec.log.LGM_cur
                                                                                                                 ),
                                                                                                               by = "ID"
                                                                                                     ),
                                                                                                   var_x = "diff.log.AnnPrec.log.LGM_cur",
                                                                                                   var_y = "nri",
                                                                                                   lab_x = c(expression(atop(
                                                                                                     "Historical Change in Precipitation",
                                                                                                     scriptstyle(log("MAP"[Contemporary]) - log("MAP"[LGM]))
                                                                                                   ))),
                                                                                                   lab_y = c(expression("NRI"["Hemispheric"]))
)


## NTI and Annual Temperature and Precipitation ------------


fig.kde.Hemispheric.NTI.diff.AnnTemp.LGM_cur.Hpi <- ks_kde_geom_contours_phyloStr_diffClim(parent_data = MPD.MNTD.Chiroptera.Comm.AllScales.CWM.Div %>%
                                                                                             filter(SamplingPool == "Hemispheric sampling") %>%
                                                                                             left_join(MPD.MNTD.LatLong.AllScales.worldClimate.diff %>%
                                                                                                         filter(SamplingPool == "Hemispheric sampling") %>%
                                                                                                         select(
                                                                                                           ID,
                                                                                                           diff.AnnTemp.LGM_cur:diff.log.AnnPrec.log.LGM_cur
                                                                                                         ),
                                                                                                       by = "ID"
                                                                                             ),
                                                                                           var_x = "diff.AnnTemp.LGM_cur",
                                                                                           var_y = "nti",
                                                                                           lab_x = c(expression(atop(
                                                                                             "Historical Change in Temperature (°C)",
                                                                                             scriptstyle("MAT"[Contemporary] - "MAT"[LGM])
                                                                                           ))),
                                                                                           lab_y = c(expression("NTI"["Hemispheric"]))
)


fig.kde.Hemispheric.NTI.diff.log.AnnPrec.log.LGM_cur.Hpi <- ks_kde_geom_contours_phyloStr_diffClim(parent_data = MPD.MNTD.Chiroptera.Comm.AllScales.CWM.Div %>%
                                                                                                     filter(SamplingPool == "Hemispheric sampling") %>%
                                                                                                     left_join(MPD.MNTD.LatLong.AllScales.worldClimate.diff %>%
                                                                                                                 filter(SamplingPool == "Hemispheric sampling") %>%
                                                                                                                 select(
                                                                                                                   ID,
                                                                                                                   diff.AnnTemp.LGM_cur:diff.log.AnnPrec.log.LGM_cur
                                                                                                                 ),
                                                                                                               by = "ID"
                                                                                                     ),
                                                                                                   var_x = "diff.log.AnnPrec.log.LGM_cur",
                                                                                                   var_y = "nti",
                                                                                                   lab_x = c(expression(atop(
                                                                                                     "Historical Change in Precipitation",
                                                                                                     scriptstyle(log("MAP"[Contemporary]) - log("MAP"[LGM]))
                                                                                                   ))),
                                                                                                   lab_y = c(expression("NTI"["Hemispheric"]))
)

# Realm ------------

## NRI and Annual Temperature and Precipitation ------------

fig.kde.Realm.NRI.diff.AnnTemp.LGM_cur.Hpi <- ks_kde_geom_contours_phyloStr_diffClim(parent_data = MPD.MNTD.Chiroptera.Comm.AllScales.CWM.Div %>%
                                                                                       filter(SamplingPool == "Realm sampling") %>%
                                                                                       left_join(MPD.MNTD.LatLong.AllScales.worldClimate.diff %>%
                                                                                                   filter(SamplingPool == "Realm sampling") %>%
                                                                                                   select(
                                                                                                     ID,
                                                                                                     diff.AnnTemp.LGM_cur:diff.log.AnnPrec.log.LGM_cur
                                                                                                   ),
                                                                                                 by = "ID"
                                                                                       ),
                                                                                     var_x = "diff.AnnTemp.LGM_cur",
                                                                                     var_y = "nri",
                                                                                     lab_x = c(expression(atop(
                                                                                       "Historical Change in Temperature (°C)",
                                                                                       scriptstyle("MAT"[Contemporary] - "MAT"[LGM])
                                                                                     ))),
                                                                                     lab_y = c(expression("NRI"["Realm"]))
)


fig.kde.Realm.NRI.diff.log.AnnPrec.log.LGM_cur.Hpi <- ks_kde_geom_contours_phyloStr_diffClim(parent_data = MPD.MNTD.Chiroptera.Comm.AllScales.CWM.Div %>%
                                                                                               filter(SamplingPool == "Realm sampling") %>%
                                                                                               left_join(MPD.MNTD.LatLong.AllScales.worldClimate.diff %>%
                                                                                                           filter(SamplingPool == "Realm sampling") %>%
                                                                                                           select(
                                                                                                             ID,
                                                                                                             diff.AnnTemp.LGM_cur:diff.log.AnnPrec.log.LGM_cur
                                                                                                           ),
                                                                                                         by = "ID"
                                                                                               ),
                                                                                             var_x = "diff.log.AnnPrec.log.LGM_cur",
                                                                                             var_y = "nri",
                                                                                             lab_x = c(expression(atop(
                                                                                               "Historical Change in Precipitation",
                                                                                               scriptstyle(log("MAP"[Contemporary]) - log("MAP"[LGM]))
                                                                                             ))),
                                                                                             lab_y = c(expression("NRI"["Realm"]))
)


## NTI and Annual Temperature and Precipitation ------------


fig.kde.Realm.NTI.diff.AnnTemp.LGM_cur.Hpi <- ks_kde_geom_contours_phyloStr_diffClim(parent_data = MPD.MNTD.Chiroptera.Comm.AllScales.CWM.Div %>%
                                                                                       filter(SamplingPool == "Realm sampling") %>%
                                                                                       left_join(MPD.MNTD.LatLong.AllScales.worldClimate.diff %>%
                                                                                                   filter(SamplingPool == "Realm sampling") %>%
                                                                                                   select(
                                                                                                     ID,
                                                                                                     diff.AnnTemp.LGM_cur:diff.log.AnnPrec.log.LGM_cur
                                                                                                   ),
                                                                                                 by = "ID"
                                                                                       ),
                                                                                     var_x = "diff.AnnTemp.LGM_cur",
                                                                                     var_y = "nti",
                                                                                     lab_x = c(expression(atop(
                                                                                       "Historical Change in Temperature (°C)",
                                                                                       scriptstyle("MAT"[Contemporary] - "MAT"[LGM])
                                                                                     ))),
                                                                                     lab_y = c(expression("NTI"["Realm"]))
)


fig.kde.Realm.NTI.diff.log.AnnPrec.log.LGM_cur.Hpi <- ks_kde_geom_contours_phyloStr_diffClim(parent_data = MPD.MNTD.Chiroptera.Comm.AllScales.CWM.Div %>%
                                                                                               filter(SamplingPool == "Realm sampling") %>%
                                                                                               left_join(MPD.MNTD.LatLong.AllScales.worldClimate.diff %>%
                                                                                                           filter(SamplingPool == "Realm sampling") %>%
                                                                                                           select(
                                                                                                             ID,
                                                                                                             diff.AnnTemp.LGM_cur:diff.log.AnnPrec.log.LGM_cur
                                                                                                           ),
                                                                                                         by = "ID"
                                                                                               ),
                                                                                             var_x = "diff.log.AnnPrec.log.LGM_cur",
                                                                                             var_y = "nti",
                                                                                             lab_x = c(expression(atop(
                                                                                               "Historical Change in Precipitation",
                                                                                               scriptstyle(log("MAP"[Contemporary]) - log("MAP"[LGM]))
                                                                                             ))),
                                                                                             lab_y = c(expression("NTI"["Realm"]))
)

# Plate ------------

## NRI and Annual Temperature and Precipitation ------------

fig.kde.Plate.NRI.diff.AnnTemp.LGM_cur.Hpi <- ks_kde_geom_contours_phyloStr_diffClim(parent_data = MPD.MNTD.Chiroptera.Comm.AllScales.CWM.Div %>%
                                                                                       filter(SamplingPool == "Plate sampling") %>%
                                                                                       left_join(MPD.MNTD.LatLong.AllScales.worldClimate.diff %>%
                                                                                                   filter(SamplingPool == "Plate sampling") %>%
                                                                                                   select(
                                                                                                     ID,
                                                                                                     diff.AnnTemp.LGM_cur:diff.log.AnnPrec.log.LGM_cur
                                                                                                   ),
                                                                                                 by = "ID"
                                                                                       ),
                                                                                     var_x = "diff.AnnTemp.LGM_cur",
                                                                                     var_y = "nri",
                                                                                     lab_x = c(expression(atop(
                                                                                       "Historical Change in Temperature (°C)",
                                                                                       scriptstyle("MAT"[Contemporary] - "MAT"[LGM])
                                                                                     ))),
                                                                                     lab_y = c(expression("NRI"["Plate"]))
)


fig.kde.Plate.NRI.diff.log.AnnPrec.log.LGM_cur.Hpi <- ks_kde_geom_contours_phyloStr_diffClim(parent_data = MPD.MNTD.Chiroptera.Comm.AllScales.CWM.Div %>%
                                                                                               filter(SamplingPool == "Plate sampling") %>%
                                                                                               left_join(MPD.MNTD.LatLong.AllScales.worldClimate.diff %>%
                                                                                                           filter(SamplingPool == "Plate sampling") %>%
                                                                                                           select(
                                                                                                             ID,
                                                                                                             diff.AnnTemp.LGM_cur:diff.log.AnnPrec.log.LGM_cur
                                                                                                           ),
                                                                                                         by = "ID"
                                                                                               ),
                                                                                             var_x = "diff.log.AnnPrec.log.LGM_cur",
                                                                                             var_y = "nri",
                                                                                             lab_x = c(expression(atop(
                                                                                               "Historical Change in Precipitation",
                                                                                               scriptstyle(log("MAP"[Contemporary]) - log("MAP"[LGM]))
                                                                                             ))),
                                                                                             lab_y = c(expression("NRI"["Plate"]))
)


## NTI and Annual Temperature and Precipitation ------------


fig.kde.Plate.NTI.diff.AnnTemp.LGM_cur.Hpi <- ks_kde_geom_contours_phyloStr_diffClim(parent_data = MPD.MNTD.Chiroptera.Comm.AllScales.CWM.Div %>%
                                                                                       filter(SamplingPool == "Plate sampling") %>%
                                                                                       left_join(MPD.MNTD.LatLong.AllScales.worldClimate.diff %>%
                                                                                                   filter(SamplingPool == "Plate sampling") %>%
                                                                                                   select(
                                                                                                     ID,
                                                                                                     diff.AnnTemp.LGM_cur:diff.log.AnnPrec.log.LGM_cur
                                                                                                   ),
                                                                                                 by = "ID"
                                                                                       ),
                                                                                     var_x = "diff.AnnTemp.LGM_cur",
                                                                                     var_y = "nti",
                                                                                     lab_x = c(expression(atop(
                                                                                       "Historical Change in Temperature (°C)",
                                                                                       scriptstyle("MAT"[Contemporary] - "MAT"[LGM])
                                                                                     ))),
                                                                                     lab_y = c(expression("NTI"["Plate"]))
)


fig.kde.Plate.NTI.diff.log.AnnPrec.log.LGM_cur.Hpi <- ks_kde_geom_contours_phyloStr_diffClim(parent_data = MPD.MNTD.Chiroptera.Comm.AllScales.CWM.Div %>%
                                                                                               filter(SamplingPool == "Plate sampling") %>%
                                                                                               left_join(MPD.MNTD.LatLong.AllScales.worldClimate.diff %>%
                                                                                                           filter(SamplingPool == "Plate sampling") %>%
                                                                                                           select(
                                                                                                             ID,
                                                                                                             diff.AnnTemp.LGM_cur:diff.log.AnnPrec.log.LGM_cur
                                                                                                           ),
                                                                                                         by = "ID"
                                                                                               ),
                                                                                             var_x = "diff.log.AnnPrec.log.LGM_cur",
                                                                                             var_y = "nti",
                                                                                             lab_x = c(expression(atop(
                                                                                               "Historical Change in Precipitation",
                                                                                               scriptstyle(log("MAP"[Contemporary]) - log("MAP"[LGM]))
                                                                                             ))),
                                                                                             lab_y = c(expression("NTI"["Plate"]))
)
# Biome ------------

## NRI and Annual Temperature and Precipitation ------------

fig.kde.Biome.NRI.diff.AnnTemp.LGM_cur.Hpi <- ks_kde_geom_contours_phyloStr_diffClim(parent_data = MPD.MNTD.Chiroptera.Comm.AllScales.CWM.Div %>%
                                                                                       filter(SamplingPool == "Biome sampling") %>%
                                                                                       left_join(MPD.MNTD.LatLong.AllScales.worldClimate.diff %>%
                                                                                                   filter(SamplingPool == "Biome sampling") %>%
                                                                                                   select(
                                                                                                     ID,
                                                                                                     diff.AnnTemp.LGM_cur:diff.log.AnnPrec.log.LGM_cur
                                                                                                   ),
                                                                                                 by = "ID"
                                                                                       ),
                                                                                     var_x = "diff.AnnTemp.LGM_cur",
                                                                                     var_y = "nri",
                                                                                     lab_x = c(expression(atop(
                                                                                       "Historical Change in Temperature (°C)",
                                                                                       scriptstyle("MAT"[Contemporary] - "MAT"[LGM])
                                                                                     ))),
                                                                                     lab_y = c(expression("NRI"["Biome"]))
)


fig.kde.Biome.NRI.diff.log.AnnPrec.log.LGM_cur.Hpi <- ks_kde_geom_contours_phyloStr_diffClim(parent_data = MPD.MNTD.Chiroptera.Comm.AllScales.CWM.Div %>%
                                                                                               filter(SamplingPool == "Biome sampling") %>%
                                                                                               left_join(MPD.MNTD.LatLong.AllScales.worldClimate.diff %>%
                                                                                                           filter(SamplingPool == "Biome sampling") %>%
                                                                                                           select(
                                                                                                             ID,
                                                                                                             diff.AnnTemp.LGM_cur:diff.log.AnnPrec.log.LGM_cur
                                                                                                           ),
                                                                                                         by = "ID"
                                                                                               ),
                                                                                             var_x = "diff.log.AnnPrec.log.LGM_cur",
                                                                                             var_y = "nri",
                                                                                             lab_x = c(expression(atop(
                                                                                               "Historical Change in Precipitation",
                                                                                               scriptstyle(log("MAP"[Contemporary]) - log("MAP"[LGM]))
                                                                                             ))),
                                                                                             lab_y = c(expression("NRI"["Biome"]))
)


## NTI and Annual Temperature and Precipitation ------------


fig.kde.Biome.NTI.diff.AnnTemp.LGM_cur.Hpi <- ks_kde_geom_contours_phyloStr_diffClim(parent_data = MPD.MNTD.Chiroptera.Comm.AllScales.CWM.Div %>%
                                                                                       filter(SamplingPool == "Biome sampling") %>%
                                                                                       left_join(MPD.MNTD.LatLong.AllScales.worldClimate.diff %>%
                                                                                                   filter(SamplingPool == "Biome sampling") %>%
                                                                                                   select(
                                                                                                     ID,
                                                                                                     diff.AnnTemp.LGM_cur:diff.log.AnnPrec.log.LGM_cur
                                                                                                   ),
                                                                                                 by = "ID"
                                                                                       ),
                                                                                     var_x = "diff.AnnTemp.LGM_cur",
                                                                                     var_y = "nti",
                                                                                     lab_x = c(expression(atop(
                                                                                       "Historical Change in Temperature (°C)",
                                                                                       scriptstyle("MAT"[Contemporary] - "MAT"[LGM])
                                                                                     ))),
                                                                                     lab_y = c(expression("NTI"["Biome"]))
)


fig.kde.Biome.NTI.diff.log.AnnPrec.log.LGM_cur.Hpi <- ks_kde_geom_contours_phyloStr_diffClim(parent_data = MPD.MNTD.Chiroptera.Comm.AllScales.CWM.Div %>%
                                                                                               filter(SamplingPool == "Biome sampling") %>%
                                                                                               left_join(MPD.MNTD.LatLong.AllScales.worldClimate.diff %>%
                                                                                                           filter(SamplingPool == "Biome sampling") %>%
                                                                                                           select(
                                                                                                             ID,
                                                                                                             diff.AnnTemp.LGM_cur:diff.log.AnnPrec.log.LGM_cur
                                                                                                           ),
                                                                                                         by = "ID"
                                                                                               ),
                                                                                             var_x = "diff.log.AnnPrec.log.LGM_cur",
                                                                                             var_y = "nti",
                                                                                             lab_x = c(expression(atop(
                                                                                               "Historical Change in Precipitation",
                                                                                               scriptstyle(log("MAP"[Contemporary]) - log("MAP"[LGM]))
                                                                                             ))),
                                                                                             lab_y = c(expression("NTI"["Biome"]))
)

# Ecoregion ---------

## NRI and Annual Temperature and Precipitation -----------

fig.kde.Ecoregion.NRI.diff.AnnTemp.LGM_cur.Hpi <- ks_kde_geom_contours_phyloStr_diffClim(parent_data = MPD.MNTD.Chiroptera.Comm.AllScales.CWM.Div %>%
                                                                                           filter(SamplingPool == "Ecoregion sampling") %>%
                                                                                           left_join(MPD.MNTD.LatLong.AllScales.worldClimate.diff %>%
                                                                                                       filter(SamplingPool == "Ecoregion sampling") %>%
                                                                                                       select(
                                                                                                         ID,
                                                                                                         diff.AnnTemp.LGM_cur:diff.log.AnnPrec.log.LGM_cur
                                                                                                       ),
                                                                                                     by = "ID"
                                                                                           ),
                                                                                         var_x = "diff.AnnTemp.LGM_cur",
                                                                                         var_y = "nri",
                                                                                         lab_x = c(expression(atop(
                                                                                           "Historical Change in Temperature (°C)",
                                                                                           scriptstyle("MAT"[Contemporary] - "MAT"[LGM])
                                                                                         ))),
                                                                                         lab_y = c(expression("NRI"["Ecoregion"]))
)


fig.kde.Ecoregion.NRI.diff.log.AnnPrec.log.LGM_cur.Hpi <- ks_kde_geom_contours_phyloStr_diffClim(parent_data = MPD.MNTD.Chiroptera.Comm.AllScales.CWM.Div %>%
                                                                                                   filter(SamplingPool == "Ecoregion sampling") %>%
                                                                                                   left_join(MPD.MNTD.LatLong.AllScales.worldClimate.diff %>%
                                                                                                               filter(SamplingPool == "Ecoregion sampling") %>%
                                                                                                               select(
                                                                                                                 ID,
                                                                                                                 diff.AnnTemp.LGM_cur:diff.log.AnnPrec.log.LGM_cur
                                                                                                               ),
                                                                                                             by = "ID"
                                                                                                   ),
                                                                                                 var_x = "diff.log.AnnPrec.log.LGM_cur",
                                                                                                 var_y = "nri",
                                                                                                 lab_x = c(expression(atop(
                                                                                                   "Historical Change in Precipitation",
                                                                                                   scriptstyle(log("MAP"[Contemporary]) - log("MAP"[LGM]))
                                                                                                 ))),
                                                                                                 lab_y = c(expression("NRI"["Ecoregion"]))
)


## NTI and Annual Temperature and Precipitation  --------


fig.kde.Ecoregion.NTI.diff.AnnTemp.LGM_cur.Hpi <- ks_kde_geom_contours_phyloStr_diffClim(parent_data = MPD.MNTD.Chiroptera.Comm.AllScales.CWM.Div %>%
                                                                                           filter(SamplingPool == "Ecoregion sampling") %>%
                                                                                           left_join(MPD.MNTD.LatLong.AllScales.worldClimate.diff %>%
                                                                                                       filter(SamplingPool == "Ecoregion sampling") %>%
                                                                                                       select(
                                                                                                         ID,
                                                                                                         diff.AnnTemp.LGM_cur:diff.log.AnnPrec.log.LGM_cur
                                                                                                       ),
                                                                                                     by = "ID"
                                                                                           ),
                                                                                         var_x = "diff.AnnTemp.LGM_cur",
                                                                                         var_y = "nti",
                                                                                         lab_x = c(expression(atop(
                                                                                           "Historical Change in Temperature (°C)",
                                                                                           scriptstyle("MAT"[Contemporary] - "MAT"[LGM])
                                                                                         ))),
                                                                                         lab_y = c(expression("NTI"["Ecoregion"]))
)


fig.kde.Ecoregion.NTI.diff.log.AnnPrec.log.LGM_cur.Hpi <- ks_kde_geom_contours_phyloStr_diffClim(parent_data = MPD.MNTD.Chiroptera.Comm.AllScales.CWM.Div %>%
                                                                                                   filter(SamplingPool == "Ecoregion sampling") %>%
                                                                                                   left_join(MPD.MNTD.LatLong.AllScales.worldClimate.diff %>%
                                                                                                               filter(SamplingPool == "Ecoregion sampling") %>%
                                                                                                               select(
                                                                                                                 ID,
                                                                                                                 diff.AnnTemp.LGM_cur:diff.log.AnnPrec.log.LGM_cur
                                                                                                               ),
                                                                                                             by = "ID"
                                                                                                   ),
                                                                                                 var_x = "diff.log.AnnPrec.log.LGM_cur",
                                                                                                 var_y = "nti",
                                                                                                 lab_x = c(expression(atop(
                                                                                                   "Historical Change in Precipitation",
                                                                                                   scriptstyle(log("MAP"[Contemporary]) - log("MAP"[LGM]))
                                                                                                 ))),
                                                                                                 lab_y = c(expression("NTI"["Ecoregion"]))
)


# Combining plots --------------------------------------------

(fig.kde.SamplingPool.NRI.diffTemp.diffPrec <- ggarrange(
  fig.kde.Global.NRI.diff.AnnTemp.LGM_cur.Hpi +
    theme(axis.title.x = element_blank()),
  fig.kde.Global.NRI.diff.log.AnnPrec.log.LGM_cur.Hpi +
    theme(axis.title.x = element_blank()),
  fig.kde.Hemispheric.NRI.diff.AnnTemp.LGM_cur.Hpi +
    theme(axis.title.x = element_blank()),
  fig.kde.Hemispheric.NRI.diff.log.AnnPrec.log.LGM_cur.Hpi +
    theme(axis.title.x = element_blank()),
  fig.kde.Realm.NRI.diff.AnnTemp.LGM_cur.Hpi +
    theme(axis.title.x = element_blank()),
  fig.kde.Realm.NRI.diff.log.AnnPrec.log.LGM_cur.Hpi +
    theme(axis.title.x = element_blank()),
  fig.kde.Plate.NRI.diff.AnnTemp.LGM_cur.Hpi +
    theme(axis.title.x = element_blank()),
  fig.kde.Plate.NRI.diff.log.AnnPrec.log.LGM_cur.Hpi +
    theme(axis.title.x = element_blank()),
  fig.kde.Biome.NRI.diff.AnnTemp.LGM_cur.Hpi +
    theme(axis.title.x = element_blank()),
  fig.kde.Biome.NRI.diff.log.AnnPrec.log.LGM_cur.Hpi +
    theme(axis.title.x = element_blank()),  
  fig.kde.Ecoregion.NRI.diff.AnnTemp.LGM_cur.Hpi,
  fig.kde.Ecoregion.NRI.diff.log.AnnPrec.log.LGM_cur.Hpi,
  labels = c("A", "B", 
             "C", "D",
             "E", "F",
             "G", "H",
             "I", "J",
             "K", "L"
  ),
  ncol = 2, nrow = 6,
  widths = c(1, 1),
  heights = c(0.9, 1),
  align = "v",
  font.label = list(size = 28)
)
)

ggsave(
  file = "figures/fig.kde.SamplingPool.NRI.diffTemp.diffPrec.png",
  fig.kde.SamplingPool.NRI.diffTemp.diffPrec,
  width = 16.5,
  height = 7.5*6,
  units = c("in"),
  dpi = 250,
  limitsize = FALSE
)


(fig.kde.SamplingPool.NTI.diffTemp.diffPrec <- ggarrange(
  fig.kde.Global.NTI.diff.AnnTemp.LGM_cur.Hpi +
    theme(axis.title.x = element_blank()),
  fig.kde.Global.NTI.diff.log.AnnPrec.log.LGM_cur.Hpi +
    theme(axis.title.x = element_blank()),
  fig.kde.Hemispheric.NTI.diff.AnnTemp.LGM_cur.Hpi +
    theme(axis.title.x = element_blank()),
  fig.kde.Hemispheric.NTI.diff.log.AnnPrec.log.LGM_cur.Hpi +
    theme(axis.title.x = element_blank()),
  fig.kde.Realm.NTI.diff.AnnTemp.LGM_cur.Hpi +
    theme(axis.title.x = element_blank()),
  fig.kde.Realm.NTI.diff.log.AnnPrec.log.LGM_cur.Hpi +
    theme(axis.title.x = element_blank()),
  fig.kde.Plate.NTI.diff.AnnTemp.LGM_cur.Hpi +
    theme(axis.title.x = element_blank()),
  fig.kde.Plate.NTI.diff.log.AnnPrec.log.LGM_cur.Hpi +
    theme(axis.title.x = element_blank()),
  fig.kde.Biome.NTI.diff.AnnTemp.LGM_cur.Hpi +
    theme(axis.title.x = element_blank()),
  fig.kde.Biome.NTI.diff.log.AnnPrec.log.LGM_cur.Hpi +
    theme(axis.title.x = element_blank()),  
  fig.kde.Ecoregion.NTI.diff.AnnTemp.LGM_cur.Hpi,
  fig.kde.Ecoregion.NTI.diff.log.AnnPrec.log.LGM_cur.Hpi,
  labels = c("A", "B", 
             "C", "D",
             "E", "F",
             "G", "H",
             "I", "J",
             "K", "L"
  ),
  ncol = 2, nrow = 6,
  widths = c(1, 1),
  heights = c(0.9, 1),
  align = "v",
  font.label = list(size = 28)
)
)

ggsave(
  file = "figures/fig.kde.SamplingPool.NTI.diffTemp.diffPrec.png",
  fig.kde.SamplingPool.NTI.diffTemp.diffPrec,
  width = 16.5,
  height = 7.5*6,
  units = c("in"),
  dpi = 250,
  limitsize = FALSE
)


