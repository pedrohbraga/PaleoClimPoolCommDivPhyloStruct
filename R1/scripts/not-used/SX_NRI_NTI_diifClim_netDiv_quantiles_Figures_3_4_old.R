#### NRI ####

(fig.Global.mean.NRI.diff.AnnTemp.LGM_cur.quantile <- quantile.XY(
  n.classes = 100,
  Y = "nri",
  X = "diff.AnnTemp.LGM_cur",
  dataset = MPD.MNTD.LatLong.AllScales.raref.rel.worldClimate.diff.CWM.Div %>%
    filter(SamplingPool == "Global sampling") %>%
    drop_na("nri") %>%
    drop_na("diff.AnnTemp.LGM_cur"),
  lab_x = c(expression(atop(
    ~italic(Q)*"(Historical Change in Temperature) (°C)",
    scriptstyle("MAT"[Contemporary] - "MAT"[LGM])
  ))),
  lab_y = c(expression(
    bar("NRI")["Global"])
  )
)
)

(fig.Hemispheric.mean.NRI.diff.AnnTemp.LGM_cur.quantile <- quantile.XY(
  n.classes = 100,
  Y = "nri",
  X = "diff.AnnTemp.LGM_cur",
  dataset = MPD.MNTD.LatLong.AllScales.raref.rel.worldClimate.diff.CWM.Div %>%
    filter(SamplingPool == "Hemispheric sampling") %>%
    drop_na("nri") %>%
    drop_na("diff.AnnTemp.LGM_cur"),
  lab_x = c(expression(atop(
    ~italic(Q)*"(Historical Change in Temperature) (°C)",
    scriptstyle("MAT"[Contemporary] - "MAT"[LGM])
  ))),
  lab_y = c(expression(bar("NRI")["Hemispheric"]))
)
)

(fig.Realm.mean.NRI.diff.AnnTemp.LGM_cur.quantile <- quantile.XY(
  n.classes = 100,
  Y = "nri",
  X = "diff.AnnTemp.LGM_cur",
  dataset = MPD.MNTD.LatLong.AllScales.raref.rel.worldClimate.diff.CWM.Div %>%
    filter(SamplingPool == "Realm sampling") %>%
    drop_na("nri") %>%
    drop_na("diff.AnnTemp.LGM_cur"),
  lab_x = c(expression(atop(
    ~italic(Q)*"(Historical Change in Temperature) (°C)",
    scriptstyle("MAT"[Contemporary] - "MAT"[LGM])
  ))),
  lab_y = c(expression(bar("NRI")["Realm"]))
)
)

(fig.Plate.mean.NRI.diff.AnnTemp.LGM_cur.quantile <- quantile.XY(
  n.classes = 100,
  Y = "nri",
  X = "diff.AnnTemp.LGM_cur",
  dataset = MPD.MNTD.LatLong.AllScales.raref.rel.worldClimate.diff.CWM.Div %>%
    filter(SamplingPool == "Plate sampling") %>%
    drop_na("nri") %>%
    drop_na("diff.AnnTemp.LGM_cur"),
  lab_x = c(expression(atop(
    ~italic(Q)*"(Historical Change in Temperature) (°C)",
    scriptstyle("MAT"[Contemporary] - "MAT"[LGM])
  ))),
  lab_y = c(expression(bar("NRI")["Plate"]))
)
)

(fig.Biome.mean.NRI.diff.AnnTemp.LGM_cur.quantile <- quantile.XY(
  n.classes = 100,
  Y = "nri",
  X = "diff.AnnTemp.LGM_cur",
  dataset = MPD.MNTD.LatLong.AllScales.raref.rel.worldClimate.diff.CWM.Div %>%
    filter(SamplingPool == "Biome sampling") %>%
    drop_na("nri") %>%
    drop_na("diff.AnnTemp.LGM_cur"),
  lab_x = c(expression(atop(
    ~italic(Q)*"(Historical Change in Temperature) (°C)",
    scriptstyle("MAT"[Contemporary] - "MAT"[LGM])
  ))),
  lab_y = c(expression(bar("NRI")["Biome"]))
)
)

(fig.Ecoregion.mean.NRI.diff.AnnTemp.LGM_cur.quantile <- quantile.XY(
  n.classes = 100,
  Y = "nri",
  X = "diff.AnnTemp.LGM_cur",
  dataset = MPD.MNTD.LatLong.AllScales.raref.rel.worldClimate.diff.CWM.Div %>%
    filter(SamplingPool == "Ecoregion sampling") %>%
    drop_na("nri") %>%
    drop_na("diff.AnnTemp.LGM_cur"),
  lab_x = c(expression(atop(
    ~italic(Q)*"(Historical Change in Temperature) °C",
    scriptstyle("MAT"[Contemporary] - "MAT"[LGM])
  ))),
  lab_y = c(expression(bar("NRI")["Ecoregion"]))
)
)


###

(fig.Global.mean.NRI.diff.log.AnnPrec.log.LGM_cur.quantile <- quantile.XY(
  n.classes = 100,
  Y = "nri",
  X = "diff.log.AnnPrec.log.LGM_cur",
  dataset = MPD.MNTD.LatLong.AllScales.raref.rel.worldClimate.diff.CWM.Div %>%
    filter(SamplingPool == "Global sampling") %>%
    drop_na("nri") %>%
    drop_na("diff.log.AnnPrec.log.LGM_cur"),
  lab_x =  c(expression(atop(
    ~italic(Q)*"(Historical Change in Precipitation)",
    scriptstyle(log("MAP"[Contemporary]) - log("MAP"[LGM]))
  ))),
  lab_y = c(expression(bar("NRI")["Global"]))
)
)

(fig.Hemispheric.mean.NRI.diff.log.AnnPrec.log.LGM_cur.quantile <- quantile.XY(
  n.classes = 100,
  Y = "nri",
  X = "diff.log.AnnPrec.log.LGM_cur",
  dataset = MPD.MNTD.LatLong.AllScales.raref.rel.worldClimate.diff.CWM.Div %>%
    filter(SamplingPool == "Hemispheric sampling") %>%
    drop_na("nri") %>%
    drop_na("diff.log.AnnPrec.log.LGM_cur"),
  lab_x =  c(expression(atop(
    ~italic(Q)*"(Historical Change in Precipitation)",
    scriptstyle(log("MAP"[Contemporary]) - log("MAP"[LGM]))
  ))),
  lab_y = c(expression(bar("NRI")["Hemispheric"]))
)
)

(fig.Realm.mean.NRI.diff.log.AnnPrec.log.LGM_cur.quantile <- quantile.XY(
  n.classes = 100,
  Y = "nri",
  X = "diff.log.AnnPrec.log.LGM_cur",
  dataset = MPD.MNTD.LatLong.AllScales.raref.rel.worldClimate.diff.CWM.Div %>%
    filter(SamplingPool == "Realm sampling") %>%
    drop_na("nri") %>%
    drop_na("diff.log.AnnPrec.log.LGM_cur"),
  lab_x =  c(expression(atop(
    ~italic(Q)*"(Historical Change in Precipitation)",
    scriptstyle(log("MAP"[Contemporary]) - log("MAP"[LGM]))
  ))),
  lab_y = c(expression(bar("NRI")["Realm"]))
)
)

(fig.Plate.mean.NRI.diff.log.AnnPrec.log.LGM_cur.quantile <- quantile.XY(
  n.classes = 100,
  Y = "nri",
  X = "diff.log.AnnPrec.log.LGM_cur",
  dataset = MPD.MNTD.LatLong.AllScales.raref.rel.worldClimate.diff.CWM.Div %>%
    filter(SamplingPool == "Plate sampling") %>%
    drop_na("nri") %>%
    drop_na("diff.log.AnnPrec.log.LGM_cur"),
  lab_x =  c(expression(atop(
    ~italic(Q)*"(Historical Change in Precipitation)",
    scriptstyle(log("MAP"[Contemporary]) - log("MAP"[LGM]))
  ))),
  lab_y = c(expression(bar("NRI")["Plate"]))
)
)

(fig.Biome.mean.NRI.diff.log.AnnPrec.log.LGM_cur.quantile <- quantile.XY(
  n.classes = 100,
  Y = "nri",
  X = "diff.log.AnnPrec.log.LGM_cur",
  dataset = MPD.MNTD.LatLong.AllScales.raref.rel.worldClimate.diff.CWM.Div %>%
    filter(SamplingPool == "Biome sampling") %>%
    drop_na("nri") %>%
    drop_na("diff.log.AnnPrec.log.LGM_cur"),
  lab_x =  c(expression(atop(
    ~italic(Q)*"(Historical Change in Precipitation)",
    scriptstyle(log("MAP"[Contemporary]) - log("MAP"[LGM]))
  ))),
  lab_y = c(expression(bar("NRI")["Biome"]))
)
)

(fig.Ecoregion.mean.NRI.diff.log.AnnPrec.log.LGM_cur.quantile <- quantile.XY(
  n.classes = 100,
  Y = "nri",
  X = "diff.log.AnnPrec.log.LGM_cur",
  dataset = MPD.MNTD.LatLong.AllScales.raref.rel.worldClimate.diff.CWM.Div %>%
    filter(SamplingPool == "Ecoregion sampling") %>%
    drop_na("nri") %>%
    drop_na("diff.log.AnnPrec.log.LGM_cur"),
  lab_x =  c(expression(atop(
    ~italic(Q)*"(Historical Change in Precipitation)",
    scriptstyle(log("MAP"[Contemporary]) - log("MAP"[LGM]))
  ))),
  lab_y = c(expression(bar("NRI")["Ecoregion"]))
)
)


# Combining plots --------------------------------------------

# (fig.quantile.SamplingPool.NRI.diffTemp.diffPrec <- ggarrange(
#   fig.Global.mean.NRI.diff.AnnTemp.LGM_cur.quantile +
#     theme(axis.title.x = element_blank()),
#   fig.Global.mean.NRI.diff.log.AnnPrec.log.LGM_cur.quantile +
#     theme(axis.title.x = element_blank()),
#   fig.Hemispheric.mean.NRI.diff.AnnTemp.LGM_cur.quantile +
#     theme(axis.title.x = element_blank()),
#   fig.Hemispheric.mean.NRI.diff.log.AnnPrec.log.LGM_cur.quantile +
#     theme(axis.title.x = element_blank()),
#   fig.Realm.mean.NRI.diff.AnnTemp.LGM_cur.quantile +
#     theme(axis.title.x = element_blank()),
#   fig.Realm.mean.NRI.diff.log.AnnPrec.log.LGM_cur.quantile +
#     theme(axis.title.x = element_blank()),
#   fig.Plate.mean.NRI.diff.AnnTemp.LGM_cur.quantile +
#     theme(axis.title.x = element_blank()),
#   fig.Plate.mean.NRI.diff.log.AnnPrec.log.LGM_cur.quantile +
#     theme(axis.title.x = element_blank()),
#   fig.Biome.mean.NRI.diff.AnnTemp.LGM_cur.quantile +
#     theme(axis.title.x = element_blank()),
#   fig.Biome.mean.NRI.diff.log.AnnPrec.log.LGM_cur.quantile +
#     theme(axis.title.x = element_blank()),  
#   fig.Ecoregion.mean.NRI.diff.AnnTemp.LGM_cur.quantile,
#   fig.Ecoregion.mean.NRI.diff.log.AnnPrec.log.LGM_cur.quantile,
#   labels = c("A", "B", 
#              "C", "D",
#              "E", "F",
#              "G", "H",
#              "I", "J",
#              "K", "L"
#   ),
#   ncol = 2, nrow = 6,
#   widths = c(1, 1),
#   heights = c(0.9, 1),
#   align = "v",
#   font.label = list(size = 28)
# )
# )

# ggsave(
#   file = "figures/fig.quantile.SamplingPool.NRI.diffTemp.diffPrec.png",
#   fig.quantile.SamplingPool.NRI.diffTemp.diffPrec,
#   width = 16.5,
#   height = 7.5*6,
#   units = c("in"),
#   dpi = 250,
#   limitsize = FALSE
# )
# 

###


(fig.Global.mean.NRI.netDiv_CWM_std_tw_rs_1.quantile <- quantile.XY(
  n.classes = 100,
  Y = "nri",
  X = "netDiv_CWM_std_tw_rs_1",
  dataset = MPD.MNTD.LatLong.AllScales.raref.rel.worldClimate.diff.CWM.Div %>%
    filter(SamplingPool == "Global sampling") %>%
    drop_na("nri") %>%
    drop_na("netDiv_CWM_std_tw_rs_1"),
  lab_x = c(expression(atop(
    ~italic(Q)*"(Local Net Diversification Rate)",
    lambda[CWM[tw]] - mu[CWM[tw]]))),
  lab_y = c(expression(bar("NRI")["Global"]))
)
)

(fig.Hemispheric.mean.NRI.netDiv_CWM_std_tw_rs_1.quantile <- quantile.XY(
  n.classes = 100,
  Y = "nri",
  X = "netDiv_CWM_std_tw_rs_1",
  dataset = MPD.MNTD.LatLong.AllScales.raref.rel.worldClimate.diff.CWM.Div %>%
    filter(SamplingPool == "Hemispheric sampling") %>%
    drop_na("nri") %>%
    drop_na("netDiv_CWM_std_tw_rs_1"),
  lab_x = c(expression(
    atop(
    ~italic(Q)*"(Local Net Diversification Rate)",
    lambda[CWM[tw]] - mu[CWM[tw]]))),
  lab_y = c(expression(bar("NRI")["Hemispheric"]))
)
)

(fig.Realm.mean.NRI.netDiv_CWM_std_tw_rs_1.quantile <- quantile.XY(
  n.classes = 100,
  Y = "nri",
  X = "netDiv_CWM_std_tw_rs_1",
  dataset = MPD.MNTD.LatLong.AllScales.raref.rel.worldClimate.diff.CWM.Div %>%
    filter(SamplingPool == "Realm sampling") %>%
    drop_na("nri") %>%
    drop_na("netDiv_CWM_std_tw_rs_1"),
  lab_x = c(expression(atop(
    ~italic(Q)*"(Local Net Diversification Rate)",
    lambda[CWM[tw]] - mu[CWM[tw]]))),
  lab_y = c(expression(bar("NRI")["Realm"]))
)
)

(fig.Plate.mean.NRI.netDiv_CWM_std_tw_rs_1.quantile <- quantile.XY(
  n.classes = 100,
  Y = "nri",
  X = "netDiv_CWM_std_tw_rs_1",
  dataset = MPD.MNTD.LatLong.AllScales.raref.rel.worldClimate.diff.CWM.Div %>%
    filter(SamplingPool == "Plate sampling") %>%
    drop_na("nri") %>%
    drop_na("netDiv_CWM_std_tw_rs_1"),
  lab_x = c(expression(atop(
    ~italic(Q)*"(Local Net Diversification Rate)",
    lambda[CWM[tw]] - mu[CWM[tw]]))),
  lab_y = c(expression(bar("NRI")["Plate"]))
)
)

(fig.Biome.mean.NRI.netDiv_CWM_std_tw_rs_1.quantile <- quantile.XY(
  n.classes = 100,
  Y = "nri",
  X = "netDiv_CWM_std_tw_rs_1",
  dataset = MPD.MNTD.LatLong.AllScales.raref.rel.worldClimate.diff.CWM.Div %>%
    filter(SamplingPool == "Biome sampling") %>%
    drop_na("nri") %>%
    drop_na("netDiv_CWM_std_tw_rs_1"),
  lab_x = c(expression(atop(
    ~italic(Q)*"(Local Net Diversification Rate)",
    lambda[CWM[tw]] - mu[CWM[tw]]))),
  lab_y = c(expression(bar("NRI")["Biome"]))
)
)

(fig.Ecoregion.mean.NRI.netDiv_CWM_std_tw_rs_1.quantile <- quantile.XY(
  n.classes = 100,
  Y = "nri",
  X = "netDiv_CWM_std_tw_rs_1",
  dataset = MPD.MNTD.LatLong.AllScales.raref.rel.worldClimate.diff.CWM.Div %>%
    filter(SamplingPool == "Ecoregion sampling") %>%
    drop_na("nri") %>%
    drop_na("netDiv_CWM_std_tw_rs_1"),
  lab_x = c(expression(atop(
    ~italic(Q)*"(Local Net Diversification Rate)",
    lambda[CWM[tw]] - mu[CWM[tw]]))),
  lab_y = c(expression(bar("NRI")["Ecoregion"]))
)
)

# Combining plots --------------------------------------------

(fig.quantile.SamplingPool.NRI.diffTemp.diffPrec.netDiv <- ggpubr::ggarrange(
  fig.Global.mean.NRI.diff.AnnTemp.LGM_cur.quantile +
    theme(axis.title.x = element_blank()),
  fig.Global.mean.NRI.diff.log.AnnPrec.log.LGM_cur.quantile +
    theme(axis.title.x = element_blank()),
  fig.Global.mean.NRI.netDiv_CWM_std_tw_rs_1.quantile +
    theme(axis.title.x = element_blank()),
  fig.Hemispheric.mean.NRI.diff.AnnTemp.LGM_cur.quantile +
    theme(axis.title.x = element_blank()),
  fig.Hemispheric.mean.NRI.diff.log.AnnPrec.log.LGM_cur.quantile +
    theme(axis.title.x = element_blank()),
  fig.Hemispheric.mean.NRI.netDiv_CWM_std_tw_rs_1.quantile +
    theme(axis.title.x = element_blank()),
  fig.Realm.mean.NRI.diff.AnnTemp.LGM_cur.quantile +
    theme(axis.title.x = element_blank()),
  fig.Realm.mean.NRI.diff.log.AnnPrec.log.LGM_cur.quantile +
    theme(axis.title.x = element_blank()),
  fig.Realm.mean.NRI.netDiv_CWM_std_tw_rs_1.quantile +
    theme(axis.title.x = element_blank()),
  fig.Plate.mean.NRI.diff.AnnTemp.LGM_cur.quantile +
    theme(axis.title.x = element_blank()),
  fig.Plate.mean.NRI.diff.log.AnnPrec.log.LGM_cur.quantile +
    theme(axis.title.x = element_blank()),
  fig.Plate.mean.NRI.netDiv_CWM_std_tw_rs_1.quantile +
    theme(axis.title.x = element_blank()),
  fig.Biome.mean.NRI.diff.AnnTemp.LGM_cur.quantile +
    theme(axis.title.x = element_blank()),
  fig.Biome.mean.NRI.diff.log.AnnPrec.log.LGM_cur.quantile +
    theme(axis.title.x = element_blank()),
  fig.Biome.mean.NRI.netDiv_CWM_std_tw_rs_1.quantile +
    theme(axis.title.x = element_blank()),  
  fig.Ecoregion.mean.NRI.diff.AnnTemp.LGM_cur.quantile,
  fig.Ecoregion.mean.NRI.diff.log.AnnPrec.log.LGM_cur.quantile,
  fig.Ecoregion.mean.NRI.netDiv_CWM_std_tw_rs_1.quantile,
  labels = c("A", "B", "C", 
             "D", "E", "F",
             "G", "H", "I", 
             "J", "K", "L",
             "M", "N", "O",
             "P", "Q", "R"
  ),
  ncol = 3, nrow = 6,
  widths = c(1, 1),
  heights = c(0.9, 1),
  align = "v",
  font.label = list(size = 28)
)
)




ggsave(
  file = "figures/fig.quantile.SamplingPool.NRI.diffTemp.diffPrec.netDiv.svg",
  fig.quantile.SamplingPool.NRI.diffTemp.diffPrec.netDiv,
  width = 8*3.2,
  height = 7*6,
  units = c("in"),
  dpi = 400,
  limitsize = FALSE
)



#### NTI ####

(fig.Global.mean.NTI.diff.AnnTemp.LGM_cur.quantile <- quantile.XY(
  n.classes = 100,
  Y = "nti",
  X = "diff.AnnTemp.LGM_cur",
  dataset = MPD.MNTD.LatLong.AllScales.raref.rel.worldClimate.diff.CWM.Div %>%
    filter(SamplingPool == "Global sampling") %>%
    drop_na("nti") %>%
    drop_na("diff.AnnTemp.LGM_cur"),
  lab_x = c(expression(atop(
    ~italic(Q)*"(Historical Change in Temperature) (°C)",
    scriptstyle("MAT"[Contemporary] - "MAT"[LGM])
  ))),
  lab_y = c(expression(
    bar("NTI")["Global"])
  )
)
)

(fig.Hemispheric.mean.NTI.diff.AnnTemp.LGM_cur.quantile <- quantile.XY(
  n.classes = 100,
  Y = "nti",
  X = "diff.AnnTemp.LGM_cur",
  dataset = MPD.MNTD.LatLong.AllScales.raref.rel.worldClimate.diff.CWM.Div %>%
    filter(SamplingPool == "Hemispheric sampling") %>%
    drop_na("nti") %>%
    drop_na("diff.AnnTemp.LGM_cur"),
  lab_x = c(expression(atop(
    ~italic(Q)*"(Historical Change in Temperature) (°C)",
    scriptstyle("MAT"[Contemporary] - "MAT"[LGM])
  ))),
  lab_y = c(expression(bar("NTI")["Hemispheric"]))
)
)

(fig.Realm.mean.NTI.diff.AnnTemp.LGM_cur.quantile <- quantile.XY(
  n.classes = 100,
  Y = "nti",
  X = "diff.AnnTemp.LGM_cur",
  dataset = MPD.MNTD.LatLong.AllScales.raref.rel.worldClimate.diff.CWM.Div %>%
    filter(SamplingPool == "Realm sampling") %>%
    drop_na("nti") %>%
    drop_na("diff.AnnTemp.LGM_cur"),
  lab_x = c(expression(atop(
    ~italic(Q)*"(Historical Change in Temperature) (°C)",
    scriptstyle("MAT"[Contemporary] - "MAT"[LGM])
  ))),
  lab_y = c(expression(bar("NTI")["Realm"]))
)
)

(fig.Plate.mean.NTI.diff.AnnTemp.LGM_cur.quantile <- quantile.XY(
  n.classes = 100,
  Y = "nti",
  X = "diff.AnnTemp.LGM_cur",
  dataset = MPD.MNTD.LatLong.AllScales.raref.rel.worldClimate.diff.CWM.Div %>%
    filter(SamplingPool == "Plate sampling") %>%
    drop_na("nti") %>%
    drop_na("diff.AnnTemp.LGM_cur"),
  lab_x = c(expression(atop(
    ~italic(Q)*"(Historical Change in Temperature) (°C)",
    scriptstyle("MAT"[Contemporary] - "MAT"[LGM])
  ))),
  lab_y = c(expression(bar("NTI")["Plate"]))
)
)

(fig.Biome.mean.NTI.diff.AnnTemp.LGM_cur.quantile <- quantile.XY(
  n.classes = 100,
  Y = "nti",
  X = "diff.AnnTemp.LGM_cur",
  dataset = MPD.MNTD.LatLong.AllScales.raref.rel.worldClimate.diff.CWM.Div %>%
    filter(SamplingPool == "Biome sampling") %>%
    drop_na("nti") %>%
    drop_na("diff.AnnTemp.LGM_cur"),
  lab_x = c(expression(atop(
    ~italic(Q)*"(Historical Change in Temperature) (°C)",
    scriptstyle("MAT"[Contemporary] - "MAT"[LGM])
  ))),
  lab_y = c(expression(bar("NTI")["Biome"]))
)
)

(fig.Ecoregion.mean.NTI.diff.AnnTemp.LGM_cur.quantile <- quantile.XY(
  n.classes = 100,
  Y = "nti",
  X = "diff.AnnTemp.LGM_cur",
  dataset = MPD.MNTD.LatLong.AllScales.raref.rel.worldClimate.diff.CWM.Div %>%
    filter(SamplingPool == "Ecoregion sampling") %>%
    drop_na("nti") %>%
    drop_na("diff.AnnTemp.LGM_cur"),
  lab_x = c(expression(atop(
    ~italic(Q)*"(Historical Change in Temperature) °C",
    scriptstyle("MAT"[Contemporary] - "MAT"[LGM])
  ))),
  lab_y = c(expression(bar("NTI")["Ecoregion"]))
)
)


###

(fig.Global.mean.NTI.diff.log.AnnPrec.log.LGM_cur.quantile <- quantile.XY(
  n.classes = 100,
  Y = "nti",
  X = "diff.log.AnnPrec.log.LGM_cur",
  dataset = MPD.MNTD.LatLong.AllScales.raref.rel.worldClimate.diff.CWM.Div %>%
    filter(SamplingPool == "Global sampling") %>%
    drop_na("nti") %>%
    drop_na("diff.log.AnnPrec.log.LGM_cur"),
  lab_x =  c(expression(atop(
    ~italic(Q)*"(Historical Change in Precipitation)",
    scriptstyle(log("MAP"[Contemporary]) - log("MAP"[LGM]))
  ))),
  lab_y = c(expression(bar("NTI")["Global"]))
)
)

(fig.Hemispheric.mean.NTI.diff.log.AnnPrec.log.LGM_cur.quantile <- quantile.XY(
  n.classes = 100,
  Y = "nti",
  X = "diff.log.AnnPrec.log.LGM_cur",
  dataset = MPD.MNTD.LatLong.AllScales.raref.rel.worldClimate.diff.CWM.Div %>%
    filter(SamplingPool == "Hemispheric sampling") %>%
    drop_na("nti") %>%
    drop_na("diff.log.AnnPrec.log.LGM_cur"),
  lab_x =  c(expression(atop(
    ~italic(Q)*"(Historical Change in Precipitation)",
    scriptstyle(log("MAP"[Contemporary]) - log("MAP"[LGM]))
  ))),
  lab_y = c(expression(bar("NTI")["Hemispheric"]))
)
)

(fig.Realm.mean.NTI.diff.log.AnnPrec.log.LGM_cur.quantile <- quantile.XY(
  n.classes = 100,
  Y = "nti",
  X = "diff.log.AnnPrec.log.LGM_cur",
  dataset = MPD.MNTD.LatLong.AllScales.raref.rel.worldClimate.diff.CWM.Div %>%
    filter(SamplingPool == "Realm sampling") %>%
    drop_na("nti") %>%
    drop_na("diff.log.AnnPrec.log.LGM_cur"),
  lab_x =  c(expression(atop(
    ~italic(Q)*"(Historical Change in Precipitation)",
    scriptstyle(log("MAP"[Contemporary]) - log("MAP"[LGM]))
  ))),
  lab_y = c(expression(bar("NTI")["Realm"]))
)
)

(fig.Plate.mean.NTI.diff.log.AnnPrec.log.LGM_cur.quantile <- quantile.XY(
  n.classes = 100,
  Y = "nti",
  X = "diff.log.AnnPrec.log.LGM_cur",
  dataset = MPD.MNTD.LatLong.AllScales.raref.rel.worldClimate.diff.CWM.Div %>%
    filter(SamplingPool == "Plate sampling") %>%
    drop_na("nti") %>%
    drop_na("diff.log.AnnPrec.log.LGM_cur"),
  lab_x =  c(expression(atop(
    ~italic(Q)*"(Historical Change in Precipitation)",
    scriptstyle(log("MAP"[Contemporary]) - log("MAP"[LGM]))
  ))),
  lab_y = c(expression(bar("NTI")["Plate"]))
)
)

(fig.Biome.mean.NTI.diff.log.AnnPrec.log.LGM_cur.quantile <- quantile.XY(
  n.classes = 100,
  Y = "nti",
  X = "diff.log.AnnPrec.log.LGM_cur",
  dataset = MPD.MNTD.LatLong.AllScales.raref.rel.worldClimate.diff.CWM.Div %>%
    filter(SamplingPool == "Biome sampling") %>%
    drop_na("nti") %>%
    drop_na("diff.log.AnnPrec.log.LGM_cur"),
  lab_x =  c(expression(atop(
    ~italic(Q)*"(Historical Change in Precipitation)",
    scriptstyle(log("MAP"[Contemporary]) - log("MAP"[LGM]))
  ))),
  lab_y = c(expression(bar("NTI")["Biome"]))
)
)

(fig.Ecoregion.mean.NTI.diff.log.AnnPrec.log.LGM_cur.quantile <- quantile.XY(
  n.classes = 100,
  Y = "nti",
  X = "diff.log.AnnPrec.log.LGM_cur",
  dataset = MPD.MNTD.LatLong.AllScales.raref.rel.worldClimate.diff.CWM.Div %>%
    filter(SamplingPool == "Ecoregion sampling") %>%
    drop_na("nti") %>%
    drop_na("diff.log.AnnPrec.log.LGM_cur"),
  lab_x =  c(expression(atop(
    ~italic(Q)*"(Historical Change in Precipitation)",
    scriptstyle(log("MAP"[Contemporary]) - log("MAP"[LGM]))
  ))),
  lab_y = c(expression(bar("NTI")["Ecoregion"]))
)
)


# Combining plots --------------------------------------------

# (fig.quantile.SamplingPool.NTI.diffTemp.diffPrec <- ggarrange(
#   fig.Global.mean.NTI.diff.AnnTemp.LGM_cur.quantile +
#     theme(axis.title.x = element_blank()),
#   fig.Global.mean.NTI.diff.log.AnnPrec.log.LGM_cur.quantile +
#     theme(axis.title.x = element_blank()),
#   fig.Hemispheric.mean.NTI.diff.AnnTemp.LGM_cur.quantile +
#     theme(axis.title.x = element_blank()),
#   fig.Hemispheric.mean.NTI.diff.log.AnnPrec.log.LGM_cur.quantile +
#     theme(axis.title.x = element_blank()),
#   fig.Realm.mean.NTI.diff.AnnTemp.LGM_cur.quantile +
#     theme(axis.title.x = element_blank()),
#   fig.Realm.mean.NTI.diff.log.AnnPrec.log.LGM_cur.quantile +
#     theme(axis.title.x = element_blank()),
#   fig.Plate.mean.NTI.diff.AnnTemp.LGM_cur.quantile +
#     theme(axis.title.x = element_blank()),
#   fig.Plate.mean.NTI.diff.log.AnnPrec.log.LGM_cur.quantile +
#     theme(axis.title.x = element_blank()),
#   fig.Biome.mean.NTI.diff.AnnTemp.LGM_cur.quantile +
#     theme(axis.title.x = element_blank()),
#   fig.Biome.mean.NTI.diff.log.AnnPrec.log.LGM_cur.quantile +
#     theme(axis.title.x = element_blank()),  
#   fig.Ecoregion.mean.NTI.diff.AnnTemp.LGM_cur.quantile,
#   fig.Ecoregion.mean.NTI.diff.log.AnnPrec.log.LGM_cur.quantile,
#   labels = c("A", "B", 
#              "C", "D",
#              "E", "F",
#              "G", "H",
#              "I", "J",
#              "K", "L"
#   ),
#   ncol = 2, nrow = 6,
#   widths = c(1, 1),
#   heights = c(0.9, 1),
#   align = "v",
#   font.label = list(size = 28)
# )
# )

# ggsave(
#   file = "figures/fig.quantile.SamplingPool.NTI.diffTemp.diffPrec.png",
#   fig.quantile.SamplingPool.NTI.diffTemp.diffPrec,
#   width = 16.5,
#   height = 7.5*6,
#   units = c("in"),
#   dpi = 250,
#   limitsize = FALSE
# )
# 

###


(fig.Global.mean.NTI.netDiv_CWM_std_tw_rs_1.quantile <- quantile.XY(
  n.classes = 100,
  Y = "nti",
  X = "netDiv_CWM_std_tw_rs_1",
  dataset = MPD.MNTD.LatLong.AllScales.raref.rel.worldClimate.diff.CWM.Div %>%
    filter(SamplingPool == "Global sampling") %>%
    drop_na("nti") %>%
    drop_na("netDiv_CWM_std_tw_rs_1"),
  lab_x = c(expression(atop(
    ~italic(Q)*"(Local Net Diversification Rate)",
    lambda[CWM[tw]] - mu[CWM[tw]]))),
  lab_y = c(expression(bar("NTI")["Global"]))
)
)

(fig.Hemispheric.mean.NTI.netDiv_CWM_std_tw_rs_1.quantile <- quantile.XY(
  n.classes = 100,
  Y = "nti",
  X = "netDiv_CWM_std_tw_rs_1",
  dataset = MPD.MNTD.LatLong.AllScales.raref.rel.worldClimate.diff.CWM.Div %>%
    filter(SamplingPool == "Hemispheric sampling") %>%
    drop_na("nti") %>%
    drop_na("netDiv_CWM_std_tw_rs_1"),
  lab_x = c(expression(atop(
    ~italic(Q)*"(Local Net Diversification Rate)",
    lambda[CWM[tw]] - mu[CWM[tw]]))),
  lab_y = c(expression(bar("NTI")["Hemispheric"]))
)
)

(fig.Realm.mean.NTI.netDiv_CWM_std_tw_rs_1.quantile <- quantile.XY(
  n.classes = 100,
  Y = "nti",
  X = "netDiv_CWM_std_tw_rs_1",
  dataset = MPD.MNTD.LatLong.AllScales.raref.rel.worldClimate.diff.CWM.Div %>%
    filter(SamplingPool == "Realm sampling") %>%
    drop_na("nti") %>%
    drop_na("netDiv_CWM_std_tw_rs_1"),
  lab_x = c(expression(atop(
    ~italic(Q)*"(Local Net Diversification Rate)",
    lambda[CWM[tw]] - mu[CWM[tw]]))),
  lab_y = c(expression(bar("NTI")["Realm"]))
)
)

(fig.Plate.mean.NTI.netDiv_CWM_std_tw_rs_1.quantile <- quantile.XY(
  n.classes = 100,
  Y = "nti",
  X = "netDiv_CWM_std_tw_rs_1",
  dataset = MPD.MNTD.LatLong.AllScales.raref.rel.worldClimate.diff.CWM.Div %>%
    filter(SamplingPool == "Plate sampling") %>%
    drop_na("nti") %>%
    drop_na("netDiv_CWM_std_tw_rs_1"),
  lab_x = c(expression(atop(
    ~italic(Q)*"(Local Net Diversification Rate)",
    lambda[CWM[tw]] - mu[CWM[tw]]))),
  lab_y = c(expression(bar("NTI")["Plate"]))
)
)

(fig.Biome.mean.NTI.netDiv_CWM_std_tw_rs_1.quantile <- quantile.XY(
  n.classes = 100,
  Y = "nti",
  X = "netDiv_CWM_std_tw_rs_1",
  dataset = MPD.MNTD.LatLong.AllScales.raref.rel.worldClimate.diff.CWM.Div %>%
    filter(SamplingPool == "Biome sampling") %>%
    drop_na("nti") %>%
    drop_na("netDiv_CWM_std_tw_rs_1"),
  lab_x = c(expression(atop(
    ~italic(Q)*"(Local Net Diversification Rate)",
    lambda[CWM[tw]] - mu[CWM[tw]]))),
  lab_y = c(expression(bar("NTI")["Biome"]))
)
)

(fig.Ecoregion.mean.NTI.netDiv_CWM_std_tw_rs_1.quantile <- quantile.XY(
  n.classes = 100,
  Y = "nti",
  X = "netDiv_CWM_std_tw_rs_1",
  dataset = MPD.MNTD.LatLong.AllScales.raref.rel.worldClimate.diff.CWM.Div %>%
    filter(SamplingPool == "Ecoregion sampling") %>%
    drop_na("nti") %>%
    drop_na("netDiv_CWM_std_tw_rs_1"),
  lab_x = c(expression(atop(
    ~italic(Q)*"(Local Net Diversification Rate)",
    lambda[CWM[tw]] - mu[CWM[tw]]))),
  lab_y = c(expression(bar("NTI")["Ecoregion"]))
)
)

# Combining plots --------------------------------------------

(fig.quantile.SamplingPool.NTI.diffTemp.diffPrec.netDiv <- ggpubr::ggarrange(
  fig.Global.mean.NTI.diff.AnnTemp.LGM_cur.quantile +
    theme(axis.title.x = element_blank()),
  fig.Global.mean.NTI.diff.log.AnnPrec.log.LGM_cur.quantile +
    theme(axis.title.x = element_blank()),
  fig.Global.mean.NTI.netDiv_CWM_std_tw_rs_1.quantile +
    theme(axis.title.x = element_blank()),
  fig.Hemispheric.mean.NTI.diff.AnnTemp.LGM_cur.quantile +
    theme(axis.title.x = element_blank()),
  fig.Hemispheric.mean.NTI.diff.log.AnnPrec.log.LGM_cur.quantile +
    theme(axis.title.x = element_blank()),
  fig.Hemispheric.mean.NTI.netDiv_CWM_std_tw_rs_1.quantile +
    theme(axis.title.x = element_blank()),
  fig.Realm.mean.NTI.diff.AnnTemp.LGM_cur.quantile +
    theme(axis.title.x = element_blank()),
  fig.Realm.mean.NTI.diff.log.AnnPrec.log.LGM_cur.quantile +
    theme(axis.title.x = element_blank()),
  fig.Realm.mean.NTI.netDiv_CWM_std_tw_rs_1.quantile +
    theme(axis.title.x = element_blank()),
  fig.Plate.mean.NTI.diff.AnnTemp.LGM_cur.quantile +
    theme(axis.title.x = element_blank()),
  fig.Plate.mean.NTI.diff.log.AnnPrec.log.LGM_cur.quantile +
    theme(axis.title.x = element_blank()),
  fig.Plate.mean.NTI.netDiv_CWM_std_tw_rs_1.quantile +
    theme(axis.title.x = element_blank()),
  fig.Biome.mean.NTI.diff.AnnTemp.LGM_cur.quantile +
    theme(axis.title.x = element_blank()),
  fig.Biome.mean.NTI.diff.log.AnnPrec.log.LGM_cur.quantile +
    theme(axis.title.x = element_blank()),
  fig.Biome.mean.NTI.netDiv_CWM_std_tw_rs_1.quantile +
    theme(axis.title.x = element_blank()),  
  fig.Ecoregion.mean.NTI.diff.AnnTemp.LGM_cur.quantile,
  fig.Ecoregion.mean.NTI.diff.log.AnnPrec.log.LGM_cur.quantile,
  fig.Ecoregion.mean.NTI.netDiv_CWM_std_tw_rs_1.quantile,
  labels = c("A", "B", "C", 
             "D", "E", "F",
             "G", "H", "I", 
             "J", "K", "L",
             "M", "N", "O",
             "P", "Q", "R"
  ),
  ncol = 3, nrow = 6,
  widths = c(1, 1),
  heights = c(0.9, 1),
  align = "v",
  font.label = list(size = 28)
)
)



ggsave(
  file = "figures/fig.quantile.SamplingPool.NTI.diffTemp.diffPrec.netDiv.svg",
  fig.quantile.SamplingPool.NTI.diffTemp.diffPrec.netDiv,
  width = 8*3.2,
  height = 7*6,
  units = c("in"),
  dpi = 400,
  limitsize = FALSE
)


ragg::agg_png("figures/fig.quantile.SamplingPool.NTI.diffTemp.diffPrec.netDiv.ragg.png", 
              width = 8*3.2, 
              height = 7*6, 
              units = "in", 
              res = 400,
              scaling = 1)
fig.quantile.SamplingPool.NTI.diffTemp.diffPrec.netDiv
dev.off()

