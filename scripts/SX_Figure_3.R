NRI.diff.AnnTemp.LGM_cur.Global.plot
NRI.diff.AnnTemp.LGM_cur.Hemispheric.plot
NRI.diff.AnnTemp.LGM_cur.Realm.plot
NRI.diff.AnnTemp.LGM_cur.Plate.plot
NRI.diff.AnnTemp.LGM_cur.Biome.plot
NRI.diff.AnnTemp.LGM_cur.Ecoregion.plot


ggarrange(NRI.diff.AnnTemp.LGM_cur.Global.plot + theme(axis.title.x = element_blank()),
          NRI.diff.AnnTemp.LGM_cur.Hemispheric.plot + theme(axis.title.x = element_blank()),
          NRI.diff.AnnTemp.LGM_cur.Realm.plot + theme(axis.title.x = element_blank()),
          NRI.diff.AnnTemp.LGM_cur.Plate.plot + theme(axis.title.x = element_blank()),
          NRI.diff.AnnTemp.LGM_cur.Biome.plot + theme(axis.title.x = element_blank()),
          NRI.diff.AnnTemp.LGM_cur.Ecoregion.plot,
          NTI.diff.AnnTemp.LGM_cur.Global.plot + theme(axis.title.x = element_blank()),
          NTI.diff.AnnTemp.LGM_cur.Hemispheric.plot + theme(axis.title.x = element_blank()),
          NTI.diff.AnnTemp.LGM_cur.Realm.plot + theme(axis.title.x = element_blank()),
          NTI.diff.AnnTemp.LGM_cur.Plate.plot + theme(axis.title.x = element_blank()),
          NTI.diff.AnnTemp.LGM_cur.Biome.plot + theme(axis.title.x = element_blank()),
          NTI.diff.AnnTemp.LGM_cur.Ecoregion.plot,
          ncol = 2, nrow = 6,
          align = "hv",
          common.legend = TRUE)


ggarrange(NRI.diff.AnnTemp.LGM_cur.Global.plot + theme(axis.title.x = element_blank()),
          NTI.diff.AnnTemp.LGM_cur.Global.plot + theme(axis.title.x = element_blank()),
          NRI.diff.AnnTemp.LGM_cur.Hemispheric.plot + theme(axis.title.x = element_blank()),
          NTI.diff.AnnTemp.LGM_cur.Hemispheric.plot + theme(axis.title.x = element_blank()),
          NRI.diff.AnnTemp.LGM_cur.Realm.plot + theme(axis.title.x = element_blank()),
          NTI.diff.AnnTemp.LGM_cur.Realm.plot + theme(axis.title.x = element_blank()),
          NRI.diff.AnnTemp.LGM_cur.Plate.plot + theme(axis.title.x = element_blank()),
          NTI.diff.AnnTemp.LGM_cur.Plate.plot + theme(axis.title.x = element_blank()),
          NTI.diff.AnnTemp.LGM_cur.Biome.plot + theme(axis.title.x = element_blank()),
          NRI.diff.AnnTemp.LGM_cur.Biome.plot + theme(axis.title.x = element_blank()),
          NRI.diff.AnnTemp.LGM_cur.Ecoregion.plot,
          NTI.diff.AnnTemp.LGM_cur.Ecoregion.plot,
          ncol = 2, nrow = 6,
          align = "hv",
          common.legend = TRUE)


ggsave(filename = "NRI.diff.AnnTemp.LGM_cur.ALL.plot.png", 
       dpi = 300,
       width = 5*2, height = 5.25*6, 
       units = "in")
