S2_MergingPhyloStructure_MPD_MNTD 

MPD.MNTD.LatLong.AllScales <- MPD.LatLong.Env.AllScales %>%
  dplyr::select(-c(bio_1:bio_19)) %>%
  dplyr::rename(mpd.runs = runs,
                mpd.ntaxa = ntaxa) %>%
  unite(ID_SamplingPool, ID, SamplingPool, remove = FALSE) %>%
  left_join(MNTD.LatLong.Env.AllScales %>%
              dplyr::select(-c(bio_1:bio_19)) %>%
              dplyr::rename(mntd.runs = runs,
                     mntd.ntaxa = ntaxa) %>% 
              unite(ID_SamplingPool, ID, SamplingPool, remove = FALSE),
            by = "ID_SamplingPool")

# verify

MPD.MNTD.LatLong.AllScales <- MPD.LatLong.Env.AllScales %>%
  dplyr::select(-c(bio_1:bio_19)) %>%
  dplyr::rename(mpd.runs = runs,
         mpd.ntaxa = ntaxa) %>%
  unite(ID_SamplingPool, ID, SamplingPool, remove = FALSE) %>%
  left_join(MNTD.LatLong.Env.AllScales %>%
              dplyr::select(-c(Longitude:bio_19)) %>%
              dplyr::rename(mntd.runs = runs,
                     mntd.ntaxa = ntaxa) %>% 
              unite(ID_SamplingPool, ID, SamplingPool, 
                    remove = TRUE),
            by = "ID_SamplingPool")


MPD.MNTD.LatLong.AllScales$ID_Realm <- factor(MPD.MNTD.LatLong.AllScales$ID_Realm, 
                                              levels=c('Neotropical',
                                                       'Nearctic',
                                                       'Afrotropical',
                                                       'Palearctic', 
                                                       'Indomalay', 
                                                       'Australasian'))

MPD.MNTD.LatLong.AllScales$SamplingPool <-  factor(MPD.MNTD.LatLong.AllScales$SamplingPool,  
                                                   levels = c("Global sampling",
                                                              "Hemispheric sampling",
                                                              "Realm sampling",
                                                              "Plate sampling",
                                                              "Biome sampling",
                                                              "Ecoregion sampling"))