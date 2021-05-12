# S1b_ClimateDatasetPreparation

# Date: 2018-11-28

# Current (1979 – 2013): Anthropocene                         # cur_*
# Pleistocene: late-Holocene, Meghalayan (4.2-0.3 ka)         # LH_*
# Pleistocene: mid-Holocene, Northgrippian (8.326-4.2 ka)     # MH_*
# Pleistocene: early-Holocene, Greenlandian (11.7-8.326 ka)   # EH_*
# Pleistocene: Last Glacial Maximum (ca. 21 ka)               # LGM_*
# Pleistocene: MIS19 (ca. 787 ka)  Brunhes–Matuyama reversal  # MIS19_*
# Pliocene: mid-Pliocene warm period (3.264-3.025 Ma)         # mPWP_*
# Pliocene: M2 (ca. 3.3 Ma)                                   # mP2_*

# Brown, Hill, Dolan, Carnaval, Haywood (2018) 
# PaleoClim, high spatial resolution paleoclimate surfaces for global land areas.  
# Nature – Scientific Data. 5:180254

# Bio_1=Annual Mean Temperature [°C*10]
# Bio_2=Mean Diurnal Range [°C]
# Bio_3=Isothermality [Bio_2/Bio_7]
# Bio_4=Temperature Seasonality [standard deviation*100]
# Bio_5=Max Temperature of Warmest Month [°C*10]
# Bio_6=Min Temperature of Coldest Month [°C*10]
# Bio_7=Temperature Annual Range [°C*10]
# Bio_8=Mean Temperature of Wettest Quarter [°C*10]
# Bio_9=Mean Temperature of Driest Quarter [°C*10]
# Bio_10=Mean Temperature of Warmest Quarter [°C*10]
# Bio_11=Mean Temperature of Coldest Quarter [°C*10]
# Bio_12=Annual Precipitation [mm/year]
# Bio_13=Precipitation of Wettest Month [mm/month]
# Bio_14=Precipitation of Driest Month [mm/month]
# Bio_15=Precipitation Seasonality [coefficient of variation]
# Bio_16=Precipitation of Wettest Quarter [mm/quarter]
# Bio_17=Precipitation of Driest Quarter [mm/quarter]
# Bio_18=Precipitation of Warmest Quarter [mm/quarter]
# Bio_19=Precipitation of Coldest Quarter [mm/quarter]

period.kYA <- data.frame(Period = c("Current", 
                                    "Late-Holocene",
                                    "Mid-Holocene",
                                    "Early-Holocene",
                                    "LGM",
                                    "Brunhes–Matuyama reversal",
                                    "Mid-Pliocene warm period",
                                    "Pliocene (M2)"),
                         Code = c("cur",
                                  "LH", 
                                  "MH",
                                  "EH", 
                                  "LGM", 
                                  "MIS19", 
                                  "mPWP", 
                                  "mP2"),
                         Time = c(1,
                                  (4.2 - 0.3) / 2,
                                  (8.326 - 4.2) / 2,
                                  (11.7-8.326) / 2,
                                  21,
                                  787,
                                  (3264-3025)/2,
                                  3300))


# Load datasets ####

current.env.raw <- read.csv("data/matrices//world_CHELSA_cur_V1_2B_r2_5m_50km.csv" , row.names = 1, header=TRUE)

lateHolo.env.raw <- read.csv("data/matrices//world_LH_v1_2_5m_50km.csv" , row.names = 1, header=TRUE)
midHolo.env.raw <- read.csv("data/matrices//world_MH_v1_2_5m_50km.csv" , row.names = 1, header=TRUE)
earlyHolo.env.raw <- read.csv("data/matrices//world_EH_v1_2_5m_50km.csv" , row.names = 1, header=TRUE)

Pleisto.LGM.env.raw <- read.csv("data/matrices//world_chelsa_LGM_v1_2B_r2_5m_50km.csv" , row.names = 1, header=TRUE)

Pleisto.BrunMatRev.env.raw <- read.csv("data/matrices//world_MIS19_v1_r2_5m_50km.csv" , row.names = 1, header=TRUE)

midPliocene.env.raw <- read.csv("data/matrices//world_mPWP_v1_r2_5m_50km.csv" , row.names = 1, header=TRUE)

Pliocene.M2.env.raw <- read.csv("data/matrices//world_M2_v1_r2_5m_50km.csv" , 
                                row.names = 1, header=TRUE)

# Add prefix to all variables of the objects

colnames(current.env.raw) <- paste("cur", colnames(current.env.raw), sep = "_")

colnames(lateHolo.env.raw) <- paste("LH", colnames(lateHolo.env.raw), sep = "_")
colnames(midHolo.env.raw) <- paste("MH", colnames(midHolo.env.raw), sep = "_")
colnames(earlyHolo.env.raw) <- paste("EH", colnames(earlyHolo.env.raw), sep = "_")

colnames(Pleisto.LGM.env.raw) <- paste("LGM", colnames(Pleisto.LGM.env.raw), sep = "_")
colnames(Pleisto.BrunMatRev.env.raw) <- paste("MIS19", colnames(Pleisto.BrunMatRev.env.raw), sep = "_")

colnames(midPliocene.env.raw) <- paste("mPWP", colnames(midPliocene.env.raw), sep = "_")
colnames(Pliocene.M2.env.raw) <- paste("mP2", colnames(Pliocene.M2.env.raw), sep = "_")

# Bind all columns to create a big data.frame

worldClimate <- cbind(current.env.raw,
                      lateHolo.env.raw,
                      midHolo.env.raw,
                      earlyHolo.env.raw,
                      Pleisto.LGM.env.raw,
                      Pleisto.BrunMatRev.env.raw,
                      midPliocene.env.raw,
                      Pliocene.M2.env.raw)

str(worldClimate); dim(worldClimate)

# Check for NAs
(rows.NA <- c(row.names(which(is.na(worldClimate), 
                              arr.ind=TRUE))))

# Solve NAs using zoo::na.spline

worldClimate <- as.data.frame(zoo::na.spline(worldClimate))

# Identify rows that have NA in all environmental datasets
(rows.NA <- c(row.names(which(is.na(worldClimate), 
                              arr.ind=TRUE))))

head(worldClimate); str(worldClimate)

# Subset variables and explore the mean values through time

## bio_1

(meanAnnTemp <- worldClimate %>% 
    select(ends_with("bio_1")) %>%
    bind_cols(regions.LatLong[ ,6, drop = FALSE]) %>%
    group_by(ID_Realm) %>% 
    summarise_all(funs(mean)) %>%
    gather(Period, bio_1,-ID_Realm) %>%
    transform(Period = str_replace(Period, "_bio_1", "")) %>%
    full_join(period.kYA, by = c("Period" = "Code")))

(meanAnnPrec <- worldClimate %>% 
    select(ends_with("bio_12")) %>%
    bind_cols(regions.LatLong[ ,6, drop = FALSE]) %>%
    group_by(ID_Realm) %>% 
    summarise_all(funs(mean)) %>%
    gather(Period, bio_12,-ID_Realm) %>%
    transform(Period = str_replace(Period, "_bio_12", "")) %>%
    full_join(period.kYA, by = c("Period" = "Code")))


# Build a plot

# meanAnnTemp$ID_Realm<-as.factor(meanAnnTemp$ID_Realm)
levels(meanAnnTemp$Period.y) <- period.kYA$Period

ggplot(data=meanAnnTemp, aes(x=Period.y, y=bio_1, 
                             group=ID_Realm, 
                             color=ID_Realm)) +
  geom_line() + 
  geom_point() +
  scale_color_brewer(palette="Paired") +
  theme_minimal()

## Difference between now and LGM

head(worldClimate)

diff.AnnTemp.LGM_cur <- worldClimate$LGM_bio_1 - worldClimate$cur_bio_1
diff.AnnTemp.LGM_MIS19 <- worldClimate$LGM_bio_1 - worldClimate$MIS19_bio_1
diff.AnnPrec.LGM_cur <- worldClimate$LGM_bio_12 - worldClimate$cur_bio_12
diff.AnnPrec.LGM_MIS19 <- worldClimate$LGM_bio_12 - worldClimate$MIS19_bio_12

MPD.MNTD.LatLong.diff.Env.Global <- cbind(filter(MPD.LatLong.Env.AllScales, 
                                                 SamplingPool == "Global sampling")[, c(1, 4:7, 27:37)],
                                          filter(MNTD.LatLong.Env.AllScales, 
                                                 SamplingPool == "Global sampling")[, c(28:33, 36)],
                                          diff.AnnTemp.LGM_cur,
                                          diff.AnnTemp.LGM_MIS19,
                                          diff.AnnPrec.LGM_cur,
                                          diff.AnnPrec.LGM_MIS19,
                                          worldClimate$cur_bio_1)

MPD.LatLong.diffTemp.Global
gather(MPD.LatLong.diffTemp.Global, NRI.NTI, -1:)

MPD.MNTD.LatLong.diff.Env.Global <- MPD.LatLong.diffTemp.Global  %>%
  filter(!is.na(ID_Realm.))

summary(lm(filter(MPD.LatLong.Env.AllScales, SamplingPool == "Global sampling")$mpd.obs.z ~ diff.AnnTemp.cur_LGM))


## AnnTemp plot #####

NRI.diff.AnnTemp.cur_LGM.plot <- ggplot(MPD.MNTD.LatLong.diff.Env.Global, 
                                        aes(x = diff.AnnTemp.LGM_cur, y = NRI)) + 
  geom_point(alpha = 0.4, size = 0.9) +
  labs(x = c(expression(Delta~"Mean Annual Temperature (mm)"))) +
  scale_y_continuous(breaks = pretty(MPD.MNTD.LatLong.diff.Env.Global$NRI, n = 7)) +
  geom_hline(yintercept = 0, alpha = 0.5) +
  theme_classic()

NRI.diff.AnnTemp.cur_LGM.plot

# Divide by biogeographical realm, in the horizontal direction
NRI.diff.AnnTemp.cur_LGM.plot + 
  facet_grid(. ~ ID_Realm., scales = "free")
theme(strip.background = element_rect(color="white", linetype = NULL))


# Regions that underwent high losses of temperature during the Last Glacial Maximum are mainly composed of bat communities less
# phylognetically structured, in relation to regions that underwent less drastic changes in climate.

# AnnPrec plot ####

NRI.diff.AnnPrec.cur_LGM.plot <- ggplot(MPD.MNTD.LatLong.diff.Env.Global, 
                                        aes(x=diff.AnnPrec.LGM_cur, y=NRI)) + 
  geom_point(alpha = 0.5) +
  labs(x = c(expression(Delta~"Annual Precipitation (mm)"))) +
  scale_y_continuous(breaks = pretty(MPD.MNTD.LatLong.diff.Env.Global$NRI, n = 7)) +
  geom_hline(yintercept = 0, alpha = 0.5) +
  theme_classic()

NRI.diff.AnnPrec.cur_LGM.plot

# Divide by biogeographical realm, in the horizontal direction
NRI.diff.AnnPrec.cur_LGM.plot + 
  facet_grid(. ~ ID_Realm., scales = "free") +
  theme(strip.background = element_rect(color="white", linetype = NULL))

# NTI

NTI.diff.AnnPrec.cur_LGM.plot <- ggplot(MPD.MNTD.LatLong.diff.Env.Global, 
                                        aes(x=diff.AnnPrec.LGM_cur, y=NTI)) + 
  geom_point(alpha = 0.5) +
  #  geom_smooth(method="glm" , color="black") +
  theme_classic()  

NTI.diff.AnnPrec.cur_LGM.plot

# Divide by biogeographical realm, in the horizontal direction
NTI.diff.AnnPrec.cur_LGM.plot + 
  facet_grid(. ~ ID_Realm., scales = "free")

