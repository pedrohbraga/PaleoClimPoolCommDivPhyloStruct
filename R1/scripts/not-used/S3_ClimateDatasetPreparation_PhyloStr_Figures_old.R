######################################################################################
### Climate dataset preparation and complementary analyses                      ######
#                                                                                    #
# Author: Pedro Henrique Pereira Braga                                               #
# Last Update: "2022-03-28"                                                          #
#                                                                                    #
######################################################################################

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

period.kYA <- data.frame(
  Period = c(
    "Current",
    "Late-Holocene",
    "Mid-Holocene",
    "Early-Holocene",
    "LGM",
    "Brunhes–Matuyama reversal",
    "Mid-Pliocene warm period",
    "Pliocene (M2)"
  ),
  Code = c(
    "cur",
    "LH",
    "MH",
    "EH",
    "LGM",
    "MIS19",
    "mPWP",
    "mP2"
  ),
  Time = c(
    0.01,
    (4.2 + 0.3) / 2,
    (8.326 + 4.2) / 2,
    (11.7 + 8.326) / 2,
    21,
    787,
    (3264 + 3025) / 2,
    3300
  )
)


# Load data sets ####

current.raw <- read.csv("data/matrices/world_CHELSA_cur_V1_2B_r2_5m_50km.csv",
                        row.names = 1,
                        header = TRUE
)

Pleisto.LGM.raw <- read.csv("data/matrices/world_chelsa_LGM_v1_2B_r2_5m_50km.csv",
                            row.names = 1,
                            header = TRUE
)

# Add prefix to all variables of the objects

colnames(current.raw) <- paste("cur",
                               colnames(current.raw),
                               sep = "_"
)

colnames(Pleisto.LGM.raw) <- paste("LGM",
                                   colnames(Pleisto.LGM.raw),
                                   sep = "_"
)

# Bind all columns to create a big data.frame

worldClimate <- cbind(
  current.raw,
  Pleisto.LGM.raw
)

str(worldClimate)
dim(worldClimate)

# Check for NAs
(rows.NA <- c(row.names(which(is.na(worldClimate),
                              arr.ind = TRUE
))))

# Solve NAs using zoo::na.spline

worldClimate <- as.data.frame(zoo::na.spline(worldClimate))

# Identify rows that have NA in all environmental datasets
c(row.names(which(is.na(worldClimate),
                  arr.ind = TRUE
)))

head(worldClimate)
str(worldClimate)

# Calculate the differences in climate between now and LGM ####

worldClimate.diff <- data.frame(
  ID = 1:nrow(worldClimate),
  diff.AnnTemp.LGM_cur = worldClimate$LGM_bio_1 - worldClimate$cur_bio_1,
  diff.AnnPrec.LGM_cur = worldClimate$LGM_bio_12 - worldClimate$cur_bio_12,
  LGM_bio_1 = worldClimate$LGM_bio_1,
  cur_bio_1 = worldClimate$cur_bio_1,
  LGM_bio_12 = worldClimate$LGM_bio_12,
  cur_bio_12 = worldClimate$cur_bio_12
) %>%
  mutate(
    log.LGM_bio_12 = log1p(LGM_bio_12),
    log.cur_bio_12 = log1p(cur_bio_12),
    diff.log.AnnPrec.log.LGM_cur = log.LGM_bio_12 - log.cur_bio_12,
    diff.AnnTemp.LGM_cur = diff.AnnTemp.LGM_cur / 10
  )

MPD.MNTD.LatLong.AllScales.raref.rel.worldClimate.diff <- MPD.MNTD.LatLong.AllScales.rarefaction.relative %>%
  right_join(worldClimate.diff,
             by = "ID"
  )

MPD.MNTD.LatLong.AllScales.raref.rel.worldClimate.diff$ID_Realm <- factor(MPD.MNTD.LatLong.AllScales.raref.rel.worldClimate.diff$ID_Realm,
                                                                levels = c(
                                                                  "Nearctic",
                                                                  "Neotropical",
                                                                  "Afrotropical",
                                                                  "Palearctic",
                                                                  "Indomalay",
                                                                  "Australasian",
                                                                  "Oceanic"
                                                                )
)

# Explore the outcomes

head(MPD.MNTD.LatLong.AllScales.raref.rel.worldClimate.diff)
dim(MPD.MNTD.LatLong.AllScales.raref.rel.worldClimate.diff)
