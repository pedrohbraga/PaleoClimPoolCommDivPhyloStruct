######################################################################################
### Climate dataset preparation and complementary analyses                      ######
#                                                                                    #
# Author: Pedro Henrique Pereira Braga                                               #
# Last Update: "2019-01-28"                                                          #
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
                         Time = c(0.01,
                                  (4.2 + 0.3) / 2,
                                  (8.326 + 4.2) / 2,
                                  (11.7 + 8.326) / 2,
                                  21,
                                  787,
                                  (3264 + 3025)/2,
                                  3300))


# Load data sets ####

current.raw <- read.csv("data/matrices/world_CHELSA_cur_V1_2B_r2_5m_50km.csv", 
                        row.names = 1, 
                        header = TRUE)

Pleisto.LGM.raw <- read.csv("data/matrices/world_chelsa_LGM_v1_2B_r2_5m_50km.csv", 
                            row.names = 1, 
                            header = TRUE)

# Add prefix to all variables of the objects

colnames(current.raw) <- paste("cur", 
                               colnames(current.raw),
                               sep = "_")

colnames(Pleisto.LGM.raw) <- paste("LGM", 
                                   colnames(Pleisto.LGM.raw), 
                                   sep = "_")

# Bind all columns to create a big data.frame

worldClimate <- cbind(current.raw,
                      Pleisto.LGM.raw)

str(worldClimate); dim(worldClimate)

# Check for NAs
(rows.NA <- c(row.names(which(is.na(worldClimate), 
                              arr.ind = TRUE))))

# Solve NAs using zoo::na.spline

worldClimate <- as.data.frame(zoo::na.spline(worldClimate))

# Identify rows that have NA in all environmental datasets
c(row.names(which(is.na(worldClimate), 
                  arr.ind = TRUE)
            )
  )

head(worldClimate); str(worldClimate)

# Calculate the differences in climate between now and LGM ####

worldClimate.diff <- data.frame(ID = 1:nrow(worldClimate),
                                diff.AnnTemp.LGM_cur = worldClimate$LGM_bio_1 - worldClimate$cur_bio_1,
                                diff.AnnPrec.LGM_cur = worldClimate$LGM_bio_12 - worldClimate$cur_bio_12,
                                LGM_bio_1 = worldClimate$LGM_bio_1, 
                                cur_bio_1 = worldClimate$cur_bio_1,   
                                LGM_bio_12 = worldClimate$LGM_bio_12, 
                                cur_bio_12 = worldClimate$cur_bio_12)

MPD.MNTD.LatLong.AllScales.worldClimate.diff <- MPD.LatLong.AllScales %>%
  right_join(worldClimate.diff,
             by = "ID")

MPD.MNTD.LatLong.AllScales.rarefaction.relative.worldClimate.diff <- MPD.LatLong.AllScales.rarefaction.relative %>%
  right_join(worldClimate.diff,
             by = "ID")

MPD.MNTD.LatLong.AllScales.worldClimate.diff$ID_Realm <- factor(MPD.MNTD.LatLong.AllScales.worldClimate.diff$ID_Realm, 
                                                                levels = c("Nearctic", 
                                                                           "Neotropical", 
                                                                           "Afrotropical", 
                                                                           "Palearctic",
                                                                           "Indomalay",
                                                                           "Australasian",
                                                                           "Oceanic")
                                                                )




################################################################################ 
# Subset variables and explore the mean values through time and represent them #
################################################################################

## bio_1

(meanAnnTemp <- worldClimate %>% 
   select(ends_with("bio_1")) %>%
   bind_cols(world_grid_50km_cat_df[ , "ID_Realm", drop = FALSE]) %>%
   group_by(ID_Realm) %>% 
   summarise_all(funs(mean)) %>%
   gather(Period, bio_1,-ID_Realm) %>%
   transform(Period = str_replace(Period, "_bio_1", "")) %>%
   full_join(period.kYA, by = c("Period" = "Code")))

## bio_12

(meanAnnPrec <- worldClimate %>% 
    select(ends_with("bio_12")) %>%
    bind_cols(regions.LatLong[ ,6, drop = FALSE]) %>%
    group_by(ID_Realm) %>% 
    summarise_all(funs(mean)) %>%
    gather(Period, bio_12,-ID_Realm) %>%
    transform(Period = str_replace(Period, "_bio_12", "")) %>%
    full_join(period.kYA, by = c("Period" = "Code"))) 

### Represent mean temperature and precipitation change through time ####

meanAnnTemp$ID_Realm <- factor(meanAnnTemp$ID_Realm, 
                               levels = c("Nearctic", 
                                          "Neotropical", 
                                          "Afrotropical", 
                                          "Palearctic",
                                          "Indomalay",
                                          "Australasian"))

meanAnnTemp$factorTime <- as.factor(meanAnnTemp$Time)
# meanAnnTemp$ID_Realm<-as.factor(meanAnnTemp$ID_Realm)

levels(meanAnnTemp$Period.y) <- period.kYA$Period

# Plot the mean annual temperature across the time steps to visualize how it
# changed from the past to now

ggplot(data=filter(meanAnnTemp,
                   is.na(ID_Realm) == FALSE), 
       aes(x=factorTime, y=bio_1/10, 
           group=ID_Realm, 
           color=ID_Realm)) +
  stat_smooth(aes(x = factorTime, y = bio_1/10), method = "gam",
              formula = y ~ poly(x, 7), se = FALSE) +
  geom_line() + 
  geom_point() +
  scale_y_continuous(position = "right") +
  scale_x_discrete(limits = rev(levels(meanAnnTemp$factorTime))) +
  #  geom_text(aes(label=ifelse(ID_Realm == "Palearctic",
  #                             as.character(Period.y),'')),
  #            hjust=0, vjust=0) +
  #  scale_x_reverse() +
  #  geom_vline(data=filter(meanAnnTemp,
  #                         is.na(ID_Realm) == FALSE), mapping=aes(xintercept=Time), color="blue") +
  #  geom_text(mapping=aes(x=factorTime, y=min(bio_1/10), label=Period.y), size=4, angle=90, vjust=-0.4, hjust=0) +
  #  opts(title="geom_vline", plot.title=theme_text(size=40, vjust=1.5)) +
  scale_color_viridis(option="D", 
                      begin = 0,
                      end = 1,
                      direction = 1, 
                      discrete=TRUE,
                      name="Realm") +
  labs(x = "kYA",
       y = c(expression("Mean Annual Temperature (°C)"))) +
  geom_hline(yintercept = 0, alpha = 0.5) +
  theme_classic() +
  theme(legend.position = "bottom", 
        legend.direction = "horizontal",
        legend.key.size = unit(0.5, "cm"),
        legend.text = element_text(size = 9),
        legend.title = element_text(face = "bold", 
                                    size = 9),
        legend.box = "horizontal",
        axis.text.x = element_text(angle = 0, #hjust = 1, 
                                   face = "bold", size = 14),
        axis.text.y = element_text(face = "bold", size = 14),
        axis.title = element_text(size = 16, face="bold")) +
  guides(colour = guide_legend(nrow = 1))

ggsave(filename = "meanAnnTemp.Time.png", 
       dpi = 300, 
       width =7.5, height = 7.5, 
       units = "in")

###########################
#### For precipitation ####
###########################

# Plot the mean annual precipitation across the time steps to visualize how it
# changed from the past to now

meanAnnPrec$ID_Realm <- factor(meanAnnPrec$ID_Realm, 
                               levels = c("Nearctic", 
                                          "Neotropical", 
                                          "Afrotropical", 
                                          "Palearctic",
                                          "Indomalay",
                                          "Australasian"))

meanAnnPrec$factorTime <- as.factor(meanAnnPrec$Time)
# meanAnnPrec$ID_Realm<-as.factor(meanAnnPrec$ID_Realm)
levels(meanAnnPrec$Period.y) <- period.kYA$Period

ggplot(data=filter(meanAnnPrec,
                   is.na(ID_Realm) == FALSE), 
       aes(x=factorTime, y=bio_12, 
           group=ID_Realm, 
           color=ID_Realm)) +
  stat_smooth(aes(x = factorTime, y = bio_12), method = "gam",
              formula = y ~ poly(x, 7), se = FALSE) +
  geom_line() + 
  geom_point() +
  scale_y_continuous(position = "right") +
  scale_x_discrete(limits = rev(levels(meanAnnPrec$factorTime))) +
  #  geom_text(aes(label=ifelse(ID_Realm == "Palearctic",
  #                             as.character(Period.y),'')),
  #            hjust=0, vjust=0) +
  #  scale_x_reverse() +
  #  geom_vline(data=filter(meanAnnPrec,
  #                         is.na(ID_Realm) == FALSE), mapping=aes(xintercept=Time), color="blue") +
  #  geom_text(mapping=aes(x=factorTime, y=min(bio_12/10), label=Period.y), size=4, angle=90, vjust=-0.4, hjust=0) +
  #  opts(title="geom_vline", plot.title=theme_text(size=40, vjust=1.5)) +
  scale_color_viridis(option="D", 
                      begin = 0,
                      end = 1,
                      direction = 1, 
                      discrete=TRUE,
                      name="Realm") +
  labs(x = "kYA",
       y = c(expression("Mean Annual Precipitation (mm)"))) +
  geom_hline(yintercept = 0, alpha = 0.5) +
  theme_classic() +
  theme(legend.position = "bottom", 
        legend.direction = "horizontal",
        legend.key.size = unit(0.5, "cm"),
        legend.text = element_text(size = 9),
        legend.title = element_text(face = "bold", 
                                    size = 9),
        legend.box = "horizontal",
        axis.text.x = element_text(angle = 0, #hjust = 1, 
                                   face = "bold", size = 14),
        axis.text.y = element_text(face = "bold", size = 14),
        axis.title = element_text(size = 16, face="bold")) +
  guides(colour = guide_legend(nrow = 1))

ggsave(filename = "meanAnnPrec.Time.png", 
       dpi = 300, 
       width =7.5, height = 7.5, 
       units = "in")



summary(lm(filter(MPD.MNTD.LatLong.AllScales.worldClimate.diff, 
                  SamplingPool == "Global sampling")$ses.mpd.z.query ~ diff.AnnTemp.LGM_cur))

# Richness and NTI
plot(MPD.MNTD.LatLong.diff.Global$ntaxa, MPD.MNTD.LatLong.diff.Global$NTI, na.rm = TRUE)

# Richness and NRI
plot(MPD.MNTD.LatLong.diff.Global$ntaxa, MPD.MNTD.LatLong.diff.Global$NRI, na.rm = TRUE)

## AnnTemp plot #####

## Boxplot vs. Quantil. Regression
## SD

(NRI.diff.AnnTemp.LGM_cur.plot <- ggplot(MPD.MNTD.LatLong.AllScales.rarefaction.relative.worldClimate.diff %>%
                                           filter(SamplingPool == "Global sampling"), 
                                         aes(x = diff.AnnTemp.LGM_cur/10, 
                                             y = nri.rarefac.mean,
                                             group = factor(ID_Realm),
                                             colour = factor(ID_Realm))) +
   geom_point(alpha = 0.4, size = 1.75) +
   scale_color_viridis(option="D", 
                       begin = 0,
                       end = 1,
                       direction = 1, 
                       discrete = TRUE,
                       name="Realm") +
   #   geom_smooth(method = "loess", span = 0.8, se = TRUE, fullrange = FALSE) +
   labs(x = c(expression(atop("Historical Temperature Change (°C)", 
                              scriptstyle("MAT"[Contemporary]-"MAT"[LGM])))),
        y = c(expression("NRI"["Global"]))) +
   scale_y_continuous(breaks = pretty(MPD.MNTD.LatLong.AllScales.worldClimate.diff$nri, n = 7)) +
   geom_hline(yintercept = 0, alpha = 0.5) +
   theme_classic() +
   theme(legend.position = "bottom", 
         legend.direction = "horizontal",
         legend.key.size = unit(0.5, "cm"),
         legend.text = element_text(size = 9),
         legend.title = element_text(face = "bold", size = 9),
         legend.box = "horizontal",
         axis.text.x = element_text(face = "bold", size = 14),
         axis.text.y = element_text(face = "bold", size = 14),
         axis.title = element_text(size = 16, face="bold")) + 
   guides(colour = guide_legend(nrow = 1)))

NRI.diff.AnnTemp.LGM_cur.plot

ggsave(filename = "NRI.diff.AnnTemp.LGM_cur.Global.D.plot.png", 
       dpi = 300,
       width = 7.5, height = 7.5, 
       units = "in")


# Divide by biogeographical realm, in the horizontal direction
NRI.diff.AnnTemp.LGM_cur.plot + 
  facet_grid(. ~ ID_Realm., scales = "free") +
  theme(strip.background = element_rect(color="white", linetype = NULL),
        strip.text=element_text(color = "black",
                                face = "bold",
                                size = 15),
        legend.position = "")


ggsave(filename = "NRI.diff.AnnTemp.LGM_cur.Global.D.Split.plot.png", 
       dpi = 300, 
       width =22.5, height = 5, 
       units = "in")


# Regions that underwent high losses of temperature during the Last Glacial
# Maximum are mainly composed of bat communities less phylognetically
# structured, in relation to regions that underwent less drastic changes in
# climate.

#### AnnPrec plot ####

(NRI.diff.AnnPrec.LGM_cur.plot <- ggplot(MPD.MNTD.LatLong.diff.Global, 
                                         aes(x = diff.AnnPrec.LGM_cur/10, 
                                             y = NRI,
                                             group = factor(ID_Realm),
                                             colour = factor(ID_Realm))) +
   geom_point(alpha = 0.4, size = 1.75) +
   scale_color_viridis(option="D", 
                       begin = 0,
                       end = 1,
                       direction = 1, 
                       discrete = TRUE,
                       name="Realm") +
   #   geom_smooth(method = "loess", span = 0.8, se = TRUE, fullrange = FALSE) +
   labs(x = c(expression(atop("Historical Precipitation Change (mm)", 
                              scriptstyle("MAP"[Contemporary]-"MAP"[LGM])))),
        y = c(expression("NRI"["Global"]))) +
   scale_y_continuous(breaks = pretty(MPD.MNTD.LatLong.diff.Global$NRI, n = 7)) +
   #  scale_x_continuous(trans = "pseudo_log") +
   geom_hline(yintercept = 0, alpha = 0.5) +
   theme_classic() +
   theme(legend.position = "bottom", 
         legend.direction = "horizontal",
         legend.key.size = unit(0.5, "cm"),
         legend.text = element_text(size = 9),
         legend.title = element_text(face = "bold", size = 9),
         legend.box = "horizontal",
         axis.text.x = element_text(face = "bold", size = 14),
         axis.text.y = element_text(face = "bold", size = 14),
         axis.title = element_text(size = 16, face="bold")) + 
   guides(colour = guide_legend(nrow = 1)))

NRI.diff.AnnPrec.LGM_cur.plot

ggsave(filename = "NRI.diff.AnnPrec.LGM_cur.Global.D.plot.png", 
       dpi = 300, 
       width = 7.5, height = 7.5, 
       units = "in")

# Divide by biogeographical realm, in the horizontal direction
NRI.diff.AnnPrec.LGM_cur.plot + 
  facet_grid(. ~ ID_Realm., scales = "free") +
  theme(strip.background = element_rect(color="white", linetype = NULL),
        strip.text=element_text(color = "black",
                                face = "bold",
                                size = 15),
        legend.position = "")

ggsave(filename = "NRI.diff.AnnPrec.LGM_cur.Global.D.Split.plot.png", 
       dpi = 300, 
       width =22.5, height = 5, 
       units = "in")


# NTI

## AnnTemp plot #####

(NTI.diff.AnnTemp.LGM_cur.plot <- ggplot(MPD.MNTD.LatLong.diff.Global, 
                                         aes(x = diff.AnnTemp.LGM_cur/10, 
                                             y = NTI,
                                             group = factor(ID_Realm),
                                             colour = factor(ID_Realm))) +
   geom_point(alpha = 0.4, size = 1.75) +
   scale_color_viridis(option="D", 
                       begin = 0,
                       end = 1,
                       direction = 1, 
                       discrete = TRUE,
                       name="Realm") +
   # geom_smooth(method = "loess", span = 0.8, se = TRUE, fullrange = FALSE) +
   labs(x = c(expression(atop("Historical Temperature Change (°C)", 
                              scriptstyle("MAT"[Contemporary]-"MAT"[LGM])))),
        y = c(expression("NTI"["Global"]))) +
   scale_y_continuous(breaks = pretty(MPD.MNTD.LatLong.diff.Global$NTI, n = 7)) +
   geom_hline(yintercept = 0, alpha = 0.5) +
   theme_classic() +
   theme(legend.position = "bottom", 
         legend.direction = "horizontal",
         legend.key.size = unit(0.5, "cm"),
         legend.text = element_text(size = 9),
         legend.title = element_text(face = "bold", size = 9),
         legend.box = "horizontal",
         axis.text.x = element_text(face = "bold", size = 14),
         axis.text.y = element_text(face = "bold", size = 14),
         axis.title = element_text(size = 16, face="bold")) + 
   guides(colour = guide_legend(nrow = 1)))

NTI.diff.AnnTemp.LGM_cur.plot

ggsave(filename = "NTI.diff.AnnTemp.LGM_cur.Global.D.plot.png", 
       dpi = 300,
       width = 7.5, height = 7.5, 
       units = "in")


# Divide by biogeographical realm, in the horizontal direction
NTI.diff.AnnTemp.LGM_cur.plot + 
  facet_grid(. ~ ID_Realm., scales = "free") +
  theme(strip.background = element_rect(color="white", linetype = NULL),
        strip.text=element_text(color = "black",
                                face = "bold",
                                size = 15),
        legend.position = "")


ggsave(filename = "NTI.diff.AnnTemp.LGM_cur.Global.D.Split.plot.png", 
       dpi = 300, 
       width =22.5, height = 5, 
       units = "in")


# Regions that underwent high losses of temperature during the Last Glacial Maximum are mainly composed of bat communities less
# phylognetically structured, in relation to regions that underwent less drastic changes in climate.

#### AnnPrec plot ####


(NTI.diff.AnnPrec.LGM_cur.plot <- ggplot(MPD.MNTD.LatLong.diff.Global, 
                                         aes(x = diff.AnnPrec.LGM_cur/10, 
                                             y = NTI,
                                             group = factor(ID_Realm),
                                             colour = factor(ID_Realm))) +
   geom_point(alpha = 0.4, size = 1.75) +
   scale_color_viridis(option="D", 
                       begin = 0,
                       end = 1,
                       direction = 1, 
                       discrete = TRUE,
                       name="Realm") +
   #   geom_smooth(method = "loess", span = 0.8, se = TRUE, fullrange = FALSE) +
   labs(x = c(expression(atop("Historical Precipitation Change (mm)", 
                              scriptstyle("MAP"[Contemporary]-"MAP"[LGM])))),
        y = c(expression("NTI"["Global"]))) +
   scale_y_continuous(breaks = pretty(MPD.MNTD.LatLong.diff.Global$NTI, n = 7)) +
   #  scale_x_continuous(trans = "pseudo_log") +
   geom_hline(yintercept = 0, alpha = 0.5) +
   theme_classic() +
   theme(legend.position = "bottom", 
         legend.direction = "horizontal",
         legend.key.size = unit(0.5, "cm"),
         legend.text = element_text(size = 9),
         legend.title = element_text(face = "bold", size = 9),
         legend.box = "horizontal",
         axis.text.x = element_text(face = "bold", size = 14),
         axis.text.y = element_text(face = "bold", size = 14),
         axis.title = element_text(size = 16, face="bold")) + 
   guides(colour = guide_legend(nrow = 1)))

NTI.diff.AnnPrec.LGM_cur.plot

ggsave(filename = "NTI.diff.AnnPrec.LGM_cur.Global.D.plot.png", 
       dpi = 300, 
       width = 7.5, height = 7.5, 
       units = "in")

# Divide by biogeographical realm, in the horizontal direction

NTI.diff.AnnPrec.LGM_cur.plot + 
  facet_grid(. ~ ID_Realm., scales = "free") +
  theme(strip.background = element_rect(color="white", linetype = NULL),
        strip.text=element_text(color = "black",
                                face = "bold",
                                size = 15),
        legend.position = "")

ggsave(filename = "NTI.diff.AnnPrec.LGM_cur.Global.D.Split.plot.png", 
       dpi = 300, 
       width =22.5, height = 5, 
       units = "in")

