# S2d_ClimateDiff_PhyloStr_Plate_Representation

MPD.MNTD.diff.worldClimate.Plate <- diff.worldClimate %>%
  left_join(cbind(filter(MPD.LatLong.Env.AllScales, 
                         SamplingPool == "Plate sampling")[, c(1, 4:9, 29:35, 38)],
                  filter(MNTD.LatLong.Env.AllScales, 
                         SamplingPool == "Plate sampling")[, c(30:35, 37, 38)]), 
            by = c("ID" = "ID")) %>%
  filter(!is.na(ID_Realm))

# head(MPD.MNTD.LatLong.Plate)

MPD.MNTD.diff.worldClimate.Plate$ID_Realm <- factor(MPD.MNTD.diff.worldClimate.Plate$ID_Realm, 
                                                    levels = c("Nearctic", 
                                                               "Neotropical", 
                                                               "Afrotropical", 
                                                               "Palearctic",
                                                               "Indomalay",
                                                               "Australasian",
                                                               "Oceanic"))

# MPD.LatLong.Env.AllScales %>% count(SamplingPool)
# MNTD.LatLong.Env.AllScales %>% count(SamplingPool)

plot(MPD.MNTD.diff.worldClimate.Global$ntaxa, MPD.MNTD.diff.worldClimate.Global$NRI)
plot(MPD.MNTD.diff.worldClimate.Hemispheric$ntaxa, MPD.MNTD.diff.worldClimate.Hemispheric$NRI)
plot(MPD.MNTD.diff.worldClimate.Plate$ntaxa, MPD.MNTD.diff.worldClimate.Plate$NRI)

###############################
#### Statistical analyses #####
###############################

###############################
#### Representation  #####
###############################

(NRI.diff.AnnTemp.LGM_cur.Plate.plot <- ggplot(MPD.MNTD.diff.worldClimate.Plate, 
                                               aes(x = diff.AnnTemp.LGM_cur/10, 
                                                   y = NRI,
                                                   group = factor(ID_Realm),
                                                   colour = factor(ID_Realm))) +
   geom_point(alpha = 0.4, size = 1.75) +
   scale_color_viridis(option="D", 
                       begin = 0,
                       end = 1,
                       direction = 1, 
                       discrete = TRUE,
                       name="Plate") +
   #   geom_smooth(method = "loess", span = 0.8, se = TRUE, fullrange = FALSE) +
   labs(x = c(expression(atop("Historical Temperature Change (°C)", 
                              scriptstyle("MAT"[Contemporary]-"MAT"[LGM])))),
        y = c(expression("NRI"["Plate"]))) +
   scale_y_continuous(breaks = pretty(MPD.MNTD.diff.worldClimate.Plate$NRI, n = 7)) +
   geom_hline(yintercept = 0, alpha = 0.5) +
   theme_classic() +
   theme(legend.position = "bottom", 
         legend.direction = "horizontal",
         legend.key.size = unit(0.5, "cm"),
         legend.text = element_text(size = 9),
         legend.title = element_blank(),
         legend.box = "horizontal",
         axis.text.x = element_text(face = "bold", size = 14),
         axis.text.y = element_text(face = "bold", size = 14),
         axis.title = element_text(size = 16, face="bold")) + 
   guides(colour = guide_legend(nrow = 1)))

# NRI.diff.AnnTemp.LGM_cur.Plate.plot

ggsave(filename = "NRI.diff.AnnTemp.LGM_cur.Plate.D.plot.png", 
       dpi = 300,
       width = 7.5, height = 7.5, 
       units = "in")

# Divide by biogeographical Plate, in the horizontal direction
NRI.diff.AnnTemp.LGM_cur.Plate.plot + 
  facet_grid(. ~ ID_Realm, scales = "free") +
  theme(strip.background = element_rect(color="white", linetype = NULL),
        strip.text=element_text(color = "black",
                                face = "bold",
                                size = 15),
        legend.position = "")

ggsave(filename = "NRI.diff.AnnTemp.LGM_cur.Plate.D.Split.plot.png", 
       dpi = 300, 
       width = 22.5, height = 5, 
       units = "in")

#### AnnPrec plot ####

(NRI.diff.AnnPrec.LGM_cur.Plate.plot <- ggplot(MPD.MNTD.diff.worldClimate.Plate, 
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
                       name="Plate") +
   #   geom_smooth(method = "loess", span = 0.8, se = TRUE, fullrange = FALSE) +
   labs(x = c(expression(atop("Historical Precipitation Change (mm)", 
                              scriptstyle("MAP"[Contemporary]-"MAP"[LGM])))),
        y = c(expression("NRI"["Plate"]))) +
   scale_y_continuous(breaks = pretty(MPD.MNTD.diff.worldClimate.Plate$NRI, n = 7)) +
   #  scale_x_continuous(trans = "pseudo_log") +
   geom_hline(yintercept = 0, alpha = 0.5) +
   theme_classic() +
   theme(legend.position = "bottom", 
         legend.direction = "horizontal",
         legend.key.size = unit(0.5, "cm"),
         legend.text = element_text(size = 9),
         legend.title = element_blank(),
         legend.box = "horizontal",
         axis.text.x = element_text(face = "bold", size = 14),
         axis.text.y = element_text(face = "bold", size = 14),
         axis.title = element_text(size = 16, face="bold")) + 
   guides(colour = guide_legend(nrow = 1)))

# NRI.diff.AnnPrec.LGM_cur.Plate.plot

ggsave(filename = "NRI.diff.AnnPrec.LGM_cur.Plate.D.plot.png", 
       dpi = 300, 
       width = 7.5, height = 7.5, 
       units = "in")

# Divide by biogeographical Plate, in the horizontal direction
NRI.diff.AnnPrec.LGM_cur.Plate.plot + 
  facet_grid(. ~ ID_Realm, scales = "free") +
  theme(strip.background = element_rect(color="white", linetype = NULL),
        strip.text=element_text(color = "black",
                                face = "bold",
                                size = 15),
        legend.position = "")

ggsave(filename = "NRI.diff.AnnPrec.LGM_cur.Plate.D.Split.plot.png", 
       dpi = 300, 
       width =22.5, height = 5, 
       units = "in")


# NTI

## AnnTemp plot #####

(NTI.diff.AnnTemp.LGM_cur.Plate.plot <- ggplot(MPD.MNTD.diff.worldClimate.Plate, 
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
                       name="Plate") +
   # geom_smooth(method = "loess", span = 0.8, se = TRUE, fullrange = FALSE) +
   labs(x = c(expression(atop("Historical Temperature Change (°C)", 
                              scriptstyle("MAT"[Contemporary]-"MAT"[LGM])))),
        y = c(expression("NTI"["Plate"]))) +
   scale_y_continuous(breaks = pretty(MPD.MNTD.diff.worldClimate.Plate$NTI, n = 7)) +
   geom_hline(yintercept = 0, alpha = 0.5) +
   theme_classic() +
   theme(legend.position = "bottom", 
         legend.direction = "horizontal",
         legend.key.size = unit(0.5, "cm"),
         legend.text = element_text(size = 9),
         legend.title = element_blank(),
         legend.box = "horizontal",
         axis.text.x = element_text(face = "bold", size = 14),
         axis.text.y = element_text(face = "bold", size = 14),
         axis.title = element_text(size = 16, face="bold")) + 
   guides(colour = guide_legend(nrow = 1)))

NTI.diff.AnnTemp.LGM_cur.Plate.plot

ggsave(filename = "NTI.diff.AnnTemp.LGM_cur.Plate.D.plot.png", 
       dpi = 300,
       width = 7.5, height = 7.5, 
       units = "in")


# Divide by biogeographical Plate, in the horizontal direction
NTI.diff.AnnTemp.LGM_cur.Plate.plot + 
  facet_grid(. ~ ID_Realm, scales = "free") +
  theme(strip.background = element_rect(color="white", linetype = NULL),
        strip.text=element_text(color = "black",
                                face = "bold",
                                size = 15),
        legend.position = "")


ggsave(filename = "NTI.diff.AnnTemp.LGM_cur.Plate.D.Split.plot.png", 
       dpi = 300, 
       width =22.5, height = 5, 
       units = "in")


# Regions that underwent high losses of temperature during the Last Glacial Maximum are mainly composed of bat communities less
# phylognetically structured, in relation to regions that underwent less drastic changes in climate.

#### AnnPrec plot ####


(NTI.diff.AnnPrec.LGM_cur.Plate.plot <- ggplot(MPD.MNTD.diff.worldClimate.Plate, 
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
                       name="Plate") +
   #   geom_smooth(method = "loess", span = 0.8, se = TRUE, fullrange = FALSE) +
   labs(x = c(expression(atop("Historical Precipitation Change (mm)", 
                              scriptstyle("MAP"[Contemporary]-"MAP"[LGM])))),
        y = c(expression("NTI"["Plate"]))) +
   scale_y_continuous(breaks = pretty(MPD.MNTD.diff.worldClimate.Plate$NTI, n = 7)) +
   #  scale_x_continuous(trans = "pseudo_log") +
   geom_hline(yintercept = 0, alpha = 0.5) +
   theme_classic() +
   theme(legend.position = "bottom", 
         legend.direction = "horizontal",
         legend.key.size = unit(0.5, "cm"),
         legend.text = element_text(size = 9),
         legend.title = element_blank(),
         legend.box = "horizontal",
         axis.text.x = element_text(face = "bold", size = 14),
         axis.text.y = element_text(face = "bold", size = 14),
         axis.title = element_text(size = 16, face="bold")) + 
   guides(colour = guide_legend(nrow = 1)))

NTI.diff.AnnPrec.LGM_cur.Plate.plot

ggsave(filename = "NTI.diff.AnnPrec.LGM_cur.Plate.D.plot.png", 
       dpi = 300, 
       width = 7.5, height = 7.5, 
       units = "in")

# Divide by biogeographical Plate, in the horizontal direction
NTI.diff.AnnPrec.LGM_cur.Plate.plot + 
  facet_grid(. ~ ID_Realm, scales = "free") +
  theme(strip.background = element_rect(color="white", linetype = NULL),
        strip.text=element_text(color = "black",
                                face = "bold",
                                size = 15),
        legend.position = "")

ggsave(filename = "NTI.diff.AnnPrec.LGM_cur.Plate.D.Split.plot.png", 
       dpi = 300, 
       width =22.5, height = 5, 
       units = "in")

