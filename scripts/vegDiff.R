MPD.LatLong.Env.AllScales.P.05$VegCoverPropDiff <- (MPD.LatLong.Env.AllScales.P.05$VegCoverMOD - MPD.LatLong.Env.AllScales.P.05$VegCoverLGM)/(MPD.LatLong.Env.AllScales.P.05$VegCoverLGM)


#######################################################################################################

Bats.MPD.VegCoverPropDiff <- lm(NRI ~ VegCoverPropDiff, data = MPD.LatLong.Env.AllScales.P.05)
summary(Bats.MPD.VegCoverPropDiff)
plot(Bats.MPD.VegCoverPropDiff)

Bats.MPD.VegCoverPropDiff.LMM <- lmer(NRI ~ VegCoverPropDiff + (1|ID_Realm) + (1|ID_Biome), 
                                 data = subset(MPD.LatLong.Env.AllScales.P.05, SamplingPool = "Global sampling"), REML=TRUE)
summary(Bats.MPD.VegCoverPropDiff.LMM)

Bats.MPD.Latitude <- lm(NRI ~ Latitude, data = subset(MPD.LatLong.Env.AllScales.P.05, SamplingPool = "Global sampling"))
summary(Bats.MPD.Latitude)

##

NRI.VegCoverPropDiff.Global.plot <- ggplot(subset(MPD.LatLong.Env.AllScales.P.05, SamplingPool = "Global sampling"), aes(VegCoverPropDiff/10, NRI)) +
  facet_wrap(~ID_Realm, nrow = 1, strip.position = "top") +
  geom_jitter(width = .5, size=1) +
  geom_smooth(method="loess", se=F) +
  # possible values for trans : 'log2', 'log10','sqrt'
  scale_x_continuous() +
  theme_bw() + 
  theme(strip.background=element_rect(fill = "black")) + 
  theme(strip.text=element_text(color = "white", face = "bold", size = 15), 
        #        axis.text.x = element_blank(),
        axis.text.y = element_text(face = "bold", 
                                   size = 14),
        axis.title=element_text(size = 16, face="bold")) +
  
  labs(y="NRI",
       x="Annual Mean Temperature", 
       title="Global sampling")

NRI.VegCoverPropDiff.Global.plot

##

NRI.VegCoverPropDiff.Realm.plot <- ggplot(subset(MPD.LatLong.Env.AllScales.P.05, SamplingPool = "Realm sampling"), aes(VegCoverPropDiff/10, NRI)) +
  facet_wrap(~ID_Realm, nrow = 1, strip.position = "top") +
  geom_jitter(width = .5, size=1) +
  geom_smooth(method="loess", se=F) +
  # possible values for trans : 'log2', 'log10','sqrt'
  scale_x_continuous() +
  theme_bw() + 
  theme(strip.background=element_rect(fill = "black")) + 
  theme(strip.text=element_text(color = "white", face = "bold", size = 15), 
        #        axis.text.x = element_blank(),
        axis.text.y = element_text(face = "bold", 
                                   size = 14),
        axis.title=element_text(size = 16, face="bold")) +
  
  labs(y="NRI",
       x="Annual Mean Temperature", 
       title="Realm sampling")

NRI.VegCoverPropDiff.Realm.plot

##

NRI.VegCoverPropDiff.boxplot <- ggplot(subset(MPD.LatLong.Env.AllScales.P.05, SamplingPool = "Realm sampling"), aes(VegCoverPropDiff/10, NRI)) + 
  geom_jitter(width = .5, size=1) +
  geom_smooth(method="loess", se=F) +
  # possible values for trans : 'log2', 'log10','sqrt'
  scale_x_continuous() +
  
  labs(y="NRI", 
       x="Annual Mean Temperature", 
       title="Global scale")

NRI.VegCoverPropDiff.boxplot

##

NRI.Latitude.boxplot <- ggplot(subset(MPD.LatLong.Env.AllScales.P.05, SamplingPool = "Global sampling"), aes(Latitude, NRI)) + 
  geom_jitter(width = .5, size=1) +
  geom_smooth(method="loess", se=F) +
  # possible values for trans : 'log2', 'log10','sqrt'
  scale_x_continuous() +
  
  labs(y="NRI", 
       x="Latitude", 
       title="Global scale")

#  scale_y_continuous(trans='log10') +
NRI.Latitude.boxplot

