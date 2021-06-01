length(Chiroptera.FaurbySven.Phylo.ultra$tip.label); ncol(Chiroptera.FaurbySven.comm)

########################
### Global sampling ####
########################

# Calculate mpd values for each community
SES.MPD.Chiroptera.World <- ses.mpd(Chiroptera.FaurbySven.comm,
                                  cophenetic(Chiroptera.FaurbySven.Phylo.ultra), # Convert phylogeny to distance matrix
                                  null.model = "taxa.labels",
                                  abundance.weighted = FALSE, 
                                  runs = 99)


MPD.LatLong.Env.World <- data.frame(cbind(regions.LatLong, 
                                          Env,
                                          SES.MPD.Chiroptera.World,
                                          SamplingPool = rep("Global sampling", nrow(regions.LatLong))))

head(MPD.LatLong.Env.World)


##########################################
#### New World and Old World Sampling ####
##########################################

NW.LatLong <- subset(regions.LatLong, ID_Realm == "Nearctic" | ID_Realm == "Neotropical")
OW.LatLong <- subset(regions.LatLong, !(ID_Realm %in% c("Nearctic", "Neotropical")))
nrow(NW.LatLong) + nrow(OW.LatLong)

NW.Comm <- Chiroptera.FaurbySven.comm[rownames(Chiroptera.FaurbySven.comm) %in% row.names(NW.LatLong), ]
OW.Comm <- Chiroptera.FaurbySven.comm[rownames(Chiroptera.FaurbySven.comm) %in% row.names(OW.LatLong), ]

dim(NW.Comm) ; nrow(NW.LatLong)
dim(OW.Comm) ; nrow(OW.LatLong)

### Subset the datasets

Chiroptera.NW.Comm <- match.phylo.comm(Mammal.FaurbySven.tree.1, NW.Comm)$comm
Chiroptera.NW.Tree <- match.phylo.comm(Mammal.FaurbySven.tree.1, NW.Comm)$phy

Chiroptera.OW.Comm <- match.phylo.comm(Mammal.FaurbySven.tree.1, OW.Comm)$comm
Chiroptera.OW.Tree <- match.phylo.comm(Mammal.FaurbySven.tree.1, OW.Comm)$phy

#############################################################

#### New World ####
# Calculate mpd values for each community
SES.MPD.Chiroptera.NW <- ses.mpd(Chiroptera.NW.Comm,
                                    cophenetic(Chiroptera.NW.Tree), # Convert phylogeny to distance matrix
                                    null.model = "taxa.labels",
                                    abundance.weighted = FALSE, 
                                    runs = 99)


#### Old World ####

SES.MPD.Chiroptera.OW <- ses.mpd(Chiroptera.OW.Comm,
                                 cophenetic(Chiroptera.OW.Tree), # Convert phylogeny to distance matrix
                                 null.model = "taxa.labels",
                                 abundance.weighted = FALSE, 
                                 runs = 99)

MPD.LatLong.Env.NW <- data.frame(cbind(regions.LatLong[(rownames(regions.LatLong) %in% row.names(Chiroptera.NW.Comm)), ], 
                                       Env[(rownames(Env) %in% row.names(Chiroptera.NW.Comm)), ],
                                       SES.MPD.Chiroptera.NW,
                                       SamplingPool = rep("Hemispheric sampling", length(MPD.Chiroptera.NW))))
dim(MPD.LatLong.Env.NW)

MPD.LatLong.Env.OW <- data.frame(cbind(regions.LatLong[(rownames(regions.LatLong) %in% row.names(Chiroptera.OW.Comm)), ], 
                                       Env[(rownames(Env) %in% row.names(Chiroptera.OW.Comm)), ],
                                       SES.MPD.Chiroptera.OW,
                                       SamplingPool = rep("Hemispheric sampling", length(MPD.Chiroptera.OW))))
dim(MPD.LatLong.Env.OW)

MPD.LatLong.Env.NW.OW <- rbind(MPD.LatLong.Env.NW, MPD.LatLong.Env.OW)

dim(MPD.LatLong.Env.NW.OW)

##########################################
############# Realm Sampling #############
##########################################

# Create a list
SubRegion.Realm.Comm.Phylo <- list()

for(i in droplevels(unique(regions.LatLong$ID_Realm))) {
  
  SubRegion.Realm.LatLong <- subset(regions.LatLong, ID_Realm == toString(i))
  SubRegion.Realm.Comm. <- Chiroptera.FaurbySven.comm[rownames(Chiroptera.FaurbySven.comm) %in% row.names(SubRegion.Realm.LatLong), ]
  
  ### Subset the datasets
  SubRegion.Realm.CommPhylo <- match.phylo.comm(Mammal.FaurbySven.tree.1, SubRegion.Realm.Comm.)
  SubRegion.Realm.Comm.Phylo[[paste0(i)]] <- list(Comm = SubRegion.Realm.CommPhylo$comm,
                                                  Phylo = SubRegion.Realm.CommPhylo$phy)
}

MPD.LatLong.Env.Realm = NULL

for(j in as.factor(ls(SubRegion.Realm.Comm.Phylo))){
  # Calculate mpd values for each community
  SES.MPD.Chiroptera.Realm <- ses.mpd(SubRegion.Realm.Comm.Phylo[[j]]$Comm,
                                    cophenetic(SubRegion.Realm.Comm.Phylo[[j]]$Phylo), # Convert phylogeny to distance matrix
                                    null.model = "taxa.labels",
                                    abundance.weighted = FALSE, 
                                    runs = 99)

  MPD.LatLong.Env.SubRegion.Realm <- data.frame(cbind(regions.LatLong[(rownames(regions.LatLong) %in% row.names(SubRegion.Realm.Comm.Phylo[[j]]$Comm)), ], 
                                                      Env[(rownames(Env) %in% row.names(SubRegion.Realm.Comm.Phylo[[j]]$Comm)), ],
                                                      SES.MPD.Chiroptera.Realm,
                                                      SamplingPool = rep("Realm sampling", nrow(SES.MPD.Chiroptera.Realm))))
  
  MPD.LatLong.Env.Realm <- rbind(MPD.LatLong.Env.Realm, MPD.LatLong.Env.SubRegion.Realm)
}

##########################################
############# Biome Sampling #############
##########################################

# Create a list
SubRegion.Biome.Comm.Phylo <- list()

for(i in droplevels(unique(regions.LatLong$ID_Biome))) {
  
  SubRegion.Biome.LatLong <- subset(regions.LatLong, ID_Biome == toString(i))
  SubRegion.Biome.Comm. <- Chiroptera.FaurbySven.comm[rownames(Chiroptera.FaurbySven.comm) %in% row.names(SubRegion.Biome.LatLong), ]
  
  ### Subset the datasets
  SubRegion.Biome.CommPhylo <- match.phylo.comm(Mammal.FaurbySven.tree.1, SubRegion.Biome.Comm.)
  SubRegion.Biome.Comm.Phylo[[paste0(i)]] <- list(Comm = SubRegion.Biome.CommPhylo$comm,
                                                  Phylo = SubRegion.Biome.CommPhylo$phy)
}

MPD.LatLong.Env.Biome = NULL

for(j in as.factor(ls(SubRegion.Biome.Comm.Phylo))){
  # Calculate standardized versions under the uniform model
  SES.MPD.Chiroptera.Biome <-   ses.mpd(SubRegion.Biome.Comm.Phylo[[j]]$Comm,
                                        cophenetic(SubRegion.Biome.Comm.Phylo[[j]]$Phylo), # Convert phylogeny to distance matrix
                                        null.model = "taxa.labels",
                                        abundance.weighted = FALSE, 
                                        runs = 99)

  MPD.LatLong.Env.SubRegion.Biome <- data.frame(cbind(regions.LatLong[(rownames(regions.LatLong) %in% row.names(SubRegion.Biome.Comm.Phylo[[j]]$Comm)), ], 
                                                      Env[(rownames(Env) %in% row.names(SubRegion.Biome.Comm.Phylo[[j]]$Comm)), ],
                                                      SES.MPD.Chiroptera.Biome,
                                                      SamplingPool = rep("Biome sampling", nrow(SES.MPD.Chiroptera.Biome))))
  
  MPD.LatLong.Env.Biome <- rbind(MPD.LatLong.Env.Biome, MPD.LatLong.Env.SubRegion.Biome)
}


#################################
#### Combine all data frames ####
#################################

MPD.LatLong.Env.AllScales <- rbind(MPD.LatLong.Env.World,
                                   MPD.LatLong.Env.NW.OW,
                                   MPD.LatLong.Env.Realm,
                                   MPD.LatLong.Env.Biome)

head(MPD.LatLong.Env.AllScales)



# Calculate NRI
MPD.LatLong.Env.AllScales$NRI <- -1 * ((MPD.LatLong.Env.AllScales$mpd.obs - MPD.LatLong.Env.AllScales$mpd.rand.mean) / MPD.LatLong.Env.AllScales$mpd.rand.sd)

# Rename dataset

MPD.LatLong.Env.AllScales$ID_Realm. <- factor(MPD.LatLong.Env.AllScales$ID_Realm, levels=c('Neotropical','Nearctic','Afrotropical','Palearctic', 'Indomalay', 'Australasian'))


MPD.LatLong.Env.AllScales.P.05 <- subset(MPD.LatLong.Env.AllScales, mpd.obs.p <= 0.05) 

#########################################################

NRI.Realm.boxplot <- ggplot(MPD.LatLong.Env.AllScales.P.05, aes(x=SamplingPool, y=NRI)) + 
  geom_boxplot(aes(fill=SamplingPool)) + 
  facet_wrap(~ID_Realm., nrow = 1, strip.position = "top") +
  #scale_y_continuous(limit = c(0, 15)) +
  geom_hline(yintercept=0) +
  labs(x="") +
  theme_bw() + 
  theme(strip.background=element_rect(fill = "black")) + 
  theme(strip.text=element_text(color = "white", face = "bold", size = 15), 
        axis.text.x = element_blank(),
        axis.text.y = element_text(face = "bold", 
                                   size = 14),
        axis.title=element_text(size = 16, face="bold"))

NRI.Realm.boxplot

##

NRI.Biome.boxplot <- ggplot(MPD.LatLong.Env.AllScales.P.05, aes(x=SamplingPool, y=NRI)) + 
  geom_boxplot(aes(fill=SamplingPool)) + 
  facet_wrap(~ID_Biome, nrow = 1, strip.position = "top") +
  #scale_y_continuous(limit = c(0, 15)) +
  geom_hline(yintercept=0) +
  labs(x="") +
  theme_bw() + 
  theme(strip.background=element_rect(fill = "black")) + 
  theme(strip.text=element_text(color = "white", face = "bold", size = 15), 
        axis.text.x = element_blank(),
        axis.text.y = element_text(face = "bold", 
                                   size = 14),
        axis.title=element_text(size = 16, face="bold"))

NRI.Biome.boxplot


########################################################

Bats.NRI.Realm <- lm(NRI ~ ID_Realm., data = MPD.LatLong.Env.AllScales.P.05)
plot(Bats.NRI.Realm)
summary(Bats.NRI.Realm)

TukeyHSD(x=Bats.NRI.Realm, 'ID_Realm.', conf.level=0.95)


Bats.NRI.Realm.LMM <- lmer(NRI ~ ID_Realm. + (1|SamplingPool), data = MPD.LatLong.Env.AllScales.P.05)
plot(Bats.NRI.Realm)
summary(Bats.NRI.Realm.LMM)

TukeyHSD(x=Bats.NRI.Realm, 'ID_Realm.', conf.level=0.95)


#######################################################################################################

Bats.MPD.AnMnTemp_Bio1 <- lm(NRI ~ AnMnTemp_Bio1, data = MPD.LatLong.Env.AllScales.P.05)
summary(Bats.MPD.AnMnTemp_Bio1)
plot(Bats.MPD.AnMnTemp_Bio1)

Bats.MPD.AnMnTemp_Bio1.LMM <- lmer(NRI ~ AnMnTemp_Bio1 + (1|ID_Realm) + (1|ID_Biome), 
                                   data = subset(MPD.LatLong.Env.AllScales.P.05, SamplingPool = "Global sampling"), REML=TRUE)
summary(Bats.MPD.AnMnTemp_Bio1.LMM)

Bats.MPD.Latitude <- lm(NRI ~ Latitude, data = subset(MPD.LatLong.Env.AllScales.P.05, SamplingPool = "Global sampling"))
summary(Bats.MPD.Latitude)

theme_set(theme_bw())  # pre-set the bw theme.

NRI.AnMnTemp_Bio1.boxplot <- ggplot(subset(MPD.LatLong.Env.AllScales.P.05, SamplingPool = "Global sampling"), aes(AnMnTemp_Bio1/10, NRI)) +
  facet_wrap(~ID_Biome, nrow = 1, strip.position = "top") +
  geom_jitter(width = .5, size=1) +
  geom_smooth(method="loess", se=F) +
  # possible values for trans : 'log2', 'log10','sqrt'
  scale_x_continuous() +
  
  labs(y="NRI", 
       x="Annual Mean Temperature", 
       title="Global scale")

NRI.AnMnTemp_Bio1.boxplot

NRI.AnMnTemp_Bio1.boxplot <- ggplot(subset(MPD.LatLong.Env.AllScales.P.05, SamplingPool = "Realm sampling"), aes(AnMnTemp_Bio1/10, NRI)) + 
  geom_jitter(width = .5, size=1) +
  geom_smooth(method="loess", se=F) +
  # possible values for trans : 'log2', 'log10','sqrt'
  scale_x_continuous() +
  
  labs(y="NRI", 
       x="Annual Mean Temperature", 
       title="Global scale")

NRI.AnMnTemp_Bio1.boxplot

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
